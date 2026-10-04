//! Splice graph, change-point trimming, and exon reconstruction.
//!
//! Node ids follow the order transfrags are added. Successor ids stay
//! sorted, so later graph walks do not depend on hash order.

use std::collections::{BTreeSet, HashSet};

use crate::changepoint::{find_threshold_points, run_changepoint, ChangePoint};
use crate::graph::Graph;
use crate::types::{Strand, Transfrag};
use crate::util::{bisect_left, bisect_right};

#[derive(Clone, Debug)]
pub(crate) struct FragPath {
    pub nodes: Vec<usize>,
    pub expr: f64,
    pub has_source: bool,
    pub has_sink: bool,
}

#[derive(Clone, Debug)]
pub(crate) struct SpliceGraph {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    guided_ends: bool,
    guided_assembly: bool,
    transfrags: Vec<Transfrag>,
    ref_transfrags: Vec<Transfrag>,
    graph: Graph<(i64, i64)>,
    is_start: Vec<bool>,
    is_stop: Vec<bool>,
    node_bounds: Vec<i64>,
    expr: Vec<f32>,
    ref_starts: Vec<i64>,
    ref_stops: Vec<i64>,
    /// Original `[start, start+1)` intervals. Trimming does not update them.
    tree_starts: Vec<(i64, usize)>,
    tree_ends: Vec<(i64, usize)>,
    start_sites: BTreeSet<i64>,
    stop_sites: BTreeSet<i64>,
}

impl SpliceGraph {
    pub(crate) fn create(
        transfrags: Vec<Transfrag>,
        guided_ends: bool,
        guided_assembly: bool,
    ) -> Result<Self, String> {
        let mut graph = Self {
            chrom: String::new(),
            start: 0,
            end: 0,
            strand: Strand::Na,
            guided_ends,
            guided_assembly,
            transfrags: Vec::new(),
            ref_transfrags: Vec::new(),
            graph: Graph::new((-1, -1), (-2, -2)),
            is_start: vec![false, false],
            is_stop: vec![false, false],
            node_bounds: Vec::new(),
            expr: Vec::new(),
            ref_starts: Vec::new(),
            ref_stops: Vec::new(),
            tree_starts: Vec::new(),
            tree_ends: Vec::new(),
            start_sites: BTreeSet::new(),
            stop_sites: BTreeSet::new(),
        };
        graph.rebuild(transfrags)?;
        Ok(graph)
    }

    fn rebuild(&mut self, transfrags: Vec<Transfrag>) -> Result<(), String> {
        self.transfrags.clear();
        self.ref_transfrags.clear();
        self.tree_starts.clear();
        self.tree_ends.clear();
        self.ref_starts.clear();
        self.ref_stops.clear();
        self.chrom.clear();
        let mut start = None;
        let mut end = None;
        let mut strand = None;
        let mut ref_starts = BTreeSet::new();
        let mut ref_stops = BTreeSet::new();
        for transfrag in transfrags {
            if transfrag.exons.is_empty() {
                return Err(format!("Transcript '{}' has no exons", transfrag.id));
            }
            if self.chrom.is_empty() {
                self.chrom = transfrag.chrom.clone();
            } else if self.chrom != transfrag.chrom {
                return Err("chrom mismatch".to_string());
            }
            if let Some(seen) = strand {
                if seen != transfrag.strand {
                    return Err("strand mismatch".to_string());
                }
            } else {
                strand = Some(transfrag.strand);
            }
            let transfrag_start = transfrag.start();
            let transfrag_end = transfrag.end();
            start = Some(start.map_or(transfrag_start, |value: i64| value.min(transfrag_start)));
            end = Some(end.map_or(transfrag_end, |value: i64| value.max(transfrag_end)));
            if transfrag.is_ref {
                ref_starts.insert(transfrag.tx_start());
                ref_stops.insert(transfrag.tx_stop());
                self.ref_transfrags.push(transfrag);
            } else {
                let index = self.transfrags.len();
                self.tree_starts.push((transfrag.start(), index));
                self.tree_ends.push((transfrag.end(), index));
                self.transfrags.push(transfrag);
            }
        }
        self.start = start.ok_or_else(|| "splice graph has no transfrags".to_string())?;
        self.end = end.unwrap_or(self.start);
        self.strand = strand.unwrap_or(Strand::Na);
        self.ref_starts = ref_starts.into_iter().collect();
        self.ref_stops = ref_stops.into_iter().collect();
        self.expr = vec![0.0; (self.end - self.start) as usize];
        let expressed: Vec<(f64, Vec<(i64, i64)>)> = self
            .transfrags
            .iter()
            .map(|transfrag| (transfrag.expr, transfrag.exons.clone()))
            .collect();
        for (expr, exons) in expressed {
            self.add_expression(expr, &exons);
        }
        self.node_bounds = self.boundaries();
        self.build_graph()?;
        self.mark_start_stop();
        Ok(())
    }

    fn add_expression(&mut self, expr: f64, exons: &[(i64, i64)]) {
        for &(start, end) in exons {
            let from = (start - self.start) as usize;
            let to = (end - self.start) as usize;
            for slot in self.expr.iter_mut().take(to).skip(from) {
                *slot += expr as f32;
            }
        }
    }

    fn boundaries(&self) -> Vec<i64> {
        let mut bounds = BTreeSet::from([self.start, self.end]);
        bounds.extend(self.start_sites.iter().copied());
        bounds.extend(self.stop_sites.iter().copied());
        bounds.extend(find_threshold_points(&self.expr, self.start));
        for transfrag in &self.transfrags {
            bounds.extend(transfrag.splices());
        }
        if self.guided_ends || self.guided_assembly {
            for transfrag in &self.ref_transfrags {
                if self.guided_ends {
                    bounds.insert(transfrag.start());
                    bounds.insert(transfrag.end());
                }
                if self.guided_assembly {
                    bounds.extend(transfrag.splices());
                }
            }
        }
        bounds.into_iter().collect()
    }

    fn build_graph(&mut self) -> Result<(), String> {
        let mut paths = Vec::new();
        for transfrag in self.iter_transfrags() {
            let mut nodes = split_exons(&transfrag.exons, &self.node_bounds)?;
            if self.strand == Strand::Neg {
                nodes.reverse();
            }
            paths.push(nodes);
        }
        self.graph = Graph::new((-1, -1), (-2, -2));
        self.is_start = vec![false, false];
        self.is_stop = vec![false, false];
        for nodes in paths {
            for node in &nodes {
                let id = self.graph.add_node(*node);
                if id >= self.is_start.len() {
                    self.is_start.push(false);
                    self.is_stop.push(false);
                }
            }
            self.graph.add_path(&nodes);
        }
        for id in self.graph.node_ids() {
            if self.graph.preds[id].is_empty() {
                self.is_start[id] = true;
            }
            if self.graph.succs[id].is_empty() {
                self.is_stop[id] = true;
            }
        }
        Ok(())
    }

    fn mark_start_stop(&mut self) {
        let mut starts: Vec<i64> = self.start_sites.iter().copied().collect();
        let mut stops: Vec<i64> = self.stop_sites.iter().copied().collect();
        if self.guided_ends {
            starts.extend(self.ref_starts.iter().copied());
            stops.extend(self.ref_stops.iter().copied());
        }
        let neg = self.strand == Strand::Neg;
        for pos in starts {
            // Positive-strand starts use bisect_right; negative-strand starts use bisect_left.
            self.mark_bound(pos, !neg, true);
        }
        for pos in stops {
            self.mark_bound(pos, neg, false);
        }
    }

    fn mark_bound(&mut self, pos: i64, right: bool, start: bool) {
        let index = if right {
            bisect_right(&self.node_bounds, pos, 0)
        } else {
            bisect_left(&self.node_bounds, pos, 0)
        };
        if index == 0 || index >= self.node_bounds.len() {
            return;
        }
        let node = (self.node_bounds[index - 1], self.node_bounds[index]);
        let Some(id) = self.graph.id_of(&node) else {
            return;
        };
        if start {
            self.is_start[id] = true;
        } else {
            self.is_stop[id] = true;
        }
    }

    fn iter_transfrags(&self) -> impl Iterator<Item = &Transfrag> {
        self.transfrags.iter().chain(
            self.guided_assembly
                .then_some(self.ref_transfrags.iter())
                .into_iter()
                .flatten(),
        )
    }

    pub(crate) fn transfrag_paths(&self) -> Result<Vec<FragPath>, String> {
        let mut starts = HashSet::new();
        let mut stops = HashSet::new();
        for id in self.graph.node_ids() {
            if self.is_start[id] {
                starts.insert(id);
            }
            if self.is_stop[id] {
                stops.insert(id);
            }
        }
        let mut paths = Vec::new();
        for transfrag in self.iter_transfrags() {
            let mut nodes = split_exons(&transfrag.exons, &self.node_bounds)?;
            if self.strand == Strand::Neg {
                nodes.reverse();
            }
            let ids: Vec<usize> = nodes
                .iter()
                .map(|node| {
                    self.graph
                        .id_of(node)
                        .ok_or_else(|| format!("missing splice node {node:?}"))
                })
                .collect::<Result<_, _>>()?;
            paths.push(FragPath {
                has_source: ids.first().is_some_and(|id| starts.contains(id)),
                has_sink: ids.last().is_some_and(|id| stops.contains(id)),
                nodes: ids,
                expr: transfrag.expr,
            });
        }
        Ok(paths)
    }

    pub(crate) fn detect_change_points(&self, pvalue: f64, fold_change: f64) -> Vec<ChangePoint> {
        let mut points = Vec::new();
        for id in self.graph.node_ids() {
            let Some(&(start, end)) = self.graph.get(id) else {
                continue;
            };
            let slice = self.expr_slice(start, end);
            for mut point in run_changepoint(slice, pvalue, fold_change, 20, 11) {
                point.pos += start;
                point.start += start;
                point.end += start;
                points.push(point);
            }
        }
        points
    }

    pub(crate) fn apply_change_point(&mut self, point: &ChangePoint, trim: bool) {
        let stop = (self.strand != Strand::Neg && point.sign < 0.0)
            || (self.strand == Strand::Neg && point.sign > 0.0);
        if stop {
            self.stop_sites.insert(point.pos);
        } else {
            self.start_sites.insert(point.pos);
        }
        if trim {
            self.trim(point);
        }
    }

    fn trim(&mut self, point: &ChangePoint) {
        if point.sign < 0.0 {
            let hits: Vec<usize> = self
                .tree_ends
                .iter()
                .filter(|(end, _)| overlaps(*end - 1, *end, point.pos, point.end))
                .map(|(_, index)| *index)
                .collect();
            for index in hits {
                let exons = &mut self.transfrags[index].exons;
                let Some(&(exon_start, _)) = exons.last() else {
                    continue;
                };
                if point.pos <= exon_start && exon_start <= point.end {
                    continue;
                }
                if let Some(exon) = exons.last_mut() {
                    exon.1 = point.pos;
                }
            }
        } else {
            let hits: Vec<usize> = self
                .tree_starts
                .iter()
                .filter(|(start, _)| overlaps(*start, *start + 1, point.start, point.pos))
                .map(|(_, index)| *index)
                .collect();
            for index in hits {
                let exons = &mut self.transfrags[index].exons;
                let Some(&(_, exon_end)) = exons.first() else {
                    continue;
                };
                if point.start <= exon_end && exon_end <= point.pos {
                    continue;
                }
                if let Some(exon) = exons.first_mut() {
                    exon.0 = point.pos;
                }
            }
        }
    }

    pub(crate) fn recreate(&mut self) -> Result<(), String> {
        let transfrags = self
            .transfrags
            .iter()
            .cloned()
            .chain(self.ref_transfrags.iter().cloned())
            .collect();
        self.rebuild(transfrags)
    }

    pub(crate) fn split(self) -> Result<Vec<Self>, String> {
        let component = self.graph.components();
        let mut ids: Vec<i32> = component.iter().copied().filter(|id| *id >= 0).collect();
        ids.sort_unstable();
        ids.dedup();
        if ids.len() <= 1 {
            return Ok(vec![self]);
        }
        let mut groups: Vec<Vec<Transfrag>> = vec![Vec::new(); ids.len()];
        for transfrag in self.iter_transfrags() {
            let nodes = split_exons(&transfrag.exons, &self.node_bounds)?;
            let Some(node) = nodes.first() else {
                continue;
            };
            let Some(node_id) = self.graph.id_of(node) else {
                continue;
            };
            let component_id = component[node_id];
            if component_id < 0 {
                continue;
            }
            let slot = ids.iter().position(|id| *id == component_id).unwrap();
            groups[slot].push(transfrag.clone());
        }
        // Yield components in the order their first transfrag was seen.
        let mut ordered = Vec::new();
        let mut emitted = vec![false; groups.len()];
        for transfrag in self.iter_transfrags() {
            let nodes = split_exons(&transfrag.exons, &self.node_bounds)?;
            let Some(node) = nodes.first() else {
                continue;
            };
            let Some(node_id) = self.graph.id_of(node) else {
                continue;
            };
            let component_id = component[node_id];
            let Some(slot) = ids.iter().position(|id| *id == component_id) else {
                continue;
            };
            if emitted[slot] {
                continue;
            }
            emitted[slot] = true;
            ordered.push(Self::create(
                groups[slot].clone(),
                self.guided_ends,
                self.guided_assembly,
            )?);
        }
        Ok(ordered)
    }

    pub(crate) fn node_gtf(&self) -> Vec<String> {
        let graph_id = format!(
            "G_{}_{}_{}_{}",
            self.chrom,
            self.start,
            self.end,
            self.strand.as_gtf()
        );
        let mut lines = Vec::new();
        for id in self.graph.node_ids() {
            let Some(&(start, end)) = self.graph.get(id) else {
                continue;
            };
            let slice = self.expr_slice(start, end);
            let (min, max, mean) = if slice.is_empty() {
                (0.0, 0.0, 0.0)
            } else {
                let min = slice.iter().copied().fold(f32::INFINITY, f32::min);
                let max = slice.iter().copied().fold(f32::NEG_INFINITY, f32::max);
                let mean = slice.iter().sum::<f32>() / slice.len() as f32;
                (min, max, mean)
            };
            let ref_starts = subset(&self.ref_starts, start, end);
            let ref_stops = subset(&self.ref_stops, start, end);
            lines.push(gtf_line(
                &self.chrom,
                "node",
                start,
                end,
                "0",
                self.strand.as_gtf(),
                &[
                    ("graph_id", graph_id.clone()),
                    ("expr_min", crate::changepoint::format_f32(min)),
                    ("expr_max", crate::changepoint::format_f32(max)),
                    ("expr_mean", crate::changepoint::format_f32(mean)),
                    ("ref_starts", join_ints(&ref_starts)),
                    ("ref_stops", join_ints(&ref_stops)),
                ],
            ));
        }
        lines
    }

    pub(crate) fn change_gtf(&self, point: &ChangePoint) -> Vec<String> {
        let graph_id = format!(
            "G_{}_{}_{}_{}",
            self.chrom,
            self.start,
            self.end,
            self.strand.as_gtf()
        );
        let attrs = [
            ("graph_id", graph_id),
            ("sign", py_sign(point.sign)),
            ("pvalue", crate::util::py_float(point.pvalue)),
            (
                "foldchange",
                crate::changepoint::format_f32(point.foldchange as f32),
            ),
        ];
        vec![
            gtf_line(
                &self.chrom,
                "changept",
                point.pos,
                point.pos + 1,
                "0",
                self.strand.as_gtf(),
                &attrs,
            ),
            gtf_line(
                &self.chrom,
                "changeinterval",
                point.start,
                point.end,
                "0",
                self.strand.as_gtf(),
                &attrs,
            ),
        ]
    }

    pub(crate) fn reconstruct_exons(&self, path: &[usize]) -> Vec<(i64, i64)> {
        if path.is_empty() {
            return Vec::new();
        }
        let mut ids: Vec<usize> = path.to_vec();
        if self.strand == Strand::Neg {
            ids.reverse();
        }
        let nodes: Vec<(i64, i64)> = ids
            .iter()
            .filter_map(|id| self.graph.get(*id).copied())
            .collect();
        if nodes.is_empty() {
            return Vec::new();
        }
        let mut exons = Vec::new();
        let mut chain_start = nodes[0].0;
        let mut chain_end = nodes[0].1;
        for &(start, end) in &nodes[1..] {
            if chain_end != start {
                exons.push((chain_start, chain_end));
                chain_start = start;
            }
            chain_end = end;
        }
        exons.push((chain_start, chain_end));
        exons
    }

    fn expr_slice(&self, start: i64, end: i64) -> &[f32] {
        let from = (start - self.start) as usize;
        let to = (end - self.start) as usize;
        &self.expr[from..to]
    }
}

fn overlaps(start: i64, end: i64, other_start: i64, other_end: i64) -> bool {
    start < other_end && other_start < end
}

fn subset(values: &[i64], start: i64, end: i64) -> Vec<i64> {
    let from = bisect_right(values, start, 0);
    let to = bisect_right(values, end, 0);
    values[from..to].to_vec()
}

fn join_ints(values: &[i64]) -> String {
    values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn py_sign(sign: f64) -> String {
    crate::util::py_float(sign)
}

pub(crate) fn gtf_line(
    chrom: &str,
    feature: &str,
    start: i64,
    end: i64,
    score: &str,
    strand: &str,
    attrs: &[(&str, String)],
) -> String {
    let text = attrs
        .iter()
        .map(|(key, value)| format!("{key} \"{value}\";"))
        .collect::<Vec<_>>()
        .join(" ");
    format!(
        "{chrom}\ttaco\t{feature}\t{}\t{end}\t{score}\t{strand}\t.\t{text}",
        start + 1
    )
}

fn split_exons(exons: &[(i64, i64)], bounds: &[i64]) -> Result<Vec<(i64, i64)>, String> {
    let mut nodes = Vec::new();
    let mut end_index = 0usize;
    for &(start, end) in exons {
        let start_index = bisect_right(bounds, start, end_index);
        end_index = bisect_left(bounds, end, start_index);
        if start_index == end_index {
            if start_index == 0 || start_index >= bounds.len() {
                return Err(format!("exon {start}-{end} falls outside node boundaries"));
            }
            nodes.push((bounds[start_index - 1], bounds[start_index]));
        } else {
            if start_index == 0 {
                return Err(format!("exon {start}-{end} falls outside node boundaries"));
            }
            for index in (start_index - 1)..end_index {
                if index + 1 >= bounds.len() {
                    return Err(format!("exon {start}-{end} falls outside node boundaries"));
                }
                nodes.push((bounds[index], bounds[index + 1]));
            }
        }
    }
    Ok(nodes)
}
