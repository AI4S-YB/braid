//! K-mer path graph, suffix rescue, and bottleneck paths.
//!
//! `find_paths` keeps expression in `f32` from the moment it starts, matching
//! the Cython finder. Smoothing amounts are also stored as `f32`.

use std::collections::HashMap;

use crate::graph::Graph;
use crate::splice::SpliceGraph;
use crate::util::{bisect_right, py_float};

const MIN_SCORE: f32 = 1e-10;
const SPACER: i64 = 1;
const TERMINATOR: i64 = 0;

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
enum Node {
    Source,
    Sink,
    Kmer(Vec<usize>),
}

pub(crate) struct PathGraph {
    graph: Graph<Node>,
    exprs: Vec<f64>,
    smooth_rev: Vec<f32>,
    smooth_fwd: Vec<f32>,
    k: usize,
    lost_kmers: usize,
    graph_expr: f64,
    short_expr: f64,
    short_count: usize,
    valid: bool,
}

struct Factory {
    paths: Vec<Vec<usize>>,
    exprs: Vec<f64>,
    has_source: Vec<bool>,
    has_sink: Vec<bool>,
    total_expr: f64,
    longest: usize,
    chrom: String,
    start: i64,
    end: i64,
    strand: String,
}

impl Factory {
    fn from_graph(splice: &SpliceGraph) -> Result<Self, String> {
        let mut factory = Self {
            paths: Vec::new(),
            exprs: Vec::new(),
            has_source: Vec::new(),
            has_sink: Vec::new(),
            total_expr: 0.0,
            longest: 0,
            chrom: splice.chrom.clone(),
            start: splice.start,
            end: splice.end,
            strand: splice.strand.as_gtf().to_string(),
        };
        for path in splice.transfrag_paths()? {
            factory.longest = factory.longest.max(path.nodes.len());
            factory.total_expr += path.expr;
            factory.has_source.push(path.has_source);
            factory.has_sink.push(path.has_sink);
            factory.paths.push(path.nodes);
            factory.exprs.push(path.expr);
        }
        Ok(factory)
    }

    fn create(&self, k: usize) -> Result<PathGraph, String> {
        let mut graph = PathGraph {
            graph: Graph::new(Node::Source, Node::Sink),
            exprs: vec![0.0, 0.0],
            smooth_rev: vec![0.0, 0.0],
            smooth_fwd: vec![0.0, 0.0],
            k,
            lost_kmers: 0,
            graph_expr: 0.0,
            short_expr: 0.0,
            short_count: 0,
            valid: false,
        };
        let mut short = Vec::new();
        for (index, path) in self.paths.iter().enumerate() {
            let full = self.has_source[index] && self.has_sink[index];
            let is_short = path.len() < k;
            if is_short && !full {
                short.push(index);
                graph.short_expr += self.exprs[index];
                continue;
            }
            let mut kmers = Vec::new();
            if self.has_source[index] {
                kmers.push(Node::Source);
            }
            if is_short {
                kmers.push(Node::Kmer(path.clone()));
            } else {
                for kmer in path.windows(k) {
                    kmers.push(Node::Kmer(kmer.to_vec()));
                }
            }
            if self.has_sink[index] {
                kmers.push(Node::Sink);
            }
            if kmers.is_empty() {
                continue;
            }
            graph.add_path(&kmers, self.exprs[index]);
            graph.graph_expr += self.exprs[index];
        }
        for id in graph.graph.unreachable() {
            graph.lost_kmers += 1;
            graph.graph.remove_node_id(id);
        }
        graph.short_count = short.len();
        graph.valid = graph.graph.is_valid();
        // Stash the short indexes in lost_kmers? No, rescue needs them.
        // They are recomputed by the caller from the same rule.
        let _ = short;
        Ok(graph)
    }

    fn stats(
        &self,
        graph: &PathGraph,
        kmax: usize,
        lost_short: usize,
        lost_expr: f64,
        optimal: bool,
    ) -> String {
        let frac = if self.total_expr == 0.0 {
            0.0
        } else {
            graph.graph_expr / self.total_expr
        };
        let score = (frac * graph.graph.len() as f64).round() as i64;
        [
            self.chrom.clone(),
            self.start.to_string(),
            self.end.to_string(),
            self.strand.clone(),
            graph.k.to_string(),
            kmax.to_string(),
            self.paths.len().to_string(),
            graph.short_count.to_string(),
            py_float(graph.short_expr),
            lost_short.to_string(),
            py_float(lost_expr),
            graph.graph.len().to_string(),
            graph.lost_kmers.to_string(),
            py_float(self.total_expr),
            py_float(graph.graph_expr),
            py_float(frac),
            i32::from(graph.valid).to_string(),
            score.to_string(),
            i32::from(optimal).to_string(),
        ]
        .join("\t")
    }

    fn create_optimal(&self, user_kmax: usize) -> Result<(Option<PathGraph>, Vec<String>), String> {
        if self.paths.is_empty() {
            return Ok((None, Vec::new()));
        }
        let mut kmax = self.longest.max(1);
        if user_kmax > 0 {
            kmax = user_kmax.min(kmax);
        }
        let mut lines = Vec::new();
        let objective = |k: i64| -> Result<i64, String> {
            let graph = self.create(k as usize)?;
            lines.push(self.stats(&graph, kmax, 0, 0.0, false));
            if !graph.valid {
                return Ok(-k);
            }
            let frac = if self.total_expr == 0.0 {
                0.0
            } else {
                graph.graph_expr / self.total_expr
            };
            Ok((frac * graph.graph.len() as f64).round() as i64)
        };
        let Some((k, _)) = maximize_bisect(1, kmax as i64, objective)? else {
            return Err("path graph k selection failed".to_string());
        };
        let mut graph = self.create(k as usize)?;
        let (lost_short, lost_expr) = self.rescue(&mut graph)?;
        lines.push(self.stats(&graph, kmax, lost_short, lost_expr, true));
        Ok((Some(graph), lines))
    }

    fn rescue(&self, graph: &mut PathGraph) -> Result<(usize, f64), String> {
        let mut indexes = Vec::new();
        for (index, path) in self.paths.iter().enumerate() {
            let full = self.has_source[index] && self.has_sink[index];
            if path.len() < graph.k && !full {
                indexes.push(index);
            }
        }
        if indexes.is_empty() {
            return Ok((0, 0.0));
        }
        // Score every short transfrag against the pre-rescue k-mer expression,
        // then add the whole batch. Later short transfrags do not see earlier ones.
        let index = SuffixIndex::build(graph)?;
        let mut staged = Vec::new();
        let mut lost = 0usize;
        let mut lost_expr = 0.0;
        for i in indexes {
            let mut total = 0.0;
            let mut hits = Vec::new();
            for kmer in index.search(&self.paths[i]) {
                let Some(id) = graph.graph.id_of(&Node::Kmer(kmer)) else {
                    continue;
                };
                let expr = graph.exprs[id];
                total += expr;
                hits.push((id, expr));
            }
            if hits.is_empty() {
                lost += 1;
                lost_expr += self.exprs[i];
                continue;
            }
            let count = hits.len() as f64;
            for (id, expr) in hits {
                let amount = if total == 0.0 {
                    self.exprs[i] / count
                } else {
                    self.exprs[i] * (expr / total)
                };
                staged.push((id, amount));
            }
        }
        for (id, amount) in staged {
            graph.exprs[id] += amount;
            graph.smooth_rev[id] += amount as f32;
            graph.smooth_fwd[id] += amount as f32;
        }
        Ok((lost, lost_expr))
    }
}

impl PathGraph {
    fn ensure(&mut self, id: usize) {
        while self.exprs.len() <= id {
            self.exprs.push(0.0);
            self.smooth_rev.push(0.0);
            self.smooth_fwd.push(0.0);
        }
    }

    fn add_path(&mut self, nodes: &[Node], expr: f64) {
        let ids = self.graph.add_path(nodes);
        for id in &ids {
            self.ensure(*id);
            self.exprs[*id] += expr;
        }
        if let Some(&first) = ids.first() {
            self.smooth_rev[first] += expr as f32;
        }
        if let Some(&last) = ids.last() {
            self.smooth_fwd[last] += expr as f32;
        }
    }

    fn apply_smoothing(&mut self) -> Result<(), String> {
        let order = self.graph.topo()?;
        let forward = smooth_pass(&order, &self.graph.succs, &self.exprs, &mut self.smooth_fwd);
        let mut reverse_order = order;
        reverse_order.reverse();
        let backward = smooth_pass(
            &reverse_order,
            &self.graph.preds,
            &self.exprs,
            &mut self.smooth_rev,
        );
        for i in 0..self.exprs.len() {
            self.exprs[i] += f64::from(forward[i]) + f64::from(backward[i]);
        }
        Ok(())
    }

    fn reconstruct(&self, path: &[usize]) -> Vec<usize> {
        if path.len() < 2 {
            return Vec::new();
        }
        let Some(Node::Kmer(first)) = self.graph.get(path[1]) else {
            return Vec::new();
        };
        let mut nodes = first.clone();
        for id in path.iter().take(path.len() - 1).skip(2) {
            if let Some(Node::Kmer(kmer)) = self.graph.get(*id) {
                if let Some(last) = kmer.last() {
                    nodes.push(*last);
                }
            }
        }
        nodes
    }
}

fn smooth_pass(
    order: &[usize],
    neighbors: &[Vec<usize>],
    exprs: &[f64],
    smooth: &mut [f32],
) -> Vec<f32> {
    let mut delta = vec![0.0f32; smooth.len()];
    for &node in order {
        let amount = smooth[node];
        let nbrs = &neighbors[node];
        if nbrs.is_empty() || amount == 0.0 {
            continue;
        }
        let total: f64 = nbrs.iter().map(|&next| exprs[next]).sum();
        if total == 0.0 {
            let share = f64::from(amount) / nbrs.len() as f64;
            let share_f = share as f32;
            for &next in nbrs {
                delta[next] += share_f;
                smooth[next] += share_f;
            }
        } else {
            for &next in nbrs {
                let share = (exprs[next] / total) * f64::from(amount);
                let share_f = share as f32;
                delta[next] += share_f;
                smooth[next] += share_f;
            }
        }
    }
    delta
}

struct SuffixIndex {
    text: Vec<i64>,
    sa: Vec<usize>,
    starts: Vec<i64>,
    lengths: Vec<usize>,
}

impl SuffixIndex {
    fn build(graph: &PathGraph) -> Result<Self, String> {
        let mut text = Vec::new();
        let mut starts = Vec::new();
        let mut lengths = Vec::new();
        for id in graph.graph.node_ids() {
            let Some(Node::Kmer(kmer)) = graph.graph.get(id) else {
                continue;
            };
            starts.push(text.len() as i64);
            lengths.push(kmer.len());
            text.extend(kmer.iter().map(|node| *node as i64));
            text.push(SPACER);
        }
        text.push(TERMINATOR);
        let sa = suffix_array(&text);
        Ok(Self {
            text,
            sa,
            starts,
            lengths,
        })
    }

    fn search(&self, pattern: &[usize]) -> Vec<Vec<usize>> {
        if pattern.is_empty() || self.text.is_empty() {
            return Vec::new();
        }
        let pattern: Vec<i64> = pattern.iter().map(|node| *node as i64).collect();
        let (left, right) = suffix_range(&pattern, &self.text, &self.sa);
        let mut hits = Vec::new();
        for slot in left..right {
            let pos = self.sa[slot];
            let owner = bisect_right(&self.starts, pos as i64, 0) as isize - 1;
            if owner < 0 {
                continue;
            }
            let owner = owner as usize;
            if owner >= self.starts.len() {
                continue;
            }
            let start = self.starts[owner] as usize;
            let end = start + self.lengths[owner];
            if end > self.text.len() {
                continue;
            }
            hits.push(
                self.text[start..end]
                    .iter()
                    .map(|node| *node as usize)
                    .collect(),
            );
        }
        hits
    }
}

fn suffix_cmp(pattern: &[i64], text: &[i64], pos: usize) -> i32 {
    let mut i = 0;
    while i < pattern.len() && pos + i < text.len() {
        if pattern[i] > text[pos + i] {
            return 1;
        }
        if pattern[i] < text[pos + i] {
            return -1;
        }
        i += 1;
    }
    0
}

fn suffix_range(pattern: &[i64], text: &[i64], sa: &[usize]) -> (usize, usize) {
    let n = sa.len();
    let mut left = 0;
    let mut right = n;
    while left < right {
        let mid = (left + right) / 2;
        if suffix_cmp(pattern, text, sa[mid]) > 0 {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    let start = left;
    right = n;
    while left < right {
        let mid = (left + right) / 2;
        if suffix_cmp(pattern, text, sa[mid]) < 0 {
            right = mid;
        } else {
            left = mid + 1;
        }
    }
    (start, right)
}

fn suffix_array(text: &[i64]) -> Vec<usize> {
    let n = text.len();
    if n == 0 {
        return Vec::new();
    }
    let mut rank = vec![0i32; n];
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by_key(|&i| text[i]);
    let mut class = 0i32;
    rank[order[0]] = 0;
    for pair in order.windows(2) {
        if text[pair[1]] != text[pair[0]] {
            class += 1;
        }
        rank[pair[1]] = class;
    }
    let mut sa = order;
    let mut tmp = vec![0i32; n];
    let mut width = 1usize;
    while width < n && class < n as i32 - 1 {
        sa.sort_by(|&a, &b| {
            let a2 = if a + width < n { rank[a + width] } else { -1 };
            let b2 = if b + width < n { rank[b + width] } else { -1 };
            rank[a].cmp(&rank[b]).then(a2.cmp(&b2))
        });
        tmp[sa[0]] = 0;
        class = 0;
        for i in 1..n {
            let prev = sa[i - 1];
            let cur = sa[i];
            let prev2 = if prev + width < n {
                rank[prev + width]
            } else {
                -1
            };
            let cur2 = if cur + width < n {
                rank[cur + width]
            } else {
                -1
            };
            if rank[prev] != rank[cur] || prev2 != cur2 {
                class += 1;
            }
            tmp[cur] = class;
        }
        rank.copy_from_slice(&tmp);
        width *= 2;
    }
    sa
}

/// Python 2 floor division. Ties keep the smaller `x`. `None` means the
/// iteration budget ran out.
fn maximize_bisect<F>(
    mut xmin: i64,
    mut xmax: i64,
    mut function: F,
) -> Result<Option<(i64, i64)>, String>
where
    F: FnMut(i64) -> Result<i64, String>,
{
    let max_iterations = xmax;
    let mut x = (xmin + xmax) / 2;
    let mut cache: HashMap<i64, i64> = HashMap::new();
    let mut iterations = 0i64;
    let value = |k: i64, cache: &mut HashMap<i64, i64>, function: &mut F| -> Result<i64, String> {
        if let Some(found) = cache.get(&k) {
            return Ok(*found);
        }
        let score = function(k)?;
        cache.insert(k, score);
        Ok(score)
    };
    while iterations < max_iterations {
        if xmax - xmin <= 3 {
            let mut best_y = None;
            let mut best_x = xmin;
            for x1 in xmin..=xmax {
                let y1 = value(x1, &mut cache, &mut function)?;
                if best_y.is_none_or(|best: i64| y1 > best) {
                    best_y = Some(y1);
                    best_x = x1;
                }
            }
            return Ok(Some((best_x, best_y.unwrap_or(0))));
        }
        let probes = [(x + xmin) / 2, x, (x + xmax) / 2];
        let mut best_y = None;
        let mut best_i = 0;
        for (i, probe) in probes.iter().enumerate() {
            let y1 = value(*probe, &mut cache, &mut function)?;
            if best_y.is_none_or(|best: i64| y1 > best) {
                best_y = Some(y1);
                best_i = i;
            }
        }
        if best_i == 0 {
            xmax = x;
            x = probes[0];
        } else if best_i == 1 {
            xmin = probes[0];
            xmax = probes[2];
        } else {
            xmin = x;
            x = probes[2];
        }
        iterations += 1;
    }
    Ok(None)
}

pub(crate) struct FoundPath {
    pub nodes: Vec<usize>,
    pub expr: f64,
}

pub(crate) fn assemble_paths(
    splice: &SpliceGraph,
    kmax: usize,
    path_frac: f64,
    max_paths: usize,
) -> Result<(Vec<FoundPath>, Vec<String>), String> {
    let factory = Factory::from_graph(splice)?;
    let (graph, lines) = factory.create_optimal(kmax)?;
    let Some(mut graph) = graph else {
        return Ok((Vec::new(), lines));
    };
    if graph.graph.len() == 0 || !graph.valid {
        return Ok((Vec::new(), lines));
    }
    graph.apply_smoothing()?;
    let mut found = Vec::new();
    for (kmers, expr) in find_paths(&graph, path_frac, max_paths)? {
        found.push(FoundPath {
            nodes: graph.reconstruct(&kmers),
            expr,
        });
    }
    Ok((found, lines))
}

fn find_paths(
    graph: &PathGraph,
    path_frac: f64,
    max_paths: usize,
) -> Result<Vec<(Vec<usize>, f64)>, String> {
    if graph.exprs[0] < f64::from(MIN_SCORE) {
        return Ok(Vec::new());
    }
    let order = graph.graph.topo()?;
    let mut exprs: Vec<f32> = graph.exprs.iter().map(|value| *value as f32).collect();
    let (mut path, mut expr) = find_one(&order, &mut exprs, &graph.graph.succs)?;
    let mut lowest = expr * path_frac as f32;
    if MIN_SCORE > lowest {
        lowest = MIN_SCORE;
    }
    let mut results = vec![(path, f64::from(expr))];
    let mut iterations = 1usize;
    loop {
        if max_paths > 0 && iterations >= max_paths {
            break;
        }
        (path, expr) = find_one(&order, &mut exprs, &graph.graph.succs)?;
        if expr <= lowest {
            break;
        }
        results.push((path, f64::from(expr)));
        iterations += 1;
    }
    Ok(results)
}

fn find_one(
    order: &[usize],
    exprs: &mut [f32],
    succs: &[Vec<usize>],
) -> Result<(Vec<usize>, f32), String> {
    let n = exprs.len();
    let mut min_exprs = vec![MIN_SCORE; n];
    let mut prevs = vec![1usize; n];
    min_exprs[0] = exprs[0];
    for &node in order {
        let min_expr = min_exprs[node];
        for &next in &succs[node] {
            let candidate = min_expr.min(exprs[next]);
            if prevs[next] == 1 || candidate > min_exprs[next] {
                min_exprs[next] = candidate;
                prevs[next] = node;
            }
        }
    }
    let expr = min_exprs[1];
    let mut path = vec![1usize];
    let mut prev = 1usize;
    for _ in 0..=n {
        prev = prevs[prev];
        path.push(prev);
        if prev == 0 {
            path.reverse();
            for &node in &path {
                let updated = exprs[node] - expr;
                exprs[node] = if MIN_SCORE >= updated {
                    MIN_SCORE
                } else {
                    updated
                };
            }
            return Ok((path, expr));
        }
    }
    Err("path graph traceback did not reach the source".to_string())
}

#[derive(Clone, Debug)]
pub(crate) struct Isoform {
    pub exons: Vec<(i64, i64)>,
    pub expr: f64,
    pub rel_frac: f64,
    pub abs_frac: f64,
    pub gene_id: u64,
    pub tss_id: u64,
}

struct Cluster {
    expr: f64,
    nodes: HashMap<usize, ()>,
    paths: Vec<(Vec<usize>, f64)>,
}

impl Cluster {
    fn new() -> Self {
        Self {
            expr: 1e-10,
            nodes: HashMap::new(),
            paths: Vec::new(),
        }
    }

    fn overlaps(&self, path: &[usize]) -> bool {
        path.iter().any(|node| self.nodes.contains_key(node))
    }

    fn add(&mut self, path: &[usize], expr: f64) {
        self.expr += expr;
        for node in path {
            self.nodes.insert(*node, ());
        }
        self.paths.push((path.to_vec(), expr));
    }

    fn merge(&mut self, others: Vec<Cluster>) {
        for other in others {
            self.expr += other.expr;
            self.nodes.extend(other.nodes);
            self.paths.extend(other.paths);
        }
        self.paths.sort_by(|a, b| b.1.total_cmp(&a.1));
    }
}

pub(crate) fn build_isoforms(
    paths: &[(Vec<usize>, f64)],
    min_frac: f64,
    exons_of: impl Fn(&[usize]) -> Vec<(i64, i64)>,
) -> Vec<Vec<Isoform>> {
    let mut clusters: Vec<Cluster> = Vec::new();
    for (path, expr) in paths {
        let matches: Vec<usize> = clusters
            .iter()
            .enumerate()
            .filter(|(_, cluster)| cluster.overlaps(path))
            .map(|(index, _)| index)
            .collect();
        if matches.is_empty() {
            let mut cluster = Cluster::new();
            cluster.add(path, *expr);
            clusters.push(cluster);
            continue;
        }
        let mut discard = false;
        for cluster in &clusters {
            let best = cluster.paths[0].1;
            let rel = if best == 0.0 { 0.0 } else { expr / best };
            if rel < min_frac {
                discard = true;
                break;
            }
        }
        if discard {
            continue;
        }
        let first = matches[0];
        let mut rest = Vec::new();
        for index in matches.into_iter().skip(1).rev() {
            rest.push(clusters.remove(index));
        }
        rest.reverse();
        clusters[first].merge(rest);
        clusters[first].add(path, *expr);
    }
    let mut genes = Vec::new();
    for cluster in clusters {
        let best = cluster.paths.first().map(|path| path.1).unwrap_or(0.0);
        let mut isoforms = Vec::new();
        for (path, expr) in cluster.paths {
            isoforms.push(Isoform {
                exons: exons_of(&path),
                rel_frac: if best == 0.0 { 0.0 } else { expr / best },
                abs_frac: if cluster.expr == 0.0 {
                    0.0
                } else {
                    expr / cluster.expr
                },
                expr,
                gene_id: 0,
                tss_id: 0,
            });
        }
        genes.push(isoforms);
    }
    genes
}

pub(crate) fn assign_ids(
    isoforms: &mut [Isoform],
    strand: crate::types::Strand,
    gene_id: &mut u64,
    tss_id: &mut u64,
) {
    *gene_id += 1;
    let gene = *gene_id;
    let mut tss_at: HashMap<i64, u64> = HashMap::new();
    for isoform in isoforms {
        let start = isoform.exons[0].0;
        let end = isoform.exons.last().map(|exon| exon.1).unwrap_or(start);
        let pos = if strand == crate::types::Strand::Neg {
            end
        } else {
            start
        };
        let id = *tss_at.entry(pos).or_insert_with(|| {
            *tss_id += 1;
            *tss_id
        });
        isoform.gene_id = gene;
        isoform.tss_id = id;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bisect_keeps_the_smaller_x_on_a_tie_and_uses_floor_division() {
        let (x, y) = maximize_bisect(1, 10, |k| Ok(-(k - 3) * (k - 3)))
            .unwrap()
            .unwrap();
        assert_eq!((x, y), (3, 0));
        let (x, y) = maximize_bisect(1, 4, |k| Ok(if k == 1 || k == 2 { 5 } else { 1 }))
            .unwrap()
            .unwrap();
        assert_eq!((x, y), (1, 5));
    }

    #[test]
    fn suffix_search_returns_every_matching_full_kmer() {
        let mut graph = PathGraph {
            graph: Graph::new(Node::Source, Node::Sink),
            exprs: vec![0.0, 0.0],
            smooth_rev: vec![0.0, 0.0],
            smooth_fwd: vec![0.0, 0.0],
            k: 2,
            lost_kmers: 0,
            graph_expr: 0.0,
            short_expr: 0.0,
            short_count: 0,
            valid: false,
        };
        graph.add_path(&[Node::Kmer(vec![2, 3, 2, 3])], 1.0);
        graph.add_path(&[Node::Kmer(vec![3, 4])], 1.0);
        let index = SuffixIndex::build(&graph).unwrap();
        let hits = index.search(&[2, 3]);
        assert_eq!(hits, vec![vec![2, 3, 2, 3], vec![2, 3, 2, 3]]);
        let hits = index.search(&[3, 4]);
        assert_eq!(hits, vec![vec![3, 4]]);
    }

    #[test]
    fn cluster_discards_a_weak_overlap_against_any_existing_cluster() {
        let paths = vec![
            (vec![1], 10.0),
            (vec![1], 0.1),
            (vec![2], 0.1),
            (vec![1], 0.5),
        ];
        let genes = build_isoforms(&paths, 0.05, |_| Vec::new());
        let kept: Vec<f64> = genes
            .iter()
            .flat_map(|gene| gene.iter().map(|iso| iso.expr))
            .collect();
        assert_eq!(kept, vec![10.0, 0.5, 0.1]);
    }

    #[test]
    fn negative_strand_tss_uses_the_exon_end() {
        use crate::types::Strand;
        let mut isoforms = vec![
            Isoform {
                exons: vec![(10, 20)],
                expr: 1.0,
                rel_frac: 1.0,
                abs_frac: 1.0,
                gene_id: 0,
                tss_id: 0,
            },
            Isoform {
                exons: vec![(10, 30)],
                expr: 1.0,
                rel_frac: 1.0,
                abs_frac: 1.0,
                gene_id: 0,
                tss_id: 0,
            },
            Isoform {
                exons: vec![(4, 20)],
                expr: 1.0,
                rel_frac: 1.0,
                abs_frac: 1.0,
                gene_id: 0,
                tss_id: 0,
            },
        ];
        let mut gene = 0;
        let mut tss = 0;
        assign_ids(&mut isoforms, Strand::Neg, &mut gene, &mut tss);
        assert_eq!(isoforms[0].tss_id, isoforms[2].tss_id);
        assert_ne!(isoforms[0].tss_id, isoforms[1].tss_id);
        let mut forward = isoforms.clone();
        let mut gene = 0;
        let mut tss = 0;
        assign_ids(&mut forward, Strand::Pos, &mut gene, &mut tss);
        assert_eq!(forward[0].tss_id, forward[1].tss_id);
        assert_ne!(forward[0].tss_id, forward[2].tss_id);
    }
}
