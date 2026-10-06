//! Capped and no-cap gene grouping, collapse voting, and transcript sort from
//! `tama_collapse.py`.
//!
//! Capped hunter ids are sorted as strings; no-cap hunters start with the
//! longest 5' end at each exon-count level. Prey walks, gene membership, merged-read
//! walks, and splice-error joins follow `Py27Dict` (CPython 2.7, seed 0).
//! A coordinate vote is the highest count, then the extreme end: smallest
//! start, largest end. That choice does not depend on coordinate-key order.

use plena_model::{Exon, Strand, Transcript};
use std::collections::{BTreeMap, HashMap};

use crate::ids::{keys, reinsert, GeneGroups};
use crate::py27::Py27Dict;
use crate::py2float::py2_str_round;
use crate::{coverage_percent, identity_percent};

#[derive(Clone, Debug)]
pub(crate) struct ReadModel {
    pub cluster_id: String,
    pub scaff: String,
    pub strand: String,
    pub start_pos: i64,
    pub end_pos: i64,
    pub exon_starts: Vec<i64>,
    pub exon_ends: Vec<i64>,
    pub sj_pre: Vec<String>,
    pub sj_post: Vec<String>,
    pub seq_length: i64,
    pub h_count: i64,
    pub s_count: i64,
    pub i_count: i64,
    pub d_count: i64,
    pub mis_count: i64,
    pub polya_seq: String,
    pub a_count: i64,
    /// Fraction of A, not yet multiplied by 100.
    pub a_percent: f64,
}

impl ReadModel {
    pub(crate) fn num_exons(&self) -> usize {
        self.exon_starts.len()
    }

    fn coverage(&self) -> f64 {
        coverage_percent(self.seq_length, self.h_count, self.s_count)
    }

    fn identity(&self, method: &str) -> Result<f64, String> {
        identity_percent(
            self.seq_length,
            self.h_count,
            self.s_count,
            self.i_count,
            self.d_count,
            self.mis_count,
            method,
        )
    }
}

#[derive(Clone, Debug)]
pub(crate) struct GroupParams {
    pub capped: bool,
    pub five_prime: i64,
    pub exon_diff: i64,
    pub three_prime: i64,
    pub ends: String,
    pub duplicates: String,
    pub sj_priority: bool,
    pub ident_method: String,
}

#[derive(Clone, Debug)]
struct CollapseResult {
    starts: Vec<i64>,
    ends: Vec<i64>,
    start_wobble: Vec<i64>,
    end_wobble: Vec<i64>,
    sj_start: Vec<i64>,
    sj_end: Vec<i64>,
    start_nuc: Vec<String>,
    end_nuc: Vec<String>,
}

#[derive(Clone, Debug)]
struct Merged {
    trans_id: String,
    scaff: String,
    strand: String,
    reads: Py27Dict<String, ReadModel>,
    num_exons: usize,
    start_pos: i64,
    end_pos: i64,
    collapse_starts: Vec<i64>,
    collapse_ends: Vec<i64>,
    start_wobble: Vec<i64>,
    end_wobble: Vec<i64>,
    sj_start: Vec<i64>,
    sj_end: Vec<i64>,
    start_nuc: Vec<String>,
    end_nuc: Vec<String>,
}

impl Merged {
    fn new(trans_id: String) -> Self {
        Self {
            trans_id,
            scaff: String::new(),
            strand: "none".to_string(),
            reads: Py27Dict::new(),
            num_exons: 0,
            start_pos: 0,
            end_pos: 0,
            collapse_starts: Vec::new(),
            collapse_ends: Vec::new(),
            start_wobble: Vec::new(),
            end_wobble: Vec::new(),
            sj_start: Vec::new(),
            sj_end: Vec::new(),
            start_nuc: Vec::new(),
            end_nuc: Vec::new(),
        }
    }

    fn add_read(&mut self, read: ReadModel) -> Result<(), String> {
        if self.num_exons < read.num_exons() {
            self.num_exons = read.num_exons();
        }
        self.scaff = read.scaff.clone();
        if self.strand == "none" {
            self.strand = read.strand.clone();
        } else if self.strand != read.strand {
            return Err(format!(
                "merged transcripts are not on the same strand ({})",
                read.cluster_id
            ));
        }
        self.reads.insert(read.cluster_id.clone(), read);
        Ok(())
    }

    fn apply(&mut self, collapsed: CollapseResult) -> Result<(), String> {
        if collapsed.starts.is_empty() {
            return Err(format!("collapse produced no exons for {}", self.trans_id));
        }
        self.start_pos = collapsed.starts[0];
        self.end_pos = *collapsed.ends.last().unwrap();
        self.collapse_starts = collapsed.starts;
        self.collapse_ends = collapsed.ends;
        self.start_wobble = collapsed.start_wobble;
        self.end_wobble = collapsed.end_wobble;
        self.sj_start = collapsed.sj_start;
        self.sj_end = collapsed.sj_end;
        self.start_nuc = collapsed.start_nuc;
        self.end_nuc = collapsed.end_nuc;
        Ok(())
    }

    fn bed_line(&self) -> Result<String, String> {
        Self::bed_blocks(
            &self.scaff,
            &self.strand,
            &format!(
                "{};{}",
                self.trans_id.split('.').next().unwrap_or(""),
                self.trans_id
            ),
            (self.start_pos, self.end_pos),
            &self.collapse_starts,
            &self.collapse_ends,
            self.num_exons,
        )
    }

    fn read_bed_line(read: &ReadModel, final_trans_id: &str) -> Result<String, String> {
        Self::bed_blocks(
            &read.scaff,
            &read.strand,
            &format!("{final_trans_id};{}", read.cluster_id),
            (read.start_pos, read.end_pos),
            &read.exon_starts,
            &read.exon_ends,
            read.num_exons(),
        )
    }

    fn bed_blocks(
        scaff: &str,
        strand: &str,
        name: &str,
        (start_pos, end_pos): (i64, i64),
        starts: &[i64],
        ends: &[i64],
        num_exons: usize,
    ) -> Result<String, String> {
        if starts.len() < num_exons || ends.len() < num_exons {
            return Err(format!("exon list shorter than num_exons for {name}"));
        }
        let mut lengths = Vec::with_capacity(num_exons);
        let mut relatives = Vec::with_capacity(num_exons);
        for i in 0..num_exons {
            let exon_start = starts[i];
            let exon_end = ends[i];
            let relative = exon_start - start_pos;
            if relative < 0 {
                return Err(format!("negative relative exon start for {name}"));
            }
            lengths.push((exon_end - exon_start).to_string());
            relatives.push(relative.to_string());
        }
        Ok(format!(
            "{scaff}\t{}\t{}\t{name}\t40\t{strand}\t{}\t{}\t255,0,0\t{num_exons}\t{}\t{}",
            start_pos - 1,
            end_pos - 1,
            start_pos - 1,
            end_pos - 1,
            lengths.join(","),
            relatives.join(",")
        ))
    }

    fn trans_report_line(&self, ident_method: &str) -> Result<String, String> {
        let mut high_q = None;
        let mut low_q = None;
        let mut high_c = None;
        let mut low_c = None;
        for (_, read) in self.reads.iter() {
            let quality = read.identity(ident_method)?;
            let coverage = read.coverage();
            high_q = Some(high_q.map_or(quality, |v: f64| v.max(quality)));
            low_q = Some(low_q.map_or(quality, |v: f64| v.min(quality)));
            high_c = Some(high_c.map_or(coverage, |v: f64| v.max(coverage)));
            low_c = Some(low_c.map_or(coverage, |v: f64| v.min(coverage)));
        }
        let (Some(high_c), Some(low_c), Some(high_q), Some(low_q)) = (high_c, low_c, high_q, low_q)
        else {
            return Err(format!("transcript {} has no reads", self.trans_id));
        };
        let start_wobble = join_ints(&self.start_wobble);
        let end_wobble = join_ints(&self.end_wobble);
        let mut sj_start: Vec<String> = self.sj_start.iter().map(|n| n.to_string()).collect();
        let mut sj_end: Vec<String> = self.sj_end.iter().map(|n| n.to_string()).collect();
        if self.strand == "+" {
            sj_start.reverse();
            sj_end.reverse();
        }
        let mut nucs = self.start_nuc.clone();
        if nucs.len() > 1 {
            if self.strand == "+" {
                if nucs.last().map(String::as_str) == Some("na") {
                    nucs.pop();
                    nucs.reverse();
                } else {
                    return Err(format!(
                        "plus collapse_error_nuc_list does not end in na ({})",
                        self.trans_id
                    ));
                }
            } else if self.strand == "-" {
                if nucs.first().map(String::as_str) == Some("na") {
                    nucs.remove(0);
                } else {
                    return Err(format!(
                        "minus collapse_error_nuc_list does not start with na ({})",
                        self.trans_id
                    ));
                }
            } else {
                return Err(format!("bad strand {} in trans report", self.strand));
            }
        }
        Ok(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            self.trans_id,
            self.reads.len(),
            py2_str_round(high_c),
            py2_str_round(low_c),
            py2_str_round(high_q),
            py2_str_round(low_q),
            start_wobble,
            end_wobble,
            sj_start.join(","),
            sj_end.join(","),
            nucs.join(";")
        ))
    }
}

fn join_ints(values: &[i64]) -> String {
    values
        .iter()
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn py_get<T: Clone>(values: &[T], index: i64) -> Result<T, String> {
    let resolved = if index < 0 {
        values.len() as i64 + index
    } else {
        index
    };
    if resolved < 0 || resolved >= values.len() as i64 {
        return Err(format!(
            "index {index} outside list of length {}",
            values.len()
        ));
    }
    Ok(values[resolved as usize].clone())
}

fn fuzzy(coord1: i64, coord2: i64, threshold: i64) -> (&'static str, i64) {
    if coord1 == coord2 {
        ("perfect_match", 0)
    } else {
        let diff = coord1 - coord2;
        if diff.abs() <= threshold {
            ("wobbly_match", diff)
        } else {
            ("no_match", diff)
        }
    }
}

/// Capped comparison. Different exon counts and non-overlapping exons are
/// different transcripts. `no_cap` is not this path.
fn same_transcript(a: &ReadModel, b: &ReadModel, params: &GroupParams) -> Result<bool, String> {
    if !params.capped {
        return Ok(same_nocap(a, b, params));
    }
    if a.num_exons() != b.num_exons() {
        return Ok(false);
    }
    let n = a.num_exons();
    for i in 0..n {
        let j = if a.strand == "+" {
            -1 - i as i64
        } else if a.strand == "-" {
            i as i64
        } else {
            return Err(format!("bad strand {}", a.strand));
        };
        let a_start = py_get(&a.exon_starts, j)?;
        let b_start = py_get(&b.exon_starts, j)?;
        let a_end = py_get(&a.exon_ends, j)?;
        let b_end = py_get(&b.exon_ends, j)?;
        if a_start >= b_end || b_start >= a_end {
            return Ok(false);
        }
        let mut start_threshold = params.exon_diff;
        let mut end_threshold = params.exon_diff;
        if a.strand == "+" {
            if i == 0 {
                end_threshold = params.three_prime;
            }
            if i + 1 == n {
                start_threshold = params.five_prime;
            }
        } else if i == 0 {
            start_threshold = params.three_prime;
            if i + 1 == n {
                end_threshold = params.five_prime;
            }
        } else if i + 1 == n {
            end_threshold = params.five_prime;
        }
        let (start_flag, _) = fuzzy(a_start, b_start, start_threshold);
        let (end_flag, _) = fuzzy(a_end, b_end, end_threshold);
        if start_flag == "no_match" || end_flag == "no_match" {
            return Ok(false);
        }
    }
    Ok(true)
}

fn priority_start(pre: &str, post: &str) -> i64 {
    match (post == "0", pre == "0") {
        (true, true) => 0,
        (true, false) => 1,
        (false, true) => 2,
        (false, false) => 3,
    }
}

fn priority_end(pre: &str, post: &str) -> i64 {
    match (post == "0", pre == "0") {
        (true, true) => 0,
        (false, true) => 1,
        (true, false) => 2,
        (false, false) => 3,
    }
}

fn sj_priority(
    read: &ReadModel,
    i: usize,
    max_exon_num: usize,
    use_priority: bool,
) -> Result<(i64, i64, String, String), String> {
    let mut start_pri = 0i64;
    let mut end_pri = 0i64;
    let mut start_err = "na".to_string();
    let mut end_err = "na".to_string();
    if read.num_exons() > 1 {
        if read.strand == "+" {
            if i == 0 {
                let pre = py_get(&read.sj_pre, -1)?;
                let post = py_get(&read.sj_post, -1)?;
                start_pri = priority_start(&pre, &post);
                start_err = format!("{pre}>{post}");
            } else if i + 1 < max_exon_num {
                let pre_s = py_get(&read.sj_pre, -1 - i as i64)?;
                let post_s = py_get(&read.sj_post, -1 - i as i64)?;
                let pre_e = py_get(&read.sj_pre, -(i as i64))?;
                let post_e = py_get(&read.sj_post, -(i as i64))?;
                start_pri = priority_start(&pre_s, &post_s);
                end_pri = priority_end(&pre_e, &post_e);
                start_err = format!("{pre_s}>{post_s}");
                end_err = format!("{pre_e}>{post_e}");
            } else if i + 1 == max_exon_num {
                let pre_e = py_get(&read.sj_pre, -(i as i64))?;
                let post_e = py_get(&read.sj_post, -(i as i64))?;
                end_pri = priority_end(&pre_e, &post_e);
                end_err = format!("{pre_e}>{post_e}");
            } else {
                return Err(format!(
                    "plus priority index {i} of {max_exon_num} for {}",
                    read.cluster_id
                ));
            }
        } else if read.strand == "-" {
            if i == 0 {
                let pre = py_get(&read.sj_pre, 0)?;
                let post = py_get(&read.sj_post, 0)?;
                end_pri = priority_end(&pre, &post);
                end_err = format!("{pre}>{post}");
            } else if i + 1 < max_exon_num {
                let pre_s = py_get(&read.sj_pre, i as i64 - 1)?;
                let post_s = py_get(&read.sj_post, i as i64 - 1)?;
                let pre_e = py_get(&read.sj_pre, i as i64)?;
                let post_e = py_get(&read.sj_post, i as i64)?;
                start_pri = priority_start(&pre_s, &post_s);
                end_pri = priority_end(&pre_e, &post_e);
                start_err = format!("{pre_s}>{post_s}");
                end_err = format!("{pre_e}>{post_e}");
            } else if i + 1 == max_exon_num {
                let pre_s = py_get(&read.sj_pre, i as i64 - 1)?;
                let post_s = py_get(&read.sj_post, i as i64 - 1)?;
                start_pri = priority_start(&pre_s, &post_s);
                start_err = format!("{pre_s}>{post_s}");
            } else {
                return Err(format!(
                    "minus priority index {i} of {max_exon_num} for {}",
                    read.cluster_id
                ));
            }
        } else {
            return Err(format!("bad strand {}", read.strand));
        }
    }
    if !use_priority {
        start_pri = 0;
        end_pri = 0;
    }
    Ok((start_pri, end_pri, start_err, end_err))
}

fn collapse_transcripts(
    reads: &[ReadModel],
    params: &GroupParams,
) -> Result<CollapseResult, String> {
    if reads.is_empty() {
        return Err("collapse called with no transcripts".to_string());
    }
    let strand = &reads[0].strand;
    let max_exon_num = reads.iter().map(ReadModel::num_exons).max().unwrap();
    let mut collapse_starts = Vec::new();
    let mut collapse_ends = Vec::new();
    let mut start_wobble = Vec::new();
    let mut end_wobble = Vec::new();
    let mut sj_start = Vec::new();
    let mut sj_end = Vec::new();
    let mut start_nuc = Vec::new();
    let mut end_nuc = Vec::new();

    for i in 0..max_exon_num {
        let j = if strand == "+" {
            -1 - i as i64
        } else if strand == "-" {
            i as i64
        } else {
            return Err(format!("bad strand {strand}"));
        };
        let mut start_counts: BTreeMap<i64, BTreeMap<i64, i64>> = BTreeMap::new();
        let mut end_counts: BTreeMap<i64, BTreeMap<i64, i64>> = BTreeMap::new();
        let mut start_errors: BTreeMap<i64, BTreeMap<i64, Py27Dict<String, i32>>> = BTreeMap::new();
        let mut end_errors: BTreeMap<i64, BTreeMap<i64, Py27Dict<String, i32>>> = BTreeMap::new();
        let mut start_range = Vec::new();
        let mut end_range = Vec::new();

        for read in reads {
            if read.strand != *strand {
                return Err("collapse strand mismatch".to_string());
            }
            if i >= read.num_exons() {
                continue;
            }
            let this_max = read.num_exons();
            let (start_pri, end_pri, start_err, end_err) =
                sj_priority(read, i, this_max, params.sj_priority)?;
            let e_start = py_get(&read.exon_starts, j)?;
            let e_end = py_get(&read.exon_ends, j)?;
            let truncated_five = !params.capped && i + 1 == this_max && i + 1 < max_exon_num;
            if !(truncated_five && strand == "+") {
                vote(
                    &mut start_counts,
                    &mut start_errors,
                    &mut start_range,
                    start_pri,
                    e_start,
                    start_err,
                );
            }
            if !(truncated_five && strand == "-") {
                vote(
                    &mut end_counts,
                    &mut end_errors,
                    &mut end_range,
                    end_pri,
                    e_end,
                    end_err,
                );
            }
        }

        let best_start_pri = *start_counts
            .keys()
            .next()
            .ok_or("empty start priority list")?;
        let best_end_pri = *end_counts.keys().next().ok_or("empty end priority list")?;
        sj_start.push(best_start_pri);
        sj_end.push(best_end_pri);

        let start_bucket = &start_counts[&best_start_pri];
        let end_bucket = &end_counts[&best_end_pri];
        let (mut best_start, long_start) = choose_start(start_bucket);
        let (mut best_end, long_end) = choose_end(end_bucket);

        if !params.capped && i + 1 == max_exon_num {
            if strand == "+" {
                best_start = long_start;
            } else {
                best_end = long_end;
            }
        }
        if params.ends == "longest_ends" {
            if i + 1 == max_exon_num {
                if strand == "+" {
                    best_start = long_start;
                } else {
                    best_end = long_end;
                }
            }
            if i == 0 {
                if strand == "+" {
                    best_end = long_end;
                } else {
                    best_start = long_start;
                }
            }
        }
        if best_start > best_end {
            if long_start < long_end {
                best_start = long_start;
                best_end = long_end;
            } else {
                return Err(format!(
                    "collapse start {best_start} is past end {best_end}"
                ));
            }
        }
        collapse_starts.push(best_start);
        collapse_ends.push(best_end);
        start_wobble.push(wobble(&mut start_range)?);
        end_wobble.push(wobble(&mut end_range)?);
        start_nuc.push(join_error_keys(
            start_errors
                .get(&best_start_pri)
                .and_then(|m| m.get(&best_start))
                .ok_or("missing start error strings")?,
        ));
        end_nuc.push(join_error_keys(
            end_errors
                .get(&best_end_pri)
                .and_then(|m| m.get(&best_end))
                .ok_or("missing end error strings")?,
        ));
    }

    let (collapse_starts, start_wobble) = zip_sort(collapse_starts, start_wobble);
    let (collapse_ends, end_wobble) = zip_sort(collapse_ends, end_wobble);
    check_exon_order(&collapse_starts, &collapse_ends)?;
    Ok(CollapseResult {
        starts: collapse_starts,
        ends: collapse_ends,
        start_wobble,
        end_wobble,
        sj_start,
        sj_end,
        start_nuc,
        end_nuc,
    })
}

fn vote(
    counts: &mut BTreeMap<i64, BTreeMap<i64, i64>>,
    errors: &mut BTreeMap<i64, BTreeMap<i64, Py27Dict<String, i32>>>,
    range: &mut Vec<i64>,
    priority: i64,
    coord: i64,
    error: String,
) {
    *counts
        .entry(priority)
        .or_default()
        .entry(coord)
        .or_insert(0) += 1;
    errors
        .entry(priority)
        .or_default()
        .entry(coord)
        .or_default()
        .insert(error, 1);
    range.push(coord);
}

fn choose_start(bucket: &BTreeMap<i64, i64>) -> (i64, i64) {
    let most = bucket.values().copied().max().unwrap_or(0);
    let best = bucket
        .iter()
        .filter(|(_, count)| **count == most)
        .map(|(coord, _)| *coord)
        .min()
        .unwrap_or(-1);
    let long = bucket.keys().copied().min().unwrap_or(-1);
    (best, long)
}

fn choose_end(bucket: &BTreeMap<i64, i64>) -> (i64, i64) {
    let most = bucket.values().copied().max().unwrap_or(0);
    let best = bucket
        .iter()
        .filter(|(_, count)| **count == most)
        .map(|(coord, _)| *coord)
        .max()
        .unwrap_or(-1);
    let long = bucket.keys().copied().max().unwrap_or(-1);
    (best, long)
}

fn wobble(range: &mut [i64]) -> Result<i64, String> {
    if range.is_empty() {
        return Err("empty wobble range".to_string());
    }
    range.sort();
    Ok(range[range.len() - 1] - range[0])
}

fn join_error_keys(errors: &Py27Dict<String, i32>) -> String {
    errors
        .iter()
        .map(|(key, _)| key.as_str())
        .collect::<Vec<_>>()
        .join("-")
}

fn zip_sort(coords: Vec<i64>, wobbles: Vec<i64>) -> (Vec<i64>, Vec<i64>) {
    let mut pairs: Vec<(i64, i64)> = coords.into_iter().zip(wobbles).collect();
    pairs.sort_by_key(|pair| pair.0);
    (
        pairs.iter().map(|pair| pair.0).collect(),
        pairs.iter().map(|pair| pair.1).collect(),
    )
}

fn check_exon_order(starts: &[i64], ends: &[i64]) -> Result<(), String> {
    let mut prev_start = -1i64;
    let mut prev_end = -1i64;
    for i in 0..starts.len() {
        if ends[i] <= starts[i] {
            return Err(format!(
                "exon end {} is not after start {}",
                ends[i], starts[i]
            ));
        }
        if starts[i] <= prev_start || ends[i] <= prev_end {
            return Err("collapsed exons are not strictly increasing".to_string());
        }
        prev_start = starts[i];
        prev_end = ends[i];
    }
    Ok(())
}

fn solo_collapse(read: &ReadModel, params: &GroupParams) -> Result<CollapseResult, String> {
    let mut starts = read.exon_starts.clone();
    let mut ends = read.exon_ends.clone();
    starts.sort();
    ends.sort();
    let n = starts.len();
    let zeros = vec![0; n];
    let mut sj_start = Vec::with_capacity(n);
    let mut sj_end = Vec::with_capacity(n);
    let mut start_nuc = Vec::with_capacity(n);
    let mut end_nuc = Vec::with_capacity(n);
    for i in 0..n {
        let (start_pri, end_pri, start_err, end_err) = sj_priority(read, i, n, params.sj_priority)?;
        sj_start.push(start_pri);
        sj_end.push(end_pri);
        start_nuc.push(start_err);
        end_nuc.push(end_err);
    }
    Ok(CollapseResult {
        starts,
        ends,
        start_wobble: zeros.clone(),
        end_wobble: zeros,
        sj_start,
        sj_end,
        start_nuc,
        end_nuc,
    })
}

struct TransGroup {
    group_count: i64,
    trans_group: Py27Dict<String, Py27Dict<i64, i32>>,
    group_trans: Py27Dict<i64, Py27Dict<String, i32>>,
}

impl TransGroup {
    fn add_to(&mut self, a: &str, b: &str) -> Result<(), String> {
        let groups: Vec<i64> = self
            .trans_group
            .get_str(a)
            .ok_or("missing short transcript")?
            .iter()
            .map(|(g, _)| *g)
            .collect();
        if groups.len() == 1
            && self
                .group_trans
                .get(&groups[0])
                .is_some_and(|g| g.len() == 1)
        {
            self.group_trans.remove(&groups[0]);
            self.trans_group.remove(&a.to_string());
        }
        let dest: Vec<i64> = self
            .trans_group
            .get_str(b)
            .ok_or("missing long transcript")?
            .iter()
            .map(|(g, _)| *g)
            .collect();
        for g in dest {
            self.trans_group
                .or_insert_with(a.to_string(), Py27Dict::new)
                .insert(g, 1);
            self.group_trans
                .get_mut(&g)
                .ok_or("missing destination group")?
                .insert(a.to_string(), 1);
        }
        Ok(())
    }
    fn new() -> Self {
        Self {
            group_count: 0,
            trans_group: Py27Dict::new(),
            group_trans: Py27Dict::new(),
        }
    }

    fn has_trans(&self, trans_id: &str) -> bool {
        self.trans_group.get_str(trans_id).is_some()
    }

    fn same_group(&self, a: &str, b: &str) -> bool {
        let Some(a_groups) = self.trans_group.get_str(a) else {
            return false;
        };
        let Some(b_groups) = self.trans_group.get_str(b) else {
            return false;
        };
        a_groups.iter().any(|(group, _)| b_groups.contains(group))
    }

    fn new_group(&mut self, trans_id: &str) -> Result<(), String> {
        if self.has_trans(trans_id) {
            return Err(format!("transcript {trans_id} is already grouped"));
        }
        self.group_count += 1;
        let group = self.group_count;
        self.trans_group
            .or_insert_with(trans_id.to_string(), Py27Dict::new)
            .insert(group, 1);
        self.group_trans
            .or_insert_with(group, Py27Dict::new)
            .insert(trans_id.to_string(), 1);
        Ok(())
    }

    fn only_group(&self, trans_id: &str) -> Result<i64, String> {
        let groups = self
            .trans_group
            .get_str(trans_id)
            .ok_or_else(|| format!("transcript {trans_id} has no group"))?;
        let mut found = None;
        for (group, _) in groups.iter() {
            if found.is_some() {
                return Err(format!(
                    "capped transcript {trans_id} is in multiple groups"
                ));
            }
            found = Some(*group);
        }
        found.ok_or_else(|| format!("transcript {trans_id} has an empty group map"))
    }

    fn merge(&mut self, trans_a: &str, trans_b: &str) -> Result<(), String> {
        let a_group = self.only_group(trans_a)?;
        let b_group = self.only_group(trans_b)?;
        if a_group == b_group {
            return Err(format!(
                "groups for {trans_a} and {trans_b} are already the same"
            ));
        }
        let a_members = keys(self.group_trans.get(&a_group).ok_or("missing group A")?);
        let b_members = keys(self.group_trans.get(&b_group).ok_or("missing group B")?);
        let (keep, drop, incoming) = if a_members.len() > b_members.len() {
            (a_group, b_group, b_members)
        } else {
            (b_group, a_group, a_members)
        };
        for member in incoming {
            self.group_trans
                .get_mut(&keep)
                .ok_or("missing survivor group")?
                .insert(member.clone(), 1);
            if let Some(groups) = self.trans_group.get_str_mut(&member) {
                groups.remove(&drop);
                groups.insert(keep, 1);
            }
        }
        self.group_trans.remove(&drop);
        Ok(())
    }
}

fn simplify_capped(
    reads: &[ReadModel],
    params: &GroupParams,
) -> Result<Py27Dict<i64, Py27Dict<String, i32>>, String> {
    let mut index = HashMap::new();
    let mut ungrouped: Py27Dict<String, i32> = Py27Dict::new();
    for (i, read) in reads.iter().enumerate() {
        index.insert(read.cluster_id.clone(), i);
        ungrouped.insert(read.cluster_id.clone(), 1);
    }
    let mut unsearched: Py27Dict<String, i32> = Py27Dict::new();
    let mut ungrouped_count = reads.len();
    let mut unsearched_count = 0usize;
    let mut hunter_id = String::new();
    let mut groups = TransGroup::new();

    while ungrouped_count > 0 {
        if unsearched_count == 0 {
            let mut ids = keys(&ungrouped);
            ids.sort();
            hunter_id = ids
                .into_iter()
                .next()
                .ok_or("ungrouped transcript list is empty")?
                .clone();
            ungrouped.remove(&hunter_id);
            unsearched_count = 1;
        }
        while unsearched_count > 0 && ungrouped_count > 0 {
            if hunter_id == "new_hunter" {
                let mut ids = keys(&unsearched);
                ids.sort();
                hunter_id = ids
                    .into_iter()
                    .next()
                    .ok_or("unsearched transcript list is empty")?
                    .clone();
                unsearched.remove(&hunter_id);
            }
            let hunter = &reads[*index
                .get(&hunter_id)
                .ok_or_else(|| format!("missing hunter {hunter_id}"))?];
            if !groups.has_trans(&hunter_id) {
                groups.new_group(&hunter_id)?;
            }
            let prey_ids = keys(&ungrouped);
            for prey_id in prey_ids {
                if prey_id == hunter_id {
                    continue;
                }
                let prey = &reads[*index
                    .get(&prey_id)
                    .ok_or_else(|| format!("missing prey {prey_id}"))?];
                if hunter.strand != prey.strand {
                    return Err("transcripts inside one gene have different strands".to_string());
                }
                if !groups.has_trans(&prey_id) {
                    groups.new_group(&prey_id)?;
                }
                if groups.same_group(&hunter_id, &prey_id) {
                    continue;
                }
                if same_transcript(hunter, prey, params)? {
                    groups.merge(&hunter_id, &prey_id)?;
                    unsearched.insert(prey_id, 1);
                }
            }
            for id in keys(&unsearched) {
                ungrouped.remove(&id);
            }
            unsearched_count = unsearched.len();
            hunter_id = "new_hunter".to_string();
            ungrouped_count = ungrouped.len();
        }
    }
    Ok(groups.group_trans)
}

fn same_nocap(a: &ReadModel, b: &ReadModel, p: &GroupParams) -> bool {
    let na = a.num_exons();
    let nb = b.num_exons();
    let n = na.min(nb);
    for i in 0..n {
        let ia = if a.strand == "+" { na - i - 1 } else { i };
        let ib = if a.strand == "+" { nb - i - 1 } else { i };
        let (sa, ea) = (a.exon_starts[ia], a.exon_ends[ia]);
        let (sb, eb) = (b.exon_starts[ib], b.exon_ends[ib]);
        if sa >= eb || sb >= ea {
            return false;
        }
        let st = if a.strand == "-" && i == 0 {
            p.three_prime
        } else if a.strand == "+" && na == nb && i + 1 == n {
            p.five_prime
        } else {
            p.exon_diff
        };
        let et = if a.strand == "+" && i == 0 {
            p.three_prime
        } else if a.strand == "-" && na == nb && i + 1 == n {
            p.five_prime
        } else {
            p.exon_diff
        };
        let sm = (sa - sb).abs() <= st;
        let em = (ea - eb).abs() <= et;
        if i + 1 == n {
            if a.strand == "+" {
                if !em {
                    return false;
                }
                if !sm && na != nb && !((na < nb && sa > sb) || (nb < na && sb > sa)) {
                    return false;
                }
            } else {
                if !sm {
                    return false;
                }
                if !em && na != nb && !((na < nb && ea < eb) || (nb < na && eb < ea)) {
                    return false;
                }
            }
        } else if !sm || !em {
            return false;
        }
    }
    true
}

fn simplify_nocap(
    reads: &[ReadModel],
    p: &GroupParams,
) -> Result<Py27Dict<i64, Py27Dict<String, i32>>, String> {
    let lookup: HashMap<&str, &ReadModel> =
        reads.iter().map(|r| (r.cluster_id.as_str(), r)).collect();
    let mut levels: BTreeMap<usize, Py27Dict<String, i32>> = BTreeMap::new();
    for r in reads {
        levels
            .entry(r.num_exons())
            .or_default()
            .insert(r.cluster_id.clone(), 1);
    }
    let mut groups = TransGroup::new();
    let mut degraded = std::collections::HashSet::new();
    let mut pending: Py27Dict<String, i32> = Py27Dict::new();
    for (&level, ids) in levels.iter().rev() {
        let mut ungrouped = reinsert(ids);
        while !ungrouped.is_empty() {
            let mut by_five: BTreeMap<i64, Py27Dict<String, i32>> = BTreeMap::new();
            for (id, _) in ungrouped.iter() {
                let r = lookup[id.as_str()];
                let coord = if r.strand == "+" {
                    r.start_pos
                } else {
                    -r.end_pos
                };
                by_five.entry(coord).or_default().insert(id.clone(), 1);
            }
            let mut hunter = by_five
                .first_key_value()
                .unwrap()
                .1
                .iter()
                .next()
                .unwrap()
                .0
                .clone();
            ungrouped.remove(&hunter);
            loop {
                if !groups.has_trans(&hunter) {
                    groups.new_group(&hunter)?;
                }
                let h = lookup[hunter.as_str()];
                if !degraded.contains(&hunter) {
                    for (prey, _) in ungrouped.iter() {
                        if prey == &hunter {
                            continue;
                        }
                        if !groups.has_trans(prey) {
                            groups.new_group(prey)?;
                        }
                        if !groups.same_group(&hunter, prey)
                            && same_nocap(h, lookup[prey.as_str()], p)
                        {
                            groups.add_to(prey, &hunter)?;
                            pending.insert(prey.clone(), 1);
                        }
                    }
                }
                for (_, smaller) in levels.range(..level).rev() {
                    let candidates = reinsert(smaller);
                    for (prey, _) in candidates.iter() {
                        if !groups.has_trans(prey) {
                            groups.new_group(prey)?;
                        }
                        if !groups.same_group(&hunter, prey)
                            && same_nocap(h, lookup[prey.as_str()], p)
                        {
                            groups.add_to(prey, &hunter)?;
                            degraded.insert(prey.clone());
                        }
                    }
                }
                let mut queued = keys(&pending);
                queued.sort();
                let Some(next) = queued.first() else {
                    break;
                };
                hunter = next.clone();
                pending.remove(&hunter);
            }
        }
    }
    Ok(groups.group_trans)
}

/// Only the geometry needed for exon-overlap grouping, borrowed from either
/// a collapsed read or a merge model. Coordinates remain algorithm-local.
pub(crate) struct TranscriptGeometry<'a> {
    pub id: &'a str,
    pub starts: &'a [i64],
    pub ends: &'a [i64],
}

pub(crate) trait GroupedTranscript {
    fn geometry(&self) -> TranscriptGeometry<'_>;
}

impl GroupedTranscript for ReadModel {
    fn geometry(&self) -> TranscriptGeometry<'_> {
        TranscriptGeometry {
            id: &self.cluster_id,
            starts: &self.exon_starts,
            ends: &self.exon_ends,
        }
    }
}

/// Group one strand. Returned transcript ids under each start are in the
/// dict order `process_loci` later walks.
pub(crate) fn gene_group<T: GroupedTranscript>(
    trans_ids: &[String],
    reads: &HashMap<String, T>,
) -> Result<(GeneGroups, Vec<i64>), String> {
    let models: Vec<TranscriptGeometry<'_>> = trans_ids
        .iter()
        .map(|id| {
            reads
                .get(id)
                .map(GroupedTranscript::geometry)
                .ok_or_else(|| format!("missing transcript {id}"))
        })
        .collect::<Result<_, _>>()?;
    let mut gene_trans: Py27Dict<i64, Py27Dict<String, i32>> = Py27Dict::new();
    let mut trans_gene: HashMap<String, i64> = HashMap::new();
    let mut gene_start: HashMap<i64, i64> = HashMap::new();
    let mut gene_count = 0i64;

    if models.len() == 1 {
        gene_count += 1;
        let id = models[0].id.to_string();
        let start = models[0].starts[0];
        gene_start.insert(gene_count, start);
        gene_trans
            .or_insert_with(gene_count, Py27Dict::new)
            .insert(id.clone(), 1);
        trans_gene.insert(id, gene_count);
    }

    for i in 0..models.len() {
        for j in (i + 1)..models.len() {
            let left = &models[i];
            let right = &models[j];
            let left_id = left.id;
            let right_id = right.id;
            if let (Some(gene), Some(other)) = (trans_gene.get(left_id), trans_gene.get(right_id)) {
                if gene == other {
                    continue;
                }
            }
            let overlap = exons_overlap(left, right);
            if !overlap {
                if !trans_gene.contains_key(left_id) {
                    gene_count += 1;
                    trans_gene.insert(left_id.to_string(), gene_count);
                    gene_trans
                        .or_insert_with(gene_count, Py27Dict::new)
                        .insert(left_id.to_string(), 1);
                    gene_start.insert(gene_count, left.starts[0]);
                }
                if !trans_gene.contains_key(right_id) {
                    gene_count += 1;
                    trans_gene.insert(right_id.to_string(), gene_count);
                    gene_trans
                        .or_insert_with(gene_count, Py27Dict::new)
                        .insert(right_id.to_string(), 1);
                    gene_start.insert(gene_count, right.starts[0]);
                }
            } else if !trans_gene.contains_key(left_id) && !trans_gene.contains_key(right_id) {
                gene_count += 1;
                trans_gene.insert(left_id.to_string(), gene_count);
                trans_gene.insert(right_id.to_string(), gene_count);
                let members = gene_trans.or_insert_with(gene_count, Py27Dict::new);
                members.insert(left_id.to_string(), 1);
                members.insert(right_id.to_string(), 1);
                gene_start.insert(gene_count, left.starts[0].min(right.starts[0]));
            } else if !trans_gene.contains_key(left_id) {
                let gene_num = *trans_gene.get(right_id).unwrap();
                trans_gene.insert(left_id.to_string(), gene_num);
                gene_trans
                    .get_mut(&gene_num)
                    .ok_or("missing gene while adding left transcript")?
                    .insert(left_id.to_string(), 1);
                gene_start.insert(gene_num, left.starts[0].min(right.starts[0]));
            } else if !trans_gene.contains_key(right_id) {
                let gene_num = *trans_gene.get(left_id).unwrap();
                trans_gene.insert(right_id.to_string(), gene_num);
                gene_trans
                    .get_mut(&gene_num)
                    .ok_or("missing gene while adding right transcript")?
                    .insert(right_id.to_string(), 1);
                gene_start.insert(gene_num, left.starts[0].min(right.starts[0]));
            } else {
                let gene_num = *trans_gene.get(left_id).unwrap();
                let other_num = *trans_gene.get(right_id).unwrap();
                if gene_num != other_num {
                    let moving = keys(
                        gene_trans
                            .get(&other_num)
                            .ok_or("missing gene during merge")?,
                    );
                    for member in moving {
                        trans_gene.insert(member.clone(), gene_num);
                        gene_trans
                            .get_mut(&gene_num)
                            .ok_or("missing survivor gene")?
                            .insert(member, 1);
                    }
                    gene_trans.remove(&other_num);
                    let kept = *gene_start.get(&gene_num).ok_or("missing gene start")?;
                    let other = *gene_start
                        .get(&other_num)
                        .ok_or("missing other gene start")?;
                    gene_start.insert(gene_num, kept.min(other));
                    gene_start.remove(&other_num);
                }
            }
        }
    }

    let mut start_gene: HashMap<i64, i64> = HashMap::new();
    for (gene_num, start) in &gene_start {
        if start_gene.contains_key(start) {
            return Err(format!("multiple genes start at {start}"));
        }
        start_gene.insert(*start, *gene_num);
    }
    let mut starts: Vec<i64> = start_gene.keys().copied().collect();
    starts.sort();
    let mut by_start: Py27Dict<i64, Py27Dict<String, i32>> = Py27Dict::new();
    for start in &starts {
        let gene_num = start_gene[start];
        let members = gene_trans
            .get(&gene_num)
            .ok_or("missing gene while building start map")?;
        let dest = by_start.or_insert_with(*start, Py27Dict::new);
        for (trans_id, _) in members.iter() {
            dest.insert(trans_id.clone(), 1);
        }
    }
    Ok((by_start, starts))
}

fn exons_overlap(left: &TranscriptGeometry<'_>, right: &TranscriptGeometry<'_>) -> bool {
    for (a_start, a_end) in left.starts.iter().zip(left.ends) {
        for (b_start, b_end) in right.starts.iter().zip(right.ends) {
            if *a_start <= *b_end && *a_end >= *b_start {
                return true;
            }
        }
    }
    false
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum Tok {
    Text(String),
    Pad,
}

impl Tok {
    fn as_int(&self) -> Result<i64, String> {
        match self {
            Tok::Pad => Ok(0),
            Tok::Text(text) => text
                .parse()
                .map_err(|_| format!("bad transcript sort field {text}")),
        }
    }
}

fn position_line(starts: &[i64], ends: &[i64]) -> Result<String, String> {
    if starts.is_empty() || ends.is_empty() {
        return Err("transcript has no exons to sort".to_string());
    }
    let mut starts = starts.to_vec();
    let mut ends = ends.to_vec();
    starts.sort();
    ends.sort();
    let mut fields = vec![starts[0].to_string(), ends[ends.len() - 1].to_string()];
    for i in 0..starts.len() {
        fields.push(starts[i].to_string());
        fields.push(ends[i].to_string());
    }
    Ok(fields.join(","))
}

fn sort_transcripts(models: Vec<Merged>, params: &GroupParams) -> Result<Vec<Merged>, String> {
    let mut by_key: HashMap<String, Merged> = HashMap::new();
    let mut keys = Vec::new();
    for mut model in models {
        model.collapse_starts.sort();
        model.collapse_ends.sort();
        let key = position_line(&model.collapse_starts, &model.collapse_ends)?;
        match by_key.entry(key) {
            std::collections::hash_map::Entry::Occupied(mut entry) => {
                if params.duplicates != "merge_dup" {
                    return Err("duplicate collapsed models and -d is not merge_dup".to_string());
                }
                let incoming: Vec<ReadModel> =
                    model.reads.iter().map(|(_, read)| read.clone()).collect();
                let survivor = entry.get_mut();
                for read in incoming {
                    survivor.add_read(read)?;
                }
                let members: Vec<ReadModel> = survivor
                    .reads
                    .iter()
                    .map(|(_, read)| read.clone())
                    .collect();
                let collapsed = collapse_transcripts(&members, params)?;
                survivor.apply(collapsed)?;
            }
            std::collections::hash_map::Entry::Vacant(entry) => {
                keys.push(entry.key().clone());
                entry.insert(model);
            }
        }
    }
    let order = sort_position_keys(&keys)?;
    let mut sorted = Vec::with_capacity(order.len());
    for key in order {
        sorted.push(
            by_key
                .remove(&key)
                .ok_or_else(|| format!("missing sorted transcript {key}"))?,
        );
    }
    Ok(sorted)
}

pub(crate) fn sort_position_keys(keys: &[String]) -> Result<Vec<String>, String> {
    let width = keys
        .iter()
        .map(|key| key.split(',').count())
        .max()
        .unwrap_or(0);
    let mut rows: Vec<(String, Vec<Tok>)> = Vec::with_capacity(keys.len());
    for key in keys {
        let mut toks: Vec<Tok> = key
            .split(',')
            .map(|field| Tok::Text(field.to_string()))
            .collect();
        while toks.len() < width {
            toks.push(Tok::Pad);
        }
        rows.push((key.clone(), toks));
    }
    rows.sort_by(|left, right| {
        left.1[0]
            .as_int()
            .unwrap_or(0)
            .cmp(&right.1[0].as_int().unwrap_or(0))
    });
    iterate_sort(&mut rows, 0)?;
    Ok(rows.into_iter().map(|(key, _)| key).collect())
}

fn iterate_sort(rows: &mut [(String, Vec<Tok>)], col: usize) -> Result<(), String> {
    let mut groups: Vec<Vec<usize>> = Vec::new();
    let mut open: Option<usize> = None;
    for j in 1..rows.len() {
        if rows[j].1[col] == rows[j - 1].1[col] {
            if let Some(index) = open {
                groups[index].push(j);
            } else {
                open = Some(groups.len());
                groups.push(vec![j - 1, j]);
            }
        } else {
            open = None;
        }
    }
    for indices in groups {
        if indices.is_empty() {
            continue;
        }
        let mut previous = indices[0];
        for index in indices.iter().skip(1) {
            if *index != previous + 1 {
                return Err("transcript sort ties are not consecutive".to_string());
            }
            previous = *index;
        }
        let min = indices[0];
        let max = indices[indices.len() - 1] + 1;
        if col + 1 >= rows[min].1.len() {
            return Err("transcript sort ran out of columns".to_string());
        }
        rows[min..max].sort_by(|left, right| {
            left.1[col + 1]
                .as_int()
                .unwrap_or(0)
                .cmp(&right.1[col + 1].as_int().unwrap_or(0))
        });
        iterate_sort(&mut rows[min..max], col + 1)?;
    }
    Ok(())
}

pub(crate) struct WrittenGene {
    pub transcripts: Vec<Transcript>,
    pub bed: Vec<String>,
    pub trans_report: Vec<String>,
    pub trans_read: Vec<String>,
    pub polya: Vec<String>,
}

pub(crate) fn collapse_locus(
    ordered_ids: &[String],
    reads: &HashMap<String, ReadModel>,
    params: &GroupParams,
    mut gene_count: i64,
) -> Result<(i64, Vec<WrittenGene>), String> {
    let mut forward = Vec::new();
    let mut reverse = Vec::new();
    for id in ordered_ids {
        let read = reads
            .get(id)
            .ok_or_else(|| format!("missing grouped read {id}"))?;
        if read.strand == "+" {
            forward.push(id.clone());
        } else if read.strand == "-" {
            reverse.push(id.clone());
        } else {
            return Err(format!("bad strand on {id}"));
        }
    }
    let (forward_genes, forward_starts) = gene_group(&forward, reads)?;
    let (reverse_genes, _reverse_starts) = gene_group(&reverse, reads)?;
    let mut starts: Vec<i64> = forward_starts;
    for (start, _) in reverse_genes.iter() {
        if !starts.contains(start) {
            starts.push(*start);
        }
    }
    starts.sort();
    let mut written = Vec::new();
    for start in starts {
        let mut bundles: Vec<Vec<ReadModel>> = Vec::new();
        if let Some(members) = forward_genes.get(&start) {
            bundles.push(models_in_dict_order(members, reads)?);
        }
        if let Some(members) = reverse_genes.get(&start) {
            bundles.push(models_in_dict_order(members, reads)?);
        }
        for bundle in bundles {
            gene_count += 1;
            written.push(collapse_gene(&bundle, params, gene_count)?);
        }
    }
    Ok((gene_count, written))
}

fn models_in_dict_order(
    members: &Py27Dict<String, i32>,
    reads: &HashMap<String, ReadModel>,
) -> Result<Vec<ReadModel>, String> {
    let mut models = Vec::new();
    for (id, _) in members.iter() {
        models.push(
            reads
                .get(id)
                .ok_or_else(|| format!("missing gene member {id}"))?
                .clone(),
        );
    }
    Ok(models)
}

fn collapse_gene(
    reads: &[ReadModel],
    params: &GroupParams,
    gene_count: i64,
) -> Result<WrittenGene, String> {
    let groups = if params.capped {
        simplify_capped(reads, params)?
    } else {
        simplify_nocap(reads, params)?
    };
    let mut merged_models = Vec::new();
    let mut tmp_count = 0i64;
    for (_, members) in groups.iter() {
        tmp_count += 1;
        let mut merged = Merged::new(format!("G{gene_count}.tmp.{tmp_count}"));
        let mut member_reads = Vec::new();
        for (id, _) in members.iter() {
            let read = reads
                .iter()
                .find(|read| read.cluster_id == *id)
                .ok_or_else(|| format!("missing collapsed member {id}"))?
                .clone();
            member_reads.push(read.clone());
            merged.add_read(read)?;
        }
        let collapsed = if member_reads.len() > 1 {
            collapse_transcripts(&member_reads, params)?
        } else {
            solo_collapse(&member_reads[0], params)?
        };
        merged.apply(collapsed)?;
        merged_models.push(merged);
    }
    let sorted = sort_transcripts(merged_models, params)?;
    let mut transcripts = Vec::new();
    let mut bed = Vec::new();
    let mut trans_report = Vec::new();
    let mut trans_read = Vec::new();
    let mut polya = Vec::new();
    for (index, merged) in sorted.into_iter().enumerate() {
        let final_id = format!("G{gene_count}.{}", index + 1);
        let mut merged = merged;
        merged.trans_id = final_id.clone();
        bed.push(merged.bed_line()?);
        transcripts.push(Transcript {
            chrom: merged.scaff.clone(),
            strand: if merged.strand == "+" {
                Strand::Forward
            } else {
                Strand::Reverse
            },
            exons: merged
                .collapse_starts
                .iter()
                .zip(&merged.collapse_ends)
                .map(|(&start, &end)| Exon { start, end })
                .collect(),
            gene_id: format!("G{gene_count}"),
            transcript_id: final_id.clone(),
            source: None,
            score: None,
            cds_start: None,
            cds_end: None,
        });
        trans_report.push(merged.trans_report_line(&params.ident_method)?);
        for (_, read) in merged.reads.iter() {
            trans_read.push(Merged::read_bed_line(read, &final_id)?);
            let percent = read.a_percent * 100.0;
            if percent > 70.0 {
                polya.push(format!(
                    "{}\t{}\t{}\t{}\t{}\t{}",
                    read.cluster_id,
                    final_id,
                    read.strand,
                    py2_str_round(percent),
                    read.a_count,
                    read.polya_seq
                ));
            }
        }
    }
    Ok(WrittenGene {
        transcripts,
        bed,
        trans_report,
        trans_read,
        polya,
    })
}
