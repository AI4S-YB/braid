//! Sub-exon update from read coverage and splice junctions.
//!
//! This follows `get_non_overlapping_exons` and `update_exons`. Coverage has
//! to clear both a fraction of the peak and an absolute floor of 20. Junction
//! ends within 10 bp are merged toward the better-supported site. Sharp
//! coverage changes, after a Gaussian smooth, can split an exon further.

use crate::gtf::{partition_exons, Gene, Isoform};

#[derive(Clone, Debug)]
pub struct Coverage {
    pub origin: i64,
    pub diff: Vec<i32>,
    pub junctions: Vec<(i64, i64)>,
}

impl Coverage {
    pub fn new(start: i64, end: i64) -> Self {
        let len = (end - start).max(0) as usize;
        Self {
            origin: start,
            diff: vec![0; len + 1],
            junctions: Vec::new(),
        }
    }

    pub fn observe(
        &mut self,
        read_start: i64,
        read_end: i64,
        blocks: &[(i64, i64)],
        junctions: &[(i64, i64)],
    ) {
        let gene_end = self.origin + self.diff.len() as i64 - 1;
        if read_start < self.origin || read_end >= gene_end {
            return;
        }
        for &(start, end) in blocks {
            let from = start - self.origin;
            let to = end - self.origin;
            if from < 0
                || to < 0
                || to as usize >= self.diff.len()
                || from as usize >= self.diff.len()
            {
                continue;
            }
            self.diff[from as usize] += 1;
            self.diff[to as usize] -= 1;
        }
        self.junctions.extend(junctions.iter().copied());
    }
}

pub fn refine_gene(
    gene: &mut Gene,
    coverage: &Coverage,
    exon_fraction: f64,
    splice_fraction: f64,
    z_score: f64,
) {
    let called = call_exons(coverage, exon_fraction, splice_fraction, z_score);
    if called.is_empty() {
        return;
    }
    let updated = update_exons(&called, &gene.exons, 20);
    if updated.is_empty() {
        return;
    }
    let old = gene.exons.clone();
    let isoforms = std::mem::take(&mut gene.isoforms);
    gene.exons = updated;
    gene.isoforms = remap_isoforms(&old, &gene.exons, isoforms);
}

fn call_exons(
    coverage: &Coverage,
    exon_fraction: f64,
    splice_fraction: f64,
    z_score: f64,
) -> Vec<(i64, i64)> {
    let size = coverage.diff.len().saturating_sub(1);
    if size == 0 {
        return Vec::new();
    }
    let mut depth = Vec::with_capacity(size);
    let mut acc = 0i32;
    let mut peak = 0i32;
    for &delta in coverage.diff.iter().take(size) {
        acc += delta;
        depth.push(acc);
        peak = peak.max(acc);
    }
    if peak <= 0 {
        return Vec::new();
    }
    let floor = ((peak as f64) * exon_fraction).max(20.0);
    let mut positions = Vec::new();
    for (offset, &cov) in depth.iter().enumerate() {
        if cov as f64 > floor {
            positions.push(coverage.origin + offset as i64);
        }
    }
    let mut blocks = continuous(&positions);
    blocks.retain(|(start, end)| end - start >= 20);
    blocks = fill_holes(blocks, 20);
    if blocks.is_empty() {
        return Vec::new();
    }
    let exons = split_on_junctions(&blocks, &coverage.junctions, splice_fraction);
    cut_by_derivative(&exons, &depth, coverage.origin, z_score)
}

fn continuous(positions: &[i64]) -> Vec<(i64, i64)> {
    if positions.is_empty() {
        return Vec::new();
    }
    let mut blocks = Vec::new();
    let mut start = positions[0];
    let mut end = positions[0] + 1;
    for &pos in &positions[1..] {
        if pos == end {
            end = pos + 1;
        } else {
            blocks.push((start, end));
            start = pos;
            end = pos + 1;
        }
    }
    blocks.push((start, end));
    blocks
}

fn fill_holes(blocks: Vec<(i64, i64)>, hole: i64) -> Vec<(i64, i64)> {
    if blocks.len() < 2 {
        return blocks;
    }
    let mut out = Vec::new();
    let (mut start, mut end) = blocks[0];
    for &(next, next_end) in &blocks[1..] {
        if next - end > hole {
            out.push((start, end));
            start = next;
            end = next_end;
        } else {
            end = next_end;
        }
    }
    out.push((start, end));
    out
}

fn split_on_junctions(
    blocks: &[(i64, i64)],
    junctions: &[(i64, i64)],
    fraction: f64,
) -> Vec<(i64, i64)> {
    let mut freq: Vec<(i64, i64)> = Vec::new();
    for &(start, end) in junctions {
        add_count(&mut freq, start);
        add_count(&mut freq, end);
    }
    if freq.len() <= 1 {
        return blocks.to_vec();
    }
    let peak = freq.iter().map(|item| item.1).max().unwrap_or(0);
    let floor = (peak as f64 * fraction).max(20.0);
    let mut exons = Vec::new();
    for &(start, end) in blocks {
        let mut cuts: Vec<(i64, i64)> = freq
            .iter()
            .copied()
            .filter(|(pos, count)| *count as f64 >= floor && *pos > start && *pos < end)
            .collect();
        cuts = merge_boundaries(cuts, 10);
        cuts.retain(|(pos, _)| (*pos - start).abs() > 10 && (end - *pos).abs() > 10);
        cuts.sort_by_key(|item| item.0);
        let mut points = vec![start];
        points.extend(cuts.into_iter().map(|item| item.0));
        points.push(end);
        for pair in points.windows(2) {
            if pair[1] > pair[0] {
                exons.push((pair[0], pair[1]));
            }
        }
    }
    exons
}

fn add_count(freq: &mut Vec<(i64, i64)>, pos: i64) {
    if let Some(slot) = freq.iter_mut().find(|item| item.0 == pos) {
        slot.1 += 1;
    } else {
        freq.push((pos, 1));
    }
}

fn merge_boundaries(mut boundaries: Vec<(i64, i64)>, distance: i64) -> Vec<(i64, i64)> {
    boundaries.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    let mut merged = Vec::new();
    while let Some((pos, mut value)) = boundaries.first().copied() {
        boundaries.remove(0);
        let mut drop = Vec::new();
        for (index, &(other, other_value)) in boundaries.iter().enumerate() {
            if (other - pos).abs() <= distance {
                value += other_value;
                drop.push(index);
            }
        }
        for index in drop.into_iter().rev() {
            boundaries.remove(index);
        }
        merged.push((pos, value));
    }
    merged
}

pub(crate) fn cut_by_derivative(
    exons: &[(i64, i64)],
    depth: &[i32],
    origin: i64,
    z_score: f64,
) -> Vec<(i64, i64)> {
    let mut values = Vec::new();
    let mut positions = Vec::new();
    for &(start, end) in exons {
        if end - start < 2 {
            continue;
        }
        let series: Vec<f64> = (start..end)
            .map(|pos| {
                let index = (pos - origin) as usize;
                if index < depth.len() {
                    depth[index] as f64
                } else {
                    0.0
                }
            })
            .collect();
        let smooth = gaussian(&series, 1.0);
        for (offset, pair) in smooth.windows(2).enumerate() {
            values.push(pair[1] - pair[0]);
            positions.push(start + 1 + offset as i64);
        }
    }
    if values.is_empty() {
        return exons.to_vec();
    }
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let var = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / values.len() as f64;
    let std = var.sqrt();
    if std == 0.0 {
        return exons.to_vec();
    }
    let mut cuts: Vec<i64> = values
        .iter()
        .zip(&positions)
        .filter_map(|(value, pos)| (((value - mean) / std).abs() > z_score).then_some(*pos))
        .collect();
    if cuts.is_empty() {
        return exons.to_vec();
    }
    let mut kept = vec![cuts[0]];
    for pos in cuts.drain(1..) {
        if pos - *kept.last().unwrap() > 10 {
            kept.push(pos);
        }
    }
    let mut next = Vec::new();
    for &(start, end) in exons {
        let local: Vec<i64> = kept
            .iter()
            .copied()
            .filter(|pos| pos - start > 10 && end - pos > 10)
            .collect();
        if local.is_empty() {
            next.push((start, end));
            continue;
        }
        let mut points = vec![start];
        points.extend(local);
        points.push(end);
        for pair in points.windows(2) {
            next.push((pair[0], pair[1]));
        }
    }
    next
}

fn gaussian(values: &[f64], sigma: f64) -> Vec<f64> {
    if values.is_empty() {
        return Vec::new();
    }
    let radius = (4.0 * sigma).ceil() as i64;
    let mut kernel = Vec::new();
    let mut sum = 0.0;
    for offset in -radius..=radius {
        let weight = (-0.5 * (offset as f64 / sigma).powi(2)).exp();
        kernel.push(weight);
        sum += weight;
    }
    for weight in &mut kernel {
        *weight /= sum;
    }
    let n = values.len() as i64;
    (0..n)
        .map(|center| {
            let mut acc = 0.0;
            for (slot, &weight) in kernel.iter().enumerate() {
                let at = center + slot as i64 - radius;
                acc += weight * values[reflect(at, n)];
            }
            acc
        })
        .collect()
}

fn reflect(mut index: i64, n: i64) -> usize {
    if n <= 1 {
        return 0;
    }
    let period = 2 * (n - 1);
    index = index.rem_euclid(period);
    if index >= n {
        index = period - index;
    }
    index as usize
}

fn update_exons(
    discovered: &[(i64, i64)],
    reference: &[(i64, i64)],
    distance: i64,
) -> Vec<(i64, i64)> {
    if discovered.is_empty() {
        return reference.to_vec();
    }
    if reference.len() == 1 {
        return reference.to_vec();
    }
    let points: Vec<i64> = reference.iter().flat_map(|exon| [exon.0, exon.1]).collect();
    let corrected: Vec<(i64, i64)> = discovered
        .iter()
        .filter_map(|&(start, end)| {
            let start = snap(start, &points, distance);
            let end = snap(end, &points, distance);
            (end > start).then_some((start, end))
        })
        .collect();
    partition_union(&corrected, reference)
}

fn snap(point: i64, reference: &[i64], distance: i64) -> i64 {
    let mut best: Option<(i64, i64)> = None;
    for &candidate in reference {
        let gap = (point - candidate).abs();
        if gap < distance && best.is_none_or(|(have, _)| gap < have) {
            best = Some((gap, candidate));
        }
    }
    best.map(|(_, candidate)| candidate).unwrap_or(point)
}

fn partition_union(discovered: &[(i64, i64)], reference: &[(i64, i64)]) -> Vec<(i64, i64)> {
    let mut all = Vec::new();
    all.extend(discovered.iter().copied());
    all.extend(reference.iter().copied());
    all.sort_unstable();
    if all.is_empty() {
        return Vec::new();
    }
    let mut parts = Vec::new();
    let (mut current_start, mut current_end) = all[0];
    for &(start, end) in &all[1..] {
        if start >= current_end {
            if current_start < current_end {
                parts.push((current_start, current_end));
            }
            current_start = start;
            current_end = end;
        } else {
            let mut temp = [current_start, current_end, start, end];
            temp.sort_unstable();
            if temp[0] < temp[1] {
                parts.push((temp[0], temp[1]));
            }
            if temp[1] < temp[2] {
                parts.push((temp[1], temp[2]));
            }
            current_start = temp[2];
            current_end = temp[3];
        }
    }
    if current_start < current_end {
        parts.push((current_start, current_end));
    }
    let covered = partition_exons(&parts);
    if covered.is_empty() {
        parts
    } else {
        covered
    }
}

fn remap_isoforms(old: &[(i64, i64)], new: &[(i64, i64)], isoforms: Vec<Isoform>) -> Vec<Isoform> {
    let mut mapped = Vec::new();
    for mut isoform in isoforms {
        let mut indexes = Vec::new();
        for &old_index in &isoform.exons {
            let Some(&(start, end)) = old.get(old_index) else {
                continue;
            };
            for (new_index, &(next, next_end)) in new.iter().enumerate() {
                if start <= next && end >= next_end {
                    indexes.push(new_index);
                }
            }
        }
        indexes.sort_unstable();
        indexes.dedup();
        if indexes.is_empty() {
            continue;
        }
        isoform.exons = indexes;
        mapped.push(isoform);
    }
    mapped.sort_by(|a, b| {
        b.exons
            .len()
            .cmp(&a.exons.len())
            .then(a.order.cmp(&b.order))
    });
    mapped
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flat_coverage_keeps_one_exon() {
        let mut coverage = Coverage::new(0, 100);
        for _ in 0..30 {
            coverage.observe(0, 99, &[(0, 100)], &[]);
        }
        let exons = call_exons(&coverage, 0.02, 0.02, 10.0);
        assert_eq!(exons, vec![(0, 100)]);
    }
}
