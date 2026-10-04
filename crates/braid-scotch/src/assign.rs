//! Read-to-exon and read-to-isoform compatibility.
//!
//! A sub-exon is matched when the read covers at least `high` of its length,
//! skipped when it covers at most `low`, and ignored in between. The outer
//! end whose covered bases are greater is treated as a possible truncation.
//! Known isoforms are one-hot over sub-exons (`1` included, `-1` excluded).
//! A read conflicts with an isoform when one requires an exon the other skips.

use crate::gtf::Gene;

#[derive(Clone, Debug)]
pub struct MatchParams {
    pub low: f64,
    pub high: f64,
    pub small_exon: i64,
    pub small_exon_high: i64,
    pub truncation_match: f64,
}

impl Default for MatchParams {
    fn default() -> Self {
        Self {
            low: 0.1,
            high: 0.6,
            small_exon: 0,
            small_exon_high: 80,
            truncation_match: 0.4,
        }
    }
}

#[derive(Clone, Debug)]
pub enum Outcome {
    Known(usize),
    Novel(Vec<i8>),
    Uncategorized,
}

#[derive(Clone, Debug)]
pub struct Placement {
    pub gene: usize,
    pub outcome: Outcome,
    pub score: f64,
    pub priority: f64,
}

pub fn place_in_metagene(
    positions: &[i64],
    read_start: i64,
    read_end: i64,
    genes: &[Gene],
    indexes: &[usize],
    params: &MatchParams,
    poly: &[bool],
) -> Option<Placement> {
    let mut scored = Vec::new();
    for (slot, &gene_index) in indexes.iter().enumerate() {
        let gene = &genes[gene_index];
        if gene.exons.is_empty() {
            continue;
        }
        let counts = exon_base_counts(positions, &gene.exons);
        let bases: u32 = counts.iter().sum();
        let span = read_end.max(gene.end) - read_start.min(gene.start);
        scored.push((
            gene_index,
            counts,
            bases,
            span,
            gene.end - gene.start,
            poly[slot],
        ));
    }
    if scored.is_empty() {
        return None;
    }
    let any_exon = scored.iter().any(|item| item.2 > 0);
    if !any_exon {
        scored.sort_by(|a, b| b.3.cmp(&a.3).then(a.0.cmp(&b.0)));
        let (gene_index, counts, _, _, _, poly_flag) = scored[0].clone();
        return Some(finish(
            gene_index,
            &genes[gene_index],
            &counts,
            params,
            poly_flag,
        ));
    }
    scored.retain(|item| item.2 > 0);
    scored.sort_by(|a, b| {
        b.2.cmp(&a.2)
            .then(b.3.cmp(&a.3))
            .then(b.4.cmp(&a.4))
            .then(a.0.cmp(&b.0))
    });
    let total_bases: u32 = scored.iter().map(|item| item.2).sum();
    if total_bases < 2 || scored.len() == 1 {
        let (gene_index, counts, _, _, _, poly_flag) = scored[0].clone();
        return Some(finish(
            gene_index,
            &genes[gene_index],
            &counts,
            params,
            poly_flag,
        ));
    }
    let mut ranked = Vec::new();
    for (gene_index, counts, bases, span, length, poly_flag) in scored {
        let placement = finish(gene_index, &genes[gene_index], &counts, params, poly_flag);
        let rank = match placement.outcome {
            Outcome::Known(_) => 0,
            Outcome::Novel(_) => 1,
            Outcome::Uncategorized => 2,
        };
        ranked.push((rank, bases, span, length, placement));
    }
    ranked.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then(b.1.cmp(&a.1))
            .then(b.2.cmp(&a.2))
            .then(b.3.cmp(&a.3))
            .then(a.4.gene.cmp(&b.4.gene))
    });
    ranked.into_iter().next().map(|item| item.4)
}

fn finish(
    gene_index: usize,
    gene: &Gene,
    counts: &[u32],
    params: &MatchParams,
    poly: bool,
) -> Placement {
    let (outcome, score) = assign_gene(gene, counts, params, poly);
    let priority = match &outcome {
        Outcome::Known(_) => score,
        Outcome::Novel(_) => -2.0,
        Outcome::Uncategorized => -3.0,
    };
    Placement {
        gene: gene_index,
        outcome,
        score,
        priority,
    }
}

fn assign_gene(gene: &Gene, counts: &[u32], params: &MatchParams, poly: bool) -> (Outcome, f64) {
    let bases: u32 = counts.iter().sum();
    if bases == 0 {
        return (Outcome::Uncategorized, -1.0);
    }
    let lengths: Vec<i64> = gene.exons.iter().map(|(a, b)| b - a).collect();
    let mut pct: Vec<f64> = counts
        .iter()
        .zip(&lengths)
        .map(|(&count, &len)| {
            if count == 0 || len <= 0 {
                0.0
            } else {
                round2(count as f64 / len as f64)
            }
        })
        .collect();
    let adjusted = adjust_truncation(&pct, &lengths, gene.strand, poly, params.truncation_match);
    pct = adjusted;
    let mut vector = exon_vector(&pct, params.high, params.low, poly, gene.strand);
    if vector.iter().all(|value| *value == 0) {
        if let Some((index, _)) = pct
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1).then(b.0.cmp(&a.0)))
        {
            pct[index] = 1.0;
            vector = exon_vector(&pct, params.high, params.low, poly, gene.strand);
        }
    }
    let novel_vector = vector.clone();
    let mut compatible = known_compatible(gene, &vector, &lengths, params.small_exon);
    if compatible.iter().filter(|&&on| on).count() == 0 {
        let mean = mean(&lengths);
        let cap = params
            .small_exon
            .max(params.small_exon_high.min(mean.floor() as i64));
        let mut threshold = params.small_exon + 10;
        while threshold <= cap {
            compatible = known_compatible(gene, &vector, &lengths, threshold);
            if compatible.iter().any(|&on| on) {
                break;
            }
            threshold += 10;
        }
    }
    let chosen = resolve_compatible(gene, &compatible, &pct, params.low, params.high);
    if let Some(index) = chosen {
        let score = isoform_score(gene, index, bases);
        (Outcome::Known(index), score)
    } else if bases > 0 {
        (Outcome::Novel(novel_vector), -1.0)
    } else {
        (Outcome::Uncategorized, -1.0)
    }
}

fn known_compatible(gene: &Gene, vector: &[i8], lengths: &[i64], small: i64) -> Vec<bool> {
    let masked: Vec<i8> = vector
        .iter()
        .zip(lengths)
        .map(|(&value, &len)| if len >= small { value } else { 0 })
        .collect();
    gene.isoforms
        .iter()
        .map(|isoform| !conflicts(&isoform.exons, gene.exons.len(), &masked))
        .collect()
}

fn conflicts(exons: &[usize], n_exons: usize, read: &[i8]) -> bool {
    let mut included = vec![false; n_exons];
    for &exon in exons {
        if exon < n_exons {
            included[exon] = true;
        }
    }
    for (exon, &state) in read.iter().enumerate() {
        let isoform = if included[exon] { 1 } else { -1 };
        if (isoform == 1 && state == -1) || (isoform == -1 && state == 1) {
            return true;
        }
    }
    false
}

fn resolve_compatible(
    gene: &Gene,
    compatible: &[bool],
    pct: &[f64],
    low: f64,
    high: f64,
) -> Option<usize> {
    let hits: Vec<usize> = compatible
        .iter()
        .enumerate()
        .filter_map(|(index, on)| on.then_some(index))
        .collect();
    if hits.is_empty() {
        return None;
    }
    if hits.len() == 1 {
        return Some(hits[0]);
    }
    let transformed: Vec<f64> = pct
        .iter()
        .map(|&value| {
            if value > high {
                1.0
            } else if value < low {
                0.0
            } else if high > low {
                (value - low) / (high - low)
            } else {
                0.0
            }
        })
        .collect();
    let mut best = hits[0];
    let mut best_dist = f64::MAX;
    for index in hits {
        let mut distance = 0.0;
        for (exon, &value) in transformed.iter().enumerate() {
            let isoform = if gene.isoforms[index].exons.contains(&exon) {
                1.0
            } else {
                0.0
            };
            distance += (isoform - value).abs();
        }
        if distance < best_dist {
            best_dist = distance;
            best = index;
        }
    }
    Some(best)
}

fn isoform_score(gene: &Gene, index: usize, bases: u32) -> f64 {
    let length: i64 = gene.isoforms[index]
        .exons
        .iter()
        .map(|&exon| gene.exons[exon].1 - gene.exons[exon].0)
        .sum();
    if length <= 0 {
        -1.0
    } else {
        bases as f64 / length as f64
    }
}

fn adjust_truncation(
    pct: &[f64],
    lengths: &[i64],
    strand: char,
    poly: bool,
    truncation: f64,
) -> Vec<f64> {
    if poly {
        let mut out = pct.to_vec();
        if strand == '-' {
            for index in (0..out.len()).rev() {
                if out[index] > 0.0 {
                    if out[index] >= truncation || lengths[index] as f64 * out[index] >= 100.0 {
                        out[index] = 1.0;
                    }
                    break;
                }
            }
        } else {
            for index in 0..out.len() {
                if out[index] > 0.0 {
                    if out[index] >= truncation || lengths[index] as f64 * out[index] >= 100.0 {
                        out[index] = 1.0;
                    }
                    break;
                }
            }
        }
        return out;
    }
    let mut left = pct.to_vec();
    let mut left_bases = 0.0;
    for index in 0..left.len() {
        if left[index] > 0.0 {
            if left[index] >= truncation || lengths[index] as f64 * left[index] > 100.0 {
                left_bases = left[index] * lengths[index] as f64;
                left[index] = 1.0;
            }
            break;
        }
    }
    let mut right = pct.to_vec();
    let mut right_bases = 0.0;
    for index in (0..right.len()).rev() {
        if right[index] > 0.0 {
            if right[index] >= truncation || lengths[index] as f64 * right[index] > 100.0 {
                right_bases = right[index] * lengths[index] as f64;
                right[index] = 1.0;
            }
            break;
        }
    }
    if left_bases >= right_bases {
        left
    } else {
        right
    }
}

fn exon_vector(pct: &[f64], high: f64, low: f64, poly: bool, strand: char) -> Vec<i8> {
    let mut raw = Vec::with_capacity(pct.len());
    for &value in pct {
        let mapped = if value >= high { 1 } else { 0 };
        let skipped = if value <= low { -1 } else { 0 };
        raw.push(mapped + skipped);
    }
    if poly {
        if strand == '-' {
            let last = raw.iter().rposition(|&value| value == 1);
            return raw
                .into_iter()
                .enumerate()
                .map(|(index, value)| match last {
                    Some(last) if index <= last => value,
                    Some(_) => 0,
                    None => 0,
                })
                .collect();
        }
        let first = raw.iter().position(|&value| value == 1);
        return raw
            .into_iter()
            .enumerate()
            .map(|(index, value)| match first {
                Some(first) if index >= first => value,
                Some(_) => 0,
                None => 0,
            })
            .collect();
    }
    let first = raw.iter().position(|&value| value == 1);
    let last = raw.iter().rposition(|&value| value == 1);
    match (first, last) {
        (Some(first), Some(last)) => raw
            .into_iter()
            .enumerate()
            .map(|(index, value)| {
                if (first..=last).contains(&index) {
                    value
                } else {
                    0
                }
            })
            .collect(),
        _ => vec![0; raw.len()],
    }
}

pub fn exon_base_counts(positions: &[i64], exons: &[(i64, i64)]) -> Vec<u32> {
    let mut counts = vec![0u32; exons.len()];
    let mut exon = 0usize;
    for &pos in positions {
        while exon < exons.len() && pos >= exons[exon].1 {
            exon += 1;
        }
        if exon < exons.len() && pos >= exons[exon].0 {
            counts[exon] += 1;
        }
    }
    counts
}

fn round2(value: f64) -> f64 {
    let scaled = value * 100.0;
    let floor = scaled.floor();
    let fraction = scaled - floor;
    let up = if (fraction - 0.5).abs() < 1e-8 {
        (floor as i64).rem_euclid(2) == 1
    } else {
        fraction > 0.5
    };
    (if up { floor + 1.0 } else { floor }) / 100.0
}

fn mean(values: &[i64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<i64>() as f64 / values.len() as f64
    }
}

/// Longest run of `base` in `seq[start, end)`.
pub fn homopolymer(seq: &str, start: i64, end: i64, base: u8) -> i64 {
    let bytes = seq.as_bytes();
    let from = start.max(0) as usize;
    let to = (end.max(0) as usize).min(bytes.len());
    if from >= to {
        return 0;
    }
    let mut best = 0i64;
    let mut run = 0i64;
    for &byte in &bytes[from..to] {
        if byte == base || byte == base.to_ascii_lowercase() {
            run += 1;
            best = best.max(run);
        } else {
            run = 0;
        }
    }
    best
}

pub fn query_poly(seq: &str, query_start: i64, query_end: i64, umi: Option<&str>) -> bool {
    if seq.is_empty() {
        return false;
    }
    if let Some(umi) = umi.filter(|umi| !umi.is_empty() && *umi != "NA") {
        let rc = revcomp(umi);
        if let Some(pos) = seq.find(umi) {
            let from = pos + umi.len();
            let to = (from + 15).min(seq.len());
            if seq.as_bytes()[from..to]
                .iter()
                .filter(|&&b| b == b'T' || b == b't')
                .count()
                >= 10
            {
                return true;
            }
        }
        if let Some(pos) = seq.find(&rc) {
            let from = pos.saturating_sub(15);
            if seq.as_bytes()[from..pos]
                .iter()
                .filter(|&&b| b == b'A' || b == b'a')
                .count()
                >= 10
            {
                return true;
            }
        } else if seq.find(umi).is_none() {
            return detect_edge_poly(seq, query_start, query_end);
        }
        return false;
    }
    detect_edge_poly(seq, query_start, query_end)
}

pub fn detect_edge_poly(seq: &str, query_start: i64, query_end: i64) -> bool {
    let window = 15i64;
    let n = seq.len() as i64;
    let head_from = (query_start - 2 * window).max(0);
    let head_to = (query_start + window).min(query_end).min(n - 1);
    let tail_from = query_start.max(query_end - window);
    let tail_to = (query_end + 2 * window).min(n - 1);
    window_poly(seq, head_from, head_to + 1, b'T') || window_poly(seq, tail_from, tail_to + 1, b'A')
}

fn window_poly(seq: &str, start: i64, end: i64, base: u8) -> bool {
    let bytes = seq.as_bytes();
    let from = start.max(0) as usize;
    let to = (end.max(0) as usize).min(bytes.len());
    if to < from + 15 {
        return false;
    }
    bytes[from..to].windows(15).any(|window| {
        window
            .iter()
            .filter(|&&byte| byte == base || byte == base.to_ascii_lowercase())
            .count()
            >= 10
    })
}

fn revcomp(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'a' => 't',
            't' => 'a',
            'g' => 'c',
            'c' => 'g',
            _ => 'N',
        })
        .collect()
}

#[cfg(test)]
#[derive(Clone, Debug)]
pub(crate) struct Diagnosis {
    pub kind: &'static str,
    pub isoform: String,
    pub scores: Vec<f64>,
    pub exon_vector: Vec<i8>,
    pub pct0: Vec<f64>,
    pub compatible: Vec<i8>,
}

/// Same decision as read assignment, plus the exon vector and per-isoform scores
/// SCOTCH returns from `map_read_to_gene`.
#[cfg(test)]
pub(crate) fn diagnose(gene: &Gene, counts: &[u32], params: &MatchParams, poly: bool) -> Diagnosis {
    let bases: u32 = counts.iter().sum();
    let n_iso = gene.isoforms.len();
    if bases == 0 {
        return Diagnosis {
            kind: "uncategorized",
            isoform: "-".to_string(),
            scores: vec![-1.0; n_iso],
            exon_vector: vec![0; gene.exons.len()],
            pct0: vec![0.0; gene.exons.len()],
            compatible: vec![0; n_iso],
        };
    }
    let lengths: Vec<i64> = gene.exons.iter().map(|(a, b)| b - a).collect();
    let pct0: Vec<f64> = counts
        .iter()
        .zip(&lengths)
        .map(|(&count, &len)| {
            if count == 0 || len <= 0 {
                0.0
            } else {
                round2(count as f64 / len as f64)
            }
        })
        .collect();
    let mut pct = adjust_truncation(&pct0, &lengths, gene.strand, poly, params.truncation_match);
    let mut vector = exon_vector(&pct, params.high, params.low, poly, gene.strand);
    if vector.iter().all(|value| *value == 0) {
        if let Some((index, _)) = pct
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1).then(b.0.cmp(&a.0)))
        {
            pct[index] = 1.0;
            vector = exon_vector(&pct, params.high, params.low, poly, gene.strand);
        }
    }
    let exon_vector = vector.clone();
    let mut compatible = known_compatible(gene, &vector, &lengths, params.small_exon);
    if compatible.iter().filter(|&&on| on).count() == 0 {
        let mean = mean(&lengths);
        let cap = params
            .small_exon
            .max(params.small_exon_high.min(mean.floor() as i64));
        let mut threshold = params.small_exon + 10;
        while threshold <= cap {
            compatible = known_compatible(gene, &vector, &lengths, threshold);
            if compatible.iter().any(|&on| on) {
                break;
            }
            threshold += 10;
        }
    }
    let chosen = resolve_compatible(gene, &compatible, &pct, params.low, params.high);
    let mut flags = vec![0i8; n_iso];
    let mut scores = vec![-1.0; n_iso];
    if let Some(index) = chosen {
        flags[index] = 1;
        scores[index] = isoform_score(gene, index, bases);
        Diagnosis {
            kind: "known",
            isoform: gene.isoforms[index].name.clone(),
            scores,
            exon_vector,
            pct0,
            compatible: flags,
        }
    } else {
        Diagnosis {
            kind: "novel",
            isoform: "-".to_string(),
            scores,
            exon_vector,
            pct0,
            compatible: flags,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gtf::parse_gtf;

    fn gene() -> Gene {
        parse_gtf(
            "\
chr1\t.\tgene\t101\t600\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\";
chr1\t.\ttranscript\t101\t600\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ALL\";
chr1\t.\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ALL\";
chr1\t.\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ALL\";
chr1\t.\texon\t501\t600\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ALL\";
chr1\t.\ttranscript\t101\t400\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"HEAD\";
chr1\t.\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"HEAD\";
chr1\t.\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"HEAD\";
",
        )
        .unwrap()
        .genes
        .into_iter()
        .next()
        .unwrap()
    }

    #[test]
    fn full_and_skipped_and_novel_patterns_pick_different_outcomes() {
        let gene = gene();
        let params = MatchParams::default();
        let full: Vec<i64> = (100..200).chain(300..400).chain(500..600).collect();
        let place = place_in_metagene(
            &full,
            100,
            600,
            std::slice::from_ref(&gene),
            &[0],
            &params,
            &[false],
        )
        .unwrap();
        assert!(matches!(place.outcome, Outcome::Known(0)));
        let head: Vec<i64> = (100..200).chain(300..400).collect();
        let place = place_in_metagene(
            &head,
            100,
            400,
            std::slice::from_ref(&gene),
            &[0],
            &params,
            &[false],
        )
        .unwrap();
        match place.outcome {
            Outcome::Known(index) => assert_eq!(gene.isoforms[index].name, "HEAD"),
            other => panic!("expected HEAD, got {other:?}"),
        }
        let novel: Vec<i64> = (100..200).chain(500..600).collect();
        let place = place_in_metagene(
            &novel,
            100,
            600,
            std::slice::from_ref(&gene),
            &[0],
            &params,
            &[false],
        )
        .unwrap();
        match place.outcome {
            Outcome::Novel(vector) => assert_eq!(vector, vec![1, -1, 1]),
            other => panic!("expected novel, got {other:?}"),
        }
    }
}
