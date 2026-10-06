//! Differential transcript usage.
//!
//! Gene counts use a two-sided Wilcoxon rank-sum test with the normal
//! approximation, tie correction, and continuity correction. Transcript
//! usage uses a Dirichlet-multinomial likelihood-ratio test. Gene p-values
//! are Holm-adjusted. Transcript and DTU gene p-values are Benjamini-Hochberg
//! adjusted. The chi-squared degrees of freedom are the number of isoform
//! columns, including a collapsed `other` column, which matches `LRT_test`.

use std::collections::BTreeMap;
use std::path::Path;

use statrs::distribution::{ChiSquared, ContinuousCDF, Normal};
use statrs::function::gamma::{digamma, ln_gamma};

#[derive(Clone, Debug)]
pub struct DtuSettings {
    pub epsilon: f64,
    pub group_novel: bool,
    pub rare: f64,
    pub min_cells: usize,
    pub min_total: f64,
}

impl Default for DtuSettings {
    fn default() -> Self {
        Self {
            epsilon: 0.01,
            group_novel: false,
            rare: 0.05,
            min_cells: 20,
            min_total: 20.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct DtuFailure {
    pub message: String,
    pub usage: bool,
}

impl std::fmt::Display for DtuFailure {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.message)
    }
}

#[derive(Clone, Debug)]
pub struct TranscriptStat {
    pub alpha1: f64,
    pub alpha2: f64,
    pub tu1: f64,
    pub tu2: f64,
    pub tu_diff: f64,
    pub tu_var1: f64,
    pub tu_var2: f64,
    pub dispersion1: f64,
    pub dispersion2: f64,
    pub isoform_switch: bool,
    pub isoform_switch_es: f64,
    pub p_dtu_gene: f64,
    pub p_transcript: f64,
    pub p_transcript_adj: f64,
    pub p_dtu_gene_adj: f64,
}

#[derive(Clone, Debug)]
pub struct DtuRow {
    pub gene: String,
    pub isoform: String,
    pub pct1: Option<f64>,
    pub pct2: Option<f64>,
    pub log_fc: Option<f64>,
    pub p_gene: Option<f64>,
    pub pct_diff: Option<f64>,
    pub p_gene_adj: Option<f64>,
    pub transcript: Option<TranscriptStat>,
}

#[derive(Clone, Debug)]
pub struct DtuTables {
    pub rows: Vec<DtuRow>,
}

#[derive(Clone, Debug)]
struct Matrix {
    cells: Vec<String>,
    features: Vec<String>,
    values: Vec<Vec<f64>>,
}

pub fn differential_usage(
    genes_a: &Path,
    transcripts_a: &Path,
    map_a: &Path,
    genes_b: &Path,
    transcripts_b: &Path,
    map_b: &Path,
    settings: &DtuSettings,
) -> Result<DtuTables, DtuFailure> {
    validate(settings)?;
    let gene_a = load_matrix(genes_a)?;
    let gene_b = load_matrix(genes_b)?;
    let mut tx_a = load_matrix(transcripts_a)?;
    let mut tx_b = load_matrix(transcripts_b)?;
    let map_a = load_map(map_a)?;
    let map_b = load_map(map_b)?;
    if settings.group_novel {
        collapse_novel_columns(&mut tx_a);
        collapse_novel_columns(&mut tx_b);
    }
    let mut gene_stats = wilcoxon_genes(&gene_a, &gene_b, settings.epsilon);
    let holm = holm_adjust(&gene_stats.values().map(|stat| stat.p).collect::<Vec<_>>());
    for (stat, adj) in gene_stats.values_mut().zip(holm) {
        stat.p_adj = adj;
    }
    let mut tested: Vec<(String, Vec<String>, Vec<TranscriptStat>)> = Vec::new();
    let mut gene_names: Vec<String> = gene_stats.keys().cloned().collect();
    for gene in map_a.keys().chain(map_b.keys()) {
        if !gene_names.iter().any(|have| have == gene) {
            gene_names.push(gene.clone());
        }
    }
    gene_names.sort();
    let mut dtu_gene_p = Vec::new();
    let mut dtu_tx_p = Vec::new();
    for gene in &gene_names {
        let cols_a = map_a.get(gene).map(Vec::as_slice).unwrap_or(&[]);
        let cols_b = map_b.get(gene).map(Vec::as_slice).unwrap_or(&[]);
        if cols_a.is_empty() && cols_b.is_empty() {
            continue;
        }
        if let Some((names, stats)) = transcript_lrt(&tx_a, cols_a, &tx_b, cols_b, settings) {
            dtu_gene_p.push(stats[0].p_dtu_gene);
            for stat in &stats {
                dtu_tx_p.push(stat.p_transcript);
            }
            tested.push((gene.clone(), names, stats));
        }
    }
    let gene_adj = bh_adjust(&dtu_gene_p);
    let mut tx_adj = bh_adjust(&dtu_tx_p).into_iter();
    for (gene_index, (_, _, stats)) in tested.iter_mut().enumerate() {
        let gene_p = gene_adj[gene_index];
        for stat in stats.iter_mut() {
            stat.p_dtu_gene_adj = gene_p;
            stat.p_transcript_adj = tx_adj.next().unwrap_or(1.0);
        }
    }
    let tested_map: BTreeMap<_, _> = tested
        .into_iter()
        .map(|(gene, names, stats)| (gene, (names, stats)))
        .collect();
    let mut rows = Vec::new();
    for gene in gene_names {
        let gene_stat = gene_stats.get(&gene);
        if let Some((names, stats)) = tested_map.get(&gene) {
            for (name, stat) in names.iter().zip(stats) {
                rows.push(row_from(&gene, name, gene_stat, Some(stat.clone())));
            }
        } else if gene_stat.is_some() {
            rows.push(row_from(&gene, "", gene_stat, None));
        }
    }
    Ok(DtuTables { rows })
}

pub fn write_dtu(path: &Path, tables: &DtuTables) -> Result<(), DtuFailure> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|err| DtuFailure {
                message: format!("create {}: {err}", parent.display()),
                usage: false,
            })?;
        }
    }
    let mut text = String::from(
        "gene\tisoform\tpct1\tpct2\tlogFC\tp_gene\tpct_diff\tp_gene_adj\talpha1\talpha2\tTU1\tTU2\tTU_diff\tTU_var1\tTU_var2\tdispersion1\tdispersion2\tisoform_switch\tisoform_switch_ES\tp_DTU_gene\tp_transcript\tp_transcript_adj\tp_DTU_gene_adj\n",
    );
    for row in &tables.rows {
        text.push_str(&tsv(row));
        text.push('\n');
    }
    std::fs::write(path, text).map_err(|err| DtuFailure {
        message: format!("write {}: {err}", path.display()),
        usage: false,
    })
}

fn validate(settings: &DtuSettings) -> Result<(), DtuFailure> {
    if settings.epsilon < 0.0 || settings.rare < 0.0 || settings.min_total < 0.0 {
        return Err(DtuFailure {
            message: "DTU epsilon, rare-isoform threshold, and minimum total must be non-negative"
                .to_string(),
            usage: true,
        });
    }
    Ok(())
}

fn row_from(
    gene: &str,
    isoform: &str,
    gene_stat: Option<&GeneStat>,
    transcript: Option<TranscriptStat>,
) -> DtuRow {
    DtuRow {
        gene: gene.to_string(),
        isoform: isoform.to_string(),
        pct1: gene_stat.map(|stat| stat.pct1),
        pct2: gene_stat.map(|stat| stat.pct2),
        log_fc: gene_stat.map(|stat| stat.log_fc),
        p_gene: gene_stat.map(|stat| stat.p),
        pct_diff: gene_stat.map(|stat| stat.pct_diff),
        p_gene_adj: gene_stat.map(|stat| stat.p_adj),
        transcript,
    }
}

#[derive(Clone, Debug)]
struct GeneStat {
    pct1: f64,
    pct2: f64,
    log_fc: f64,
    p: f64,
    pct_diff: f64,
    p_adj: f64,
}

fn wilcoxon_genes(a: &Matrix, b: &Matrix, epsilon: f64) -> BTreeMap<String, GeneStat> {
    let mut out = BTreeMap::new();
    for feature in &a.features {
        let Some(col_b) = b.features.iter().position(|name| name == feature) else {
            continue;
        };
        let col_a = a.features.iter().position(|name| name == feature).unwrap();
        let left: Vec<f64> = a.values.iter().map(|row| row[col_a] + epsilon).collect();
        let right: Vec<f64> = b.values.iter().map(|row| row[col_b] + epsilon).collect();
        let pct1 = fraction_above(&left, epsilon);
        let pct2 = fraction_above(&right, epsilon);
        let mean1 = mean(&left);
        let mean2 = mean(&right);
        let log_fc = if mean1 > 0.0 && mean2 > 0.0 {
            (mean1 / mean2).log2()
        } else {
            f64::NAN
        };
        out.insert(
            feature.clone(),
            GeneStat {
                pct1,
                pct2,
                log_fc,
                p: wilcoxon(&left, &right),
                pct_diff: pct1 - pct2,
                p_adj: 1.0,
            },
        );
    }
    out
}

fn fraction_above(values: &[f64], epsilon: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().filter(|value| **value > epsilon).count() as f64 / values.len() as f64
}

fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        0.0
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

fn wilcoxon(left: &[f64], right: &[f64]) -> f64 {
    let n1 = left.len();
    let n2 = right.len();
    if n1 == 0 || n2 == 0 {
        return 1.0;
    }
    let mut values: Vec<(f64, u8)> = left
        .iter()
        .map(|&value| (value, 0))
        .chain(right.iter().map(|&value| (value, 1)))
        .collect();
    values.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut ranks = vec![0.0; values.len()];
    let mut tie = 0.0;
    let mut index = 0;
    while index < values.len() {
        let mut next = index + 1;
        while next < values.len() && values[next].0 == values[index].0 {
            next += 1;
        }
        let rank = (index + 1 + next) as f64 / 2.0;
        for slot in ranks.iter_mut().take(next).skip(index) {
            *slot = rank;
        }
        let size = (next - index) as f64;
        if size > 1.0 {
            tie += size.powi(3) - size;
        }
        index = next;
    }
    let w: f64 = ranks
        .iter()
        .zip(&values)
        .filter(|(_, value)| value.1 == 0)
        .map(|(rank, _)| *rank)
        .sum();
    let n = (n1 + n2) as f64;
    let n1f = n1 as f64;
    let n2f = n2 as f64;
    let expected = n1f * (n + 1.0) / 2.0;
    let mut variance = n1f * n2f * (n + 1.0) / 12.0;
    if n > 1.0 {
        variance -= n1f * n2f * tie / (12.0 * n * (n - 1.0));
    }
    if variance <= 0.0 {
        return 1.0;
    }
    let diff = w - expected;
    let z = if diff == 0.0 {
        0.0
    } else {
        (diff - diff.signum() * 0.5) / variance.sqrt()
    };
    let normal = Normal::new(0.0, 1.0).expect("standard normal");
    (2.0 * (1.0 - normal.cdf(z.abs()))).clamp(0.0, 1.0)
}

fn transcript_lrt(
    a: &Matrix,
    cols_a: &[String],
    b: &Matrix,
    cols_b: &[String],
    settings: &DtuSettings,
) -> Option<(Vec<String>, Vec<TranscriptStat>)> {
    if a.cells.len() < settings.min_cells || b.cells.len() < settings.min_cells {
        return None;
    }
    let (names, left, right) = align_isoforms(a, cols_a, b, cols_b, settings.rare);
    if names.is_empty() {
        return None;
    }
    let sum_left: f64 = left.iter().flatten().sum();
    let sum_right: f64 = right.iter().flatten().sum();
    if sum_left < settings.min_total
        || sum_right < settings.min_total
        || sum_left == 0.0
        || sum_right == 0.0
    {
        return None;
    }
    let left_pos = positive(&left);
    let right_pos = positive(&right);
    if left_pos.is_empty() || right_pos.is_empty() {
        return None;
    }
    let mut pooled = left_pos.clone();
    pooled.extend(right_pos.iter().cloned());
    let h0 = mle(&pooled)?;
    let h1 = mle(&left_pos)?;
    let h2 = mle(&right_pos)?;
    let k = names.len();
    let p_gene = chi_p(2.0 * (h0.0 - h1.0 - h2.0), k);
    let mut p_tx = Vec::new();
    if k > 2 {
        for column in 0..k {
            let left_c = collapse_column(&left_pos, column);
            let right_c = collapse_column(&right_pos, column);
            let mut pooled_c = left_c.clone();
            pooled_c.extend(right_c.iter().cloned());
            let t0 = mle(&pooled_c)?;
            let t1 = mle(&left_c)?;
            let t2 = mle(&right_c)?;
            p_tx.push(chi_p(2.0 * (t0.0 - t1.0 - t2.0), k));
        }
    } else {
        p_tx = vec![p_gene; k];
    }
    let alpha0_1: f64 = h1.1.iter().sum();
    let alpha0_2: f64 = h2.1.iter().sum();
    let tu1: Vec<f64> = h1.1.iter().map(|alpha| alpha / alpha0_1).collect();
    let tu2: Vec<f64> = h2.1.iter().map(|alpha| alpha / alpha0_2).collect();
    let (switch, effect) = isoform_switch(&names, &tu1, &tu2);
    let stats = names
        .iter()
        .enumerate()
        .map(|(index, _)| TranscriptStat {
            alpha1: h1.1[index],
            alpha2: h2.1[index],
            tu1: tu1[index],
            tu2: tu2[index],
            tu_diff: (tu1[index] - tu2[index]).abs(),
            tu_var1: tu1[index] * (1.0 - tu1[index]) / (1.0 + alpha0_1),
            tu_var2: tu2[index] * (1.0 - tu2[index]) / (1.0 + alpha0_2),
            dispersion1: 1.0 / (alpha0_1 + 1.0),
            dispersion2: 1.0 / (alpha0_2 + 1.0),
            isoform_switch: switch,
            isoform_switch_es: effect,
            p_dtu_gene: p_gene,
            p_transcript: p_tx[index],
            p_transcript_adj: 1.0,
            p_dtu_gene_adj: 1.0,
        })
        .collect();
    Some((names, stats))
}

fn isoform_switch(names: &[String], tu1: &[f64], tu2: &[f64]) -> (bool, f64) {
    let keep: Vec<usize> = names
        .iter()
        .enumerate()
        .filter_map(|(index, name)| (name != "other").then_some(index))
        .collect();
    if keep.is_empty() {
        return (false, 0.0);
    }
    let max1 = keep
        .iter()
        .map(|&index| tu1[index])
        .fold(f64::MIN, f64::max);
    let max2 = keep
        .iter()
        .map(|&index| tu2[index])
        .fold(f64::MIN, f64::max);
    let tops: Vec<usize> = keep
        .into_iter()
        .filter(|&index| tu1[index] == max1 || tu2[index] == max2)
        .collect();
    let switch = tops.len() != 1;
    let effect = if tops.len() >= 2 && switch {
        ((tu1[tops[0]] - tu2[tops[0]]) + (tu2[tops[1]] - tu1[tops[1]])).abs()
    } else {
        (tu1[tops[0]] - tu2[tops[0]]).abs()
    };
    (switch, effect)
}

fn align_isoforms(
    a: &Matrix,
    cols_a: &[String],
    b: &Matrix,
    cols_b: &[String],
    rare: f64,
) -> (Vec<String>, Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let left = select_columns(a, cols_a);
    let right = select_columns(b, cols_b);
    let mut order: Vec<String> = cols_a
        .iter()
        .filter(|name| usage(&left, cols_a, name) >= rare || usage(&right, cols_b, name) >= rare)
        .cloned()
        .collect();
    for name in cols_b {
        if (usage(&left, cols_a, name) >= rare || usage(&right, cols_b, name) >= rare)
            && !order.iter().any(|have| have == name)
        {
            order.push(name.clone());
        }
    }
    if order.is_empty() {
        order.extend(cols_a.iter().cloned());
        for name in cols_b {
            if !order.iter().any(|have| have == name) {
                order.push(name.clone());
            }
        }
        return pack(&order, &left, cols_a, &right, cols_b, false);
    }
    pack(&order, &left, cols_a, &right, cols_b, true)
}

fn pack(
    order: &[String],
    left: &[Vec<f64>],
    cols_a: &[String],
    right: &[Vec<f64>],
    cols_b: &[String],
    collapse_rest: bool,
) -> (Vec<String>, Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mut names = order.to_vec();
    let mut out_left = columns_or_zero(left, cols_a, order);
    let mut out_right = columns_or_zero(right, cols_b, order);
    if collapse_rest {
        let rest_left = rest_sum(left, cols_a, order);
        let rest_right = rest_sum(right, cols_b, order);
        if rest_left.iter().any(|&value| value != 0.0)
            || rest_right.iter().any(|&value| value != 0.0)
        {
            names.push("other".to_string());
            push_column(&mut out_left, &rest_left);
            push_column(&mut out_right, &rest_right);
        }
    }
    if out_left.first().map(Vec::len) != out_right.first().map(Vec::len) {
        let width_left = out_left.first().map_or(0, Vec::len);
        let width_right = out_right.first().map_or(0, Vec::len);
        if width_left < width_right {
            names.push("other".to_string());
            push_column(&mut out_left, &vec![0.0; left.len()]);
        } else if width_right < width_left {
            names.push("other".to_string());
            push_column(&mut out_right, &vec![0.0; right.len()]);
        }
    }
    (names, out_left, out_right)
}

fn select_columns(matrix: &Matrix, columns: &[String]) -> Vec<Vec<f64>> {
    if matrix.cells.is_empty() {
        return Vec::new();
    }
    matrix
        .values
        .iter()
        .map(|row| {
            columns
                .iter()
                .map(|name| {
                    matrix
                        .features
                        .iter()
                        .position(|feature| feature == name)
                        .map_or(0.0, |index| row[index])
                })
                .collect()
        })
        .collect()
}

fn columns_or_zero(rows: &[Vec<f64>], have: &[String], want: &[String]) -> Vec<Vec<f64>> {
    rows.iter()
        .map(|row| {
            want.iter()
                .map(|name| {
                    have.iter()
                        .position(|feature| feature == name)
                        .map_or(0.0, |index| row[index])
                })
                .collect()
        })
        .collect()
}

fn rest_sum(rows: &[Vec<f64>], have: &[String], keep: &[String]) -> Vec<f64> {
    rows.iter()
        .map(|row| {
            have.iter()
                .enumerate()
                .filter(|(_, name)| !keep.iter().any(|kept| kept == *name))
                .map(|(index, _)| row[index])
                .sum()
        })
        .collect()
}

fn push_column(rows: &mut Vec<Vec<f64>>, column: &[f64]) {
    if rows.is_empty() {
        *rows = column.iter().map(|&value| vec![value]).collect();
        return;
    }
    for (row, value) in rows.iter_mut().zip(column) {
        row.push(*value);
    }
}

fn usage(rows: &[Vec<f64>], columns: &[String], name: &str) -> f64 {
    let Some(index) = columns.iter().position(|column| column == name) else {
        return 0.0;
    };
    let total: f64 = rows.iter().flatten().sum();
    if total == 0.0 {
        return 0.0;
    }
    rows.iter().map(|row| row[index]).sum::<f64>() / total
}

fn positive(rows: &[Vec<f64>]) -> Vec<Vec<f64>> {
    rows.iter()
        .filter(|row| row.iter().sum::<f64>() > 0.0)
        .cloned()
        .collect()
}

fn collapse_column(rows: &[Vec<f64>], column: usize) -> Vec<Vec<f64>> {
    rows.iter()
        .map(|row| {
            let rest: f64 = row
                .iter()
                .enumerate()
                .filter(|(index, _)| *index != column)
                .map(|(_, value)| *value)
                .sum();
            vec![row[column], rest]
        })
        .collect()
}

fn chi_p(stat: f64, df: usize) -> f64 {
    let stat = if stat.is_finite() { stat.max(0.0) } else { 0.0 };
    if df == 0 {
        return 1.0;
    }
    let chi = ChiSquared::new(df as f64).expect("positive degrees of freedom");
    (1.0 - chi.cdf(stat)).clamp(0.0, 1.0)
}

fn mle(rows: &[Vec<f64>]) -> Option<(f64, Vec<f64>)> {
    if rows.is_empty() || rows[0].is_empty() {
        return None;
    }
    let width = rows[0].len();
    let mut starts = vec![vec![1.0; width]];
    let total: f64 = rows.iter().flatten().sum();
    if total > 0.0 {
        starts.push(
            (0..width)
                .map(|column| {
                    let sum: f64 = rows.iter().map(|row| row[column]).sum();
                    (sum / total).max(1e-4)
                })
                .collect(),
        );
    }
    let mut best: Option<(f64, Vec<f64>)> = None;
    for start in starts {
        if let Some(fitted) = newton(&start, rows) {
            if best.as_ref().is_none_or(|(nll, _)| fitted.0 < *nll) {
                best = Some(fitted);
            }
        }
    }
    best
}

fn newton(start: &[f64], rows: &[Vec<f64>]) -> Option<(f64, Vec<f64>)> {
    let mut alpha = project(start);
    let (mut nll, _, _) = dirichlet_nll(&alpha, rows)?;
    for _ in 0..40 {
        let (current, grad, hess) = dirichlet_nll(&alpha, rows)?;
        nll = current;
        let grad_norm: f64 = grad.iter().map(|value| value * value).sum::<f64>().sqrt();
        if grad_norm < 1e-5 {
            break;
        }
        let mut step = solve(hess, grad.iter().map(|value| -value).collect())
            .unwrap_or_else(|| grad.iter().map(|value| -value).collect());
        let slope: f64 = grad.iter().zip(&step).map(|(g, s)| g * s).sum();
        if slope >= 0.0 {
            step = grad.iter().map(|value| -value).collect();
        }
        let slope: f64 = grad.iter().zip(&step).map(|(g, s)| g * s).sum();
        let mut scale = 1.0;
        let mut improved = false;
        while scale >= 1e-8 {
            let trial = project(
                &alpha
                    .iter()
                    .zip(&step)
                    .map(|(value, delta)| value + scale * delta)
                    .collect::<Vec<_>>(),
            );
            if let Some((trial_nll, _, _)) = dirichlet_nll(&trial, rows) {
                if trial_nll <= nll + 1e-4 * scale * slope {
                    alpha = trial;
                    nll = trial_nll;
                    improved = true;
                    break;
                }
            }
            scale *= 0.5;
        }
        if !improved {
            break;
        }
    }
    if nll.is_finite() {
        Some((nll, alpha))
    } else {
        None
    }
}

/// Trigamma ψ₁(x) for x > 0. statrs 0.18 exposes digamma but not trigamma.
fn trigamma(mut x: f64) -> f64 {
    if x <= 0.0 || !x.is_finite() {
        return f64::NAN;
    }
    let mut extra = 0.0;
    while x < 8.0 {
        extra += 1.0 / (x * x);
        x += 1.0;
    }
    let inv = 1.0 / x;
    let inv2 = inv * inv;
    extra
        + inv
            * (1.0
                + inv
                    * (0.5
                        + inv
                            * (1.0 / 6.0
                                + inv2
                                    * (-1.0 / 30.0
                                        + inv2
                                            * (1.0 / 42.0
                                                + inv2 * (-1.0 / 30.0 + inv2 * (5.0 / 66.0)))))))
}

fn project(alpha: &[f64]) -> Vec<f64> {
    alpha
        .iter()
        .map(|value| {
            if value.is_finite() {
                value.max(1e-4)
            } else {
                1e-4
            }
        })
        .collect()
}

fn dirichlet_nll(alpha: &[f64], rows: &[Vec<f64>]) -> Option<(f64, Vec<f64>, Vec<Vec<f64>>)> {
    let width = alpha.len();
    let alpha0: f64 = alpha.iter().sum();
    if !alpha0.is_finite() || alpha.iter().any(|value| *value <= 0.0) {
        return None;
    }
    let mut nll = 0.0;
    let mut grad = vec![0.0; width];
    let mut hess = vec![vec![0.0; width]; width];
    for row in rows {
        let total: f64 = row.iter().sum();
        let mut log_lik = ln_gamma(alpha0) - ln_gamma(alpha0 + total);
        for column in 0..width {
            log_lik += ln_gamma(alpha[column] + row[column]) - ln_gamma(alpha[column]);
        }
        if !log_lik.is_finite() {
            return None;
        }
        nll -= log_lik;
        let d_alpha = digamma(alpha0) - digamma(alpha0 + total);
        let t_alpha = trigamma(alpha0) - trigamma(alpha0 + total);
        for column in 0..width {
            grad[column] -= d_alpha + digamma(alpha[column] + row[column]) - digamma(alpha[column]);
            for entry in &mut hess[column] {
                *entry -= t_alpha;
            }
            hess[column][column] -= trigamma(alpha[column] + row[column]) - trigamma(alpha[column]);
        }
    }
    if nll.is_finite() {
        Some((nll, grad, hess))
    } else {
        None
    }
}

fn solve(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Option<Vec<f64>> {
    let n = rhs.len();
    for col in 0..n {
        let mut pivot = col;
        for row in (col + 1)..n {
            if matrix[row][col].abs() > matrix[pivot][col].abs() {
                pivot = row;
            }
        }
        if matrix[pivot][col].abs() < 1e-10 {
            return None;
        }
        matrix.swap(col, pivot);
        rhs.swap(col, pivot);
        let div = matrix[col][col];
        for value in matrix[col].iter_mut().skip(col) {
            *value /= div;
        }
        rhs[col] /= div;
        for row in 0..n {
            if row == col {
                continue;
            }
            let factor = matrix[row][col];
            if factor == 0.0 {
                continue;
            }
            let pivot: Vec<f64> = matrix[col][col..].to_vec();
            for (value, pivot_value) in matrix[row][col..].iter_mut().zip(&pivot) {
                *value -= factor * pivot_value;
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    rhs.iter().all(|value| value.is_finite()).then_some(rhs)
}

fn holm_adjust(pvalues: &[f64]) -> Vec<f64> {
    adjust(pvalues, true)
}

fn bh_adjust(pvalues: &[f64]) -> Vec<f64> {
    adjust(pvalues, false)
}

fn adjust(pvalues: &[f64], holm: bool) -> Vec<f64> {
    let m = pvalues.len();
    if m == 0 {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&left, &right| pvalues[left].total_cmp(&pvalues[right]));
    let mut out = vec![1.0; m];
    if holm {
        let mut running = 0.0f64;
        for (rank, &index) in order.iter().enumerate() {
            let value = (pvalues[index] * (m - rank) as f64).min(1.0);
            running = running.max(value);
            out[index] = running;
        }
    } else {
        let mut running = 1.0f64;
        let total = m as f64;
        for rank in (0..m).rev() {
            let value = (pvalues[order[rank]] * total / (rank as f64 + 1.0)).min(1.0);
            running = running.min(value);
            out[order[rank]] = running;
        }
    }
    out
}

fn collapse_novel_columns(matrix: &mut Matrix) {
    let novel: Vec<usize> = matrix
        .features
        .iter()
        .enumerate()
        .filter_map(|(index, name)| name.contains("novel").then_some(index))
        .collect();
    if novel.is_empty() {
        return;
    }
    let mut name = strip_novel_suffix(&matrix.features[novel[0]]);
    if name.is_empty() {
        name = "novel".to_string();
    }
    let keep_existing = matrix
        .features
        .iter()
        .enumerate()
        .find(|(index, feature)| *feature == &name && !novel.contains(index))
        .map(|(index, _)| index);
    for row in &mut matrix.values {
        let sum: f64 = novel.iter().map(|&index| row[index]).sum();
        if let Some(index) = keep_existing {
            row[index] += sum;
        }
    }
    if keep_existing.is_none() {
        for row in &mut matrix.values {
            let sum: f64 = novel.iter().map(|&index| row[index]).sum();
            row.push(sum);
        }
        matrix.features.push(name);
    }
    let drop: std::collections::BTreeSet<usize> = novel.into_iter().collect();
    for row in &mut matrix.values {
        let mut kept = Vec::new();
        for (index, value) in row.iter().enumerate() {
            if !drop.contains(&index) {
                kept.push(*value);
            }
        }
        *row = kept;
    }
    matrix.features = matrix
        .features
        .iter()
        .enumerate()
        .filter(|(index, _)| !drop.contains(index))
        .map(|(_, name)| name.clone())
        .collect();
}

fn strip_novel_suffix(name: &str) -> String {
    let bytes = name.as_bytes();
    let mut index = 0;
    while index + 1 < bytes.len() {
        if bytes[index] == b'_' || bytes[index] == b'-' {
            let start = index + 1;
            if start < bytes.len() && bytes[start].is_ascii_digit() {
                let mut end = start;
                while end < bytes.len() && bytes[end].is_ascii_digit() {
                    end += 1;
                }
                let mut out = String::new();
                out.push_str(&name[..index]);
                out.push_str(&name[end..]);
                return out;
            }
        }
        index += 1;
    }
    name.to_string()
}

fn load_matrix(path: &Path) -> Result<Matrix, DtuFailure> {
    let text = std::fs::read_to_string(path).map_err(|err| DtuFailure {
        message: format!("open {}: {err}", path.display()),
        usage: false,
    })?;
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let Some(header) = lines.next() else {
        return Err(DtuFailure {
            message: format!("{} is empty", path.display()),
            usage: true,
        });
    };
    let header = parse_csv(header);
    if header.is_empty() {
        return Err(DtuFailure {
            message: format!("{} has no columns", path.display()),
            usage: true,
        });
    }
    let features = header[1..].to_vec();
    let mut cells = Vec::new();
    let mut values = Vec::new();
    for (lineno, line) in lines.enumerate() {
        let fields = parse_csv(line);
        if fields.len() != header.len() {
            return Err(DtuFailure {
                message: format!(
                    "{}:{} has {} columns, expected {}",
                    path.display(),
                    lineno + 2,
                    fields.len(),
                    header.len()
                ),
                usage: true,
            });
        }
        let mut row = Vec::with_capacity(features.len());
        for field in &fields[1..] {
            let value: f64 = field.parse().map_err(|_| DtuFailure {
                message: format!("{}:{} has a non-numeric count", path.display(), lineno + 2),
                usage: true,
            })?;
            row.push(value);
        }
        cells.push(fields[0].clone());
        values.push(row);
    }
    Ok(Matrix {
        cells,
        features,
        values,
    })
}

fn load_map(path: &Path) -> Result<BTreeMap<String, Vec<String>>, DtuFailure> {
    let text = std::fs::read_to_string(path).map_err(|err| DtuFailure {
        message: format!("open {}: {err}", path.display()),
        usage: false,
    })?;
    let mut lines = text.lines().filter(|line| !line.trim().is_empty());
    let Some(header) = lines.next() else {
        return Err(DtuFailure {
            message: format!("{} is empty", path.display()),
            usage: true,
        });
    };
    let header: Vec<&str> = header.split('\t').collect();
    let gene_col = header.iter().position(|name| *name == "genes");
    let tx_col = header.iter().position(|name| *name == "transcripts");
    let (Some(gene_col), Some(tx_col)) = (gene_col, tx_col) else {
        return Err(DtuFailure {
            message: format!("{} needs genes and transcripts columns", path.display()),
            usage: true,
        });
    };
    let mut map: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (lineno, line) in lines.enumerate() {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() <= gene_col.max(tx_col) {
            return Err(DtuFailure {
                message: format!("{}:{} is missing a column", path.display(), lineno + 2),
                usage: true,
            });
        }
        map.entry(fields[gene_col].to_string())
            .or_default()
            .push(fields[tx_col].to_string());
    }
    Ok(map)
}

fn parse_csv(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut quoted = false;
    let mut chars = line.chars().peekable();
    while let Some(ch) = chars.next() {
        if quoted {
            if ch == '"' {
                if chars.peek() == Some(&'"') {
                    field.push('"');
                    chars.next();
                } else {
                    quoted = false;
                }
            } else {
                field.push(ch);
            }
        } else if ch == '"' && field.is_empty() {
            quoted = true;
        } else if ch == ',' {
            fields.push(std::mem::take(&mut field));
        } else {
            field.push(ch);
        }
    }
    fields.push(field);
    fields
}

fn tsv(row: &DtuRow) -> String {
    let blank = String::new();
    let tx = row.transcript.as_ref();
    [
        row.gene.clone(),
        row.isoform.clone(),
        opt_num(row.pct1),
        opt_num(row.pct2),
        opt_num(row.log_fc),
        opt_p(row.p_gene),
        opt_num(row.pct_diff),
        opt_p(row.p_gene_adj),
        tx.map(|stat| fmt_num(stat.alpha1))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_num(stat.alpha2))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_num(stat.tu1))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_num(stat.tu2))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_num(stat.tu_diff))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_num(stat.tu_var1))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_num(stat.tu_var2))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_num(stat.dispersion1))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_num(stat.dispersion2))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| stat.isoform_switch.to_string())
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_num(stat.isoform_switch_es))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_p(stat.p_dtu_gene))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_p(stat.p_transcript))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_p(stat.p_transcript_adj))
            .unwrap_or_else(|| blank.clone()),
        tx.map(|stat| fmt_p(stat.p_dtu_gene_adj))
            .unwrap_or_else(|| blank.clone()),
    ]
    .join("\t")
}

fn opt_num(value: Option<f64>) -> String {
    value.map(fmt_num).unwrap_or_default()
}

fn opt_p(value: Option<f64>) -> String {
    value.map(fmt_p).unwrap_or_default()
}

fn fmt_num(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6}")
    } else {
        String::new()
    }
}

fn fmt_p(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        String::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn write_pair(dir: &std::path::Path, name: &str, gene: &str, tx: &str, map: &str) {
        std::fs::create_dir_all(dir).unwrap();
        std::fs::write(dir.join("gene_counts.csv"), gene).unwrap();
        std::fs::write(dir.join("transcript_counts.csv"), tx).unwrap();
        std::fs::write(dir.join("gene_transcript.tsv"), map).unwrap();
        let _ = name;
    }

    fn table(header: &str, count: usize, prefix: &str, values: &str) -> String {
        let mut text = format!("{header}\n");
        for index in 0..count {
            text.push_str(&format!("{prefix}{index},{values}\n"));
        }
        text
    }

    #[test]
    fn switched_usage_is_significant_and_matched_usage_is_not() {
        let root = std::env::temp_dir().join(format!("plena-scotch-dtu-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let strong_a = root.join("strong-a");
        let strong_b = root.join("strong-b");
        let map = "genes\ttranscripts\nG\tG_A\nG\tG_B\n";
        write_pair(
            &strong_a,
            "a",
            &table("cell,G", 22, "a", "20"),
            &table("cell,G_A,G_B", 22, "a", "20,0"),
            map,
        );
        write_pair(
            &strong_b,
            "b",
            &table("cell,G", 22, "b", "2"),
            &table("cell,G_A,G_B", 22, "b", "0,2"),
            map,
        );
        let tables = differential_usage(
            &strong_a.join("gene_counts.csv"),
            &strong_a.join("transcript_counts.csv"),
            &strong_a.join("gene_transcript.tsv"),
            &strong_b.join("gene_counts.csv"),
            &strong_b.join("transcript_counts.csv"),
            &strong_b.join("gene_transcript.tsv"),
            &DtuSettings::default(),
        )
        .unwrap();
        let row = tables
            .rows
            .iter()
            .find(|row| row.gene == "G" && row.transcript.is_some())
            .unwrap();
        let stat = row.transcript.as_ref().unwrap();
        assert!(stat.p_dtu_gene < 1e-4, "p {}", stat.p_dtu_gene);
        assert!(stat.isoform_switch);
        assert!(row.p_gene.unwrap() < 1e-6);

        let same_a = root.join("same-a");
        let same_b = root.join("same-b");
        write_pair(
            &same_a,
            "a",
            &table("cell,G", 22, "a", "24"),
            &table("cell,G_A,G_B", 22, "a", "18,6"),
            map,
        );
        write_pair(
            &same_b,
            "b",
            &table("cell,G", 22, "b", "24"),
            &table("cell,G_A,G_B", 22, "b", "18,6"),
            map,
        );
        let tables = differential_usage(
            &same_a.join("gene_counts.csv"),
            &same_a.join("transcript_counts.csv"),
            &same_a.join("gene_transcript.tsv"),
            &same_b.join("gene_counts.csv"),
            &same_b.join("transcript_counts.csv"),
            &same_b.join("gene_transcript.tsv"),
            &DtuSettings::default(),
        )
        .unwrap();
        let stat = tables.rows[0].transcript.as_ref().unwrap();
        assert!(stat.p_dtu_gene > 0.05, "p {}", stat.p_dtu_gene);
        assert!(!stat.isoform_switch);
        let _ = std::fs::remove_dir_all(&root);
    }
}
