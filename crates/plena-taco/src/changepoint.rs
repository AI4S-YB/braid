//! Change-point detection.
//!
//! The Mann-Whitney statistic and the mean-squared-error split follow
//! TACO's copies, not current SciPy. Expression arrays stay `f32`.
//! Smoothing and the slope use `f64`.

use crate::util::py_float;

#[derive(Clone, Debug)]
pub(crate) struct ChangePoint {
    pub pos: i64,
    pub start: i64,
    pub end: i64,
    pub pvalue: f64,
    pub foldchange: f64,
    pub sign: f64,
}

pub(crate) fn rankdata(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| values[i].total_cmp(&values[j]));
    let mut ranks = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && values[order[j]] == values[order[i]] {
            j += 1;
        }
        let average = ((i + 1) + j) as f64 / 2.0;
        for slot in order.iter().take(j).skip(i) {
            ranks[*slot] = average;
        }
        i = j;
    }
    ranks
}

pub(crate) fn tiecorrect(ranks: &[f64]) -> f64 {
    let n = ranks.len();
    if n < 2 {
        return 1.0;
    }
    let mut sorted = ranks.to_vec();
    sorted.sort_by(f64::total_cmp);
    let mut sum = 0.0;
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && sorted[j] == sorted[i] {
            j += 1;
        }
        let count = (j - i) as f64;
        sum += count * count * count - count;
        i = j;
    }
    let size = n as f64;
    1.0 - sum / (size * size * size - size)
}

/// Two-sided Mann-Whitney p-value with the continuity correction.
fn mannwhitneyu(x: &[f64], y: &[f64]) -> Result<f64, String> {
    let n1 = x.len();
    let n2 = y.len();
    let mut both = Vec::with_capacity(n1 + n2);
    both.extend_from_slice(x);
    both.extend_from_slice(y);
    let ranked = rankdata(&both);
    let rank_sum: f64 = ranked[..n1].iter().sum();
    let u1 = (n1 * n2) as f64 + (n1 * (n1 + 1)) as f64 / 2.0 - rank_sum;
    let u2 = (n1 * n2) as f64 - u1;
    let correction = tiecorrect(&ranked);
    if correction == 0.0 {
        return Err("All numbers are identical in amannwhitneyu".to_string());
    }
    let sd = (correction * (n1 * n2) as f64 * (n1 + n2 + 1) as f64 / 12.0).sqrt();
    let mean_rank = (n1 * n2) as f64 / 2.0 + 0.5;
    let z = (u1.max(u2) - mean_rank).abs() / sd;
    Ok(norm_sf(z) * 2.0)
}

fn norm_sf(z: f64) -> f64 {
    0.5 * erfc(z / std::f64::consts::SQRT_2)
}

/// Abramowitz and Stegun 7.1.26, maximum error about 1.5e-7.
fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let p = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    p * (-x * x).exp()
}

fn mwu(left: &[f64], right: &[f64]) -> f64 {
    if left.is_empty() || right.is_empty() {
        return 1.0;
    }
    if constant(left) && constant(right) && left[0] == right[0] {
        return 1.0;
    }
    mannwhitneyu(left, right).unwrap_or(1.0)
}

fn constant(values: &[f64]) -> bool {
    values.iter().all(|value| *value == values[0])
}

/// Index of the minimum within-cluster sum of squares. `-1` when the
/// series is too short to split. Arithmetic stays in `f32`.
fn mse_index(values: &[f32]) -> isize {
    let n = values.len();
    if n < 2 {
        return -1;
    }
    let mut right_sum: f32 = values.iter().sum();
    let mut left_sum = 0.0f32;
    let mut best_i = -1isize;
    let mut best = 0.0f32;
    let mut prev = 0usize;
    for i in 1..n {
        for value in values.iter().take(i).skip(prev) {
            left_sum += *value;
            right_sum -= *value;
        }
        let left_mean = left_sum / i as f32;
        let right_mean = right_sum / (n - i) as f32;
        let mut score = 0.0f32;
        for value in values.iter().take(i) {
            let delta = *value - left_mean;
            score += delta * delta;
        }
        for value in values.iter().skip(i) {
            let delta = *value - right_mean;
            score += delta * delta;
        }
        if best_i < 0 || score < best {
            best = score;
            best_i = i as isize;
        }
        prev = i;
    }
    best_i
}

fn mean_f32(values: &[f32]) -> f32 {
    if values.is_empty() {
        return 0.0;
    }
    let sum: f32 = values.iter().sum();
    sum / values.len() as f32
}

fn hanning(len: usize) -> Vec<f64> {
    if len <= 1 {
        return vec![1.0; len];
    }
    (0..len)
        .map(|n| 0.5 - 0.5 * (2.0 * std::f64::consts::PI * n as f64 / (len - 1) as f64).cos())
        .collect()
}

/// Reflect, convolve with a Hanning window, and return a series one
/// element longer than the input. That is TACO's slice of the valid
/// convolution.
fn smooth(values: &[f32], window_len: usize) -> Vec<f64> {
    let window = hanning(window_len);
    let weight: f64 = window.iter().sum();
    let mut signal = Vec::with_capacity(values.len() + 2 * (window_len - 1));
    // x[window_len-1:0:-1] and x[-1:-window_len:-1], each window_len-1 samples.
    for index in (1..window_len).rev() {
        signal.push(f64::from(values[index]));
    }
    signal.extend(values.iter().map(|value| f64::from(*value)));
    let tail = values.len() - (window_len - 1);
    for index in (tail..values.len()).rev() {
        signal.push(f64::from(values[index]));
    }
    let mut full = vec![0.0; signal.len() - window_len + 1];
    for (out_index, slot) in full.iter_mut().enumerate() {
        let origin = out_index + window_len - 1;
        let mut sum = 0.0;
        for (m, coeff) in window.iter().enumerate() {
            sum += (coeff / weight) * signal[origin - m];
        }
        *slot = sum;
    }
    let from = window_len / 2 - 1;
    let to = full.len() - window_len / 2;
    full[from..to].to_vec()
}

fn gradient(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![0.0];
    }
    let mut slope = vec![0.0; n];
    slope[0] = values[1] - values[0];
    slope[n - 1] = values[n - 1] - values[n - 2];
    for i in 1..n - 1 {
        slope[i] = (values[i + 1] - values[i - 1]) / 2.0;
    }
    slope
}

fn signum(value: f64) -> f64 {
    if value > 0.0 {
        1.0
    } else if value < 0.0 {
        -1.0
    } else {
        0.0
    }
}

/// `(steps left including the point, steps right including the point, sign)`.
/// A zero slope returns `(0, 0, 0)` and the caller rejects it.
fn slope_extract(slope: &[f64], index: usize) -> (usize, usize, f64) {
    let signs: Vec<f64> = slope.iter().copied().map(signum).collect();
    let sign = signs[index];
    if sign == 0.0 {
        return (0, 0, 0.0);
    }
    let mut left = 0;
    while index >= left && signs[index - left] == sign {
        left += 1;
        if index < left {
            break;
        }
    }
    let mut right = 0;
    while index + right < signs.len() && signs[index + right] == sign {
        right += 1;
    }
    (left, right, sign)
}

fn bin_seg(
    values: &[f32],
    slope: &[f64],
    pvalue: f64,
    fold_change: f64,
    size_cutoff: usize,
    points: &mut Vec<ChangePoint>,
    offset: usize,
) {
    if values.len() < size_cutoff {
        return;
    }
    let mut diff_indexes = Vec::new();
    for i in 0..values.len() - 1 {
        if values[i + 1] - values[i] != 0.0 {
            diff_indexes.push(i);
        }
    }
    let condensed: Vec<f32> = diff_indexes.iter().map(|&i| values[i]).collect();
    let diff_i = mse_index(&condensed);
    if diff_i < 0 {
        return;
    }
    let diff_i = diff_i as usize;
    let i = diff_indexes[diff_i];
    if i <= 1 {
        return;
    }
    let left: Vec<f64> = condensed[..diff_i].iter().map(|v| f64::from(*v)).collect();
    let right: Vec<f64> = condensed[diff_i..].iter().map(|v| f64::from(*v)).collect();
    let p = mwu(&left, &right);
    if p >= pvalue {
        return;
    }
    let mut mean_left = mean_f32(&values[..i]);
    let mut mean_right = mean_f32(&values[i..]);
    if mean_left > mean_right {
        std::mem::swap(&mut mean_left, &mut mean_right);
    }
    let ratio = f64::from(mean_left / mean_right);
    if ratio >= fold_change {
        return;
    }
    let (left_run, right_run, sign) = slope_extract(slope, offset + i);
    if left_run == 0 || right_run == 0 {
        return;
    }
    points.push(ChangePoint {
        pos: (offset + i) as i64,
        start: (offset + i - left_run) as i64,
        end: (offset + i + right_run) as i64,
        pvalue: p,
        foldchange: ratio,
        sign,
    });
    if i > left_run {
        bin_seg(
            &values[..i - left_run],
            slope,
            pvalue,
            fold_change,
            size_cutoff,
            points,
            offset,
        );
    }
    if offset + i + right_run < offset + values.len() {
        bin_seg(
            &values[i + right_run..],
            slope,
            pvalue,
            fold_change,
            size_cutoff,
            points,
            offset + i + right_run,
        );
    }
}

pub(crate) fn run_changepoint(
    values: &[f32],
    pvalue: f64,
    fold_change: f64,
    size_cutoff: usize,
    window_len: usize,
) -> Vec<ChangePoint> {
    if values.len() < window_len {
        return Vec::new();
    }
    let slope = gradient(&smooth(values, window_len));
    let mut points = Vec::new();
    bin_seg(
        values,
        &slope,
        pvalue,
        fold_change,
        size_cutoff,
        &mut points,
        0,
    );
    points
}

pub(crate) fn find_threshold_points(values: &[f32], start: i64) -> Vec<i64> {
    if values.len() < 2 {
        return Vec::new();
    }
    let mut points = Vec::new();
    for i in 1..values.len() {
        let fell = values[i - 1] > 0.0 && values[i] <= 0.0;
        let rose = values[i - 1] <= 0.0 && values[i] > 0.0;
        if fell || rose {
            points.push(start + i as i64);
        }
    }
    points
}

pub(crate) fn format_f32(value: f32) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value.is_sign_negative() {
            "-inf".to_string()
        } else {
            "inf".to_string()
        };
    }
    for places in 1..12 {
        let text = format!("{value:.places$}");
        if text.parse::<f32>().ok() == Some(value) {
            return if text.contains(['.', 'e', 'E']) {
                text
            } else {
                format!("{text}.0")
            };
        }
    }
    py_float(f64::from(value))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranks_and_tie_correction_match_scipy_examples() {
        let ranks = rankdata(&[0.0, 2.0, 3.0, 2.0]);
        assert_eq!(ranks, vec![1.0, 2.5, 4.0, 2.5]);
        assert!((tiecorrect(&[1.0, 2.5, 2.5, 4.0]) - 0.9).abs() < 1e-12);
    }

    #[test]
    fn ramp_change_point_matches_the_taco_oracle() {
        // 50 full transcripts over [1000, 1220) plus 20 partials ending at
        // 1100..1119. Positions are relative to 1000.
        let mut expr = vec![50.0f32; 220];
        for end in 100..120 {
            for slot in expr.iter_mut().take(end) {
                *slot += 10.0;
            }
        }
        assert_eq!(expr[0], 250.0);
        assert_eq!(expr[99], 250.0);
        assert_eq!(expr[100], 240.0);
        assert_eq!(expr[110], 140.0);
        assert_eq!(*expr.last().unwrap(), 50.0);
        let points = run_changepoint(&expr, 0.05, 0.80, 20, 11);
        assert_eq!(points.len(), 1, "{points:?}");
        let point = &points[0];
        assert_eq!(point.pos, 109);
        assert_eq!(point.start, 95);
        assert_eq!(point.end, 125);
        assert_eq!(point.sign, -1.0);
        assert!(
            (point.pvalue - 0.0001826717911095504).abs() < 1e-6,
            "{}",
            point.pvalue
        );
        assert!(
            (point.foldchange - 0.22351081669330597).abs() < 1e-5,
            "{}",
            point.foldchange
        );
    }
}
