//! Python 2 `str(round(x, 2))`, the text TAMA writes for coverage, identity,
//! and polyA percentages.
//!
//! Scaling by 100 in binary does not match CPython here: `0.015 * 100` lands
//! above 1.5 even though the double is below the decimal halfway point.
//! Non-halfway values use Rust's correctly rounded `{:.2}` (half-even agrees
//! with half-away when the value is not exactly halfway). Exact halfway
//! values are binary floats whose 2-valuation is -3 (`_Py_double_round` in
//! CPython 2.7). Those round half away from zero, which is what Python 2 does
//! and what `{:.2}` does not (`85.125` prints `85.13`, not `85.12`).

/// `str(round(x, 2))` for a finite percentage-scale `x`.
pub fn py2_str_round(x: f64) -> String {
    if x == 0.0 || !x.is_finite() {
        return "0.0".to_string();
    }
    if two_valuation(x) == -3 {
        let y = x.abs() * 100.0;
        let cents = (y.floor() as i64) + 1;
        let cents = if x < 0.0 { -cents } else { cents };
        return format_cents(cents);
    }
    strip_fixed2(&format!("{x:.2}"))
}

/// Power of two in `x` when written as `odd * 2^k`. Zero is not passed in.
fn two_valuation(x: f64) -> i32 {
    let bits = x.abs().to_bits();
    let raw_exp = ((bits >> 52) & 0x7ff) as i32;
    let frac = bits & ((1u64 << 52) - 1);
    if raw_exp == 0 {
        if frac == 0 {
            return 0;
        }
        return frac.trailing_zeros() as i32 - 1074;
    }
    let significand = (1u64 << 52) | frac;
    significand.trailing_zeros() as i32 + (raw_exp - 1023 - 52)
}

fn strip_fixed2(text: &str) -> String {
    let (neg, rest) = match text.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, text),
    };
    let (whole, frac) = rest.split_once('.').unwrap_or((rest, "0"));
    let frac = frac.trim_end_matches('0');
    let frac = if frac.is_empty() { "0" } else { frac };
    if neg {
        format!("-{whole}.{frac}")
    } else {
        format!("{whole}.{frac}")
    }
}

fn format_cents(cents: i64) -> String {
    let neg = cents < 0;
    let cents = cents.unsigned_abs();
    let whole = cents / 100;
    let frac = cents % 100;
    let body = if frac == 0 {
        format!("{whole}.0")
    } else if frac % 10 == 0 {
        format!("{whole}.{}", frac / 10)
    } else {
        format!("{whole}.{frac:02}")
    };
    if neg {
        format!("-{body}")
    } else {
        body
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn matches_python2_round_strings() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/parity/py2_round_str.tsv");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|err| panic!("read {}: {err}", path.display()));
        let mut rows = 0;
        let mut mismatches = Vec::new();
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            rows += 1;
            let mut cols = line.split('\t');
            let _kind = cols.next().unwrap();
            let expr = cols.next().unwrap();
            let value_repr = cols.next().unwrap();
            let _round_repr = cols.next().unwrap();
            let round_str = cols.next().unwrap();
            let value: f64 = value_repr
                .parse()
                .unwrap_or_else(|_| panic!("parse {value_repr}"));
            let got = py2_str_round(value);
            if got != round_str {
                mismatches.push(format!(
                    "{expr}: value {value_repr} got {got} want {round_str}"
                ));
            }
        }
        assert!(rows > 0);
        assert!(
            mismatches.is_empty(),
            "{} mismatches:\n{}",
            mismatches.len(),
            mismatches.join("\n")
        );
    }
}
