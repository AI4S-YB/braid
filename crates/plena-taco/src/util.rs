pub(crate) fn bisect_left(a: &[i64], x: i64, mut lo: usize) -> usize {
    let mut hi = a.len();
    if lo > hi {
        lo = hi;
    }
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] < x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

pub(crate) fn bisect_right(a: &[i64], x: i64, mut lo: usize) -> usize {
    let mut hi = a.len();
    if lo > hi {
        lo = hi;
    }
    while lo < hi {
        let mid = (lo + hi) / 2;
        if a[mid] <= x {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    lo
}

/// Python `str` for finite values that Rust would print without a decimal point.
pub(crate) fn py_float(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x.is_sign_negative() {
            "-inf".to_string()
        } else {
            "inf".to_string()
        };
    }
    let s = format!("{x}");
    if s.contains(['.', 'e', 'E']) {
        s
    } else {
        format!("{s}.0")
    }
}

pub(crate) fn insert_sorted(values: &mut Vec<usize>, value: usize) {
    if let Err(index) = values.binary_search(&value) {
        values.insert(index, value);
    }
}
