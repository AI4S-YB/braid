//! `cigar_list`, `mapped_seq_length`, and `trans_coordinates` from
//! `tama_collapse.py`. `=` and `X` become `M`. Operators other than
//! H, S, M, I, D, N are ignored, including `P`.

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CigarParts {
    pub lengths: Vec<u32>,
    pub ops: Vec<String>,
}

pub fn cigar_parts(cigar: &str) -> Result<CigarParts, String> {
    let cigar = cigar.replace('=', "M").replace('X', "M");
    let op_blob: String = cigar
        .chars()
        .map(|c| if c.is_ascii_digit() { ' ' } else { c })
        .collect();
    let len_blob: String = replace_letter_runs(&cigar);
    let ops: Vec<String> = op_blob.split_whitespace().map(str::to_string).collect();
    let mut lengths = Vec::new();
    for token in len_blob.split_whitespace() {
        let len: u32 = token
            .parse()
            .map_err(|_| format!("bad CIGAR length '{token}' in {cigar}"))?;
        lengths.push(len);
    }
    if lengths.len() > ops.len() {
        return Err(format!("CIGAR has more lengths than operators: {cigar}"));
    }
    Ok(CigarParts { lengths, ops })
}

fn replace_letter_runs(cigar: &str) -> String {
    let mut out = String::with_capacity(cigar.len());
    let mut in_run = false;
    for c in cigar.chars() {
        if c.is_ascii_alphabetic() {
            if !in_run {
                out.push(' ');
                in_run = true;
            }
        } else {
            in_run = false;
            out.push(c);
        }
    }
    out
}

pub fn mapped_seq_length(cigar: &str) -> Result<u32, String> {
    let parts = cigar_parts(cigar)?;
    let mut total = 0u32;
    for (i, len) in parts.lengths.iter().enumerate() {
        match parts.ops[i].as_str() {
            "M" | "D" => total += *len,
            _ => {}
        }
    }
    Ok(total)
}

pub struct ExonCoords {
    pub end: i64,
    pub starts: Vec<i64>,
    pub ends: Vec<i64>,
}

pub fn trans_coordinates(start_pos: i64, cigar: &str) -> Result<ExonCoords, String> {
    let parts = cigar_parts(cigar)?;
    let mut end_pos = start_pos;
    let mut starts = vec![start_pos];
    let mut ends = Vec::new();
    for (i, len) in parts.lengths.iter().enumerate() {
        let len = i64::from(*len);
        match parts.ops[i].as_str() {
            "M" | "D" => end_pos += len,
            "N" => {
                ends.push(end_pos);
                end_pos += len;
                starts.push(end_pos);
            }
            _ => {}
        }
    }
    ends.push(end_pos);
    Ok(ExonCoords {
        end: end_pos,
        starts,
        ends,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn parity(name: &str) -> String {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/parity")
            .join(name);
        std::fs::read_to_string(&path)
            .unwrap_or_else(|err| panic!("read {}: {err}", path.display()))
    }

    fn assert_golden(text: &str) {
        let mut rows = 0;
        for line in text.lines() {
            if line.is_empty() {
                continue;
            }
            rows += 1;
            let mut cols = line.split('\t');
            let start: i64 = cols.next().unwrap().parse().unwrap();
            let cigar = cols.next().unwrap();
            let map_len: u32 = cols.next().unwrap().parse().unwrap();
            let end: i64 = cols.next().unwrap().parse().unwrap();
            let starts: Vec<i64> = cols
                .next()
                .unwrap()
                .split(',')
                .map(|n| n.parse().unwrap())
                .collect();
            let ends: Vec<i64> = cols
                .next()
                .unwrap()
                .split(',')
                .map(|n| n.parse().unwrap())
                .collect();
            assert_eq!(mapped_seq_length(cigar).unwrap(), map_len, "{cigar}");
            let coords = trans_coordinates(start, cigar).unwrap();
            assert_eq!(coords.end, end, "{cigar}");
            assert_eq!(coords.starts, starts, "{cigar}");
            assert_eq!(coords.ends, ends, "{cigar}");
        }
        assert!(rows > 0);
    }

    #[test]
    fn matches_python2_synthetic_cigars() {
        assert_golden(&parity("cigar_synthetic.tsv"));
    }

    #[test]
    fn matches_python2_gmap_test_sam() {
        assert_golden(&parity("cigar_gmap.tsv"));
    }
}
