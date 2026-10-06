//! `detect_polya` and `reverse_complement` from `tama_collapse.py`.
//!
//! The window is taken with the transcript coordinate used directly as a
//! Python list index. On `+` that is `fasta[end:end+window]`. On `-` it is
//! `fasta[start-window:start]`, then the reverse complement. An empty window
//! is still reported with length 1 so the percentage division matches.

pub struct PolyA {
    pub sequence: String,
    /// Length used in the percentage. Empty sequence becomes 1.
    pub length: i64,
    pub a_count: i64,
    pub n_count: i64,
    /// Fraction of A, not yet multiplied by 100.
    pub a_percent: f64,
    pub n_percent: f64,
}

pub fn reverse_complement(seq: &str) -> Result<String, String> {
    let mut out = String::with_capacity(seq.len());
    for base in seq.chars().rev() {
        let comp = match base {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'N' | 'R' | 'Y' | 'K' | 'M' | 'S' | 'W' | 'B' | 'D' | 'H' | 'V' => 'N',
            other => return Err(format!("unknown base in reverse complement {other}")),
        };
        out.push(comp);
    }
    Ok(out)
}

pub fn detect_polya(
    strand: &str,
    genome: &str,
    start_pos: i64,
    end_pos: i64,
    window: i64,
) -> Result<PolyA, String> {
    let sequence = if strand == "+" {
        slice_clamped(genome, end_pos, end_pos + window).to_string()
    } else if strand == "-" {
        let mut window_start = start_pos - window;
        if window_start < 0 {
            window_start = 0;
        }
        let raw = slice_clamped(genome, window_start, start_pos);
        reverse_complement(raw)?
    } else {
        return Err(format!(
            "Error with strand information for poly a detection {strand}"
        ));
    };
    let mut length = sequence.len() as i64;
    if length == 0 {
        length = 1;
    }
    let a_count = sequence.bytes().filter(|base| *base == b'A').count() as i64;
    let n_count = sequence.bytes().filter(|base| *base == b'N').count() as i64;
    Ok(PolyA {
        sequence,
        length,
        a_count,
        n_count,
        a_percent: a_count as f64 / length as f64,
        n_percent: n_count as f64 / length as f64,
    })
}

fn slice_clamped(text: &str, start: i64, end: i64) -> &str {
    let len = text.len() as i64;
    let start = start.max(0).min(len) as usize;
    let end = end.max(start as i64).min(len) as usize;
    &text[start..end]
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::trans_coordinates;
    use plena_io::{read_fasta, read_sam, SamClass};
    use std::collections::HashMap;
    use std::path::Path;

    fn repo(rel: &str) -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
    }

    #[test]
    fn plus_window_uses_end_as_a_list_index() {
        let genome = "NNNAAAAAAAAAAAAAAAAACCC";
        let polya = detect_polya("+", genome, 1, 3, 20).unwrap();
        assert_eq!(polya.sequence, "AAAAAAAAAAAAAAAAACCC");
        assert_eq!(polya.a_count, 17);
        assert_eq!(polya.length, 20);
    }

    #[test]
    fn gmap_polya_read_matches_oracle_sequence() {
        let genomes: HashMap<_, _> =
            read_fasta(&repo("../../tests/parity/gmap_collapse/test_genome.fa"))
                .unwrap()
                .into_iter()
                .map(|rec| (rec.id, rec.seq))
                .collect();
        let sam = read_sam(&repo("../../tests/parity/gmap_collapse/gmap_test.sam")).unwrap();
        let rec = sam
            .iter()
            .find(|rec| rec.qname == "11_c69636/1/797")
            .unwrap();
        assert_eq!(rec.class, SamClass::Reverse);
        let coords = trans_coordinates(rec.pos, &rec.cigar).unwrap();
        let polya = detect_polya("-", &genomes[&rec.rname], rec.pos, coords.end, 20).unwrap();
        assert_eq!(polya.sequence, "AAGGAAAGAGGAAAAAAAAA");
        assert_eq!(polya.a_count, 15);
        assert_eq!(crate::py2_str_round(polya.a_percent * 100.0), "75.0");
    }
}
