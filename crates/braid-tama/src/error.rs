//! `mismatch_seq` and `calc_error_rate` from `tama_collapse.py`.
//!
//! The IUPAC table is copied as written. Several letters do not match the
//! standard code (`K` is A/C, `M` is G/T, and so on). Splice-window scans
//! keep the same off-by-threshold steps: an operator that pushes the running
//! length past `sj_err_threshold` is still inspected, and the pre-junction
//! token list is reversed before it is joined.

use crate::py27::Py27Dict;

const NO_MISMATCH: &str = "0";

/// Nonstandard on purpose. Do not "correct" it.
fn iupac(base: u8) -> Option<&'static [u8]> {
    Some(match base {
        b'A' => b"A",
        b'T' => b"T",
        b'C' => b"C",
        b'G' => b"G",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"AC",
        b'M' => b"GT",
        b'Y' => b"AG",
        b'R' => b"CT",
        b'V' => b"CGT",
        b'H' => b"AGT",
        b'D' => b"ACT",
        b'B' => b"ACG",
        b'N' => b"ATCG",
        _ => return None,
    })
}

fn bases_match(genome: u8, read: u8) -> Result<bool, String> {
    let read = (read as char).to_ascii_uppercase() as u8;
    let genome_set =
        iupac(genome).ok_or_else(|| format!("unknown genome base {}", genome as char))?;
    let read_set = iupac(read).ok_or_else(|| format!("unknown read base {}", read as char))?;
    Ok(genome_set.iter().any(|base| read_set.contains(base)))
}

fn slice_ascii(text: &str, start: i64, end: i64) -> Result<&str, String> {
    if start < 0 || end < start || end > text.len() as i64 {
        return Err(format!(
            "slice {start}:{end} outside sequence of length {}",
            text.len()
        ));
    }
    Ok(&text[start as usize..end as usize])
}

fn mismatch_seq(
    genome: &str,
    query: &str,
    genome_pos: i64,
    seq_pos: i64,
) -> Result<(Vec<i64>, Vec<i64>, Vec<String>), String> {
    if genome.len() != query.len() {
        return Err(format!(
            "Genome seq is not the same length as query seq ({genome_pos} {seq_pos})"
        ));
    }
    let mut genome_mismatches = Vec::new();
    let mut seq_mismatches = Vec::new();
    let mut nuc_mismatches = Vec::new();
    let genome = genome.as_bytes();
    let query_bytes = query.as_bytes();
    for i in 0..genome.len() {
        if !bases_match(genome[i], query_bytes[i])? {
            genome_mismatches.push(genome_pos + i as i64);
            seq_mismatches.push(seq_pos + i as i64);
            nuc_mismatches.push(format!("{}.{}", query_bytes[i] as char, genome[i] as char));
        }
    }
    Ok((genome_mismatches, seq_mismatches, nuc_mismatches))
}

/// `variation_dict` and `var_coverage_dict`. Nested maps use CPython 2.7 order.
#[derive(Clone, Debug)]
pub struct VariationBook {
    pub sites:
        Py27Dict<String, Py27Dict<i64, Py27Dict<String, Py27Dict<String, Py27Dict<String, i32>>>>>,
    pub coverage: Py27Dict<String, Py27Dict<i64, Py27Dict<String, i32>>>,
}

impl VariationBook {
    pub fn new() -> Self {
        Self {
            sites: Py27Dict::new(),
            coverage: Py27Dict::new(),
        }
    }

    pub fn update(&mut self, scaffold: &str, pos: i64, kind: &str, seq: &str, read_id: &str) {
        self.sites
            .or_insert_with(scaffold.to_string(), Py27Dict::new)
            .or_insert_with(pos, Py27Dict::new)
            .or_insert_with(kind.to_string(), Py27Dict::new)
            .or_insert_with(seq.to_string(), Py27Dict::new)
            .insert(read_id.to_string(), 1);
        self.coverage
            .or_insert_with(scaffold.to_string(), Py27Dict::new)
            .or_insert_with(pos, Py27Dict::new)
            .insert(read_id.to_string(), 1);
    }
}

impl Default for VariationBook {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Clone, Debug)]
pub struct ErrorRate {
    pub h_count: i64,
    pub s_count: i64,
    pub i_count: i64,
    pub d_count: i64,
    pub mis_count: i64,
    pub insertions: Vec<i64>,
    pub insertion_lengths: Vec<i64>,
    pub deletions: Vec<i64>,
    pub deletion_lengths: Vec<i64>,
    pub mismatches: Vec<i64>,
    pub nuc_mismatches: Vec<String>,
    pub sj_pre: Vec<String>,
    pub sj_post: Vec<String>,
}

impl ErrorRate {
    pub fn error_line(&self) -> String {
        format!(
            "{};{};{};{};{}",
            self.h_count, self.s_count, self.i_count, self.d_count, self.mis_count
        )
    }

    /// `pre>post` tokens joined the way `_local_density_error.txt` writes them.
    pub fn sj_nuc_line(&self) -> String {
        if self.sj_pre.is_empty() {
            return "na".to_string();
        }
        self.sj_pre
            .iter()
            .zip(&self.sj_post)
            .map(|(pre, post)| format!("{pre}>{post}"))
            .collect::<Vec<_>>()
            .join(";")
    }
}

pub fn coverage_percent(seq_length: i64, h_count: i64, s_count: i64) -> f64 {
    (seq_length - h_count - s_count) as f64 / seq_length as f64 * 100.0
}

/// `ident_cov` counts clips against the full read. `ident_map` ignores them
/// and divides by the mapped length.
pub fn identity_percent(
    seq_length: i64,
    h_count: i64,
    s_count: i64,
    i_count: i64,
    d_count: i64,
    mis_count: i64,
    method: &str,
) -> Result<f64, String> {
    if method == "ident_cov" {
        let nonmatch = h_count + s_count + i_count + d_count + mis_count;
        Ok((seq_length - nonmatch) as f64 / seq_length as f64 * 100.0)
    } else if method == "ident_map" {
        let mapped = seq_length - h_count - s_count;
        let nonmatch = i_count + d_count + mis_count;
        Ok((mapped - nonmatch) as f64 / mapped as f64 * 100.0)
    } else {
        Err(format!("unknown identity method {method}"))
    }
}

pub fn calc_error_rate(
    start_pos: i64,
    cigar: &str,
    seq: &str,
    scaffold: &str,
    read_id: &str,
    genome: &str,
    sj_err_threshold: i64,
    variations: &mut VariationBook,
) -> Result<ErrorRate, String> {
    let parts = crate::cigar_parts(cigar)?;
    let lens: Vec<i64> = parts.lengths.iter().map(|n| i64::from(*n)).collect();
    let ops = parts.ops;
    let mut h_count = 0i64;
    let mut s_count = 0i64;
    let mut i_count = 0i64;
    let mut d_count = 0i64;
    let mut mis_count = 0i64;
    let mut mismatches = Vec::new();
    let mut nuc_mismatches = Vec::new();
    let mut insertions = Vec::new();
    let mut insertion_lengths = Vec::new();
    let mut deletions = Vec::new();
    let mut deletion_lengths = Vec::new();
    let mut sj_pre = Vec::new();
    let mut sj_post = Vec::new();
    let mut genome_pos = start_pos - 1;
    let mut seq_pos = 0i64;

    for i in 0..lens.len() {
        let flag = ops[i].as_str();
        let len = lens[i];
        match flag {
            "H" => {
                h_count += len;
                let var_pos = genome_pos - h_count;
                let var_seq = "N".repeat(h_count as usize);
                variations.update(scaffold, var_pos, "H", &var_seq, read_id);
            }
            "S" => {
                s_count += len;
                let seq_start = seq_pos;
                let seq_end = seq_pos + len;
                let mut var_pos = genome_pos - s_count;
                if i + 1 == lens.len() {
                    var_pos = genome_pos + 1 - s_count;
                }
                let var_seq = slice_ascii(seq, seq_start, seq_end)?;
                variations.update(scaffold, var_pos, "S", var_seq, read_id);
                seq_pos += len;
            }
            "M" => {
                let seq_slice = slice_ascii(seq, seq_pos, seq_pos + len)?;
                let genome_slice = slice_ascii(genome, genome_pos, genome_pos + len)?;
                let (genome_mm, seq_mm, nuc_mm) =
                    mismatch_seq(genome_slice, seq_slice, genome_pos, seq_pos)?;
                mis_count += genome_mm.len() as i64;
                for index in 0..seq_mm.len() {
                    let var_seq = (seq.as_bytes()[seq_mm[index] as usize] as char).to_string();
                    variations.update(scaffold, genome_mm[index], "M", &var_seq, read_id);
                }
                mismatches.extend(genome_mm);
                nuc_mismatches.extend(nuc_mm);
                seq_pos += len;
                genome_pos += len;
            }
            "I" => {
                let var_seq = slice_ascii(seq, seq_pos, seq_pos + len)?;
                variations.update(scaffold, genome_pos, "I", var_seq, read_id);
                seq_pos += len;
                i_count += len;
                insertions.push(genome_pos);
                insertion_lengths.push(len);
            }
            "D" => {
                deletions.push(genome_pos);
                deletion_lengths.push(len);
                variations.update(scaffold, genome_pos, "D", &len.to_string(), read_id);
                genome_pos += len;
                d_count += len;
            }
            "N" => {
                if i == 0 || i + 1 >= ops.len() {
                    return Err(format!("CIGAR N at the edge: {cigar}"));
                }
                let pre = scan_pre(
                    &ops,
                    &lens,
                    i,
                    genome_pos,
                    &mismatches,
                    &nuc_mismatches,
                    sj_err_threshold,
                );
                genome_pos += len;
                let post = scan_post(
                    &ops,
                    &lens,
                    i,
                    genome_pos,
                    seq,
                    seq_pos,
                    genome,
                    sj_err_threshold,
                )?;
                sj_pre.push(pre);
                sj_post.push(post);
            }
            other => {
                return Err(format!("Error with cigar flag {len}{other}"));
            }
        }
    }

    Ok(ErrorRate {
        h_count,
        s_count,
        i_count,
        d_count,
        mis_count,
        insertions,
        insertion_lengths,
        deletions,
        deletion_lengths,
        mismatches,
        nuc_mismatches,
        sj_pre,
        sj_post,
    })
}

fn scan_pre(
    ops: &[String],
    lens: &[i64],
    n_index: usize,
    genome_pos: i64,
    mismatches: &[i64],
    nucs: &[String],
    sj_err_threshold: i64,
) -> String {
    let mut prev_index = n_index as isize - 1;
    let mut prev_flag = ops[prev_index as usize].clone();
    let mut prev_len = lens[prev_index as usize];
    let mut prev_total = 0i64;
    let mut window = 0i64;
    let mut mismatch_index = mismatches.len() as isize - 1;
    let mut this_genome_pos = genome_pos;
    let mut hit_intron = false;
    let mut builder = Vec::new();

    while prev_total <= sj_err_threshold && !hit_intron {
        window += prev_len;
        prev_total = window;
        if prev_flag == "M" {
            if mismatches.is_empty() {
                if prev_total <= sj_err_threshold {
                    builder.push(format!("{prev_len}{prev_flag}"));
                }
            } else if mismatch_index >= 0 {
                let mut last_pos = mismatches[mismatch_index as usize];
                let mut last_nuc = nucs[mismatch_index as usize].clone();
                let mut dist = genome_pos - last_pos;
                let mut added = 0;
                while last_pos >= this_genome_pos - prev_len
                    && last_pos <= this_genome_pos
                    && dist <= sj_err_threshold
                {
                    builder.push(format!("{dist}.{last_nuc}"));
                    added += 1;
                    mismatch_index -= 1;
                    if mismatch_index < 0 {
                        break;
                    }
                    last_pos = mismatches[mismatch_index as usize];
                    last_nuc = nucs[mismatch_index as usize].clone();
                    dist = genome_pos - last_pos;
                }
                if added == 0 && prev_total <= sj_err_threshold {
                    builder.push(format!("{prev_len}{prev_flag}"));
                }
            } else if prev_total <= sj_err_threshold {
                builder.push(format!("{prev_len}{prev_flag}"));
            }
        } else if prev_flag != "N" {
            builder.push(format!("{prev_len}{prev_flag}"));
        } else {
            hit_intron = true;
        }

        prev_index -= 1;
        if prev_index < 0 {
            break;
        }
        this_genome_pos -= prev_len;
        prev_flag = ops[prev_index as usize].clone();
        prev_len = lens[prev_index as usize];
    }

    if builder.is_empty() {
        builder.push(NO_MISMATCH.to_string());
    }
    builder.reverse();
    builder.join("_")
}

fn scan_post(
    ops: &[String],
    lens: &[i64],
    n_index: usize,
    genome_pos: i64,
    seq: &str,
    seq_pos: i64,
    genome: &str,
    sj_err_threshold: i64,
) -> Result<String, String> {
    let mut next_index = n_index + 1;
    let mut next_flag = ops[next_index].clone();
    let mut next_len = lens[next_index];
    let mut next_total = 0i64;
    let mut window = 0i64;
    let mut seq_cursor = seq_pos;
    let mut genome_cursor = genome_pos;
    let mut hit_intron = false;
    let mut builder = Vec::new();

    while next_total <= sj_err_threshold && !hit_intron {
        window += next_len;
        next_total = window;
        if next_flag == "M" {
            let seq_slice = slice_ascii(seq, seq_cursor, seq_cursor + next_len)?;
            let genome_slice = slice_ascii(genome, genome_cursor, genome_cursor + next_len)?;
            let (genome_mm, _seq_mm, nuc_mm) =
                mismatch_seq(genome_slice, seq_slice, genome_cursor, seq_cursor)?;
            if genome_mm.is_empty() {
                if next_total <= sj_err_threshold {
                    builder.push(format!("{next_len}{next_flag}"));
                }
            } else {
                let mut index = 0usize;
                let mut dist = genome_mm[0] - genome_pos;
                let mut added = 0;
                while dist <= window && dist < sj_err_threshold {
                    let pos = genome_mm[index];
                    dist = pos - genome_pos;
                    if dist <= sj_err_threshold {
                        builder.push(format!("{dist}.{}", nuc_mm[index]));
                        added += 1;
                    }
                    index += 1;
                    if index >= genome_mm.len() {
                        break;
                    }
                }
                if added == 0 && next_total <= sj_err_threshold {
                    builder.push(format!("{next_len}{next_flag}"));
                }
            }
        } else if next_flag != "N" {
            builder.push(format!("{next_len}{next_flag}"));
        } else {
            hit_intron = true;
        }

        if next_flag != "D" {
            seq_cursor += next_len;
        }
        if next_flag != "I" {
            genome_cursor += next_len;
        }
        next_index += 1;
        if next_index >= ops.len() {
            break;
        }
        next_flag = ops[next_index].clone();
        next_len = lens[next_index];
    }

    if builder.is_empty() {
        builder.push(NO_MISMATCH.to_string());
    }
    Ok(builder.join("_"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use braid_io::{read_fasta, read_sam, SamClass};
    use std::collections::HashMap;
    use std::path::Path;

    fn repo(rel: &str) -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
    }

    #[test]
    fn iupac_table_is_the_tama_table() {
        assert!(bases_match(b'A', b'a').unwrap());
        assert!(!bases_match(b'A', b'T').unwrap());
        assert!(bases_match(b'N', b'C').unwrap());
        assert!(bases_match(b'K', b'A').unwrap());
        assert!(!bases_match(b'K', b'G').unwrap());
    }

    #[test]
    fn gmap_error_counts_and_splice_windows_match_oracle() {
        let genome_records =
            read_fasta(&repo("../../tests/parity/gmap_collapse/test_genome.fa")).unwrap();
        let mut genomes = HashMap::new();
        for record in genome_records {
            genomes.insert(record.id, record.seq);
        }
        let sam = read_sam(&repo("../../tests/parity/gmap_collapse/gmap_test.sam")).unwrap();
        let read_txt =
            std::fs::read_to_string(repo("../../tests/parity/gmap_collapse/gmap_read.txt"))
                .unwrap();
        let lde_txt = std::fs::read_to_string(repo(
            "../../tests/parity/gmap_collapse/gmap_local_density_error.txt",
        ))
        .unwrap();
        let mut lde_nuc: HashMap<String, String> = HashMap::new();
        for line in lde_txt.lines().skip(1) {
            let mut cols = line.split('\t');
            let id = cols.next().unwrap().to_string();
            for _ in 0..9 {
                cols.next().unwrap();
            }
            lde_nuc.insert(id, cols.next().unwrap().to_string());
        }

        let mut read_lines = read_txt.lines().skip(1);
        let mut checked = 0usize;
        let mut splice_checked = 0usize;
        for rec in &sam {
            let line = read_lines.next().expect("read.txt shorter than SAM");
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols[0], rec.qname);
            let mapped = match rec.class {
                SamClass::Forward => "forward_strand",
                SamClass::Reverse => "reverse_strand",
                SamClass::Unmapped => "unmapped",
                SamClass::Chimeric => "chimeric",
                SamClass::NotPrimary => "not_primary",
            };
            assert_eq!(cols[1], mapped, "{}", rec.qname);
            if !rec.class.keeps_model() {
                assert_eq!(cols[5], "NA", "{}", rec.qname);
                continue;
            }
            let genome = genomes.get(&rec.rname).expect("scaffold");
            let mut book = VariationBook::new();
            let rate = calc_error_rate(
                rec.pos, &rec.cigar, &rec.seq, &rec.rname, &rec.qname, genome, 10, &mut book,
            )
            .unwrap_or_else(|err| panic!("{}: {err}", rec.qname));
            let seq_length = rec.seq.len() as i64 + rate.h_count;
            assert_eq!(rate.error_line(), cols[5], "{}", rec.qname);
            assert_eq!(seq_length.to_string(), cols[6], "{}", rec.qname);
            assert_eq!(cols[7], rec.cigar, "{}", rec.qname);
            let cov = coverage_percent(seq_length, rate.h_count, rate.s_count);
            let ident = identity_percent(
                seq_length,
                rate.h_count,
                rate.s_count,
                rate.i_count,
                rate.d_count,
                rate.mis_count,
                "ident_cov",
            )
            .unwrap();
            assert_eq!(crate::py2_str_round(cov), cols[3], "{}", rec.qname);
            assert_eq!(crate::py2_str_round(ident), cols[4], "{}", rec.qname);
            if cov < 99.0 || ident < 85.0 {
                assert_eq!(cols[2], "discarded", "{}", rec.qname);
            } else {
                assert!(
                    cols[2] == "accepted" || cols[2] == "local_density_error",
                    "{} {}",
                    rec.qname,
                    cols[2]
                );
                assert_eq!(rate.sj_nuc_line(), lde_nuc[&rec.qname], "{}", rec.qname);
                splice_checked += 1;
            }
            checked += 1;
        }
        assert!(read_lines.next().is_none());
        assert!(checked > 0);
        assert!(splice_checked > 0);
    }
}
