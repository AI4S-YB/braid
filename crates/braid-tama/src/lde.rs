//! Splice-junction local density from `length_error_type`, `simple_sj_error`,
//! and `sj_error_local_density`.
//!
//! `simple_sj_error` reverses the pre-junction tokens, pads out to
//! `sj_err_threshold` with the match character, then reverses that string.
//! A lone `"0"` expands to `sj_err_threshold` match characters and does not
//! advance the running position counter.

use crate::cigar_parts;

pub struct LdeRead<'a> {
    pub cluster_id: &'a str,
    pub scaffold: &'a str,
    pub start_pos: i64,
    pub end_pos: i64,
    pub strand: &'a str,
    pub exon_count: usize,
    pub cigar: &'a str,
    pub sj_pre: &'a [String],
    pub sj_post: &'a [String],
}

pub struct LdeResult {
    pub bad_sj_flag: i64,
    pub line: String,
}

pub fn local_density(
    read: &LdeRead<'_>,
    sj_err_threshold: i64,
    lde_threshold: i64,
    ses_match: char,
) -> Result<LdeResult, String> {
    let mut bad_sj_flag = 0i64;
    let mut bad_sj_nums = Vec::new();
    let mut bad_counts = Vec::new();
    let mut profiles = Vec::new();
    let mut nucs = Vec::new();
    let mut simples = Vec::new();
    let junctions = read.exon_count.saturating_sub(1);

    for i in 0..junctions {
        let pre = &read.sj_pre[i];
        let post = &read.sj_post[i];
        let pre_tokens: Vec<&str> = pre.split('_').collect();
        let post_tokens: Vec<&str> = post.split('_').collect();
        let pre_simple = simple_side(&pre_tokens, true, sj_err_threshold, ses_match)?;
        let post_simple = simple_side(&post_tokens, false, sj_err_threshold, ses_match)?;
        let (pre_i, pre_d, pre_m, pre_s, pre_h) = error_counts(&pre_tokens)?;
        let (post_i, post_d, post_m, post_s, post_h) = error_counts(&post_tokens)?;
        let pre_all = pre_i + pre_d + pre_m + pre_s + pre_h;
        let post_all = post_i + post_d + post_m + post_s + post_h;

        nucs.push(format!("{pre}>{post}"));
        simples.push(format!("{pre_simple}>{post_simple}"));
        profiles.push(format!(
            "{pre_i},{pre_d},{pre_m},{pre_s},{pre_h}>{post_i},{post_d},{post_m},{post_s},{post_h}"
        ));
        bad_counts.push(format!("{pre_all}>{post_all}"));

        let sj_num = if read.strand == "+" {
            i as i64 + 1
        } else if read.strand == "-" {
            read.exon_count as i64 - i as i64 - 1
        } else {
            return Err(format!("bad strand {}", read.strand));
        };
        let mut this_bad = false;
        if pre_all > lde_threshold {
            bad_sj_flag += 1;
            this_bad = true;
        }
        if post_all > lde_threshold {
            bad_sj_flag += 1;
            this_bad = true;
        }
        if this_bad {
            bad_sj_nums.push(sj_num.to_string());
        }
    }

    let flag = if bad_sj_flag == 0 {
        "lde_pass"
    } else {
        "lde_fail"
    };
    let line = [
        read.cluster_id,
        flag,
        read.scaffold,
        &read.start_pos.to_string(),
        &read.end_pos.to_string(),
        read.strand,
        &read.exon_count.to_string(),
        &if bad_sj_nums.is_empty() {
            "0".to_string()
        } else {
            bad_sj_nums.join(",")
        },
        &if bad_counts.is_empty() {
            "na".to_string()
        } else {
            bad_counts.join(",")
        },
        &if profiles.is_empty() {
            "na".to_string()
        } else {
            profiles.join(";")
        },
        &if nucs.is_empty() {
            "na".to_string()
        } else {
            nucs.join(";")
        },
        &if simples.is_empty() {
            "na".to_string()
        } else {
            simples.join(";")
        },
        read.cigar,
    ]
    .join("\t");

    Ok(LdeResult { bad_sj_flag, line })
}

fn error_counts(tokens: &[&str]) -> Result<(i64, i64, i64, i64, i64), String> {
    let mut indel = [0i64; 5];
    for token in tokens {
        let (length, kind) = length_error_type(token)?;
        let slot = match kind {
            "mismatch" => 2,
            "I" => 0,
            "D" => 1,
            "S" => 3,
            "H" => 4,
            "0" | "M" => continue,
            other => return Err(format!("unexpected error type {other}")),
        };
        indel[slot] += length;
    }
    Ok((indel[0], indel[1], indel[2], indel[3], indel[4]))
}

fn length_error_type(token: &str) -> Result<(i64, &'static str), String> {
    let kind = if token.split('.').count() == 3 {
        "mismatch"
    } else if token == "0" {
        "0"
    } else if token.contains('I') {
        "I"
    } else if token.contains('D') {
        "D"
    } else if token.contains('M') {
        "M"
    } else if token.contains('S') {
        "S"
    } else if token.contains('H') {
        "H"
    } else {
        return Err(format!("Error with error_string: {token}"));
    };
    let length = match kind {
        "0" => 0,
        "mismatch" => 1,
        letter => {
            let bits: Vec<&str> = token.split(letter).collect();
            if bits.len() != 2 {
                return Err(format!("Error with error_string: {token}"));
            }
            bits[0]
                .parse()
                .map_err(|_| format!("bad error length in {token}"))?
        }
    };
    Ok((length, kind))
}

fn simple_side(
    tokens: &[&str],
    reverse: bool,
    threshold: i64,
    ses: char,
) -> Result<String, String> {
    let mut ordered: Vec<&str> = tokens.to_vec();
    if reverse {
        ordered.reverse();
    }
    let mut simple = Vec::new();
    let mut count = 1i64;
    for token in ordered {
        if token.split('.').count() == 3 {
            push_mismatch(&mut simple, &mut count, token, ses);
        } else if token.split('.').count() == 1 {
            push_op(&mut simple, &mut count, token, threshold, ses)?;
        } else {
            return Err(format!("Error with LDE error char {token}"));
        }
    }
    while (simple.len() as i64) < threshold {
        simple.push(ses);
    }
    if reverse {
        simple.reverse();
    }
    Ok(simple.into_iter().collect())
}

fn push_mismatch(simple: &mut Vec<char>, count: &mut i64, token: &str, ses: char) {
    let mismatch_position = token.split('.').next().unwrap().parse::<i64>().unwrap_or(0) + 1;
    if *count == 1 {
        for _ in 0..(mismatch_position - 1) {
            simple.push(ses);
            *count += 1;
        }
        simple.push('X');
        *count += 1;
    } else if *count > 1 {
        let diff = mismatch_position - *count;
        for _ in 0..diff {
            simple.push(ses);
            *count += 1;
        }
        simple.push('X');
        *count += 1;
    }
}

fn push_op(
    simple: &mut Vec<char>,
    count: &mut i64,
    token: &str,
    threshold: i64,
    ses: char,
) -> Result<(), String> {
    if token == "0" {
        for _ in 0..threshold {
            simple.push(ses);
        }
        return Ok(());
    }
    let parts = cigar_parts(token)?;
    if parts.lengths.is_empty() {
        return Err(format!("Error with LDE error char {token}"));
    }
    let dig = i64::from(parts.lengths[0]);
    let symbol = match parts.ops[0].as_str() {
        "M" => ses,
        "I" => 'I',
        "D" => 'D',
        "S" => 'S',
        "H" => 'H',
        other => return Err(format!("Error with LDE error char {other}")),
    };
    for _ in 0..dig {
        simple.push(symbol);
        *count += 1;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::trans_coordinates;
    use crate::{calc_error_rate, coverage_percent, identity_percent, VariationBook};
    use braid_io::{read_fasta, read_sam, SamClass};
    use std::collections::HashMap;
    use std::path::Path;

    fn repo(rel: &str) -> std::path::PathBuf {
        Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)
    }

    #[test]
    fn gmap_lde_lines_match_oracle() {
        let genomes: HashMap<_, _> =
            read_fasta(&repo("../../tests/parity/gmap_collapse/test_genome.fa"))
                .unwrap()
                .into_iter()
                .map(|rec| (rec.id, rec.seq))
                .collect();
        let sam = read_sam(&repo("../../tests/parity/gmap_collapse/gmap_test.sam")).unwrap();
        let lde_txt = std::fs::read_to_string(repo(
            "../../tests/parity/gmap_collapse/gmap_local_density_error.txt",
        ))
        .unwrap();
        let read_txt =
            std::fs::read_to_string(repo("../../tests/parity/gmap_collapse/gmap_read.txt"))
                .unwrap();
        let mut expected_lde = lde_txt.lines().skip(1);
        let mut read_lines = read_txt.lines().skip(1);
        let mut produced = 0usize;
        for rec in &sam {
            let read_line = read_lines.next().unwrap();
            let accept = read_line.split('\t').nth(2).unwrap();
            if !rec.class.keeps_model() {
                continue;
            }
            let genome = &genomes[&rec.rname];
            let mut book = VariationBook::new();
            let rate = calc_error_rate(
                rec.pos, &rec.cigar, &rec.seq, &rec.rname, &rec.qname, genome, 10, &mut book,
            )
            .unwrap();
            let seq_length = rec.seq.len() as i64 + rate.h_count;
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
            if cov < 99.0 || ident < 85.0 {
                assert_eq!(accept, "discarded", "{}", rec.qname);
                continue;
            }
            let coords = trans_coordinates(rec.pos, &rec.cigar).unwrap();
            let strand = match rec.class {
                SamClass::Forward => "+",
                SamClass::Reverse => "-",
                _ => unreachable!(),
            };
            let lde = local_density(
                &LdeRead {
                    cluster_id: &rec.qname,
                    scaffold: &rec.rname,
                    start_pos: rec.pos,
                    end_pos: coords.end,
                    strand,
                    exon_count: coords.starts.len(),
                    cigar: &rec.cigar,
                    sj_pre: &rate.sj_pre,
                    sj_post: &rate.sj_post,
                },
                10,
                1000,
                '_',
            )
            .unwrap();
            let expected = expected_lde
                .next()
                .unwrap_or_else(|| panic!("missing {}", rec.qname));
            assert_eq!(lde.line, expected, "{}", rec.qname);
            let want = if lde.bad_sj_flag > 0 {
                "local_density_error"
            } else {
                "accepted"
            };
            assert_eq!(accept, want, "{}", rec.qname);
            produced += 1;
        }
        assert!(expected_lde.next().is_none());
        assert!(produced > 0);
    }
}
