//! Original-mode TAMA collapse.
//!
//! `low_mem` and `no_cap` are refused. Variant support is forced back to 5
//! before the variant files and the report are written, which is what
//! `tama_collapse.py` does in original mode.

use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

use braid_io::{read_fasta, read_sam, SamClass};

use crate::group::{collapse_locus, GroupParams, ReadModel};
use crate::lde::{local_density, LdeRead};
use crate::py2float::py2_str_round;
use crate::{
    calc_error_rate, coverage_percent, detect_polya, identity_percent, mapped_seq_length,
    trans_coordinates, VariationBook,
};

#[derive(Clone, Debug)]
pub struct CollapseSettings {
    pub cap: String,
    pub ends: String,
    pub coverage: f64,
    pub identity: f64,
    pub ident_method: String,
    pub five_prime: i64,
    pub exon_diff: i64,
    pub three_prime: i64,
    pub duplicates: String,
    pub sj_priority: String,
    pub sj_threshold: i64,
    pub lde: i64,
    pub ses: String,
    pub bam: String,
    pub log: String,
    pub run_mode: String,
    /// Value parsed from `-vc`. Original mode overwrites this with 5.
    pub var_support: i64,
    pub sam_label: String,
    pub fasta_label: String,
    pub prefix_label: String,
}

#[derive(Clone, Debug)]
pub struct CollapseTexts {
    pub bed: String,
    pub read_txt: String,
    pub trans_report: String,
    pub trans_read: String,
    pub polya: String,
    pub strand_check: String,
    pub local_density: String,
    pub variants: String,
    pub varcov: String,
    pub report: String,
}

pub fn original_collapse(
    sam_path: &Path,
    fasta_path: &Path,
    settings: &CollapseSettings,
) -> Result<CollapseTexts, String> {
    if settings.run_mode != "original" {
        return Err("only -rm original is ported".to_string());
    }
    if settings.bam != "SAM" {
        return Err("BAM input is not ported".to_string());
    }
    if settings.cap != "capped" {
        return Err("only -x capped is ported".to_string());
    }
    if settings.ends != "common_ends" && settings.ends != "longest_ends" {
        return Err(format!("unknown -e value {}", settings.ends));
    }
    if settings.duplicates != "merge_dup" && settings.duplicates != "no_merge" {
        return Err(format!("unknown -d value {}", settings.duplicates));
    }
    if settings.sj_priority != "no_priority" && settings.sj_priority != "sj_priority" {
        return Err(format!("unknown -sj value {}", settings.sj_priority));
    }
    if settings.ident_method != "ident_cov" && settings.ident_method != "ident_map" {
        return Err(format!("unknown -icm value {}", settings.ident_method));
    }
    let ses = settings
        .ses
        .chars()
        .next()
        .filter(|_| settings.ses.chars().count() == 1)
        .ok_or_else(|| format!("-ses must be one character, got {}", settings.ses))?;

    let fasta = read_fasta(fasta_path)?;
    let mut genome = HashMap::new();
    for record in fasta {
        genome.insert(record.id, record.seq);
    }
    let records = read_sam(sam_path)?;
    let params = GroupParams {
        five_prime: settings.five_prime,
        exon_diff: settings.exon_diff,
        three_prime: settings.three_prime,
        ends: settings.ends.clone(),
        duplicates: settings.duplicates.clone(),
        sj_priority: settings.sj_priority == "sj_priority",
        ident_method: settings.ident_method.clone(),
    };

    let mut book = VariationBook::new();
    let mut read_lines = vec![
        "read_id\tmapped_flag\taccept_flag\tpercent_coverage\tpercent_identity\terror_line<h;s;i;d;m>\tlength\tcigar"
            .to_string(),
    ];
    let mut lde_lines = vec![
        "cluster_id\tlde_flag\tscaff_name\tstart_pos\tend_pos\tstrand\tnum_exons\tbad_sj_num_line\tbad_sj_error_count_line\tsj_error_profile_idmsh\tsj_error_nuc\tsj_error_simple\tcigar"
            .to_string(),
    ];
    let mut strand_lines = vec!["read_id\tscaff_name\tstart_pos\tcigar\tstrands".to_string()];
    let mut accepted = 0i64;
    let mut discarded = 0i64;
    let mut models: HashMap<String, ReadModel> = HashMap::new();
    let mut groups: HashMap<usize, Vec<String>> = HashMap::new();
    let mut trans_group: HashMap<String, usize> = HashMap::new();
    let mut this_scaffold: Option<String> = None;
    let mut group_start = 0i64;
    let mut group_end = 0i64;
    let mut group_count = 0usize;
    let mut scaffolds = Vec::new();

    for record in &records {
        let mapped = mapped_label(record.class);
        if record.class == SamClass::Forward && record.xs == Some('-') {
            strand_lines.push(format!(
                "{}\t{}\t{}\t{}\t+-",
                record.qname, record.rname, record.pos, record.cigar
            ));
        } else if record.class == SamClass::Reverse && record.xs == Some('+') {
            strand_lines.push(format!(
                "{}\t{}\t{}\t{}\t-+",
                record.qname, record.rname, record.pos, record.cigar
            ));
        }
        if !record.class.keeps_model() {
            read_lines.push(format!(
                "{}\t{mapped}\t{mapped}\tNA\tNA\tNA\tNA\tNA",
                record.qname
            ));
            discarded += 1;
            continue;
        }
        let sequence = genome
            .get(&record.rname)
            .ok_or_else(|| format!("scaffold {} is not in the fasta", record.rname))?;
        let _mapped_len = mapped_seq_length(&record.cigar)?;
        let coords = trans_coordinates(record.pos, &record.cigar)?;
        let rate = calc_error_rate(
            record.pos,
            &record.cigar,
            &record.seq,
            &record.rname,
            &record.qname,
            sequence,
            settings.sj_threshold,
            &mut book,
        )?;
        let seq_length = record.seq.len() as i64 + rate.h_count;
        let coverage = coverage_percent(seq_length, rate.h_count, rate.s_count);
        let identity = identity_percent(
            seq_length,
            rate.h_count,
            rate.s_count,
            rate.i_count,
            rate.d_count,
            rate.mis_count,
            &settings.ident_method,
        )?;
        let error_line = format!(
            "{};{};{};{};{}",
            rate.h_count, rate.s_count, rate.i_count, rate.d_count, rate.mis_count
        );
        if coverage < settings.coverage || identity < settings.identity {
            read_lines.push(format!(
                "{}\t{mapped}\tdiscarded\t{}\t{}\t{error_line}\t{seq_length}\t{}",
                record.qname,
                py2_str_round(coverage),
                py2_str_round(identity),
                record.cigar
            ));
            discarded += 1;
            continue;
        }
        let strand = if record.class == SamClass::Forward {
            "+"
        } else {
            "-"
        };
        let lde = local_density(
            &LdeRead {
                cluster_id: &record.qname,
                scaffold: &record.rname,
                start_pos: record.pos,
                end_pos: coords.end,
                strand,
                exon_count: coords.starts.len(),
                cigar: &record.cigar,
                sj_pre: &rate.sj_pre,
                sj_post: &rate.sj_post,
            },
            settings.sj_threshold,
            settings.lde,
            ses,
        )?;
        lde_lines.push(lde.line);
        if lde.bad_sj_flag > 0 {
            read_lines.push(format!(
                "{}\t{mapped}\tlocal_density_error\t{}\t{}\t{error_line}\t{seq_length}\t{}",
                record.qname,
                py2_str_round(coverage),
                py2_str_round(identity),
                record.cigar
            ));
            discarded += 1;
            continue;
        }
        read_lines.push(format!(
            "{}\t{mapped}\taccepted\t{}\t{}\t{error_line}\t{seq_length}\t{}",
            record.qname,
            py2_str_round(coverage),
            py2_str_round(identity),
            record.cigar
        ));
        accepted += 1;
        let polya = detect_polya(strand, sequence, record.pos, coords.end, 20)?;
        let model = ReadModel {
            cluster_id: record.qname.clone(),
            scaff: record.rname.clone(),
            strand: strand.to_string(),
            start_pos: record.pos,
            end_pos: coords.end,
            exon_starts: coords.starts,
            exon_ends: coords.ends,
            sj_pre: rate.sj_pre,
            sj_post: rate.sj_post,
            seq_length,
            h_count: rate.h_count,
            s_count: rate.s_count,
            i_count: rate.i_count,
            d_count: rate.d_count,
            mis_count: rate.mis_count,
            polya_seq: polya.sequence,
            a_count: polya.a_count,
            a_percent: polya.a_percent,
        };
        if models.contains_key(&record.qname) {
            let old = &models[&record.qname];
            let old_cov = coverage_percent(old.seq_length, old.h_count, old.s_count);
            let old_ident = identity_percent(
                old.seq_length,
                old.h_count,
                old.s_count,
                old.i_count,
                old.d_count,
                old.mis_count,
                &settings.ident_method,
            )?;
            if !prefer_new(
                old_cov,
                old_ident,
                coverage,
                identity,
                settings.coverage,
                settings.identity,
            ) {
                continue;
            }
            detach(&record.qname, &mut trans_group, &mut groups)?;
            models.insert(record.qname.clone(), model);
        } else {
            models.insert(record.qname.clone(), model);
        }

        if this_scaffold.is_none() {
            this_scaffold = Some(record.rname.clone());
            group_start = record.pos;
            group_end = coords_end_of(&models[&record.qname]);
            groups.insert(0, vec![record.qname.clone()]);
            trans_group.insert(record.qname.clone(), 0);
            scaffolds.push(record.rname.clone());
            continue;
        }
        let end_pos = models[&record.qname].end_pos;
        if this_scaffold.as_deref() == Some(record.rname.as_str()) {
            if record.pos >= group_start && record.pos <= group_end {
                groups
                    .get_mut(&group_count)
                    .ok_or_else(|| format!("missing locus {group_count}"))?
                    .push(record.qname.clone());
                trans_group.insert(record.qname.clone(), group_count);
                if end_pos > group_end {
                    group_end = end_pos;
                }
            } else if record.pos > group_end {
                group_count += 1;
                group_start = record.pos;
                group_end = end_pos;
                groups.insert(group_count, vec![record.qname.clone()]);
                trans_group.insert(record.qname.clone(), group_count);
            } else {
                return Err(format!("Sam file not sorted! {}", record.qname));
            }
        } else {
            this_scaffold = Some(record.rname.clone());
            group_start = record.pos;
            group_end = end_pos;
            group_count += 1;
            groups.insert(group_count, vec![record.qname.clone()]);
            trans_group.insert(record.qname.clone(), group_count);
            scaffolds.push(record.rname.clone());
        }
    }

    if groups.is_empty() {
        return Err("no groups found".to_string());
    }

    let mut gene_count = 0i64;
    let mut bed = Vec::new();
    let mut trans_report = vec![
        "transcript_id\tnum_clusters\thigh_coverage\tlow_coverage\thigh_quality_percent\tlow_quality_percent\tstart_wobble_list\tend_wobble_list\tcollapse_sj_start_err\tcollapse_sj_end_err\tcollapse_error_nuc"
            .to_string(),
    ];
    let mut trans_read = Vec::new();
    let mut polya = vec!["cluster_id\ttrans_id\tstrand\ta_percent\ta_count\tsequence".to_string()];
    for index in 0..=group_count {
        let Some(ids) = groups.get(&index).cloned() else {
            continue;
        };
        for id in &ids {
            if let Some(read) = models.get(id) {
                stamp_coverage(&mut book, read);
            }
        }
        let (next, genes) = collapse_locus(&ids, &models, &params, gene_count)?;
        gene_count = next;
        for gene in genes {
            bed.extend(gene.bed);
            trans_report.extend(gene.trans_report);
            trans_read.extend(gene.trans_read);
            polya.extend(gene.polya);
        }
    }

    let (variants, varcov) = write_variants(&book, &genome, &scaffolds)?;
    let report = report_text(settings, gene_count, bed.len() as i64, accepted, discarded);
    Ok(CollapseTexts {
        bed: join_lines(&bed),
        read_txt: join_lines(&read_lines),
        trans_report: join_lines(&trans_report),
        trans_read: join_lines(&trans_read),
        polya: join_lines(&polya),
        strand_check: join_lines(&strand_lines),
        local_density: join_lines(&lde_lines),
        variants,
        varcov,
        report,
    })
}

fn coords_end_of(model: &ReadModel) -> i64 {
    model.end_pos
}

fn mapped_label(class: SamClass) -> &'static str {
    match class {
        SamClass::Forward => "forward_strand",
        SamClass::Reverse => "reverse_strand",
        SamClass::Unmapped => "unmapped",
        SamClass::Chimeric => "chimeric",
        SamClass::NotPrimary => "not_primary",
    }
}

fn prefer_new(
    old_cov: f64,
    old_ident: f64,
    new_cov: f64,
    new_ident: f64,
    cov_threshold: f64,
    ident_threshold: f64,
) -> bool {
    let old_pass = old_cov > cov_threshold && old_ident > ident_threshold;
    let new_pass = new_cov > cov_threshold && new_ident > ident_threshold;
    if new_pass && !old_pass {
        return true;
    }
    if old_pass && new_pass {
        return new_cov > old_cov;
    }
    false
}

fn detach(
    read_id: &str,
    trans_group: &mut HashMap<String, usize>,
    groups: &mut HashMap<usize, Vec<String>>,
) -> Result<(), String> {
    let old = trans_group
        .remove(read_id)
        .ok_or_else(|| format!("multimap {read_id} has no locus"))?;
    let members = groups
        .get_mut(&old)
        .ok_or_else(|| format!("missing locus {old}"))?;
    if members.len() > 1 {
        members.retain(|id| id != read_id);
        Ok(())
    } else if members.len() == 1 {
        groups.remove(&old);
        Ok(())
    } else {
        Err(format!("empty locus while moving multimap {read_id}"))
    }
}

fn stamp_coverage(book: &mut VariationBook, read: &ReadModel) {
    let Some(by_pos) = book.coverage.get_mut(&read.scaff) else {
        return;
    };
    for (start, end) in read.exon_starts.iter().zip(&read.exon_ends) {
        for coord in *start..*end {
            if let Some(ids) = by_pos.get_mut(&coord) {
                ids.insert(read.cluster_id.clone(), 1);
            }
        }
    }
}

fn write_variants(
    book: &VariationBook,
    genome: &HashMap<String, String>,
    scaffolds: &[String],
) -> Result<(String, String), String> {
    let mut lines = vec![
        "scaffold\tposition\ttype\tref_allele\talt_allele\tcount\tcov_count\tcluster_list"
            .to_string(),
    ];
    let mut group_order = Vec::new();
    let mut group_pos: HashMap<String, Vec<String>> = HashMap::new();
    for scaffold in scaffolds {
        let Some(by_pos) = book.sites.get(scaffold) else {
            continue;
        };
        let Some(cov_by_pos) = book.coverage.get(scaffold) else {
            return Err(format!("variant coverage is missing {scaffold}"));
        };
        let Some(sequence) = genome.get(scaffold) else {
            return Err(format!("fasta is missing {scaffold}"));
        };
        let mut positions: Vec<i64> = by_pos.iter().map(|(pos, _)| *pos).collect();
        positions.sort();
        for pos in positions {
            if pos < 0 || pos as usize >= sequence.len() {
                continue;
            }
            let Some(cov_reads) = cov_by_pos.get(&pos) else {
                return Err(format!("variant coverage is missing {scaffold}:{pos}"));
            };
            let mut cov_ids: Vec<&str> = cov_reads.iter().map(|(id, _)| id.as_str()).collect();
            cov_ids.sort();
            let cov_count = cov_ids.len();
            let Some(types) = by_pos.get(&pos) else {
                continue;
            };
            let mut ref_allele = (sequence.as_bytes()[pos as usize] as char).to_string();
            let mut accepted = false;
            for kind in ["H", "S", "M", "I", "D"] {
                let Some(alts) = types.get_str(kind) else {
                    continue;
                };
                if kind != "M" {
                    ref_allele = "NA".to_string();
                }
                for (alt, reads) in alts.iter() {
                    let ids: Vec<&str> = reads.iter().map(|(id, _)| id.as_str()).collect();
                    if ids.len() >= 5 {
                        accepted = true;
                        lines.push(format!(
                            "{scaffold}\t{pos}\t{kind}\t{ref_allele}\t{alt}\t{}\t{cov_count}\t{}",
                            ids.len(),
                            ids.join(",")
                        ));
                    }
                }
            }
            if accepted {
                let cov_line = cov_ids.join(",");
                if !group_pos.contains_key(&cov_line) {
                    group_order.push(cov_line.clone());
                    group_pos.insert(cov_line.clone(), Vec::new());
                }
                group_pos
                    .get_mut(&cov_line)
                    .unwrap()
                    .push(format!("{scaffold}_{pos}"));
            }
        }
    }
    let mut varcov = vec!["positions\toverlap_clusters".to_string()];
    for cov_line in group_order {
        let mut positions = group_pos.remove(&cov_line).unwrap_or_default();
        positions.sort();
        varcov.push(format!("{}\t{cov_line}", positions.join(",")));
    }
    Ok((join_lines(&lines), join_lines(&varcov)))
}

fn report_text(
    settings: &CollapseSettings,
    genes: i64,
    transcripts: i64,
    accepted: i64,
    discarded: i64,
) -> String {
    let params = format!(
        "-s {} -f {} -p {} -x {} -e {} -c {} -i {} -icm {} -a {} -m {} -z {} -d {} -sj {} -sjt {} -lde {} -ses {} -b {} -log {} -rm {} -vc 5",
        settings.sam_label,
        settings.fasta_label,
        settings.prefix_label,
        settings.cap,
        settings.ends,
        py2_report_float(settings.coverage),
        py2_report_float(settings.identity),
        settings.ident_method,
        settings.five_prime,
        settings.exon_diff,
        settings.three_prime,
        settings.duplicates,
        settings.sj_priority,
        settings.sj_threshold,
        settings.lde,
        settings.ses,
        settings.bam,
        settings.log,
        settings.run_mode
    );
    format!(
        "TAMA Collapse has run successfully!\nParameters used:\t{params}\nTotal Gene Count:\t{genes}\nTotal Transcript Count:\t{transcripts}\nTotal Accepted Reads:\t{accepted}\nTotal Discarded Reads:\t{discarded}\n"
    )
}

fn py2_report_float(value: f64) -> String {
    if value.fract() == 0.0 {
        format!("{value:.1}")
    } else {
        format!("{value}")
    }
}

fn join_lines(lines: &[String]) -> String {
    if lines.is_empty() {
        return String::new();
    }
    let mut text = lines.join("\n");
    text.push('\n');
    text
}

pub fn write_collapse_texts(prefix: &Path, texts: &CollapseTexts) -> Result<(), String> {
    let files = [
        (".bed", &texts.bed),
        ("_read.txt", &texts.read_txt),
        ("_trans_report.txt", &texts.trans_report),
        ("_trans_read.bed", &texts.trans_read),
        ("_polya.txt", &texts.polya),
        ("_strand_check.txt", &texts.strand_check),
        ("_local_density_error.txt", &texts.local_density),
        ("_variants.txt", &texts.variants),
        ("_varcov.txt", &texts.varcov),
        ("_report.txt", &texts.report),
    ];
    for (suffix, body) in files {
        let path = suffix_path(prefix, suffix);
        fs::write(&path, body).map_err(|err| format!("write {}: {err}", path.display()))?;
    }
    Ok(())
}

fn suffix_path(prefix: &Path, suffix: &str) -> PathBuf {
    let mut name = prefix.as_os_str().to_os_string();
    name.push(suffix);
    PathBuf::from(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gmap_collapse_matches_oracle() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let texts = original_collapse(
            &root.join("../../tests/parity/gmap_collapse/gmap_test.sam"),
            &root.join("../../tests/parity/gmap_collapse/test_genome.fa"),
            &CollapseSettings {
                cap: "capped".to_string(),
                ends: "common_ends".to_string(),
                coverage: 99.0,
                identity: 85.0,
                ident_method: "ident_cov".to_string(),
                five_prime: 10,
                exon_diff: 10,
                three_prime: 10,
                duplicates: "merge_dup".to_string(),
                sj_priority: "no_priority".to_string(),
                sj_threshold: 10,
                lde: 1000,
                ses: "_".to_string(),
                bam: "SAM".to_string(),
                log: "log_off".to_string(),
                run_mode: "original".to_string(),
                var_support: 5,
                sam_label: "test_files/gmap_test.sam".to_string(),
                fasta_label: "test_files/test_genome.fa".to_string(),
                prefix_label: "/work/braid/tests/parity/gmap_collapse/gmap".to_string(),
            },
        )
        .unwrap_or_else(|err| panic!("{err}"));
        let golden = root.join("../../tests/parity/gmap_collapse");
        expect_file(&golden.join("gmap.bed"), &texts.bed);
        expect_file(&golden.join("gmap_read.txt"), &texts.read_txt);
        expect_file(&golden.join("gmap_trans_report.txt"), &texts.trans_report);
        expect_file(&golden.join("gmap_trans_read.bed"), &texts.trans_read);
        expect_file(&golden.join("gmap_polya.txt"), &texts.polya);
        expect_file(&golden.join("gmap_strand_check.txt"), &texts.strand_check);
        expect_file(
            &golden.join("gmap_local_density_error.txt"),
            &texts.local_density,
        );
        expect_file(&golden.join("gmap_variants.txt"), &texts.variants);
        expect_file(&golden.join("gmap_varcov.txt"), &texts.varcov);
        expect_file(&golden.join("gmap_report.txt"), &texts.report);
    }

    fn expect_file(path: &Path, actual: &str) {
        let expected = fs::read_to_string(path).unwrap_or_else(|err| panic!("{path:?}: {err}"));
        if expected == actual {
            return;
        }
        let exp_lines: Vec<&str> = expected.lines().collect();
        let got_lines: Vec<&str> = actual.lines().collect();
        let mut message = format!(
            "{} differs (expected {} lines, got {})",
            path.display(),
            exp_lines.len(),
            got_lines.len()
        );
        let limit = exp_lines.len().max(got_lines.len()).min(8);
        for i in 0..limit {
            let left = exp_lines.get(i).copied().unwrap_or("<missing>");
            let right = got_lines.get(i).copied().unwrap_or("<missing>");
            if left != right {
                message.push_str(&format!(
                    "\nline {}: expected {left}\n      got {right}",
                    i + 1
                ));
                break;
            }
        }
        if exp_lines.len() == got_lines.len() {
            for (i, (left, right)) in exp_lines.iter().zip(got_lines.iter()).enumerate() {
                if left != right {
                    message.push_str(&format!(
                        "\nline {}: expected {left}\n      got {right}",
                        i + 1
                    ));
                    break;
                }
            }
        }
        panic!("{message}");
    }
}
