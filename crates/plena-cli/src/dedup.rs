//! One primary alignment per cell barcode and UMI.
//!
//! Isosceles counts reads and does not read a UMI tag. Bambu-Clump
//! deduplicates UMIs only after a chromosome has more than 100 distinct UMIs.
//! Both commands keep the longest reference span for each barcode and UMI,
//! then the lexicographically smaller read name, then the earlier alignment,
//! so a repeated cell and UMI is one molecule at any input size. IsoQuant
//! and SCOTCH do not use this step.

use std::fs;
use std::path::Path;
use std::process::Command;

use crate::tools::{self, run};

struct Kept {
    span: i64,
    qname: String,
    order: usize,
    line: String,
}

pub fn dedup_bam(
    input: &Path,
    output_bam: &Path,
    cell_tag: &str,
    umi_tag: &str,
) -> Result<usize, String> {
    let samtools = tools::samtools()?;
    let view = Command::new(&samtools)
        .args(["view", "-h", "--"])
        .arg(input)
        .output()
        .map_err(|err| format!("samtools view {}: {err}", input.display()))?;
    if !view.status.success() {
        return Err(format!(
            "samtools view {}: {}",
            input.display(),
            String::from_utf8_lossy(&view.stderr)
        ));
    }
    let text = String::from_utf8(view.stdout)
        .map_err(|_| "samtools view returned non-UTF8".to_string())?;
    let mut header = String::new();
    let mut best: Vec<(String, Kept)> = Vec::new();
    let mut order = 0usize;
    for raw in text.lines() {
        if raw.is_empty() {
            continue;
        }
        if raw.starts_with('@') {
            header.push_str(raw);
            header.push('\n');
            continue;
        }
        let fields: Vec<&str> = raw.split('\t').collect();
        if fields.len() < 11 {
            return Err(format!(
                "SAM record {} has {} columns",
                fields[0],
                fields.len()
            ));
        }
        let flag: i32 = fields[1]
            .parse()
            .map_err(|_| format!("bad flag {} on {}", fields[1], fields[0]))?;
        if flag & 0x4 != 0 || flag & 0x100 != 0 || flag & 0x800 != 0 {
            continue;
        }
        let pos: i64 = fields[3]
            .parse()
            .map_err(|_| format!("bad POS {} on {}", fields[3], fields[0]))?;
        if pos <= 0 || fields[5] == "*" {
            continue;
        }
        let cell = tag_value(&fields, cell_tag)
            .ok_or_else(|| format!("read {} has no {cell_tag} tag", fields[0]))?;
        let umi = tag_value(&fields, umi_tag)
            .ok_or_else(|| format!("read {} has no {umi_tag} tag", fields[0]))?;
        let span = reference_span(fields[5])?;
        let key = format!("{cell}\t{umi}");
        let candidate = Kept {
            span,
            qname: fields[0].to_string(),
            order,
            line: raw.to_string(),
        };
        order += 1;
        if let Some((_, current)) = best.iter_mut().find(|(existing, _)| existing == &key) {
            let better = candidate.span > current.span
                || (candidate.span == current.span && candidate.qname < current.qname);
            if better {
                *current = candidate;
            }
        } else {
            best.push((key, candidate));
        }
    }
    if best.is_empty() {
        return Err(format!(
            "{} has no primary alignments with {cell_tag} and {umi_tag}",
            input.display()
        ));
    }
    best.sort_by(|left, right| {
        left.1
            .line
            .split('\t')
            .nth(2)
            .cmp(&right.1.line.split('\t').nth(2))
            .then(left.1.order.cmp(&right.1.order))
    });
    let parent = output_bam.parent().unwrap_or(Path::new("."));
    fs::create_dir_all(parent).map_err(|err| format!("create {}: {err}", parent.display()))?;
    let sam_path = output_bam.with_extension("dedup.sam");
    let mut sam = header;
    for (_, kept) in &best {
        sam.push_str(&kept.line);
        sam.push('\n');
    }
    fs::write(&sam_path, sam).map_err(|err| format!("write {}: {err}", sam_path.display()))?;
    let log = parent.join("dedup-sort.log");
    run(
        &samtools,
        &[
            "sort".to_string(),
            "-o".to_string(),
            output_bam.display().to_string(),
            sam_path.display().to_string(),
        ],
        &log,
    )?;
    tools::ensure_bam_index(output_bam)?;
    let _ = fs::remove_file(&sam_path);
    Ok(best.len())
}

fn tag_value(fields: &[&str], tag: &str) -> Option<String> {
    let prefix = format!("{tag}:");
    fields.iter().skip(11).find_map(|field| {
        let rest = field.strip_prefix(&prefix)?;
        let (_kind, value) = rest.split_once(':')?;
        Some(value.to_string())
    })
}

fn reference_span(cigar: &str) -> Result<i64, String> {
    let mut span = 0i64;
    let mut number = String::new();
    for character in cigar.chars() {
        if character.is_ascii_digit() {
            number.push(character);
            continue;
        }
        let len: i64 = number
            .parse()
            .map_err(|_| format!("bad CIGAR length in {cigar}"))?;
        number.clear();
        if matches!(character, 'M' | 'D' | 'N' | '=' | 'X') {
            span += len;
        } else if !matches!(character, 'I' | 'S' | 'H' | 'P') {
            return Err(format!("bad CIGAR operator {character} in {cigar}"));
        }
    }
    if !number.is_empty() || span <= 0 {
        return Err(format!("CIGAR {cigar} does not cover the reference"));
    }
    Ok(span)
}
