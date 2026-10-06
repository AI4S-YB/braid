//! Reference GTF to non-overlapping sub-exons and metagenes.
//!
//! GTF coordinates are 1-based and inclusive. Internally a feature becomes
//! the half-open interval `[start - 1, end)`, which is what SCOTCH stores
//! after converting GTF to BED.

use std::collections::BTreeMap;
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Isoform {
    pub name: String,
    /// Sub-exon indexes, in genomic order.
    pub exons: Vec<usize>,
    pub order: usize,
}

#[derive(Clone, Debug)]
pub struct Gene {
    pub name: String,
    pub id: String,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub strand: char,
    pub exons: Vec<(i64, i64)>,
    pub isoforms: Vec<Isoform>,
    /// Column name in the count matrix. Disambiguated when names collide.
    pub label: String,
}

#[derive(Clone, Debug)]
pub struct Annotation {
    pub genes: Vec<Gene>,
    /// Each metagene is a list of gene indexes. Overlapping and book-ended
    /// genes on the same chromosome share a metagene.
    pub metagenes: Vec<Vec<usize>>,
    pub gtf_text: String,
}

#[derive(Clone, Debug)]
struct Tx {
    id: String,
    order: usize,
    exons: Vec<(i64, i64)>,
}

#[derive(Clone, Debug)]
struct GeneBuild {
    name: String,
    id: String,
    chrom: String,
    start: Option<i64>,
    end: Option<i64>,
    strand: char,
    transcripts: Vec<Tx>,
}

/// Keep genes whose name, id, or label is listed. Metagenes are rebuilt.
pub fn retain_genes(annotation: &mut Annotation, names: &[String]) -> Result<(), String> {
    let mut missing = Vec::new();
    for name in names {
        let found = annotation
            .genes
            .iter()
            .any(|gene| gene.name == *name || gene.id == *name || gene.label == *name);
        if !found {
            missing.push(name.as_str());
        }
    }
    if !missing.is_empty() {
        return Err(format!(
            "gene subset not in the GTF: {}",
            missing.join(", ")
        ));
    }
    let want: std::collections::BTreeSet<&str> = names.iter().map(String::as_str).collect();
    annotation.genes.retain(|gene| {
        want.contains(gene.name.as_str())
            || want.contains(gene.id.as_str())
            || want.contains(gene.label.as_str())
    });
    if annotation.genes.is_empty() {
        return Err("gene subset matches no gene in the GTF".to_string());
    }
    annotation.metagenes = cluster_metagenes(&annotation.genes);
    Ok(())
}

pub fn load_gtf(path: &Path) -> Result<Annotation, String> {
    let text =
        std::fs::read_to_string(path).map_err(|err| format!("open {}: {err}", path.display()))?;
    parse_gtf(&text)
}

pub fn parse_gtf(text: &str) -> Result<Annotation, String> {
    let mut genes: BTreeMap<String, GeneBuild> = BTreeMap::new();
    let mut order = 0usize;
    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim_end_matches(['\r']);
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 9 {
            return Err(format!(
                "GTF line {} has {} columns",
                lineno + 1,
                cols.len()
            ));
        }
        let kind = cols[2];
        if kind != "gene" && kind != "transcript" && kind != "exon" {
            continue;
        }
        let start: i64 = cols[3]
            .parse()
            .map_err(|_| format!("GTF line {}: bad start", lineno + 1))?;
        let end: i64 = cols[4]
            .parse()
            .map_err(|_| format!("GTF line {}: bad end", lineno + 1))?;
        if start <= 0 || end < start {
            return Err(format!("GTF line {}: invalid coordinates", lineno + 1));
        }
        let strand = match cols[6] {
            "+" | "-" | "." => cols[6].chars().next().unwrap(),
            other => {
                return Err(format!("GTF line {}: bad strand {other}", lineno + 1));
            }
        };
        let attr = parse_attrs(cols[8]);
        let gene_id = attr
            .get("gene_id")
            .cloned()
            .ok_or_else(|| format!("GTF line {}: missing gene_id", lineno + 1))?;
        let gene_name = attr
            .get("gene_name")
            .cloned()
            .unwrap_or_else(|| gene_id.clone());
        let entry = genes.entry(gene_id.clone()).or_insert_with(|| GeneBuild {
            name: gene_name.clone(),
            id: gene_id.clone(),
            chrom: cols[0].to_string(),
            start: None,
            end: None,
            strand,
            transcripts: Vec::new(),
        });
        if entry.chrom != cols[0] {
            return Err(format!(
                "gene {gene_id} is on both {} and {}",
                entry.chrom, cols[0]
            ));
        }
        if kind == "gene" {
            entry.start = Some(start - 1);
            entry.end = Some(end);
            entry.strand = strand;
            entry.name = gene_name;
            continue;
        }
        let tx_id = attr
            .get("transcript_id")
            .cloned()
            .ok_or_else(|| format!("GTF line {}: missing transcript_id", lineno + 1))?;
        if entry.transcripts.iter().all(|tx| tx.id != tx_id) {
            entry.transcripts.push(Tx {
                id: tx_id.clone(),
                order,
                exons: Vec::new(),
            });
            order += 1;
        }
        if kind == "exon" {
            let tx = entry
                .transcripts
                .iter_mut()
                .find(|tx| tx.id == tx_id)
                .expect("transcript just inserted");
            tx.exons.push((start - 1, end));
        }
    }

    let mut built = Vec::new();
    for build in genes.into_values() {
        let mut exons = Vec::new();
        for tx in &build.transcripts {
            exons.extend(tx.exons.iter().copied());
        }
        if exons.is_empty() {
            continue;
        }
        exons.sort_unstable();
        let (span_start, span_end) = exons.iter().fold((exons[0].0, exons[0].1), |acc, iv| {
            (acc.0.min(iv.0), acc.1.max(iv.1))
        });
        let sub = partition_exons(&exons);
        if sub.is_empty() {
            continue;
        }
        let mut isoforms = Vec::new();
        for tx in &build.transcripts {
            if tx.exons.is_empty() {
                continue;
            }
            let mut idx = Vec::new();
            for (i, subexon) in sub.iter().enumerate() {
                if tx
                    .exons
                    .iter()
                    .any(|exon| exon.0 <= subexon.0 && exon.1 >= subexon.1)
                {
                    idx.push(i);
                }
            }
            if idx.is_empty() {
                continue;
            }
            isoforms.push(Isoform {
                name: tx.id.clone(),
                exons: idx,
                order: tx.order,
            });
        }
        isoforms.sort_by(|a, b| {
            b.exons
                .len()
                .cmp(&a.exons.len())
                .then(a.order.cmp(&b.order))
                .then(a.name.cmp(&b.name))
        });
        built.push(Gene {
            name: build.name,
            id: build.id,
            chrom: build.chrom,
            start: build.start.unwrap_or(span_start),
            end: build.end.unwrap_or(span_end),
            strand: build.strand,
            exons: sub,
            isoforms,
            label: String::new(),
        });
    }
    if built.is_empty() {
        return Err("GTF contains no gene with exons".to_string());
    }
    built.sort_by(|a, b| {
        a.chrom
            .cmp(&b.chrom)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.id.cmp(&b.id))
    });
    let mut name_count: BTreeMap<String, usize> = BTreeMap::new();
    for gene in &built {
        *name_count.entry(gene.name.clone()).or_default() += 1;
    }
    for gene in &mut built {
        gene.label = if name_count[&gene.name] > 1 {
            format!("{}_{}", gene.name, gene.id)
        } else {
            gene.name.clone()
        };
    }
    let metagenes = cluster_metagenes(&built);
    Ok(Annotation {
        genes: built,
        metagenes,
        gtf_text: text.to_string(),
    })
}

/// Disjoint half-open pieces of the union of `exons`.
pub fn partition_exons(exons: &[(i64, i64)]) -> Vec<(i64, i64)> {
    let mut points = Vec::new();
    for &(start, end) in exons {
        if end > start {
            points.push(start);
            points.push(end);
        }
    }
    points.sort_unstable();
    points.dedup();
    let mut out = Vec::new();
    for pair in points.windows(2) {
        let (start, end) = (pair[0], pair[1]);
        if exons.iter().any(|&(a, b)| a <= start && b >= end) {
            out.push((start, end));
        }
    }
    out
}

fn cluster_metagenes(genes: &[Gene]) -> Vec<Vec<usize>> {
    let mut groups = Vec::new();
    let mut current: Vec<usize> = Vec::new();
    let mut chrom = String::new();
    let mut end = 0i64;
    for (index, gene) in genes.iter().enumerate() {
        let fresh = current.is_empty() || gene.chrom != chrom || gene.start > end;
        if fresh {
            if !current.is_empty() {
                groups.push(std::mem::take(&mut current));
            }
            chrom = gene.chrom.clone();
            end = gene.end;
        } else {
            end = end.max(gene.end);
        }
        current.push(index);
    }
    if !current.is_empty() {
        groups.push(current);
    }
    groups
}

fn parse_attrs(text: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for part in text.split(';') {
        let part = part.trim();
        if part.is_empty() {
            continue;
        }
        let mut pieces = part.splitn(2, char::is_whitespace);
        let Some(key) = pieces.next() else { continue };
        let Some(value) = pieces.next() else { continue };
        let value = value.trim().trim_matches('"').to_string();
        if matches!(key, "gene_id" | "gene_name" | "transcript_id" | "gene_type") {
            out.insert(key.to_string(), value);
        }
    }
    out
}

pub fn subexon_indexes(bits: &[bool]) -> Vec<usize> {
    bits.iter()
        .enumerate()
        .filter_map(|(i, on)| on.then_some(i))
        .collect()
}

/// Decimal of `sum(2^i for bits[i])`, the SCOTCH novel-isoform id.
pub fn bitset_decimal(bits: &[bool]) -> String {
    let mut digits = vec![0u8];
    for &bit in bits.iter().rev() {
        let mut carry = 0u8;
        for digit in &mut digits {
            let value = *digit * 2 + carry;
            *digit = value % 10;
            carry = value / 10;
        }
        if carry > 0 {
            digits.push(carry);
        }
        if bit {
            let mut add = 1u8;
            for digit in &mut digits {
                let value = *digit + add;
                *digit = value % 10;
                add = value / 10;
                if add == 0 {
                    break;
                }
            }
            if add > 0 {
                digits.push(add);
            }
        }
    }
    let mut text: String = digits
        .into_iter()
        .rev()
        .map(|d| char::from(b'0' + d))
        .collect();
    while text.len() > 1 && text.starts_with('0') {
        text.remove(0);
    }
    text
}

pub fn novel_name(bits: &[bool]) -> String {
    format!("novelIsoform_{}", bitset_decimal(bits))
}

/// Merge sub-exons that touch and convert them to 1-based inclusive GTF bounds.
pub fn gtf_exons(gene: &Gene, indexes: &[usize]) -> Vec<(i64, i64)> {
    let mut pieces: Vec<(i64, i64)> = indexes
        .iter()
        .filter_map(|&index| gene.exons.get(index).copied())
        .collect();
    pieces.sort_unstable();
    if pieces.is_empty() {
        return Vec::new();
    }
    let mut merged = Vec::new();
    let (mut start, mut end) = pieces[0];
    for &(next, next_end) in &pieces[1..] {
        if next == end {
            end = end.max(next_end);
        } else {
            merged.push((start + 1, end));
            start = next;
            end = next_end;
        }
    }
    merged.push((start + 1, end));
    merged
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overlapping_exons_become_subexons_and_bookends_share_a_metagene() {
        let text = "\
chr1\t.\tgene\t1\t150\t.\t+\t.\tgene_id \"G\"; gene_name \"G\";
chr1\t.\ttranscript\t1\t100\t.\t+\t.\tgene_id \"G\"; gene_name \"G\"; transcript_id \"A\";
chr1\t.\texon\t1\t100\t.\t+\t.\tgene_id \"G\"; gene_name \"G\"; transcript_id \"A\";
chr1\t.\ttranscript\t51\t150\t.\t+\t.\tgene_id \"G\"; gene_name \"G\"; transcript_id \"B\";
chr1\t.\texon\t51\t150\t.\t+\t.\tgene_id \"G\"; gene_name \"G\"; transcript_id \"B\";
chr1\t.\tgene\t151\t160\t.\t+\t.\tgene_id \"H\"; gene_name \"H\";
chr1\t.\ttranscript\t151\t160\t.\t+\t.\tgene_id \"H\"; gene_name \"H\"; transcript_id \"C\";
chr1\t.\texon\t151\t160\t.\t+\t.\tgene_id \"H\"; gene_name \"H\"; transcript_id \"C\";
chr1\t.\tgene\t200\t210\t.\t+\t.\tgene_id \"I\"; gene_name \"I\";
chr1\t.\ttranscript\t200\t210\t.\t+\t.\tgene_id \"I\"; gene_name \"I\"; transcript_id \"D\";
chr1\t.\texon\t200\t210\t.\t+\t.\tgene_id \"I\"; gene_name \"I\"; transcript_id \"D\";
";
        let ann = parse_gtf(text).unwrap();
        assert_eq!(ann.genes[0].exons, vec![(0, 50), (50, 100), (100, 150)]);
        assert_eq!(ann.genes[0].isoforms[0].exons, vec![0, 1]);
        assert_eq!(ann.genes[0].isoforms[1].exons, vec![1, 2]);
        assert_eq!(ann.metagenes.len(), 2);
        assert_eq!(ann.metagenes[0].len(), 2);
        assert_eq!(bitset_decimal(&[true, true, false]), "3");
        assert_eq!(novel_name(&[true, false, true]), "novelIsoform_5");
    }
}
