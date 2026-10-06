//! FLAIR `collapse_isoforms_precise.py` from BrooksLabUCSC/flair commit
//! `573414c551332bf6348a9d04ba7cf562f67416cb`.
//!
//! `flair collapse` invokes that script with isoform-level ends (`-i`), the
//! script's default support proportion `0.25`, and `--nosplice chrM`. This
//! module matches that call. It does not reassign reads with minimap2; that
//! later filter is a separate FLAIR step and is not part of this collapse.

use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs;
use std::path::Path;
use std::process::{Command, Stdio};

type ChromSites = HashMap<String, HashSet<i64>>;

struct SiteSearch<'a> {
    settings: &'a PreciseSettings,
    annot_tss: &'a ChromSites,
    annot_tes: &'a ChromSites,
    chrom: Option<&'a str>,
    junc: (f64, f64),
    max_results: usize,
}

struct Junctions {
    pairs: Vec<(i64, i64)>,
    span: (i64, i64),
}

#[derive(Clone, Debug)]
pub struct PreciseSettings {
    pub window: i64,
    pub max_results: usize,
    /// Proportion when `< 1`, otherwise a read count. The script default is `0.25`.
    pub minsupport: f64,
    /// `-i` / isoform-level TSS. `flair collapse` sets this unless `--gene_tss`.
    pub isoform_tss: bool,
    /// `none`, `longest`, or `best_only`.
    pub no_redundant: String,
    pub nosplice: Vec<String>,
    pub clean: bool,
}

impl Default for PreciseSettings {
    fn default() -> Self {
        Self {
            window: 100,
            max_results: 2,
            minsupport: 0.25,
            isoform_tss: true,
            no_redundant: "none".to_string(),
            nosplice: vec!["chrM".to_string()],
            clean: false,
        }
    }
}

/// Collapse a BED12 text the way `collapse_isoforms_precise.py` does.
/// `gtf` is the optional annotation used to snap TSS/TES.
pub fn collapse_precise_bed(
    bed: &str,
    gtf: Option<&str>,
    settings: &PreciseSettings,
) -> Result<String, String> {
    let (annot_tss, annot_tes) = match gtf {
        Some(text) => load_annot_ends(text)?,
        None => (HashMap::new(), HashMap::new()),
    };
    let mut isoforms: OrderMap<String, OrderMap<String, Chain>> = OrderMap::new();
    let mut single: OrderMap<String, BTreeMap<(i64, i64), SeRead>> = OrderMap::new();
    let nosplice: HashSet<&str> = settings.nosplice.iter().map(String::as_str).collect();

    for (lineno, raw) in bed.lines().enumerate() {
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let fields: Vec<String> = raw.split('\t').map(str::to_string).collect();
        if fields.len() < 12 {
            return Err(format!(
                "BED line {} has {} columns",
                lineno + 1,
                fields.len()
            ));
        }
        let chrom = fields[0].clone();
        let tss: i64 = fields[1]
            .parse()
            .map_err(|_| format!("bad BED start on line {}", lineno + 1))?;
        let tes: i64 = fields[2]
            .parse()
            .map_err(|_| format!("bad BED end on line {}", lineno + 1))?;
        let Junctions {
            pairs: junctions,
            span: junc,
        } = junctions_bed12(&fields)?;
        if junctions.is_empty() {
            let chrom_se = single.get_or_insert_with(chrom, BTreeMap::new);
            let entry = chrom_se.entry((tss, tes)).or_insert_with(|| SeRead {
                count: 0,
                line: fields,
            });
            entry.count += 1;
            continue;
        }
        if nosplice.contains(chrom.as_str()) {
            continue;
        }
        let key = junction_key(&junctions);
        let chrom_iso = isoforms.get_or_insert_with(chrom, OrderMap::new);
        let chain = chrom_iso.get_or_insert_with(key, || Chain {
            tss: BTreeMap::new(),
            tss_tes: BTreeMap::new(),
            line: fields,
            junc,
        });
        *chain.tss.entry(tss).or_insert(0) += 1;
        *chain
            .tss_tes
            .entry(tss)
            .or_default()
            .entry(tes)
            .or_insert(0) += 1;
    }

    let mut se_chroms: Vec<(String, usize)> = single
        .keys()
        .map(|chrom| {
            let n = single.get(chrom).map(|m| m.len()).unwrap_or(0);
            (chrom.clone(), n)
        })
        .collect();
    se_chroms.sort_by_key(|chrom| Reverse(chrom.1));

    let mut out_lines: Vec<Vec<String>> = Vec::new();
    for (chrom, _) in &se_chroms {
        let groups = group_single_exon(single.get(chrom).unwrap(), settings.window);
        out_lines.extend(collapse_single_exon(
            &groups, settings, &annot_tss, &annot_tes,
        ));
    }
    for chrom in isoforms.keys() {
        let chains = isoforms.get(chrom).unwrap();
        for chain in chains.values() {
            out_lines.extend(collapse_chain(chain, settings, &annot_tss, &annot_tes));
        }
    }
    if out_lines.is_empty() {
        return Ok(String::new());
    }
    Ok(out_lines
        .into_iter()
        .map(|fields| fields.join("\t"))
        .collect::<Vec<_>>()
        .join("\n")
        + "\n")
}

pub fn collapse_precise_files(
    query: &Path,
    gtf: Option<&Path>,
    settings: &PreciseSettings,
) -> Result<String, String> {
    let bed = fs::read_to_string(query).map_err(|e| format!("read {}: {e}", query.display()))?;
    let gtf_text = match gtf {
        Some(path) => {
            Some(fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?)
        }
        None => None,
    };
    collapse_precise_bed(&bed, gtf_text.as_deref(), settings)
}

/// Primary alignments in a coordinate-sorted SAM or BAM, as BED12.
/// BAM is decoded with `samtools view`.
pub fn alignments_to_bed12(path: &Path, bam: bool) -> Result<String, String> {
    let text = if bam {
        let output = Command::new("samtools")
            .args(["view", "--"])
            .arg(path)
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .output()
            .map_err(|e| format!("BAM input requires samtools on PATH: {e}"))?;
        if !output.status.success() {
            return Err(format!(
                "samtools view {}: {}",
                path.display(),
                String::from_utf8_lossy(&output.stderr)
            ));
        }
        String::from_utf8(output.stdout)
            .map_err(|_| "samtools view returned non-UTF8".to_string())?
    } else {
        fs::read_to_string(path).map_err(|e| format!("read {}: {e}", path.display()))?
    };
    let mut lines = Vec::new();
    for (lineno, raw) in text.lines().enumerate() {
        if raw.is_empty() || raw.starts_with('@') {
            continue;
        }
        if let Some(bed) =
            sam_to_bed12(raw).map_err(|e| format!("{} line {}: {e}", path.display(), lineno + 1))?
        {
            lines.push(bed);
        }
    }
    Ok(lines.join("\n") + if lines.is_empty() { "" } else { "\n" })
}

/// BED12 rows to a GTF. Names that contain `_` are split the way FLAIR's
/// `bed_to_gtf.py` splits them. Other names are used as both ids so a read
/// name still becomes a transcript model.
pub fn bed12_to_gtf(bed: &str, source: &str) -> Result<String, String> {
    let mut out = String::new();
    for (lineno, raw) in bed.lines().enumerate() {
        if raw.is_empty() {
            continue;
        }
        let line: Vec<&str> = raw.split('\t').collect();
        if line.len() < 12 {
            return Err(format!(
                "BED line {} has {} columns",
                lineno + 1,
                line.len()
            ));
        }
        let chrom = line[0];
        let start: i64 = line[1]
            .parse()
            .map_err(|_| format!("bad start on BED line {}", lineno + 1))?;
        let strand = line[5];
        let name = line[3].replace(';', ":");
        let (transcript_id, gene_id) = split_iso_gene(&name);
        let rel: Vec<i64> = line[11]
            .trim_end_matches(',')
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| {
                s.parse::<i64>()
                    .map_err(|_| format!("bad blockStarts on BED line {}", lineno + 1))
            })
            .collect::<Result<_, _>>()?;
        let sizes: Vec<i64> = line[10]
            .trim_end_matches(',')
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| {
                s.parse::<i64>()
                    .map_err(|_| format!("bad blockSizes on BED line {}", lineno + 1))
            })
            .collect::<Result<_, _>>()?;
        if rel.len() != sizes.len() || rel.is_empty() {
            return Err(format!("BED line {} has inconsistent blocks", lineno + 1));
        }
        let abs: Vec<i64> = rel.iter().map(|r| start + r).collect();
        let transcript_end = abs.last().copied().unwrap() + sizes.last().copied().unwrap();
        let attrs = format!("gene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";");
        out.push_str(&format!(
            "{chrom}\t{source}\ttranscript\t{}\t{transcript_end}\t.\t{strand}\t.\t{attrs}\n",
            start + 1
        ));
        for (b, (tstart, size)) in abs.iter().zip(sizes.iter()).enumerate() {
            out.push_str(&format!(
                "{chrom}\t{source}\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"{gene_id}\"; transcript_id \"{transcript_id}\"; exon_number \"{b}\";\n",
                tstart + 1,
                tstart + size
            ));
        }
    }
    Ok(out)
}

/// `bed_to_gtf.py` splits on the last accession marker, and that marker's
/// underscore stays out of both ids (`iso_chr1` → transcript `iso`, gene `chr1`).
/// Any other underscore is a plain last-underscore split. Names without one
/// are used for both ids.
fn split_iso_gene(name: &str) -> (String, String) {
    for marker in ["_chr", "_XM", "_XR", "_NM", "_NR", "_R2_"] {
        if let Some(i) = name.rfind(marker) {
            return (name[..i].to_string(), name[i + 1..].to_string());
        }
    }
    if let Some(i) = name.rfind('_') {
        (name[..i].to_string(), name[i + 1..].to_string())
    } else {
        (name.to_string(), name.to_string())
    }
}

fn load_annot_ends(gtf: &str) -> Result<(ChromSites, ChromSites), String> {
    let mut tss: ChromSites = HashMap::new();
    let mut tes: ChromSites = HashMap::new();
    for (lineno, raw) in gtf.lines().enumerate() {
        if raw.is_empty() || raw.starts_with('#') {
            continue;
        }
        let line: Vec<&str> = raw.split('\t').collect();
        if line.len() < 5 {
            return Err(format!(
                "GTF line {} has {} columns",
                lineno + 1,
                line.len()
            ));
        }
        if line[2] != "transcript" {
            continue;
        }
        let start: i64 = line[3]
            .parse::<i64>()
            .map_err(|_| format!("bad GTF start on line {}", lineno + 1))?
            - 1;
        let end: i64 = line[4]
            .parse()
            .map_err(|_| format!("bad GTF end on line {}", lineno + 1))?;
        tss.entry(line[0].to_string()).or_default().insert(start);
        tes.entry(line[0].to_string()).or_default().insert(end);
    }
    Ok((tss, tes))
}

struct SeRead {
    count: i64,
    line: Vec<String>,
}

struct Chain {
    tss: BTreeMap<i64, i64>,
    tss_tes: BTreeMap<i64, BTreeMap<i64, i64>>,
    line: Vec<String>,
    junc: (i64, i64),
}

struct SeGroup {
    tss: BTreeMap<i64, i64>,
    tss_tes: BTreeMap<i64, BTreeMap<i64, i64>>,
    bounds: [f64; 3],
    line: Vec<String>,
    strand: OrderMap<String, i64>,
}

fn group_single_exon(reads: &BTreeMap<(i64, i64), SeRead>, window: i64) -> Vec<SeGroup> {
    let mut groups: Vec<SeGroup> = Vec::new();
    for (&(tss, tes), read) in reads {
        let mut overlapped: Vec<(usize, f64)> = Vec::new();
        let start = groups.len().saturating_sub(7);
        for g in (start..groups.len()).rev() {
            let (is_overlap, coverage) = overlap(
                (tss as f64, tes as f64),
                (groups[g].bounds[0], groups[g].bounds[1]),
            );
            if is_overlap {
                overlapped.push((g, coverage));
            } else {
                break;
            }
        }
        overlapped.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));
        let mut added = false;
        for (g, _) in overlapped {
            if add_se(&mut groups[g], tss, tes, read, window) {
                added = true;
                break;
            }
        }
        if !added {
            let mut group = SeGroup {
                tss: BTreeMap::new(),
                tss_tes: BTreeMap::new(),
                bounds: [tss as f64, tes as f64, read.count as f64],
                line: read.line.clone(),
                strand: OrderMap::new(),
            };
            add_se_new(&mut group, tss, tes, read);
            groups.push(group);
        }
    }
    groups
}

fn add_se(group: &mut SeGroup, tss: i64, tes: i64, read: &SeRead, window: i64) -> bool {
    if tss as f64 > group.bounds[1] || tes as f64 > group.bounds[1] + window as f64 {
        return false;
    }
    let support = read.count as f64;
    group.bounds[2] += support;
    group.bounds[0] += support * (tss as f64 - group.bounds[0]) / group.bounds[2];
    group.bounds[1] += support * (tes as f64 - group.bounds[1]) / group.bounds[2];
    add_se_new(group, tss, tes, read);
    true
}

fn add_se_new(group: &mut SeGroup, tss: i64, tes: i64, read: &SeRead) {
    let strand = read.line.get(5).cloned().unwrap_or_else(|| ".".to_string());
    *group.strand.get_or_insert_with(strand, || 0) += 1;
    *group.tss.entry(tss).or_insert(0) += read.count;
    *group
        .tss_tes
        .entry(tss)
        .or_default()
        .entry(tes)
        .or_insert(0) += read.count;
}

fn overlap(mut a: (f64, f64), mut b: (f64, f64)) -> (bool, f64) {
    if b.0 < a.0 {
        std::mem::swap(&mut a, &mut b);
    }
    let is_overlap = a.0 <= b.0 && b.0 <= a.1;
    if !is_overlap || a.1 == a.0 {
        return (is_overlap, 0.0);
    }
    let coverage = (a.1.min(b.1) - b.0) / (a.1 - a.0);
    (true, coverage)
}

fn collapse_single_exon(
    groups: &[SeGroup],
    settings: &PreciseSettings,
    annot_tss: &ChromSites,
    annot_tes: &ChromSites,
) -> Vec<Vec<String>> {
    let mut towrite: Vec<Vec<String>> = Vec::new();
    let mut used: HashMap<(i64, i64), usize> = HashMap::new();
    let mut names: HashMap<String, usize> = HashMap::new();
    for group in groups {
        let chrom = group.line[0].clone();
        let ends = find_best_sites(
            &group.tss,
            &group.tss_tes,
            &SiteSearch {
                settings,
                annot_tss,
                annot_tes,
                chrom: Some(chrom.as_str()),
                junc: (group.bounds[0], group.bounds[1]),
                max_results: 1,
            },
        );
        let name = group.line[3].clone();
        let strand_count = group.strand.values().copied().max().unwrap_or(0);
        for (tss, tes, support, _, _) in ends {
            if let Some(index) = used.get(&(tss, tes)) {
                if !settings.clean {
                    let last = towrite[*index].last_mut().unwrap();
                    let value: i64 = last.parse().unwrap_or(0);
                    *last = (value + support).to_string();
                }
                continue;
            }
            used.insert((tss, tes), towrite.len());
            let mut edited = edit_line_bed12(&group.line, tss, tes, Some(tes - tss));
            edited[8] = strand_count.to_string();
            if !settings.clean {
                edited.push(support.to_string());
            }
            let seen = names.entry(name.clone()).or_insert(0);
            if *seen >= 1 {
                edited[3] = rename_isoform(&name, *seen);
            }
            *seen += 1;
            towrite.push(edited);
        }
    }
    if settings.isoform_tss {
        return towrite;
    }
    homogenize_single_exon(towrite, settings)
}

fn homogenize_single_exon(
    towrite: Vec<Vec<String>>,
    settings: &PreciseSettings,
) -> Vec<Vec<String>> {
    let mut all_starts: BTreeMap<i64, i64> = BTreeMap::new();
    let mut all_ends: BTreeMap<i64, i64> = BTreeMap::new();
    for line in &towrite {
        let tss: i64 = line[1].parse().unwrap_or(0);
        let tes: i64 = line[2].parse().unwrap_or(0);
        let support: i64 = line.last().and_then(|s| s.parse().ok()).unwrap_or(0);
        *all_starts.entry(tss).or_insert(0) += support;
        *all_ends.entry(tes).or_insert(0) += support;
    }
    let mut rewritten: Vec<Vec<String>> = Vec::new();
    let mut used: HashMap<(i64, i64), usize> = HashMap::new();
    for line in towrite {
        let mut tss: i64 = line[1].parse().unwrap_or(0);
        let mut tes: i64 = line[2].parse().unwrap_or(0);
        let mut tss_support = all_starts.get(&tss).copied().unwrap_or(0);
        for t in ((tss - settings.window)..(tss + settings.window)).rev() {
            if let Some(other) = all_starts.get(&t).copied() {
                if other > tss_support {
                    tss = t;
                    tss_support = other;
                }
            }
        }
        let mut tes_support = all_ends.get(&tes).copied().unwrap_or(0);
        for t in (tes - settings.window)..(tes + settings.window) {
            if let Some(other) = all_ends.get(&t).copied() {
                if other > tes_support {
                    tes = t;
                    tes_support = other;
                }
            }
        }
        if let Some(index) = used.get(&(tss, tes)) {
            if !settings.clean {
                let last = rewritten[*index].last_mut().unwrap();
                let value: i64 = last.parse().unwrap_or(0);
                let extra: i64 = line.last().and_then(|s| s.parse().ok()).unwrap_or(0);
                *last = (value + extra).to_string();
            }
            continue;
        }
        used.insert((tss, tes), rewritten.len());
        rewritten.push(edit_line_bed12(&line, tss, tes, Some(tes - tss)));
    }
    rewritten
}

fn collapse_chain(
    chain: &Chain,
    settings: &PreciseSettings,
    annot_tss: &ChromSites,
    annot_tes: &ChromSites,
) -> Vec<Vec<String>> {
    let chrom = chain.line[0].as_str();
    let ends = find_best_sites(
        &chain.tss,
        &chain.tss_tes,
        &SiteSearch {
            settings,
            annot_tss,
            annot_tes,
            chrom: Some(chrom),
            junc: (chain.junc.0 as f64, chain.junc.1 as f64),
            max_results: settings.max_results,
        },
    );
    if ends.is_empty() {
        return Vec::new();
    }
    if settings.isoform_tss && settings.no_redundant == "longest" {
        let tss = ends.iter().map(|e| e.0).min().unwrap();
        let tes = ends.iter().map(|e| e.1).max().unwrap();
        let mut line = edit_line_bed12(&chain.line, tss, tes, None);
        if !settings.clean {
            let left = ends.iter().find(|e| e.0 == tss).map(|e| e.3).unwrap_or(0);
            let right = ends.iter().find(|e| e.1 == tes).map(|e| e.4).unwrap_or(0);
            line.push(format_support_average(left, right));
        }
        return vec![line];
    }
    if settings.isoform_tss && settings.no_redundant == "best_only" {
        let mut line = edit_line_bed12(&chain.line, ends[0].0, ends[0].1, None);
        if !settings.clean {
            line.push(ends[0].3.to_string());
        }
        return vec![line];
    }
    let name = chain.line.get(9).cloned().unwrap_or_default();
    let mut rows = Vec::new();
    for (i, (tss, tes, support, _, _)) in ends.into_iter().enumerate() {
        let mut edited = edit_line_bed12(&chain.line, tss, tes, None);
        if !settings.clean {
            edited.push(support.to_string());
        }
        if i >= 1 {
            edited[3] = rename_isoform(&name, i);
        }
        rows.push(edited);
        if settings.no_redundant == "longest" {
            break;
        }
    }
    rows
}

fn format_support_average(left: i64, right: i64) -> String {
    let value = (left as f64 + right as f64) / 2.0;
    if value.fract() == 0.0 {
        format!("{value:.1}")
    } else {
        // Python `str` of a binary float. The default path does not use this.
        format!("{value}")
    }
}

fn rename_isoform(name: &str, index: usize) -> String {
    if let Some(pos) = name.rfind('_') {
        format!("{}-{index}{}", &name[..pos], &name[pos..])
    } else {
        format!("{name}-{index}")
    }
}

/// `(tss, tes, support, tss_count, tes_count)`.
fn find_best_sites(
    sites_tss: &BTreeMap<i64, i64>,
    sites_tes: &BTreeMap<i64, BTreeMap<i64, i64>>,
    search: &SiteSearch<'_>,
) -> Vec<(i64, i64, i64, i64, i64)> {
    let total: i64 = sites_tss.values().sum();
    if total == 0 {
        return Vec::new();
    }
    let found_tss = find_tsss(sites_tss.clone(), total, true, search);
    if found_tss.is_empty() {
        return Vec::new();
    }
    let mut ends = Vec::new();
    for tss in found_tss {
        let mut specific: BTreeMap<i64, i64> = BTreeMap::new();
        for (&tss2, tes_map) in sites_tes {
            if (tss2 - tss.pos).abs() <= search.settings.window {
                for (&tes, &count) in tes_map {
                    *specific.entry(tes).or_insert(0) += count;
                }
            }
        }
        let found_tes = find_tsss(specific, total, false, search);
        for tes in found_tes {
            ends.push((
                tss.pos,
                tes.pos,
                tss.count.max(tes.count),
                tss.count,
                tes.count,
            ));
        }
    }
    ends
}

struct Site {
    pos: i64,
    count: i64,
}

fn find_tsss(
    mut sites: BTreeMap<i64, i64>,
    total: i64,
    finding_tss: bool,
    search: &SiteSearch<'_>,
) -> Vec<Site> {
    let settings = search.settings;
    let mut remaining = sites.values().sum::<i64>() as f64;
    let total_f = total as f64;
    let n_sites = sites.len().max(1) as f64;
    let avg = remaining / n_sites;
    let mut found = Vec::new();
    while ((settings.minsupport < 1.0 && remaining / total_f > settings.minsupport)
        || remaining >= settings.minsupport)
        && found.len() < search.max_results
    {
        let (next, best_pos, best_count) = find_best_tss(&sites, finding_tss, settings.window);
        let new_remaining: i64 = next.values().sum();
        let used = remaining - new_remaining as f64;
        remaining = new_remaining as f64;
        if !found.is_empty()
            && (settings.no_redundant == "best_only"
                || (settings.minsupport < 1.0
                    && (used / total_f < settings.minsupport || best_count < 4 || used < avg)))
        {
            break;
        }
        let mut closest = None;
        let annot = if finding_tss {
            search.annot_tss
        } else {
            search.annot_tes
        };
        if let Some(chrom) = search.chrom {
            if let Some(set) = annot.get(chrom) {
                let mut best_dist = i64::MAX;
                for t in (best_pos - settings.window)..(best_pos + settings.window) {
                    if !set.contains(&t) {
                        continue;
                    }
                    let on_side = if finding_tss {
                        (t as f64) < search.junc.0
                    } else {
                        (t as f64) > search.junc.1
                    };
                    if !on_side {
                        continue;
                    }
                    let dist = (t - best_pos).abs();
                    if dist < best_dist {
                        best_dist = dist;
                        closest = Some(t);
                    }
                }
            }
        }
        found.push(Site {
            pos: closest.unwrap_or(best_pos),
            count: best_count,
        });
        sites = next;
        if sites.is_empty() {
            break;
        }
    }
    found
}

fn find_best_tss(
    sites: &BTreeMap<i64, i64>,
    finding_tss: bool,
    window: i64,
) -> (BTreeMap<i64, i64>, i64, i64) {
    let order: Vec<i64> = if finding_tss {
        sites.keys().copied().collect()
    } else {
        sites.keys().copied().rev().collect()
    };
    let mut best_pos = 0;
    let mut best_weight = 0.0;
    let mut best_count = 0;
    let mut have = false;
    for s in order {
        let mut weight = 0.0;
        for (&s2, &count) in sites {
            if s == s2 {
                weight += count as f64;
            } else if (s - s2).abs() < window {
                weight += ((window - (s - s2).abs()) as f64 / window as f64) * count as f64;
            }
        }
        if !have || weight > best_weight {
            best_pos = s;
            best_weight = weight;
            best_count = sites[&s];
            have = true;
        }
    }
    let mut left = sites.clone();
    left.retain(|s, _| (s - best_pos).abs() > window);
    (left, best_pos, best_count)
}

fn edit_line_bed12(line: &[String], tss: i64, tes: i64, blocksize: Option<i64>) -> Vec<String> {
    let mut line = line.to_vec();
    if let Some(size) = blocksize {
        line[10] = format!("{size},");
        line[11] = "0,".to_string();
        line[1] = tss.to_string();
        line[6] = tss.to_string();
        line[2] = tes.to_string();
        line[7] = tes.to_string();
        return line;
    }
    let mut sizes: Vec<i64> = split_csv_nums(&line[10]);
    let mut starts: Vec<i64> = split_csv_nums(&line[11]);
    if sizes.is_empty() || starts.is_empty() {
        return line;
    }
    let tstart: i64 = line[1].parse().unwrap_or(0);
    let tend: i64 = line[2].parse().unwrap_or(0);
    let n = sizes.len().min(starts.len());
    sizes.truncate(n);
    starts.truncate(n);
    sizes[0] += tstart - tss;
    let last = sizes.len() - 1;
    sizes[last] += tes - tend;
    for start in &mut starts {
        *start += tstart;
    }
    starts[0] = tss;
    line[1] = tss.to_string();
    line[6] = tss.to_string();
    line[2] = tes.to_string();
    line[7] = tes.to_string();
    line[10] = format!(
        "{},",
        sizes
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );
    line[11] = format!(
        "{},",
        starts
            .iter()
            .map(|v| (v - tss).to_string())
            .collect::<Vec<_>>()
            .join(",")
    );
    line
}

fn split_csv_nums(text: &str) -> Vec<i64> {
    let mut parts: Vec<&str> = text.split(',').collect();
    if parts.last().is_some_and(|s| s.is_empty()) {
        parts.pop();
    }
    parts.iter().filter_map(|s| s.parse().ok()).collect()
}

fn junctions_bed12(line: &[String]) -> Result<Junctions, String> {
    let chrstart: i64 = line[1]
        .parse()
        .map_err(|_| format!("bad start {}", line[1]))?;
    let mut rel = split_csv_nums(&line[11]);
    let mut sizes = split_csv_nums(&line[10]);
    // `split(',')[:-1]` drops the empty field after a trailing comma.
    if rel.len() != sizes.len() {
        let n = rel.len().min(sizes.len());
        rel.truncate(n);
        sizes.truncate(n);
    }
    if rel.len() <= 1 {
        return Ok(Junctions {
            pairs: Vec::new(),
            span: (0, 0),
        });
    }
    let starts: Vec<i64> = rel.iter().map(|n| n + chrstart + 1).collect();
    let sizes: Vec<i64> = sizes.iter().map(|n| n - 1).collect();
    let mut pairs = Vec::new();
    for b in 0..starts.len() - 1 {
        pairs.push((starts[b] + sizes[b], starts[b + 1]));
    }
    let span = (starts[0] + sizes[0], *starts.last().unwrap());
    Ok(Junctions { pairs, span })
}

fn junction_key(junctions: &[(i64, i64)]) -> String {
    let mut junctions = junctions.to_vec();
    junctions.sort();
    let body = junctions
        .iter()
        .map(|(a, b)| format!("({a}, {b})"))
        .collect::<Vec<_>>()
        .join(", ");
    format!("[{body}]")
}

fn sam_to_bed12(line: &str) -> Result<Option<String>, String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return Err(format!("SAM line has {} columns", fields.len()));
    }
    let flag: i32 = fields[1]
        .parse()
        .map_err(|_| format!("bad flag {}", fields[1]))?;
    if flag & 0x4 != 0 || flag & 0x100 != 0 || flag & 0x800 != 0 {
        return Ok(None);
    }
    let pos: i64 = fields[3]
        .parse()
        .map_err(|_| format!("bad POS {}", fields[3]))?;
    if pos <= 0 || fields[5] == "*" {
        return Ok(None);
    }
    let mapq = fields[4];
    let strand = if flag & 0x10 != 0 { "-" } else { "+" };
    let exons = exons_from_cigar(pos, fields[5])?;
    if exons.is_empty() {
        return Ok(None);
    }
    let bed_start = exons[0].0;
    let bed_end = exons.last().unwrap().1;
    let blocks = exons.len();
    let sizes = exons
        .iter()
        .map(|(s, e)| (e - s).to_string())
        .collect::<Vec<_>>()
        .join(",");
    let rels = exons
        .iter()
        .map(|(s, _)| (s - bed_start).to_string())
        .collect::<Vec<_>>()
        .join(",");
    Ok(Some(format!(
        "{chrom}\t{bed_start}\t{bed_end}\t{name}\t{mapq}\t{strand}\t{bed_start}\t{bed_end}\t0\t{blocks}\t{sizes},\t{rels},",
        chrom = fields[2],
        name = fields[0],
    )))
}

fn exons_from_cigar(pos: i64, cigar: &str) -> Result<Vec<(i64, i64)>, String> {
    let mut exons = Vec::new();
    let mut cursor = pos - 1;
    let mut exon_start = cursor;
    let mut number = String::new();
    let mut saw = false;
    for c in cigar.chars() {
        if c.is_ascii_digit() {
            number.push(c);
            continue;
        }
        let len: i64 = number
            .parse()
            .map_err(|_| format!("bad CIGAR length in {cigar}"))?;
        number.clear();
        match c {
            'M' | 'D' | '=' | 'X' => {
                if !saw {
                    exon_start = cursor;
                    saw = true;
                }
                cursor += len;
            }
            'N' => {
                if saw {
                    exons.push((exon_start, cursor));
                    saw = false;
                }
                cursor += len;
            }
            'I' | 'S' | 'H' | 'P' => {}
            other => return Err(format!("bad CIGAR operator {other} in {cigar}")),
        }
    }
    if !number.is_empty() {
        return Err(format!("CIGAR ends in a length: {cigar}"));
    }
    if saw {
        exons.push((exon_start, cursor));
    }
    Ok(exons)
}

/// Insertion-ordered map. Iteration follows first insertion, matching CPython 3.7+.
struct OrderMap<K, V> {
    keys: Vec<K>,
    values: Vec<V>,
}

impl<K: PartialEq + Clone, V> OrderMap<K, V> {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            values: Vec::new(),
        }
    }

    fn keys(&self) -> impl Iterator<Item = &K> {
        self.keys.iter()
    }

    fn values(&self) -> impl Iterator<Item = &V> {
        self.values.iter()
    }

    fn get(&self, key: &K) -> Option<&V> {
        self.keys
            .iter()
            .position(|k| k == key)
            .map(|i| &self.values[i])
    }

    fn get_or_insert_with(&mut self, key: K, make: impl FnOnce() -> V) -> &mut V {
        if let Some(i) = self.keys.iter().position(|k| k == &key) {
            return &mut self.values[i];
        }
        self.keys.push(key);
        self.values.push(make());
        let last = self.values.len() - 1;
        &mut self.values[last]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/parity/flair_collapse")
    }

    #[test]
    fn precise_collapse_matches_upstream_bed() {
        let dir = root();
        let bed = fs::read_to_string(dir.join("in.bed")).unwrap();
        let settings = PreciseSettings::default();
        let got = collapse_precise_bed(&bed, None, &settings).unwrap();
        let gold = fs::read_to_string(dir.join("expected.bed")).unwrap();
        assert_eq!(got, gold, "precise collapse without an annotation");
        let gtf = fs::read_to_string(dir.join("annot.gtf")).unwrap();
        let got = collapse_precise_bed(&bed, Some(&gtf), &settings).unwrap();
        let gold = fs::read_to_string(dir.join("expected-gtf.bed")).unwrap();
        assert_eq!(got, gold, "precise collapse snapping to annotated ends");
    }

    #[test]
    fn bed12_names_follow_the_flair_split() {
        let bed = "\
chr1\t0\t10\tiso_chr1\t60\t+\t0\t10\t0\t1\t10,\t0,\n\
chr1\t0\t10\tread_NM_001\t60\t+\t0\t10\t0\t1\t10,\t0,\n\
chr1\t0\t10\talone\t60\t+\t0\t10\t0\t1\t10,\t0,\n";
        let gtf = bed12_to_gtf(bed, "FLAIR").unwrap();
        assert!(gtf.contains("transcript_id \"iso\""));
        assert!(gtf.contains("gene_id \"chr1\""));
        assert!(gtf.contains("transcript_id \"read\""));
        assert!(gtf.contains("gene_id \"NM_001\""));
        assert!(gtf.contains("transcript_id \"alone\""));
        assert!(gtf.contains("gene_id \"alone\""));
        assert!(gtf.contains("exon_number \"0\""));
    }
}
