//! Combine FLAIR transcriptomes.
//!
//! Dict order is Python 3 insertion order. Two upstream accidents are kept
//! because they change the files:
//! `biggestdiff` is initialized to 0 and never assigned again, so the "longest"
//! end group is the last group whose length is greater than 0.
//! A chain that fails the filter still names its map entry from `theseisos`
//! left behind by the previous chain, including a filtered chain.

use braid_io::suffix_path;
use std::collections::HashMap;
use std::fs;
use std::path::Path;

#[derive(Clone, Debug)]
pub struct CombineSettings {
    pub end_window: i64,
    /// Integer percent. 10 means usage must be greater than 0.10.
    pub min_percent: i64,
    pub include_single_exon: bool,
    /// `usageandlongest`, `usageonly`, `none`, or a read-count threshold.
    pub filter: String,
    pub convert_gtf: bool,
}

impl Default for CombineSettings {
    fn default() -> Self {
        Self {
            end_window: 200,
            min_percent: 10,
            include_single_exon: false,
            filter: "usageandlongest".to_string(),
            convert_gtf: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct CombineTexts {
    pub bed: String,
    pub counts: String,
    pub isoform_map: String,
    pub fasta: Option<String>,
    pub gtf: Option<String>,
}

#[derive(Clone, Debug)]
struct Iso {
    start: i64,
    end: i64,
    sample: String,
    name: String,
    usage: f64,
    counts: i64,
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
enum ChainBody {
    Introns(Vec<(i64, i64)>),
    Single(String),
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
struct FusionLocus {
    chrom: String,
    strand: String,
    start: Option<i64>,
    end: Option<i64>,
    body: ChainBody,
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
enum ChainId {
    Isoform {
        chrom: String,
        strand: String,
        body: ChainBody,
    },
    Fusion(Vec<FusionLocus>),
}

struct Group {
    rep: Iso,
    members: Vec<Iso>,
}

struct SampleData {
    label: String,
    isos: Vec<(ChainId, Iso)>,
    seqs: HashMap<String, String>,
}

pub fn combine(manifest: &Path, settings: &CombineSettings) -> Result<CombineTexts, String> {
    let count_threshold = parse_count_filter(&settings.filter)?;
    let min_usage = settings.min_percent as f64 / 100.0;
    // Upstream ignores every FASTA unless all samples supply one.
    let generate_fa = samples_have_fasta(manifest)?;
    let samples = load_manifest(manifest, settings.include_single_exon, generate_fa)?;
    if samples.is_empty() {
        return Err(format!(
            "no samples found in manifest file {}",
            manifest.display()
        ));
    }
    let sample_labels: Vec<String> = samples.iter().map(|s| s.label.clone()).collect();

    let mut chains: Vec<ChainId> = Vec::new();
    let mut chain_isos: HashMap<ChainId, Vec<Iso>> = HashMap::new();
    let mut seqs: HashMap<String, HashMap<String, String>> = HashMap::new();
    for sample in &samples {
        if generate_fa {
            seqs.insert(sample.label.clone(), sample.seqs.clone());
        }
        for (chain, iso) in &sample.isos {
            if let Some(bucket) = chain_isos.get_mut(chain) {
                bucket.push(iso.clone());
            } else {
                chains.push(chain.clone());
                chain_isos.insert(chain.clone(), vec![iso.clone()]);
            }
        }
    }

    let mut bed = String::new();
    let mut fa = String::new();
    let mut names: Vec<String> = Vec::new();
    let mut isomap: HashMap<String, Vec<String>> = HashMap::new();
    let mut support_order: Vec<String> = Vec::new();
    let mut support: HashMap<String, HashMap<String, i64>> = HashMap::new();
    let mut last_members: Option<Vec<Iso>> = None;
    let mut iso_count: i64 = 1;

    for chain in &chains {
        let groups = combine_ends(
            chain_isos.remove(chain).unwrap_or_default(),
            settings.end_window,
        )?;
        let single = is_single_exon(chain);
        let fusion = matches!(chain, ChainId::Fusion(_));
        // Last group with a non-zero length. See the module comment.
        let mut longest: Option<(i64, i64)> = None;
        let mut max_chain_usage = 0.0;
        let mut chain_counts = 0i64;
        for group in &groups {
            if (group.rep.end - group.rep.start).abs() > 0 {
                longest = Some((group.rep.start, group.rep.end));
            }
            for iso in &group.members {
                if iso.usage > max_chain_usage {
                    max_chain_usage = iso.usage;
                }
                chain_counts += iso.counts;
            }
        }
        let chain_kept = settings.filter == "none"
            || max_chain_usage > min_usage
            || count_threshold.is_some_and(|n| chain_counts > n);
        if !chain_kept {
            let previous = last_members.as_ref().ok_or_else(|| {
                "flair combine: a filtered intron chain ran before any kept isoform (upstream UnboundLocalError on theseisos)".to_string()
            })?;
            let outname = lowexp_name(previous)?;
            let mut ids = Vec::new();
            for group in &groups {
                ids.extend(group.members.iter().map(source_id));
                last_members = Some(group.members.clone());
            }
            push_map(&mut names, &mut isomap, &outname, ids);
            continue;
        }

        for (end_count, group) in (1i64..).zip(&groups) {
            let mut members = group.members.clone();
            members.sort_by_key(|iso| std::cmp::Reverse(iso.end - iso.start));
            let max_usage = members
                .iter()
                .map(|iso| iso.usage)
                .fold(f64::NEG_INFINITY, f64::max);
            let group_counts: i64 = members.iter().map(|iso| iso.counts).sum();
            let kept = settings.filter == "none"
                || max_usage > min_usage
                || (longest == Some((group.rep.start, group.rep.end))
                    && !single
                    && settings.filter == "usageandlongest")
                || count_threshold.is_some_and(|n| group_counts > n);
            let outname = if kept {
                let name = if fusion {
                    let gene = mode(&member_genes(&members))?;
                    format!("flairiso{iso_count}-{end_count}_{gene}")
                } else if let Some(annotated) = members.iter().find(|iso| annotated_enst(&iso.name))
                {
                    format!("{iso_count}-{end_count}_{}", annotated.name)
                } else {
                    let gene = gene_from_members(&members)?;
                    format!("flairiso{iso_count}-{end_count}_{gene}")
                };
                write_bed(&mut bed, chain, &group.rep, &name)?;
                if generate_fa {
                    let seq = seqs
                        .get(&group.rep.sample)
                        .and_then(|table| table.get(&group.rep.name))
                        .ok_or_else(|| {
                            format!(
                                "sequence for {} is not in the {} fasta",
                                group.rep.name, group.rep.sample
                            )
                        })?;
                    fa.push('>');
                    fa.push_str(&name);
                    fa.push('\n');
                    fa.push_str(seq);
                    fa.push('\n');
                }
                name
            } else {
                lowexp_name(&members)?
            };
            push_map(
                &mut names,
                &mut isomap,
                &outname,
                members.iter().map(source_id).collect(),
            );
            add_support(
                &mut support_order,
                &mut support,
                &sample_labels,
                &outname,
                &members,
            );
            last_members = Some(members);
        }
        iso_count += 1;
    }

    let mut isoform_map = String::new();
    for name in &names {
        let ids = isomap.get(name).map(Vec::as_slice).unwrap_or(&[]);
        isoform_map.push_str(name);
        isoform_map.push('\t');
        isoform_map.push_str(&ids.join("\t"));
        isoform_map.push('\n');
    }
    let mut counts = String::new();
    counts.push_str("ids");
    for label in &sample_labels {
        counts.push('\t');
        counts.push_str(label);
    }
    counts.push('\n');
    for name in &support_order {
        counts.push_str(name);
        let row = &support[name];
        for label in &sample_labels {
            counts.push('\t');
            counts.push_str(&row.get(label).copied().unwrap_or(0).to_string());
        }
        counts.push('\n');
    }

    let gtf = if settings.convert_gtf {
        Some(bed_to_gtf(&bed)?)
    } else {
        None
    };
    Ok(CombineTexts {
        bed,
        counts,
        isoform_map,
        fasta: if generate_fa { Some(fa) } else { None },
        gtf,
    })
}

pub fn write_combine_texts(prefix: &Path, texts: &CombineTexts) -> Result<(), String> {
    write_one(&suffix_path(prefix, ".bed"), &texts.bed)?;
    write_one(&suffix_path(prefix, ".counts.tsv"), &texts.counts)?;
    write_one(&suffix_path(prefix, ".isoform.map.txt"), &texts.isoform_map)?;
    if let Some(fa) = &texts.fasta {
        write_one(&suffix_path(prefix, ".fa"), fa)?;
    }
    if let Some(gtf) = &texts.gtf {
        write_one(&suffix_path(prefix, ".gtf"), gtf)?;
    }
    Ok(())
}

fn write_one(path: &Path, text: &str) -> Result<(), String> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)
                .map_err(|err| format!("failed to create {}: {err}", parent.display()))?;
        }
    }
    fs::write(path, text).map_err(|err| format!("failed to write {}: {err}", path.display()))
}

fn samples_have_fasta(manifest: &Path) -> Result<bool, String> {
    let text = read_text(manifest)?;
    let mut saw = false;
    let mut all = true;
    for line in text.split_terminator('\n') {
        let line = py_rstrip(line.trim_end_matches('\r'));
        if line.is_empty() && text_is_blank_row(line) {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if !(3..=5).contains(&cols.len()) {
            return Err(format!(
                "Expected between 3 to 5 columns in manifest, got {} in {}",
                cols.len(),
                manifest.display()
            ));
        }
        saw = true;
        if cols.len() <= 3 || cols[3].is_empty() {
            all = false;
        }
    }
    Ok(saw && all)
}

fn text_is_blank_row(line: &str) -> bool {
    line.is_empty()
}

fn load_manifest(
    manifest: &Path,
    include_single_exon: bool,
    generate_fa: bool,
) -> Result<Vec<SampleData>, String> {
    let text = read_text(manifest)?;
    let mut samples = Vec::new();
    for raw in text.split_terminator('\n') {
        let line = py_rstrip(raw.trim_end_matches('\r'));
        if line.is_empty() {
            return Err(format!(
                "Expected between 3 to 5 columns in manifest, got 1 in {}",
                manifest.display()
            ));
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if !(3..=5).contains(&cols.len()) {
            return Err(format!(
                "Expected between 3 to 5 columns in manifest, got {} in {}",
                cols.len(),
                manifest.display()
            ));
        }
        let label = format!("{}__{}", cols[0], cols[1]);
        let fa_path = if cols.len() > 3 { cols[3] } else { "" };
        let map_path = if cols.len() > 4 { cols[4] } else { "" };
        let support = if map_path.is_empty() {
            None
        } else {
            Some(read_support(Path::new(map_path))?)
        };
        let fusion = cols[1] == "fusionisoform";
        let isos = if fusion {
            read_fusion_bed(Path::new(cols[2]), &label, support.as_ref())?
        } else {
            read_isoform_bed(
                Path::new(cols[2]),
                &label,
                support.as_ref(),
                include_single_exon,
            )?
        };
        let seqs = if !generate_fa {
            HashMap::new()
        } else {
            read_fasta_table(Path::new(fa_path))?
        };
        samples.push(SampleData { label, isos, seqs });
    }
    Ok(samples)
}

struct Support {
    iso_reads: HashMap<String, i64>,
    gene_reads: HashMap<String, i64>,
}

fn read_support(path: &Path) -> Result<Support, String> {
    let text = read_text(path)?;
    let mut iso_reads = HashMap::new();
    let mut gene_reads: HashMap<String, i64> = HashMap::new();
    for raw in text.split_terminator('\n') {
        let line = raw.trim_end_matches('\r');
        if line.is_empty() {
            continue;
        }
        let (iso, reads) = line.split_once('\t').ok_or_else(|| {
            format!(
                "read map {} is missing a tab in line {line}",
                path.display()
            )
        })?;
        let gene = last_token(iso);
        let n = reads.split(',').count() as i64;
        *gene_reads.entry(gene.to_string()).or_insert(0) += n;
        iso_reads.insert(iso.to_string(), n);
    }
    Ok(Support {
        iso_reads,
        gene_reads,
    })
}

fn usage_of(support: Option<&Support>, name: &str) -> Result<(f64, i64), String> {
    let Some(support) = support else {
        return Ok((1.0, 0));
    };
    let counts = *support
        .iso_reads
        .get(name)
        .ok_or_else(|| format!("isoform {name} is not in the read map"))?;
    let gene = last_token(name);
    let total = *support.gene_reads.get(gene).unwrap_or(&0);
    if total == 0 {
        return Err(format!("gene {gene} has no reads in the read map"));
    }
    Ok((counts as f64 / total as f64, counts))
}

fn read_isoform_bed(
    path: &Path,
    sample: &str,
    support: Option<&Support>,
    include_single_exon: bool,
) -> Result<Vec<(ChainId, Iso)>, String> {
    let text = read_text(path)?;
    let mut out = Vec::new();
    for raw in text.split_terminator('\n') {
        let line = py_rstrip(raw.trim_end_matches('\r'));
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let BedFields {
            chrom,
            start,
            end,
            name,
            strand,
            exons,
            sizes,
            rels,
        } = bed_fields(&cols, path)?;
        let mut keep_name = name.to_string();
        let body = if exons > 1 {
            ChainBody::Introns(intron_chain(start, &sizes, &rels))
        } else if include_single_exon {
            ChainBody::Single(format!("{chrom}-{}", round_to_10k(start)))
        } else {
            continue;
        };
        let (usage, counts) = usage_of(support, &keep_name)?;
        keep_name = clean_iso_name(&keep_name);
        let chain = ChainId::Isoform {
            chrom: chrom.to_string(),
            strand: strand.to_string(),
            body,
        };
        out.push((
            chain,
            Iso {
                start,
                end,
                sample: sample.to_string(),
                name: keep_name,
                usage,
                counts,
            },
        ));
    }
    Ok(out)
}

fn read_fusion_bed(
    path: &Path,
    sample: &str,
    support: Option<&Support>,
) -> Result<Vec<(ChainId, Iso)>, String> {
    let text = read_text(path)?;
    let mut order: Vec<String> = Vec::new();
    let mut loci: HashMap<String, Vec<FusionLocus>> = HashMap::new();
    for raw in text.split_terminator('\n') {
        let line = py_rstrip(raw.trim_end_matches('\r'));
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        let BedFields {
            chrom,
            start,
            end,
            name,
            strand,
            exons,
            sizes,
            rels,
        } = bed_fields(&cols, path)?;
        let stripped = name.split('_').skip(1).collect::<Vec<_>>().join("_");
        let body = if exons > 1 {
            ChainBody::Introns(intron_chain(start, &sizes, &rels))
        } else {
            ChainBody::Single(format!("{chrom}-{}", round_to_10k(start)))
        };
        let locus = FusionLocus {
            chrom: chrom.to_string(),
            strand: strand.to_string(),
            start: Some(start),
            end: Some(end),
            body,
        };
        if let Some(bucket) = loci.get_mut(&stripped) {
            bucket.push(locus);
        } else {
            order.push(stripped.clone());
            loci.insert(stripped, vec![locus]);
        }
    }
    let mut out = Vec::new();
    for name in order {
        let mut loci = loci.remove(&name).unwrap_or_default();
        if loci.is_empty() {
            continue;
        }
        let fusion_start = if loci[0].strand == "+" {
            loci[0]
                .start
                .take()
                .ok_or("fusion locus is missing a start")?
        } else {
            loci[0].end.take().ok_or("fusion locus is missing an end")?
        };
        let last = loci.len() - 1;
        let fusion_end = if loci[last].strand == "+" {
            loci[last]
                .end
                .take()
                .ok_or("fusion locus is missing an end")?
        } else {
            loci[last]
                .start
                .take()
                .ok_or("fusion locus is missing a start")?
        };
        let (usage, counts) = usage_of(support, &name)?;
        let stored = clean_iso_name(&name);
        out.push((
            ChainId::Fusion(loci),
            Iso {
                start: fusion_start,
                end: fusion_end,
                sample: sample.to_string(),
                name: stored,
                usage,
                counts,
            },
        ));
    }
    Ok(out)
}

struct BedFields<'a> {
    chrom: &'a str,
    start: i64,
    end: i64,
    name: &'a str,
    strand: &'a str,
    exons: i64,
    sizes: Vec<i64>,
    rels: Vec<i64>,
}

fn bed_fields<'a>(cols: &'a [&str], path: &Path) -> Result<BedFields<'a>, String> {
    if cols.len() < 12 {
        return Err(format!(
            "{} has a BED row with {} columns; FLAIR combine needs 12",
            path.display(),
            cols.len()
        ));
    }
    let start = parse_i64(cols[1], path)?;
    let end = parse_i64(cols[2], path)?;
    let exons = parse_i64(cols[9], path)?;
    let sizes = parse_csv_i64(cols[10], path)?;
    let rels = parse_csv_i64(cols[11], path)?;
    Ok(BedFields {
        chrom: cols[0],
        start,
        end,
        name: cols[3],
        strand: cols[5],
        exons,
        sizes,
        rels,
    })
}

fn intron_chain(start: i64, sizes: &[i64], rels: &[i64]) -> Vec<(i64, i64)> {
    let mut introns = Vec::new();
    for i in 0..sizes.len().saturating_sub(1) {
        introns.push((start + rels[i] + sizes[i], start + rels[i + 1]));
    }
    introns
}

fn read_fasta_table(path: &Path) -> Result<HashMap<String, String>, String> {
    let text = read_text(path)?;
    let mut table = HashMap::new();
    let mut last: Option<String> = None;
    for raw in text.split_terminator('\n') {
        let line = raw.trim_end_matches('\r');
        // Like upstream, even an empty sequence line replaces the previous
        // line for this header; FASTA records are not concatenated here.
        if let Some(header) = line.strip_prefix('>') {
            last = Some(py_rstrip(header).to_string());
        } else {
            let header = last.as_ref().ok_or_else(|| {
                format!(
                    "fasta {} has a sequence line before a header",
                    path.display()
                )
            })?;
            let key = clean_iso_name(header);
            table.insert(key.clone(), py_rstrip(line).to_string());
            last = Some(key);
        }
    }
    Ok(table)
}

fn combine_ends(mut isos: Vec<Iso>, window: i64) -> Result<Vec<Group>, String> {
    isos.sort_by(|a, b| {
        a.start
            .cmp(&b.start)
            .then(a.end.cmp(&b.end))
            .then(a.sample.cmp(&b.sample))
            .then(a.name.cmp(&b.name))
            .then(cmp_usage(a.usage, b.usage))
            .then(a.counts.cmp(&b.counts))
    });
    let mut groups = Vec::new();
    let mut current: Vec<Iso> = Vec::new();
    let mut last_start = 0i64;
    let mut last_end = 0i64;
    for iso in isos {
        let start = iso.start;
        let end = iso.end;
        if start - last_start <= window && end - last_end <= window {
            current.push(iso);
        } else {
            if !current.is_empty() {
                push_group(&mut groups, &current)?;
            }
            current = vec![iso];
        }
        last_start = start;
        last_end = end;
    }
    if !current.is_empty() {
        push_group(&mut groups, &current)?;
    }
    Ok(groups)
}

fn push_group(groups: &mut Vec<Group>, members: &[Iso]) -> Result<(), String> {
    let rep = best_ends(members)?;
    if let Some(existing) = groups.iter_mut().find(|group| same_iso(&group.rep, &rep)) {
        existing.members = members.to_vec();
    } else {
        groups.push(Group {
            rep,
            members: members.to_vec(),
        });
    }
    Ok(())
}

fn best_ends(members: &[Iso]) -> Result<Iso, String> {
    let mut best: Option<Iso> = None;
    for iso in members {
        match &best {
            None => {
                if iso.usage > 0.0 {
                    best = Some(iso.clone());
                } else if iso.usage == 0.0 {
                    return Err(
                        "flair combine: isoform usage is 0, which raises TypeError in getbestends"
                            .to_string(),
                    );
                }
            }
            Some(chosen) => {
                let longer = iso.end - iso.start > chosen.end - chosen.start;
                if iso.usage > chosen.usage || (iso.usage == chosen.usage && longer) {
                    best = Some(iso.clone());
                }
            }
        }
    }
    best.ok_or_else(|| "flair combine: no isoform end with usage > 0".to_string())
}

fn same_iso(left: &Iso, right: &Iso) -> bool {
    left.start == right.start
        && left.end == right.end
        && left.sample == right.sample
        && left.name == right.name
        && left.usage == right.usage
        && left.counts == right.counts
}

fn cmp_usage(left: f64, right: f64) -> std::cmp::Ordering {
    if left < right {
        std::cmp::Ordering::Less
    } else if left > right {
        std::cmp::Ordering::Greater
    } else {
        std::cmp::Ordering::Equal
    }
}

fn is_single_exon(chain: &ChainId) -> bool {
    matches!(
        chain,
        ChainId::Isoform {
            body: ChainBody::Single(_),
            ..
        }
    )
}

fn write_bed(out: &mut String, chain: &ChainId, rep: &Iso, name: &str) -> Result<(), String> {
    match chain {
        ChainId::Fusion(loci) => {
            let mut loci = loci.clone();
            if loci.is_empty() {
                return Err("fusion isoform has no loci".to_string());
            }
            if loci[0].strand == "+" {
                loci[0].start = Some(rep.start);
            } else {
                loci[0].end = Some(rep.start);
            }
            let last = loci.len() - 1;
            if loci[last].strand == "+" {
                loci[last].end = Some(rep.end);
            } else {
                loci[last].start = Some(rep.end);
            }
            for (index, locus) in loci.iter().enumerate() {
                let start = locus
                    .start
                    .ok_or_else(|| format!("fusion locus {} is missing a start", index + 1))?;
                let end = locus
                    .end
                    .ok_or_else(|| format!("fusion locus {} is missing an end", index + 1))?;
                let (sizes, rels) = match &locus.body {
                    ChainBody::Single(_) => (vec![end - start], vec![0]),
                    ChainBody::Introns(introns) => exon_blocks(introns, start, end),
                };
                let locus_name = format!("fusiongene{}_{name}", index + 1);
                push_bed_line(
                    out,
                    &locus.chrom,
                    (start, end),
                    &locus_name,
                    &locus.strand,
                    &sizes,
                    &rels,
                );
            }
        }
        ChainId::Isoform {
            chrom,
            strand,
            body,
        } => {
            let (sizes, rels) = match body {
                ChainBody::Single(_) => (vec![rep.end - rep.start], vec![0]),
                ChainBody::Introns(introns) => exon_blocks(introns, rep.start, rep.end),
            };
            push_bed_line(
                out,
                chrom,
                (rep.start, rep.end),
                name,
                strand,
                &sizes,
                &rels,
            );
        }
    }
    Ok(())
}

fn exon_blocks(introns: &[(i64, i64)], start: i64, end: i64) -> (Vec<i64>, Vec<i64>) {
    let mut sizes = Vec::new();
    let mut rels = vec![0i64];
    for (intron_start, intron_end) in introns {
        sizes.push(intron_start - (start + rels[rels.len() - 1]));
        rels.push(intron_end - start);
    }
    sizes.push(end - (start + rels[rels.len() - 1]));
    (sizes, rels)
}

fn push_bed_line(
    out: &mut String,
    chrom: &str,
    (start, end): (i64, i64),
    name: &str,
    strand: &str,
    sizes: &[i64],
    rels: &[i64],
) {
    out.push_str(&format!(
        "{chrom}\t{start}\t{end}\t{name}\t1000\t{strand}\t{start}\t{end}\t0\t{}\t{},\t{},\n",
        sizes.len(),
        join_i64(sizes),
        join_i64(rels),
    ));
}

fn join_i64(values: &[i64]) -> String {
    values
        .iter()
        .map(i64::to_string)
        .collect::<Vec<_>>()
        .join(",")
}

fn annotated_enst(name: &str) -> bool {
    let mut parts = name.split("ENSG");
    let prefix = parts.next().unwrap_or("");
    let one_ensg = parts.next().is_some() && parts.next().is_none();
    name.starts_with("ENST") && prefix.len() < 25 && one_ensg
}

fn gene_from_members(members: &[Iso]) -> Result<String, String> {
    let mut gene: Option<String> = None;
    for iso in members {
        if let Some(at) = iso.name.rfind("ENSG") {
            gene = Some(format!("ENSG{}", &iso.name[at + 4..]));
        }
        let is_ensg = gene.as_ref().is_some_and(|value| value.starts_with("ENSG"));
        if gene.is_none() || !is_ensg {
            if let Some(at) = iso.name.rfind("chr") {
                gene = Some(format!("chr{}", &iso.name[at + 3..]));
            }
        }
    }
    if let Some(gene) = gene {
        Ok(gene)
    } else {
        mode(&member_genes(members))
    }
}

fn lowexp_name(members: &[Iso]) -> Result<String, String> {
    let mut gene: Option<String> = None;
    for iso in members {
        let token = last_token(&iso.name);
        if token.starts_with("ENSG") {
            gene = Some(token.to_string());
        }
    }
    let gene = match gene {
        Some(gene) => gene,
        None => mode(&member_genes(members))?,
    };
    Ok(format!("lowexpiso_{gene}"))
}

fn member_genes(members: &[Iso]) -> Vec<String> {
    members
        .iter()
        .map(|iso| last_token(&iso.name).to_string())
        .collect()
}

fn mode(values: &[String]) -> Result<String, String> {
    if values.is_empty() {
        return Err("flair combine: statistics.mode() of empty data".to_string());
    }
    let mut counts: HashMap<&str, usize> = HashMap::new();
    let mut order = Vec::new();
    for value in values {
        let count = counts.entry(value.as_str()).or_insert(0);
        if *count == 0 {
            order.push(value.as_str());
        }
        *count += 1;
    }
    let max = counts.values().copied().max().unwrap_or(0);
    order
        .into_iter()
        .find(|value| counts[value] == max)
        .map(str::to_string)
        .ok_or_else(|| "flair combine: statistics.mode() failed".to_string())
}

fn source_id(iso: &Iso) -> String {
    format!("{}..{}", iso.sample, iso.name)
}

fn push_map(
    order: &mut Vec<String>,
    map: &mut HashMap<String, Vec<String>>,
    name: &str,
    ids: Vec<String>,
) {
    if let Some(bucket) = map.get_mut(name) {
        bucket.extend(ids);
    } else {
        order.push(name.to_string());
        map.insert(name.to_string(), ids);
    }
}

fn add_support(
    order: &mut Vec<String>,
    support: &mut HashMap<String, HashMap<String, i64>>,
    samples: &[String],
    name: &str,
    members: &[Iso],
) {
    if !support.contains_key(name) {
        order.push(name.to_string());
        let mut row = HashMap::new();
        for sample in samples {
            row.entry(sample.clone()).or_insert(0);
        }
        support.insert(name.to_string(), row);
    }
    let row = support.get_mut(name).expect("row was just inserted");
    for iso in members {
        *row.entry(iso.sample.clone()).or_insert(0) += iso.counts;
    }
}

fn clean_iso_name(name: &str) -> String {
    name.replace("_PAR_Y", "")
}

fn last_token(name: &str) -> &str {
    name.rsplit('_').next().unwrap_or(name)
}

/// Python 3 `round(start, -4)`: half to even, then `int`.
fn round_to_10k(start: i64) -> i64 {
    let unit = 10_000;
    let quot = start.div_euclid(unit);
    let rem = start.rem_euclid(unit);
    if rem < 5_000 {
        quot * unit
    } else if rem > 5_000 {
        (quot + 1) * unit
    } else if quot % 2 == 0 {
        quot * unit
    } else {
        (quot + 1) * unit
    }
}

fn parse_count_filter(filter: &str) -> Result<Option<i64>, String> {
    if filter.is_empty() || !filter.bytes().all(|byte| byte.is_ascii_digit()) {
        return Ok(None);
    }
    filter
        .parse::<i64>()
        .map(Some)
        .map_err(|_| format!("filter count does not fit an i64: {filter}"))
}

fn parse_i64(text: &str, path: &Path) -> Result<i64, String> {
    text.parse::<i64>()
        .map_err(|_| format!("bad integer {text} in {}", path.display()))
}

fn parse_csv_i64(text: &str, path: &Path) -> Result<Vec<i64>, String> {
    let text = text.trim_end_matches(',');
    if text.is_empty() {
        return Ok(Vec::new());
    }
    text.split(',')
        .map(|piece| parse_i64(piece, path))
        .collect()
}

fn py_rstrip(text: &str) -> &str {
    text.trim_end_matches([' ', '\t', '\n', '\r', '\u{000b}', '\u{000c}'])
}

fn read_text(path: &Path) -> Result<String, String> {
    fs::read_to_string(path).map_err(|err| format!("failed to read {}: {err}", path.display()))
}

fn bed_to_gtf(bed: &str) -> Result<String, String> {
    let mut out = String::new();
    for raw in bed.split_terminator('\n') {
        let line = py_rstrip(raw.trim_end_matches('\r'));
        if line.is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 12 {
            return Err("combined bed is missing columns for gtf conversion".to_string());
        }
        let start: i64 = cols[1]
            .parse()
            .map_err(|_| format!("bad bed start {}", cols[1]))?;
        let thick_start: i64 = cols[6]
            .parse()
            .map_err(|_| format!("bad thick start {}", cols[6]))?;
        let thick_end: i64 = cols[7]
            .parse()
            .map_err(|_| format!("bad thick end {}", cols[7]))?;
        let name = cols[3].replace(';', ":");
        if !name.contains('_') {
            return Err(
                "Entry name should contain underscore-delimited transcriptid and geneid; no GTF conversion was done"
                    .to_string(),
            );
        }
        let (transcript_id, gene_id) = split_iso_gene(&name);
        let attributes = format!("gene_id \"{gene_id}\"; transcript_id \"{transcript_id}\";");
        let rels = parse_csv_i64(cols[11], Path::new("combined.bed"))?;
        let sizes = parse_csv_i64(cols[10], Path::new("combined.bed"))?;
        if rels.is_empty() || sizes.is_empty() || rels.len() != sizes.len() {
            return Err("combined bed block columns do not match".to_string());
        }
        let tstarts: Vec<i64> = rels.iter().map(|rel| rel + start).collect();
        let transcript_end = tstarts[tstarts.len() - 1] + sizes[sizes.len() - 1];
        push_gtf(
            &mut out,
            cols[0],
            "transcript",
            start + 1,
            transcript_end,
            cols[5],
            &attributes,
        );
        let transcript_bed_end: i64 = cols[2]
            .parse()
            .map_err(|_| format!("bad bed end {}", cols[2]))?;
        if thick_start != thick_end && (thick_start != start || thick_end != transcript_bed_end) {
            push_gtf(
                &mut out,
                cols[0],
                "CDS",
                thick_start + 1,
                thick_end,
                cols[5],
                &attributes,
            );
            if cols[5] == "+" {
                push_gtf(
                    &mut out,
                    cols[0],
                    "start_codon",
                    thick_start + 1,
                    thick_start + 3,
                    "+",
                    &attributes,
                );
                push_gtf(
                    &mut out,
                    cols[0],
                    "5UTR",
                    start + 1,
                    thick_start + 1,
                    "+",
                    &attributes,
                );
                push_gtf(
                    &mut out,
                    cols[0],
                    "3UTR",
                    thick_end,
                    transcript_end,
                    "+",
                    &attributes,
                );
            } else if cols[5] == "-" {
                push_gtf(
                    &mut out,
                    cols[0],
                    "start_codon",
                    thick_end - 2,
                    thick_end,
                    "-",
                    &attributes,
                );
                push_gtf(
                    &mut out,
                    cols[0],
                    "3UTR",
                    start + 1,
                    thick_start + 1,
                    "-",
                    &attributes,
                );
                push_gtf(
                    &mut out,
                    cols[0],
                    "5UTR",
                    thick_end,
                    transcript_end,
                    "-",
                    &attributes,
                );
            }
        }
        for (index, (tstart, size)) in tstarts.iter().zip(&sizes).enumerate() {
            let attributes = format!(
                "gene_id \"{gene_id}\"; transcript_id \"{transcript_id}\"; exon_number \"{index}\";"
            );
            push_gtf(
                &mut out,
                cols[0],
                "exon",
                tstart + 1,
                tstart + size,
                cols[5],
                &attributes,
            );
        }
    }
    Ok(out)
}

fn push_gtf(
    out: &mut String,
    chrom: &str,
    kind: &str,
    start: i64,
    end: i64,
    strand: &str,
    attributes: &str,
) {
    out.push_str(&format!(
        "{chrom}\tFLAIR\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n"
    ));
}

fn split_iso_gene(name: &str) -> (String, String) {
    // The marker checks are the upstream order. The gene is everything after
    // the underscore that starts the matched marker.
    let from = ["_chr", "_XM", "_XR", "_NM", "_NR", "_R2_"]
        .iter()
        .find_map(|mark| name.rfind(mark))
        .or_else(|| name.rfind('_'));
    match from {
        Some(at) => (name[..at].to_string(), name[at + 1..].to_string()),
        None => (name.to_string(), name.to_string()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn root() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../tests/parity/flair_combine")
    }

    fn write_manifest(name: &str, lines: &[String]) -> PathBuf {
        let path = std::env::temp_dir().join(format!("braid-flair-{name}.txt"));
        fs::write(&path, lines.join("\n") + "\n").unwrap();
        path
    }

    fn assert_text(got: &str, path: &Path) {
        let gold =
            fs::read_to_string(path).unwrap_or_else(|err| panic!("read {}: {err}", path.display()));
        if got != gold {
            let mut shown = 0;
            for (index, (left, right)) in gold.lines().zip(got.lines()).enumerate() {
                if left != right {
                    panic!(
                        "{} line {} differs\n gold: {left}\n  got: {right}",
                        path.display(),
                        index + 1
                    );
                }
                shown = index + 1;
            }
            panic!(
                "{} length differs after {shown} shared lines (gold {}, got {})",
                path.display(),
                gold.lines().count(),
                got.lines().count()
            );
        }
    }

    fn assert_outputs(texts: &CombineTexts, dir: &Path, fasta: bool, gtf: bool) {
        assert_text(&texts.bed, &dir.join("combined.bed"));
        assert_text(&texts.counts, &dir.join("combined.counts.tsv"));
        assert_text(&texts.isoform_map, &dir.join("combined.isoform.map.txt"));
        if fasta {
            assert_text(
                texts.fasta.as_deref().expect("fasta"),
                &dir.join("combined.fa"),
            );
        } else {
            assert!(texts.fasta.is_none());
        }
        if gtf {
            assert_text(
                texts.gtf.as_deref().expect("gtf"),
                &dir.join("combined.gtf"),
            );
        } else {
            assert!(texts.gtf.is_none());
        }
    }

    #[test]
    fn official_combine_matches_flair() {
        let dir = root().join("official");
        let bed = dir.join("in/collapse.isoforms.bed");
        let fa = dir.join("in/collapse.isoforms.fa");
        let map = dir.join("in/collapse.isoform.read.map.txt");
        let cases = [
            (
                "bed-only",
                vec![format!("A1\tisoforms\t{}", bed.display())],
                false,
            ),
            (
                "bed-fa",
                vec![format!("A1\tisoforms\t{}\t{}", bed.display(), fa.display())],
                true,
            ),
            (
                "bed-fa-map",
                vec![format!(
                    "A1\tisoforms\t{}\t{}\t{}",
                    bed.display(),
                    fa.display(),
                    map.display()
                )],
                true,
            ),
        ];
        for (name, lines, fasta) in cases {
            let manifest = write_manifest(name, &lines);
            let texts = combine(&manifest, &CombineSettings::default()).unwrap();
            assert_outputs(&texts, &dir.join(name), fasta, false);
        }
    }

    #[test]
    fn synthetic_combine_matches_flair() {
        let dir = root().join("synthetic");
        let inn = dir.join("in");
        let main = vec![
            format!(
                "S1\tisoforms\t{}\t{}\t{}",
                inn.join("s1.bed").display(),
                inn.join("s1.fa").display(),
                inn.join("s1.map").display()
            ),
            format!(
                "S2\tisoforms\t{}\t{}\t{}",
                inn.join("s2.bed").display(),
                inn.join("s2.fa").display(),
                inn.join("s2.map").display()
            ),
        ];
        let cases = [
            (
                "usageandlongest",
                main.clone(),
                CombineSettings {
                    convert_gtf: true,
                    ..CombineSettings::default()
                },
                true,
                true,
            ),
            (
                "usageonly",
                main.clone(),
                CombineSettings {
                    filter: "usageonly".to_string(),
                    ..CombineSettings::default()
                },
                true,
                false,
            ),
            (
                "filternone",
                main,
                CombineSettings {
                    filter: "none".to_string(),
                    ..CombineSettings::default()
                },
                true,
                false,
            ),
            (
                "dropped",
                vec![format!(
                    "D\tisoforms\t{}\t{}\t{}",
                    inn.join("drop.bed").display(),
                    inn.join("drop.fa").display(),
                    inn.join("drop.map").display()
                )],
                CombineSettings::default(),
                true,
                false,
            ),
            (
                "dropped_none",
                vec![format!(
                    "D\tisoforms\t{}\t{}\t{}",
                    inn.join("drop.bed").display(),
                    inn.join("drop.fa").display(),
                    inn.join("drop.map").display()
                )],
                CombineSettings {
                    filter: "none".to_string(),
                    ..CombineSettings::default()
                },
                true,
                false,
            ),
            (
                "numeric",
                vec![format!(
                    "N\tisoforms\t{}\t{}\t{}",
                    inn.join("num.bed").display(),
                    inn.join("num.fa").display(),
                    inn.join("num.map").display()
                )],
                CombineSettings {
                    filter: "3".to_string(),
                    ..CombineSettings::default()
                },
                true,
                false,
            ),
            (
                "se_off",
                vec![format!(
                    "E\tisoforms\t{}\t{}\t{}",
                    inn.join("se.bed").display(),
                    inn.join("se.fa").display(),
                    inn.join("se.map").display()
                )],
                CombineSettings::default(),
                true,
                false,
            ),
            (
                "se_on",
                vec![format!(
                    "E\tisoforms\t{}\t{}\t{}",
                    inn.join("se.bed").display(),
                    inn.join("se.fa").display(),
                    inn.join("se.map").display()
                )],
                CombineSettings {
                    include_single_exon: true,
                    ..CombineSettings::default()
                },
                true,
                false,
            ),
            (
                "fusion",
                vec![format!(
                    "F\tfusionisoform\t{}\t{}\t{}",
                    inn.join("fus.bed").display(),
                    inn.join("fus.fa").display(),
                    inn.join("fus.map").display()
                )],
                CombineSettings {
                    convert_gtf: true,
                    ..CombineSettings::default()
                },
                true,
                true,
            ),
        ];
        for (name, lines, settings, fasta, gtf) in cases {
            let manifest = write_manifest(name, &lines);
            let texts = combine(&manifest, &settings).unwrap_or_else(|err| panic!("{name}: {err}"));
            assert_outputs(&texts, &dir.join(name), fasta, gtf);
        }
    }

    #[test]
    fn filtered_chain_before_any_kept_isoform_errors() {
        let dir = std::env::temp_dir().join("braid-flair-first-drop");
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::write(
            dir.join("one.bed"),
            "chr1\t100\t400\ttiny_ENSG1\t40\t+\t100\t400\t0\t2\t100,100,\t0,200,\n",
        )
        .unwrap();
        fs::write(
            dir.join("one.map"),
            format!("tiny_ENSG1\tt1\nbulk_ENSG1\t{}", vec!["b"; 50].join(",")),
        )
        .unwrap();
        let manifest = write_manifest(
            "first-drop",
            &[format!(
                "D\tisoforms\t{}\t\t{}",
                dir.join("one.bed").display(),
                dir.join("one.map").display()
            )],
        );
        let err = combine(&manifest, &CombineSettings::default()).unwrap_err();
        assert!(err.contains("UnboundLocalError"), "{err}");
    }

    #[test]
    fn consecutive_filtered_chains_match_flair() {
        let dir = root().join("regressions/consecutive-filtered-chains");
        let manifest = write_manifest(
            "consecutive-filtered-chains",
            &[format!(
                "S\tisoforms\t{}\t\t{}",
                dir.join("in.bed").display(),
                dir.join("in.map").display()
            )],
        );
        let texts = combine(&manifest, &CombineSettings::default()).unwrap();
        assert_outputs(&texts, &dir, false, false);
    }

    #[test]
    fn partial_fasta_manifest_ignores_all_fasta_paths() {
        let dir = root().join("regressions/mixed-fasta-missing-unused-file");
        let manifest = write_manifest(
            "mixed-fasta-missing-unused-file",
            &[
                format!(
                    "A\tisoforms\t{}\t{}",
                    dir.join("in.bed").display(),
                    dir.join("absent.fa").display()
                ),
                format!("B\tisoforms\t{}", dir.join("in.bed").display()),
            ],
        );
        let texts = combine(&manifest, &CombineSettings::default()).unwrap();
        assert_outputs(&texts, &dir, false, false);
    }

    #[test]
    fn fasta_blank_line_matches_flair() {
        let dir = root().join("regressions/fasta-blank-line");
        let manifest = write_manifest(
            "fasta-blank-line",
            &[format!(
                "A\tisoforms\t{}\t{}",
                dir.join("in.bed").display(),
                dir.join("in.fa").display()
            )],
        );
        let texts = combine(&manifest, &CombineSettings::default()).unwrap();
        assert_outputs(&texts, &dir, true, false);
    }

    #[test]
    fn round_to_10k_matches_python3() {
        assert_eq!(round_to_10k(0), 0);
        assert_eq!(round_to_10k(4999), 0);
        assert_eq!(round_to_10k(5000), 0);
        assert_eq!(round_to_10k(5001), 10_000);
        assert_eq!(round_to_10k(14999), 10_000);
        assert_eq!(round_to_10k(15000), 20_000);
        assert_eq!(round_to_10k(25000), 20_000);
    }
}
