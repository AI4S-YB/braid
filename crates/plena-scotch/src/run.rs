//! Quantify one sample group per BAM or SAM input.
//!
//! Reads that share a cell barcode and UMI keep the longest alignment.
//! Samples share the reference annotation and discover novel isoforms
//! jointly. Counts stay per sample. A novel isoform is reported only when a
//! discovery chunk assigns at least 10 reads.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

use plena_io::read_fasta;

use crate::assign::{
    detect_edge_poly, homopolymer, place_in_metagene, query_poly, MatchParams, Outcome, Placement,
};
use crate::bam::{is_unmapped, open_alignments, Read};
use crate::coverage::{refine_gene, Coverage};
use crate::gtf::{gtf_exons, load_gtf, retain_genes, Annotation, Gene};
use crate::novel::{discover, group_novel, NovelCall};

const BUCKET: i64 = 100_000;
const NOVEL_MIN_READS: usize = 10;
const NOVEL_MAX_ROUNDS: usize = 50;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Platform {
    TenxOnt,
    TenxPacbio,
    ParseOnt,
    Bulk,
}

#[derive(Clone, Debug)]
pub struct ScotchSettings {
    pub bams: Vec<PathBuf>,
    pub gtf: PathBuf,
    pub out: PathBuf,
    pub fasta: Option<PathBuf>,
    pub update_gtf: bool,
    pub platform: Platform,
    pub barcode_cell: Option<String>,
    pub barcode_umi: Option<String>,
    pub genes: Vec<String>,
    pub params: MatchParams,
    pub novel_read_n: usize,
    pub novel_read_pct: f64,
    pub group_novel: bool,
    pub exon_fraction: f64,
    pub splice_fraction: f64,
    pub z_score: f64,
    pub workers: usize,
}

impl Default for ScotchSettings {
    fn default() -> Self {
        Self {
            bams: Vec::new(),
            gtf: PathBuf::new(),
            out: PathBuf::new(),
            fasta: None,
            update_gtf: false,
            platform: Platform::TenxOnt,
            barcode_cell: None,
            barcode_umi: None,
            genes: Vec::new(),
            params: MatchParams::default(),
            novel_read_n: 0,
            novel_read_pct: 0.0,
            group_novel: true,
            exon_fraction: 0.02,
            splice_fraction: 0.02,
            z_score: 10.0,
            workers: 1,
        }
    }
}

#[derive(Clone, Debug)]
pub struct ScotchFailure {
    pub message: String,
    pub usage: bool,
}

impl std::fmt::Display for ScotchFailure {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.message)
    }
}

struct Sample {
    name: String,
    files: Vec<PathBuf>,
}

struct Meta {
    start: i64,
    end: i64,
    genes: Vec<usize>,
}

struct Hit {
    sample: usize,
    qname: String,
    cell: String,
    gene: usize,
    outcome: Outcome,
    score: f64,
}

struct Winner {
    length: i64,
    ordinal: u64,
}

struct Identity {
    key: String,
    cell: String,
    length: i64,
}

#[derive(Default)]
struct BucketIndex {
    by_chrom: HashMap<String, BTreeMap<i64, Vec<usize>>>,
}

impl BucketIndex {
    fn insert(&mut self, index: usize, chrom: &str, start: i64, end: i64) {
        if end <= start {
            return;
        }
        let map = self.by_chrom.entry(chrom.to_string()).or_default();
        let mut bucket = start.div_euclid(BUCKET) * BUCKET;
        let last = (end - 1).div_euclid(BUCKET) * BUCKET;
        while bucket <= last {
            map.entry(bucket).or_default().push(index);
            bucket += BUCKET;
        }
    }

    fn overlap(&self, chrom: &str, start: i64, end: i64) -> Vec<usize> {
        if end <= start {
            return Vec::new();
        }
        let Some(map) = self.by_chrom.get(chrom) else {
            return Vec::new();
        };
        let mut bucket = start.div_euclid(BUCKET) * BUCKET;
        let last = (end - 1).div_euclid(BUCKET) * BUCKET;
        let mut found = Vec::new();
        while bucket <= last {
            if let Some(list) = map.get(&bucket) {
                found.extend(list.iter().copied());
            }
            bucket += BUCKET;
        }
        found.sort_unstable();
        found.dedup();
        found
    }
}

pub fn quantify(settings: &ScotchSettings) -> Result<(), ScotchFailure> {
    let _ = settings.workers;
    validate(settings)?;
    let mut annotation = load_gtf(&settings.gtf).map_err(|message| fail(&message, false))?;
    if !settings.genes.is_empty() {
        retain_genes(&mut annotation, &settings.genes).map_err(|message| fail(&message, true))?;
    }
    let samples = expand_samples(&settings.bams)?;
    let genome = load_genome(settings.fasta.as_deref())?;
    let (cell_tag, umi_tag) = barcode_tags(settings);
    let mut coverage = if settings.update_gtf {
        annotation
            .genes
            .iter()
            .map(|gene| Coverage::new(gene.start, gene.end))
            .collect()
    } else {
        Vec::new()
    };
    let mut winners = Vec::new();
    let mut aligned = 0usize;
    let mut tagged = 0usize;
    for sample in &samples {
        let (winner, sample_aligned, sample_tagged) = scan_sample(
            sample,
            settings,
            &annotation,
            cell_tag,
            umi_tag,
            coverage.as_mut_slice(),
        )?;
        aligned += sample_aligned;
        tagged += sample_tagged;
        winners.push(winner);
    }
    if aligned == 0 {
        return Err(fail("no aligned reads in the BAM or SAM input", false));
    }
    if settings.platform != Platform::Bulk && tagged == 0 {
        return Err(fail(
            &format!("no reads carried barcode tags {cell_tag} and {umi_tag}"),
            false,
        ));
    }
    if settings.update_gtf {
        for (gene, depth) in annotation.genes.iter_mut().zip(&coverage) {
            refine_gene(
                gene,
                depth,
                settings.exon_fraction,
                settings.splice_fraction,
                settings.z_score,
            );
        }
    }
    let (metas, meta_index) = index_metagenes(&annotation);
    let mut hits = Vec::new();
    let ctx = AssignCtx {
        settings,
        annotation: &annotation,
        genome: genome.as_ref(),
        metas: &metas,
        meta_index: &meta_index,
        cell_tag,
        umi_tag,
    };
    for (index, sample) in samples.iter().enumerate() {
        assign_sample(sample, index, &winners[index], &ctx, &mut hits)?;
    }
    let (counts, novels, assignments) = summarize(&hits, &annotation, settings);
    write_outputs(
        settings,
        &samples,
        &annotation,
        &counts,
        &novels,
        &assignments,
    )
}

fn validate(settings: &ScotchSettings) -> Result<(), ScotchFailure> {
    if settings.bams.is_empty() {
        return Err(fail(
            "provide at least one BAM or SAM path with --bam",
            true,
        ));
    }
    let params = &settings.params;
    if !(0.0..=1.0).contains(&params.low)
        || !(0.0..=1.0).contains(&params.high)
        || params.high < params.low
        || !(0.0..=1.0).contains(&params.truncation_match)
        || params.small_exon < 0
        || params.small_exon_high < params.small_exon
    {
        return Err(fail(
            "match thresholds must sit in [0, 1], with match-high at least match-low",
            true,
        ));
    }
    if !(0.0..=1.0).contains(&settings.novel_read_pct)
        || settings.exon_fraction < 0.0
        || settings.splice_fraction < 0.0
        || settings.z_score < 0.0
    {
        return Err(fail(
            "novel fraction must sit in [0, 1] and coverage thresholds must be non-negative",
            true,
        ));
    }
    Ok(())
}

fn fail(message: &str, usage: bool) -> ScotchFailure {
    ScotchFailure {
        message: message.to_string(),
        usage,
    }
}

fn barcode_tags(settings: &ScotchSettings) -> (&str, &str) {
    let cell = settings.barcode_cell.as_deref().unwrap_or("CB");
    let umi = settings
        .barcode_umi
        .as_deref()
        .unwrap_or(match settings.platform {
            Platform::TenxPacbio => "XM",
            _ => "UB",
        });
    (cell, umi)
}

fn load_genome(path: Option<&Path>) -> Result<Option<HashMap<String, String>>, ScotchFailure> {
    let Some(path) = path else {
        return Ok(None);
    };
    let records = read_fasta(path).map_err(|message| fail(&message, false))?;
    let mut genome = HashMap::new();
    for record in records {
        genome.insert(record.id, record.seq);
    }
    Ok(Some(genome))
}

fn expand_samples(paths: &[PathBuf]) -> Result<Vec<Sample>, ScotchFailure> {
    let mut used = BTreeSet::new();
    let mut samples = Vec::new();
    for path in paths {
        if path.is_dir() {
            let mut files = Vec::new();
            let entries = std::fs::read_dir(path)
                .map_err(|err| fail(&format!("read {}: {err}", path.display()), false))?;
            for entry in entries {
                let entry =
                    entry.map_err(|err| fail(&format!("read {}: {err}", path.display()), false))?;
                let file = entry.path();
                if file.is_file() && is_alignment(&file) {
                    files.push(file);
                }
            }
            files.sort();
            if files.is_empty() {
                return Err(fail(
                    &format!("{} contains no BAM or SAM files", path.display()),
                    true,
                ));
            }
            let stem = path
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("sample")
                .to_string();
            samples.push(Sample {
                name: unique_name(&mut used, stem),
                files,
            });
        } else if path.is_file() {
            let stem = path
                .file_stem()
                .and_then(|name| name.to_str())
                .unwrap_or("sample")
                .to_string();
            samples.push(Sample {
                name: unique_name(&mut used, stem),
                files: vec![path.clone()],
            });
        } else {
            return Err(fail(
                &format!("alignment path {} does not exist", path.display()),
                true,
            ));
        }
    }
    Ok(samples)
}

fn is_alignment(path: &Path) -> bool {
    path.extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("bam") || ext.eq_ignore_ascii_case("sam"))
}

fn unique_name(used: &mut BTreeSet<String>, stem: String) -> String {
    if used.insert(stem.clone()) {
        return stem;
    }
    let mut suffix = 2usize;
    loop {
        let name = format!("{stem}_{suffix}");
        if used.insert(name.clone()) {
            return name;
        }
        suffix += 1;
    }
}

fn scan_sample(
    sample: &Sample,
    settings: &ScotchSettings,
    annotation: &Annotation,
    cell_tag: &str,
    umi_tag: &str,
    coverage: &mut [Coverage],
) -> Result<(HashMap<String, Winner>, usize, usize), ScotchFailure> {
    let mut winners = HashMap::new();
    let mut ordinal = 0u64;
    let mut aligned = 0usize;
    let mut tagged = 0usize;
    let gene_index = gene_buckets(annotation);
    for path in &sample.files {
        for read in open_alignments(path).map_err(|message| fail(&message, false))? {
            let read = read.map_err(|message| fail(&message, false))?;
            if is_unmapped(&read) {
                continue;
            }
            aligned += 1;
            if settings.update_gtf {
                observe_read(&read, annotation, &gene_index, coverage);
            }
            let Some(identity) = identify(&read, path, &sample.name, settings, cell_tag, umi_tag)
            else {
                continue;
            };
            tagged += 1;
            let winner = winners.entry(identity.key).or_insert(Winner {
                length: -1,
                ordinal: 0,
            });
            if identity.length > winner.length {
                *winner = Winner {
                    length: identity.length,
                    ordinal,
                };
            }
            ordinal += 1;
        }
    }
    Ok((winners, aligned, tagged))
}

fn gene_buckets(annotation: &Annotation) -> BucketIndex {
    let mut index = BucketIndex::default();
    for (gene_index, gene) in annotation.genes.iter().enumerate() {
        index.insert(gene_index, &gene.chrom, gene.start, gene.end);
    }
    index
}

fn observe_read(
    read: &Read,
    annotation: &Annotation,
    index: &BucketIndex,
    coverage: &mut [Coverage],
) {
    let blocks = read.blocks();
    let junctions = read.junctions();
    for gene_index in index.overlap(&read.chrom, read.start, read.end) {
        let gene = &annotation.genes[gene_index];
        if gene.start <= read.start && read.end < gene.end {
            coverage[gene_index].observe(read.start, read.end, &blocks, &junctions);
        }
    }
}

fn identify(
    read: &Read,
    path: &Path,
    sample: &str,
    settings: &ScotchSettings,
    cell_tag: &str,
    umi_tag: &str,
) -> Option<Identity> {
    let length = read.query_alignment_length();
    match settings.platform {
        Platform::Bulk => Some(Identity {
            key: read.qname.clone(),
            cell: sample.replace(':', "."),
            length,
        }),
        Platform::ParseOnt => {
            let parts: Vec<&str> = read.qname.split('_').collect();
            if parts.len() < 5 {
                return None;
            }
            let n = parts.len();
            let cell = format!("{}_{}_{}", parts[n - 5], parts[n - 4], parts[n - 3]);
            let sublib = sublibrary(path);
            Some(Identity {
                key: format!("{cell}_{}_{sublib}", parts[n - 1]),
                cell,
                length,
            })
        }
        Platform::TenxOnt | Platform::TenxPacbio => {
            let cell = read.tag(cell_tag)?.to_string();
            let umi = read.tag(umi_tag)?;
            Some(Identity {
                key: format!("{cell}_{umi}"),
                cell,
                length,
            })
        }
    }
}

fn sublibrary(path: &Path) -> String {
    let name = path
        .file_name()
        .and_then(|text| text.to_str())
        .unwrap_or("");
    let Some(rest) = name.split("sublibrary").nth(1) else {
        return "1".to_string();
    };
    let digits: String = rest.chars().take_while(|ch| ch.is_ascii_digit()).collect();
    if digits.is_empty() {
        "1".to_string()
    } else {
        digits
    }
}

fn index_metagenes(annotation: &Annotation) -> (Vec<Meta>, BucketIndex) {
    let mut metas = Vec::new();
    let mut index = BucketIndex::default();
    for genes in &annotation.metagenes {
        if genes.is_empty() {
            continue;
        }
        let chrom = annotation.genes[genes[0]].chrom.clone();
        let start = genes
            .iter()
            .map(|&gene| annotation.genes[gene].start)
            .min()
            .unwrap();
        let end = genes
            .iter()
            .map(|&gene| annotation.genes[gene].end)
            .max()
            .unwrap();
        index.insert(metas.len(), &chrom, start, end);
        metas.push(Meta {
            start,
            end,
            genes: genes.clone(),
        });
    }
    (metas, index)
}

struct AssignCtx<'a> {
    settings: &'a ScotchSettings,
    annotation: &'a Annotation,
    genome: Option<&'a HashMap<String, String>>,
    metas: &'a [Meta],
    meta_index: &'a BucketIndex,
    cell_tag: &'a str,
    umi_tag: &'a str,
}

fn assign_sample(
    sample: &Sample,
    sample_index: usize,
    winners: &HashMap<String, Winner>,
    ctx: &AssignCtx<'_>,
    hits: &mut Vec<Hit>,
) -> Result<(), ScotchFailure> {
    let mut ordinal = 0u64;
    for path in &sample.files {
        for read in open_alignments(path).map_err(|message| fail(&message, false))? {
            let read = read.map_err(|message| fail(&message, false))?;
            if is_unmapped(&read) {
                continue;
            }
            let Some(identity) = identify(
                &read,
                path,
                &sample.name,
                ctx.settings,
                ctx.cell_tag,
                ctx.umi_tag,
            ) else {
                continue;
            };
            let keep = winners
                .get(&identity.key)
                .is_some_and(|winner| winner.ordinal == ordinal);
            ordinal += 1;
            if !keep {
                continue;
            }
            if let Some(place) = place_read(&read, ctx) {
                hits.push(Hit {
                    sample: sample_index,
                    qname: read.qname,
                    cell: identity.cell,
                    gene: place.gene,
                    outcome: place.outcome,
                    score: place.score,
                });
            }
        }
    }
    Ok(())
}

fn place_read(read: &Read, ctx: &AssignCtx<'_>) -> Option<Placement> {
    let positions = read.ref_positions();
    let mut best: Option<Placement> = None;
    for meta_index_value in ctx.meta_index.overlap(&read.chrom, read.start, read.end) {
        let meta = &ctx.metas[meta_index_value];
        if meta.end <= read.start || meta.start >= read.end {
            continue;
        }
        let polys: Vec<bool> = meta
            .genes
            .iter()
            .map(|&gene| poly_for(read, &ctx.annotation.genes[gene], ctx.settings, ctx.genome))
            .collect();
        let place = place_in_metagene(
            &positions,
            read.start,
            read.end,
            &ctx.annotation.genes,
            &meta.genes,
            &ctx.settings.params,
            &polys,
        );
        if let Some(place) = place {
            best = Some(match best {
                Some(current) if prefer(&current, &place, &ctx.annotation.genes) => current,
                _ => place,
            });
        }
    }
    best
}

fn prefer(current: &Placement, other: &Placement, genes: &[Gene]) -> bool {
    if current.priority != other.priority {
        return current.priority > other.priority;
    }
    let current_len = genes[current.gene].end - genes[current.gene].start;
    let other_len = genes[other.gene].end - genes[other.gene].start;
    if current_len != other_len {
        return current_len > other_len;
    }
    genes[current.gene].id <= genes[other.gene].id
}

fn poly_for(
    read: &Read,
    gene: &Gene,
    settings: &ScotchSettings,
    genome: Option<&HashMap<String, String>>,
) -> bool {
    if settings.platform == Platform::TenxPacbio {
        return genome
            .and_then(|genome| genomic_tail(read, gene, genome))
            .unwrap_or(true);
    }
    let detected = match settings.platform {
        Platform::ParseOnt | Platform::Bulk => {
            detect_edge_poly(&read.seq, read.query_start, read.query_end)
        }
        Platform::TenxOnt => {
            let umi = read.tag("UR").or_else(|| {
                settings
                    .barcode_umi
                    .as_deref()
                    .or(Some("UB"))
                    .and_then(|tag| read.tag(tag))
            });
            query_poly(&read.seq, read.query_start, read.query_end, umi)
        }
        Platform::TenxPacbio => true,
    };
    if detected {
        if let Some(genome) = genome {
            return genomic_tail(read, gene, genome).unwrap_or(true);
        }
    }
    detected
}

fn genomic_tail(read: &Read, gene: &Gene, genome: &HashMap<String, String>) -> Option<bool> {
    let sequence = genome.get(&gene.chrom)?;
    let (pos, base) = if gene.strand == '+' {
        (read.end, b'A')
    } else {
        (read.start, b'T')
    };
    let start = pos.saturating_sub(15);
    let end = pos.saturating_add(15);
    Some(homopolymer(sequence, start, end, base) < 8)
}

type CountMap = BTreeMap<(usize, String, String, String), u64>;
type NovelMap = BTreeMap<(usize, String), Vec<usize>>;

struct SavedAssignment {
    sample: usize,
    qname: String,
    cell: String,
    gene: usize,
    isoform: String,
    score: f64,
}

fn summarize(
    hits: &[Hit],
    annotation: &Annotation,
    settings: &ScotchSettings,
) -> (CountMap, NovelMap, Vec<SavedAssignment>) {
    let mut isoform_name: Vec<Option<String>> = vec![None; hits.len()];
    let mut novel_exons: HashMap<(usize, String), Vec<usize>> = HashMap::new();
    let mut by_gene: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (index, hit) in hits.iter().enumerate() {
        by_gene.entry(hit.gene).or_default().push(index);
    }
    for (gene_index, indexes) in &by_gene {
        let gene = &annotation.genes[*gene_index];
        let novel_hits: Vec<usize> = indexes
            .iter()
            .copied()
            .filter(|index| matches!(hits[*index].outcome, Outcome::Novel(_)))
            .collect();
        let reads: Vec<(String, Vec<i8>)> = novel_hits
            .iter()
            .filter_map(|index| match &hits[*index].outcome {
                Outcome::Novel(vector) => Some((hits[*index].qname.clone(), vector.clone())),
                _ => None,
            })
            .collect();
        let (isoforms, calls) = discover(&reads, NOVEL_MIN_READS, NOVEL_MAX_ROUNDS);
        let mut exons_of: HashMap<String, Vec<usize>> = isoforms
            .iter()
            .map(|isoform| (isoform.name.clone(), isoform.exons.clone()))
            .collect();
        let strand = if settings.platform == Platform::ParseOnt {
            '.'
        } else {
            gene.strand
        };
        let grouped = if settings.group_novel {
            group_novel(
                &isoforms
                    .iter()
                    .map(|isoform| (isoform.name.clone(), isoform.exons.clone()))
                    .collect::<Vec<_>>(),
                strand,
            )
        } else {
            Vec::new()
        };
        let rename = |name: &str| {
            grouped
                .iter()
                .find(|(from, _)| from == name)
                .map(|(_, to)| to.clone())
                .unwrap_or_else(|| name.to_string())
        };
        let mut support: BTreeMap<String, usize> = BTreeMap::new();
        let mut known = 0usize;
        let mut assigned = Vec::new();
        for call in &calls {
            if let NovelCall::Isoform(name) = call {
                let name = rename(name);
                *support.entry(name.clone()).or_default() += 1;
                assigned.push(Some(name));
            } else {
                assigned.push(None);
            }
        }
        for index in indexes {
            if matches!(hits[*index].outcome, Outcome::Known(_)) {
                known += 1;
            }
        }
        let novel_total: usize = support.values().sum();
        let denom = known + novel_total;
        let mut dropped = BTreeSet::new();
        for (name, count) in &support {
            let fraction = if denom == 0 {
                0.0
            } else {
                *count as f64 / denom as f64
            };
            if *count < settings.novel_read_n
                || (settings.novel_read_pct > 0.0 && fraction < settings.novel_read_pct)
            {
                dropped.insert(name.clone());
            }
        }
        for (slot, hit_index) in novel_hits.iter().enumerate() {
            isoform_name[*hit_index] = match assigned.get(slot).and_then(Option::as_ref) {
                Some(name) if dropped.contains(name) => Some("uncategorized_novel".to_string()),
                Some(name) => {
                    if let Some(exons) = exons_of.remove(name) {
                        novel_exons.insert((*gene_index, name.clone()), exons);
                    }
                    Some(name.clone())
                }
                None => Some("uncategorized".to_string()),
            };
        }
    }
    let mut counts = CountMap::new();
    let mut assignments = Vec::new();
    let mut novels = NovelMap::new();
    for (index, hit) in hits.iter().enumerate() {
        let gene = &annotation.genes[hit.gene];
        let isoform = match &hit.outcome {
            Outcome::Known(isoform) => gene.isoforms[*isoform].name.clone(),
            Outcome::Novel(_) => isoform_name[index]
                .clone()
                .unwrap_or_else(|| "uncategorized".to_string()),
            Outcome::Uncategorized => "uncategorized".to_string(),
        };
        if let Some(exons) = novel_exons.get(&(hit.gene, isoform.clone())) {
            novels
                .entry((hit.gene, isoform.clone()))
                .or_insert_with(|| exons.clone());
        }
        *counts
            .entry((
                hit.sample,
                hit.cell.clone(),
                gene.label.clone(),
                isoform.clone(),
            ))
            .or_default() += 1;
        assignments.push(SavedAssignment {
            sample: hit.sample,
            qname: hit.qname.clone(),
            cell: hit.cell.clone(),
            gene: hit.gene,
            isoform,
            score: hit.score,
        });
    }
    (counts, novels, assignments)
}

fn write_outputs(
    settings: &ScotchSettings,
    samples: &[Sample],
    annotation: &Annotation,
    counts: &CountMap,
    novels: &NovelMap,
    assignments: &[SavedAssignment],
) -> Result<(), ScotchFailure> {
    std::fs::create_dir_all(&settings.out)
        .map_err(|err| fail(&format!("create {}: {err}", settings.out.display()), false))?;
    let mut gtf = annotation.gtf_text.clone();
    if !gtf.is_empty() && !gtf.ends_with('\n') {
        gtf.push('\n');
    }
    for ((gene_index, name), exons) in novels {
        append_novel(&mut gtf, &annotation.genes[*gene_index], name, exons);
    }
    let annotation_path = settings.out.join("annotation.gtf");
    std::fs::write(&annotation_path, gtf).map_err(|err| {
        fail(
            &format!("write {}: {err}", annotation_path.display()),
            false,
        )
    })?;
    for (sample_index, sample) in samples.iter().enumerate() {
        write_sample(
            settings,
            sample_index,
            sample,
            annotation,
            counts,
            assignments,
        )?;
    }
    Ok(())
}

fn append_novel(text: &mut String, gene: &Gene, name: &str, exons: &[usize]) {
    let pieces = gtf_exons(gene, exons);
    if pieces.is_empty() {
        return;
    }
    let start = pieces.iter().map(|piece| piece.0).min().unwrap();
    let end = pieces.iter().map(|piece| piece.1).max().unwrap();
    let transcript = format!(
        "gene_id \"{}\"; transcript_id \"{name}\"; gene_name \"{}\";",
        gene.id, gene.name
    );
    text.push_str(&format!(
        "{}\tSCOTCH\ttranscript\t{start}\t{end}\t.\t{}\t.\t{transcript}\n",
        gene.chrom, gene.strand
    ));
    for (exon_start, exon_end) in pieces {
        text.push_str(&format!(
            "{}\tSCOTCH\texon\t{exon_start}\t{exon_end}\t.\t{}\t.\t{transcript}\n",
            gene.chrom, gene.strand
        ));
    }
}

fn write_sample(
    settings: &ScotchSettings,
    sample_index: usize,
    sample: &Sample,
    annotation: &Annotation,
    counts: &CountMap,
    assignments: &[SavedAssignment],
) -> Result<(), ScotchFailure> {
    let count_dir = settings.out.join(&sample.name).join("count_matrix");
    let aux_dir = settings.out.join(&sample.name).join("auxiliary");
    std::fs::create_dir_all(&count_dir)
        .and_then(|_| std::fs::create_dir_all(&aux_dir))
        .map_err(|err| fail(&format!("create {}: {err}", count_dir.display()), false))?;
    let mut cells = BTreeSet::new();
    let mut isoforms: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for (sample_id, cell, gene, isoform) in counts.keys() {
        if *sample_id == sample_index {
            cells.insert(format!("{cell}:sample{sample_index}"));
            isoforms
                .entry(gene.clone())
                .or_default()
                .insert(isoform.clone());
        }
    }
    let mut gene_header = vec!["cell".to_string()];
    let mut tx_header = vec!["cell".to_string()];
    let mut map = String::from("genes\ttranscripts\n");
    let mut gene_names: Vec<String> = isoforms.keys().cloned().collect();
    gene_names.sort();
    for gene in &gene_names {
        gene_header.push(gene.clone());
        let mut names: Vec<&String> = isoforms[gene].iter().collect();
        names.sort();
        for isoform in names {
            let column = format!("{gene}_{isoform}");
            tx_header.push(column.clone());
            map.push_str(gene);
            map.push('\t');
            map.push_str(&column);
            map.push('\n');
        }
    }
    let mut gene_csv = csv_line(&gene_header);
    let mut tx_csv = csv_line(&tx_header);
    for cell in &cells {
        let barcode = cell
            .rsplit_once(":sample")
            .map(|(barcode, _)| barcode)
            .unwrap_or(cell.as_str());
        let mut gene_row = vec![cell.clone()];
        let mut tx_row = vec![cell.clone()];
        for gene in &gene_names {
            let mut gene_count = 0u64;
            let mut names: Vec<&String> = isoforms[gene].iter().collect();
            names.sort();
            for isoform in names {
                let count = counts
                    .get(&(
                        sample_index,
                        barcode.to_string(),
                        gene.clone(),
                        isoform.clone(),
                    ))
                    .copied()
                    .unwrap_or(0);
                gene_count += count;
                tx_row.push(count.to_string());
            }
            gene_row.push(gene_count.to_string());
        }
        gene_csv.push_str(&csv_line(&gene_row));
        tx_csv.push_str(&csv_line(&tx_row));
    }
    std::fs::write(count_dir.join("gene_counts.csv"), gene_csv)
        .and_then(|_| std::fs::write(count_dir.join("transcript_counts.csv"), tx_csv))
        .and_then(|_| std::fs::write(count_dir.join("gene_transcript.tsv"), map))
        .map_err(|err| fail(&format!("write counts for {}: {err}", sample.name), false))?;
    let mut assignment_text = String::from("sample\tqname\tcell\tgene\tgene_id\tisoform\tscore\n");
    for SavedAssignment {
        sample: sample_id,
        qname,
        cell,
        gene,
        isoform,
        score,
    } in assignments
    {
        if *sample_id != sample_index {
            continue;
        }
        assignment_text.push_str(&format!(
            "{}\t{qname}\t{cell}:sample{sample_index}\t{}\t{}\t{isoform}\t{score:.6}\n",
            sample.name, annotation.genes[*gene].label, annotation.genes[*gene].id
        ));
    }
    let assignment_path = aux_dir.join("assignments.tsv");
    std::fs::write(&assignment_path, assignment_text).map_err(|err| {
        fail(
            &format!("write {}: {err}", assignment_path.display()),
            false,
        )
    })
}

fn csv_line(fields: &[String]) -> String {
    let mut line = String::new();
    for (index, field) in fields.iter().enumerate() {
        if index > 0 {
            line.push(',');
        }
        if field.contains([',', '"', '\n', '\r']) {
            line.push('"');
            line.push_str(&field.replace('"', "\"\""));
            line.push('"');
        } else {
            line.push_str(field);
        }
    }
    line.push('\n');
    line
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seq(bases: usize) -> String {
        "ACGT".repeat(bases / 4)
    }

    fn record(qname: &str, cigar: &str, pos: i64, cb: &str, ub: &str) -> String {
        let length = cigar_query_len(cigar);
        format!(
            "{qname}\t0\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{}\t*\tCB:Z:{cb}\tUB:Z:{ub}\n",
            seq(length)
        )
    }

    fn cigar_query_len(cigar: &str) -> usize {
        let mut number = String::new();
        let mut total = 0usize;
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                number.push(ch);
                continue;
            }
            let len: usize = number.parse().unwrap_or(0);
            number.clear();
            if matches!(ch, 'M' | 'I' | 'S' | '=' | 'X') {
                total += len;
            }
        }
        total
    }

    #[test]
    fn quantifies_known_isoforms_and_a_ten_read_novel() {
        let root = std::env::temp_dir().join(format!("plena-scotch-quant-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(&root).unwrap();
        let gtf = root.join("genes.gtf");
        std::fs::write(
            &gtf,
            "\
chr1\t.\tgene\t101\t600\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\";
chr1\t.\ttranscript\t101\t600\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ALL\";
chr1\t.\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ALL\";
chr1\t.\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ALL\";
chr1\t.\texon\t501\t600\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ALL\";
chr1\t.\ttranscript\t101\t400\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"HEAD\";
chr1\t.\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"HEAD\";
chr1\t.\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"HEAD\";
",
        )
        .unwrap();
        let mut sam = String::new();
        sam.push_str(&record("SHORT", "100M", 101, "FULL", "U1"));
        sam.push_str(&record("FULL", "100M100N100M100N100M", 101, "FULL", "U1"));
        sam.push_str(&record("HEAD", "100M100N100M", 101, "HEAD", "U1"));
        for index in 0..10 {
            sam.push_str(&record(
                &format!("N{index}"),
                "100M300N100M",
                101,
                &format!("N{index}"),
                "U",
            ));
        }
        let sam_path = root.join("cells.sam");
        std::fs::write(&sam_path, sam).unwrap();
        let out = root.join("out");
        quantify(&ScotchSettings {
            bams: vec![sam_path],
            gtf,
            out: out.clone(),
            ..ScotchSettings::default()
        })
        .unwrap();
        let genes =
            std::fs::read_to_string(out.join("cells/count_matrix/gene_counts.csv")).unwrap();
        let transcripts =
            std::fs::read_to_string(out.join("cells/count_matrix/transcript_counts.csv")).unwrap();
        let annotation = std::fs::read_to_string(out.join("annotation.gtf")).unwrap();
        let assigned =
            std::fs::read_to_string(out.join("cells/auxiliary/assignments.tsv")).unwrap();
        assert!(genes.contains("FULL:sample0"));
        assert!(genes.contains("HEAD:sample0"));
        assert!(transcripts.contains("GENE1_ALL"));
        assert!(transcripts.contains("GENE1_HEAD"));
        assert!(transcripts.contains("GENE1_novelIsoform_5"));
        assert!(annotation.contains("transcript_id \"novelIsoform_5\""));
        assert!(!assigned.lines().any(|line| line.contains("\tSHORT\t")));
        assert_eq!(
            assigned
                .lines()
                .filter(|line| line.split('\t').nth(1) == Some("FULL"))
                .count(),
            1
        );
        let full = transcripts
            .lines()
            .find(|line| line.starts_with("FULL:sample0"))
            .unwrap();
        let header: Vec<&str> = transcripts.lines().next().unwrap().split(',').collect();
        let all = header.iter().position(|name| *name == "GENE1_ALL").unwrap();
        let novel = header
            .iter()
            .position(|name| *name == "GENE1_novelIsoform_5")
            .unwrap();
        let fields: Vec<&str> = full.split(',').collect();
        assert_eq!(fields[all], "1");
        assert_eq!(fields[novel], "0");
        let _ = std::fs::remove_dir_all(&root);
    }
}
