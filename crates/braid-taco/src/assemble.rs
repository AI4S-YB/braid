//! Aggregate sample GTFs and assemble one transcriptome.
//!
//! Sample ids, locus ids, and gene ids follow input order. Neighbor iteration
//! elsewhere is sorted, so this does not reproduce one Python 2.7 hash seed.
//! The isoform sequence matches TACO when the winning bottleneck is unique.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::path::{Path, PathBuf};

use braid_model::{Exon, Strand as ModelStrand, Transcript};

use crate::path::{assemble_paths, assign_ids, build_isoforms};
use crate::splice::{gtf_line, SpliceGraph};
use crate::types::{Strand, Transfrag};
use crate::util::py_float;

#[derive(Debug)]
pub struct TacoFailure {
    pub message: String,
    pub usage: bool,
}

impl std::fmt::Display for TacoFailure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.message)
    }
}

fn usage(message: impl Into<String>) -> TacoFailure {
    TacoFailure {
        message: message.into(),
        usage: true,
    }
}

fn runtime(message: impl Into<String>) -> TacoFailure {
    TacoFailure {
        message: message.into(),
        usage: false,
    }
}

#[derive(Clone, Debug)]
pub struct TacoSettings {
    pub sample_file: PathBuf,
    pub gtf_expr_attr: String,
    pub filter_min_length: i64,
    pub filter_min_expr: f64,
    pub isoform_frac: f64,
    pub max_isoforms: i64,
    pub assemble_unstranded: bool,
    pub change_point: bool,
    pub change_point_pvalue: f64,
    pub change_point_fold_change: f64,
    pub change_point_trim: bool,
    pub path_kmax: i64,
    pub path_frac: f64,
    pub max_paths: i64,
    pub ref_gtf: Option<PathBuf>,
    pub guided_strand: bool,
    pub guided_ends: bool,
    pub guided_assembly: bool,
    pub filter_splice_juncs: bool,
    pub splice_motifs: Vec<String>,
    pub genome_fasta: Option<PathBuf>,
}

impl Default for TacoSettings {
    fn default() -> Self {
        Self {
            sample_file: PathBuf::new(),
            gtf_expr_attr: "FPKM".to_string(),
            filter_min_length: 200,
            filter_min_expr: 0.5,
            isoform_frac: 0.05,
            max_isoforms: 0,
            assemble_unstranded: false,
            change_point: true,
            change_point_pvalue: 0.01,
            change_point_fold_change: 0.85,
            change_point_trim: true,
            path_kmax: 0,
            path_frac: 0.0,
            max_paths: 0,
            ref_gtf: None,
            guided_strand: false,
            guided_ends: false,
            guided_assembly: false,
            filter_splice_juncs: false,
            splice_motifs: Vec::new(),
            genome_fasta: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct TacoTexts {
    pub transcripts: Vec<Transcript>,
    pub samples: String,
    pub sample_stats: String,
    pub transfrags_bed: String,
    pub transfrags_filtered_bed: String,
    pub loci: String,
    pub bedgraph_pos: String,
    pub bedgraph_neg: String,
    pub bedgraph_none: String,
    pub splice_junctions: String,
    pub splice_graph: String,
    pub change_points: String,
    pub path_graph_stats: String,
    pub assembly_gtf: String,
    pub assembly_bed: String,
}

struct Sample {
    path: String,
    id: String,
    is_ref: bool,
}

struct SampleStats {
    id: String,
    total: usize,
    length: usize,
    expr: usize,
    splice: usize,
}

pub fn assemble(settings: &TacoSettings) -> Result<TacoTexts, TacoFailure> {
    validate(settings)?;
    let samples = read_samples(&settings.sample_file)?;
    let mut work = Vec::new();
    if let Some(reference) = &settings.ref_gtf {
        let id = "R".to_string();
        if samples.iter().any(|sample| sample.id == id) {
            return Err(runtime(format!("sample_id '{id}' is not unique")));
        }
        work.push(Sample {
            path: reference.display().to_string(),
            id,
            is_ref: true,
        });
    }
    work.extend(samples);

    let genome = if settings.filter_splice_juncs {
        let path = settings.genome_fasta.as_ref().unwrap();
        Some(load_fasta(path).map_err(runtime)?)
    } else {
        None
    };
    let motifs = allowed_motifs(&settings.splice_motifs);
    aggregate(settings, &work, genome.as_ref(), &motifs)
}

#[derive(Default)]
struct AssemblyIds {
    gene: u64,
    tss: u64,
    transcript: u64,
}

#[derive(Default)]
struct AssemblyOutput {
    splice_lines: Vec<String>,
    change_lines: Vec<String>,
    path_lines: Vec<String>,
    assembly_gtf: Vec<String>,
    assembly_bed: Vec<String>,
    transcripts: Vec<Transcript>,
}

fn aggregate(
    settings: &TacoSettings,
    work: &[Sample],
    genome: Option<&HashMap<String, String>>,
    motifs: &HashSet<String>,
) -> Result<TacoTexts, TacoFailure> {
    let mut kept: Vec<Transfrag> = Vec::new();
    let mut filtered_lines = Vec::new();
    let mut stats = Vec::new();
    for sample in work {
        let (transfrags, total_expr) = parse_gtf(
            Path::new(&sample.path),
            &sample.id,
            &settings.gtf_expr_attr,
            sample.is_ref,
        )
        .map_err(runtime)?;
        let total_count = transfrags.len();
        let mut n_length = 0usize;
        let mut n_expr = 0usize;
        let mut n_splice = 0usize;
        for mut transfrag in transfrags {
            if total_expr > 0.0 {
                transfrag.expr = 1.0e6 * transfrag.expr / total_expr;
            }
            let mut keep = true;
            if transfrag.length() < settings.filter_min_length {
                keep = false;
                n_length += 1;
            }
            if transfrag.expr < settings.filter_min_expr {
                keep = false;
                n_expr += 1;
            }
            if settings.filter_splice_juncs {
                let bad = bad_introns(&transfrag, genome.unwrap(), motifs).map_err(runtime)?;
                if bad > 0 {
                    keep = false;
                    n_splice += bad;
                }
            }
            if keep {
                kept.push(transfrag);
            } else {
                filtered_lines.push(transfrag_bed(&transfrag));
            }
        }
        stats.push(SampleStats {
            id: sample.id.clone(),
            total: total_count,
            length: n_length,
            expr: n_expr,
            splice: n_splice,
        });
    }

    let mut order: Vec<usize> = (0..kept.len()).collect();
    order.sort_by(|&a, &b| {
        kept[a]
            .chrom
            .cmp(&kept[b].chrom)
            .then(kept[a].start().cmp(&kept[b].start()))
            .then(a.cmp(&b))
    });
    let kept: Vec<Transfrag> = order.into_iter().map(|index| kept[index].clone()).collect();

    let mut bed = String::new();
    let mut placed = Vec::new();
    for transfrag in &kept {
        let offset = bed.len();
        let line = transfrag_bed(transfrag);
        bed.push_str(&line);
        bed.push('\n');
        placed.push((offset, transfrag));
    }

    let loci = index_loci(&placed);
    let num_samples = work.iter().filter(|sample| !sample.is_ref).count();
    if num_samples == 0 {
        return Err(runtime("no input samples"));
    }

    let mut graphs = [Vec::new(), Vec::new(), Vec::new()];
    let mut junctions = Vec::new();
    let mut output = AssemblyOutput::default();
    let mut ids = AssemblyIds::default();

    for locus in &loci {
        let group: Vec<Transfrag> = placed[locus.first..locus.first + locus.count]
            .iter()
            .map(|(_, transfrag)| (*transfrag).clone())
            .collect();
        let mut locus_obj = Locus::create(
            group,
            settings.guided_strand,
            settings.guided_ends,
            settings.guided_assembly,
        )
        .map_err(runtime)?;
        locus_obj.impute();
        for (strand_index, lines) in locus_obj.bedgraph().into_iter().enumerate() {
            graphs[strand_index].extend(lines);
        }
        junctions.extend(locus_obj.junctions());
        for strand in Strand::ALL {
            if !locus_obj.strands.contains(&strand) {
                continue;
            }
            if strand == Strand::Na && !settings.assemble_unstranded {
                continue;
            }
            let frags = locus_obj.transfrags[strand as usize].clone();
            if frags.is_empty() {
                continue;
            }
            let graph = SpliceGraph::create(frags, settings.guided_ends, settings.guided_assembly)
                .map_err(runtime)?;
            for part in graph.split().map_err(runtime)? {
                assemble_gene(
                    part,
                    &locus.name,
                    settings,
                    num_samples,
                    &mut ids,
                    &mut output,
                )?;
            }
        }
    }

    stats.sort_by(|a, b| a.id.cmp(&b.id));
    let mut sample_stats = String::from(
        "sample_id\tnum_transfrags\tfiltered_length\tfiltered_expr\tfiltered_splice\n",
    );
    for stat in &stats {
        sample_stats.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\n",
            stat.id, stat.total, stat.length, stat.expr, stat.splice
        ));
    }
    let mut samples_txt = String::from("gtf\tsample_id\n");
    for sample in work.iter().filter(|sample| !sample.is_ref) {
        samples_txt.push_str(&format!("{}\t{}\n", sample.path, sample.id));
    }
    let mut loci_txt = String::new();
    for locus in &loci {
        loci_txt.push_str(&format!(
            "L{}\t{}\t{}\t{}\t{}\t{}\n",
            locus.id, locus.chrom, locus.start, locus.end, locus.filepos, locus.count
        ));
    }

    let AssemblyOutput {
        splice_lines,
        change_lines,
        path_lines,
        assembly_gtf,
        assembly_bed,
        transcripts,
    } = output;
    Ok(TacoTexts {
        transcripts,
        samples: samples_txt,
        sample_stats,
        transfrags_bed: bed,
        transfrags_filtered_bed: join_lines(filtered_lines),
        loci: loci_txt,
        bedgraph_pos: sort_bed(graphs[0].clone()),
        bedgraph_neg: sort_bed(graphs[1].clone()),
        bedgraph_none: sort_bed(graphs[2].clone()),
        splice_junctions: with_header(
            "track name=junctions description=\"Splice Junctions\" graphType=junctions",
            sort_bed(junctions),
        ),
        splice_graph: sort_gtf(splice_lines),
        change_points: sort_gtf(change_lines),
        path_graph_stats: with_header(
            "chrom\tstart\tend\tstrand\tk\tkmax\ttransfrags\tshort_transfrags\tshort_expr\tlost_short\tlost_short_expr\tkmers\tlost_kmers\ttot_expr\tgraph_expr\texpr_frac\tvalid\topt\tis_opt",
            sort_bed(path_lines),
        ),
        assembly_gtf: sort_gtf(assembly_gtf),
        assembly_bed: sort_bed(assembly_bed),
    })
}

fn assemble_gene(
    mut graph: SpliceGraph,
    locus_id: &str,
    settings: &TacoSettings,
    num_samples: usize,
    ids: &mut AssemblyIds,
    output: &mut AssemblyOutput,
) -> Result<(), TacoFailure> {
    let AssemblyIds {
        gene: gene_id,
        tss: tss_id,
        transcript: transcript_id,
    } = ids;
    let AssemblyOutput {
        splice_lines,
        change_lines,
        path_lines,
        assembly_gtf,
        assembly_bed,
        transcripts,
    } = output;
    if settings.change_point {
        let points = graph.detect_change_points(
            settings.change_point_pvalue,
            settings.change_point_fold_change,
        );
        let any = !points.is_empty();
        for point in &points {
            graph.apply_change_point(point, settings.change_point_trim);
            change_lines.extend(graph.change_gtf(point));
        }
        if any {
            graph.recreate().map_err(runtime)?;
        }
    }
    splice_lines.extend(graph.node_gtf());
    let (found, stats) = assemble_paths(
        &graph,
        settings.path_kmax.max(0) as usize,
        settings.path_frac,
        settings.max_paths.max(0) as usize,
    )
    .map_err(runtime)?;
    path_lines.extend(stats);
    let paths: Vec<(Vec<usize>, f64)> = found
        .into_iter()
        .map(|path| (path.nodes, path.expr))
        .collect();
    let mut genes = build_isoforms(&paths, settings.isoform_frac, |path| {
        graph.reconstruct_exons(path)
    });
    for isoforms in &mut genes {
        if settings.max_isoforms > 0 {
            isoforms.truncate(settings.max_isoforms as usize);
        }
        if isoforms.is_empty() {
            continue;
        }
        assign_ids(isoforms, graph.strand, gene_id, tss_id);
        for isoform in isoforms {
            *transcript_id += 1;
            let gene = format!("G{}", isoform.gene_id);
            let tss = format!("TSS{}", isoform.tss_id);
            let tid = format!("TU{transcript_id}");
            let expr = isoform.expr / num_samples as f64;
            let score = (1000.0 * isoform.rel_frac).round() as i64;
            let attrs = [
                ("expr", format!("{expr:.3}")),
                ("rel_frac", format!("{:.5}", isoform.rel_frac)),
                ("abs_frac", format!("{:.5}", isoform.abs_frac)),
                ("locus_id", locus_id.to_string()),
                ("gene_id", gene.clone()),
                ("tss_id", tss.clone()),
                ("transcript_id", tid.clone()),
            ];
            let tx_start = isoform.exons[0].0;
            let tx_end = isoform.exons.last().unwrap().1;
            assembly_gtf.push(gtf_line(
                &graph.chrom,
                "transcript",
                tx_start,
                tx_end,
                &score.to_string(),
                graph.strand.as_gtf(),
                &attrs,
            ));
            let exon_attrs = [
                ("locus_id", locus_id.to_string()),
                ("gene_id", gene.clone()),
                ("tss_id", tss),
                ("transcript_id", tid.clone()),
            ];
            for &(start, end) in &isoform.exons {
                assembly_gtf.push(gtf_line(
                    &graph.chrom,
                    "exon",
                    start,
                    end,
                    &score.to_string(),
                    graph.strand.as_gtf(),
                    &exon_attrs,
                ));
            }
            let name = format!("{gene}|{tid}({expr:.1})");
            assembly_bed.push(assembly_bed_line(
                &graph.chrom,
                &name,
                graph.strand,
                score,
                &isoform.exons,
            ));
            if let Some(strand) = model_strand(graph.strand) {
                transcripts.push(Transcript {
                    chrom: graph.chrom.clone(),
                    strand,
                    exons: isoform
                        .exons
                        .iter()
                        .map(|&(start, end)| Exon {
                            start: start + 1,
                            end: end + 1,
                        })
                        .collect(),
                    gene_id: gene,
                    transcript_id: tid,
                    source: Some("taco".to_string()),
                    score: Some(expr),
                    cds_start: None,
                    cds_end: None,
                });
            }
        }
    }
    Ok(())
}

fn model_strand(strand: Strand) -> Option<ModelStrand> {
    match strand {
        Strand::Pos => Some(ModelStrand::Forward),
        Strand::Neg => Some(ModelStrand::Reverse),
        Strand::Na => None,
    }
}

struct Locus {
    strands: HashSet<Strand>,
    expr: [Vec<f32>; 3],
    mask: [Vec<u8>; 3],
    transfrags: [Vec<Transfrag>; 3],
    guided_strand: bool,
    chrom: String,
    start: i64,
}

impl Locus {
    fn create(
        transfrags: Vec<Transfrag>,
        guided_strand: bool,
        guided_ends: bool,
        guided_assembly: bool,
    ) -> Result<Self, String> {
        let chrom = transfrags[0].chrom.clone();
        let start = transfrags.iter().map(Transfrag::start).min().unwrap();
        let end = transfrags.iter().map(Transfrag::end).max().unwrap();
        let len = (end - start) as usize;
        let mut locus = Self {
            strands: HashSet::new(),
            expr: [vec![0.0; len], vec![0.0; len], vec![0.0; len]],
            mask: [vec![0; len], vec![0; len], vec![0; len]],
            transfrags: [Vec::new(), Vec::new(), Vec::new()],
            guided_strand,
            chrom,
            start,
        };
        let _ = (guided_ends, guided_assembly);
        for transfrag in transfrags {
            if transfrag.chrom != locus.chrom {
                return Err("chrom mismatch".to_string());
            }
            locus.add(transfrag);
        }
        Ok(locus)
    }

    fn add(&mut self, transfrag: Transfrag) {
        let strand = transfrag.strand as usize;
        self.strands.insert(transfrag.strand);
        for &(start, end) in &transfrag.exons {
            let from = (start - self.start) as usize;
            let to = (end - self.start) as usize;
            if !(transfrag.is_ref && self.guided_strand) {
                for slot in self.expr[strand].iter_mut().take(to).skip(from) {
                    *slot += transfrag.expr as f32;
                }
            }
            for slot in self.mask[strand].iter_mut().take(to).skip(from) {
                *slot = 1;
            }
        }
        self.transfrags[strand].push(transfrag);
    }

    fn impute(&mut self) {
        let mut pending = std::mem::take(&mut self.transfrags[Strand::Na as usize]);
        let mut resolved_count = 0usize;
        loop {
            if pending.is_empty() {
                break;
            }
            let mut still = Vec::new();
            let mut resolved = Vec::new();
            for transfrag in pending.drain(..) {
                match self.predict(&transfrag) {
                    Strand::Na => still.push(transfrag),
                    strand => resolved.push((transfrag, strand)),
                }
            }
            if resolved.is_empty() {
                pending = still;
                break;
            }
            for (mut transfrag, strand) in resolved {
                transfrag.strand = strand;
                self.add(transfrag);
                resolved_count += 1;
            }
            pending = still;
        }
        if resolved_count > 0 {
            self.expr[Strand::Na as usize].fill(0.0);
            self.mask[Strand::Na as usize].fill(0);
            self.strands.remove(&Strand::Na);
            for transfrag in pending {
                self.add(transfrag);
            }
        } else {
            self.transfrags[Strand::Na as usize] = pending;
        }
    }

    fn predict(&self, transfrag: &Transfrag) -> Strand {
        let mut hit = [false, false];
        for &(start, end) in &transfrag.exons {
            let from = (start - self.start) as usize;
            let to = (end - self.start) as usize;
            for (strand, hit) in hit.iter_mut().enumerate() {
                let covered = self.mask[strand][from..to].iter().any(|bit| *bit != 0)
                    || self.expr[strand][from..to].iter().any(|value| *value > 0.0);
                if covered {
                    *hit = true;
                }
            }
        }
        match (hit[0], hit[1]) {
            (true, false) => Strand::Pos,
            (false, true) => Strand::Neg,
            _ => Strand::Na,
        }
    }

    fn bedgraph(&self) -> [Vec<String>; 3] {
        let mut out = [Vec::new(), Vec::new(), Vec::new()];
        for strand in Strand::ALL {
            if !self.strands.contains(&strand) {
                continue;
            }
            out[strand as usize] =
                bedgraph_array(&self.expr[strand as usize], &self.chrom, self.start);
        }
        out
    }

    fn junctions(&self) -> Vec<String> {
        let mut totals: std::collections::BTreeMap<(i64, i64, Strand), f64> =
            std::collections::BTreeMap::new();
        for strand in [Strand::Pos, Strand::Neg] {
            for transfrag in &self.transfrags[strand as usize] {
                if transfrag.is_ref {
                    continue;
                }
                for (start, end) in transfrag.introns() {
                    *totals.entry((start, end, strand)).or_default() += transfrag.expr;
                }
            }
        }
        totals
            .into_iter()
            .map(|((start, end, strand), expr)| {
                format!(
                    "{}\t{}\t{}\tJUNC\t{}\t{}\t{}\t{}\t255,0,0\t2\t1,1\t0,{}",
                    self.chrom,
                    start - 1,
                    end + 1,
                    py_float(expr),
                    strand.as_gtf(),
                    start - 1,
                    end + 1,
                    end + 1 - start
                )
            })
            .collect()
    }
}

struct LocusSpan {
    id: usize,
    name: String,
    chrom: String,
    start: i64,
    end: i64,
    filepos: usize,
    count: usize,
    first: usize,
}

fn index_loci(placed: &[(usize, &Transfrag)]) -> Vec<LocusSpan> {
    if placed.is_empty() {
        return Vec::new();
    }
    let mut loci = Vec::new();
    let mut id = 1usize;
    let mut first = 0usize;
    let mut chrom = placed[0].1.chrom.clone();
    let mut start = placed[0].1.start();
    let mut end = placed[0].1.end();
    let mut filepos = placed[0].0;
    let mut count = 1usize;
    for (index, (offset, transfrag)) in placed.iter().enumerate().skip(1) {
        if !(chrom == transfrag.chrom && start < transfrag.end() && transfrag.start() < end) {
            loci.push(LocusSpan {
                name: format!("L{id}"),
                id,
                chrom: chrom.clone(),
                start,
                end,
                filepos,
                count,
                first,
            });
            id += 1;
            first = index;
            chrom = transfrag.chrom.clone();
            start = transfrag.start();
            end = transfrag.end();
            filepos = *offset;
            count = 1;
        } else {
            if transfrag.end() > end {
                end = transfrag.end();
            }
            count += 1;
        }
    }
    loci.push(LocusSpan {
        name: format!("L{id}"),
        id,
        chrom,
        start,
        end,
        filepos,
        count,
        first,
    });
    loci
}

fn validate(settings: &TacoSettings) -> Result<(), TacoFailure> {
    if !settings.sample_file.exists() {
        return Err(usage(format!(
            "Sample file {} not found",
            settings.sample_file.display()
        )));
    }
    if settings.filter_min_length < 0 {
        return Err(usage("filter_min_length < 0"));
    }
    if settings.filter_min_expr < 0.0 {
        return Err(usage("filter_min_expr < 0"));
    }
    if !(0.0..=1.0).contains(&settings.isoform_frac) {
        return Err(usage("isoform_frac out of range (0.0-1.0)"));
    }
    if settings.max_isoforms < 0 {
        return Err(usage("max_isoforms < 0"));
    }
    if settings.path_kmax < 0 {
        return Err(usage("kmax must be >= 0"));
    }
    if settings.max_paths < 0 {
        return Err(usage("max_paths must be >= 0"));
    }
    if !(0.0..=1.0).contains(&settings.path_frac) {
        return Err(usage("path_frac not in range (0.0-1.0)"));
    }
    if settings.change_point {
        if !(0.0..=1.0).contains(&settings.change_point_pvalue) {
            return Err(usage("change point pvalue invalid"));
        }
        if !(0.0..=1.0).contains(&settings.change_point_fold_change) {
            return Err(usage("change point fold change invalid"));
        }
    }
    if settings.filter_splice_juncs && settings.genome_fasta.is_none() {
        return Err(usage("Filtering of splice junctions enabled (--filter-splice-juncs) but reference genome FASTA file not specified"));
    }
    if let Some(fasta) = &settings.genome_fasta {
        if settings.filter_splice_juncs && !fasta.exists() {
            return Err(usage("Reference genome FASTA file not found"));
        }
    }
    if !settings.splice_motifs.is_empty() && !settings.filter_splice_juncs {
        return Err(usage(
            "Additional splice motifs require --filter-splice-juncs",
        ));
    }
    for motif in &settings.splice_motifs {
        if motif.len() != 4 {
            return Err(usage("Only motifs of length 4 are allowed"));
        }
    }
    if settings.ref_gtf.is_none()
        && (settings.guided_strand || settings.guided_ends || settings.guided_assembly)
    {
        return Err(usage(
            "Guided assembly modes require a reference GTF (--ref-gtf)",
        ));
    }
    if let Some(reference) = &settings.ref_gtf {
        if !reference.exists() {
            return Err(usage(format!(
                "reference GTF file {} not found",
                reference.display()
            )));
        }
    }
    Ok(())
}

fn read_samples(path: &Path) -> Result<Vec<Sample>, TacoFailure> {
    let text = fs::read_to_string(path)
        .map_err(|err| runtime(format!("read {}: {err}", path.display())))?;
    let mut samples = Vec::new();
    let mut files = HashSet::new();
    let mut ids = HashSet::new();
    let mut next_id = 1u64;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut fields = line.split('\t');
        let gtf = fields.next().unwrap_or("").to_string();
        if gtf.is_empty() || !Path::new(&gtf).exists() {
            return Err(runtime(format!("GTF file '{gtf}' is not valid")));
        }
        if !files.insert(gtf.clone()) {
            return Err(runtime(format!("GTF file '{gtf}' is not unique")));
        }
        let id = if let Some(id) = fields.next() {
            let id = id.to_string();
            if id.is_empty() {
                return Err(runtime(format!("sample_id '{id}' is not unique")));
            }
            id
        } else {
            let id = next_id.to_string();
            next_id += 1;
            id
        };
        if !ids.insert(id.clone()) {
            return Err(runtime(format!("sample_id '{id}' is not unique")));
        }
        samples.push(Sample {
            path: gtf,
            id,
            is_ref: false,
        });
    }
    if samples.is_empty() {
        return Err(runtime("no input samples"));
    }
    Ok(samples)
}

fn parse_gtf(
    path: &Path,
    sample_id: &str,
    expr_attr: &str,
    is_ref: bool,
) -> Result<(Vec<Transfrag>, f64), String> {
    let text = fs::read_to_string(path).map_err(|err| format!("read {}: {err}", path.display()))?;
    let mut order = Vec::new();
    let mut map: HashMap<String, Transfrag> = HashMap::new();
    let mut total = 0.0;
    let mut next_id = 1usize;
    for raw in text.lines() {
        if raw.is_empty() || raw.starts_with('#') || raw.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = raw.split('\t').collect();
        if fields.len() < 9 {
            return Err(format!(
                "GTF line has fewer than 9 fields in {}",
                path.display()
            ));
        }
        let feature = fields[2];
        if feature != "transcript" && feature != "exon" {
            continue;
        }
        let start: i64 = fields[3]
            .parse::<i64>()
            .map_err(|_| format!("bad GTF start in {}", path.display()))?
            - 1;
        let end: i64 = fields[4]
            .parse()
            .map_err(|_| format!("bad GTF end in {}", path.display()))?;
        let strand = Strand::from_gtf(fields[6])
            .ok_or_else(|| format!("invalid strand '{}' in {}", fields[6], path.display()))?;
        let attrs = parse_attrs(fields[8]);
        let transcript_id = attrs
            .get("transcript_id")
            .cloned()
            .ok_or_else(|| format!("transcript_id attribute not found in {}", path.display()))?;
        if feature == "transcript" {
            if map.contains_key(&transcript_id) {
                return Err(format!("Transcript '{transcript_id}' duplicate detected"));
            }
            let expr = if is_ref {
                0.0
            } else {
                let text = attrs
                    .get(expr_attr)
                    .ok_or_else(|| format!("GTF expression attribute '{expr_attr}' not found"))?;
                text.parse::<f64>().map_err(|_| {
                    format!("bad expression '{text}' for transcript '{transcript_id}'")
                })?
            };
            total += expr;
            let id = format!("{sample_id}.{next_id}");
            next_id += 1;
            order.push(transcript_id.clone());
            map.insert(
                transcript_id,
                Transfrag {
                    chrom: fields[0].to_string(),
                    strand,
                    exons: Vec::new(),
                    id,
                    expr,
                    is_ref,
                },
            );
        } else if let Some(transfrag) = map.get_mut(&transcript_id) {
            transfrag.exons.push((start, end));
        } else {
            return Err(format!("Transcript '{transcript_id}' exon feature appeared in gtf file prior to transcript feature"));
        }
    }
    let mut transfrags = Vec::new();
    for id in order {
        let transfrag = map.remove(&id).unwrap();
        if transfrag.exons.is_empty() {
            return Err(format!("Transcript '{}' has no exons", transfrag.id));
        }
        transfrags.push(transfrag);
    }
    Ok((transfrags, total))
}

fn parse_attrs(text: &str) -> HashMap<String, String> {
    let mut attrs = HashMap::new();
    if text == "." {
        return attrs;
    }
    for piece in text.split(';') {
        let piece = piece.trim();
        if piece.is_empty() {
            continue;
        }
        let parts: Vec<&str> = piece.split(' ').collect();
        if parts.is_empty() || parts[0].is_empty() {
            continue;
        }
        let value = parts[1..].join(" ").trim_matches('"').to_string();
        attrs.insert(parts[0].to_string(), value);
    }
    attrs
}

fn transfrag_bed(transfrag: &Transfrag) -> String {
    let start = transfrag.start();
    let end = transfrag.end();
    let sizes = transfrag
        .exons
        .iter()
        .map(|(exon_start, exon_end)| (exon_end - exon_start).to_string())
        .collect::<Vec<_>>()
        .join(",");
    let starts = transfrag
        .exons
        .iter()
        .map(|(exon_start, _)| (exon_start - start).to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t0\t0\t0\t{}\t{}\t{}",
        transfrag.chrom,
        start,
        end,
        transfrag.id,
        py_float(transfrag.expr),
        transfrag.strand.as_gtf(),
        transfrag.exons.len(),
        sizes,
        starts
    )
}

fn assembly_bed_line(
    chrom: &str,
    name: &str,
    strand: Strand,
    score: i64,
    exons: &[(i64, i64)],
) -> String {
    let start = exons[0].0;
    let end = exons.last().unwrap().1;
    let sizes = exons
        .iter()
        .map(|(exon_start, exon_end)| (exon_end - exon_start).to_string())
        .collect::<Vec<_>>()
        .join(",");
    let starts = exons
        .iter()
        .map(|(exon_start, _)| (exon_start - start).to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!("{chrom}\t{start}\t{end}\t{name}\t{score}\t{}\t{start}\t{start}\t0\t{}\t{sizes},\t{starts},", strand.as_gtf(), exons.len())
}

fn bedgraph_array(values: &[f32], chrom: &str, origin: i64) -> Vec<String> {
    let mut lines = Vec::new();
    let mut index = 0;
    while index < values.len() {
        let value = values[index];
        let mut next = index + 1;
        while next < values.len() && values[next] == value {
            next += 1;
        }
        if value != 0.0 {
            lines.push(format!(
                "{chrom}\t{}\t{}\t{:.6}",
                origin + index as i64,
                origin + next as i64,
                f64::from(value)
            ));
        }
        index = next;
    }
    lines
}

fn allowed_motifs(extra: &[String]) -> HashSet<String> {
    let mut motifs = HashSet::from(["GTAG".to_string(), "GCAG".to_string(), "ATAC".to_string()]);
    for motif in extra {
        motifs.insert(motif.to_ascii_uppercase());
    }
    motifs
}

fn bad_introns(
    transfrag: &Transfrag,
    genome: &HashMap<String, String>,
    motifs: &HashSet<String>,
) -> Result<usize, String> {
    let mut bad = 0usize;
    for (start, end) in transfrag.introns() {
        if end < start + 2 {
            return Err(format!(
                "splice junction {start}-{end} on {} is shorter than 2",
                transfrag.chrom
            ));
        }
        let mut motif = fetch(genome, &transfrag.chrom, start, start + 2)?;
        motif.push_str(&fetch(genome, &transfrag.chrom, end - 2, end)?);
        if transfrag.strand == Strand::Neg {
            motif = reverse_complement(&motif)?;
        }
        if !motifs.contains(&motif) {
            bad += 1;
        }
    }
    Ok(bad)
}

fn fetch(
    genome: &HashMap<String, String>,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<String, String> {
    let Some(seq) = genome.get(chrom) else {
        return Err(format!("reference sequence '{chrom}' not found"));
    };
    if start < 0 || end > seq.len() as i64 || start > end {
        return Err(format!(
            "splice junction {start}-{end} on {chrom} is outside the reference sequence"
        ));
    }
    Ok(seq[start as usize..end as usize].to_string())
}

fn reverse_complement(seq: &str) -> Result<String, String> {
    seq.chars()
        .rev()
        .map(|base| match base {
            'A' => Ok('T'),
            'T' => Ok('A'),
            'G' => Ok('C'),
            'C' => Ok('G'),
            'N' => Ok('N'),
            other => Err(format!("cannot complement base '{other}'")),
        })
        .collect()
}

fn load_fasta(path: &Path) -> Result<HashMap<String, String>, String> {
    let mut genome = HashMap::new();
    for record in braid_io::read_fasta(path)? {
        if genome.insert(record.id.clone(), record.seq).is_some() {
            return Err(format!("duplicate fasta sequence {}", record.id));
        }
    }
    Ok(genome)
}

fn join_lines(lines: Vec<String>) -> String {
    if lines.is_empty() {
        String::new()
    } else {
        let mut text = lines.join("\n");
        text.push('\n');
        text
    }
}

fn with_header(header: &str, body: String) -> String {
    if body.is_empty() {
        format!("{header}\n")
    } else {
        format!("{header}\n{body}")
    }
}

fn sort_bed(lines: Vec<String>) -> String {
    let mut rows: Vec<(String, i64, usize, String)> = lines
        .into_iter()
        .enumerate()
        .map(|(index, line)| {
            let mut fields = line.split('\t');
            let chrom = fields.next().unwrap_or("").to_string();
            let start = fields
                .next()
                .and_then(|text| text.parse().ok())
                .unwrap_or(0);
            (chrom, start, index, line)
        })
        .collect();
    rows.sort();
    join_lines(rows.into_iter().map(|row| row.3).collect())
}

fn sort_gtf(lines: Vec<String>) -> String {
    let mut rows: Vec<(String, i64, i64, usize, String)> = lines
        .into_iter()
        .enumerate()
        .map(|(index, line)| {
            let fields: Vec<&str> = line.split('\t').collect();
            let chrom = fields.first().copied().unwrap_or("").to_string();
            let start = fields
                .get(3)
                .and_then(|text| text.parse().ok())
                .unwrap_or(0);
            let feature = if fields.get(2).copied() == Some("transcript") {
                0
            } else {
                1
            };
            (chrom, start, feature, index, line)
        })
        .collect();
    rows.sort();
    join_lines(rows.into_iter().map(|row| row.4).collect())
}

pub fn write_taco_texts(dir: &Path, texts: &TacoTexts) -> Result<(), TacoFailure> {
    if dir.exists() {
        return Err(usage(format!(
            "Output directory '{}' already exists",
            dir.display()
        )));
    }
    fs::create_dir_all(dir)
        .map_err(|err| runtime(format!("failed to create {}: {err}", dir.display())))?;
    let files = [
        ("samples.txt", &texts.samples),
        ("sample_stats.txt", &texts.sample_stats),
        ("transfrags.bed", &texts.transfrags_bed),
        ("transfrags.filtered.bed", &texts.transfrags_filtered_bed),
        ("loci.txt", &texts.loci),
        ("expr.pos.bedgraph", &texts.bedgraph_pos),
        ("expr.neg.bedgraph", &texts.bedgraph_neg),
        ("expr.none.bedgraph", &texts.bedgraph_none),
        ("splice_junctions.bed", &texts.splice_junctions),
        ("splice_graph.gtf", &texts.splice_graph),
        ("change_points.gtf", &texts.change_points),
        ("path_graph_stats.txt", &texts.path_graph_stats),
        ("assembly.gtf", &texts.assembly_gtf),
        ("assembly.bed", &texts.assembly_bed),
    ];
    for (name, text) in files {
        let path = dir.join(name);
        fs::write(&path, text)
            .map_err(|err| runtime(format!("failed to write {}: {err}", path.display())))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn case(name: &str, gtfs: &[(&str, &str)]) -> (PathBuf, TacoSettings) {
        let dir = std::env::temp_dir().join(format!("braid-taco-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        let mut rows = String::new();
        for (file, body) in gtfs {
            let path = dir.join(file);
            fs::write(&path, body).unwrap();
            rows.push_str(&format!("{}\n", path.display()));
        }
        let sample = dir.join("samples.tsv");
        fs::write(&sample, rows).unwrap();
        let settings = TacoSettings {
            sample_file: sample,
            ..TacoSettings::default()
        };
        (dir, settings)
    }

    const THREE: &str = "\
chr1\tA\ttranscript\t1000\t2000\t.\t+\t.\ttranscript_id \"t1\"; FPKM \"4\";\n\
chr1\tA\texon\t1000\t1200\t.\t+\t.\ttranscript_id \"t1\";\n\
chr1\tA\texon\t1400\t1600\t.\t+\t.\ttranscript_id \"t1\";\n\
chr1\tA\texon\t1800\t2000\t.\t+\t.\ttranscript_id \"t1\";\n";

    #[test]
    fn three_exon_isoform_matches_taco_outputs() {
        let (dir, settings) = case("three", &[("a.gtf", THREE)]);
        let texts = assemble(&settings).unwrap();
        assert_eq!(texts.assembly_gtf, "\
chr1\ttaco\ttranscript\t1000\t2000\t1000\t+\t.\texpr \"1000000.000\"; rel_frac \"1.00000\"; abs_frac \"1.00000\"; locus_id \"L1\"; gene_id \"G1\"; tss_id \"TSS1\"; transcript_id \"TU1\";\n\
chr1\ttaco\texon\t1000\t1200\t1000\t+\t.\tlocus_id \"L1\"; gene_id \"G1\"; tss_id \"TSS1\"; transcript_id \"TU1\";\n\
chr1\ttaco\texon\t1400\t1600\t1000\t+\t.\tlocus_id \"L1\"; gene_id \"G1\"; tss_id \"TSS1\"; transcript_id \"TU1\";\n\
chr1\ttaco\texon\t1800\t2000\t1000\t+\t.\tlocus_id \"L1\"; gene_id \"G1\"; tss_id \"TSS1\"; transcript_id \"TU1\";\n");
        assert_eq!(texts.assembly_bed, "chr1\t999\t2000\tG1|TU1(1000000.0)\t1000\t+\t999\t999\t0\t3\t201,201,201,\t0,400,800,\n");
        assert_eq!(
            texts.transfrags_bed,
            "chr1\t999\t2000\t1.1\t1000000.0\t+\t0\t0\t0\t3\t201,201,201\t0,400,800\n"
        );
        assert_eq!(texts.transfrags_filtered_bed, "");
        assert_eq!(texts.loci, "L1\tchr1\t999\t2000\t0\t1\n");
        assert_eq!(
            texts.bedgraph_pos,
            "\
chr1\t999\t1200\t1000000.000000\n\
chr1\t1399\t1600\t1000000.000000\n\
chr1\t1799\t2000\t1000000.000000\n"
        );
        assert_eq!(texts.bedgraph_neg, "");
        assert_eq!(texts.bedgraph_none, "");
        assert_eq!(texts.change_points, "");
        assert_eq!(
            texts.splice_junctions,
            "\
track name=junctions description=\"Splice Junctions\" graphType=junctions\n\
chr1\t1199\t1400\tJUNC\t1000000.0\t+\t1199\t1400\t255,0,0\t2\t1,1\t0,200\n\
chr1\t1599\t1800\tJUNC\t1000000.0\t+\t1599\t1800\t255,0,0\t2\t1,1\t0,200\n"
        );
        assert!(texts
            .splice_graph
            .contains("graph_id \"G_chr1_999_2000_+\""));
        assert!(texts.splice_graph.contains("expr_mean \"1000000.0\""));
        assert_eq!(
            texts.path_graph_stats,
            "\
chrom\tstart\tend\tstrand\tk\tkmax\ttransfrags\tshort_transfrags\tshort_expr\tlost_short\tlost_short_expr\tkmers\tlost_kmers\ttot_expr\tgraph_expr\texpr_frac\tvalid\topt\tis_opt\n\
chr1\t999\t2000\t+\t1\t3\t1\t0\t0.0\t0\t0.0\t3\t0\t1000000.0\t1000000.0\t1.0\t1\t3\t0\n\
chr1\t999\t2000\t+\t2\t3\t1\t0\t0.0\t0\t0.0\t2\t0\t1000000.0\t1000000.0\t1.0\t1\t2\t0\n\
chr1\t999\t2000\t+\t3\t3\t1\t0\t0.0\t0\t0.0\t1\t0\t1000000.0\t1000000.0\t1.0\t1\t1\t0\n\
chr1\t999\t2000\t+\t1\t3\t1\t0\t0.0\t0\t0.0\t3\t0\t1000000.0\t1000000.0\t1.0\t1\t3\t1\n"
        );
        assert_eq!(texts.transcripts.len(), 1);
        assert_eq!(texts.transcripts[0].score, Some(1_000_000.0));
        assert_eq!(
            texts.transcripts[0].exons,
            vec![
                Exon {
                    start: 1000,
                    end: 1201
                },
                Exon {
                    start: 1400,
                    end: 1601
                },
                Exon {
                    start: 1800,
                    end: 2001
                }
            ]
        );
        assert_eq!(texts.transcripts[0].gene_id, "G1");
        assert_eq!(texts.transcripts[0].transcript_id, "TU1");
        assert_eq!(texts.sample_stats.lines().nth(1).unwrap(), "1\t1\t0\t0\t0");
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn two_samples_divide_expression_by_the_sample_count() {
        let (dir, settings) = case("two", &[("a.gtf", THREE), ("b.gtf", THREE)]);
        let texts = assemble(&settings).unwrap();
        assert!(texts.assembly_gtf.contains("expr \"1000000.000\""));
        assert!(texts.bedgraph_pos.contains("2000000.000000"));
        assert!(texts.splice_junctions.contains("2000000.0"));
        assert_eq!(texts.sample_stats.lines().nth(1).unwrap(), "1\t1\t0\t0\t0");
        assert_eq!(texts.sample_stats.lines().nth(2).unwrap(), "2\t1\t0\t0\t0");
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn unstranded_transcript_takes_the_only_overlapping_strand() {
        let gtf = "\
chr1\tA\ttranscript\t1000\t1300\t.\t+\t.\ttranscript_id \"a\"; FPKM \"4\";\n\
chr1\tA\texon\t1000\t1300\t.\t+\t.\ttranscript_id \"a\";\n\
chr1\tA\ttranscript\t1000\t1300\t.\t.\t.\ttranscript_id \"b\"; FPKM \"4\";\n\
chr1\tA\texon\t1000\t1300\t.\t.\t.\ttranscript_id \"b\";\n";
        let (dir, settings) = case("impute", &[("a.gtf", gtf)]);
        let texts = assemble(&settings).unwrap();
        assert!(texts.assembly_gtf.contains("expr \"1000000.000\""));
        assert_eq!(texts.transcripts.len(), 1);
        assert_eq!(texts.transcripts[0].strand, ModelStrand::Forward);
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn short_transcript_is_filtered_and_still_normalizes_expression() {
        let gtf = "\
chr1\tA\ttranscript\t1000\t2000\t.\t+\t.\ttranscript_id \"t1\"; FPKM \"4\";\n\
chr1\tA\texon\t1000\t1200\t.\t+\t.\ttranscript_id \"t1\";\n\
chr1\tA\texon\t1400\t1600\t.\t+\t.\ttranscript_id \"t1\";\n\
chr1\tA\texon\t1800\t2000\t.\t+\t.\ttranscript_id \"t1\";\n\
chr1\tA\ttranscript\t3000\t3050\t.\t+\t.\ttranscript_id \"s\"; FPKM \"4\";\n\
chr1\tA\texon\t3000\t3050\t.\t+\t.\ttranscript_id \"s\";\n";
        let (dir, settings) = case("short", &[("a.gtf", gtf)]);
        let texts = assemble(&settings).unwrap();
        assert!(texts.assembly_gtf.contains("expr \"500000.000\""));
        assert!(texts.transfrags_filtered_bed.contains("1.2"));
        assert!(!texts.transfrags_bed.contains("1.2"));
        assert_eq!(texts.sample_stats.lines().nth(1).unwrap(), "1\t2\t1\t0\t0");
        let _ = fs::remove_dir_all(dir);
    }

    #[test]
    fn negative_strand_exon_is_reported_on_the_minus_strand() {
        let gtf = "\
chr1\tA\ttranscript\t1000\t1200\t.\t-\t.\ttranscript_id \"t1\"; FPKM \"4\";\n\
chr1\tA\texon\t1000\t1200\t.\t-\t.\ttranscript_id \"t1\";\n";
        let (dir, settings) = case("neg", &[("a.gtf", gtf)]);
        let texts = assemble(&settings).unwrap();
        assert!(texts.assembly_gtf.contains("\t-\t"));
        assert!(texts.assembly_gtf.contains("exon\t1000\t1200\t"));
        assert_eq!(texts.transcripts[0].strand, ModelStrand::Reverse);
        assert_eq!(
            texts.transcripts[0].exons,
            vec![Exon {
                start: 1000,
                end: 1201
            }]
        );
        let _ = fs::remove_dir_all(dir);
    }
}
