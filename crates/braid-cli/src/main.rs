use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;

use braid_flair::{combine, write_combine_texts, CombineSettings};
use braid_scotch::{
    differential_usage, quantify, write_dtu, DtuSettings, Platform, ScotchSettings,
};
use braid_taco::{assemble, write_taco_texts, TacoSettings};
use braid_tama::{
    collapse_to_files, merge_sources, read_merge_sources, write_merge_texts, CollapseSettings,
    MergeSettings, TAMA_COLLAPSE_DATE,
};
use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "braid",
    version,
    about = "Reconcile transcript models with TAMA collapse, TAMA merge, FLAIR combine, TACO, and SCOTCH.",
    after_help = "TAMA supports capped/no_cap and original/low_mem modes. BAM input requires samtools on PATH. Legacy TAMA single-dash options are accepted; console logs may differ from upstream."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Collapse alignments into transcript models.
    Collapse(CollapseArgs),
    /// Merge annotation sets with TAMA.
    Merge(MergeArgs),
    /// Combine transcriptomes. FLAIR combine is available.
    Combine(CombineArgs),
    /// TAMA collapse with original parameter spellings.
    TamaCollapse(TamaCollapseArgs),
    /// TAMA merge with original parameter spellings.
    TamaMerge(TamaMergeArgs),
    /// Same flags as `flair combine`.
    FlairCombine(FlairCombineArgs),
    /// Meta-assemble sample GTFs with TACO.
    Taco(TacoArgs),
    /// Quantify single-cell full-length isoforms with SCOTCH.
    Scotch(ScotchArgs),
    /// Test differential transcript usage between two SCOTCH count directories.
    ScotchDtu(ScotchDtuArgs),
}

#[derive(Args)]
struct CollapseArgs {
    /// Algorithm. Only `tama` is available.
    #[arg(long, default_value = "tama")]
    algo: String,
    #[command(flatten)]
    tama: TamaCollapseArgs,
}

#[derive(Args)]
struct MergeArgs {
    /// Algorithm. Only `tama` is available.
    #[arg(long, default_value = "tama")]
    algo: String,
    #[command(flatten)]
    tama: TamaMergeArgs,
}

#[derive(Args)]
struct CombineArgs {
    /// Algorithm. Only `flair` is available.
    #[arg(long, default_value = "flair")]
    algo: String,
    #[command(flatten)]
    flair: FlairCombineArgs,
}

#[derive(Args, Debug)]
struct FlairCombineArgs {
    /// Manifest. Each row is sample, type, bed, and optional fasta and read map.
    #[arg(short = 'm', long = "manifest")]
    manifest: PathBuf,
    /// Output prefix.
    #[arg(
        short = 'o',
        long = "output_prefix",
        default_value = "flair.combined.isoforms"
    )]
    output_prefix: PathBuf,
    /// Window for comparing ends of isoforms with the same intron chain.
    #[arg(short = 'w', long = "endwindow", default_value_t = 200)]
    endwindow: i64,
    /// Minimum percent usage required in one sample. 10 means usage greater than 0.10.
    #[arg(short = 'p', long = "minpercentusage", default_value_t = 10)]
    minpercentusage: i64,
    /// Also write a GTF from the combined BED.
    #[arg(short = 'c', long = "convert_gtf")]
    convert_gtf: bool,
    /// Keep single-exon isoforms. Off by default.
    #[arg(short = 's', long = "include_se")]
    include_se: bool,
    /// `usageandlongest`, `usageonly`, `none`, or a read-count threshold.
    #[arg(short = 'f', long = "filter", default_value = "usageandlongest")]
    filter: String,
}

#[derive(Args, Debug)]
struct TamaCollapseArgs {
    /// Sorted SAM file.
    #[arg(short = 's', required_unless_present = "version_date")]
    sam: Option<PathBuf>,
    /// Genome FASTA.
    #[arg(short = 'f', required_unless_present = "version_date")]
    fasta: Option<PathBuf>,
    /// Output prefix.
    #[arg(short = 'p', required_unless_present = "version_date")]
    prefix: Option<PathBuf>,
    /// `capped` or `no_cap` (allows 5' degradation).
    #[arg(short = 'x', default_value = "capped")]
    cap: String,
    /// `common_ends` or `longest_ends`.
    #[arg(short = 'e', default_value = "common_ends")]
    ends: String,
    /// Coverage threshold.
    #[arg(short = 'c', default_value_t = 99.0)]
    coverage: f64,
    /// Identity threshold.
    #[arg(short = 'i', default_value_t = 85.0)]
    identity: f64,
    /// `ident_cov` or `ident_map`.
    #[arg(long = "icm", default_value = "ident_cov")]
    ident_method: String,
    /// 5' threshold.
    #[arg(short = 'a', default_value_t = 10)]
    five_prime: i64,
    /// Exon / splice-junction threshold.
    #[arg(short = 'm', default_value_t = 10)]
    exon: i64,
    /// 3' threshold.
    #[arg(short = 'z', default_value_t = 10)]
    three_prime: i64,
    /// `merge_dup` or `no_merge`.
    #[arg(short = 'd', default_value = "merge_dup")]
    duplicates: String,
    /// `no_priority` or `sj_priority`.
    #[arg(long = "sj", default_value = "no_priority")]
    sj_priority: String,
    /// Splice-junction error threshold.
    #[arg(long = "sjt", default_value_t = 10)]
    sj_threshold: i64,
    /// Local density error threshold.
    #[arg(long = "lde", default_value_t = 1000)]
    lde: i64,
    /// Match symbol in the simple error string. Default `_`.
    #[arg(long = "ses", default_value = "_")]
    ses: String,
    /// `SAM` or `BAM`. BAM requires samtools on PATH.
    #[arg(short = 'b')]
    bam: Option<String>,
    /// Recorded in the report; upstream console logs are not reproduced.
    #[arg(long = "log", default_value = "log_on")]
    log: String,
    /// `original` or `low_mem`. low_mem processes one locus at a time.
    #[arg(long = "rm", default_value = "original")]
    run_mode: String,
    /// Variant support threshold.
    #[arg(long = "vc", default_value_t = 5)]
    var_support: i64,
    /// Print the TAMA date and exit.
    #[arg(short = 'v')]
    version_date: Option<String>,
}

#[derive(Args, Debug)]
struct TamaMergeArgs {
    /// File list.
    #[arg(short = 'f')]
    filelist: PathBuf,
    /// Output prefix.
    #[arg(short = 'p')]
    prefix: PathBuf,
    /// `common_ends` or `longest_ends`.
    #[arg(short = 'e', default_value = "common_ends")]
    ends: String,
    /// 5' threshold. TAMA merge's code default is 20.
    #[arg(short = 'a', default_value_t = 20)]
    five_prime: i64,
    /// Exon / splice-junction threshold.
    #[arg(short = 'm', default_value_t = 10)]
    exon: i64,
    /// 3' threshold. TAMA merge's code default is 20.
    #[arg(short = 'z', default_value_t = 20)]
    three_prime: i64,
    /// `no_merge` or `merge_dup`. Merge's default is `no_merge`.
    #[arg(short = 'd', default_value = "no_merge")]
    duplicates: String,
    /// Source name whose gene and transcript ids are kept.
    #[arg(short = 's', default_value = "no_source_id")]
    source_id: String,
    /// Source name whose CDS is kept.
    #[arg(long = "cds", default_value = "no_cds")]
    cds: String,
}

#[derive(Args, Debug)]
struct TacoArgs {
    /// Sample table. Each row is a GTF path and an optional sample id.
    sample_file: PathBuf,
    /// Directory where TACO writes its output files.
    #[arg(short = 'o', long = "output-dir", default_value = "output")]
    output_dir: PathBuf,
    /// Accepted for compatibility. This port always uses one process so transcript ids stay reproducible.
    #[arg(short = 'p', long = "num-processes", default_value_t = 1)]
    num_processes: i64,
    /// Accepted for compatibility.
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,
    /// Not supported. This port does not resume a partial run.
    #[arg(long = "resume")]
    resume: bool,
    /// Not supported. Pass sample GTFs instead of a preexisting BED file.
    #[arg(long = "assemble")]
    assemble_bed: Option<PathBuf>,
    /// GTF attribute that holds expression.
    #[arg(long = "gtf-expr-attr", default_value = "FPKM")]
    gtf_expr_attr: String,
    /// Drop input transcripts shorter than this.
    #[arg(long = "filter-min-length", default_value_t = 200)]
    filter_min_length: i64,
    /// Drop input transcripts whose normalized expression is below this.
    #[arg(long = "filter-min-expr", default_value_t = 0.5)]
    filter_min_expr: f64,
    /// Keep isoforms at least this fraction of the major isoform in the gene.
    #[arg(long = "isoform-frac", default_value_t = 0.05)]
    isoform_frac: f64,
    /// Maximum isoforms to report for each gene. Zero means no limit.
    #[arg(long = "max-isoforms", default_value_t = 0)]
    max_isoforms: i64,
    /// Assemble transcripts that stay unstranded.
    #[arg(
        long = "assemble-unstranded",
        conflicts_with = "no_assemble_unstranded"
    )]
    assemble_unstranded: bool,
    /// Leave unstranded transcripts out of the assembly.
    #[arg(long = "no-assemble-unstranded")]
    no_assemble_unstranded: bool,
    /// Run change-point detection. On unless `--no-change-point` is set.
    #[arg(long = "change-point", conflicts_with = "no_change_point")]
    change_point: bool,
    /// Skip change-point detection.
    #[arg(long = "no-change-point")]
    no_change_point: bool,
    /// Mann-Whitney p-value threshold for a change point.
    #[arg(long = "change-point-pvalue", default_value_t = 0.01)]
    change_point_pvalue: f64,
    /// Fold-change threshold for a change point.
    #[arg(long = "change-point-fold-change", default_value_t = 0.85)]
    change_point_fold_change: f64,
    /// Trim transcript ends around change points. On unless `--no-change-point-trim` is set.
    #[arg(long = "change-point-trim", conflicts_with = "no_change_point_trim")]
    change_point_trim: bool,
    /// Leave transcript ends untrimmed.
    #[arg(long = "no-change-point-trim")]
    no_change_point_trim: bool,
    /// Largest k considered for the path graph. Zero uses the longest transcript.
    #[arg(long = "path-kmax", default_value_t = 0)]
    path_kmax: i64,
    /// Stop the path search below this fraction of the best path.
    #[arg(long = "path-frac", default_value_t = 0.0)]
    path_frac: f64,
    /// Stop after this many paths. Zero means no limit.
    #[arg(long = "max-paths", default_value_t = 0)]
    max_paths: i64,
    /// Reference GTF used by the guided modes.
    #[arg(long = "ref-gtf")]
    ref_gtf: Option<PathBuf>,
    /// Resolve unstranded transcripts with the reference strand.
    #[arg(long = "guided-strand")]
    guided_strand: bool,
    /// Use reference transcript ends as graph boundaries.
    #[arg(long = "guided-ends")]
    guided_ends: bool,
    /// Include reference transcripts in the path graph.
    #[arg(long = "guided-assembly")]
    guided_assembly: bool,
    /// Drop transcripts with splice motifs other than GTAG, GCAG, and ATAC.
    #[arg(long = "filter-splice-juncs")]
    filter_splice_juncs: bool,
    /// Extra 4-base splice motif allowed by `--filter-splice-juncs`. Repeatable.
    #[arg(long = "add-splice-motif")]
    splice_motifs: Vec<String>,
    /// Genome FASTA used by `--filter-splice-juncs`.
    #[arg(long = "ref-genome-fasta")]
    genome_fasta: Option<PathBuf>,
}

#[derive(Args, Debug)]
struct ScotchArgs {
    /// BAM or SAM file, or a directory of them. Repeat for one sample per path.
    #[arg(long = "bam", required = true)]
    bam: Vec<PathBuf>,
    /// Reference GTF.
    #[arg(long)]
    gtf: PathBuf,
    /// Output directory. Each sample is written underneath it.
    #[arg(long)]
    out: PathBuf,
    /// Genome FASTA. Used to reject internal priming when a poly(A) or poly(T) is detected.
    #[arg(long)]
    fasta: Option<PathBuf>,
    /// Refine sub-exons from read coverage before assignment.
    #[arg(long)]
    update_gtf: bool,
    /// `10x-ont`, `10x-pacbio`, `parse-ont`, or `bulk`.
    #[arg(long, default_value = "10x-ont")]
    platform: String,
    /// Cell-barcode tag. Defaults are CB for 10x and the read name for Parse.
    #[arg(long)]
    barcode_cell: Option<String>,
    /// UMI tag. Defaults are UB for 10x-ont and XM for 10x-pacbio.
    #[arg(long)]
    barcode_umi: Option<String>,
    /// Limit quantification to these gene names or ids. Repeat the flag for more than one gene.
    #[arg(long = "gene")]
    gene: Vec<String>,
    /// Exon coverage at or below this fraction is a skip.
    #[arg(long, default_value_t = 0.1)]
    match_low: f64,
    /// Exon coverage at or above this fraction is a match.
    #[arg(long, default_value_t = 0.6)]
    match_high: f64,
    /// Ignore exons shorter than this when matching known isoforms.
    #[arg(long, default_value_t = 0)]
    small_exon: i64,
    /// Upper bound of the small-exon threshold.
    #[arg(long, default_value_t = 80)]
    small_exon_high: i64,
    /// Truncated ends at or above this fraction are treated as fully covered.
    #[arg(long, default_value_t = 0.4)]
    truncation_match: f64,
    /// Drop novel isoforms supported by fewer reads than this.
    #[arg(long, default_value_t = 0)]
    novel_read_n: i64,
    /// Drop novel isoforms whose share of known plus novel reads is below this.
    #[arg(long, default_value_t = 0.0)]
    novel_read_pct: f64,
    /// Keep each novel isoform separate instead of grouping truncated copies.
    #[arg(long)]
    no_group_novel: bool,
    /// Minimum coverage as a fraction of the peak when refining exons.
    #[arg(long, default_value_t = 0.02)]
    coverage_exon: f64,
    /// Minimum junction support as a fraction of the strongest junction.
    #[arg(long, default_value_t = 0.02)]
    coverage_splice: f64,
    /// Z-score for splitting an exon at a coverage change.
    #[arg(long, default_value_t = 10.0)]
    z_score: f64,
    /// Accepted for compatibility. Quantification runs in one process.
    #[arg(long, default_value_t = 1)]
    workers: i64,
}

#[derive(Args, Debug)]
struct ScotchDtuArgs {
    /// Sample directory, or its `count_matrix` directory, for group A.
    #[arg(long)]
    a: Option<PathBuf>,
    /// Sample directory, or its `count_matrix` directory, for group B.
    #[arg(long)]
    b: Option<PathBuf>,
    /// Gene count CSV for group A.
    #[arg(long)]
    gene_a: Option<PathBuf>,
    /// Transcript count CSV for group A.
    #[arg(long)]
    transcript_a: Option<PathBuf>,
    /// Gene-to-transcript map for group A.
    #[arg(long)]
    map_a: Option<PathBuf>,
    /// Gene count CSV for group B.
    #[arg(long)]
    gene_b: Option<PathBuf>,
    /// Transcript count CSV for group B.
    #[arg(long)]
    transcript_b: Option<PathBuf>,
    /// Gene-to-transcript map for group B.
    #[arg(long)]
    map_b: Option<PathBuf>,
    /// Output TSV.
    #[arg(long)]
    out: PathBuf,
    /// Added to gene counts before the Wilcoxon test.
    #[arg(long, default_value_t = 0.01)]
    epsilon: f64,
    /// Sum transcript columns whose names contain `novel` before the test.
    #[arg(long)]
    group_novel: bool,
    /// Isoforms below this usage in both groups collapse into `other`.
    #[arg(long, default_value_t = 0.05)]
    rare: f64,
}

// clap's short options are single characters. Translate only TAMA option
// tokens, preserving values and non-UTF-8 paths, before handing them to clap.
fn normalize_tama_args(mut args: Vec<OsString>) -> Vec<OsString> {
    let Some(command) = args.get(1).and_then(|s| s.to_str()) else {
        return args;
    };
    let legacy: &[&str] = match command {
        "collapse" | "tama-collapse" => &["icm", "sj", "sjt", "lde", "ses", "log", "rm", "vc"],
        "merge" | "tama-merge" => &["cds"],
        _ => return args,
    };
    let mut value_next = false;
    for arg in args.iter_mut().skip(2) {
        if value_next {
            value_next = false;
            continue;
        }
        let Some(text) = arg.to_str() else { continue };
        if text == "--" {
            break;
        }
        let (option, attached) = text
            .split_once('=')
            .map_or((text, false), |(key, _)| (key, true));
        let old = option
            .strip_prefix('-')
            .is_some_and(|key| legacy.contains(&key));
        value_next = !attached
            && (old
                || matches!(
                    option,
                    "-s" | "-f"
                        | "-p"
                        | "-x"
                        | "-e"
                        | "-c"
                        | "-i"
                        | "-a"
                        | "-m"
                        | "-z"
                        | "-d"
                        | "-b"
                        | "-v"
                        | "--algo"
                )
                || option
                    .strip_prefix("--")
                    .is_some_and(|key| legacy.contains(&key)));
        if old {
            *arg = OsString::from(format!("-{text}"));
        }
    }
    args
}

fn main() -> ExitCode {
    let cli = Cli::parse_from(normalize_tama_args(std::env::args_os().collect()));
    match cli.command {
        Command::Collapse(args) => {
            if args.algo != "tama" {
                eprintln!("unknown collapse algorithm '{}'", args.algo);
                return ExitCode::from(2);
            }
            run_collapse(args.tama)
        }
        Command::TamaCollapse(args) => run_collapse(args),
        Command::Merge(args) => {
            if args.algo != "tama" {
                eprintln!("unknown merge algorithm '{}'", args.algo);
                return ExitCode::from(2);
            }
            run_merge(args.tama)
        }
        Command::TamaMerge(args) => run_merge(args),
        Command::Combine(args) => {
            if args.algo != "flair" {
                eprintln!("unknown combine algorithm '{}'", args.algo);
                return ExitCode::from(2);
            }
            run_flair_combine(args.flair)
        }
        Command::FlairCombine(args) => run_flair_combine(args),
        Command::Taco(args) => run_taco(args),
        Command::Scotch(args) => run_scotch(args),
        Command::ScotchDtu(args) => run_scotch_dtu(args),
    }
}

fn run_collapse(args: TamaCollapseArgs) -> ExitCode {
    if args.version_date.is_some() {
        println!("{TAMA_COLLAPSE_DATE}");
        println!("Program did not run");
        return ExitCode::SUCCESS;
    }
    let (Some(sam), Some(fasta), Some(prefix)) = (args.sam, args.fasta, args.prefix) else {
        eprintln!("TAMA collapse requires -s, -f and -p");
        return ExitCode::from(2);
    };
    let settings = CollapseSettings {
        cap: args.cap,
        ends: args.ends,
        coverage: args.coverage,
        identity: args.identity,
        ident_method: args.ident_method,
        five_prime: args.five_prime,
        exon_diff: args.exon,
        three_prime: args.three_prime,
        duplicates: args.duplicates,
        sj_priority: args.sj_priority,
        sj_threshold: args.sj_threshold,
        lde: args.lde,
        ses: args.ses,
        bam: args.bam.unwrap_or_else(|| "SAM".to_string()),
        log: args.log,
        run_mode: args.run_mode,
        var_support: args.var_support,
        sam_label: sam.display().to_string(),
        fasta_label: fasta.display().to_string(),
        prefix_label: prefix.display().to_string(),
    };
    match collapse_to_files(&sam, &fasta, &prefix, &settings) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("{err}");
            ExitCode::from(1)
        }
    }
}

fn run_flair_combine(args: FlairCombineArgs) -> ExitCode {
    let settings = CombineSettings {
        end_window: args.endwindow,
        min_percent: args.minpercentusage,
        include_single_exon: args.include_se,
        filter: args.filter,
        convert_gtf: args.convert_gtf,
    };
    let texts = match combine(&args.manifest, &settings) {
        Ok(texts) => texts,
        Err(err) => {
            eprintln!("{err}");
            return ExitCode::from(1);
        }
    };
    match write_combine_texts(&args.output_prefix, &texts) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("{err}");
            ExitCode::from(1)
        }
    }
}

fn run_taco(args: TacoArgs) -> ExitCode {
    if args.resume {
        eprintln!("--resume is not supported");
        return ExitCode::from(2);
    }
    if args.assemble_bed.is_some() {
        eprintln!("--assemble is not supported");
        return ExitCode::from(2);
    }
    let _ = (
        args.num_processes,
        args.verbose,
        args.change_point,
        args.change_point_trim,
        args.no_assemble_unstranded,
    );
    if args.output_dir.exists() {
        eprintln!(
            "Output directory '{}' already exists",
            args.output_dir.display()
        );
        return ExitCode::from(2);
    }
    let settings = TacoSettings {
        sample_file: args.sample_file,
        gtf_expr_attr: args.gtf_expr_attr,
        filter_min_length: args.filter_min_length,
        filter_min_expr: args.filter_min_expr,
        isoform_frac: args.isoform_frac,
        max_isoforms: args.max_isoforms,
        assemble_unstranded: args.assemble_unstranded,
        change_point: !args.no_change_point,
        change_point_pvalue: args.change_point_pvalue,
        change_point_fold_change: args.change_point_fold_change,
        change_point_trim: !args.no_change_point_trim,
        path_kmax: args.path_kmax,
        path_frac: args.path_frac,
        max_paths: args.max_paths,
        ref_gtf: args.ref_gtf,
        guided_strand: args.guided_strand,
        guided_ends: args.guided_ends,
        guided_assembly: args.guided_assembly,
        filter_splice_juncs: args.filter_splice_juncs,
        splice_motifs: args.splice_motifs,
        genome_fasta: args.genome_fasta,
    };
    match assemble(&settings).and_then(|texts| write_taco_texts(&args.output_dir, &texts)) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("{err}");
            ExitCode::from(if err.usage { 2 } else { 1 })
        }
    }
}

fn run_scotch(args: ScotchArgs) -> ExitCode {
    let platform = match args.platform.as_str() {
        "10x-ont" => Platform::TenxOnt,
        "10x-pacbio" => Platform::TenxPacbio,
        "parse-ont" => Platform::ParseOnt,
        "bulk" => Platform::Bulk,
        other => {
            eprintln!("unknown platform '{other}'");
            return ExitCode::from(2);
        }
    };
    if args.novel_read_n < 0 || args.workers < 0 {
        eprintln!("--novel-read-n and --workers must be non-negative");
        return ExitCode::from(2);
    }
    let mut settings = ScotchSettings {
        bams: args.bam,
        gtf: args.gtf,
        out: args.out,
        fasta: args.fasta,
        update_gtf: args.update_gtf,
        platform,
        barcode_cell: args.barcode_cell,
        barcode_umi: args.barcode_umi,
        genes: args.gene,
        novel_read_n: args.novel_read_n as usize,
        novel_read_pct: args.novel_read_pct,
        group_novel: !args.no_group_novel,
        exon_fraction: args.coverage_exon,
        splice_fraction: args.coverage_splice,
        z_score: args.z_score,
        workers: args.workers as usize,
        ..ScotchSettings::default()
    };
    settings.params.low = args.match_low;
    settings.params.high = args.match_high;
    settings.params.small_exon = args.small_exon;
    settings.params.small_exon_high = args.small_exon_high;
    settings.params.truncation_match = args.truncation_match;
    match quantify(&settings) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("{err}");
            ExitCode::from(if err.usage { 2 } else { 1 })
        }
    }
}

fn count_table(
    dir: Option<PathBuf>,
    gene: Option<PathBuf>,
    transcript: Option<PathBuf>,
    map: Option<PathBuf>,
) -> Result<(PathBuf, PathBuf, PathBuf), String> {
    let base = dir.as_ref().map(|path| {
        if path.join("count_matrix").is_dir() {
            path.join("count_matrix")
        } else {
            path.clone()
        }
    });
    let pick = |explicit: Option<PathBuf>, name: &str| -> Result<PathBuf, String> {
        if let Some(path) = explicit {
            return Ok(path);
        }
        base.as_ref()
            .map(|dir| dir.join(name))
            .ok_or_else(|| format!("provide --a/--b or an explicit {name} path"))
    };
    Ok((
        pick(gene, "gene_counts.csv")?,
        pick(transcript, "transcript_counts.csv")?,
        pick(map, "gene_transcript.tsv")?,
    ))
}

fn run_scotch_dtu(args: ScotchDtuArgs) -> ExitCode {
    let group_a = match count_table(args.a, args.gene_a, args.transcript_a, args.map_a) {
        Ok(paths) => paths,
        Err(err) => {
            eprintln!("{err}");
            return ExitCode::from(2);
        }
    };
    let group_b = match count_table(args.b, args.gene_b, args.transcript_b, args.map_b) {
        Ok(paths) => paths,
        Err(err) => {
            eprintln!("{err}");
            return ExitCode::from(2);
        }
    };
    let settings = DtuSettings {
        epsilon: args.epsilon,
        group_novel: args.group_novel,
        rare: args.rare,
        ..DtuSettings::default()
    };
    match differential_usage(
        &group_a.0, &group_a.1, &group_a.2, &group_b.0, &group_b.1, &group_b.2, &settings,
    )
    .and_then(|tables| write_dtu(&args.out, &tables))
    {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("{err}");
            ExitCode::from(if err.usage { 2 } else { 1 })
        }
    }
}

fn run_merge(args: TamaMergeArgs) -> ExitCode {
    let settings = MergeSettings {
        ends: args.ends,
        five_prime: args.five_prime,
        exon_diff: args.exon,
        three_prime: args.three_prime,
        duplicates: args.duplicates,
        source_id: args.source_id,
        cds: args.cds,
    };
    let result = read_merge_sources(&args.filelist)
        .and_then(|sources| merge_sources(&sources, &settings))
        .and_then(|texts| write_merge_texts(&args.prefix, &texts));
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("{e}");
            ExitCode::from(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn normalized(args: &[&str]) -> Vec<OsString> {
        normalize_tama_args(args.iter().map(OsString::from).collect())
    }

    #[test]
    fn legacy_options_do_not_rewrite_values_or_other_algorithms() {
        assert_eq!(
            normalized(&[
                "braid",
                "tama-collapse",
                "-p",
                "-rm",
                "-icm=ident_map",
                "--",
                "-log"
            ]),
            [
                "braid",
                "tama-collapse",
                "-p",
                "-rm",
                "--icm=ident_map",
                "--",
                "-log"
            ]
            .map(OsString::from),
        );
        let flair = ["braid", "flair-combine", "-m", "-rm", "-icm"];
        assert_eq!(normalized(&flair), flair.map(OsString::from));
    }

    #[cfg(unix)]
    #[test]
    fn legacy_options_preserve_non_utf8_paths() {
        use std::os::unix::ffi::OsStringExt;
        let path = OsString::from_vec(b"input-\xff.sam".to_vec());
        let args = vec![
            "braid".into(),
            "collapse".into(),
            "-s".into(),
            path.clone(),
            "-rm".into(),
            "original".into(),
        ];
        let got = normalize_tama_args(args);
        assert_eq!(got[3], path);
        assert_eq!(got[4], "--rm");
    }
}
