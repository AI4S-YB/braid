use std::ffi::OsString;
use std::path::PathBuf;
use std::process::ExitCode;

use braid_flair::{combine, write_combine_texts, CombineSettings};
use braid_tama::{
    collapse_to_files, merge_sources, read_merge_sources, write_merge_texts, CollapseSettings,
    MergeSettings, TAMA_COLLAPSE_DATE,
};
use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "braid",
    version,
    about = "Reconcile transcript models with TAMA collapse, TAMA merge and FLAIR combine.",
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
