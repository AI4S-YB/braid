//! Pooled transcriptome construction.
//!
//! `tama` and `flair` run in process. IsoQuant, Bambu, and StringTie are
//! pinned external programs. Every algorithm writes `transcripts.gtf`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::Args;
use plena_flair::{alignments_to_bed12, bed12_to_gtf, collapse_precise_bed, PreciseSettings};
use plena_tama::{collapse_to_files, CollapseSettings};

use crate::tools::{self, append_log, publish};

const BAMBU_DISCOVER: &str = include_str!("scripts/bambu_discover.R");

#[derive(Args, Debug)]
pub struct DiscoverArgs {
    /// `tama`, `flair`, `isoquant`, `bambu`, or `stringtie`.
    #[arg(long)]
    pub algo: String,
    /// Coordinate-sorted genome BAM. Cells are pooled.
    #[arg(long)]
    pub bam: PathBuf,
    /// Genome FASTA matching the BAM.
    #[arg(long)]
    pub genome: PathBuf,
    /// Guide annotation. Required for IsoQuant and Bambu.
    #[arg(long)]
    pub gtf: Option<PathBuf>,
    /// Output directory. `transcripts.gtf` is written here.
    #[arg(long)]
    pub out: PathBuf,
    /// Threads for external programs.
    #[arg(long, default_value_t = 1)]
    pub threads: i64,
}

pub fn run(args: DiscoverArgs) -> ExitCode {
    if args.threads < 1 {
        eprintln!("--threads must be at least 1");
        return ExitCode::from(2);
    }
    let result = match args.algo.as_str() {
        "tama" => discover_tama(&args),
        "flair" => discover_flair(&args),
        "isoquant" => discover_isoquant(&args),
        "bambu" => discover_bambu(&args),
        "stringtie" => discover_stringtie(&args),
        other => {
            eprintln!("unknown discover algorithm '{other}'");
            return ExitCode::from(2);
        }
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("{err}");
            ExitCode::from(1)
        }
    }
}

fn discover_tama(args: &DiscoverArgs) -> Result<(), String> {
    fs::create_dir_all(&args.out).map_err(|err| format!("create {}: {err}", args.out.display()))?;
    let raw = args.out.join("raw");
    fs::create_dir_all(&raw).map_err(|err| format!("create {}: {err}", raw.display()))?;
    let prefix = raw.join("tama");
    let settings = CollapseSettings {
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
        bam: "BAM".to_string(),
        log: "log_off".to_string(),
        run_mode: "original".to_string(),
        var_support: 5,
        sam_label: args.bam.display().to_string(),
        fasta_label: args.genome.display().to_string(),
        prefix_label: prefix.display().to_string(),
    };
    collapse_to_files(&args.bam, &args.genome, &prefix, &settings)?;
    let bed_path = raw.join("tama.bed");
    let bed = fs::read_to_string(&bed_path)
        .map_err(|err| format!("read {}: {err}", bed_path.display()))?;
    let gtf = bed12_to_gtf(&normalize_tama_names(&bed), "TAMA")?;
    finish(args, &gtf, "tama")
}

fn discover_flair(args: &DiscoverArgs) -> Result<(), String> {
    let raw = args.out.join("raw");
    fs::create_dir_all(&raw).map_err(|err| format!("create {}: {err}", raw.display()))?;
    let bed = alignments_to_bed12(&args.bam, true)?;
    fs::write(raw.join("flair.bed"), &bed).map_err(|err| format!("write BED: {err}"))?;
    let guide = match &args.gtf {
        Some(path) => Some(
            fs::read_to_string(path).map_err(|err| format!("read {}: {err}", path.display()))?,
        ),
        None => None,
    };
    let collapsed = collapse_precise_bed(&bed, guide.as_deref(), &PreciseSettings::default())?;
    fs::write(raw.join("flair.collapsed.bed"), &collapsed)
        .map_err(|err| format!("write collapsed BED: {err}"))?;
    let gtf = bed12_to_gtf(&collapsed, "FLAIR")?;
    finish(args, &gtf, "flair")
}

fn discover_isoquant(args: &DiscoverArgs) -> Result<(), String> {
    let gtf = args
        .gtf
        .as_ref()
        .ok_or_else(|| "isoquant discovery requires --gtf".to_string())?;
    prepare_indexes(args)?;
    let program = tools::isoquant()?;
    let raw = args.out.join("raw/isoquant");
    let log = args.out.join("raw/isoquant.log");
    tools::run_exclusive(
        &program,
        &[
            "--reference".to_string(),
            args.genome.display().to_string(),
            "--genedb".to_string(),
            gtf.display().to_string(),
            "--complete_genedb".to_string(),
            "--bam".to_string(),
            args.bam.display().to_string(),
            "--data_type".to_string(),
            "nanopore".to_string(),
            "--output".to_string(),
            raw.display().to_string(),
            "--prefix".to_string(),
            "plena".to_string(),
            "--threads".to_string(),
            args.threads.to_string(),
            "--force".to_string(),
            "--analysis".to_string(),
            "transcript_discovery".to_string(),
            "--model_construction_strategy".to_string(),
            "all".to_string(),
        ],
        &log,
        "plena-isoquant.lock",
    )?;
    let model = pick_tool_gtf(&raw)?;
    let text =
        fs::read_to_string(&model).map_err(|err| format!("read {}: {err}", model.display()))?;
    finish(args, &text, "isoquant")
}

fn discover_bambu(args: &DiscoverArgs) -> Result<(), String> {
    let gtf = args
        .gtf
        .as_ref()
        .ok_or_else(|| "bambu discovery requires --gtf".to_string())?;
    prepare_indexes(args)?;
    let rscript = tools::rscript()?;
    let raw = args.out.join("raw");
    fs::create_dir_all(&raw).map_err(|err| format!("create {}: {err}", raw.display()))?;
    let script = raw.join("bambu_discover.R");
    fs::write(&script, BAMBU_DISCOVER)
        .map_err(|err| format!("write {}: {err}", script.display()))?;
    let produced = raw.join("bambu.gtf");
    let log = raw.join("bambu.log");
    tools::run(
        &rscript,
        &[
            "--vanilla".to_string(),
            script.display().to_string(),
            args.bam.display().to_string(),
            gtf.display().to_string(),
            args.genome.display().to_string(),
            produced.display().to_string(),
        ],
        &log,
    )?;
    let text = fs::read_to_string(&produced)
        .map_err(|err| format!("read {}: {err}", produced.display()))?;
    finish(args, &text, "bambu")
}

fn discover_stringtie(args: &DiscoverArgs) -> Result<(), String> {
    tools::ensure_bam_index(&args.bam)?;
    let program = tools::stringtie()?;
    let raw = args.out.join("raw");
    fs::create_dir_all(&raw).map_err(|err| format!("create {}: {err}", raw.display()))?;
    let produced = raw.join("stringtie.gtf");
    let mut command = vec![
        "-L".to_string(),
        "-p".to_string(),
        args.threads.to_string(),
        "-m".to_string(),
        "50".to_string(),
        "-c".to_string(),
        "0.001".to_string(),
        "-s".to_string(),
        "0.01".to_string(),
        "-f".to_string(),
        "0".to_string(),
        "-j".to_string(),
        "1".to_string(),
        "-a".to_string(),
        "1".to_string(),
        "-o".to_string(),
        produced.display().to_string(),
    ];
    if let Some(gtf) = &args.gtf {
        command.push("-G".to_string());
        command.push(gtf.display().to_string());
    }
    command.push(args.bam.display().to_string());
    tools::run(&program, &command, &raw.join("stringtie.log"))?;
    let text = fs::read_to_string(&produced)
        .map_err(|err| format!("read {}: {err}", produced.display()))?;
    finish(args, &text, "stringtie")
}

pub fn stringtie_merge(
    output: &Path,
    guide: Option<&Path>,
    gtfs: &[PathBuf],
) -> Result<(), String> {
    if gtfs.is_empty() {
        return Err("stringtie merge needs at least one GTF".to_string());
    }
    let program = tools::stringtie()?;
    let mut command = vec![
        "--merge".to_string(),
        "-o".to_string(),
        output.display().to_string(),
        "-m".to_string(),
        "50".to_string(),
    ];
    if let Some(guide) = guide {
        command.push("-G".to_string());
        command.push(guide.display().to_string());
    }
    for gtf in gtfs {
        command.push(gtf.display().to_string());
    }
    let log = output.with_extension("log");
    tools::run(&program, &command, &log)?;
    let text =
        fs::read_to_string(output).map_err(|err| format!("read {}: {err}", output.display()))?;
    if !tools::gtf_has_transcript(&text) {
        return Err(format!(
            "StringTie merge wrote no transcript features to {}",
            output.display()
        ));
    }
    publish("merge", "stringtie-merge.gtf", output);
    publish("merge", "stringtie-merge.log", &log);
    Ok(())
}

fn prepare_indexes(args: &DiscoverArgs) -> Result<(), String> {
    tools::ensure_bam_index(&args.bam)?;
    tools::ensure_fasta_index(&args.genome)
}

fn finish(args: &DiscoverArgs, gtf: &str, algo: &str) -> Result<(), String> {
    if !tools::gtf_has_transcript(gtf) {
        return Err(format!(
            "{algo} produced no transcript features in {}",
            args.out.display()
        ));
    }
    fs::create_dir_all(&args.out).map_err(|err| format!("create {}: {err}", args.out.display()))?;
    let path = args.out.join("transcripts.gtf");
    fs::write(&path, gtf).map_err(|err| format!("write {}: {err}", path.display()))?;
    let log = args.out.join("command.log");
    append_log(
        &log,
        &format!(
            "discover --algo {algo} --bam {} --genome {}",
            args.bam.display(),
            args.genome.display()
        ),
    );
    publish("discover", &format!("{algo}-transcripts.gtf"), &path);
    publish("discover", &format!("{algo}-command.log"), &log);
    Ok(())
}

fn normalize_tama_names(bed: &str) -> String {
    let mut lines = Vec::new();
    for line in bed.lines() {
        if line.is_empty() {
            continue;
        }
        let mut columns: Vec<String> = line.split('\t').map(str::to_string).collect();
        if columns.len() > 3 {
            if let Some(index) = columns[3].find(';') {
                columns[3].replace_range(index..=index, "_");
            }
        }
        lines.push(columns.join("\t"));
    }
    if lines.is_empty() {
        String::new()
    } else {
        lines.join("\n") + "\n"
    }
}

fn pick_tool_gtf(dir: &Path) -> Result<PathBuf, String> {
    let mut found = tools::find_files(dir, |name| name.ends_with(".gtf"));
    found.sort_by_key(|path| {
        let name = path
            .file_name()
            .and_then(|text| text.to_str())
            .unwrap_or("");
        if name.contains("transcript_models") {
            0
        } else if name.contains("extended_annotation") {
            1
        } else {
            2
        }
    });
    for path in &found {
        let text = fs::read_to_string(path).unwrap_or_default();
        if tools::gtf_has_transcript(&text) {
            return Ok(path.clone());
        }
    }
    Err(format!(
        "IsoQuant wrote no transcript features under {}",
        dir.display()
    ))
}
