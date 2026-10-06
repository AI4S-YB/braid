//! Per-cell quantification from a genome BAM that already carries barcode and
//! UMI tags. Counts are written to `counts.tsv` as barcode, transcript, count.
//! Tool directories stay under `raw/`.

use std::fs;
use std::path::PathBuf;
use std::process::ExitCode;

use braid_scotch::{quantify, Platform, ScotchSettings};
use clap::Args;

use crate::dedup;
use crate::tools::{self, append_log, publish};

const BAMBU_QUANT: &str = include_str!("scripts/bambu_quant.R");
const ISOSCELES_QUANT: &str = include_str!("scripts/isosceles_quant.R");

#[derive(Args, Debug)]
pub struct QuantArgs {
    /// `scotch`, `isoquant`, `isosceles`, or `bambu`.
    #[arg(long)]
    pub algo: String,
    /// Genome BAM or, for SCOTCH, SAM. External tools require a sorted BAM.
    #[arg(long)]
    pub bam: PathBuf,
    /// Genome FASTA matching the alignments.
    #[arg(long)]
    pub genome: PathBuf,
    /// Transcript annotation used for assignment.
    #[arg(long)]
    pub gtf: PathBuf,
    /// Output directory. `counts.tsv` is written here.
    #[arg(long)]
    pub out: PathBuf,
    /// Cell-barcode BAM tag.
    #[arg(long, default_value = "CB")]
    pub barcode_tag: String,
    /// UMI BAM tag.
    #[arg(long, default_value = "UB")]
    pub umi_tag: String,
    /// Threads for external programs.
    #[arg(long, default_value_t = 1)]
    pub threads: i64,
}

pub fn run(args: QuantArgs) -> ExitCode {
    if args.threads < 1 {
        eprintln!("--threads must be at least 1");
        return ExitCode::from(2);
    }
    let result = match args.algo.as_str() {
        "scotch" => quant_scotch(&args),
        "isoquant" => quant_isoquant(&args),
        "isosceles" => quant_isosceles(&args),
        "bambu" => quant_bambu(&args),
        other => {
            eprintln!("unknown quant algorithm '{other}'");
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

fn quant_scotch(args: &QuantArgs) -> Result<(), String> {
    let raw = args.out.join("raw/scotch");
    let settings = ScotchSettings {
        bams: vec![args.bam.clone()],
        gtf: args.gtf.clone(),
        out: raw.clone(),
        fasta: Some(args.genome.clone()),
        platform: Platform::TenxOnt,
        barcode_cell: Some(args.barcode_tag.clone()),
        barcode_umi: Some(args.umi_tag.clone()),
        ..ScotchSettings::default()
    };
    quantify(&settings).map_err(|err| err.to_string())?;
    let tables = tools::find_files(&raw, |name| name == "transcript_counts.csv");
    let table = tables.first().ok_or_else(|| {
        format!(
            "SCOTCH wrote no transcript_counts.csv under {}",
            raw.display()
        )
    })?;
    let rows = scotch_rows(table)?;
    let annotation = raw.join("annotation.gtf");
    if annotation.is_file() {
        fs::copy(&annotation, args.out.join("annotation.gtf"))
            .map_err(|err| format!("copy {}: {err}", annotation.display()))?;
        publish(
            "quant",
            "scotch-annotation.gtf",
            &args.out.join("annotation.gtf"),
        );
    }
    finish(args, rows, "scotch")
}

fn scotch_rows(path: &std::path::Path) -> Result<Vec<(String, String, f64)>, String> {
    let text = fs::read_to_string(path).map_err(|err| format!("read {}: {err}", path.display()))?;
    let mut lines = text.lines();
    let header = lines
        .next()
        .ok_or_else(|| format!("{} is empty", path.display()))?;
    let columns: Vec<&str> = header.split(',').collect();
    if columns.len() < 2 {
        return Err(format!("{} has no transcript columns", path.display()));
    }
    let mut rows = Vec::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() != columns.len() {
            return Err(format!(
                "{} has a row that does not match its header",
                path.display()
            ));
        }
        let barcode = fields[0]
            .split_once(":sample")
            .map(|(barcode, _)| barcode)
            .unwrap_or(fields[0]);
        for (index, transcript) in columns.iter().enumerate().skip(1) {
            let count: f64 = fields[index].parse().map_err(|_| {
                format!(
                    "{} has a non-numeric count for {transcript}",
                    path.display()
                )
            })?;
            if count != 0.0 {
                rows.push((barcode.to_string(), (*transcript).to_string(), count));
            }
        }
    }
    Ok(rows)
}

fn quant_isoquant(args: &QuantArgs) -> Result<(), String> {
    tools::ensure_bam_index(&args.bam)?;
    tools::ensure_fasta_index(&args.genome)?;
    let program = tools::isoquant()?;
    let raw = args.out.join("raw/isoquant");
    let log = args.out.join("raw/isoquant.log");
    tools::run_exclusive(
        &program,
        &[
            "--reference".to_string(),
            args.genome.display().to_string(),
            "--genedb".to_string(),
            args.gtf.display().to_string(),
            "--complete_genedb".to_string(),
            "--bam".to_string(),
            args.bam.display().to_string(),
            "--data_type".to_string(),
            "nanopore".to_string(),
            "--output".to_string(),
            raw.display().to_string(),
            "--prefix".to_string(),
            "braid".to_string(),
            "--threads".to_string(),
            args.threads.to_string(),
            "--force".to_string(),
            "--mode".to_string(),
            "tenX_v3".to_string(),
            "--barcoded_bam".to_string(),
            "--barcode_tag".to_string(),
            args.barcode_tag.clone(),
            "--umi_tag".to_string(),
            args.umi_tag.clone(),
            "--analysis".to_string(),
            "quantification".to_string(),
        ],
        &log,
        "braid-isoquant.lock",
    )?;
    let mut tables = tools::find_files(&raw, |name| {
        // Barcoded runs write `*.transcript_grouped_barcode_counts.linear.tsv`.
        // Bulk runs write `*.transcript_grouped_counts.linear.tsv`. The linear
        // file is feature, barcode, count. The wide matrix is not.
        name.contains("transcript_grouped")
            && name.contains("counts")
            && name.ends_with(".linear.tsv")
            && !name.contains("tpm")
    });
    tables.sort_by_key(|path| {
        let name = path
            .file_name()
            .and_then(|text| text.to_str())
            .unwrap_or("");
        if name.contains("linear") {
            0
        } else {
            1
        }
    });
    let table = tables.first().ok_or_else(|| {
        format!(
            "IsoQuant wrote no transcript_grouped_counts table under {}",
            raw.display()
        )
    })?;
    let rows = tools::read_count_table(table)?;
    finish(args, rows, "isoquant")
}

fn quant_isosceles(args: &QuantArgs) -> Result<(), String> {
    tools::ensure_fasta_index(&args.genome)?;
    let raw = args.out.join("raw");
    fs::create_dir_all(&raw).map_err(|err| format!("create {}: {err}", raw.display()))?;
    let deduped = raw.join("isosceles.dedup.bam");
    let kept = dedup::dedup_bam(&args.bam, &deduped, &args.barcode_tag, &args.umi_tag)?;
    append_log(
        &args.out.join("command.log"),
        &format!("isosceles kept {kept} alignments after CB+UMI deduplication"),
    );
    let rscript = tools::rscript()?;
    let script = raw.join("isosceles_quant.R");
    fs::write(&script, ISOSCELES_QUANT)
        .map_err(|err| format!("write {}: {err}", script.display()))?;
    let produced = raw.join("isosceles.counts.tsv");
    tools::run(
        &rscript,
        &[
            "--vanilla".to_string(),
            script.display().to_string(),
            deduped.display().to_string(),
            args.gtf.display().to_string(),
            args.genome.display().to_string(),
            produced.display().to_string(),
            args.barcode_tag.clone(),
        ],
        &raw.join("isosceles.log"),
    )?;
    let rows = tools::read_count_table(&produced)?;
    finish(args, rows, "isosceles")
}

fn quant_bambu(args: &QuantArgs) -> Result<(), String> {
    tools::ensure_fasta_index(&args.genome)?;
    let rscript = tools::rscript()?;
    let raw = args.out.join("raw");
    fs::create_dir_all(&raw).map_err(|err| format!("create {}: {err}", raw.display()))?;
    // bambu.singlecell skips UMI dedup unless a chromosome has more than 100
    // distinct UMIs. Deduplicate first so a repeated CB+UMI is one molecule
    // on a small BAM as well.
    let deduped = raw.join("bambu.dedup.bam");
    let kept = dedup::dedup_bam(&args.bam, &deduped, &args.barcode_tag, &args.umi_tag)?;
    append_log(
        &args.out.join("command.log"),
        &format!("bambu kept {kept} alignments after CB+UMI deduplication"),
    );
    let script = raw.join("bambu_quant.R");
    fs::write(&script, BAMBU_QUANT).map_err(|err| format!("write {}: {err}", script.display()))?;
    let produced = raw.join("bambu.counts.tsv");
    tools::run(
        &rscript,
        &[
            "--vanilla".to_string(),
            script.display().to_string(),
            deduped.display().to_string(),
            args.gtf.display().to_string(),
            args.genome.display().to_string(),
            produced.display().to_string(),
            args.barcode_tag.clone(),
            args.umi_tag.clone(),
        ],
        &raw.join("bambu.log"),
    )?;
    let rows = tools::read_count_table(&produced)?;
    finish(args, rows, "bambu")
}

fn finish(args: &QuantArgs, rows: Vec<(String, String, f64)>, algo: &str) -> Result<(), String> {
    fs::create_dir_all(&args.out).map_err(|err| format!("create {}: {err}", args.out.display()))?;
    let path = args.out.join("counts.tsv");
    tools::write_counts(&path, &rows)?;
    let log = args.out.join("command.log");
    append_log(
        &log,
        &format!(
            "quant --algo {algo} --bam {} --barcode-tag {} --umi-tag {}",
            args.bam.display(),
            args.barcode_tag,
            args.umi_tag
        ),
    );
    publish("quant", &format!("{algo}-counts.tsv"), &path);
    publish("quant", &format!("{algo}-command.log"), &log);
    Ok(())
}
