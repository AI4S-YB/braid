//! Resolve pinned external programs and shared output helpers.
//!
//! Resolution order is an explicit environment variable, then `~/.local/bin`,
//! then the `braid-lr` conda environment. `Rscript` does not fall through to
//! `PATH`, because the system R is older than Bambu and Isosceles require.

use std::fs::{self, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;

pub fn rscript() -> Result<PathBuf, String> {
    resolve("BRAID_RSCRIPT", "Rscript", false)
}

pub fn isoquant() -> Result<PathBuf, String> {
    resolve("BRAID_ISOQUANT", "isoquant", true)
}

pub fn stringtie() -> Result<PathBuf, String> {
    resolve("BRAID_STRINGTIE", "stringtie", true)
}

pub fn samtools() -> Result<PathBuf, String> {
    resolve("BRAID_SAMTOOLS", "samtools", true)
}

fn resolve(env_key: &str, name: &str, allow_path: bool) -> Result<PathBuf, String> {
    if let Some(value) = std::env::var_os(env_key) {
        let path = PathBuf::from(&value);
        if path.is_file() {
            return Ok(path);
        }
        return Err(format!("{env_key}={} is not a file", path.display()));
    }
    let home = std::env::var_os("HOME").map(PathBuf::from);
    if let Some(home) = home {
        for path in [
            home.join(".local/bin").join(name),
            home.join("miniforge3/envs/braid-lr/bin").join(name),
        ] {
            if path.is_file() {
                return Ok(path);
            }
        }
    }
    if allow_path {
        if let Some(found) = which(name) {
            return Ok(found);
        }
    }
    Err(format!(
        "could not find {name}. Set {env_key}, or install it under ~/.local/bin or the braid-lr conda environment"
    ))
}

fn which(name: &str) -> Option<PathBuf> {
    let output = Command::new("which").arg(name).output().ok()?;
    if !output.status.success() {
        return None;
    }
    let text = String::from_utf8(output.stdout).ok()?;
    let path = PathBuf::from(text.trim());
    path.is_file().then_some(path)
}

pub fn run(program: &Path, args: &[String], log: &Path) -> Result<(), String> {
    run_command(program, args, log, &[])
}

/// Run `program` while holding an exclusive `flock` so two IsoQuant processes
/// do not rewrite `~/.config/IsoQuant/db_config.json` at the same time.
pub fn run_exclusive(
    program: &Path,
    args: &[String],
    log: &Path,
    lock_name: &str,
) -> Result<(), String> {
    let Some(flock) = which("flock") else {
        return run(program, args, log);
    };
    let lock = std::env::temp_dir().join(lock_name);
    let mut prefixed = vec![lock.display().to_string(), program.display().to_string()];
    prefixed.extend(args.iter().cloned());
    run_command(&flock, &prefixed, log, args)
}

fn run_command(
    program: &Path,
    args: &[String],
    log: &Path,
    logged_tail: &[String],
) -> Result<(), String> {
    let output = Command::new(program)
        .args(args)
        .output()
        .map_err(|err| format!("run {}: {err}", program.display()))?;
    let (shown_program, shown_args): (&Path, &[String]) = if logged_tail.is_empty() {
        (program, args)
    } else {
        (
            Path::new(args.get(1).map(String::as_str).unwrap_or("")),
            logged_tail,
        )
    };
    let mut record = format!("$ {}", shown_program.display());
    for arg in shown_args {
        record.push(' ');
        if arg.contains(' ') {
            record.push('"');
            record.push_str(arg);
            record.push('"');
        } else {
            record.push_str(arg);
        }
    }
    record.push_str(&format!("\nstatus: {}\n", output.status));
    record.push_str("--- stdout ---\n");
    record.push_str(&String::from_utf8_lossy(&output.stdout));
    record.push_str("\n--- stderr ---\n");
    record.push_str(&String::from_utf8_lossy(&output.stderr));
    record.push('\n');
    if let Some(parent) = log.parent() {
        fs::create_dir_all(parent).ok();
    }
    fs::write(log, &record).ok();
    if output.status.success() {
        Ok(())
    } else {
        Err(format!(
            "{} failed. See {}",
            shown_program.display(),
            log.display()
        ))
    }
}

pub fn ensure_bam_index(bam: &Path) -> Result<(), String> {
    let index = PathBuf::from(format!("{}.bai", bam.display()));
    if index.is_file() {
        return Ok(());
    }
    let program = samtools()?;
    let log = bam.with_file_name("samtools-index.log");
    run(
        &program,
        &["index".to_string(), bam.display().to_string()],
        &log,
    )
}

pub fn ensure_fasta_index(fasta: &Path) -> Result<(), String> {
    let index = PathBuf::from(format!("{}.fai", fasta.display()));
    if index.is_file() {
        return Ok(());
    }
    let program = samtools()?;
    let log = fasta.with_file_name("samtools-faidx.log");
    run(
        &program,
        &["faidx".to_string(), fasta.display().to_string()],
        &log,
    )
}

pub fn find_files(dir: &Path, pred: impl Fn(&str) -> bool) -> Vec<PathBuf> {
    let mut out = Vec::new();
    let mut stack = vec![dir.to_path_buf()];
    while let Some(current) = stack.pop() {
        let Ok(entries) = fs::read_dir(&current) else {
            continue;
        };
        for entry in entries.flatten() {
            let path = entry.path();
            if path.is_dir() {
                stack.push(path);
            } else if path
                .file_name()
                .and_then(|name| name.to_str())
                .is_some_and(&pred)
            {
                out.push(path);
            }
        }
    }
    out.sort();
    out
}

pub fn gtf_has_transcript(text: &str) -> bool {
    text.lines().any(|line| {
        let mut columns = line.split('\t');
        columns.next();
        columns.next();
        columns.next() == Some("transcript")
    })
}

pub fn format_count(value: f64) -> String {
    if (value - value.round()).abs() < 1e-4 {
        format!("{}", value.round() as i64)
    } else {
        let text = format!("{value:.6}");
        text.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}

pub fn write_counts(path: &Path, rows: &[(String, String, f64)]) -> Result<(), String> {
    let mut text = String::from("barcode\ttranscript_id\tcount\n");
    for (barcode, transcript, count) in rows {
        text.push_str(barcode);
        text.push('\t');
        text.push_str(transcript);
        text.push('\t');
        text.push_str(&format_count(*count));
        text.push('\n');
    }
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).map_err(|err| format!("create {}: {err}", parent.display()))?;
    }
    fs::write(path, text).map_err(|err| format!("write {}: {err}", path.display()))
}

pub fn read_count_table(path: &Path) -> Result<Vec<(String, String, f64)>, String> {
    let text = fs::read_to_string(path).map_err(|err| format!("read {}: {err}", path.display()))?;
    let mut rows = Vec::new();
    let mut header_checked = false;
    // IsoQuant grouped counts are `feature_id group_id count` (transcript, barcode).
    // The R writers and a plain table are `barcode transcript_id count`.
    let mut transcript_first = false;
    for (index, line) in text.lines().enumerate() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            return Err(format!(
                "{} line {} has {} columns",
                path.display(),
                index + 1,
                fields.len()
            ));
        }
        if !header_checked {
            header_checked = true;
            let name = fields[0].to_ascii_lowercase();
            if name == "barcode" || name == "feature_id" || name == "cell" {
                transcript_first =
                    name == "feature_id" && fields[1].eq_ignore_ascii_case("group_id");
                continue;
            }
        }
        let count: f64 = fields[2].parse().map_err(|_| {
            format!(
                "{} line {} has a non-numeric count",
                path.display(),
                index + 1
            )
        })?;
        if count == 0.0 {
            continue;
        }
        let (barcode, transcript) = if transcript_first {
            (fields[1], fields[0])
        } else {
            (fields[0], fields[1])
        };
        rows.push((barcode.to_string(), transcript.to_string(), count));
    }
    Ok(rows)
}

pub fn publish(kind: &str, name: &str, path: &Path) {
    let Ok(root) = std::env::var("BRAID_EVIDENCE_DIR") else {
        return;
    };
    let dir = PathBuf::from(root).join(kind);
    if fs::create_dir_all(&dir).is_err() || !path.is_file() {
        return;
    }
    let dest = dir.join(name);
    let dest = if dest.exists() {
        dir.join(format!("repeat-{name}"))
    } else {
        dest
    };
    fs::copy(path, dest).ok();
}

pub fn append_log(path: &Path, line: &str) {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent).ok();
    }
    if let Ok(mut file) = OpenOptions::new().create(true).append(true).open(path) {
        writeln!(file, "{line}").ok();
    }
}
