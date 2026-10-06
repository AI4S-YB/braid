use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

struct WorkDir(PathBuf);

impl WorkDir {
    fn new() -> Self {
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "plena-cli-compat-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for WorkDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn parity() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/parity")
}

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_plena"))
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn tama_legacy_and_long_flags_write_all_oracle_outputs() {
    let golden = parity().join("gmap_collapse");
    for dash in ["-", "--"] {
        let dir = WorkDir::new();
        let prefix = dir.0.join("gmap");
        let sam = golden.join("gmap_test.sam");
        let fasta = golden.join("test_genome.fa");
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_plena"));
        cmd.args(["tama", "collapse"])
            .arg("-s")
            .arg(&sam)
            .arg("-f")
            .arg(&fasta)
            .arg("-p")
            .arg(&prefix);
        for (key, value) in [
            ("icm", "ident_cov"),
            ("sj", "no_priority"),
            ("sjt", "10"),
            ("lde", "1000"),
            ("ses", "_"),
            ("log", "log_off"),
            ("rm", "original"),
            ("vc", "9"),
        ] {
            cmd.arg(format!("{dash}{key}")).arg(value);
        }
        let output = cmd.output().unwrap();
        assert!(
            output.status.success(),
            "{dash}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(fs::read_dir(&dir.0).unwrap().count(), 10);
        for suffix in [
            ".bed",
            "_read.txt",
            "_trans_report.txt",
            "_trans_read.bed",
            "_polya.txt",
            "_strand_check.txt",
            "_local_density_error.txt",
            "_variants.txt",
            "_varcov.txt",
            "_report.txt",
        ] {
            let expected = fs::read_to_string(golden.join(format!("gmap{suffix}")))
                .unwrap()
                .replace("test_files/gmap_test.sam", &sam.display().to_string())
                .replace("test_files/test_genome.fa", &fasta.display().to_string())
                .replace(
                    "/work/plena/tests/parity/gmap_collapse/gmap",
                    &prefix.display().to_string(),
                );
            let actual = fs::read_to_string(dir.0.join(format!("gmap{suffix}"))).unwrap();
            assert_eq!(actual, expected, "{dash}: {suffix}");
        }
    }
}

#[test]
fn tama_version_does_not_require_input_files() {
    let output = run(&["tama", "collapse", "-v", "1"]);
    assert!(output.status.success());
    let stdout = String::from_utf8(output.stdout)
        .unwrap()
        .replace("\r\n", "\n");
    assert_eq!(stdout, "tc_version_date_2023_03_28\nProgram did not run\n");
    assert!(output.stderr.is_empty());
}

#[test]
fn tama_execution_still_requires_inputs() {
    for args in [
        vec!["tama", "collapse"],
        vec!["tama", "collapse", "-s", "missing.sam"],
    ] {
        let output = run(&args);
        assert_eq!(output.status.code(), Some(2));
        assert!(String::from_utf8_lossy(&output.stderr).contains("required arguments"));
    }
}

#[test]
fn collapse_preflight_errors_preserve_existing_outputs() {
    let golden = parity().join("gmap_collapse");
    for mode in ["original", "low_mem"] {
        for failure in ["cap", "sam", "fasta"] {
            let dir = WorkDir::new();
            let prefix = dir.0.join("out");
            let existing = [dir.0.join("out.bed"), dir.0.join("out_read.txt")];
            for path in &existing {
                fs::write(path, "previous successful result\n").unwrap();
            }
            let sam = if failure == "sam" {
                dir.0.join("missing.sam")
            } else {
                golden.join("gmap_test.sam")
            };
            let fasta = if failure == "fasta" {
                dir.0.join("missing.fa")
            } else {
                golden.join("test_genome.fa")
            };
            let output = Command::new(env!("CARGO_BIN_EXE_plena"))
                .args(["tama", "collapse", "--rm", mode, "-x"])
                .arg(if failure == "cap" {
                    "invalid"
                } else {
                    "capped"
                })
                .arg("-s")
                .arg(sam)
                .arg("-f")
                .arg(fasta)
                .arg("-p")
                .arg(prefix)
                .output()
                .unwrap();
            assert_eq!(output.status.code(), Some(1), "{mode}/{failure}");
            let stderr = String::from_utf8_lossy(&output.stderr);
            let expected = match failure {
                "cap" => "-x must be capped or no_cap",
                "sam" => "missing.sam",
                _ => "missing.fa",
            };
            assert!(stderr.contains(expected), "{mode}/{failure}: {stderr}");
            for path in &existing {
                assert_eq!(
                    fs::read_to_string(path).unwrap(),
                    "previous successful result\n",
                    "{mode}/{failure}: {}",
                    path.display(),
                );
            }
            assert_eq!(fs::read_dir(&dir.0).unwrap().count(), existing.len());
        }
    }
}

#[cfg(unix)]
#[test]
fn flair_output_preserves_non_utf8_paths() {
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt;

    let dir = WorkDir::new();
    let input = parity().join("flair_combine/official/in");
    let manifest = dir.0.join("manifest.tsv");
    fs::write(
        &manifest,
        format!(
            "A1\tisoforms\t{}\t{}\t{}\n",
            input.join("collapse.isoforms.bed").display(),
            input.join("collapse.isoforms.fa").display(),
            input.join("collapse.isoform.read.map.txt").display(),
        ),
    )
    .unwrap();
    let ascii_dir = dir.0.join("ascii");
    let native_dir = dir.0.join(OsString::from_vec(b"native-\xff".to_vec()));
    for output_dir in [&ascii_dir, &native_dir] {
        let output = Command::new(env!("CARGO_BIN_EXE_plena"))
            .args(["flair", "combine", "-c", "-m"])
            .arg(&manifest)
            .arg("-o")
            .arg(output_dir.join("combined"))
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    for suffix in [".bed", ".counts.tsv", ".isoform.map.txt", ".fa", ".gtf"] {
        let name = format!("combined{suffix}");
        assert_eq!(
            fs::read(native_dir.join(&name)).unwrap(),
            fs::read(ascii_dir.join(&name)).unwrap(),
            "{name}",
        );
    }
    assert_eq!(fs::read_dir(&native_dir).unwrap().count(), 5);
}

#[test]
fn merge_is_available_and_legacy_cds_reaches_input_validation() {
    let output = run(&[
        "tama", "merge", "-f", "unused", "-p", "unused", "-cds", "source",
    ]);
    assert_eq!(output.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&output.stderr).contains("read unused"));
    let output = run(&["--help"]);
    let help = String::from_utf8(output.stdout).unwrap();
    assert!(!help.contains("not implemented"));
    assert!(help.contains("  tama "));
    let tama_help = String::from_utf8(run(&["tama", "--help"]).stdout).unwrap();
    assert!(tama_help.contains("  merge "));
    assert!(tama_help.contains("  collapse "));
}

#[test]
fn scotch_quantifies_a_sam_group_and_help_lists_the_command() {
    let help = run(&["--help"]);
    let text = String::from_utf8(help.stdout).unwrap();
    assert!(text.contains("  scotch "));
    let scotch_help = String::from_utf8(run(&["scotch", "--help"]).stdout).unwrap();
    assert!(scotch_help.contains("  quant "));
    assert!(scotch_help.contains("  dtu "));
    assert!(!text.contains("not implemented"));
    let dir = WorkDir::new();
    let gtf = dir.0.join("genes.gtf");
    fs::write(
        &gtf,
        "\
chr1\t.\tgene\t101\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\";\n\
chr1\t.\ttranscript\t101\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ONLY\";\n\
chr1\t.\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"ONLY\";\n",
    )
    .unwrap();
    let sam = dir.0.join("cells.sam");
    fs::write(
        &sam,
        format!(
            "R1\t0\tchr1\t101\t60\t100M\t*\t0\t0\t{}\t*\tCB:Z:C1\tUB:Z:U1\n",
            "ACGT".repeat(25)
        ),
    )
    .unwrap();
    let out = dir.0.join("scotch");
    let output = Command::new(env!("CARGO_BIN_EXE_plena"))
        .args(["scotch", "quant"])
        .arg("--bam")
        .arg(&sam)
        .arg("--gtf")
        .arg(&gtf)
        .arg("--out")
        .arg(&out)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let counts = fs::read_to_string(out.join("cells/count_matrix/gene_counts.csv")).unwrap();
    assert!(counts.contains("C1:sample0"));
    assert!(counts.contains("GENE1"));
}

#[test]
fn taco_cli_assembles_one_sample_and_rejects_a_second_run() {
    let dir = WorkDir::new();
    let gtf = dir.0.join("a.gtf");
    fs::write(
        &gtf,
        "\
chr1\tA\ttranscript\t1000\t2000\t.\t+\t.\ttranscript_id \"t1\"; FPKM \"4\";\n\
chr1\tA\texon\t1000\t1200\t.\t+\t.\ttranscript_id \"t1\";\n\
chr1\tA\texon\t1400\t1600\t.\t+\t.\ttranscript_id \"t1\";\n\
chr1\tA\texon\t1800\t2000\t.\t+\t.\ttranscript_id \"t1\";\n",
    )
    .unwrap();
    let samples = dir.0.join("samples.tsv");
    fs::write(&samples, format!("{}\n", gtf.display())).unwrap();
    let out = dir.0.join("taco");
    let output = Command::new(env!("CARGO_BIN_EXE_plena"))
        .arg("taco")
        .arg("-o")
        .arg(&out)
        .arg(&samples)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read_dir(&out).unwrap().count(), 14);
    let assembly = fs::read_to_string(out.join("assembly.gtf")).unwrap();
    assert_eq!(
        assembly,
        "\
chr1\ttaco\ttranscript\t1000\t2000\t1000\t+\t.\texpr \"1000000.000\"; rel_frac \"1.00000\"; abs_frac \"1.00000\"; locus_id \"L1\"; gene_id \"G1\"; tss_id \"TSS1\"; transcript_id \"TU1\";\n\
chr1\ttaco\texon\t1000\t1200\t1000\t+\t.\tlocus_id \"L1\"; gene_id \"G1\"; tss_id \"TSS1\"; transcript_id \"TU1\";\n\
chr1\ttaco\texon\t1400\t1600\t1000\t+\t.\tlocus_id \"L1\"; gene_id \"G1\"; tss_id \"TSS1\"; transcript_id \"TU1\";\n\
chr1\ttaco\texon\t1800\t2000\t1000\t+\t.\tlocus_id \"L1\"; gene_id \"G1\"; tss_id \"TSS1\"; transcript_id \"TU1\";\n"
    );
    let again = Command::new(env!("CARGO_BIN_EXE_plena"))
        .arg("taco")
        .arg("-o")
        .arg(&out)
        .arg(&samples)
        .output()
        .unwrap();
    assert_eq!(again.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&again.stderr).contains("already exists"));

    let missing = Command::new(env!("CARGO_BIN_EXE_plena"))
        .arg("taco")
        .arg("-o")
        .arg(dir.0.join("other"))
        .arg(dir.0.join("missing.tsv"))
        .output()
        .unwrap();
    assert_eq!(missing.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&missing.stderr).contains("not found"));

    let help = run(&["--help"]);
    assert!(help.status.success());
    assert!(String::from_utf8_lossy(&help.stdout).contains("taco"));
    let taco_help = run(&["taco", "--help"]);
    assert!(taco_help.status.success());
    let text = String::from_utf8(taco_help.stdout).unwrap();
    assert!(text.contains("one process"));
    assert!(text.contains("--resume"));
    assert!(text.contains("--assemble"));
}

#[test]
fn flair_cli_writes_filtered_chain_regression_outputs() {
    let golden = parity().join("flair_combine/regressions/consecutive-filtered-chains");
    let dir = WorkDir::new();
    let manifest = dir.0.join("manifest.tsv");
    fs::write(
        &manifest,
        format!(
            "S\tisoforms\t{}\t\t{}\n",
            golden.join("in.bed").display(),
            golden.join("in.map").display()
        ),
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_plena"))
        .args(["flair", "combine"])
        .arg("-m")
        .arg(manifest)
        .arg("-o")
        .arg(dir.0.join("combined"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(fs::read_dir(&dir.0).unwrap().count(), 4);
    for suffix in [".bed", ".counts.tsv", ".isoform.map.txt"] {
        assert_eq!(
            fs::read(dir.0.join(format!("combined{suffix}"))).unwrap(),
            fs::read(golden.join(format!("combined{suffix}"))).unwrap()
        );
    }
}

fn fixture_args(case: &Path, dir: &Path, kind: &str, bam: bool) -> Vec<String> {
    let mut args: Vec<String> = fs::read_to_string(case.join("args.txt"))
        .unwrap()
        .lines()
        .map(str::to_string)
        .collect();
    if kind == "merge" {
        let manifest = fs::read_to_string(case.join("files.tsv"))
            .unwrap()
            .lines()
            .map(|line| {
                let (path, rest) = line.split_once('\t').unwrap();
                format!("{}\t{rest}\n", case.join(path).display())
            })
            .collect::<String>();
        let path = dir.join("files.tsv");
        fs::write(&path, manifest).unwrap();
        let index = args.iter().position(|s| s == "-f").unwrap();
        args[index + 1] = path.display().to_string();
    } else {
        for flag in ["-s", "-f"] {
            let index = args.iter().position(|s| s == flag).unwrap();
            args[index + 1] = if bam && flag == "-s" {
                case.parent()
                    .unwrap()
                    .join("gmap.bam")
                    .display()
                    .to_string()
            } else {
                case.join(&args[index + 1]).display().to_string()
            };
        }
        if bam {
            let index = args.iter().position(|s| s == "-b").unwrap();
            args[index + 1] = "BAM".into();
        }
    }
    args.extend(["-p".into(), dir.join("out").display().to_string()]);
    args
}

fn check_fixture(name: &str, kind: &str, success: bool, bam: bool) {
    let case = parity().join("tama_modes").join(name);
    let dir = WorkDir::new();
    let args = fixture_args(&case, &dir.0, kind, bam);
    let output = Command::new(env!("CARGO_BIN_EXE_plena"))
        .args(["tama", kind])
        .args(&args)
        .output()
        .unwrap();
    assert_eq!(
        output.status.success(),
        success,
        "{name} bam={bam}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    if !success {
        return;
    }
    let value = |flag: &str| &args[args.iter().position(|s| s == flag).unwrap() + 1];
    let mut expected_names = std::collections::BTreeSet::new();
    for file in fs::read_dir(&case).unwrap() {
        let file = file.unwrap();
        let name = file.file_name();
        if !name.to_string_lossy().starts_with("out") {
            continue;
        }
        expected_names.insert(name.clone());
        let mut expected = fs::read_to_string(file.path()).unwrap();
        if kind == "collapse" {
            expected = expected
                .replace("@SAM@", value("-s"))
                .replace("@FASTA@", value("-f"))
                .replace("@PREFIX@", value("-p"));
            if bam {
                expected = expected.replace("-b SAM", "-b BAM");
            }
        }
        let actual = fs::read_to_string(dir.0.join(&name)).unwrap();
        assert_eq!(actual, expected, "{name:?} in {} bam={bam}", case.display());
    }
    let actual_names: std::collections::BTreeSet<_> = fs::read_dir(&dir.0)
        .unwrap()
        .map(|f| f.unwrap().file_name())
        .filter(|s| s.to_string_lossy().starts_with("out"))
        .collect();
    assert_eq!(actual_names, expected_names);
}

#[test]
fn tama_modes_match_source_derived_goldens() {
    let index = fs::read_to_string(parity().join("tama_modes/cases.tsv")).unwrap();
    for line in index.lines() {
        let fields: Vec<_> = line.split('\t').collect();
        check_fixture(fields[0], fields[1], fields[2] == "match", false);
    }
}

fn have_samtools() -> bool {
    let available = Command::new("samtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success());
    assert!(
        available || std::env::var_os("PLENA_REQUIRE_SAMTOOLS").is_none(),
        "CI requires samtools for BAM tests"
    );
    if !available {
        eprintln!("BAM integration tests skipped: install samtools or set PLENA_REQUIRE_SAMTOOLS=1 to require it");
    }
    available
}

#[test]
fn actual_bam_matches_sam_in_all_collapse_modes() {
    if !have_samtools() {
        return;
    }
    for cap in ["capped", "no_cap"] {
        for mode in ["original", "low_mem"] {
            for variant in 0..4 {
                check_fixture(
                    &format!("gmap-{cap}-{mode}-SAM-{variant}"),
                    "collapse",
                    true,
                    true,
                );
            }
        }
    }
}

#[test]
fn truncated_bam_is_not_reported_as_success() {
    if !have_samtools() {
        return;
    }
    let dir = WorkDir::new();
    let bam = dir.0.join("broken.bam");
    fs::write(&bam, b"BAM\x01not a valid BAM file").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_plena"))
        .args(["tama", "collapse"])
        .arg("-s")
        .arg(bam)
        .arg("-b")
        .arg("BAM")
        .arg("-f")
        .arg(parity().join("gmap_collapse/test_genome.fa"))
        .arg("-p")
        .arg(dir.0.join("out"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("samtools view"));
    assert!(!fs::read_to_string(dir.0.join("out_report.txt"))
        .unwrap()
        .contains("successfully"));
}

#[test]
fn merge_default_rejects_duplicate_groups() {
    let case = parity().join("tama_modes/merge-071");
    let dir = WorkDir::new();
    let mut args = fixture_args(&case, &dir.0, "merge", false);
    let index = args.iter().position(|s| s == "-d").unwrap();
    args.drain(index..index + 2);
    let output = Command::new(env!("CARGO_BIN_EXE_plena"))
        .args(["tama", "merge"])
        .args(args)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("duplicate merged models"));
}
