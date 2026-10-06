//! Construction, quantification, and StringTie merge through the shipped CLI.
//!
//! Discovery and the large UMI tests share one coordinate-sorted genome BAM:
//! 101 distinct UMIs, with the first repeated. Bambu quantification is also
//! run on three reads, two of which share a cell barcode and a UMI.

use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::{Mutex, OnceLock};

struct WorkDir(PathBuf);

impl WorkDir {
    fn new() -> Self {
        static NEXT: Mutex<u64> = Mutex::new(0);
        let mut next = NEXT.lock().unwrap();
        let id = *next;
        *next += 1;
        drop(next);
        let path =
            std::env::temp_dir().join(format!("braid-cli-stage-{}-{}", std::process::id(), id));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for WorkDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

struct Fixture {
    _dir: WorkDir,
    bam: PathBuf,
    genome: PathBuf,
    gtf: PathBuf,
    unique_umis: usize,
    read_count: usize,
}

fn shared() -> &'static Fixture {
    static FIXTURE: OnceLock<Fixture> = OnceLock::new();
    FIXTURE.get_or_init(build_fixture)
}

fn build_fixture() -> Fixture {
    let dir = WorkDir::new();
    let genome_seq = patterned_genome(1000);
    let genome = dir.0.join("genome.fa");
    fs::write(&genome, format!(">chr1\n{genome_seq}\n")).unwrap();
    let gtf = dir.0.join("genes.gtf");
    fs::write(
        &gtf,
        "\
chr1\ttest\tgene\t1\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\";
chr1\ttest\ttranscript\t1\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"T1\";
chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"T1\";
chr1\ttest\texon\t201\t300\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"T1\";
chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"T1\";
",
    )
    .unwrap();
    let seq = exon_seq(&genome_seq);
    assert_eq!(
        seq.len(),
        300,
        "transcript sequence must be longer than 200 bp"
    );
    let qual = "I".repeat(seq.len());
    // 101 distinct UMIs, with the first repeated. Bambu-Clump deduplicates only
    // after it sees more than 100 UMIs. Blocks of four identical bases keep
    // every pair at Hamming distance 4 or more, outside IsoQuant tenX_v3
    // edit distance 3.
    let distinct = block_umis(101);
    let mut umis = vec![distinct[0].clone()];
    umis.extend(distinct);
    let mut molecules = HashSet::new();
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n");
    for (index, umi) in umis.iter().enumerate() {
        molecules.insert(umi.clone());
        sam.push_str(&format!(
            "r{index}\t0\tchr1\t1\t60\t100M100N100M100N100M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:CELL1\tUB:Z:{umi}\tXS:A:+\n"
        ));
    }
    assert!(
        molecules.len() > 100 && molecules.len() < umis.len(),
        "fixture must duplicate one UMI and clear Bambu's 100-UMI dedup threshold"
    );
    let sam_path = dir.0.join("reads.sam");
    fs::write(&sam_path, sam).unwrap();
    let bam = dir.0.join("reads.bam");
    assert_cmd(
        "samtools",
        &[
            "sort",
            "-o",
            &bam.to_string_lossy(),
            &sam_path.to_string_lossy(),
        ],
    );
    assert_cmd("samtools", &["index", &bam.to_string_lossy()]);
    assert_cmd("samtools", &["faidx", &genome.to_string_lossy()]);
    Fixture {
        _dir: dir,
        bam,
        genome,
        gtf,
        unique_umis: molecules.len(),
        read_count: umis.len(),
    }
}

fn block_umis(n: usize) -> Vec<String> {
    let alphabet = ['A', 'C', 'G', 'T'];
    (0..n)
        .map(|code| {
            let mut umi = String::new();
            let mut value = code;
            for _ in 0..4 {
                let digit = value % 4;
                value /= 4;
                for _ in 0..4 {
                    umi.push(alphabet[digit]);
                }
            }
            umi
        })
        .collect()
}

fn patterned_genome(length: usize) -> String {
    let pattern = b"ACGTACGTACGT";
    let mut bases = vec![0u8; length];
    for (index, base) in bases.iter_mut().enumerate() {
        *base = pattern[index % pattern.len()];
    }
    // GT-AG donors and acceptors at the two annotated introns.
    let splice = [
        (100, b'G'),
        (101, b'T'),
        (198, b'A'),
        (199, b'G'),
        (300, b'G'),
        (301, b'T'),
        (398, b'A'),
        (399, b'G'),
    ];
    for (index, base) in splice {
        bases[index] = base;
    }
    String::from_utf8(bases).unwrap()
}

fn exon_seq(genome: &str) -> String {
    let bytes = genome.as_bytes();
    let mut seq = String::new();
    seq.push_str(std::str::from_utf8(&bytes[0..100]).unwrap());
    seq.push_str(std::str::from_utf8(&bytes[200..300]).unwrap());
    seq.push_str(std::str::from_utf8(&bytes[400..500]).unwrap());
    seq
}

fn assert_cmd(program: &str, args: &[&str]) {
    let output = Command::new(program)
        .args(args)
        .output()
        .unwrap_or_else(|err| {
            panic!("{program} failed to start: {err}");
        });
    assert!(
        output.status.success(),
        "{program} {args:?}\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn braid(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_braid"))
        .args(args)
        .output()
        .unwrap()
}

fn assert_success(output: &Output, context: &Path) {
    if output.status.success() {
        return;
    }
    let mut logs = String::new();
    if context.is_dir() {
        let mut files = Vec::new();
        collect_logs(context, &mut files);
        for path in files {
            if let Ok(text) = fs::read_to_string(&path) {
                logs.push_str(&format!("\n--- {} ---\n{text}", path.display()));
            }
        }
    }
    panic!(
        "braid exited {:?}\nstdout:\n{}\nstderr:\n{}{logs}",
        output.status.code(),
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
}

fn collect_logs(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect_logs(&path, out);
        } else if path
            .file_name()
            .and_then(|name| name.to_str())
            .is_some_and(|name| name.ends_with(".log") || name.ends_with(".txt"))
        {
            out.push(path);
        }
    }
}

fn transcript_features(gtf: &str) -> Vec<String> {
    gtf.lines()
        .filter(|line| {
            let mut columns = line.split('\t');
            columns.next();
            columns.next();
            columns.next() == Some("transcript")
        })
        .map(str::to_string)
        .collect()
}

fn discover(algo: &str, out: &Path) {
    let fix = shared();
    let output = braid(&[
        "discover",
        "--algo",
        algo,
        "--bam",
        &fix.bam.to_string_lossy(),
        "--genome",
        &fix.genome.to_string_lossy(),
        "--gtf",
        &fix.gtf.to_string_lossy(),
        "--out",
        &out.to_string_lossy(),
        "--threads",
        "1",
    ]);
    assert_success(&output, out);
    let gtf = fs::read_to_string(out.join("transcripts.gtf")).unwrap();
    assert!(
        !transcript_features(&gtf).is_empty(),
        "{algo} transcripts.gtf has no transcript feature:\n{gtf}"
    );
}

fn quant(algo: &str, out: &Path) {
    let fix = shared();
    let output = braid(&[
        "quant",
        "--algo",
        algo,
        "--bam",
        &fix.bam.to_string_lossy(),
        "--genome",
        &fix.genome.to_string_lossy(),
        "--gtf",
        &fix.gtf.to_string_lossy(),
        "--out",
        &out.to_string_lossy(),
        "--threads",
        "1",
    ]);
    assert_success(&output, out);
    let counts = fs::read_to_string(out.join("counts.tsv")).unwrap();
    let mut lines = counts.lines();
    assert_eq!(
        lines.next(),
        Some("barcode\ttranscript_id\tcount"),
        "{algo} counts:\n{counts}"
    );
    let mut sum = 0.0;
    let mut saw_cell = false;
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 3, "{algo} row {line}");
        if fields[0] == "CELL1" {
            saw_cell = true;
            sum += fields[2].parse::<f64>().unwrap();
        }
    }
    assert!(saw_cell, "{algo} matrix is not keyed by CELL1:\n{counts}");
    let expected = fix.unique_umis as f64;
    assert!(
        (sum - expected).abs() < 0.05,
        "{algo} counted {sum} for CELL1; unique CB+UMI pairs are {expected} and reads are {}",
        fix.read_count
    );
}

fn fresh(name: &str) -> PathBuf {
    let path = shared()._dir.0.join(name);
    let _ = fs::remove_dir_all(&path);
    fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn help_lists_stage_commands_without_dropping_the_existing_ones() {
    let output = braid(&["--help"]);
    assert!(output.status.success());
    let help = String::from_utf8(output.stdout).unwrap();
    for command in [
        "discover",
        "quant",
        "stringtie-merge",
        "tama-merge",
        "merge",
        "scotch",
        "scotch-dtu",
        "taco",
    ] {
        assert!(
            help.contains(&format!("  {command}")),
            "help is missing {command}:\n{help}"
        );
    }
    assert!(!help.contains("not implemented"));
    let rejected = braid(&[
        "merge",
        "--algo",
        "stringtie",
        "-f",
        "unused",
        "-p",
        "unused",
    ]);
    assert_eq!(rejected.status.code(), Some(2));
    let stderr = String::from_utf8_lossy(&rejected.stderr);
    assert!(stderr.contains("unknown merge algorithm"));
    assert!(!stderr.contains("not implemented"));
}

#[test]
fn discover_tama_emits_a_transcript() {
    discover("tama", &fresh("discover-tama"));
}

#[test]
fn discover_flair_emits_a_transcript_and_the_repeat_matches() {
    let first = fresh("discover-flair");
    let second = fresh("discover-flair-repeat");
    discover("flair", &first);
    discover("flair", &second);
    let left = fs::read(first.join("transcripts.gtf")).unwrap();
    let right = fs::read(second.join("transcripts.gtf")).unwrap();
    assert_eq!(left, right);
}

#[test]
fn discover_isoquant_emits_a_transcript() {
    discover("isoquant", &fresh("discover-isoquant"));
}

#[test]
fn discover_bambu_emits_a_transcript() {
    discover("bambu", &fresh("discover-bambu"));
}

#[test]
fn discover_stringtie_emits_a_transcript() {
    discover("stringtie", &fresh("discover-stringtie"));
}

#[test]
fn quant_scotch_counts_each_umi_once() {
    quant("scotch", &fresh("quant-scotch"));
}

#[test]
fn quant_isoquant_counts_each_umi_once_and_the_repeat_matches() {
    let first = fresh("quant-isoquant");
    let second = fresh("quant-isoquant-repeat");
    quant("isoquant", &first);
    quant("isoquant", &second);
    assert_eq!(
        fs::read(first.join("counts.tsv")).unwrap(),
        fs::read(second.join("counts.tsv")).unwrap()
    );
}

#[test]
fn quant_isosceles_counts_each_umi_once() {
    quant("isosceles", &fresh("quant-isosceles"));
}

#[test]
fn quant_bambu_counts_each_umi_once() {
    quant("bambu", &fresh("quant-bambu"));
}

#[test]
fn quant_bambu_counts_a_duplicate_umi_on_three_reads() {
    let dir = WorkDir::new();
    let genome_seq = patterned_genome(1000);
    let genome = dir.0.join("genome.fa");
    fs::write(&genome, format!(">chr1\n{genome_seq}\n")).unwrap();
    let gtf = dir.0.join("genes.gtf");
    fs::write(
        &gtf,
        "\
chr1\ttest\tgene\t1\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\";
chr1\ttest\ttranscript\t1\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"T1\";
chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"T1\";
chr1\ttest\texon\t201\t300\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"T1\";
chr1\ttest\texon\t401\t500\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\"; transcript_id \"T1\";
",
    )
    .unwrap();
    let seq = exon_seq(&genome_seq);
    let qual = "I".repeat(seq.len());
    let umis = ["AAAA", "AAAA", "CCCC"];
    let mut sam = String::from("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n");
    for (index, umi) in umis.iter().enumerate() {
        sam.push_str(&format!(
            "r{index}\t0\tchr1\t1\t60\t100M100N100M100N100M\t*\t0\t0\t{seq}\t{qual}\tCB:Z:CELL1\tUB:Z:{umi}\tXS:A:+\n"
        ));
    }
    let sam_path = dir.0.join("reads.sam");
    fs::write(&sam_path, sam).unwrap();
    let bam = dir.0.join("reads.bam");
    assert_cmd(
        "samtools",
        &[
            "sort",
            "-o",
            &bam.to_string_lossy(),
            &sam_path.to_string_lossy(),
        ],
    );
    assert_cmd("samtools", &["index", &bam.to_string_lossy()]);
    assert_cmd("samtools", &["faidx", &genome.to_string_lossy()]);
    let out = dir.0.join("quant");
    let output = braid(&[
        "quant",
        "--algo",
        "bambu",
        "--bam",
        &bam.to_string_lossy(),
        "--genome",
        &genome.to_string_lossy(),
        "--gtf",
        &gtf.to_string_lossy(),
        "--out",
        &out.to_string_lossy(),
    ]);
    assert_success(&output, &out);
    let counts = fs::read_to_string(out.join("counts.tsv")).unwrap();
    let mut lines = counts.lines();
    assert_eq!(
        lines.next(),
        Some("barcode\ttranscript_id\tcount"),
        "{counts}"
    );
    let mut sum = 0.0;
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 3, "{line}");
        assert_eq!(fields[0], "CELL1", "{counts}");
        sum += fields[2].parse::<f64>().unwrap();
    }
    assert!(
        (sum - 2.0).abs() < 0.05,
        "three reads and two UMIs counted {sum}; the repeated UMI must contribute one molecule:\n{counts}"
    );
    let log = fs::read_to_string(out.join("command.log")).unwrap();
    assert!(
        log.contains("kept 2 alignments"),
        "dedup log does not record two molecules:\n{log}"
    );
}

#[test]
fn quant_scotch_writes_the_unannotated_isoform() {
    let dir = WorkDir::new();
    let genome = dir.0.join("genome.fa");
    fs::write(&genome, format!(">chr1\n{}\n", "ACGT".repeat(250))).unwrap();
    let gtf = dir.0.join("genes.gtf");
    fs::write(
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
    sam.push_str(&sam_record("SHORT", "100M", 101, "FULL", "U1"));
    sam.push_str(&sam_record(
        "FULL",
        "100M100N100M100N100M",
        101,
        "FULL",
        "U1",
    ));
    sam.push_str(&sam_record("HEAD", "100M100N100M", 101, "HEAD", "U1"));
    for index in 0..10 {
        sam.push_str(&sam_record(
            &format!("N{index}"),
            "100M300N100M",
            101,
            &format!("N{index}"),
            "U",
        ));
    }
    let sam_path = dir.0.join("cells.sam");
    fs::write(&sam_path, sam).unwrap();
    let out = dir.0.join("quant");
    let output = braid(&[
        "quant",
        "--algo",
        "scotch",
        "--bam",
        &sam_path.to_string_lossy(),
        "--genome",
        &genome.to_string_lossy(),
        "--gtf",
        &gtf.to_string_lossy(),
        "--out",
        &out.to_string_lossy(),
    ]);
    assert_success(&output, &out);
    let annotation = fs::read_to_string(out.join("annotation.gtf")).unwrap();
    assert!(
        annotation.contains("transcript_id \"novelIsoform_5\""),
        "SCOTCH annotation has no novel transcript:\n{annotation}"
    );
    let counts = fs::read_to_string(out.join("counts.tsv")).unwrap();
    assert!(
        counts.contains("novelIsoform_5"),
        "SCOTCH counts dropped the novel transcript:\n{counts}"
    );
}

fn sam_record(qname: &str, cigar: &str, pos: i64, cb: &str, ub: &str) -> String {
    let length = cigar_query_len(cigar);
    format!(
        "{qname}\t0\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{}\t*\tCB:Z:{cb}\tUB:Z:{ub}\n",
        "ACGT".repeat(length / 4)
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
fn stringtie_merge_keeps_models_from_both_inputs() {
    let dir = WorkDir::new();
    let left = dir.0.join("left.gtf");
    let right = dir.0.join("right.gtf");
    fs::write(&left, locus_gtf("chr1", 1, "A1")).unwrap();
    fs::write(&right, locus_gtf("chr1", 8001, "B1")).unwrap();
    let merged = dir.0.join("merged.gtf");
    let output = braid(&[
        "stringtie-merge",
        "-o",
        &merged.to_string_lossy(),
        &left.to_string_lossy(),
        &right.to_string_lossy(),
    ]);
    assert_success(&output, &dir.0);
    let text = fs::read_to_string(&merged).unwrap();
    let transcripts = transcript_features(&text);
    assert!(transcripts.len() >= 2, "merged annotation:\n{text}");
    let starts: Vec<i64> = transcripts
        .iter()
        .map(|line| line.split('\t').nth(3).unwrap().parse().unwrap())
        .collect();
    assert!(
        starts.iter().any(|start| *start < 1000) && starts.iter().any(|start| *start > 7000),
        "merged starts {starts:?} do not cover both input loci:\n{text}"
    );
}

fn locus_gtf(chrom: &str, start: i64, transcript: &str) -> String {
    let exon = |offset: i64| {
        format!(
            "{chrom}\tstringtie\texon\t{}\t{}\t.\t+\t.\tgene_id \"{transcript}\"; transcript_id \"{transcript}\";\n",
            start + offset,
            start + offset + 99
        )
    };
    format!(
        "{chrom}\tstringtie\ttranscript\t{start}\t{}\t1000\t+\t.\tgene_id \"{transcript}\"; transcript_id \"{transcript}\";\n{}{}{}",
        start + 499,
        exon(0),
        exon(200),
        exon(400)
    )
}
