//! Synthetic inputs for the SCOTCH functions exercised by upstream.py.
//!
//! Non-DTU lines are compared with the outputs recorded from SCOTCH commit
//! 15d6ad8. Gene choice and DTU are asserted separately: upstream's pandas
//! label lookup and R's L-BFGS-B fit are not copied.

use std::fs;

use crate::assign::{diagnose, exon_base_counts, query_poly, MatchParams, Outcome};
use crate::bam::open_alignments;
use crate::coverage::cut_by_derivative;
use crate::dtu::{differential_usage, DtuSettings};
use crate::gtf::{gtf_exons, parse_gtf, Gene, Isoform};
use crate::novel::{discover, group_novel, NovelCall};

fn fmt_pct(values: &[f64]) -> String {
    values
        .iter()
        .map(|value| {
            if value.fract() == 0.0 {
                format!("{}", *value as i64)
            } else {
                format!("{value:.6}")
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn fmt_i8(values: &[i8]) -> String {
    values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn fmt_scores(values: &[f64]) -> String {
    values
        .iter()
        .map(|value| format!("{value:.6}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn write_sam(path: &std::path::Path, body: &str) {
    let mut text = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n");
    text.push_str(body);
    fs::write(path, text).unwrap();
}

fn sam_line(qname: &str, pos: i64, cigar: &str, seq: &str) -> String {
    format!("{qname}\t0\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t*\tCB:Z:C1\tUB:Z:U1\n")
}

fn seq_for(len: usize) -> String {
    "ACGT".repeat(len / 4 + 1).chars().take(len).collect()
}

fn query_len(cigar: &str) -> usize {
    let mut total = 0usize;
    let mut number = String::new();
    for char in cigar.chars() {
        if char.is_ascii_digit() {
            number.push(char);
        } else {
            let length: usize = number.parse().unwrap();
            number.clear();
            if matches!(char, 'M' | 'I' | 'S' | '=' | 'X') {
                total += length;
            }
        }
    }
    total
}

const GENE: &str = "\
chr1\t.\tgene\t101\t600\t.\t+\t.\tgene_id \"G1\"; gene_name \"GENE1\";
chr1\t.\ttranscript\t101\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ALL\";
chr1\t.\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ALL\";
chr1\t.\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ALL\";
chr1\t.\texon\t501\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"ALL\";
chr1\t.\ttranscript\t101\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"LEFT\";
chr1\t.\texon\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"LEFT\";
chr1\t.\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"LEFT\";
chr1\t.\ttranscript\t301\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"RIGHT\";
chr1\t.\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"RIGHT\";
chr1\t.\texon\t501\t600\t.\t+\t.\tgene_id \"G1\"; transcript_id \"RIGHT\";
";

fn make_gene(
    name: &str,
    id: &str,
    start: i64,
    end: i64,
    exons: Vec<(i64, i64)>,
    isoform: &str,
    used: Vec<usize>,
) -> Gene {
    Gene {
        name: name.into(),
        id: id.into(),
        chrom: "chr1".into(),
        start,
        end,
        strand: '+',
        exons,
        isoforms: vec![Isoform {
            name: isoform.into(),
            exons: used,
            order: 0,
        }],
        label: name.into(),
    }
}

fn assign_line(name: &str, gene: &Gene, counts: &[u32], poly: bool) -> String {
    let found = diagnose(gene, counts, &MatchParams::default(), poly);
    let payload = if found.kind == "novel" {
        format!(
            "pct={} vector={}",
            fmt_pct(&found.pct0),
            fmt_i8(&found.exon_vector)
        )
    } else {
        format!("vector={}", fmt_i8(&found.compatible))
    };
    format!(
        "assign {name} kind={} isoform={} scores={} {payload}",
        found.kind,
        found.isoform,
        fmt_scores(&found.scores)
    )
}

fn novel_line(label: &str, vectors: &[Vec<i8>]) -> String {
    let reads: Vec<(String, Vec<i8>)> = vectors
        .iter()
        .enumerate()
        .map(|(index, vector)| (format!("r{index}"), vector.clone()))
        .collect();
    let (isoforms, calls) = discover(&reads, 10, 50);
    let mut names: Vec<&str> = isoforms.iter().map(|iso| iso.name.as_str()).collect();
    names.sort_unstable();
    // find_novel_by_chunk returns no table when the chunk stays under 10 reads.
    if isoforms.is_empty() {
        return format!("{label} isoforms= calls= empty=");
    }
    let mut assigned = Vec::new();
    let mut empty = Vec::new();
    for (read, call) in reads.iter().zip(&calls) {
        match call {
            NovelCall::Isoform(name) => assigned.push(format!("{}:{name}", read.0)),
            NovelCall::Uncategorized => empty.push(read.0.clone()),
        }
    }
    format!(
        "{label} isoforms={} calls={} empty={}",
        names.join(","),
        assigned.join(","),
        empty.join(",")
    )
}

fn group_line(label: &str, strand: char) -> String {
    let isoforms = vec![
        ("novelIsoform_1".to_string(), vec![0usize]),
        ("novelIsoform_2".to_string(), vec![1]),
        ("novelIsoform_3".to_string(), vec![0, 1]),
    ];
    let mut pairs: Vec<String> = group_novel(&isoforms, strand)
        .into_iter()
        .map(|(name, kept)| format!("{name}->{kept}"))
        .collect();
    pairs.sort();
    format!("{label} {}", pairs.join(" "))
}

fn write_pair(dir: &std::path::Path, gene: i64, left: i64, right: i64) {
    fs::create_dir_all(dir).unwrap();
    let mut genes = String::from("cell,G\n");
    let mut tx = String::from("cell,G_A,G_B\n");
    for index in 0..22 {
        genes.push_str(&format!("c{index},{gene}\n"));
        tx.push_str(&format!("c{index},{left},{right}\n"));
    }
    fs::write(dir.join("gene_counts.csv"), genes).unwrap();
    fs::write(dir.join("transcript_counts.csv"), tx).unwrap();
    fs::write(
        dir.join("gene_transcript.tsv"),
        "genes\ttranscripts\nG\tG_A\nG\tG_B\n",
    )
    .unwrap();
}

fn dtu_lines(
    root: &std::path::Path,
    label: &str,
    a: (i64, i64, i64),
    b: (i64, i64, i64),
) -> String {
    let left = root.join(format!("{label}-a"));
    let right = root.join(format!("{label}-b"));
    write_pair(&left, a.0, a.1, a.2);
    write_pair(&right, b.0, b.1, b.2);
    let tables = differential_usage(
        &left.join("gene_counts.csv"),
        &left.join("transcript_counts.csv"),
        &left.join("gene_transcript.tsv"),
        &right.join("gene_counts.csv"),
        &right.join("transcript_counts.csv"),
        &right.join("gene_transcript.tsv"),
        &DtuSettings::default(),
    )
    .unwrap();
    let mut lines = Vec::new();
    for row in &tables.rows {
        let gene_p = row.p_gene.map(|v| format!("{v:.6e}")).unwrap_or_default();
        let log_fc = row.log_fc.map(|v| format!("{v:.6}")).unwrap_or_default();
        if let Some(tx) = &row.transcript {
            lines.push(format!(
                "dtu {label} isoform={} p_gene={gene_p} logfc={log_fc} alpha={:.6},{:.6} tu={:.6},{:.6} p_dtu={:.6e} p_tx={:.6e} switch={} es={:.6}",
                row.isoform,
                tx.alpha1,
                tx.alpha2,
                tx.tu1,
                tx.tu2,
                tx.p_dtu_gene,
                tx.p_transcript,
                tx.isoform_switch,
                tx.isoform_switch_es
            ));
        }
    }
    lines.join("\n")
}

#[test]
fn key_functions_match_the_synthetic_scotch_inputs() {
    let root = std::env::temp_dir().join(format!("plena-scotch-parity-{}", std::process::id()));
    let _ = fs::remove_dir_all(&root);
    fs::create_dir_all(&root).unwrap();
    let ann = parse_gtf(GENE).unwrap();
    let gene = &ann.genes[0];
    assert_eq!(
        gene.isoforms
            .iter()
            .map(|iso| iso.name.as_str())
            .collect::<Vec<_>>(),
        vec!["ALL", "LEFT", "RIGHT"]
    );

    let mut lines = Vec::new();
    let touch = Gene {
        name: "M".into(),
        id: "M".into(),
        chrom: "chr1".into(),
        start: 100,
        end: 400,
        strand: '+',
        exons: vec![(100, 200), (200, 250), (300, 400)],
        isoforms: vec![Isoform {
            name: "T".into(),
            exons: vec![0, 1, 2],
            order: 0,
        }],
        label: "M".into(),
    };
    let merged = gtf_exons(&touch, &[0, 1, 2]);
    lines.push(format!(
        "merge {}",
        merged
            .iter()
            .map(|(start, end)| format!("{}-{end}", start - 1))
            .collect::<Vec<_>>()
            .join(",")
    ));
    for isoform in &gene.isoforms {
        let vector: Vec<i8> = (0..gene.exons.len())
            .map(|exon| if isoform.exons.contains(&exon) { 1 } else { -1 })
            .collect();
        lines.push(format!("onehot {} {}", isoform.name, fmt_i8(&vector)));
    }
    let positions: Vec<i64> = (100..150).chain(300..310).collect();
    let hits = exon_base_counts(&positions, &gene.exons);
    lines.push(format!(
        "exon_hit {}",
        hits.iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",")
    ));

    let mut body = String::new();
    let specs = [
        ("FULL", 101, "100M100N100M100N100M"),
        ("LEFT", 101, "100M100N100M"),
        ("SKIP", 101, "100M300N100M"),
        ("PART40", 101, "40M"),
        ("PART5", 101, "5M"),
        ("ENDS", 101, "50M350N50M"),
        ("MID", 301, "100M"),
    ];
    for (name, pos, cigar) in specs {
        body.push_str(&sam_line(name, pos, cigar, &seq_for(query_len(cigar))));
    }
    let poly_seq = format!(
        "{}AAAAACCC{}{}AC",
        "ACGT".repeat(10),
        "T".repeat(10),
        "ACGT".repeat(10)
    );
    assert_eq!(poly_seq.len(), 100);
    body.push_str(&format!(
        "POLY\t0\tchr1\t101\t60\t100M\t*\t0\t0\t{poly_seq}\t*\tCB:Z:C1\tUB:Z:AAAAACCC\n"
    ));
    let sam = root.join("reads.sam");
    write_sam(&sam, &body);
    let loaded: Vec<_> = open_alignments(&sam)
        .unwrap()
        .map(|read| read.unwrap())
        .collect();
    for (name, _, _) in specs {
        let read = loaded.iter().find(|read| read.qname == name).unwrap();
        let counts = exon_base_counts(&read.ref_positions(), &gene.exons);
        lines.push(assign_line(name, gene, &counts, false));
    }
    let part = loaded.iter().find(|read| read.qname == "PART40").unwrap();
    let counts = exon_base_counts(&part.ref_positions(), &gene.exons);
    lines.push(assign_line("PART40_POLY", gene, &counts, true));
    let poly = loaded.iter().find(|read| read.qname == "POLY").unwrap();
    lines.push(format!(
        "poly_umi {}",
        if query_poly(&poly.seq, poly.query_start, poly.query_end, poly.tag("UB")) {
            "true"
        } else {
            "false"
        }
    ));

    // Gene B's only isoform uses the second exon, matching the Python dict.
    let pair_genes = vec![
        make_gene("GENEA", "GA", 100, 200, vec![(100, 200)], "ONE", vec![0]),
        make_gene(
            "GENEB",
            "GB",
            100,
            400,
            vec![(100, 180), (180, 400)],
            "ONLY",
            vec![1],
        ),
    ];
    let first = root.join("first.sam");
    write_sam(&first, &sam_line("FIRST", 101, "100M", &seq_for(100)));
    let first_read = open_alignments(&first).unwrap().next().unwrap().unwrap();
    let placement = crate::assign::place_in_metagene(
        &first_read.ref_positions(),
        first_read.start,
        first_read.end,
        &pair_genes,
        &[0, 1],
        &MatchParams::default(),
        &[false, false],
    )
    .unwrap();
    let chosen = &pair_genes[placement.gene];
    let chosen_counts = exon_base_counts(&first_read.ref_positions(), &chosen.exons);
    let found = diagnose(chosen, &chosen_counts, &MatchParams::default(), false);
    let detail = match &placement.outcome {
        Outcome::Novel(vector) => format!("novel {}", fmt_i8(vector)),
        Outcome::Known(_) | Outcome::Uncategorized => format!(
            "known {} scores={}",
            fmt_i8(&found.compatible),
            fmt_scores(&found.scores)
        ),
    };
    lines.push(format!("choose gene={} {detail}", placement.gene));

    lines.push(novel_line("novel10", &vec![vec![1, -1, 1]; 10]));
    lines.push(novel_line("novel9", &vec![vec![1, -1, 1]; 9]));
    let mut split = vec![vec![1, 1, -1]; 6];
    split.extend(vec![vec![-1, -1, 1]; 6]);
    lines.push(novel_line("novel_split", &split));
    lines.push(group_line("group_plus", '+'));
    lines.push(group_line("group_minus", '-'));

    let depth: Vec<i32> = (0..200).map(|pos| if pos < 100 { 80 } else { 5 }).collect();
    let cut = cut_by_derivative(&[(0, 200)], &depth, 0, 10.0);
    lines.push(format!(
        "cut {}",
        cut.iter()
            .map(|(a, b)| format!("{a}-{b}"))
            .collect::<Vec<_>>()
            .join(",")
    ));

    lines.push(dtu_lines(&root, "switch", (20, 18, 2), (20, 2, 18)));
    lines.push(dtu_lines(&root, "matched", (20, 18, 2), (20, 18, 2)));
    lines.push(dtu_lines(&root, "de", (20, 16, 4), (5, 4, 1)));

    // Outputs of merge_exons, isoform one-hot, exon_hit, map_read_to_gene,
    // find_novel_by_chunk, group_novel_isoform, detect_poly, and
    // cut_exons_by_derivative on these inputs.
    const UPSTREAM: &str = "\
merge 100-250,300-400
onehot ALL 1,1,1
onehot LEFT 1,1,-1
onehot RIGHT -1,1,1
exon_hit 50,10,0
assign FULL kind=known isoform=ALL scores=1.000000,-1.000000,-1.000000 vector=1,0,0
assign LEFT kind=known isoform=LEFT scores=-1.000000,1.000000,-1.000000 vector=0,1,0
assign SKIP kind=novel isoform=- scores=-1.000000,-1.000000,-1.000000 pct=1,0,1 vector=1,-1,1
assign PART40 kind=known isoform=LEFT scores=-1.000000,0.200000,-1.000000 vector=0,1,0
assign PART5 kind=known isoform=LEFT scores=-1.000000,0.025000,-1.000000 vector=0,1,0
assign ENDS kind=known isoform=ALL scores=0.333333,-1.000000,-1.000000 vector=1,0,0
assign MID kind=known isoform=LEFT scores=-1.000000,0.500000,-1.000000 vector=0,1,0
assign PART40_POLY kind=novel isoform=- scores=-1.000000,-1.000000,-1.000000 pct=0.400000,0,0 vector=1,-1,-1
poly_umi true
novel10 isoforms=novelIsoform_5 calls=r0:novelIsoform_5,r1:novelIsoform_5,r2:novelIsoform_5,r3:novelIsoform_5,r4:novelIsoform_5,r5:novelIsoform_5,r6:novelIsoform_5,r7:novelIsoform_5,r8:novelIsoform_5,r9:novelIsoform_5 empty=
novel9 isoforms= calls= empty=
novel_split isoforms=novelIsoform_3,novelIsoform_4 calls=r0:novelIsoform_3,r1:novelIsoform_3,r2:novelIsoform_3,r3:novelIsoform_3,r4:novelIsoform_3,r5:novelIsoform_3,r6:novelIsoform_4,r7:novelIsoform_4,r8:novelIsoform_4,r9:novelIsoform_4,r10:novelIsoform_4,r11:novelIsoform_4 empty=
group_plus novelIsoform_1->novelIsoform_1 novelIsoform_2->novelIsoform_3 novelIsoform_3->novelIsoform_3
group_minus novelIsoform_1->novelIsoform_3 novelIsoform_2->novelIsoform_2 novelIsoform_3->novelIsoform_3
cut 0-100,100-200
";
    let got: Vec<&str> = lines
        .iter()
        .filter(|line| !line.starts_with("dtu ") && !line.starts_with("choose "))
        .map(String::as_str)
        .collect();
    assert_eq!(got, UPSTREAM.lines().collect::<Vec<_>>());
    // Upstream keeps gene 1 here: after sorting it reads label 0, so the
    // known-versus-novel order does not affect the selected row.
    assert!(lines
        .iter()
        .any(|line| line == "choose gene=0 known 1 scores=1.000000"));
    assert!(lines
        .iter()
        .any(|line| line.contains("dtu switch") && line.contains("switch=true")));
    assert!(lines.iter().any(|line| line.contains("dtu matched")
        && line.contains("switch=false")
        && line.contains("p_dtu=1.000000e0")));
    assert!(lines.iter().any(|line| line.contains("dtu de")
        && line.contains("p_gene=5.993894e-11")
        && line.contains("switch=false")));
}
