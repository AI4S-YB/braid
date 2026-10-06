use std::fs;
use std::path::Path;

use plena_model::{CollapseAlgo, CollapseInput, Strand};
use plena_tama::TamaCollapse;

#[test]
fn collapse_trait_returns_oracle_models_and_writes_outputs() {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/parity/gmap_collapse");
    let dir = std::env::temp_dir().join(format!("plena-collapse-trait-{}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    let input = CollapseInput {
        alignments: golden.join("gmap_test.sam"),
        genome: golden.join("test_genome.fa"),
        prefix: dir.join("gmap"),
    };
    let output = TamaCollapse.collapse(&input).unwrap();
    assert_eq!(output.transcripts.len(), 32);
    assert!(output
        .transcripts
        .iter()
        .any(|t| t.strand == Strand::Reverse));
    let mut bed = String::new();
    for transcript in output.transcripts {
        bed.push_str(
            &plena_io::tama_bed12(
                &transcript.chrom,
                transcript.strand,
                &transcript.exons,
                &format!("{};{}", transcript.gene_id, transcript.transcript_id),
            )
            .unwrap(),
        );
        bed.push('\n');
        assert!(transcript.cds_start.is_none());
        assert!(transcript.cds_end.is_none());
    }
    assert_eq!(bed, fs::read_to_string(golden.join("gmap.bed")).unwrap());
    assert_eq!(fs::read_dir(&dir).unwrap().count(), 10);
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
            .replace(
                "test_files/gmap_test.sam",
                &input.alignments.display().to_string(),
            )
            .replace(
                "test_files/test_genome.fa",
                &input.genome.display().to_string(),
            )
            .replace(
                "/work/plena/tests/parity/gmap_collapse/gmap",
                &input.prefix.display().to_string(),
            )
            .replace("-log log_off", "-log log_on");
        assert_eq!(
            fs::read_to_string(dir.join(format!("gmap{suffix}"))).unwrap(),
            expected,
            "{suffix}"
        );
    }
    fs::remove_dir_all(dir).unwrap();
}
