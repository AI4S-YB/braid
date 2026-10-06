use plena_model::{MergeAlgo, MergeInput, MergeSource};
use plena_tama::TamaMerge;
use std::fs;
use std::path::Path;

#[test]
fn merge_trait_accepts_sources_and_returns_oracle_models() {
    let golden =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/parity/tama_modes/merge-gmap");
    let dir = std::env::temp_dir().join(format!("plena-merge-trait-{}", std::process::id()));
    fs::create_dir(&dir).unwrap();
    let output = TamaMerge
        .merge(&MergeInput {
            sources: vec![MergeSource {
                path: golden.join("s0.bed"),
                cap: "capped".into(),
                priority: "1,1,1".into(),
                name: "testmerge".into(),
            }],
            prefix: dir.join("out"),
        })
        .unwrap();
    let mut bed = String::new();
    for model in &output.transcripts {
        bed.push_str(
            &plena_io::tama_bed12(
                &model.chrom,
                model.strand,
                &model.exons,
                &format!("{};{}", model.gene_id, model.transcript_id),
            )
            .unwrap(),
        );
        bed.push('\n');
    }
    assert_eq!(bed, fs::read_to_string(golden.join("out.bed")).unwrap());
    assert_eq!(fs::read_dir(&dir).unwrap().count(), 4);
    for suffix in [
        ".bed",
        "_gene_report.txt",
        "_trans_report.txt",
        "_merge.txt",
    ] {
        assert_eq!(
            fs::read(dir.join(format!("out{suffix}"))).unwrap(),
            fs::read(golden.join(format!("out{suffix}"))).unwrap()
        );
    }
    fs::remove_dir_all(dir).unwrap();
}
