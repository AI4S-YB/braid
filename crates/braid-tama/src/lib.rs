//! TAMA collapse and merge.
//!
//! Dict iteration that reaches a file or a grouping decision has to follow
//! CPython 2.7 with `PYTHONHASHSEED=0`. That order is not insertion order.
//! See `tests/parity/dict_order.txt`. `Py27Dict` lives in this crate and
//! must not be reused by TACO or other algorithms.
//!
//! Coverage, identity, and polyA percentages use Python 2 `str(round(x, 2))`.
//! The frozen table is `tests/parity/py2_round_str.tsv`.

mod cigar;
mod error;
mod group;
mod ids;
mod lde;
mod merge;
mod polya;
mod py27;
mod py2float;
mod run;
pub use merge::{merge_sources, read_merge_sources, write_merge_texts, MergeSettings, MergeTexts};

pub use cigar::{cigar_parts, mapped_seq_length, trans_coordinates, CigarParts};
pub use error::{calc_error_rate, coverage_percent, identity_percent, ErrorRate, VariationBook};
pub use lde::{local_density, LdeRead, LdeResult};
pub use polya::{detect_polya, reverse_complement, PolyA};
pub use py27::{py27_int_hash, py27_str_hash, Py27Dict, Py27Key};
pub use py2float::py2_str_round;
pub use run::{
    collapse_to_files, original_collapse, write_collapse_texts, CollapseSettings, CollapseTexts,
};

use braid_model::{
    CollapseAlgo, CollapseInput, CollapseOutput, MergeAlgo, MergeInput, MergeOutput,
};

pub const TAMA_COLLAPSE_DATE: &str = "tc_version_date_2023_03_28";

/// Default TAMA collapse (capped, original). Detects BAM by extension, requiring
/// samtools. Writes ten output files and returns the models for composition.
pub struct TamaCollapse;

impl CollapseAlgo for TamaCollapse {
    fn name(&self) -> &'static str {
        "tama"
    }

    fn collapse(&self, input: &CollapseInput) -> Result<CollapseOutput, String> {
        let settings = CollapseSettings {
            bam: if input
                .alignments
                .extension()
                .is_some_and(|ext| ext.eq_ignore_ascii_case("bam"))
            {
                "BAM".into()
            } else {
                "SAM".into()
            },
            sam_label: input.alignments.display().to_string(),
            fasta_label: input.genome.display().to_string(),
            prefix_label: input.prefix.display().to_string(),
            ..CollapseSettings::default()
        };
        let texts = original_collapse(&input.alignments, &input.genome, &settings)?;
        write_collapse_texts(&input.prefix, &texts)?;
        Ok(CollapseOutput {
            transcripts: texts.transcripts,
        })
    }
}

/// Merge BED sources with default TAMA thresholds and priorities supplied by
/// each MergeSource. Writes BED, transcript/gene reports and the source map.
pub struct TamaMerge;

impl MergeAlgo for TamaMerge {
    fn name(&self) -> &'static str {
        "tama"
    }

    fn merge(&self, input: &MergeInput) -> Result<MergeOutput, String> {
        let texts = merge_sources(&input.sources, &MergeSettings::default())?;
        write_merge_texts(&input.prefix, &texts)?;
        Ok(MergeOutput {
            transcripts: texts.transcripts,
        })
    }
}
