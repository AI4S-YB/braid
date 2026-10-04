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
mod lde;
mod polya;
mod py27;
mod py2float;
mod run;

pub use cigar::{cigar_parts, mapped_seq_length, trans_coordinates, CigarParts};
pub use error::{calc_error_rate, coverage_percent, identity_percent, ErrorRate, VariationBook};
pub use lde::{local_density, LdeRead, LdeResult};
pub use polya::{detect_polya, reverse_complement, PolyA};
pub use py27::{py27_int_hash, py27_str_hash, Py27Dict, Py27Key};
pub use py2float::py2_str_round;
pub use run::{original_collapse, write_collapse_texts, CollapseSettings, CollapseTexts};

use braid_model::{
    CollapseAlgo, CollapseInput, CollapseOutput, MergeAlgo, MergeInput, MergeOutput,
};

pub const TAMA_COLLAPSE_DATE: &str = "tc_version_date_2023_03_28";

pub struct TamaCollapse;

impl CollapseAlgo for TamaCollapse {
    fn name(&self) -> &'static str {
        "tama"
    }

    fn collapse(&self, _input: &CollapseInput) -> Result<CollapseOutput, String> {
        Err("TAMA collapse is not ported yet".to_string())
    }
}

pub struct TamaMerge;

impl MergeAlgo for TamaMerge {
    fn name(&self) -> &'static str {
        "tama"
    }

    fn merge(&self, _input: &MergeInput) -> Result<MergeOutput, String> {
        Err("TAMA merge is not ported yet".to_string())
    }
}
