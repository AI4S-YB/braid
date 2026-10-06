//! FLAIR combine and precise isoform collapse.
//!
//! Independent Rust implementation of `flair_combine.py` and
//! `collapse_isoforms_precise.py` from BrooksLabUCSC/flair commit
//! `573414c551332bf6348a9d04ba7cf562f67416cb`. Upstream is BSD-3-Clause,
//! copyright 2019-2026 The Regents of the University of California. Output
//! files match those scripts, including the quirks called out in each module.

mod collapse;
mod combine;

pub use collapse::{
    alignments_to_bed12, bed12_to_gtf, collapse_precise_bed, collapse_precise_files,
    PreciseSettings,
};
pub use combine::{combine, write_combine_texts, CombineSettings, CombineTexts};
