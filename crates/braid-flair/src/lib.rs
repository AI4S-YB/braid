//! FLAIR combine.
//!
//! Independent Rust implementation of `flair_combine.py` from
//! BrooksLabUCSC/flair commit `573414c551332bf6348a9d04ba7cf562f67416cb`.
//! Upstream is BSD-3-Clause, copyright 2019-2026 The Regents of the
//! University of California. Output files match that script, including the
//! quirks called out in `combine`.

mod combine;

pub use combine::{combine, write_combine_texts, CombineSettings, CombineTexts};
