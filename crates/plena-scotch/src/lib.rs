//! SCOTCH isoform quantification for single-cell full-length long reads.
//!
//! Independent Rust implementation of the preprocessing and DTU steps in
//! WGLab/SCOTCH commit `15d6ad8b6319cf806c36f677ffe3ebcc003cae16`
//! (MIT, copyright 2024 Wang Genomics Lab). The pipeline reads a reference
//! GTF and one BAM or SAM group per sample, assigns each cell UMI to a
//! known or novel isoform, and writes gene and transcript counts.
//!
//! Louvain community ids and tied gene choices are deterministic. Upstream
//! breaks those ties at random, and it draws a random isoform when a read
//! stays compatible with more than one model after scoring. A novel isoform
//! is kept only when a discovery chunk assigns at least 10 reads, matching
//! `find_novel_by_chunk`.

mod assign;
mod bam;
mod coverage;
mod dtu;
mod gtf;
mod novel;
mod run;

#[cfg(test)]
mod parity;

pub use assign::MatchParams;
pub use dtu::{differential_usage, write_dtu, DtuFailure, DtuSettings, DtuTables};
pub use run::{quantify, Platform, ScotchFailure, ScotchSettings};
