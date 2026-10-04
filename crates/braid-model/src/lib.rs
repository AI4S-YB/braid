//! Transcript annotations shared by every algorithm.
//!
//! TAMA-only fields (wobble, splice-error strings, polyA, variants) stay in
//! `braid-tama`. Coordinates here are TAMA's internal coordinates: 1-based
//! start, and end exclusive in that same 1-based system (`POS + reference
//! bases consumed`). The BED writer subtracts 1 from both, which is what
//! TAMA writes.

use std::path::PathBuf;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub fn as_bed(self) -> &'static str {
        match self {
            Strand::Forward => "+",
            Strand::Reverse => "-",
        }
    }
}

/// One exon in internal coordinates. `end` is exclusive.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Exon {
    pub start: i64,
    pub end: i64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Transcript {
    pub chrom: String,
    pub strand: Strand,
    pub exons: Vec<Exon>,
    pub gene_id: String,
    pub transcript_id: String,
    pub source: Option<String>,
    /// Expression, unused by TAMA. TACO will use it.
    pub score: Option<f64>,
    pub cds_start: Option<i64>,
    pub cds_end: Option<i64>,
}

#[derive(Clone, Debug)]
pub struct CollapseInput {
    pub alignments: PathBuf,
    pub genome: PathBuf,
    pub prefix: PathBuf,
}

#[derive(Clone, Debug)]
pub struct CollapseOutput {
    pub transcripts: Vec<Transcript>,
}

#[derive(Clone, Debug)]
pub struct MergeSource {
    pub path: PathBuf,
    pub cap: String,
    pub priority: String,
    pub name: String,
}

#[derive(Clone, Debug)]
pub struct MergeInput {
    pub sources: Vec<MergeSource>,
    pub prefix: PathBuf,
}

#[derive(Clone, Debug)]
pub struct MergeOutput {
    pub transcripts: Vec<Transcript>,
}

pub trait CollapseAlgo {
    fn name(&self) -> &'static str;
    fn collapse(&self, input: &CollapseInput) -> Result<CollapseOutput, String>;
}

pub trait MergeAlgo {
    fn name(&self) -> &'static str;
    fn merge(&self, input: &MergeInput) -> Result<MergeOutput, String>;
}
