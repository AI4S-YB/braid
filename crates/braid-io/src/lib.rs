//! Readers and writers. SAM parsing keeps the raw CIGAR string. TAMA's
//! `X`/`=` folding lives in `braid-tama`, because a normal SAM parser must
//! not treat mismatch operators as matches.

mod bed;
mod fasta;
mod sam;

pub use bed::tama_bed12;
pub use fasta::{read_fasta, FastaRecord};
pub use sam::{alignment_reader, read_sam, AlignmentReader, SamClass, SamRecord};

/// Append an output suffix without converting the native path to UTF-8.
pub fn suffix_path(prefix: &std::path::Path, suffix: &str) -> std::path::PathBuf {
    let mut name = prefix.as_os_str().to_os_string();
    name.push(suffix);
    name.into()
}
