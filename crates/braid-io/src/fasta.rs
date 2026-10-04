use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FastaRecord {
    /// Token before the first whitespace. This is Biopython `SeqRecord.id`.
    pub id: String,
    /// Header without the leading `>`.
    pub description: String,
    /// Sequence uppercased, matching `tama_collapse.py`.
    pub seq: String,
}

pub fn read_fasta(path: &Path) -> Result<Vec<FastaRecord>, String> {
    let file = File::open(path).map_err(|err| format!("open {}: {err}", path.display()))?;
    let reader = BufReader::new(file);
    let mut records = Vec::new();
    let mut header: Option<(String, String)> = None;
    let mut seq = String::new();

    for (lineno, line) in reader.lines().enumerate() {
        let line = line.map_err(|err| format!("read {}:{}: {err}", path.display(), lineno + 1))?;
        let line = line.trim_end_matches(['\r', '\n']);
        if let Some(rest) = line.strip_prefix('>') {
            if let Some((id, description)) = header.take() {
                records.push(FastaRecord {
                    id,
                    description,
                    seq: std::mem::take(&mut seq).to_ascii_uppercase(),
                });
            }
            let description = rest.to_string();
            let id = description
                .split_whitespace()
                .next()
                .unwrap_or("")
                .to_string();
            if id.is_empty() {
                return Err(format!(
                    "empty fasta header at {}:{}",
                    path.display(),
                    lineno + 1
                ));
            }
            header = Some((id, description));
        } else if header.is_some() {
            if !line.is_empty() {
                seq.push_str(line.trim());
            }
        } else if !line.is_empty() {
            return Err(format!(
                "fasta sequence before header at {}:{}",
                path.display(),
                lineno + 1
            ));
        }
    }
    if let Some((id, description)) = header {
        records.push(FastaRecord {
            id,
            description,
            seq: seq.to_ascii_uppercase(),
        });
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::Write;

    #[test]
    fn uppercases_and_uses_biopython_id() {
        let dir = std::env::temp_dir().join("braid-fasta-test");
        let _ = std::fs::create_dir_all(&dir);
        let path = dir.join("tiny.fa");
        let mut file = File::create(&path).unwrap();
        write!(file, ">chr1 dna_sm:scaffold extra\nacgt\nTT\n>chr2\nNa\n").unwrap();
        let records = read_fasta(&path).unwrap();
        assert_eq!(records[0].id, "chr1");
        assert_eq!(records[0].description, "chr1 dna_sm:scaffold extra");
        assert_eq!(records[0].seq, "ACGTTT");
        assert_eq!(records[1].seq, "NA");
    }
}
