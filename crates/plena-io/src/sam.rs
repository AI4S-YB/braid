use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Child, Command, Stdio};

/// Classes TAMA actually accepts. Any other flag is an error, matching the
/// `sam_flag_dict` lookup. This is not a SAM bitfield decoder.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SamClass {
    Forward,
    Reverse,
    Unmapped,
    Chimeric,
    NotPrimary,
}

impl SamClass {
    pub fn from_flag(flag: i32) -> Result<Self, String> {
        match flag {
            0 => Ok(SamClass::Forward),
            16 => Ok(SamClass::Reverse),
            4 => Ok(SamClass::Unmapped),
            2048 | 2064 => Ok(SamClass::Chimeric),
            256 | 272 => Ok(SamClass::NotPrimary),
            other => Err(format!("unrecognized SAM flag {other}")),
        }
    }

    pub fn keeps_model(self) -> bool {
        matches!(self, SamClass::Forward | SamClass::Reverse)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SamRecord {
    pub qname: String,
    pub flag: i32,
    pub class: SamClass,
    pub rname: String,
    pub pos: i64,
    pub cigar: String,
    pub seq: String,
    /// `XS:A:` from the optional fields. `None` when the tag is absent.
    pub xs: Option<char>,
}

pub fn parse_sam_line(line: &str) -> Result<SamRecord, String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return Err(format!("SAM line has {} columns", fields.len()));
    }
    let flag: i32 = fields[1]
        .parse()
        .map_err(|_| format!("bad SAM flag {}", fields[1]))?;
    let pos: i64 = fields[3]
        .parse()
        .map_err(|_| format!("bad SAM POS {}", fields[3]))?;
    let mut xs = None;
    for field in &fields[11..] {
        if let Some(rest) = field.strip_prefix("XS:A:") {
            xs = rest.chars().next();
        }
    }
    Ok(SamRecord {
        qname: fields[0].to_string(),
        flag,
        class: SamClass::from_flag(flag)?,
        rname: fields[2].to_string(),
        pos,
        cigar: fields[5].to_string(),
        seq: fields[9].to_string(),
        xs,
    })
}

pub fn read_sam(path: &Path) -> Result<Vec<SamRecord>, String> {
    alignment_reader(path, false)?.collect()
}

/// Streaming SAM records; BAM is decoded by samtools, as in upstream TAMA.
pub struct AlignmentReader {
    reader: Box<dyn BufRead>,
    child: Option<Child>,
    label: String,
    lineno: usize,
    finished: bool,
}

pub fn alignment_reader(path: &Path, bam: bool) -> Result<AlignmentReader, String> {
    let (reader, child): (Box<dyn BufRead>, Option<Child>) = if bam {
        let mut child = Command::new("samtools")
            .arg("view")
            .arg("--")
            .arg(path)
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .map_err(|e| format!("BAM input requires samtools on PATH: {e}"))?;
        let stdout = child
            .stdout
            .take()
            .ok_or("samtools stdout is unavailable")?;
        (Box::new(BufReader::new(stdout)), Some(child))
    } else {
        let file = File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
        (Box::new(BufReader::new(file)), None)
    };
    Ok(AlignmentReader {
        reader,
        child,
        label: path.display().to_string(),
        lineno: 0,
        finished: false,
    })
}

impl Iterator for AlignmentReader {
    type Item = Result<SamRecord, String>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        loop {
            let mut line = String::new();
            match self.reader.read_line(&mut line) {
                Ok(0) => {
                    self.finished = true;
                    if let Some(mut child) = self.child.take() {
                        match child.wait() {
                            Ok(status) if status.success() => {}
                            Ok(status) => {
                                return Some(Err(format!(
                                    "samtools view {} failed: {status}",
                                    self.label
                                )))
                            }
                            Err(e) => return Some(Err(format!("waiting for samtools: {e}"))),
                        }
                    }
                    return None;
                }
                Ok(_) => {
                    self.lineno += 1;
                    let line = line.trim_end_matches(['\r', '\n']);
                    if line.is_empty() || line.starts_with('@') {
                        continue;
                    }
                    return Some(
                        parse_sam_line(line)
                            .map_err(|e| format!("{}:{}: {e}", self.label, self.lineno)),
                    );
                }
                Err(e) => {
                    self.finished = true;
                    return Some(Err(format!("read {}: {e}", self.label)));
                }
            }
        }
    }
}

impl Drop for AlignmentReader {
    fn drop(&mut self) {
        if let Some(mut child) = self.child.take() {
            let _ = child.kill();
            let _ = child.wait();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_flags_and_unknown_flag() {
        assert_eq!(SamClass::from_flag(0).unwrap(), SamClass::Forward);
        assert_eq!(SamClass::from_flag(16).unwrap(), SamClass::Reverse);
        assert_eq!(SamClass::from_flag(2064).unwrap(), SamClass::Chimeric);
        assert!(SamClass::from_flag(1).is_err());
    }

    #[test]
    fn reads_xs_tag() {
        let line = "r1\t0\tchr\t10\t40\t5M\t*\t0\t0\tACGTN\t*\tXS:A:+\tNM:i:0";
        let rec = parse_sam_line(line).unwrap();
        assert_eq!(rec.pos, 10);
        assert_eq!(rec.cigar, "5M");
        assert_eq!(rec.xs, Some('+'));
        assert_eq!(rec.seq, "ACGTN");
    }
}
