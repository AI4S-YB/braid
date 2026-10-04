//! SAM records, and BAM records decoded by `samtools view`.
//!
//! Reference coordinates are 0-based and half-open. `N` opens a splice
//! junction and splits alignment blocks. `M`, `=`, and `X` contribute
//! reference positions used for exon overlap.

use std::collections::BTreeMap;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::process::{Child, Command, Stdio};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OpKind {
    Match,
    Insert,
    Delete,
    Skip,
    Soft,
    Other,
}

#[derive(Clone, Copy, Debug)]
pub struct Op {
    pub kind: OpKind,
    pub len: i64,
}

#[derive(Clone, Debug)]
pub struct Read {
    pub qname: String,
    pub flag: i32,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub ops: Vec<Op>,
    pub seq: String,
    pub tags: BTreeMap<String, String>,
    pub query_start: i64,
    pub query_end: i64,
}

impl Read {
    pub fn tag(&self, name: &str) -> Option<&str> {
        self.tags.get(name).map(String::as_str)
    }

    /// Query bases in the alignment, excluding soft clips. This is pysam
    /// `query_alignment_length` (`M`/`=`/`X`/`I`).
    pub fn query_alignment_length(&self) -> i64 {
        self.ops
            .iter()
            .filter(|op| matches!(op.kind, OpKind::Match | OpKind::Insert))
            .map(|op| op.len)
            .sum()
    }

    pub fn ref_positions(&self) -> Vec<i64> {
        let mut pos = self.start;
        let mut out = Vec::new();
        for op in &self.ops {
            match op.kind {
                OpKind::Match => {
                    for _ in 0..op.len {
                        out.push(pos);
                        pos += 1;
                    }
                }
                OpKind::Delete | OpKind::Skip => pos += op.len,
                OpKind::Insert | OpKind::Soft | OpKind::Other => {}
            }
        }
        out
    }

    pub fn blocks(&self) -> Vec<(i64, i64)> {
        let mut pos = self.start;
        let mut block: Option<(i64, i64)> = None;
        let mut out = Vec::new();
        for op in &self.ops {
            match op.kind {
                OpKind::Match | OpKind::Delete => {
                    let start = block.map_or(pos, |(start, _)| start);
                    pos += op.len;
                    block = Some((start, pos));
                }
                OpKind::Skip => {
                    if let Some(done) = block.take() {
                        out.push(done);
                    }
                    pos += op.len;
                }
                OpKind::Insert | OpKind::Soft | OpKind::Other => {}
            }
        }
        if let Some(done) = block {
            out.push(done);
        }
        out
    }

    pub fn junctions(&self) -> Vec<(i64, i64)> {
        let mut pos = self.start;
        let mut out = Vec::new();
        for op in &self.ops {
            match op.kind {
                OpKind::Skip => {
                    out.push((pos, pos + op.len));
                    pos += op.len;
                }
                OpKind::Match | OpKind::Delete => pos += op.len,
                OpKind::Insert | OpKind::Soft | OpKind::Other => {}
            }
        }
        out
    }
}

pub struct ReadIter {
    reader: Box<dyn BufRead>,
    child: Option<Child>,
    label: String,
    line_no: usize,
    finished: bool,
}

pub fn open_alignments(path: &Path) -> Result<ReadIter, String> {
    let bam = path
        .extension()
        .and_then(|ext| ext.to_str())
        .is_some_and(|ext| ext.eq_ignore_ascii_case("bam"));
    let (reader, child): (Box<dyn BufRead>, Option<Child>) = if bam {
        let mut child = Command::new("samtools")
            .args(["view", "--"])
            .arg(path)
            .stdout(Stdio::piped())
            .stderr(Stdio::inherit())
            .spawn()
            .map_err(|err| format!("BAM input requires samtools on PATH: {err}"))?;
        let stdout = child
            .stdout
            .take()
            .ok_or("samtools stdout is unavailable")?;
        (Box::new(BufReader::new(stdout)), Some(child))
    } else {
        let file =
            std::fs::File::open(path).map_err(|err| format!("open {}: {err}", path.display()))?;
        (Box::new(BufReader::new(file)), None)
    };
    Ok(ReadIter {
        reader,
        child,
        label: path.display().to_string(),
        line_no: 0,
        finished: false,
    })
}

impl Iterator for ReadIter {
    type Item = Result<Read, String>;

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
                                )));
                            }
                            Err(err) => {
                                return Some(Err(format!("waiting for samtools: {err}")));
                            }
                        }
                    }
                    return None;
                }
                Ok(_) => {
                    self.line_no += 1;
                    let line = line.trim_end_matches(['\r', '\n']);
                    if line.is_empty() || line.starts_with('@') {
                        continue;
                    }
                    return Some(
                        parse_read(line)
                            .map_err(|err| format!("{}:{}: {err}", self.label, self.line_no)),
                    );
                }
                Err(err) => {
                    self.finished = true;
                    return Some(Err(format!("read {}: {err}", self.label)));
                }
            }
        }
    }
}

fn parse_read(line: &str) -> Result<Read, String> {
    let fields: Vec<&str> = line.split('\t').collect();
    if fields.len() < 11 {
        return Err(format!("SAM line has {} columns", fields.len()));
    }
    let flag: i32 = fields[1]
        .parse()
        .map_err(|_| format!("bad flag {}", fields[1]))?;
    let pos: i64 = fields[3]
        .parse()
        .map_err(|_| format!("bad POS {}", fields[3]))?;
    let ops = parse_cigar(fields[5])?;
    let seq = if fields[9] == "*" {
        String::new()
    } else {
        fields[9].to_string()
    };
    let mut tags = BTreeMap::new();
    for field in &fields[11..] {
        let mut parts = field.splitn(3, ':');
        let Some(tag) = parts.next() else { continue };
        let Some(kind) = parts.next() else { continue };
        let Some(value) = parts.next() else { continue };
        if tag.len() != 2 {
            continue;
        }
        let text = match kind {
            "A" | "Z" | "H" | "B" => value.to_string(),
            "i" | "f" => value.to_string(),
            _ => value.to_string(),
        };
        tags.insert(tag.to_string(), text);
    }
    let start = pos - 1;
    let mut end = start;
    let mut query_start = 0i64;
    let mut seen_align = false;
    let mut query_pos = 0i64;
    for op in &ops {
        match op.kind {
            OpKind::Soft => {
                if !seen_align {
                    query_start += op.len;
                }
                query_pos += op.len;
            }
            OpKind::Match => {
                seen_align = true;
                end += op.len;
                query_pos += op.len;
            }
            OpKind::Insert => query_pos += op.len,
            OpKind::Delete | OpKind::Skip => end += op.len,
            OpKind::Other => {}
        }
    }
    if flag & 0x4 != 0 || fields[2] == "*" || ops.is_empty() {
        end = start;
    }
    Ok(Read {
        qname: fields[0].to_string(),
        flag,
        chrom: fields[2].to_string(),
        start,
        end,
        ops,
        seq,
        tags,
        query_start,
        query_end: query_pos,
    })
}

fn parse_cigar(text: &str) -> Result<Vec<Op>, String> {
    if text == "*" {
        return Ok(Vec::new());
    }
    let mut ops = Vec::new();
    let mut number = String::new();
    for ch in text.chars() {
        if ch.is_ascii_digit() {
            number.push(ch);
            continue;
        }
        let len: i64 = number.parse().map_err(|_| format!("bad CIGAR {text}"))?;
        number.clear();
        let kind = match ch {
            'M' | '=' | 'X' => OpKind::Match,
            'I' => OpKind::Insert,
            'D' => OpKind::Delete,
            'N' => OpKind::Skip,
            'S' => OpKind::Soft,
            'H' | 'P' => OpKind::Other,
            other => return Err(format!("bad CIGAR operator {other}")),
        };
        ops.push(Op { kind, len });
    }
    if !number.is_empty() {
        return Err(format!("bad CIGAR {text}"));
    }
    Ok(ops)
}

pub fn is_unmapped(read: &Read) -> bool {
    read.flag & 0x4 != 0 || read.chrom == "*" || read.end <= read.start
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spliced_cigar_emits_blocks_junctions_and_positions() {
        let line = format!(
            "R\t0\tchr1\t101\t60\t100M100N100M\t*\t0\t0\t{}\t*\tCB:Z:C\tUB:Z:U",
            "A".repeat(200)
        );
        let read = parse_read(&line).unwrap();
        assert_eq!(read.start, 100);
        assert_eq!(read.end, 400);
        assert_eq!(read.blocks(), vec![(100, 200), (300, 400)]);
        assert_eq!(read.junctions(), vec![(200, 300)]);
        assert_eq!(read.ref_positions().len(), 200);
        assert_eq!(read.tag("CB"), Some("C"));
        assert!(!is_unmapped(&read));
    }
}
