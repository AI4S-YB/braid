//! Internal transcript model. Coordinates are TACO's: 0-based, end exclusive.

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub(crate) enum Strand {
    Pos = 0,
    Neg = 1,
    Na = 2,
}

impl Strand {
    pub(crate) fn from_gtf(text: &str) -> Option<Self> {
        match text {
            "+" => Some(Strand::Pos),
            "-" => Some(Strand::Neg),
            "." => Some(Strand::Na),
            _ => None,
        }
    }

    pub(crate) fn as_gtf(self) -> &'static str {
        match self {
            Strand::Pos => "+",
            Strand::Neg => "-",
            Strand::Na => ".",
        }
    }

    pub(crate) const ALL: [Strand; 3] = [Strand::Pos, Strand::Neg, Strand::Na];
}

#[derive(Clone, Debug)]
pub(crate) struct Transfrag {
    pub chrom: String,
    pub strand: Strand,
    pub exons: Vec<(i64, i64)>,
    pub id: String,
    pub expr: f64,
    pub is_ref: bool,
}

impl Transfrag {
    pub(crate) fn length(&self) -> i64 {
        self.exons.iter().map(|(start, end)| end - start).sum()
    }

    pub(crate) fn start(&self) -> i64 {
        self.exons[0].0
    }

    pub(crate) fn end(&self) -> i64 {
        self.exons.last().map(|exon| exon.1).unwrap_or(0)
    }

    /// Transcription start. On the negative strand this is the exon end.
    pub(crate) fn tx_start(&self) -> i64 {
        if self.strand == Strand::Neg {
            self.end()
        } else {
            self.start()
        }
    }

    pub(crate) fn tx_stop(&self) -> i64 {
        if self.strand == Strand::Neg {
            self.start()
        } else {
            self.end()
        }
    }

    pub(crate) fn introns(&self) -> Vec<(i64, i64)> {
        let mut out = Vec::new();
        for pair in self.exons.windows(2) {
            out.push((pair[0].1, pair[1].0));
        }
        out
    }

    pub(crate) fn splices(&self) -> Vec<i64> {
        let mut out = Vec::new();
        for (left, right) in self.introns() {
            out.push(left);
            out.push(right);
        }
        out
    }
}
