use plena_model::{Exon, Strand};

/// TAMA's BED12 line. Both start and end are written as `coord - 1`.
/// Exon lengths are `end - start`. Relative starts are `exon.start - transcript.start`.
pub fn tama_bed12(
    chrom: &str,
    strand: Strand,
    exons: &[Exon],
    name: &str,
) -> Result<String, String> {
    if exons.is_empty() {
        return Err("BED transcript has no exons".to_string());
    }
    let start = exons[0].start;
    let end = exons[exons.len() - 1].end;
    let mut lengths = Vec::with_capacity(exons.len());
    let mut rels = Vec::with_capacity(exons.len());
    for exon in exons {
        let length = exon.end - exon.start;
        let rel = exon.start - start;
        if length < 0 || rel < 0 {
            return Err(format!("negative exon coordinate in {name}"));
        }
        lengths.push(length.to_string());
        rels.push(rel.to_string());
    }
    let bed_start = start - 1;
    let bed_end = end - 1;
    Ok(format!(
        "{chrom}\t{bed_start}\t{bed_end}\t{name}\t40\t{}\t{bed_start}\t{bed_end}\t255,0,0\t{}\t{}\t{}",
        strand.as_bed(),
        exons.len(),
        lengths.join(","),
        rels.join(","),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn subtracts_one_from_both_ends() {
        let exons = vec![Exon { start: 10, end: 20 }, Exon { start: 30, end: 35 }];
        let line = tama_bed12("chr1", Strand::Forward, &exons, "G1;G1.1").unwrap();
        assert_eq!(
            line,
            "chr1\t9\t34\tG1;G1.1\t40\t+\t9\t34\t255,0,0\t2\t10,5\t0,20"
        );
    }
}
