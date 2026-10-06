# Isosceles strict quantification of one UMI-deduplicated BAM.
# The BAM name is fixed to "cell" so column names are cell.<barcode>.
# Arguments: bam, gtf, genome fasta, output TSV, barcode tag.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
    stop("usage: isosceles_quant.R bam gtf genome out_tsv barcode_tag")
}
bam <- args[[1]]
gtf <- args[[2]]
genome <- args[[3]]
out_tsv <- args[[4]]
barcode_tag <- args[[5]]

suppressPackageStartupMessages(library(Isosceles))
tx <- prepare_transcripts(gtf, genome, bam_parsed = NULL)
bams <- setNames(bam, "cell")
se_tcc <- bam_to_tcc(
    bams,
    tx,
    run_mode = "strict",
    min_read_count = 1,
    is_single_cell = TRUE,
    barcode_tag = barcode_tag,
    ncpu = 1
)
se <- tcc_to_transcript(se_tcc, ncpu = 1)
mat <- SummarizedExperiment::assay(se, "counts")
transcripts <- rownames(mat)
columns <- colnames(mat)
con <- file(out_tsv, "wt")
writeLines("barcode\ttranscript_id\tcount", con)
if (!is.null(columns) && length(columns) > 0 && nrow(mat) > 0) {
    for (j in seq_along(columns)) {
        barcode <- columns[[j]]
        prefix <- "cell."
        if (startsWith(barcode, prefix)) {
            barcode <- substring(barcode, nchar(prefix) + 1)
        }
        column <- as.numeric(mat[, j])
        hits <- which(column > 0)
        for (i in hits) {
            writeLines(paste(barcode, transcripts[[i]], column[[i]], sep = "\t"), con)
        }
    }
}
close(con)
