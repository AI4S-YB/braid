# Bambu-Clump quantification.
# The BAM is already one primary alignment per CB+UMI. bambu.singlecell
# (GoekeLab/bambu d704164) reads CB and UB and only deduplicates UMIs after
# a chromosome has more than 100 distinct UMIs. Arguments: bam, gtf, genome
# fasta, output TSV, barcode tag, UMI tag.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
    stop("usage: bambu_quant.R bam gtf genome out_tsv barcode_tag umi_tag")
}
bam <- args[[1]]
gtf <- args[[2]]
genome <- args[[3]]
out_tsv <- args[[4]]
barcode_tag <- args[[5]]
umi_tag <- args[[6]]
if (barcode_tag != "CB" || umi_tag != "UB") {
    stop("this bambu.singlecell build reads CB and UB tags")
}

suppressPackageStartupMessages(library(bambu))
if (!"bambu.singlecell" %in% getNamespaceExports("bambu")) {
    stop("this bambu build does not export bambu.singlecell")
}
annotations <- prepareAnnotations(gtf)
fn <- get("bambu.singlecell", envir = asNamespace("bambu"))
quant_data <- fn(
    reads = bam,
    output = "quantData",
    annotations = annotations,
    genome = genome,
    ncore = 1,
    verbose = FALSE
)
se <- fn(
    reads = quant_data,
    output = "EM",
    annotations = annotations,
    genome = genome,
    ncore = 1,
    verbose = FALSE,
    opt.em = list(degradationBias = FALSE)
)
if (!methods::is(se, "SummarizedExperiment")) {
    stop("bambu.singlecell did not return a SummarizedExperiment")
}
mat <- SummarizedExperiment::assay(se, "counts")
barcodes <- colnames(mat)
coldata <- as.data.frame(SummarizedExperiment::colData(se))
if ("barcode" %in% colnames(coldata)) {
    barcodes <- as.character(coldata$barcode)
} else if (barcode_tag %in% colnames(coldata)) {
    barcodes <- as.character(coldata[[barcode_tag]])
}
transcripts <- rownames(mat)
con <- file(out_tsv, "wt")
writeLines("barcode\ttranscript_id\tcount", con)
if (ncol(mat) > 0 && nrow(mat) > 0) {
    for (j in seq_len(ncol(mat))) {
        column <- as.numeric(mat[, j])
        hits <- which(column > 0)
        for (i in hits) {
            writeLines(paste(barcodes[[j]], transcripts[[i]], column[[i]], sep = "\t"), con)
        }
    }
}
close(con)
