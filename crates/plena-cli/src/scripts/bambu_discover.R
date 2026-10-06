# Bambu transcriptome discovery from a genome BAM.
# Arguments: bam, gtf, genome fasta, output GTF.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("usage: bambu_discover.R bam gtf genome out_gtf")
}
bam <- args[[1]]
gtf <- args[[2]]
genome <- args[[3]]
out_gtf <- args[[4]]

suppressPackageStartupMessages(library(bambu))
annotations <- prepareAnnotations(gtf)
call <- list(
    reads = bam,
    annotations = annotations,
    genome = genome,
    ncore = 1L,
    verbose = FALSE
)
fmls <- names(formals(bambu))
if ("discovery" %in% fmls) call$discovery <- TRUE
if ("quant" %in% fmls) call$quant <- FALSE
if ("assignDist" %in% fmls) call$assignDist <- FALSE
if ("NDR" %in% fmls) call$NDR <- 1
if ("ndr" %in% fmls) call$ndr <- 1
se <- do.call(bambu, call)
if ("writeToGTF" %in% getNamespaceExports("bambu")) {
    writeToGTF(se, out_gtf)
} else {
    ranges <- if (methods::is(se, "SummarizedExperiment")) SummarizedExperiment::rowRanges(se) else se
    rtracklayer::export(ranges, out_gtf, format = "gtf")
}
