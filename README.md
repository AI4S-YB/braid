# braid

Reconcile transcript models with Rust implementations of **TAMA collapse**, **TAMA merge**, **FLAIR combine**, **FLAIR precise collapse**, **TACO meta-assembly**, and **SCOTCH** isoform quantification. Pooled construction and per-cell quantification also call pinned builds of **IsoQuant**, **Bambu**, **Bambu-Clump**, **Isosceles**, and **StringTie**.

The project focuses on reproducing upstream algorithm behavior and output formats. TAMA, FLAIR, TACO, and SCOTCH run natively in Rust. IsoQuant, Bambu, Isosceles, and StringTie stay external programs. Reading BAM files requires `samtools` on `PATH`.

## Commands

| Command | Input | Main output |
| --- | --- | --- |
| `braid tama-collapse` | Sorted SAM/BAM and genome FASTA | Collapsed BED12 models and read reports |
| `braid tama-merge` | Manifest of annotated BED12 files | Merged BED12 models and source reports |
| `braid flair-combine` | Manifest of FLAIR transcriptomes | Combined BED12, counts, and isoform map |
| `braid taco` | Manifest of sample GTF files | Assembled GTF/BED and diagnostic tracks |
| `braid scotch` | Reference GTF and one BAM/SAM group per sample | Per-sample gene and transcript counts |
| `braid scotch-dtu` | Two SCOTCH count directories | Gene Wilcoxon and transcript usage tests |
| `braid discover` | Coordinate-sorted genome BAM, genome FASTA, optional guide GTF | `transcripts.gtf` for one pooled transcriptome |
| `braid quant` | Barcode-tagged genome BAM, genome FASTA, and annotation GTF | `counts.tsv` keyed by cell barcode |
| `braid stringtie-merge` | Two or more transcript GTFs | One merged GTF |

The equivalent general commands are `collapse --algo tama`, `merge --algo tama`, and `combine --algo flair`. Those three commands still accept only `tama` or `flair`. StringTie merge is the separate `stringtie-merge` command.

## Install

### Build from source

CI uses Rust **1.96.1**. Install a Rust toolchain with Cargo, then:

```sh
git clone https://github.com/AI4S-YB/braid.git
cd braid
cargo build --release --locked
./target/release/braid --help
```

To install the executable into Cargo's binary directory:

```sh
cargo install --path crates/braid-cli --locked
braid --version
```

The examples below assume `braid` is on `PATH`. Alternatively, use `./target/release/braid` from the repository root.

### Release binaries

See [GitHub Releases](https://github.com/AI4S-YB/braid/releases). The release workflow builds Linux x86_64 GNU and Windows x86_64 MSVC archives, with SHA-256 checksums. Linux binaries require a compatible glibc environment.

## Quick start

Run TAMA collapse on the bundled test data from the repository root:

```sh
mkdir -p demo
braid tama-collapse \
  -s tests/parity/gmap_collapse/gmap_test.sam \
  -f tests/parity/gmap_collapse/test_genome.fa \
  -p demo/gmap
```

The main annotation is `demo/gmap.bed`. TAMA output prefixes require an existing parent directory.

Get the complete option list for any command:

```sh
braid tama-collapse --help
braid tama-merge --help
braid flair-combine --help
braid taco --help
```

## TAMA collapse

Provide alignments sorted by reference and position, together with the matching genome FASTA:

```sh
mkdir -p results
braid tama-collapse \
  -s reads.sorted.bam -b BAM \
  -f genome.fa \
  -p results/sample \
  -x no_cap --rm low_mem
```

The CLI defaults to SAM input: **use `-b BAM` for BAM files**, regardless of the filename extension. BAM decoding uses `samtools view`.

| Option | Default | Meaning |
| --- | --- | --- |
| `-x` | `capped` | `capped` or `no_cap` grouping |
| `--rm` | `original` | `original` or `low_mem` execution |
| `-e` | `common_ends` | `common_ends` or `longest_ends` |
| `-c` / `-i` | `99` / `85` | Minimum coverage / identity percentages |
| `-a` / `-m` / `-z` | `10` / `10` / `10` | 5′, splice-junction, and 3′ coordinate thresholds |
| `-d` | `merge_dup` | Merge duplicate models, or use `no_merge` |

Legacy multi-character options such as `-rm`, `-icm`, `-sj`, and `-vc` are also accepted in place of their double-dash spellings. `braid tama-collapse -v 1` prints the upstream TAMA version date; `braid --version` prints the braid version.

For prefix `sample`, output files are:

```text
sample.bed
sample_read.txt
sample_trans_report.txt
sample_trans_read.bed
sample_polya.txt
sample_strand_check.txt
sample_local_density_error.txt
sample_report.txt
sample_variants.txt                 # original mode only
sample_varcov.txt                   # original mode only
```

`low_mem` retains the genome and current locus while streaming output. The two modes preserve their upstream differences, including multimapping and duplicate-read handling; they are not interchangeable memory settings. In `original` mode, variant support is fixed to 5 by upstream behavior even when `--vc` is supplied.

## TAMA merge

Create a headerless, tab-separated manifest with four columns:

```text
BED path    capped/no_cap    priorities    source name
```

For example, this command writes real tab delimiters:

```sh
printf 'sample1.bed\tcapped\t1,1,1\tS1\nsample2.bed\tno_cap\t2,2,2\tS2\n' > merge.tsv
mkdir -p results
braid tama-merge -f merge.tsv -p results/merged -d merge_dup
```

The manifest must contain no spaces or blank rows. Each BED12 name must contain `gene_id;transcript_id`, and strand must be `+` or `-`. Priorities are three comma-separated integers used in coordinate voting.

Merge defaults differ from collapse: `-a 20`, `-m 10`, `-z 20`, and `-d no_merge`. Use `-d merge_dup` to permit duplicate-model merging. `-s S1` retains IDs from a named source; `--cds S1` (also `-cds S1`) retains its CDS coordinates. Comma-separated source names are accepted for both options.

Outputs are `<prefix>.bed`, `<prefix>_trans_report.txt`, `<prefix>_gene_report.txt`, and `<prefix>_merge.txt`.

## FLAIR combine

The headerless manifest has three to five tab-separated columns:

```text
sample    type    BED path    [FASTA path]    [read-map path]
```

Use `isoforms` for ordinary transcriptomes and `fusionisoform` for fusion models. For a BED-only example:

```sh
printf 'S1\tisoforms\tsample1.bed\nS2\tisoforms\tsample2.bed\n' > flair.tsv
braid flair-combine -m flair.tsv -o results/combined
```

Use FLAIR-style transcript/gene names in the BED input. Optional FASTA and read-map files should correspond to those isoforms. To provide a read map without FASTA, leave the fourth column empty. Manifest rows must not be blank.

- `-w 200` sets the transcript-end comparison window.
- `-p 10` requires usage greater than 10% for the usage filter.
- `-f` selects `usageandlongest` (default), `usageonly`, `none`, or an integer read-count threshold.
- `-s` includes single-exon isoforms, which are excluded by default.
- `-c` additionally converts the combined BED to GTF.

Outputs are `<prefix>.bed`, `<prefix>.counts.tsv`, and `<prefix>.isoform.map.txt`. `<prefix>.fa` is written only when **every** manifest row supplies a FASTA path; otherwise all FASTA paths are ignored. `-c` adds `<prefix>.gtf`.

The FASTA reader intentionally follows the pinned FLAIR script: sequence lines replace rather than concatenate earlier lines for a record. Supply one sequence line per record.

## TACO meta-assembly

Create a headerless manifest containing a GTF path and an optional sample ID, separated by a tab:

```sh
printf 'sample1.gtf\tS1\nsample2.gtf\tS2\n' > samples.tsv
braid taco samples.tsv -o results/taco
```

Each sample GTF must contain a `transcript` feature before its `exon` features, with matching `transcript_id` attributes. Sample transcript features need an expression attribute, `FPKM` by default; select another with `--gtf-expr-attr`. Expression is normalized within each sample before filtering. Sample paths and IDs must be unique.

The output directory **must not already exist**. Main results are `assembly.gtf` and `assembly.bed`; additional files contain sample statistics, filtered transcripts, loci, expression tracks, splice graphs, change points, and path statistics.

Useful options include `--filter-min-length`, `--filter-min-expr`, `--isoform-frac`, and `--max-isoforms`. Guided modes require `--ref-gtf`; splice-junction filtering requires both `--filter-splice-juncs` and `--ref-genome-fasta`.

TACO currently runs in one process. `--num-processes` is accepted for command-line compatibility but does not enable parallel execution. `--resume` and `--assemble` are unsupported. Change-point detection and trimming are enabled by default; unresolved unstranded transcripts require `--assemble-unstranded` to be assembled.

For all manifests above, relative input paths are resolved against the **current working directory**, not the manifest's directory.

## SCOTCH quantification

`braid scotch` assigns each cell UMI in a full-length long-read BAM or SAM to a known or novel isoform of a gene in the reference GTF.

```sh
braid scotch --bam sample.bam --gtf genes.gtf --out scotch_out
braid scotch-dtu --a scotch_out/sampleA --b scotch_out/sampleB --out dtu.tsv
```

Repeat `--bam` for more than one sample. A path may be a file or a directory of `*.bam` and `*.sam` files; a directory is one sample. Samples share the annotation and discover novel isoforms together. Counts are written per sample.

`--platform` is `10x-ont` (cell tag `CB`, UMI tag `UB`), `10x-pacbio` (`CB` and `XM`), `parse-ont` (barcode fields in the read name), or `bulk` (every read in a sample is one cell). `--barcode-cell` and `--barcode-umi` override the 10x tags. The longest alignment is kept for each cell and UMI. `--fasta` rejects a called poly(A) or poly(T) tail when the genome has a homopolymer at the alignment end. `--update-gtf` splits reference sub-exons using coverage before assignment. `--workers` is accepted and does not run the quantification in parallel.

A read that matches no annotated isoform can seed a novel isoform. SCOTCH keeps that isoform only when a discovery chunk assigns at least 10 reads. Shorter novel copies are grouped into a longer isoform when the extra exons sit on the truncated end, unless `--no-group-novel` is set. `--novel-read-n` and `--novel-read-pct` drop weakly supported novels after that grouping; their reads are counted as `uncategorized_novel`.

## Pooled discovery

`braid discover` pools every cell in one coordinate-sorted genome BAM and writes `transcripts.gtf` plus `command.log`. `--algo` is `tama`, `flair`, `isoquant`, `bambu`, or `stringtie`. `--genome` is the FASTA that matches the BAM. `--gtf` is optional for TAMA, FLAIR, and StringTie, and required for IsoQuant and Bambu. `--threads` defaults to 1.

`flair` is the native precise-collapse step (`collapse_isoforms_precise.py` behavior): alignments become BED12, then isoforms. It does not run minimap2. `tama` uses the existing collapse defaults and converts the BED12 models to GTF. `stringtie` runs long-read assembly (`stringtie -L`). IsoQuant discovery is the bulk transcript-discovery run. Bambu discovery calls `bambu()` with discovery on and quantification off.

Tool-specific files remain under `<out>/raw/`.

## Per-cell quantification

`braid quant` reads a genome BAM that already has cell-barcode and UMI tags, plus an annotation GTF, and writes `counts.tsv`:

```text
barcode    transcript_id    count
```

`--algo` is `scotch`, `isoquant`, `isosceles`, or `bambu`. `--barcode-tag` and `--umi-tag` default to `CB` and `UB`. A repeated cell-barcode and UMI pair is one molecule. Raw tool output stays under `<out>/raw/`.

`scotch` is the existing quantifier and still writes novel transcripts to `<out>/annotation.gtf`. `isoquant` uses `--barcoded_bam` and does not discover novel genes in this mode. `bambu` calls Bambu-Clump `bambu.singlecell`. That build reads `CB` and `UB`, and its own UMI deduplication runs only after a chromosome has more than 100 distinct UMIs, so the command first keeps one primary alignment per cell and UMI. Isosceles assigns reads to the reference transcripts and runs its EM step after the same per-cell UMI filter, because Isosceles counts alignments and does not read a UMI tag. IsoQuant and SCOTCH deduplicate UMIs themselves.

## StringTie merge

```sh
braid stringtie-merge -o merged.gtf sample1.gtf sample2.gtf
```

An optional `-G guide.gtf` is passed through to `stringtie --merge`. `braid merge --algo` still accepts only `tama`.

## Pinned external programs

| Tool | Pin | How braid runs it |
| --- | --- | --- |
| FLAIR precise collapse | BrooksLabUCSC/flair `573414c551332bf6348a9d04ba7cf562f67416cb` | In-process Rust, checked against `tests/parity/flair_collapse` |
| IsoQuant | v4.0.0, commit `7d8268918a770b8d0c6925e8cea99d6c969b4eae` (GPL-2.0-only) | `isoquant` subprocess. Sources are not in this tree |
| StringTie | AI4S-YB/stringtie-rust `743b9710b421f4409aaf9ad0fe7543e55dab116b` (StringTie 3.0.3, MIT) | `stringtie` subprocess |
| Bambu / Bambu-Clump | Pipeline GoekeLab/bambu-singlecell-spatial `aa17818929fb1f64029c8c6d46d112c1d2c9488e`. `bambu.singlecell` is GoekeLab/bambu `d704164f0e20c7fe3fe98ba0921e0893cbe3613f` (package 3.11.1, from `ghcr.io/goekelab/bambu-pipe-bambu:1.0.0`). The installed copy keeps a one-row equivalence class as a list so dplyr 1.2 can join it | `Rscript` in the `braid-lr` environment |
| Isosceles | Genentech/Isosceles `f8f8c0bb449ca3e55e6ad2c6e5a1b33070c1d387` (0.2.1) | `Rscript` in the `braid-lr` environment |

Resolution order is `BRAID_ISOQUANT`, `BRAID_STRINGTIE`, `BRAID_RSCRIPT`, or `BRAID_SAMTOOLS`, then `~/.local/bin`, then `~/miniforge3/envs/braid-lr/bin`. `Rscript` does not fall back to `PATH`. The other programs do.

Each sample directory contains `count_matrix/gene_counts.csv`, `count_matrix/transcript_counts.csv`, `count_matrix/gene_transcript.tsv`, and `auxiliary/assignments.tsv`. The output directory also contains `annotation.gtf`: the reference records plus novel transcripts. Gene p-values from `scotch-dtu` are Holm-adjusted. Transcript and DTU gene p-values are Benjamini-Hochberg adjusted. A transcript test requires at least 20 cells and 20 total counts in each group.

## Compatibility and validation

| Algorithm | Upstream reference |
| --- | --- |
| TAMA | Commit `2fa3c308282190c413e9bf0e0b49e63086eef7d4`; collapse date `2023_03_28` |
| FLAIR combine and precise collapse | Commit `573414c551332bf6348a9d04ba7cf562f67416cb` |
| IsoQuant | Version 4.0.0, commit `7d8268918a770b8d0c6925e8cea99d6c969b4eae` |
| StringTie | stringtie-rust commit `743b9710b421f4409aaf9ad0fe7543e55dab116b` (StringTie 3.0.3) |
| Bambu-Clump pipeline | GoekeLab/bambu-singlecell-spatial commit `aa17818929fb1f64029c8c6d46d112c1d2c9488e` |
| Bambu `bambu.singlecell` | GoekeLab/bambu commit `d704164f0e20c7fe3fe98ba0921e0893cbe3613f` (package 3.11.1) |
| Isosceles | Commit `f8f8c0bb449ca3e55e6ad2c6e5a1b33070c1d387` (package 0.2.1) |
| TACO | Version `0.7.3`, commit `eeaeb879b8622365123edbc61ebc100d84194b80` |
| SCOTCH | Commit `15d6ad8b6319cf806c36f677ffe3ebcc003cae16` |

Tests compare output files against frozen fixtures and exercise CLI compatibility, SAM/BAM equivalence, grouping modes, and regression cases. This is evidence for the covered cases, not a guarantee of equivalence for every input.

TAMA reproduces CPython 2.7 dictionary traversal with `PYTHONHASHSEED=0` and Python 2 numeric formatting. Some mode fixtures were generated through a development-only Python 2 semantic adapter; see [their provenance](tests/parity/tama_modes/ORIGIN.txt). FLAIR fixture provenance is recorded in [UPSTREAM.txt](tests/parity/flair_combine/UPSTREAM.txt).

TACO uses deterministic sorted neighbor traversal rather than reproducing a particular Python 2 hash seed. Tied path choices may therefore differ from upstream. Console logs and upstream error-path exit codes are not covered by file parity.

SCOTCH keeps the lowest community id when Louvain moves are tied, and the lowest gene id when gene assignments are tied. Upstream breaks those ties at random. A read that stays compatible with more than one known isoform is assigned by exon distance instead of a random draw. Novel discovery clusters every unassigned read in chunks of 1500 and does not draw the upstream random subsample. PacBio reads are treated as polyadenylated unless `--fasta` marks the tail as internal priming.

## Development

Run from the workspace root:

```sh
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --locked
```

BAM integration tests skip when `samtools` is unavailable. To require them, as CI does:

```sh
BRAID_REQUIRE_SAMTOOLS=1 cargo test --workspace --locked
```

| Directory | Responsibility |
| --- | --- |
| `crates/braid-cli` | CLI parsing and command dispatch |
| `crates/braid-model` | Shared transcript types and algorithm traits |
| `crates/braid-io` | FASTA/SAM readers, BAM decoding, BED formatting, and path utilities |
| `crates/braid-tama` | TAMA algorithms and Python 2 compatibility semantics |
| `crates/braid-flair` | FLAIR combine and precise collapse |
| `crates/braid-taco` | TACO assembly and graph algorithms |
| `crates/braid-scotch` | SCOTCH quantification and differential transcript usage |
| `tests/parity` | Frozen input/output fixtures and provenance |
| `oracle` | Development scripts for upstream comparisons |

Shared `Transcript` coordinates use 1-based starts and exclusive ends in that same coordinate system. TAMA BED formatting subtracts one from both; algorithm-local coordinate conventions are documented in their modules. Preserve these boundaries and upstream-specific ordering/formatting when refactoring.

## License

Workspace packages declare **GPL-3.0-only**. The bundled FLAIR fixtures retain their [BSD-3-Clause license](tests/parity/flair_combine/FLAIR-LICENSE.txt). The SCOTCH algorithm follows WGLab/SCOTCH, which is MIT licensed; this repository reimplements that behavior and does not copy the upstream sources.
