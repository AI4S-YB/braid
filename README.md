# braid

Reconcile transcript models with Rust implementations of **TAMA collapse**, **TAMA merge**, **FLAIR combine**, and **TACO meta-assembly**.

The project focuses on reproducing upstream algorithm behavior and output formats. Algorithms run natively in Rust; Python is used only by development oracle scripts. Reading BAM files requires `samtools` on `PATH`.

## Commands

| Command | Input | Main output |
| --- | --- | --- |
| `braid tama-collapse` | Sorted SAM/BAM and genome FASTA | Collapsed BED12 models and read reports |
| `braid tama-merge` | Manifest of annotated BED12 files | Merged BED12 models and source reports |
| `braid flair-combine` | Manifest of FLAIR transcriptomes | Combined BED12, counts, and isoform map |
| `braid taco` | Manifest of sample GTF files | Assembled GTF/BED and diagnostic tracks |

The equivalent general commands are `collapse --algo tama`, `merge --algo tama`, and `combine --algo flair`. These are currently the only algorithms available for those commands.

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

## Compatibility and validation

| Algorithm | Upstream reference |
| --- | --- |
| TAMA | Commit `2fa3c308282190c413e9bf0e0b49e63086eef7d4`; collapse date `2023_03_28` |
| FLAIR combine | Commit `573414c551332bf6348a9d04ba7cf562f67416cb` |
| TACO | Version `0.7.3`, commit `eeaeb879b8622365123edbc61ebc100d84194b80` |

Tests compare output files against frozen fixtures and exercise CLI compatibility, SAM/BAM equivalence, grouping modes, and regression cases. This is evidence for the covered cases, not a guarantee of equivalence for every input.

TAMA reproduces CPython 2.7 dictionary traversal with `PYTHONHASHSEED=0` and Python 2 numeric formatting. Some mode fixtures were generated through a development-only Python 2 semantic adapter; see [their provenance](tests/parity/tama_modes/ORIGIN.txt). FLAIR fixture provenance is recorded in [UPSTREAM.txt](tests/parity/flair_combine/UPSTREAM.txt).

TACO uses deterministic sorted neighbor traversal rather than reproducing a particular Python 2 hash seed. Tied path choices may therefore differ from upstream. Console logs and upstream error-path exit codes are not covered by file parity.

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
| `crates/braid-flair` | FLAIR combine |
| `crates/braid-taco` | TACO assembly and graph algorithms |
| `tests/parity` | Frozen input/output fixtures and provenance |
| `oracle` | Development scripts for upstream comparisons |

Shared `Transcript` coordinates use 1-based starts and exclusive ends in that same coordinate system. TAMA BED formatting subtracts one from both; algorithm-local coordinate conventions are documented in their modules. Preserve these boundaries and upstream-specific ordering/formatting when refactoring.

## License

Workspace packages declare **GPL-3.0-only**. The bundled FLAIR fixtures retain their [BSD-3-Clause license](tests/parity/flair_combine/FLAIR-LICENSE.txt).
