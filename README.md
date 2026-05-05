# Tandem: Tandem Amplification aNd Duplication Event Mapper

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

Detection, mechanism classification, and population quantification of bacterial tandem gene duplications.

## Overview

Tandem provides three complementary analysis modes in a single command-line tool:

- **Module 1 — Reference genome.** Detect tandem duplications via NUCmer self-alignment, optionally subcategorize by CDS content, and classify each by HR signature (flanking repeat alignment) and microhomology length at the junction.
- **Module 2 — Isolate.** A two-step workflow: first generate an interactive coverage plot of resequencing data; the user inspects it to identify candidate duplication boundaries; then re-run with those coordinates to confirm junctions and classify mechanism.
- **Module 3 — Population.** Quantify the population-level frequency of known duplication junctions from pooled sequencing. Reports a Junction Read Ratio (JRR) per junction.

The unifying design choice is to report the **two sequence-level observations** that license mechanism inference (HR signature, junction microhomology) rather than collapsing each duplication into a single pathway label. Users and downstream tools apply thresholds appropriate to their biological question.

## Installation

### From source

```bash
git clone https://github.com/PLACEHOLDER/tandem.git
cd tandem
pip install .
```

Tested on Python 3.8–3.12. Python dependencies (Biopython ≥ 1.79, NumPy ≥ 1.20) install automatically.

### External dependencies

```bash
conda install -c bioconda art -y 
conda install -c bioconda bowtie2 -y
conda install -c bioconda mummer4 minimap2 samtools
```

| Tool        | Required by | Why |
|-------------|-------------|-----|
| MUMmer4 (`nucmer`) | Module 1 | self-alignment for tandem candidate detection |
| minimap2    | Modules 2, 3 | read mapping for coverage and junction confirmation |
| samtools    | Modules 2, 3 | coverage parsing |
| bowtie2     | _not required by Tandem_ | only needed if comparing with breseq |

## Quick start

### Module 1 — reference genome

```bash
# Basic
tandem -r genome.fna -o results/ -t 8

# Linear chromosomes (e.g. Streptomyces)
tandem -r genome.fna -o results/ -t 8 --no-circular

# With CDS subcategorization (intergenic / tandem_single_gene / tandem_segmental)
tandem -r genome.fna --annotation genome.gff -o results/ -t 8

# Use coordinates from another tool (skip NUCmer)
tandem -r genome.fna --detection-input external_dups.tsv -o results/ -t 8
```

### Module 2 — isolate junction discovery (two-step workflow)

**Step 1: generate a coverage plot.** Run `-iso` with reads but without `-s`/`-e`:

```bash
tandem -r genome.fna -i reads_R1.fq.gz -I reads_R2.fq.gz -iso -o results/ -t 8
```

This maps reads, computes coverage, suggests candidate elevated regions automatically, and writes an interactive HTML coverage plot to `results/coverage_plot.html`. Tandem then exits without confirming junctions.

**Step 2: inspect the plot, then re-run with chosen coordinates.** Open `coverage_plot.html` in a browser, identify candidate duplication boundaries, then call again with `-s` (start positions) and `-e` (end positions). Multiple values for each test all start × end combinations:

```bash
tandem -r genome.fna -i reads_R1.fq.gz -I reads_R2.fq.gz -iso \
       -s 1854339 1856867 \
       -e 2866785 2860261 \
       --precise -flank 400 \
       -o results/ -t 8
```

Tandem confirms each junction by exact-match read counting and writes `confirmed_junctions.tsv` and `junctions_metadata.tsv` (the latter ready for Module 3). The `--precise` flag tells Tandem the coordinates are exact; without it, Tandem expands the search by `--window` bp on each side.

If you only want the coverage plot (e.g. to inspect coverage before deciding on `-iso`), use `--coverage-only`:

```bash
tandem -r genome.fna -i reads_R1.fq.gz -I reads_R2.fq.gz --coverage-only -o results/ -t 8
```

### Module 3 — population quantification

```bash
tandem -r genome.fna -i pop_R1.fq.gz -I pop_R2.fq.gz \
       -pop -m junctions.tsv -o results/ -t 8
```

The metadata file is typically `junctions_metadata.tsv` produced by Module 2. You can also write it by hand:

```
# Name    Start    End    Sequence(optional)
NHEJ_1    10000    25000
HR_1      50000    68000    ACGTACGT...
```

Module 3 searches a small neighbourhood (default ± 20 bp) around each metadata entry to locate the exact junction position, tolerating coordinate uncertainty from microhomology or different reporting conventions.

## Junction classification

Each duplication is classified by **two independent features**:

### 1. HR signature

For Module 1 (two real copies in the reference), Tandem tests both R-B-R-B-R bracketing scenarios — `[R-B][R-B]-R` and `R-[B-R][B-R]` — and aligns each copy's R region against the external R using a sliding-window HSP analysis. A consecutive match ≥ 35 bp at ≥ 92% identity flags HR.

For Modules 2 and 3 (single copy of the duplicated region Y in the reference), Tandem searches for an **R-Y-R structure** flanking Y, testing three scenarios and reporting the longest valid match:

| Scenario | Geometry | When |
|----------|----------|------|
| `s3` (primary) | R outside Y at both boundaries | Classical R-Y-R: rRNA operons or IS elements flanking the duplicated region |
| `s1` | R outside at start, partial overlap at end | Boundary at end cuts internally through a repeat |
| `s2` | Partial overlap at start, R outside at end | Boundary at start cuts internally through a repeat |

The outward search window scales with duplication size: 200 bp for < 10 kb, 2 kb for 10–50 kb, 5 kb for ≥ 50 kb. This catches flanking rRNA/IS that sits kilobases outside the duplication boundary.

A complexity filter rejects HR calls when both flanking windows are low-complexity (microsatellite periodicity), distinguishing replication slippage from true HR. Disable with `--no-hr-complexity-filter`.

### 2. Microhomology at the junction

For Module 1: the maximum exact overlap between the end of copy 1 and the start of copy 2.

For Modules 2 and 3: the longest exact match between `ref[Y_end-k : Y_end]` and `ref[Y_start : Y_start+k]` — the biological junction sequence created when copy 1's end is joined directly to copy 2's start.

Microhomology is reported as a raw bp count. Tandem does **not** bin into NHEJ/MMEJ/SSA categories — users apply thresholds appropriate to their biological question.

## Output

Each run produces:

- **TSV** result files (one per module's primary output)
- **JSON** with the same content for programmatic parsing
- **HTML report** (`tandem_report.html`) with a click-through summary table
- For Module 2: an interactive coverage plot (`coverage_plot.html`)

### Module 1 columns

| Column | Description |
|--------|-------------|
| `copy1_start`, `copy1_end` | First copy coordinates (1-based inclusive) |
| `copy2_start`, `copy2_end` | Second copy coordinates (overlap-resolved) |
| `distance` | Inter-copy distance, computed as circular distance if `--circular` |
| `nucmer_overlap_bp` | Original NUCmer alignment overlap before resolution |
| `is_hr_signature` | True if flanking repeat signature consistent with HR |
| `hr_match_len`, `hr_identity`, `hr_scenario` | HR alignment details |
| `microhomology_bp`, `microhomology_seq` | Junction microhomology |
| `classification_note` | Flags such as `junction_distant`, `too_small_for_HR_check` |
| `content_category` (only with `--annotation`) | `intergenic` / `tandem_single_gene` / `tandem_segmental` / `proximal_*` |

### Module 2 columns

Junction sequence, read support (spanning, HQ), HR signature status with scenario label, and microhomology. Also generates `junctions_metadata.tsv` ready for Module 3.

### Module 3 — Junction Read Ratio (JRR)

```
JRR = spanning_reads / (WT_coverage × (read_length − junction_length + 1) / read_length)
```

JRR approximates the fraction of cells in the population carrying the duplication. It is best suited for comparing relative abundance across time-series samples from the same population.

## Parameters

### Mode selection

| Flag | Description |
|------|-------------|
| (no flag) | Module 1 (reference analysis) |
| `-iso` | Module 2 (isolate) |
| `-pop` | Module 3 (population, requires `-m`) |
| `--coverage-only` | Coverage plot only — map reads, plot, exit |

### General

| Flag | Default | Description |
|------|---------|-------------|
| `-r` | required | Reference genome FASTA (plain or `.gz`) |
| `-i`, `-I` | — | Forward / reverse reads FASTQ |
| `-m` | — | Junction metadata file (Module 3) |
| `-t` | 1 | Threads |
| `-o` | `tandem_output` | Output directory |
| `-v` / `-vv` | — | Verbose / debug output |
| `--version` | — | Print version and exit |

### Module 1 (Reference)

| Flag | Default | Description |
|------|---------|-------------|
| `-d` | 30000 | Max distance between tandem copies (bp) |
| `--min-identity` | 80.0 | Minimum NUCmer alignment identity (%) |
| `--min-size` | 200 | Minimum duplication size (bp) |
| `-c` / `--circular` | on | Treat genome as circular |
| `--no-circular` | — | Treat genome as linear |
| `--flag-distance` | 2000 | Distance threshold for `junction_distant` flag |
| `--annotation` | — | GFF3 or GenBank file for CDS subcategorization (Stage 2) |
| `--detection-input` | — | TSV of pre-computed coordinates from another tool (skips NUCmer) |
| `--no-short-dup-filters` | — | Disable short-duplication (< 500 bp) filters |

### Module 2 (Isolate)

| Flag | Default | Description |
|------|---------|-------------|
| `-s` | — | Candidate start position(s), 1-based inclusive (multiple allowed) |
| `-e` | — | Candidate end position(s), 1-based inclusive (multiple allowed) |
| `--precise` | — | Treat `-s`/`-e` as exact (no window expansion) |
| `-w` / `--window` | 200 | Sliding window for coverage detection (bp) |
| `-flank` / `--flank` | 200 | Flank around junction for read-mapping reference (bp) |
| `--max-k` | 75 | Half-length of junction probe (bp). Increase for long reads (e.g. 500 for PacBio/Nanopore) |
| `--min-fold` | 1.7 | Minimum fold change above median coverage to flag elevated regions |
| `--merge-fold` | 1.3 | Fold-change threshold for merging adjacent elevated regions |
| `--merge-distance` | 3000 | Fixed merge distance (bp); regions within this are merged |

### Module 3 (Population)

| Flag | Default | Description |
|------|---------|-------------|
| `-flank` / `--flank` | 20 | Neighbourhood (± bp) around metadata coordinates to search for the exact junction. Defaults to 20 in Module 3; explicit user values are honoured |

### HR-detection sensitivity (all modules)

| Flag | Default | Description |
|------|---------|-------------|
| `--hr-min-consec` | 35 | Min consecutive match for HR signature (bp) |
| `--hr-min-identity` | 0.92 | Min window identity |
| `--hr-max-inward` | 2000 | Inward search window cap (bp) |
| `--hr-outward-ext` | 200 | Outward search window (bp); auto-scaled for ≥ 10 kb dups |
| `--no-hr-complexity-filter` | — | Disable microsatellite-rejection filter |

### `tandem --help` snapshot

```
$ tandem --help
usage: tandem [-h] -r REFERENCE [-i READ1] [-I READ2] [-iso] [-pop]
              [--coverage-only] [-m METADATA] [-s START [START ...]]
              [-e END [END ...]] [-w WINDOW] [-flank FLANK] [--precise]
              [--max-k MAX_K] [--min-fold MIN_FOLD] [--merge-fold MERGE_FOLD]
              [--merge-distance MERGE_DISTANCE] [-t THREADS]
              [-d MAX_DISTANCE] [--min-identity MIN_IDENTITY] [--min-size MIN_SIZE]
              [-c] [--no-circular] [--flag-distance FLAG_DISTANCE]
              [--detection-input DETECTION_INPUT] [--no-short-dup-filters]
              [--annotation ANNOTATION]
              [--hr-min-consec HR_MIN_CONSEC] [--hr-min-identity HR_MIN_IDENTITY]
              [--hr-max-inward HR_MAX_INWARD] [--hr-outward-ext HR_OUTWARD_EXT]
              [--no-hr-complexity-filter]
              [-o OUTPUT] [-v] [--version]
```

Full descriptions: `tandem --help`.

## Typical workflow

```bash
# 1. Characterise the reference genome
tandem -r ref.fna -o step1_reference/ -t 8

# 2a. Generate coverage plot for evolved isolate
tandem -r ref.fna -i isolate_R1.fq.gz -I isolate_R2.fq.gz \
       -iso -o step2_isolate/ -t 8
# → inspect step2_isolate/coverage_plot.html

# 2b. Re-run with chosen coordinates from the plot
tandem -r ref.fna -i isolate_R1.fq.gz -I isolate_R2.fq.gz \
       -iso -s 1854000 -e 2866000 --precise -flank 400 \
       -o step2_isolate/ -t 8

# 3. Quantify junctions in population samples
tandem -r ref.fna -i pop_R1.fq.gz -I pop_R2.fq.gz \
       -pop -m step2_isolate/junctions_metadata.tsv \
       -o step3_population/ -t 8
```

## Performance and validation

| Module | Sample | Runtime |
|--------|--------|---------|
| 1 | *M. genitalium* (0.58 Mb) | 1.2 s |
| 1 | *P. fluorescens* SBW25 (6.7 Mb) | 10.3 s |
| 1 | *M. tuberculosis* H37Rv (4.4 Mb) | 5.0 s |
| 2 | SBW25 isolate (~6 Mb, ~50× coverage) | 30–60 min on 44 threads |
| 3 | SBW25 day-28 population sample | 2.5–3.2 min on 44 threads |

Module 1 runtimes are dominated by NUCmer self-alignment, which is largely single-threaded. Modules 2 and 3 use thread parallelism for read mapping and junction confirmation.

**Module 3 validation** (artificial populations, in silico read mixing across six bacterial species):

| Species | n tests | Recall | MAE |
|---------|---------|--------|-----|
| *A. baumannii* ATCC17978 | 40 | 100% | 0.041 |
| *E. coli* K12 MG1655     | 40 | 100% | 0.052 |
| *M. pneumoniae* M129     | 40 | 100% | 0.042 |
| *M. tuberculosis* H37Rv  | 40 | 100% | 0.045 |
| *P. aeruginosa* PAO1     | 40 | 100% | 0.036 |
| *S. coelicolor* A3(2)    | 40 | 100% | 0.041 |
| **Aggregate** | **240** | **100%** | **0.043 (4.3%)** |

## Limitations

- Detects forward-orientation flanking repeats only; inverted-repeat-mediated inversion events are out of scope.
- Cannot detect non-templated insertions at junctions (LigD-class NHEJ in Ku/LigD-positive species).
- Static-genome analysis cannot detect duplications that have already reverted; only time-course data captures reversion.
- The complexity filter excludes replication-slippage tandem repeats but may also reject HR mediated by simple-repeat substrates such as REP elements; disable with `--no-hr-complexity-filter`.
- Tandem requires the R-Y-R geometry to be satisfied on both sides of the duplicated region; one-sided flanking repeats will not produce an HR call even when one boundary clearly sits in a repeat.
- Circular-genome wraparound is handled for distance metrics but not for HR detection windows at genome edges.

## Citation

If you use Tandem, please cite:

```bibtex
XXXXXXXXXXXXX

```

## License

MIT — see [LICENSE](LICENSE).
