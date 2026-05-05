# Changelog

All notable changes to Tandem are documented in this file. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — 2026-04-28

First public release accompanying manuscript submission.

### Added
- **Module 3 (population)** end-to-end implementation: coordinate refinement
  in a ± 20 bp neighbourhood, exact-match junction read counting, and
  Junction Read Ratio (JRR) calculation against minimap2-derived genome-wide
  coverage.
- **Single-copy junction classifier** (`classify_single_copy_junction`) used
  by Modules 2 and 3, with three R-Y-R scenarios:
  - `s3` — repeat outside Y on both sides (primary).
  - `s1`, `s2` — asymmetric, handling boundaries that cut internally
    through a repeat copy.
- **Adaptive outward extension** for HR detection in single-copy mode:
  200 bp default, 2 kb for duplications ≥ 10 kb, 5 kb for ≥ 50 kb. Catches
  flanking rRNA operons and IS elements that sit several kilobases outside
  the boundary.
- **HR sensitivity sweep** (CLI parameters `--hr-min-consec`,
  `--hr-min-identity`, `--hr-max-inward`, `--hr-outward-ext`) and a
  `--no-hr-complexity-filter` switch for users who want to retain
  simple-repeat HR substrates.
- **HTML reports** with click-through junction detail views.
- 23-test self-check (`test_hr_fixes.py`) covering parser correctness,
  windowed-identity behaviour, four-check geometry, complexity filtering,
  microhomology measurement up to 100 bp, and single-copy R-Y-R / no-HR
  controls.

### Fixed
- Parser layer now correctly reads Biopython 1.87 alignment outputs (full
  target + query strings).
- Microhomology detection no longer silently truncates at 60 bp.
- Single-copy MH measurement now compares end-of-Y to start-of-Y rather
  than a virtual second-copy reconstruction (pre-0.2.0 bug that produced
  wrong MH for single-copy junctions).
- HR sliding-window scoring uses strict gap penalties to prevent identity
  dilution by noise.

### Performance
- Module 1: 1.2–10.3 s on genomes from 0.58 Mb (*M. genitalium*) to
  6.26 Mb (*M. tuberculosis*) on 44 threads.
- Module 2: 30–60 min per SBW25 isolate (~6 Mb genome, ~50× coverage).
- Module 3: 2.5–3.2 min per SBW25 day-28 population sample.

### Validation
- Module 1: 8 reference genomes, 152 duplications detected, 23/23
  unit tests pass.
- Module 2: 5 SBW25 isolates (M1–M5), all confirmed junctions match
  published coordinates from Khomarbaghi et al. 2024.
- Module 3: 240 spiked junctions across 6 species, 4 time points each,
  100% recall, 4.3% mean absolute error.

## [0.1.0] — pre-release

Internal alpha. Module 1 + early Module 2 prototype. Not distributed.
