"""Command-line interface for Tandem.

Usage:
    tandem -r reference.fna                                    # Module 1: reference analysis
    tandem -r reference.fna -i R1.fq -I R2.fq -iso            # Module 2: coverage plot (step 1)
    tandem -r reference.fna -i R1.fq -I R2.fq -iso \\
           -s 1847972 1848410 -e 2856714 2866983 \\
           --precise -flank 400                                # Module 2: junction confirmation (step 2)
    tandem -r reference.fna -i R1.fq -I R2.fq -pop -m meta.tsv # Module 3: population quantification

Module 2 workflow:
    Step 1: Run with -iso (no -s/-e) to generate interactive coverage plot.
    Step 2: Examine the HTML coverage plot to identify candidate boundaries.
    Step 3: Re-run with -iso -s <starts> -e <ends> --precise -flank <flank>.

Options:
    -r          Reference genome FASTA (required)
    -i          Forward reads FASTQ (modules 2 & 3)
    -I          Reverse reads FASTQ (optional, paired-end)
    -iso        Isolate mode (module 2)
    -pop        Population mode (module 3)
    -m          Metadata file for module 3 (Name, Start, End, [Sequence])
    -s          Candidate start position(s), 1-based inclusive (module 2)
    -e          Candidate end position(s), 1-based inclusive (module 2)
    --precise   Exact position search (no window addition around -s/-e)
    -flank      Flank size for junction generation (default: 200)
    -t          Number of threads (default: 1)
    -o          Output directory (default: tandem_output)
    -v          Verbose output (-vv for debug)
"""

import argparse
import os
import sys
import time

from . import __version__
from . import utils
from .report import generate_report


def parse_args(argv=None):
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        prog="tandem",
        description=(
            "Tandem: Detection and classification of tandem duplication "
            "mechanisms in bacterial genomes"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Module 1: Detect tandem duplications in a reference genome
  tandem -r genome.fna -o results/

  # Module 2: Discover junctions from isolate resequencing
  tandem -r genome.fna -i reads_R1.fq -I reads_R2.fq -iso -t 8

  # Module 2: With manual region coordinates
  tandem -r genome.fna -i reads_R1.fq -I reads_R2.fq -iso -s 10000 -e 25000

  # Module 3: Quantify known junctions in population data
  tandem -r genome.fna -i pop_R1.fq -I pop_R2.fq -pop -m junctions.tsv -t 8

  # Module 3 metadata format (tab-separated):
  # Name    Start   End     Sequence(optional)
  # J1      10000   25000
  # J2      50000   68000   ACGTACGT...
        """,
    )

    # Required
    parser.add_argument(
        "-r", "--reference", required=True,
        help="Reference genome FASTA file (required)"
    )

    # Reads
    parser.add_argument(
        "-i", "--read1",
        help="Forward reads FASTQ (required for -iso and -pop)"
    )
    parser.add_argument(
        "-I", "--read2",
        help="Reverse reads FASTQ (optional, for paired-end data)"
    )

    # Mode flags
    parser.add_argument(
        "-iso", "--isolate", action="store_true",
        help="Isolate mode (Module 2): confirm duplication junctions in "
             "resequencing data. Without -s/-e: generates coverage plot for "
             "boundary identification. With -s/-e: tests all start×end "
             "combinations and reports confirmed junctions with mechanism "
             "classification."
    )
    parser.add_argument(
        "-pop", "--population", action="store_true",
        help="Population mode: quantify known junctions (module 3)"
    )
    parser.add_argument(
        "--coverage-only", action="store_true",
        help="Coverage plot mode: map reads, calculate coverage, generate "
             "interactive HTML coverage plot, then stop. Use this to inspect "
             "coverage before running full junction discovery. Requires -r and -i."
    )

    # Module 3 metadata
    parser.add_argument(
        "-m", "--metadata",
        help="Metadata file for -pop mode (TSV: Name, Start, End, [Sequence])"
    )

    # Module 2 manual coordinates (multiple values = test all combinations)
    # Coordinates are 1-based inclusive (paper/NCBI/GenBank convention)
    parser.add_argument(
        "-s", "--start", type=int, nargs="+",
        help="Start position(s) for junction search (1-based inclusive, "
             "matching paper/NCBI convention). Multiple values allowed: "
             "-s 1854339 1856867 (tests all start×end combinations)"
    )
    parser.add_argument(
        "-e", "--end", type=int, nargs="+",
        help="End position(s) for junction search (1-based inclusive, "
             "matching paper/NCBI convention). Multiple values allowed: "
             "-e 2866785 2860261 (tests all start×end combinations)"
    )

    # Parameters
    parser.add_argument(
        "-w", "--window", type=int, default=200,
        help="Sliding window size for coverage detection (default: 200)"
    )
    parser.add_argument(
        "-flank", "--flank", type=int, default=200,
        help="Flank region size for junction search (default: 200)"
    )
    parser.add_argument(
        "--precise", action="store_true",
        help="Precise mode: search range = flank only (no window added). "
             "Use when you know the exact junction position. "
             "Without this flag, search range = window + flank."
    )
    parser.add_argument(
        "--max-k", type=int, default=75,
        help="Maximum half-length of junction probe (default: 75, giving "
             "150bp junction reference for Illumina 150bp reads). "
             "Increase for long reads (e.g. --max-k 500 for PacBio/Nanopore). "
             "Junction reference length = 2 × max-k."
    )
    parser.add_argument(
        "--min-fold", type=float, default=1.7,
        help="Minimum fold change above median coverage to detect elevated "
             "regions (default: 1.7). Lower for sensitivity, higher to reduce noise."
    )
    parser.add_argument(
        "--merge-fold", type=float, default=1.3,
        help="Fold change threshold for merging adjacent elevated regions "
             "(default: 1.3). If coverage in the gap between two regions stays "
             "above median * merge-fold, they are merged into one block. "
             "Lower to merge more aggressively (e.g. 1.1 for large duplications "
             "with uneven coverage)."
    )
    parser.add_argument(
        "--merge-distance", type=int, default=3000,
        help="Fixed merge distance in bp (default: 3000). Regions within this "
             "distance are merged regardless of gap coverage. Applied after "
             "coverage-aware merging."
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=1,
        help="Number of threads (default: 1)"
    )

    # Module 1 specific
    parser.add_argument(
        "-d", "--max-distance", type=int, default=30000,
        help="Maximum distance between copies to consider tandem (default: 30000)"
    )
    parser.add_argument(
        "--min-identity", type=float, default=80.0,
        help="Minimum alignment identity for duplication detection (default: 80.0)"
    )
    parser.add_argument(
        "--min-size", type=int, default=200,
        help="Minimum duplication size in bp (default: 200). Duplications "
             "smaller than 500bp have additional filters applied (see "
             "--short-dup-gap-ratio). Mechanism inference confidence is "
             "reported as low (<300bp), moderate (300-500bp), or high (>=500bp)."
    )
    parser.add_argument(
        "-c", "--circular", action="store_true", default=True,
        help="Treat genomes as circular (default: on). Computes distance "
             "as min(linear, wrap-around through origin)."
    )
    parser.add_argument(
        "--no-circular", action="store_true",
        help="Treat genomes as linear (disables circular distance)"
    )
    parser.add_argument(
        "--flag-distance", type=int, default=2000,
        help="Distance threshold (bp) for 'junction_distant' flag (default: 2000). "
             "Duplications beyond this distance are flagged as having potentially "
             "degraded junction signals."
    )

    # Module 1 Stage 1 — optional external detection input
    parser.add_argument(
        "--detection-input",
        help="Optional TSV with pre-computed tandem duplication coordinates "
             "from an external tool (e.g. SegMantX, BISER, breseq). "
             "Columns: seq_id, copy1_start, copy1_end, copy2_start, copy2_end "
             "(1-based inclusive). Optional columns: size, identity, "
             "orientation. If provided, NUCmer detection is skipped and "
             "tandem only runs Stages 2-3 on the provided coordinates."
    )
    parser.add_argument(
        "--no-short-dup-filters", action="store_true",
        help="Disable short-duplication filters (gap<copy check, low-complexity "
             "masking) for duplications below 500bp. Default: filters on."
    )

    # Module 1 Stage 2 — optional annotation for CDS subcategorization
    parser.add_argument(
        "--annotation",
        help="Optional GFF3 or GenBank file for Stage 2 CDS subcategorization. "
             "When provided, each tandem duplication is annotated with CDS "
             "content (n_cds_copy1, n_cds_gap) and a content_category "
             "(intergenic / tandem_single_gene / tandem_segmental / proximal_*)."
    )

    # Module 1 Stage 3 — HR detection parameter tuning (for sensitivity analysis)
    parser.add_argument(
        "--hr-min-consec", type=int, default=None,
        help="Minimum consecutive match length for HR detection (default: 35). "
             "Lower = more sensitive, higher = more stringent."
    )
    parser.add_argument(
        "--hr-min-identity", type=float, default=None,
        help="Minimum alignment identity for HR detection (default: 0.92). "
             "Lower = accept more diverged repeats."
    )
    parser.add_argument(
        "--hr-max-inward", type=int, default=None,
        help="Cap for inward extent of HR search window (default: 2000). "
             "The window extends up to this many bp into the copy."
    )
    parser.add_argument(
        "--hr-outward-ext", type=int, default=None,
        help="Outward extension of HR search window beyond copy boundary "
             "(default: 200). Tolerates NUCmer boundary uncertainty."
    )
    parser.add_argument(
        "--no-hr-complexity-filter", action="store_true",
        help="Disable the HR complexity filter that rejects HR calls when "
             "both flanking windows are low-complexity (e.g. microsatellite "
             "periodicity). Default: filter is ON, which excludes tandem "
             "repeat expansions from replication slippage. Disable this if "
             "you suspect HR mediated by simple-repeat substrates (e.g. "
             "REP elements)."
    )

    # Output
    parser.add_argument(
        "-o", "--output", default="tandem_output",
        help="Output directory (default: tandem_output)"
    )

    # Verbosity
    parser.add_argument(
        "-v", "--verbose", action="count", default=0,
        help="Increase verbosity (-v for info, -vv for debug)"
    )

    parser.add_argument(
        "--version", action="version",
        version=f"tandem {__version__}"
    )

    args = parser.parse_args(argv)

    # Validate arguments
    _validate_args(args, parser)

    return args


def _validate_args(args, parser):
    """Validate argument combinations."""
    # --- Logical checks first (before file existence) ---
    has_reads = args.read1 is not None

    n_modes = sum([args.isolate, args.population, args.coverage_only])
    if n_modes > 1:
        parser.error("Cannot specify multiple modes. Choose one of: -iso, -pop, --coverage-only.")

    if args.isolate or args.population:
        if not has_reads:
            parser.error(
                f"{'Isolate' if args.isolate else 'Population'} mode "
                f"requires reads (-i). Provide at least forward reads."
            )

    if args.coverage_only:
        if not has_reads:
            parser.error(
                "Coverage-only mode requires reads (-i). "
                "Provide at least forward reads."
            )

    if args.population and not args.metadata:
        parser.error(
            "Population mode (-pop) requires a metadata file (-m). "
            "Run isolate mode (-iso) first to discover junctions, "
            "or provide metadata manually."
        )

    if has_reads and not (args.isolate or args.population or args.coverage_only):
        parser.error(
            "Reads provided but no mode specified. "
            "Use -iso for isolate mode, -pop for population mode, "
            "or --coverage-only for coverage inspection."
        )

    # Manual coordinates check
    if (args.start is not None) != (args.end is not None):
        parser.error("Both -s (start) and -e (end) must be specified together.")

    if args.start is not None and args.end is not None:
        if not args.isolate:
            parser.error("Manual coordinates (-s, -e) are only used in isolate mode (-iso).")
        # Validate 1-based minimum
        for s in args.start:
            if s < 1:
                parser.error(f"Start coordinate {s} is less than 1 (coordinates are 1-based).")
        for e in args.end:
            if e < 1:
                parser.error(f"End coordinate {e} is less than 1 (coordinates are 1-based).")
        # Convert 1-based starts to 0-based internal representation
        # (ends stay the same: 1-based inclusive == Python exclusive)
        args.start = [s - 1 for s in args.start]
        # args.end is unchanged (1-based inclusive = 0-based exclusive in Python slicing)

    # --- File existence checks ---
    if not os.path.isfile(args.reference):
        parser.error(f"Reference file not found: {args.reference}")

    if args.read1 and not os.path.isfile(args.read1):
        parser.error(f"Forward reads not found: {args.read1}")
    if args.read2 and not os.path.isfile(args.read2):
        parser.error(f"Reverse reads not found: {args.read2}")

    if args.metadata and not os.path.isfile(args.metadata):
        parser.error(f"Metadata file not found: {args.metadata}")


def main(argv=None):
    """Main entry point for Tandem."""
    args = parse_args(argv)

    # Setup logging
    utils.setup_logging(args.verbose)

    start_time = time.time()

    print(f"Tandem v{__version__}")
    print(f"{'=' * 60}")

    # Determine which module to run
    if args.coverage_only:
        module_num = -1  # special: coverage only
        mode_name = "Coverage Plot (inspection only)"
    elif args.isolate:
        module_num = 2
        mode_name = "Isolate Junction Discovery"
    elif args.population:
        module_num = 3
        mode_name = "Population Junction Quantification"
    else:
        module_num = 1
        mode_name = "Reference Genome Analysis"

    print(f"Mode: {mode_name}")
    print(f"Reference: {args.reference}")
    if args.read1:
        print(f"Forward reads: {args.read1}")
    if args.read2:
        print(f"Reverse reads: {args.read2}")
    print(f"Output: {args.output}")
    print(f"Threads: {args.threads}")
    if module_num == 2:
        if hasattr(args, 'precise') and args.precise:
            print(f"Flank: {args.flank} bp (precise mode)")
        else:
            print(f"Flank: {args.flank} bp")

    # Handle --no-circular overriding --circular
    is_circular = args.circular and not args.no_circular
    if module_num == 1:
        print(f"Circular genome: {'ON' if is_circular else 'OFF'}")
        print(f"Flag distance: {args.flag_distance} bp")
    print(f"{'=' * 60}\n")

    # Check dependencies
    deps_needed = [2] if module_num == -1 else [module_num]
    utils.check_dependencies(deps_needed)

    # Create output directory
    utils.ensure_dir(args.output)

    # Run appropriate module
    if module_num == -1:
        # Coverage-only mode
        from .coverage_plot import run_coverage_only
        results = run_coverage_only(
            ref_fasta=args.reference,
            read1=args.read1,
            read2=args.read2,
            output_dir=args.output,
            threads=args.threads,
            window=args.window,
            min_fold_change=args.min_fold,
            merge_fold_threshold=args.merge_fold,
            merge_distance=args.merge_distance,
        )
        module_name = None  # no report needed

    elif module_num == 1:
        from .module1_reference import run_module1
        results = run_module1(
            ref_fasta=args.reference,
            output_dir=args.output,
            threads=args.threads,
            max_distance=args.max_distance,
            min_identity=args.min_identity,
            min_size=args.min_size,
            circular=is_circular,
            flag_distance=args.flag_distance,
            annotation=args.annotation,
            apply_short_dup_filters=not args.no_short_dup_filters,
            detection_input=args.detection_input,
            hr_min_consec=args.hr_min_consec,
            hr_min_identity=args.hr_min_identity,
            hr_max_inward=args.hr_max_inward,
            hr_outward_ext=args.hr_outward_ext,
            hr_complexity_filter=not args.no_hr_complexity_filter,
        )
        module_name = "module1"

    elif module_num == 2:
        from .module2_isolate import run_module2
        # When --precise, junction search uses flank only (no window added)
        junction_window = 0 if args.precise else args.window
        results = run_module2(
            ref_fasta=args.reference,
            read1=args.read1,
            read2=args.read2,
            output_dir=args.output,
            threads=args.threads,
            window=args.window,
            flank=args.flank,
            junction_window=junction_window,
            manual_start=args.start,
            manual_end=args.end,
            max_k=args.max_k,
            hr_complexity_filter=not args.no_hr_complexity_filter,
        )
        module_name = "module2"

    elif module_num == 3:
        from .module3_population import run_module3
        # Module 3 uses a smaller flank window (default ±20bp) to tolerate
        # coordinate uncertainty from microhomology (0-25bp typical) or
        # different tools. If user sets --flank explicitly, use that value.
        # Larger values (>50) dramatically slow down the search without
        # meaningful biological benefit.
        if args.flank == 200:  # default for module 2
            module3_flank = 20
        else:
            module3_flank = args.flank
        results = run_module3(
            ref_fasta=args.reference,
            read1=args.read1,
            read2=args.read2,
            metadata_path=args.metadata,
            output_dir=args.output,
            threads=args.threads,
            flank=module3_flank,
            hr_complexity_filter=not args.no_hr_complexity_filter,
        )
        module_name = "module3"

    # Generate HTML report (skip for coverage-only mode)
    if module_name is not None:
        extra_info = {}
        if args.read1:
            extra_info["Forward Reads"] = os.path.basename(args.read1)
        if args.read2:
            extra_info["Reverse Reads"] = os.path.basename(args.read2)
        extra_info["Threads"] = str(args.threads)
        if module_num == 2:
            extra_info["Flank"] = f"{args.flank} bp"
            if hasattr(args, 'precise') and args.precise:
                extra_info["Mode"] = "Precise"

        generate_report(
            results=results,
            module_name=module_name,
            output_dir=args.output,
            ref_name=os.path.basename(args.reference),
            extra_info=extra_info,
        )

    elapsed = time.time() - start_time
    print(f"\n{'=' * 60}")
    if module_name is None:
        print(f"Coverage plot generated.")
        print(f"Time: {elapsed:.1f}s")
        print(f"Output: {args.output}/")
        print(f"Inspect the coverage plot, then run with -iso -s START -e END for targeted junction discovery.")
    else:
        print(f"Tandem complete. {len(results)} result(s).")
        print(f"Time: {elapsed:.1f}s")
        print(f"Output: {args.output}/")
        print(f"Report: {args.output}/tandem_report.html")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
