"""Coverage-only mode for Tandem.

Maps reads to reference, calculates coverage, generates an interactive
coverage plot for visual inspection. No junction discovery.

Usage: tandem -r reference.fna -i R1.fq -I R2.fq --coverage-only -o output/
"""

import logging
import os

from . import utils
from . import alignment
from . import coverage as cov_module

logger = logging.getLogger("tandem")


def run_coverage_only(ref_fasta, read1, read2=None, output_dir=".",
                      threads=1, window=200, min_fold_change=1.7,
                      merge_fold_threshold=1.3, merge_distance=3000):
    """Run coverage-only analysis: map reads, calculate coverage, plot.

    Args:
        ref_fasta: path to reference FASTA
        read1: path to forward reads
        read2: path to reverse reads (optional)
        output_dir: output directory
        threads: number of threads
        window: sliding window size for coverage detection
        min_fold_change: fold change threshold to show on plot
        merge_fold_threshold: fold change for merging adjacent regions
        merge_distance: fixed merge distance in bp

    Returns:
        list of elevated region dicts (for summary)
    """
    logger.info("=" * 60)
    logger.info("Coverage-only mode: generating coverage plot")
    logger.info("=" * 60)

    output_dir = utils.ensure_dir(output_dir)

    # Load reference
    logger.info(f"Loading reference: {ref_fasta}")
    sequences, headers = utils.load_fasta(ref_fasta)
    for seq_id, desc in headers:
        logger.info(f"  {seq_id}: {len(sequences[seq_id]):,} bp")

    # Map reads to reference
    logger.info("Mapping reads to reference")
    bam_path = os.path.join(output_dir, "reads_to_ref.sorted.bam")
    alignment.map_reads_to_reference(
        ref_fasta, read1, read2, output_bam=bam_path, threads=threads
    )

    # Calculate coverage and generate plots per contig
    all_regions = []

    for seq_id in sequences:
        logger.info(f"Calculating coverage for {seq_id}...")
        cov_data = cov_module.calculate_coverage_from_bam(bam_path, seq_id)

        if seq_id not in cov_data or len(cov_data[seq_id]) == 0:
            logger.warning(f"  No coverage data for {seq_id}")
            continue

        # Detect regions using threshold approach
        regions, genome_median = cov_module.detect_elevated_regions(
            cov_data[seq_id],
            window_size=window,
            min_fold_change=min_fold_change,
            merge_fold_threshold=merge_fold_threshold,
            merge_distance=merge_distance,
        )

        # Generate interactive coverage plot
        plot_path = cov_module.generate_coverage_plot(
            cov_data[seq_id],
            seq_id=seq_id,
            output_dir=output_dir,
            window_size=200,
            genome_median=genome_median,
            elevated_regions=regions,
            min_fold_change=min_fold_change,
        )

        logger.info(f"  Coverage plot: {plot_path}")

        if regions:
            for r in regions:
                r["seq_id"] = seq_id
            all_regions.extend(regions)
            logger.info(f"  {len(regions)} elevated regions detected")
            for r in regions:
                logger.info(
                    f"    {r['start']:,}-{r['end']:,} "
                    f"({r['size']:,} bp, {r['fold_change']:.1f}x)"
                )
        else:
            logger.info("  No elevated regions detected")

    logger.info(f"\nTotal elevated regions across all contigs: {len(all_regions)}")

    return all_regions
