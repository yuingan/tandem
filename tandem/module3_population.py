"""Module 3: how common known duplications are in a population.

Each junction in the metadata is measured the way the isolate scan measures
the junctions it finds (junction_scan.py): reads crossing the junction
against reads running straight through the same boundary of the reference,
both counted by the same rule, each corrected for the read placements the
scan could count at that spot. The ratio estimates f(n-1), the extra copies
per cell. For a duplication that adds one copy it is the share of cells that
carry it.

No read mapping is needed; the reads are placed with unique 20-mers.

Usage: tandem -r reference.fna -i R1.fq -I R2.fq -pop -m metadata.tsv
"""

import json
import logging
import os

from . import utils
from .junction import classify_single_copy_junction
from .junction_scan import MIN_ANCHOR, MIN_SIZE, scan
from .kmer_index import KmerIndex

logger = logging.getLogger("tandem")

COLUMNS = [
    "name", "seq_id", "provided_start", "provided_end",
    "dup_start", "dup_end", "dup_size",
    "junction_reads", "reference_reads", "dosage", "dosage_raw",
    "junction_visible", "reference_visible", "depth_dosage",
    "flanking_repeat", "repeat_len", "repeat_identity",
    "microhomology_bp", "microhomology_seq",
]


def run_module3(ref_fasta, read1, read2=None, metadata_path=None,
                output_dir=".", threads=1, flank=20,
                hr_complexity_filter=True, min_anchor=MIN_ANCHOR):
    """Measure every junction in the metadata file in one pass over the reads.

    Coordinates in the metadata may be a few bases off (tools disagree, and
    microhomology lets a junction sit anywhere along it): each is matched to
    the best-supported junction within `flank` bases of both ends. A junction
    no read crosses is still reported, with zero junction reads.
    """
    logger.info("=" * 60)
    logger.info("Module 3: Population junction quantification")
    logger.info("=" * 60)

    output_dir = utils.ensure_dir(output_dir)
    if metadata_path is None:
        logger.error("Module 3 needs a metadata file (-m) with the junctions to measure.")
        return []

    entries = utils.parse_metadata(metadata_path)
    if not entries:
        logger.error("No valid entries in metadata file.")
        return []
    if any(e.get("sequence") for e in entries):
        logger.info("  Junction sequences in the metadata are not needed and are ignored")

    logger.info(f"Loading reference: {ref_fasta}")
    sequences, _ = utils.load_fasta(ref_fasta)
    # the metadata carries no sequence names; junctions are on the main contig
    seq_id = max(sequences, key=lambda s: len(sequences[s]))
    ref = sequences[seq_id]
    logger.info(f"  Junctions are measured on {seq_id} ({len(ref):,} bp)")

    logger.info("Building reference index")
    index = KmerIndex(sequences)

    targets = [(e["user_start_1based"], e["user_end_1based"]) for e in entries]
    smallest = min(end - start + 1 for start, end in targets)
    logger.info(f"Measuring {len(targets)} junction(s)")
    hits = scan(ref, [read1, read2], min_anchor=min_anchor,
                min_size=min(MIN_SIZE, smallest), index=index,
                contig_start=index._seq_offsets[seq_id], threads=threads,
                targets=targets, slop=flank)

    results = []
    for entry, hit in zip(entries, hits):
        junction = classify_single_copy_junction(
            ref, hit["start"] - 1, hit["end"], hr_complexity_filter=hr_complexity_filter)
        results.append({
            "name": entry["name"],
            "seq_id": seq_id,
            "provided_start": entry["user_start_1based"],
            "provided_end": entry["user_end_1based"],
            "dup_start": hit["start"],
            "dup_end": hit["end"],
            "dup_size": hit["size"],
            "junction_reads": hit["junction_reads"],
            "reference_reads": hit["reference_reads"],
            "dosage": hit["dosage"],
            "dosage_raw": hit["dosage_raw"],
            "junction_visible": hit["junction_visible"],
            "reference_visible": hit["reference_visible"],
            "depth_dosage": hit["depth_dosage"],
            "flanking_repeat": junction["flanking_repeat"],
            "repeat_len": junction["repeat_len"],
            "repeat_identity": junction["repeat_identity"],
            "microhomology_bp": junction["microhomology_bp"],
            "microhomology_seq": junction["microhomology_seq"],
        })

    for r in results:
        moved = ""
        if (r["dup_start"], r["dup_end"]) != (r["provided_start"], r["provided_end"]):
            moved = f" (given {r['provided_start']:,}-{r['provided_end']:,})"
        logger.info(
            f"  {r['name']}: {r['dup_start']:,}-{r['dup_end']:,}{moved}  "
            f"junction={r['junction_reads']} reference={r['reference_reads']}  "
            f"dosage={r['dosage']} (raw {r['dosage_raw']})  depth={r['depth_dosage']}"
        )

    _save_results(results, output_dir)
    return results


def _save_results(results, output_dir):
    tsv_path = os.path.join(output_dir, "junction_quantification.tsv")
    with open(tsv_path, "w") as f:
        f.write("# dosage = junction reads / reference reads, each corrected for the read\n")
        f.write("# placements the scan can count there: f(n-1), extra copies per cell.\n")
        f.write("\t".join(COLUMNS) + "\n")
        for r in results:
            f.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in COLUMNS) + "\n")
    logger.info(f"  Results saved to {tsv_path}")

    json_path = os.path.join(output_dir, "junction_quantification.json")
    with open(json_path, "w") as f:
        json.dump({
            "unit": "dosage",
            "unit_description": (
                "Junction reads over reference reads at the same boundary, each "
                "corrected for the read placements that can be counted there. "
                "Estimates f(n-1), the extra copies per cell; for a duplication "
                "adding one copy, the share of cells that carry it."
            ),
            "results": results,
        }, f, indent=2)
    logger.info(f"  JSON results saved to {json_path}")
