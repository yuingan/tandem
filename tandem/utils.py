"""Utility functions for Tandem."""

import logging
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

logger = logging.getLogger("tandem")


def setup_logging(verbosity=0):
    """Configure logging based on verbosity level."""
    level = logging.WARNING
    if verbosity == 1:
        level = logging.INFO
    elif verbosity >= 2:
        level = logging.DEBUG

    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter(
        "[%(asctime)s] %(levelname)s - %(message)s", datefmt="%H:%M:%S"
    ))
    root = logging.getLogger("tandem")
    root.setLevel(level)
    root.addHandler(handler)


def check_dependency(name, test_cmd=None):
    """Check if an external tool is available."""
    if test_cmd is None:
        test_cmd = [name, "--version"]
    try:
        result = subprocess.run(test_cmd, capture_output=True, timeout=10)
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def check_dependencies(modules_needed):
    """Check all required external dependencies.

    Args:
        modules_needed: list of module numbers (1, 2, 3)
    """
    missing = []

    # Module 1 needs NUCmer
    if 1 in modules_needed:
        if not check_dependency("nucmer", ["nucmer", "--version"]):
            missing.append("nucmer (MUMmer4): conda install -c bioconda mummer4")

    # Modules 2 and 3 need a read mapper and samtools
    if 2 in modules_needed or 3 in modules_needed:
        if not check_dependency("minimap2"):
            if not check_dependency("bwa"):
                missing.append(
                    "minimap2 or bwa: conda install -c bioconda minimap2"
                )
        if not check_dependency("samtools"):
            missing.append("samtools: conda install -c bioconda samtools")

    if missing:
        logger.error("Missing dependencies:")
        for m in missing:
            logger.error(f"  - {m}")
        sys.exit(1)


def load_fasta(fasta_path):
    """Load a FASTA file into a dict of {seq_id: sequence_string}.

    Supports both plain and gzipped (.gz) FASTA files.
    Also returns a list of (seq_id, description) tuples preserving order.
    """
    import gzip

    sequences = {}
    headers = []
    current_id = None
    current_desc = ""
    current_seq = []

    open_func = gzip.open if str(fasta_path).endswith('.gz') else open
    mode = 'rt' if str(fasta_path).endswith('.gz') else 'r'

    with open_func(fasta_path, mode) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences[current_id] = "".join(current_seq).upper()
                    headers.append((current_id, current_desc))
                parts = line[1:].split(None, 1)
                current_id = parts[0]
                current_desc = parts[1] if len(parts) > 1 else ""
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences[current_id] = "".join(current_seq).upper()
            headers.append((current_id, current_desc))

    return sequences, headers


def decompress_if_gzipped(filepath, output_dir=None):
    """Decompress a gzipped file to a temporary location if needed.

    NUCmer and some other tools cannot read .gz files directly.
    This creates an uncompressed copy in output_dir.

    Args:
        filepath: path to input file (may or may not be .gz)
        output_dir: directory for decompressed file (default: same dir)

    Returns:
        (decompressed_path, needs_cleanup) tuple.
        If file was not gzipped, returns (original_path, False).
    """
    import gzip
    import shutil

    filepath = str(filepath)
    if not filepath.endswith('.gz'):
        return filepath, False

    if output_dir is None:
        output_dir = os.path.dirname(filepath) or '.'

    # Strip .gz extension for output name
    basename = os.path.basename(filepath)
    if basename.endswith('.gz'):
        basename = basename[:-3]
    out_path = os.path.join(output_dir, basename)

    logger.info(f"  Decompressing {os.path.basename(filepath)}...")
    with gzip.open(filepath, 'rb') as f_in:
        with open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    return out_path, True


# Supported file extensions
REFERENCE_EXTENSIONS = {
    '.fasta', '.fa', '.fna', '.fas', '.seq',
    '.fasta.gz', '.fa.gz', '.fna.gz', '.fas.gz', '.seq.gz',
}

READ_EXTENSIONS = {
    '.fastq', '.fq', '.fasta', '.fa', '.fna',
    '.fastq.gz', '.fq.gz', '.fasta.gz', '.fa.gz', '.fna.gz',
    '.bam', '.sam',
}


def validate_file_extension(filepath, allowed_extensions, file_type="file"):
    """Check if a file has a supported extension.

    Args:
        filepath: path to the file
        allowed_extensions: set of allowed extensions (including dot)
        file_type: description for error message

    Returns:
        True if valid, logs warning if not (does not block execution)
    """
    name = str(filepath).lower()
    for ext in sorted(allowed_extensions, key=len, reverse=True):
        if name.endswith(ext):
            return True

    logger.warning(
        f"  {file_type} '{os.path.basename(filepath)}' has an unrecognised extension. "
        f"Supported: {', '.join(sorted(allowed_extensions))}"
    )
    return False


def write_fasta(records, output_path):
    """Write sequences to FASTA file.

    Args:
        records: list of (seq_id, sequence) or (seq_id, sequence, description)
        output_path: output file path
    """
    with open(output_path, "w") as f:
        for rec in records:
            if len(rec) == 3:
                seq_id, seq, desc = rec
                f.write(f">{seq_id} {desc}\n")
            else:
                seq_id, seq = rec
                f.write(f">{seq_id}\n")
            # Write sequence in 80-char lines
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def reverse_complement(seq):
    """Return reverse complement of a DNA sequence."""
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def run_command(cmd, description="", check=True, capture=True):
    """Run a shell command with logging.

    Args:
        cmd: command as list of strings or single string
        description: human-readable description for logging
        check: raise on non-zero exit
        capture: capture stdout/stderr

    Returns:
        subprocess.CompletedProcess
    """
    if isinstance(cmd, str):
        cmd_str = cmd
        shell = True
    else:
        cmd_str = " ".join(str(c) for c in cmd)
        shell = False

    if description:
        logger.info(f"{description}")
    logger.debug(f"Running: {cmd_str}")

    result = subprocess.run(
        cmd,
        shell=shell,
        capture_output=capture,
        text=True if capture else None,
    )

    if check and result.returncode != 0:
        stderr = result.stderr if capture else ""
        logger.error(f"Command failed (exit {result.returncode}): {cmd_str}")
        if stderr:
            logger.error(f"stderr: {stderr[:500]}")
        raise RuntimeError(f"Command failed: {cmd_str}")

    return result


def parallel_map(func, items, threads=1, use_processes=False, desc=""):
    """Apply func to items in parallel.

    Args:
        func: callable taking one argument
        items: iterable of arguments
        threads: number of parallel workers
        use_processes: use ProcessPoolExecutor instead of ThreadPoolExecutor
        desc: description for logging

    Returns:
        list of results in order
    """
    items = list(items)
    if not items:
        return []

    if desc:
        logger.info(f"{desc} ({len(items)} items, {threads} threads)")

    if threads <= 1:
        return [func(item) for item in items]

    Executor = ProcessPoolExecutor if use_processes else ThreadPoolExecutor
    results = [None] * len(items)

    with Executor(max_workers=threads) as executor:
        future_to_idx = {
            executor.submit(func, item): i for i, item in enumerate(items)
        }
        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            try:
                results[idx] = future.result()
            except Exception as e:
                logger.error(f"Error processing item {idx}: {e}")
                results[idx] = None

    return results


def ensure_dir(path):
    """Create directory if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)


def parse_metadata(metadata_path):
    """Parse module 3 metadata file.

    Format: Name[TAB]Start[TAB]End[TAB]Sequence(optional)
    If Name is empty, assign J1, J2, etc.

    **Coordinate convention**: coordinates in the metadata file are
    1-based inclusive (matching paper/NCBI/GenBank/IGV convention).
    They are converted internally to 0-based start / 1-based-inclusive
    end (equivalent to Python slicing: seq[start:end]).

    Example: to represent a duplication of bases 1000 through 2000 in
    paper coordinates, write:
        M1_F1    1000    2000
    Internally this becomes start=999, end=2000 so that seq[999:2000]
    gives the duplicated region.

    Returns:
        list of dicts with keys: name, start (0-based), end (1-based
        inclusive = Python exclusive), sequence (or None),
        user_start_1based, user_end_1based (original values for display)
    """
    entries = []
    counter = 0

    with open(metadata_path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                logger.warning(
                    f"Metadata line {line_num}: expected at least 3 columns, got {len(parts)}. Skipping."
                )
                continue

            counter += 1
            name = parts[0].strip() if parts[0].strip() else f"J{counter}"

            try:
                user_start = int(parts[1])  # 1-based inclusive
                user_end = int(parts[2])    # 1-based inclusive
            except ValueError:
                logger.warning(
                    f"Metadata line {line_num}: invalid start/end coordinates. Skipping."
                )
                continue

            if user_start < 1:
                logger.warning(
                    f"Metadata line {line_num}: start={user_start} is less than 1 "
                    f"(coordinates are 1-based). Skipping."
                )
                continue

            # Convert to internal 0-based start / Python exclusive end
            internal_start = user_start - 1
            internal_end = user_end

            sequence = parts[3].strip().upper() if len(parts) >= 4 and parts[3].strip() else None

            entries.append({
                "name": name,
                "start": internal_start,
                "end": internal_end,
                "user_start_1based": user_start,
                "user_end_1based": user_end,
                "sequence": sequence,
            })

    logger.info(f"Parsed {len(entries)} entries from metadata file")
    return entries
