#!/usr/bin/env python3
"""
DNase preprocessing pipeline following Alphagenome methodology.

This script processes DNase-seq BAM files to produce normalized BigWig tracks
for machine learning training, following the preprocessing approach described
in the Alphagenome paper from DeepMind.

Pipeline steps:
1. Parse metadata to identify replicates per cell line
2. Validate input files (BAMs, indices, chromosome compatibility)
3. Convert BAM -> BigWig using ChromBPNet's reads_to_bigwig.py
4. Aggregate replicates by averaging
5. Normalize to 10^8 total counts
6. Validate output and write manifest

Author: Take-Home Assignment
"""

import argparse
import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam
import pyBigWig


# ============================================================================
# Configuration
# ============================================================================

TARGET_TOTAL = 1e8  # 100 million counts
TOLERANCE = 0.01    # 1% tolerance for normalization validation
CHROMBPNET_REPO = "https://github.com/kundajelab/chrombpnet"


# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class ReplicateInfo:
    """Information about a single replicate."""
    replicate_id: str
    bam_path: str
    bigwig_path: Optional[str] = None


@dataclass
class CellLineData:
    """Data for a single cell line."""
    cell_line: str
    replicates: List[ReplicateInfo] = field(default_factory=list)
    output_dir: Path = None
    final_bigwig: Path = None


@dataclass
class PipelineResult:
    """Result of processing a cell line."""
    cell_line: str
    success: bool
    error_message: Optional[str] = None
    manifest: Optional[Dict] = None


# ============================================================================
# Logging Setup
# ============================================================================

def setup_logging(log_file: Path, verbose: bool = False) -> logging.Logger:
    """Set up logging to both file and console."""
    logger = logging.getLogger(log_file.stem)
    logger.setLevel(logging.DEBUG)

    # Clear existing handlers
    logger.handlers = []

    # File handler - always DEBUG level
    fh = logging.FileHandler(log_file, mode='w')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    ))
    logger.addHandler(fh)

    # Console handler
    ch = logging.StreamHandler(sys.stderr)
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
    logger.addHandler(ch)

    return logger


# ============================================================================
# Validation Functions
# ============================================================================

def validate_bam_file(bam_path: str, logger: logging.Logger) -> Tuple[bool, str]:
    """Validate that BAM file exists, is readable, and has an index."""
    bam_path = Path(bam_path)

    # Check BAM exists
    if not bam_path.exists():
        return False, f"BAM file does not exist: {bam_path}"

    # Check BAM is readable
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            # Try to read header
            _ = bam.header
    except Exception as e:
        return False, f"BAM file is not readable: {bam_path} - {e}"

    # Check index exists
    bai_path = Path(str(bam_path) + ".bai")
    alt_bai_path = bam_path.with_suffix(".bai")

    if not bai_path.exists() and not alt_bai_path.exists():
        return False, f"BAM index not found for: {bam_path}"

    return True, ""


def validate_chromosome_compatibility(
    bam_path: str,
    chrom_sizes: Dict[str, int],
    logger: logging.Logger
) -> Tuple[bool, str]:
    """Check that BAM chromosomes are compatible with reference."""
    try:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            bam_chroms = set(bam.references)
            ref_chroms = set(chrom_sizes.keys())

            # Check for at least some overlap
            common = bam_chroms & ref_chroms
            if len(common) == 0:
                return False, (
                    f"No common chromosomes between BAM and reference. "
                    f"BAM has: {sorted(list(bam_chroms))[:5]}..., "
                    f"Reference has: {sorted(list(ref_chroms))[:5]}..."
                )

            logger.debug(f"Found {len(common)} common chromosomes")
            return True, ""

    except Exception as e:
        return False, f"Error checking chromosome compatibility: {e}"


def validate_bigwig(bw_path: Path, logger: logging.Logger) -> Tuple[bool, str]:
    """Validate that BigWig file exists and is readable."""
    if not bw_path.exists():
        return False, f"BigWig file does not exist: {bw_path}"

    try:
        with pyBigWig.open(str(bw_path)) as bw:
            if bw is None:
                return False, f"BigWig file is not valid: {bw_path}"
            chroms = bw.chroms()
            if not chroms:
                return False, f"BigWig has no chromosomes: {bw_path}"
        return True, ""
    except Exception as e:
        return False, f"Error reading BigWig: {bw_path} - {e}"


def check_output_conflicts(output_dir: Path, logger: logging.Logger) -> Tuple[bool, str]:
    """Check for output directory conflicts."""
    if output_dir.exists():
        # Check if it's a file (not directory)
        if output_dir.is_file():
            return False, f"Output path exists as a file: {output_dir}"

        # Check for existing output files
        final_bw = output_dir / "dnase_avg_norm100M.bw"
        if final_bw.exists():
            logger.warning(f"Output file already exists and will be overwritten: {final_bw}")

    return True, ""


# ============================================================================
# File I/O Functions
# ============================================================================

def load_chrom_sizes(chrom_sizes_path: str) -> Dict[str, int]:
    """Load chromosome sizes from file."""
    chrom_sizes = {}
    with open(chrom_sizes_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def parse_metadata(
    metadata_path: str,
    target_cell_lines: Optional[List[str]] = None
) -> Dict[str, CellLineData]:
    """
    Parse metadata TSV to extract cell line and replicate information.

    Expected columns: File accession, Biosample term name, and path info
    """
    df = pd.read_csv(metadata_path, sep='\t')

    # Identify relevant columns
    # Look for cell line column (various possible names)
    cell_line_col = None
    for col in ['Biosample term name', 'cell_line', 'Cell line', 'biosample']:
        if col in df.columns:
            cell_line_col = col
            break

    if cell_line_col is None:
        raise ValueError(f"Could not find cell line column in metadata. Columns: {df.columns.tolist()}")

    # Look for file/replicate ID column
    id_col = None
    for col in ['File accession', 'replicate_id', 'sample_id', 'accession']:
        if col in df.columns:
            id_col = col
            break

    if id_col is None:
        raise ValueError(f"Could not find ID column in metadata. Columns: {df.columns.tolist()}")

    # Look for BAM path column
    bam_col = None
    for col in ['bam_path', 'File path', 'path', 'bam']:
        if col in df.columns:
            bam_col = col
            break

    if bam_col is None:
        raise ValueError(f"Could not find BAM path column in metadata. Columns: {df.columns.tolist()}")

    # Filter for DNase if assay column exists
    assay_col = None
    for col in ['Assay', 'assay', 'Assay type']:
        if col in df.columns:
            assay_col = col
            break

    if assay_col:
        # Filter for DNase-seq
        df = df[df[assay_col].str.contains('DNase', case=False, na=False)]

    # Group by cell line
    cell_lines_data = {}

    for cell_line in df[cell_line_col].unique():
        if target_cell_lines and cell_line not in target_cell_lines:
            continue

        cell_df = df[df[cell_line_col] == cell_line]

        replicates = []
        for _, row in cell_df.iterrows():
            rep = ReplicateInfo(
                replicate_id=str(row[id_col]),
                bam_path=str(row[bam_col])
            )
            replicates.append(rep)

        if replicates:
            cell_lines_data[cell_line] = CellLineData(
                cell_line=cell_line,
                replicates=replicates
            )

    return cell_lines_data


# ============================================================================
# BigWig Processing Functions
# ============================================================================

def get_chrombpnet_info() -> Tuple[str, str]:
    """Get ChromBPNet repository URL and commit SHA."""
    # Try to get commit from installed package
    try:
        result = subprocess.run(
            ['pip', 'show', 'chrombpnet'],
            capture_output=True, text=True
        )
        # Parse version from output
        for line in result.stdout.split('\n'):
            if line.startswith('Version:'):
                version = line.split(':')[1].strip()
                break
        else:
            version = "unknown"
    except:
        version = "unknown"

    # Try to get commit SHA if installed from git
    try:
        import chrombpnet
        pkg_path = Path(chrombpnet.__file__).parent
        git_dir = pkg_path.parent / '.git'
        if git_dir.exists():
            result = subprocess.run(
                ['git', '-C', str(pkg_path.parent), 'rev-parse', 'HEAD'],
                capture_output=True, text=True
            )
            commit_sha = result.stdout.strip()
        else:
            commit_sha = version
    except:
        commit_sha = "unknown"

    return CHROMBPNET_REPO, commit_sha


def run_reads_to_bigwig(
    bam_path: str,
    output_prefix: str,
    chrom_sizes_path: str,
    genome_fasta: Optional[str],
    threads: int,
    logger: logging.Logger
) -> Tuple[bool, str, str]:
    """
    Run ChromBPNet's reads_to_bigwig.py to convert BAM to BigWig.

    Returns: (success, output_path, command)
    """
    # Build command
    cmd = [
        'python', '-m', 'chrombpnet.helpers.preprocessing.reads_to_bigwig',
        '-ibam', bam_path,
        '-c', chrom_sizes_path,
        '-op', output_prefix,
        '-d', 'DNASE'
    ]

    if genome_fasta:
        cmd.extend(['-g', genome_fasta])

    cmd_str = ' '.join(cmd)
    logger.info(f"Running: {cmd_str}")

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600  # 1 hour timeout
        )

        if result.returncode != 0:
            logger.error(f"reads_to_bigwig failed: {result.stderr}")
            return False, "", cmd_str

        # Output file is {prefix}_unstranded.bw
        output_bw = f"{output_prefix}_unstranded.bw"

        if not Path(output_bw).exists():
            logger.error(f"Expected output not found: {output_bw}")
            return False, "", cmd_str

        return True, output_bw, cmd_str

    except subprocess.TimeoutExpired:
        logger.error("reads_to_bigwig timed out")
        return False, "", cmd_str
    except Exception as e:
        logger.error(f"Error running reads_to_bigwig: {e}")
        return False, "", cmd_str


def compute_bigwig_total(bw_path: str) -> float:
    """Compute total signal sum from BigWig file."""
    total = 0.0
    with pyBigWig.open(bw_path) as bw:
        for chrom, length in bw.chroms().items():
            # Get all intervals
            intervals = bw.intervals(chrom)
            if intervals:
                for start, end, value in intervals:
                    total += value * (end - start)
    return total


def aggregate_bigwigs(
    bigwig_paths: List[str],
    output_path: str,
    chrom_sizes: Dict[str, int],
    method: str = "mean",
    logger: logging.Logger = None
) -> bool:
    """
    Aggregate multiple BigWig files by averaging.

    This creates a deterministic output by processing chromosomes in sorted order.
    """
    if not bigwig_paths:
        return False

    if len(bigwig_paths) == 1:
        # Just copy the single file
        shutil.copy(bigwig_paths[0], output_path)
        return True

    # Open all BigWigs
    bws = [pyBigWig.open(p) for p in bigwig_paths]

    try:
        # Create output BigWig
        out_bw = pyBigWig.open(output_path, "w")

        # Add header with sorted chromosomes for determinism
        sorted_chroms = sorted(chrom_sizes.keys())
        header = [(c, chrom_sizes[c]) for c in sorted_chroms if c in bws[0].chroms()]
        out_bw.addHeader(header)

        # Process each chromosome
        for chrom, length in header:
            if logger:
                logger.debug(f"Aggregating chromosome {chrom}")

            # Collect all intervals from all files
            all_positions = set()
            for bw in bws:
                intervals = bw.intervals(chrom)
                if intervals:
                    for start, end, _ in intervals:
                        all_positions.add(start)
                        all_positions.add(end)

            if not all_positions:
                continue

            # Sort positions for deterministic output
            positions = sorted(all_positions)

            # For each interval, compute average
            chroms_out = []
            starts_out = []
            ends_out = []
            values_out = []

            for i in range(len(positions) - 1):
                start = positions[i]
                end = positions[i + 1]

                # Get values from all BigWigs
                values = []
                for bw in bws:
                    val = bw.stats(chrom, start, end, type="mean")
                    if val and val[0] is not None:
                        values.append(val[0])
                    else:
                        values.append(0.0)

                if method == "mean":
                    avg_val = np.mean(values)
                else:
                    avg_val = np.sum(values)

                if avg_val > 0:
                    chroms_out.append(chrom)
                    starts_out.append(start)
                    ends_out.append(end)
                    values_out.append(float(avg_val))

            # Write in batch for efficiency
            if chroms_out:
                out_bw.addEntries(chroms_out, starts_out, ends=ends_out, values=values_out)

        out_bw.close()
        return True

    finally:
        for bw in bws:
            bw.close()


def normalize_bigwig(
    input_path: str,
    output_path: str,
    target_total: float,
    chrom_sizes: Dict[str, int],
    logger: logging.Logger
) -> Tuple[float, float, float]:
    """
    Normalize BigWig to target total signal.

    Returns: (pre_total, post_total, scaling_factor)
    """
    # Compute current total
    pre_total = compute_bigwig_total(input_path)

    if pre_total == 0:
        raise ValueError("Input BigWig has zero total signal")

    scaling_factor = target_total / pre_total
    logger.info(f"Pre-normalization total: {pre_total:.2e}")
    logger.info(f"Scaling factor: {scaling_factor:.6f}")

    # Apply scaling
    with pyBigWig.open(input_path) as in_bw:
        out_bw = pyBigWig.open(output_path, "w")

        # Add header with sorted chromosomes for determinism
        sorted_chroms = sorted(chrom_sizes.keys())
        header = [(c, chrom_sizes[c]) for c in sorted_chroms if c in in_bw.chroms()]
        out_bw.addHeader(header)

        for chrom, length in header:
            intervals = in_bw.intervals(chrom)
            if intervals:
                chroms_out = []
                starts_out = []
                ends_out = []
                values_out = []

                for start, end, value in intervals:
                    chroms_out.append(chrom)
                    starts_out.append(start)
                    ends_out.append(end)
                    values_out.append(float(value * scaling_factor))

                if chroms_out:
                    out_bw.addEntries(chroms_out, starts_out, ends=ends_out, values=values_out)

        out_bw.close()

    # Verify
    post_total = compute_bigwig_total(output_path)
    logger.info(f"Post-normalization total: {post_total:.2e}")

    return pre_total, post_total, scaling_factor


# ============================================================================
# Main Pipeline Functions
# ============================================================================

def process_cell_line(
    cell_line_data: CellLineData,
    chrom_sizes_path: str,
    chrom_sizes: Dict[str, int],
    genome_fasta: Optional[str],
    threads: int,
    logger: logging.Logger
) -> PipelineResult:
    """Process a single cell line through the complete pipeline."""

    cell_line = cell_line_data.cell_line
    output_dir = cell_line_data.output_dir
    commands = []

    logger.info(f"=" * 60)
    logger.info(f"Processing cell line: {cell_line}")
    logger.info(f"Number of replicates: {len(cell_line_data.replicates)}")
    logger.info(f"=" * 60)

    # Create temp directory for intermediate files
    temp_dir = output_dir / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    try:
        # Step 1: Validate all inputs
        logger.info("Step 1: Validating input files...")
        for rep in cell_line_data.replicates:
            # Validate BAM
            valid, msg = validate_bam_file(rep.bam_path, logger)
            if not valid:
                return PipelineResult(cell_line, False, msg)

            # Validate chromosome compatibility
            valid, msg = validate_chromosome_compatibility(rep.bam_path, chrom_sizes, logger)
            if not valid:
                return PipelineResult(cell_line, False, msg)

            logger.info(f"  Validated: {rep.replicate_id}")

        # Step 2: Convert BAMs to BigWigs
        logger.info("Step 2: Converting BAMs to BigWigs...")
        bigwig_paths = []

        for rep in cell_line_data.replicates:
            output_prefix = str(temp_dir / f"{rep.replicate_id}")

            success, bw_path, cmd = run_reads_to_bigwig(
                rep.bam_path,
                output_prefix,
                chrom_sizes_path,
                genome_fasta,
                threads,
                logger
            )
            commands.append(cmd)

            if not success:
                return PipelineResult(
                    cell_line, False,
                    f"Failed to convert BAM to BigWig for {rep.replicate_id}"
                )

            # Validate output
            valid, msg = validate_bigwig(Path(bw_path), logger)
            if not valid:
                return PipelineResult(cell_line, False, msg)

            rep.bigwig_path = bw_path
            bigwig_paths.append(bw_path)
            logger.info(f"  Created: {bw_path}")

        # Step 3: Aggregate replicates
        logger.info("Step 3: Aggregating replicates...")
        aggregated_path = str(temp_dir / "dnase_avg.bw")

        success = aggregate_bigwigs(
            bigwig_paths,
            aggregated_path,
            chrom_sizes,
            method="mean",
            logger=logger
        )

        if not success:
            return PipelineResult(cell_line, False, "Failed to aggregate BigWigs")

        logger.info(f"  Created aggregated BigWig: {aggregated_path}")

        # Step 4: Normalize to target total
        logger.info("Step 4: Normalizing to 10^8 total counts...")
        final_path = str(output_dir / "dnase_avg_norm100M.bw")

        pre_total, post_total, scaling_factor = normalize_bigwig(
            aggregated_path,
            final_path,
            TARGET_TOTAL,
            chrom_sizes,
            logger
        )

        # Step 5: Validate final output
        logger.info("Step 5: Validating final output...")
        valid, msg = validate_bigwig(Path(final_path), logger)
        if not valid:
            return PipelineResult(cell_line, False, msg)

        # Check normalization tolerance
        relative_error = abs(post_total - TARGET_TOTAL) / TARGET_TOTAL
        if relative_error > TOLERANCE:
            return PipelineResult(
                cell_line, False,
                f"Normalization failed: relative error {relative_error:.4f} > {TOLERANCE}"
            )

        logger.info(f"  Final BigWig validated: {final_path}")
        logger.info(f"  Normalization error: {relative_error:.6f}")

        # Clean up temp directory
        shutil.rmtree(temp_dir)

        # Get ChromBPNet info
        repo_url, commit_sha = get_chrombpnet_info()

        # Build manifest
        manifest = {
            "cell_line": cell_line,
            "replicates": [
                {"replicate_id": rep.replicate_id, "bam_path": rep.bam_path}
                for rep in cell_line_data.replicates
            ],
            "chrombpnet": {
                "repo_url": repo_url,
                "commit_sha": commit_sha
            },
            "commands": commands,
            "normalization": {
                "target_total": int(TARGET_TOTAL),
                "pre_total": float(pre_total),
                "post_total": float(post_total),
                "scaling_factor": float(scaling_factor),
                "tolerance": float(TOLERANCE)
            },
            "aggregation": {
                "method": "mean",
                "n_replicates": len(cell_line_data.replicates)
            }
        }

        logger.info("Pipeline completed successfully!")
        return PipelineResult(cell_line, True, manifest=manifest)

    except Exception as e:
        logger.exception(f"Pipeline failed with exception: {e}")
        return PipelineResult(cell_line, False, str(e))


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="DNase preprocessing pipeline following Alphagenome methodology"
    )
    parser.add_argument(
        '--metadata', required=True,
        help='Path to metadata TSV file'
    )
    parser.add_argument(
        '--chrom-sizes', required=True,
        help='Path to chromosome sizes file'
    )
    parser.add_argument(
        '--outdir', required=True,
        help='Output directory'
    )
    parser.add_argument(
        '--threads', type=int, default=4,
        help='Number of threads (default: 4)'
    )
    parser.add_argument(
        '--genome', default=None,
        help='Reference genome FASTA (optional, for shift estimation)'
    )
    parser.add_argument(
        '--cell-lines', nargs='+', default=None,
        help='Specific cell lines to process (default: all)'
    )
    parser.add_argument(
        '--verbose', action='store_true',
        help='Enable verbose output'
    )

    args = parser.parse_args()

    # Validate inputs exist
    if not Path(args.metadata).exists():
        print(f"ERROR: Metadata file not found: {args.metadata}", file=sys.stderr)
        sys.exit(1)

    if not Path(args.chrom_sizes).exists():
        print(f"ERROR: Chromosome sizes file not found: {args.chrom_sizes}", file=sys.stderr)
        sys.exit(1)

    # Load chromosome sizes
    try:
        chrom_sizes = load_chrom_sizes(args.chrom_sizes)
    except Exception as e:
        print(f"ERROR: Failed to load chromosome sizes: {e}", file=sys.stderr)
        sys.exit(1)

    # Parse metadata
    try:
        cell_lines_data = parse_metadata(args.metadata, args.cell_lines)
    except Exception as e:
        print(f"ERROR: Failed to parse metadata: {e}", file=sys.stderr)
        sys.exit(1)

    if not cell_lines_data:
        print("ERROR: No cell lines found in metadata", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(cell_lines_data)} cell lines to process")

    # Create output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Process each cell line
    all_success = True

    for cell_line, data in cell_lines_data.items():
        # Set up output directory for this cell line
        cell_outdir = outdir / cell_line
        cell_outdir.mkdir(parents=True, exist_ok=True)
        data.output_dir = cell_outdir

        # Check for output conflicts
        log_path = cell_outdir / "pipeline.log"
        logger = setup_logging(log_path, args.verbose)

        valid, msg = check_output_conflicts(cell_outdir, logger)
        if not valid:
            logger.error(msg)
            all_success = False
            continue

        # Process cell line
        result = process_cell_line(
            data,
            args.chrom_sizes,
            chrom_sizes,
            args.genome,
            args.threads,
            logger
        )

        if result.success:
            # Write manifest
            manifest_path = cell_outdir / "manifest.json"
            with open(manifest_path, 'w') as f:
                json.dump(result.manifest, f, indent=2)
            logger.info(f"Manifest written to: {manifest_path}")
        else:
            logger.error(f"Pipeline failed: {result.error_message}")
            all_success = False

    # Exit with appropriate code
    if all_success:
        print("All cell lines processed successfully!")
        sys.exit(0)
    else:
        print("Some cell lines failed. Check logs for details.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
