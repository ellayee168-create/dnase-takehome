#!/usr/bin/env python3
"""
Generate synthetic test data for the DNase preprocessing pipeline.

This script creates:
- Reference chromosome sizes file
- Synthetic BAM files with controlled read distributions
- Metadata TSV file

The synthetic data allows for deterministic testing with known expected outputs.
"""

import argparse
import os
import random
import subprocess
import sys
from pathlib import Path

import pysam
import pyBigWig
import numpy as np


# Seed for reproducibility
SEED = 42


def create_chrom_sizes(output_path: str, chromosomes: dict):
    """Create a chromosome sizes file."""
    with open(output_path, 'w') as f:
        for chrom, size in sorted(chromosomes.items()):
            f.write(f"{chrom}\t{size}\n")
    print(f"Created: {output_path}")


def create_synthetic_bam(
    output_path: str,
    chrom_sizes: dict,
    n_reads: int,
    read_length: int = 50,
    seed: int = SEED
):
    """
    Create a synthetic BAM file with random reads.

    Reads are distributed across chromosomes proportional to chromosome size.
    """
    random.seed(seed)
    np.random.seed(seed)

    # Create header
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'LN': size, 'SN': chrom} for chrom, size in sorted(chrom_sizes.items())]
    }

    # Calculate reads per chromosome (proportional to size)
    total_size = sum(chrom_sizes.values())
    reads_per_chrom = {
        chrom: max(1, int(n_reads * size / total_size))
        for chrom, size in chrom_sizes.items()
    }

    # Create BAM file
    with pysam.AlignmentFile(output_path, "wb", header=header) as outf:
        read_id = 0
        for chrom, n_chrom_reads in sorted(reads_per_chrom.items()):
            chrom_size = chrom_sizes[chrom]

            for _ in range(n_chrom_reads):
                # Random position (ensuring read fits)
                max_pos = max(0, chrom_size - read_length - 1)
                pos = random.randint(0, max_pos)

                # Create alignment
                a = pysam.AlignedSegment()
                a.query_name = f"read_{read_id}"
                a.query_sequence = "A" * read_length
                a.flag = 0 if random.random() > 0.5 else 16  # Random strand
                a.reference_id = outf.get_tid(chrom)
                a.reference_start = pos
                a.mapping_quality = 60
                a.cigar = [(0, read_length)]  # M operation
                a.query_qualities = pysam.qualitystring_to_array("I" * read_length)

                outf.write(a)
                read_id += 1

    # Sort and index
    sorted_path = output_path.replace('.bam', '.sorted.bam')
    pysam.sort("-o", sorted_path, output_path)
    os.rename(sorted_path, output_path)
    pysam.index(output_path)

    print(f"Created: {output_path} ({read_id} reads)")
    print(f"Created: {output_path}.bai")


def create_metadata_tsv(
    output_path: str,
    cell_lines: dict,
    bam_dir: str
):
    """
    Create metadata TSV file.

    Args:
        cell_lines: Dict mapping cell_line name to list of replicate IDs
        bam_dir: Directory containing BAM files
    """
    with open(output_path, 'w') as f:
        # Header
        f.write("File accession\tBiosample term name\tAssay\tbam_path\n")

        for cell_line, replicates in cell_lines.items():
            for rep_id in replicates:
                bam_path = os.path.join(bam_dir, f"{rep_id}.bam")
                f.write(f"{rep_id}\t{cell_line}\tDNase-seq\t{bam_path}\n")

    print(f"Created: {output_path}")


def create_expected_bigwig(
    output_path: str,
    chrom_sizes: dict,
    total_signal: float,
    seed: int = SEED
):
    """Create an expected output BigWig with known total signal."""
    np.random.seed(seed)

    with pyBigWig.open(output_path, "w") as bw:
        # Add header
        header = [(chrom, size) for chrom, size in sorted(chrom_sizes.items())]
        bw.addHeader(header)

        # Distribute signal across chromosomes
        total_size = sum(chrom_sizes.values())
        remaining_signal = total_signal

        for chrom, size in sorted(chrom_sizes.items()):
            # Create sparse signal (random intervals)
            n_intervals = min(100, size // 1000)
            if n_intervals == 0:
                continue

            # Generate random intervals
            starts = sorted(np.random.randint(0, size - 100, n_intervals))
            ends = [min(s + np.random.randint(50, 150), size) for s in starts]

            # Assign signal proportional to chromosome size
            chrom_signal = total_signal * (size / total_size)
            values = np.random.exponential(chrom_signal / n_intervals, n_intervals)

            # Normalize to exact signal
            values = values * (chrom_signal / (values.sum() * 100))  # Approx interval size

            chroms = [chrom] * len(starts)
            bw.addEntries(chroms, starts, ends=ends, values=values.tolist())

    print(f"Created expected BigWig: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate test data for DNase pipeline")
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--n-reads', type=int, default=10000, help='Reads per BAM')
    parser.add_argument('--small', action='store_true', help='Generate small test set')
    parser.add_argument('--seed', type=int, default=SEED, help='Random seed')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define chromosome sizes (small for testing)
    if args.small:
        chrom_sizes = {
            "chr1": 100000,
            "chr2": 80000,
            "chr3": 60000,
        }
        n_reads = args.n_reads
    else:
        chrom_sizes = {
            "chr1": 1000000,
            "chr2": 900000,
            "chr3": 800000,
            "chr4": 700000,
            "chr5": 600000,
        }
        n_reads = args.n_reads

    # Create chromosome sizes file
    chrom_sizes_path = output_dir / "reference.chrom.sizes"
    create_chrom_sizes(str(chrom_sizes_path), chrom_sizes)

    # Define cell lines and replicates
    cell_lines = {
        "GM12878": ["ENCFF001GM1_rep1", "ENCFF002GM1_rep2"],
        "HeLa-S3": ["ENCFF001HLS_rep1", "ENCFF002HLS_rep2", "ENCFF003HLS_rep3"],
        "SK-N-SH": ["ENCFF001SKN_rep1", "ENCFF002SKN_rep2"],
    }

    # Create BAM directory
    bam_dir = output_dir / "bams"
    bam_dir.mkdir(exist_ok=True)

    # Create BAM files
    seed_offset = 0
    for cell_line, replicates in cell_lines.items():
        for rep_id in replicates:
            bam_path = bam_dir / f"{rep_id}.bam"
            create_synthetic_bam(
                str(bam_path),
                chrom_sizes,
                n_reads,
                seed=args.seed + seed_offset
            )
            seed_offset += 1

    # Create metadata
    metadata_path = output_dir / "metadata.tsv"
    create_metadata_tsv(str(metadata_path), cell_lines, str(bam_dir))

    print("\nTest data generation complete!")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
