#!/usr/bin/env python3
"""
DNase Preprocessing Pipeline - Take-Home Assignment

Implement the DNase preprocessing pipeline following the Alphagenome methodology.

Your task is to implement this script to:
1. Download DNase-seq BAM files from ENCODE for the specified cell lines
2. Convert BAM -> BigWig using ChromBPNet's reads_to_bigwig.py
3. Aggregate replicates by averaging per cell line
4. Normalize to 10^8 total counts
5. Output a manifest.json with provenance information

Target cell lines and their ENCODE experiment accessions:
- MCF-7: ENCSR000EJD
- K562: ENCSR000EMT
- SK-N-SH: ENCSR000ENH

You can find BAM file URLs at:
https://www.encodeproject.org/experiments/{ACCESSION}/

Look for "alignments" files in BAM format with hg38 assembly.

Expected output structure:
    out/
    ├── MCF-7/
    │   ├── dnase_avg_norm100M.bw
    │   ├── manifest.json
    │   └── pipeline.log
    ├── K562/
    │   └── ...
    └── SK-N-SH/
        └── ...

The manifest.json should contain:
- cell_line: Name of the cell line
- replicates: List of {replicate_id, bam_path} for each replicate
- chrombpnet: {repo_url, commit_sha} for reproducibility
- commands: List of commands executed
- normalization: {target_total, pre_total, post_total, scaling_factor, tolerance}
- aggregation: {method: "mean", n_replicates: N}

Requirements:
- Deterministic output (same input -> same output)
- Handle missing/corrupted files gracefully
- Validate chromosome compatibility
- Log all steps to pipeline.log

Author: [Your Name]
"""

import argparse
import sys
from pathlib import Path


def main():
    """Main entry point - implement your pipeline here."""
    parser = argparse.ArgumentParser(
        description="DNase preprocessing pipeline following Alphagenome methodology"
    )
    # TODO: Add your command line arguments here
    # Suggested arguments:
    #   --outdir: Output directory (required)
    #   --chrom-sizes: Path to chromosome sizes file (required)
    #   --threads: Number of threads (default: 4)
    #   --cell-lines: Specific cell lines to process (default: all)
    #   --verbose: Enable verbose output

    parser.add_argument(
        '--outdir', required=True,
        help='Output directory'
    )
    parser.add_argument(
        '--chrom-sizes', required=True,
        help='Path to chromosome sizes file'
    )
    parser.add_argument(
        '--threads', type=int, default=4,
        help='Number of threads (default: 4)'
    )
    parser.add_argument(
        '--cell-lines', nargs='+', default=None,
        help='Specific cell lines to process (default: all)'
    )
    parser.add_argument(
        '--verbose', action='store_true',
        help='Enable verbose output'
    )

    # For testing compatibility - accepts metadata but you should
    # implement direct ENCODE download for the real submission
    parser.add_argument(
        '--metadata', default=None,
        help='Path to metadata TSV (for testing only)'
    )

    args = parser.parse_args()

    # TODO: Implement the pipeline
    #
    # Suggested steps:
    # 1. Define ENCODE URLs for each cell line's BAM files
    # 2. Download BAM files (or use local paths if --metadata provided for testing)
    # 3. Validate BAM files (existence, index, chromosome compatibility)
    # 4. Convert BAM -> BigWig using ChromBPNet
    # 5. Aggregate replicates by cell line (mean)
    # 6. Normalize to 10^8 total counts
    # 7. Write manifest.json with provenance
    # 8. Clean up intermediate files

    print("TODO: Implement the DNase preprocessing pipeline", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
