#!/usr/bin/env python
"""
Validate final BigWigs and generate manifest.json files.
"""

import os
import json
import numpy as np
import pyBigWig

CHROM_SIZES_FILE = "reference/hg38.chrom.sizes"
EXPECTED_SUM = 1e8
TOLERANCE = 0.001  # 0.1% tolerance

# UPDATE THIS with your actual commit SHA
CHROMBPNET_COMMIT = "YOUR_COMMIT_SHA_HERE"

def load_chrom_sizes():
    chrom_sizes = {}
    with open(CHROM_SIZES_FILE) as f:
        for line in f:
            chrom, size = line.strip().split("\t")
            chrom_sizes[chrom] = int(size)
    return chrom_sizes

def validate_bigwig(bw_path, chrom_sizes):
    """
    Validate BigWig file.
    
    Checks:
    1. File exists and is readable
    2. Chromosome names are valid
    3. Total signal sum is ~10^8
    """
    # Check 1: Readable
    if not os.path.exists(bw_path):
        raise ValueError(f"File does not exist: {bw_path}")
    
    try:
        bw = pyBigWig.open(bw_path)
    except Exception as e:
        raise ValueError(f"Cannot open BigWig: {e}")
    
    if bw is None:
        raise ValueError(f"BigWig is None: {bw_path}")
    
    # Check 2: Chromosome names
    bw_chroms = set(bw.chroms().keys())
    expected_chroms = set(chrom_sizes.keys())
    unexpected = bw_chroms - expected_chroms
    if unexpected:
        bw.close()
        raise ValueError(f"Unexpected chromosomes: {unexpected}")
    
    # Check 3: Total signal sum
    total = 0.0
    for chrom, size in bw.chroms().items():
        vals = bw.values(chrom, 0, size)
        vals = np.nan_to_num(vals, nan=0.0)
        total += np.sum(vals)
    
    bw.close()
    
    relative_error = abs(total - EXPECTED_SUM) / EXPECTED_SUM
    if relative_error > TOLERANCE:
        raise ValueError(
            f"Total signal {total:.6e} differs from expected {EXPECTED_SUM:.0e} "
            f"by {relative_error*100:.4f}% (tolerance: {TOLERANCE*100}%)"
        )
    
    return total

def write_manifest(cell_line, file_accessions, normalization, aggregation, output_path):
    """Write manifest.json."""
    manifest = {
        "cell_line": cell_line,
        "replicates": file_accessions,
        "chrombpnet": {
            "repo_url": "https://github.com/kundajelab/chrombpnet",
            "commit_sha": CHROMBPNET_COMMIT
        },
        "normalization": {
            "pre_total": normalization["pre_total"],
            "post_total": normalization["post_total"],
            "scaling_factor": normalization["scaling_factor"]
        },
        "aggregation": {
            "non_zero_average_mean": aggregation["non_zero_average_mean"],
            "n_replicates": aggregation["n_replicates"]
        }
    }
    
    with open(output_path, "w") as f:
        json.dump(manifest, f, indent=2)
    
    print(f"Wrote {output_path}")

def main():
    chrom_sizes = load_chrom_sizes()
    
    # Load processing results
    with open("processing_results.json") as f:
        results = json.load(f)
    
    for cell_line in ["GM12878", "HeLa-S3", "SK-N-SH"]:
        print(f"\n{'='*50}")
        print(f"Validating {cell_line}")
        print(f"{'='*50}")
        
        bw_path = f"out/{cell_line}/dnase_avg_norm100M.bw"
        
        # Validate
        try:
            total = validate_bigwig(bw_path, chrom_sizes)
            print(f"  ✓ BigWig valid (total={total:.6e})")
        except ValueError as e:
            print(f"  ✗ Validation failed: {e}")
            continue
        
        # Write manifest
        cell_results = results[cell_line]
        manifest_path = f"out/{cell_line}/manifest.json"
        write_manifest(
            cell_line,
            cell_results["file_accessions"],
            cell_results["normalization"],
            cell_results["aggregation"],
            manifest_path
        )
    
    print("\nDone!")

if __name__ == "__main__":
    main()