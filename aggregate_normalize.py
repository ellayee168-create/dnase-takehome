#!/usr/bin/env python
"""
Average replicates and normalize to 100M total counts.
PATCHED: Uses pyBigWig.addEntries() directly instead of bedGraphToBigWig
to avoid signal loss from compression/rounding.
"""

import os
import json
import numpy as np
import pyBigWig

CHROM_SIZES_FILE = "reference/hg38.chrom.sizes"
TARGET_SUM = 1e8

def load_chrom_sizes():
    """Load chromosome sizes from file."""
    chrom_sizes = {}
    with open(CHROM_SIZES_FILE) as f:
        for line in f:
            chrom, size = line.strip().split("\t")
            chrom_sizes[chrom] = int(size)
    return chrom_sizes

def get_bigwig_values(bw_path, chrom, size):
    """Get values for a chromosome, handling None/NaN."""
    bw = pyBigWig.open(bw_path)
    if chrom not in bw.chroms():
        bw.close()
        return np.zeros(size)
    vals = bw.values(chrom, 0, size)
    bw.close()
    vals = np.array(vals, dtype=np.float64)
    vals = np.nan_to_num(vals, nan=0.0)
    return vals

def write_bedgraph(values_by_chrom, output_path, chrom_sizes):
    """Write values to bedGraph format."""
    with open(output_path, "w") as f:
        for chrom in sorted(chrom_sizes.keys()):
            if chrom not in values_by_chrom:
                continue
            vals = values_by_chrom[chrom]
            
            # Write runs of identical values
            i = 0
            while i < len(vals):
                if vals[i] > 0:
                    start = i
                    val = vals[i]
                    while i < len(vals) and vals[i] == val:
                        i += 1
                    f.write(f"{chrom}\t{start}\t{i}\t{val}\n")
                else:
                    i += 1

def bedgraph_to_bigwig(bedgraph_path, bigwig_path):
    """
    Convert bedGraph to BigWig using pyBigWig directly.
    PATCHED: No longer uses external bedGraphToBigWig tool.
    """
    # Load chromosome sizes
    chrom_sizes = load_chrom_sizes()
    
    # Parse bedGraph file
    intervals_by_chrom = {}
    with open(bedgraph_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                value = float(parts[3])
                
                if chrom not in intervals_by_chrom:
                    intervals_by_chrom[chrom] = []
                intervals_by_chrom[chrom].append((start, end, value))
    
    # Create BigWig file
    bw = pyBigWig.open(bigwig_path, "w")
    
    # Add header with chromosome sizes (sorted order)
    header = [(chrom, size) for chrom, size in sorted(chrom_sizes.items())]
    bw.addHeader(header)
    
    # Write intervals for each chromosome
    for chrom, size in header:
        if chrom in intervals_by_chrom and intervals_by_chrom[chrom]:
            intervals = sorted(intervals_by_chrom[chrom], key=lambda x: x[0])
            
            starts = [iv[0] for iv in intervals]
            ends = [iv[1] for iv in intervals]
            values = [iv[2] for iv in intervals]
            
            bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)
    
    bw.close()

def average_bigwigs(bw_paths, output_path, chrom_sizes):
    """
    Average multiple BigWigs into one.
    Returns aggregation stats.
    """
    n_replicates = len(bw_paths)
    all_nonzero = []
    values_by_chrom = {}
    
    for chrom in sorted(chrom_sizes.keys()):
        size = chrom_sizes[chrom]
        
        # Get values from all replicates
        all_vals = []
        for bw_path in bw_paths:
            vals = get_bigwig_values(bw_path, chrom, size)
            all_vals.append(vals)
        
        # Average
        avg_vals = np.mean(all_vals, axis=0)
        values_by_chrom[chrom] = avg_vals
        
        # Track non-zero for stats
        nonzero = avg_vals[avg_vals > 0]
        all_nonzero.extend(nonzero)
    
    # Write to bedGraph then convert to BigWig
    bedgraph_path = output_path.replace(".bw", ".bedGraph")
    write_bedgraph(values_by_chrom, bedgraph_path, chrom_sizes)
    bedgraph_to_bigwig(bedgraph_path, output_path)
    os.remove(bedgraph_path)
    
    return {
        "non_zero_average_mean": float(np.mean(all_nonzero)) if all_nonzero else 0.0,
        "n_replicates": n_replicates
    }

def normalize_bigwig(input_path, output_path, chrom_sizes, target_sum=TARGET_SUM):
    """
    Normalize BigWig so total signal = target_sum.
    Returns normalization stats.
    """
    # Pass 1: Calculate total
    pre_total = 0.0
    values_by_chrom = {}
    
    for chrom in sorted(chrom_sizes.keys()):
        size = chrom_sizes[chrom]
        vals = get_bigwig_values(input_path, chrom, size)
        pre_total += np.sum(vals)
        values_by_chrom[chrom] = vals
    
    # Calculate scaling factor
    scaling_factor = target_sum / pre_total
    
    # Pass 2: Scale and write
    for chrom in values_by_chrom:
        values_by_chrom[chrom] = values_by_chrom[chrom] * scaling_factor
    
    bedgraph_path = output_path.replace(".bw", ".bedGraph")
    write_bedgraph(values_by_chrom, bedgraph_path, chrom_sizes)
    bedgraph_to_bigwig(bedgraph_path, output_path)
    os.remove(bedgraph_path)
    
    return {
        "pre_total": float(pre_total),
        "post_total": float(target_sum),
        "scaling_factor": float(scaling_factor)
    }

def process_cell_line(cell_line, file_accessions, chrom_sizes):
    """Process one cell line: average replicates, normalize."""
    print(f"\n{'='*50}")
    print(f"Processing {cell_line}")
    print(f"{'='*50}")
    
    out_dir = f"out/{cell_line}"
    os.makedirs(out_dir, exist_ok=True)
    
    # Get BigWig paths for this cell line
    bw_paths = [f"bigwigs/{acc}_unstranded.bw" for acc in file_accessions]
    
    # Check all exist
    for bw_path in bw_paths:
        if not os.path.exists(bw_path):
            raise FileNotFoundError(f"BigWig not found: {bw_path}")
    
    print(f"Averaging {len(bw_paths)} replicates...")
    avg_path = f"{out_dir}/dnase_avg.bw"
    aggregation_stats = average_bigwigs(bw_paths, avg_path, chrom_sizes)
    print(f"  Aggregation stats: {aggregation_stats}")
    
    print(f"Normalizing to {TARGET_SUM:.0e} counts...")
    final_path = f"{out_dir}/dnase_avg_norm100M.bw"
    normalization_stats = normalize_bigwig(avg_path, final_path, chrom_sizes)
    print(f"  Normalization stats: {normalization_stats}")
    
    # Clean up intermediate
    os.remove(avg_path)
    
    return aggregation_stats, normalization_stats

def main():
    chrom_sizes = load_chrom_sizes()
    
    # Load metadata and group by cell line
    cell_line_files = {}
    with open("metadata.tsv") as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split("\t")
            cell_line = parts[0]
            file_acc = parts[2]
            
            if cell_line not in cell_line_files:
                cell_line_files[cell_line] = []
            cell_line_files[cell_line].append(file_acc)
    
    # Process each cell line
    results = {}
    for cell_line in ["GM12878", "HeLa-S3", "SK-N-SH"]:
        if cell_line not in cell_line_files:
            print(f"WARNING: No files found for {cell_line}")
            continue
        
        file_accessions = cell_line_files[cell_line]
        agg_stats, norm_stats = process_cell_line(cell_line, file_accessions, chrom_sizes)
        
        results[cell_line] = {
            "file_accessions": file_accessions,
            "aggregation": agg_stats,
            "normalization": norm_stats
        }
    
    # Save results for manifest generation
    with open("processing_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    print("\nDone! Results saved to processing_results.json")

if __name__ == "__main__":
    main()
