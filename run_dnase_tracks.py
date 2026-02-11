#!/usr/bin/env python
"""
DNase preprocessing pipeline following AlphaGenome methodology.

1. Downloads DNase-seq BAM files from ENCODE for specified cell lines
2. Converts BAMs to base-resolution BigWigs (replicating ChromBPNet methodology)
3. Averages replicates within each cell line
4. Normalizes to 100M total counts
5. Validates output and generates manifest.json

Usage:
    python run_dnase_tracks.py

    # With pre-downloaded BAMs (skip download):
    python run_dnase_tracks.py --skip-download

    # Specify cell lines:
    python run_dnase_tracks.py --cell-lines GM12878 HeLa-S3

Output:
    out/<cell_line>/dnase_avg_norm100M.bw
    out/<cell_line>/manifest.json
    out/<cell_line>/pipeline.log

References:
    - AlphaGenome paper: https://doi.org/10.1038/s41586-025-10014-0
    - ChromBPNet: https://github.com/kundajelab/chrombpnet
"""

import os
import sys
import json
import logging
import argparse
import subprocess
import tempfile
import time
from datetime import datetime

import requests
import numpy as np
import pyBigWig
import pyfaidx

# =============================================================================
# Configuration
# =============================================================================

CELL_LINES = ["GM12878", "HeLa-S3", "SK-N-SH"]
GENOME_FA = "reference/hg38.fa"
CHROM_SIZES = "reference/hg38.chrom.sizes"
TARGET_SUM = 1e8  # 100 million

# DNase shifts from AlphaGenome paper:
# "Shifts of 0/+1 for DNase-seq were applied"
PLUS_SHIFT_DELTA = 0
MINUS_SHIFT_DELTA = 1

# ChromBPNet commit used as reference for methodology
CHROMBPNET_REPO = "https://github.com/kundajelab/chrombpnet"
CHROMBPNET_COMMIT = None  # Will be set from chrombpnet/.git or argument

# =============================================================================
# Logging Setup
# =============================================================================

def setup_logging(log_path):
    """Configure logging to both file and console."""
    # Remove any existing handlers
    root = logging.getLogger()
    for handler in root.handlers[:]:
        root.removeHandler(handler)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # File handler
    file_handler = logging.FileHandler(log_path, mode='w')
    file_handler.setLevel(logging.INFO)
    file_handler.setFormatter(formatter)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)
    
    # Configure root logger
    root.setLevel(logging.INFO)
    root.addHandler(file_handler)
    root.addHandler(console_handler)

# =============================================================================
# ENCODE Query Functions
# =============================================================================

def get_json(url):
    """Get JSON from ENCODE API with rate limiting."""
    response = requests.get(url, headers={"Accept": "application/json"})
    time.sleep(0.1)
    return response.json()

def trace_read_length(file_acc, depth=0):
    """Recursively trace derived_from to find read_length from FASTQ files."""
    if depth > 5:
        return []
    
    detail = get_json(f"https://www.encodeproject.org/files/{file_acc}/?format=json")
    file_format = detail.get("file_format")
    read_length = detail.get("read_length")
    
    if file_format == "fastq" and read_length:
        return [read_length]
    
    derived_from = detail.get("derived_from", [])
    read_lengths = []
    
    for df in derived_from:
        if isinstance(df, str):
            df_acc = df.split("/")[-2]
            read_lengths.extend(trace_read_length(df_acc, depth + 1))
    
    return read_lengths

def get_spot1_score(file_detail):
    """Get spot1_score (DNase equivalent of FRiP) from HotspotQualityMetric."""
    for qm in file_detail.get("quality_metrics", []):
        qm_type = qm.get("@type", [])
        if "HotspotQualityMetric" in qm_type:
            return qm.get("spot1_score")
    return None

def has_error_audit(exp_accession):
    """Check if experiment has ERROR-level audits."""
    exp_detail = get_json(f"https://www.encodeproject.org/experiments/{exp_accession}/?format=json")
    audits = exp_detail.get("audit", {})
    return "ERROR" in audits

def query_encode_for_cell_line(cell_line):
    """
    Query ENCODE for DNase BAM files for one cell line.
    
    Applies QC filters per AlphaGenome methodology:
    - FRiP (spot1_score) > 10%
    - Read length >= 36nt
    - No ERROR audits
    """
    logging.info(f"Querying ENCODE for {cell_line} DNase-seq BAM files...")
    
    search_url = "https://www.encodeproject.org/search/"
    params = {
        "type": "File",
        "file_format": "bam",
        "output_type": "alignments",
        "assay_title": "DNase-seq",
        "biosample_ontology.term_name": cell_line,
        "assembly": "GRCh38",
        "status": "released",
        "format": "json",
        "limit": "all"
    }
    
    response = requests.get(search_url, params=params, headers={"Accept": "application/json"})
    files = response.json().get("@graph", [])
    logging.info(f"  Found {len(files)} BAM files in ENCODE")
    
    valid_files = []
    
    for f in files:
        file_acc = f.get("accession")
        file_detail = get_json(f"https://www.encodeproject.org/files/{file_acc}/?format=json")
        
        # Get experiment accession
        dataset = file_detail.get("dataset", "")
        exp_acc = dataset.split("/")[-2] if dataset else None
        
        if not exp_acc:
            continue
        
        # Filter 1: No ERROR audits
        if has_error_audit(exp_acc):
            logging.debug(f"  Skipping {file_acc}: ERROR audit")
            continue
        
        # Filter 2: Read length >= 36
        read_lengths = trace_read_length(file_acc)
        if not read_lengths:
            logging.debug(f"  Skipping {file_acc}: no read_length found")
            continue
        
        min_read_length = min(read_lengths)
        if min_read_length < 36:
            logging.debug(f"  Skipping {file_acc}: read_length={min_read_length} < 36")
            continue
        
        # Filter 3: spot1_score > 0.10
        spot1_score = get_spot1_score(file_detail)
        if spot1_score is None or spot1_score <= 0.10:
            logging.debug(f"  Skipping {file_acc}: spot1_score={spot1_score}")
            continue
        
        # Get download URL
        href = file_detail.get("href")
        if href:
            download_url = f"https://www.encodeproject.org{href}"
        else:
            download_url = f"https://www.encodeproject.org/files/{file_acc}/@@download/{file_acc}.bam"
        
        valid_files.append({
            "file_accession": file_acc,
            "experiment_accession": exp_acc,
            "bam_url": download_url,
            "spot1_score": spot1_score,
            "read_length": min_read_length
        })
        
        logging.info(f"  Valid: {file_acc} (spot1={spot1_score:.3f}, read_len={min_read_length})")
    
    logging.info(f"  {len(valid_files)} files passed QC filters")
    return valid_files

# =============================================================================
# Download Functions
# =============================================================================

def download_bam(file_info, output_dir):
    """Download a BAM file from ENCODE."""
    file_acc = file_info["file_accession"]
    bam_url = file_info["bam_url"]
    output_path = os.path.join(output_dir, f"{file_acc}.bam")
    
    if os.path.exists(output_path):
        logging.info(f"  Already downloaded: {file_acc}.bam")
        return output_path
    
    logging.info(f"  Downloading {file_acc}.bam...")
    subprocess.run(["wget", "-q", "-O", output_path, bam_url], check=True)
    return output_path

# =============================================================================
# BAM to BigWig Conversion (ChromBPNet methodology)
# =============================================================================

def load_valid_chroms_from_fasta(genome_fa):
    """Load valid chromosome names from reference FASTA."""
    with pyfaidx.Fasta(genome_fa) as g:
        return set(g.keys())

def bam_to_bigwig(bam_path, output_prefix, valid_chroms):
    """
    Convert BAM to base-resolution BigWig for DNase-seq.
    
    Replicates ChromBPNet's reads_to_bigwig.py methodology:
    - No BAM quality filtering (matches ChromBPNet)
    - Filter chromosomes to those in reference FASTA
    - Apply DNase shifts: +0 for plus strand start, +1 for minus strand end
    - Count 5' ends with bedtools genomecov -5
    """
    output_bw = f"{output_prefix}_unstranded.bw"
    
    if os.path.exists(output_bw):
        logging.info(f"    Already exists: {os.path.basename(output_bw)}")
        return output_bw
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bedGraph', delete=False) as tmp_bg:
        bedgraph_path = tmp_bg.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_bed:
        tmp_bed_path = tmp_bed.name
    
    try:
        logging.info(f"    Converting BAM to shifted BED...")
        
        # BAM to BED (no quality filtering, matching ChromBPNet)
        p1 = subprocess.Popen(
            ["bedtools", "bamtobed", "-i", bam_path],
            stdout=subprocess.PIPE
        )
        
        # Filter chromosomes and apply shifts
        with open(tmp_bed_path, 'w') as out_f:
            for line in p1.stdout:
                fields = line.decode('utf-8').strip().split('\t')
                chrom = fields[0]
                
                if chrom not in valid_chroms:
                    continue
                
                start = int(fields[1])
                end = int(fields[2])
                name = fields[3]
                score = fields[4]
                strand = fields[5]
                
                # Apply shifts (ChromBPNet style)
                if strand == "+":
                    new_start = start + PLUS_SHIFT_DELTA
                    new_end = end
                else:
                    new_start = start
                    new_end = end + MINUS_SHIFT_DELTA
                
                out_f.write(f"{chrom}\t{new_start}\t{new_end}\t{name}\t{score}\t{strand}\n")
        
        p1.wait()
        
        logging.info(f"    Generating coverage bedGraph...")
        cmd = f"""
        sort -k1,1 {tmp_bed_path} | \
        bedtools genomecov -bg -5 -i stdin -g {CHROM_SIZES} | \
        LC_COLLATE="C" sort -k1,1 -k2,2n > {bedgraph_path}
        """
        subprocess.run(cmd, shell=True, check=True)
        
        logging.info(f"    Converting to BigWig...")
        subprocess.run([
            "bedGraphToBigWig", bedgraph_path, CHROM_SIZES, output_bw
        ], check=True)
        
        return output_bw
        
    finally:
        for f in [bedgraph_path, tmp_bed_path]:
            if os.path.exists(f):
                os.remove(f)

# =============================================================================
# BigWig Averaging and Normalization
# =============================================================================

def load_chrom_sizes():
    """Load chromosome sizes from file."""
    chrom_sizes = {}
    with open(CHROM_SIZES) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes

def average_bigwigs(bw_paths, output_path, chrom_sizes):
    """
    Average multiple BigWig files (replicates) into one.
    Returns aggregation stats for manifest.json.
    """
    n_replicates = len(bw_paths)
    logging.info(f"  Averaging {n_replicates} replicates...")
    
    bws = [pyBigWig.open(p) for p in bw_paths]
    all_nonzero = []
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bedGraph', delete=False) as tmp:
        bedgraph_path = tmp.name
    
    try:
        with open(bedgraph_path, 'w') as bg_out:
            for chrom in sorted(chrom_sizes.keys()):
                size = chrom_sizes[chrom]
                
                # Check all BigWigs have this chromosome
                if not all(chrom in bw.chroms() for bw in bws):
                    continue
                
                # Get values from all replicates
                all_vals = []
                for bw in bws:
                    vals = bw.values(chrom, 0, size)
                    vals = np.array(vals, dtype=np.float64)
                    vals = np.nan_to_num(vals, nan=0.0)
                    all_vals.append(vals)
                
                # Average
                avg_vals = np.mean(all_vals, axis=0)
                
                # Track non-zero for stats
                nonzero = avg_vals[avg_vals > 0]
                all_nonzero.extend(nonzero)
                
                # Write to bedGraph (run-length encoding for efficiency)
                i = 0
                while i < len(avg_vals):
                    if avg_vals[i] > 0:
                        start = i
                        val = avg_vals[i]
                        while i < len(avg_vals) and avg_vals[i] == val:
                            i += 1
                        bg_out.write(f"{chrom}\t{start}\t{i}\t{val}\n")
                    else:
                        i += 1
        
        for bw in bws:
            bw.close()
        
        # Convert to BigWig
        subprocess.run([
            "bedGraphToBigWig", bedgraph_path, CHROM_SIZES, output_path
        ], check=True)
        
        non_zero_average_mean = float(np.mean(all_nonzero)) if all_nonzero else 0.0
        
        return {
            "non_zero_average_mean": non_zero_average_mean,
            "n_replicates": n_replicates
        }
        
    finally:
        if os.path.exists(bedgraph_path):
            os.remove(bedgraph_path)

def normalize_bigwig(input_path, output_path, chrom_sizes, target_sum=TARGET_SUM):
    """
    Normalize BigWig so total signal = target_sum (100M).
    Returns normalization stats for manifest.json.
    """
    logging.info(f"  Normalizing to {target_sum:.0e} total counts...")
    
    bw_in = pyBigWig.open(input_path)
    
    # Calculate current total
    pre_total = 0.0
    values_by_chrom = {}
    
    for chrom in sorted(chrom_sizes.keys()):
        if chrom not in bw_in.chroms():
            continue
        size = chrom_sizes[chrom]
        vals = bw_in.values(chrom, 0, size)
        vals = np.array(vals, dtype=np.float64)
        vals = np.nan_to_num(vals, nan=0.0)
        pre_total += np.sum(vals)
        values_by_chrom[chrom] = vals
    
    bw_in.close()
    
    # Calculate scaling factor
    scaling_factor = target_sum / pre_total
    logging.info(f"    Pre-total: {pre_total:.2e}, scaling factor: {scaling_factor:.6f}")
    
    # Write normalized values
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bedGraph', delete=False) as tmp:
        bedgraph_path = tmp.name
    
    try:
        with open(bedgraph_path, 'w') as bg_out:
            for chrom in sorted(values_by_chrom.keys()):
                vals = values_by_chrom[chrom] * scaling_factor
                
                i = 0
                while i < len(vals):
                    if vals[i] > 0:
                        start = i
                        val = vals[i]
                        while i < len(vals) and vals[i] == val:
                            i += 1
                        bg_out.write(f"{chrom}\t{start}\t{i}\t{val}\n")
                    else:
                        i += 1
        
        subprocess.run([
            "bedGraphToBigWig", bedgraph_path, CHROM_SIZES, output_path
        ], check=True)
        
        return {
            "pre_total": float(pre_total),
            "post_total": float(target_sum),
            "scaling_factor": float(scaling_factor)
        }
        
    finally:
        if os.path.exists(bedgraph_path):
            os.remove(bedgraph_path)

# =============================================================================
# Validation
# =============================================================================

def validate_bigwig(bw_path, chrom_sizes, expected_sum=TARGET_SUM, tolerance=0.001):
    """
    Validate the final BigWig file.
    
    Checks:
    1. File exists and is readable
    2. Chromosome names are valid
    3. Total signal sum is ~10^8
    """
    logging.info(f"  Validating {os.path.basename(bw_path)}...")
    
    if not os.path.exists(bw_path):
        raise ValueError(f"File does not exist: {bw_path}")
    
    try:
        bw = pyBigWig.open(bw_path)
    except Exception as e:
        raise ValueError(f"Cannot open BigWig: {e}")
    
    if bw is None:
        raise ValueError(f"BigWig is None: {bw_path}")
    
    # Check chromosomes
    bw_chroms = set(bw.chroms().keys())
    expected_chroms = set(chrom_sizes.keys())
    unexpected = bw_chroms - expected_chroms
    if unexpected:
        bw.close()
        raise ValueError(f"Unexpected chromosomes: {unexpected}")
    
    # Check total signal
    total = 0.0
    for chrom, size in bw.chroms().items():
        vals = bw.values(chrom, 0, size)
        vals = np.nan_to_num(vals, nan=0.0)
        total += np.sum(vals)
    
    bw.close()
    
    relative_error = abs(total - expected_sum) / expected_sum
    if relative_error > tolerance:
        raise ValueError(
            f"Total signal {total:.6e} differs from expected {expected_sum:.0e} "
            f"by {relative_error*100:.4f}% (tolerance: {tolerance*100}%)"
        )
    
    logging.info(f"    ✓ Valid (total={total:.6e})")
    return total

# =============================================================================
# Manifest Generation
# =============================================================================

def write_manifest(output_path, cell_line, file_accessions, chrombpnet_commit,
                   normalization_stats, aggregation_stats):
    """Write manifest.json with all required fields."""
    manifest = {
        "cell_line": cell_line,
        "replicates": file_accessions,
        "chrombpnet": {
            "repo_url": CHROMBPNET_REPO,
            "commit_sha": chrombpnet_commit
        },
        "normalization": {
            "pre_total": normalization_stats["pre_total"],
            "post_total": normalization_stats["post_total"],
            "scaling_factor": normalization_stats["scaling_factor"]
        },
        "aggregation": {
            "non_zero_average_mean": aggregation_stats["non_zero_average_mean"],
            "n_replicates": aggregation_stats["n_replicates"]
        }
    }
    
    with open(output_path, "w") as f:
        json.dump(manifest, f, indent=2)
    
    logging.info(f"  Wrote {output_path}")

# =============================================================================
# Main Pipeline
# =============================================================================

def get_chrombpnet_commit():
    """Get ChromBPNet commit SHA from cloned repo."""
    chrombpnet_dir = "chrombpnet"
    if os.path.exists(os.path.join(chrombpnet_dir, ".git")):
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=chrombpnet_dir,
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            return result.stdout.strip()
    return None

def process_cell_line(cell_line, skip_download=False, chrombpnet_commit=None):
    """Process one cell line through the full pipeline."""
    out_dir = f"out/{cell_line}"
    os.makedirs(out_dir, exist_ok=True)
    
    log_path = os.path.join(out_dir, "pipeline.log")
    setup_logging(log_path)
    
    logging.info(f"{'='*60}")
    logging.info(f"Processing {cell_line}")
    logging.info(f"{'='*60}")
    logging.info(f"Started at: {datetime.now().isoformat()}")
    
    try:
        # Step 1: Query ENCODE for BAM files
        file_infos = query_encode_for_cell_line(cell_line)
        
        if not file_infos:
            raise ValueError(f"No valid BAM files found for {cell_line}")
        
        # Step 2: Download BAMs
        bam_dir = "bams"
        os.makedirs(bam_dir, exist_ok=True)
        
        if not skip_download:
            logging.info("Downloading BAM files...")
            for file_info in file_infos:
                download_bam(file_info, bam_dir)
        else:
            logging.info("Skipping download (--skip-download)")
        
        # Step 3: Convert BAMs to BigWigs
        logging.info("Converting BAMs to BigWigs...")
        valid_chroms = load_valid_chroms_from_fasta(GENOME_FA)
        logging.info(f"  Loaded {len(valid_chroms)} chromosomes from reference")
        
        bigwig_dir = "bigwigs"
        os.makedirs(bigwig_dir, exist_ok=True)
        
        replicate_bws = []
        file_accessions = []
        
        for file_info in file_infos:
            file_acc = file_info["file_accession"]
            bam_path = os.path.join(bam_dir, f"{file_acc}.bam")
            
            if not os.path.exists(bam_path):
                logging.warning(f"  BAM not found: {bam_path}, skipping")
                continue
            
            logging.info(f"  Processing {file_acc}...")
            bw_prefix = os.path.join(bigwig_dir, file_acc)
            bw_path = bam_to_bigwig(bam_path, bw_prefix, valid_chroms)
            
            replicate_bws.append(bw_path)
            file_accessions.append(file_acc)
        
        if not replicate_bws:
            raise ValueError("No BigWigs generated")
        
        # Step 4: Average replicates
        logging.info("Averaging replicates...")
        chrom_sizes = load_chrom_sizes()
        avg_bw_path = os.path.join(out_dir, "dnase_avg.bw")
        aggregation_stats = average_bigwigs(replicate_bws, avg_bw_path, chrom_sizes)
        logging.info(f"  Aggregation: {aggregation_stats}")
        
        # Step 5: Normalize to 100M
        logging.info("Normalizing to 100M counts...")
        final_bw_path = os.path.join(out_dir, "dnase_avg_norm100M.bw")
        normalization_stats = normalize_bigwig(avg_bw_path, final_bw_path, chrom_sizes)
        logging.info(f"  Normalization: {normalization_stats}")
        
        # Clean up intermediate averaged BigWig
        os.remove(avg_bw_path)
        
        # Step 6: Validate
        logging.info("Validating final BigWig...")
        validate_bigwig(final_bw_path, chrom_sizes)
        
        # Step 7: Write manifest
        logging.info("Writing manifest.json...")
        manifest_path = os.path.join(out_dir, "manifest.json")
        write_manifest(
            manifest_path,
            cell_line,
            file_accessions,
            chrombpnet_commit,
            normalization_stats,
            aggregation_stats
        )
        
        logging.info(f"{'='*60}")
        logging.info(f"SUCCESS: {cell_line} completed")
        logging.info(f"Output: {final_bw_path}")
        logging.info(f"{'='*60}")
        
        return True
        
    except Exception as e:
        logging.error(f"FAILED: {e}")
        logging.exception("Traceback:")
        return False

def main():
    parser = argparse.ArgumentParser(
        description="DNase preprocessing pipeline following AlphaGenome methodology"
    )
    parser.add_argument(
        "--cell-lines",
        nargs="+",
        default=CELL_LINES,
        help=f"Cell lines to process (default: {CELL_LINES})"
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip BAM download (use existing files in bams/)"
    )
    parser.add_argument(
        "--chrombpnet-commit",
        type=str,
        default=None,
        help="ChromBPNet commit SHA (auto-detected from chrombpnet/ if not specified)"
    )
    args = parser.parse_args()
    
    # Get ChromBPNet commit
    chrombpnet_commit = args.chrombpnet_commit or get_chrombpnet_commit()
    if not chrombpnet_commit:
        print("WARNING: Could not determine ChromBPNet commit SHA")
        print("  Clone the repo: git clone https://github.com/kundajelab/chrombpnet.git")
        print("  Or specify: --chrombpnet-commit <sha>")
        chrombpnet_commit = "unknown"
    else:
        print(f"ChromBPNet commit: {chrombpnet_commit}")
    
    # Check required files
    if not os.path.exists(GENOME_FA):
        print(f"ERROR: Reference genome not found: {GENOME_FA}")
        print("  Download: wget -O reference/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz && gunzip reference/hg38.fa.gz")
        sys.exit(1)
    
    if not os.path.exists(CHROM_SIZES):
        print(f"ERROR: Chrom sizes not found: {CHROM_SIZES}")
        print("  Download: wget -O reference/hg38.chrom.sizes https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes")
        sys.exit(1)
    
    # Process each cell line
    success = True
    for cell_line in args.cell_lines:
        result = process_cell_line(
            cell_line,
            skip_download=args.skip_download,
            chrombpnet_commit=chrombpnet_commit
        )
        if not result:
            success = False
    
    # Final status
    if success:
        print("\n" + "="*60)
        print("ALL CELL LINES PROCESSED SUCCESSFULLY")
        print("="*60)
        sys.exit(0)
    else:
        print("\n" + "="*60)
        print("SOME CELL LINES FAILED - check pipeline.log files")
        print("="*60)
        sys.exit(1)

if __name__ == "__main__":
    main()