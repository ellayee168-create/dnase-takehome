#!/usr/bin/env python
"""
DNase preprocessing pipeline following AlphaGenome methodology.

Supports two modes:
1. ENCODE mode (default): Query ENCODE API for BAM files
2. Local mode: Process local BAM files specified in metadata.tsv

QC Filters Applied (ENCODE mode only):
1. FRiP (spot1_score) > 10% for DNase-seq
2. Minimum FASTQ read_length >= 36 nt for DNase-seq

Processing Steps:
1. Get BAM files (from ENCODE or local metadata)
2. Convert BAMs to base-resolution BigWigs (ChromBPNet methodology)
3. Average replicates within each cell line
4. Normalize to 100M total counts
5. Validate output and generate manifest.json

Usage:
    # ENCODE mode (query ENCODE for data):
    python run_dnase_tracks.py --cell-lines GM12878 HeLa-S3

    # Local mode (use pre-downloaded BAMs):
    python run_dnase_tracks.py --metadata metadata.tsv --chrom-sizes ref.chrom.sizes --outdir out

Output:
    <outdir>/<cell_line>/dnase_avg_norm100M.bw
    <outdir>/<cell_line>/manifest.json
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

# =============================================================================
# Configuration
# =============================================================================

CELL_LINES = ["GM12878", "HeLa-S3", "SK-N-SH"]
GENOME_FA = "reference/hg38.fa"
CHROM_SIZES = "reference/hg38.chrom.sizes"
TARGET_SUM = 1e8

# =============================================================================
# AlphaGenome QC Filter Thresholds
# =============================================================================
MIN_SPOT1_SCORE = 0.10
MIN_FASTQ_READ_LENGTH = 36
TARGET_PIPELINE = "dnase-alignment-filtering-step-v-1"

# DNase shifts
PLUS_SHIFT_DELTA = 0
MINUS_SHIFT_DELTA = 1

# ChromBPNet reference
CHROMBPNET_REPO = "https://github.com/kundajelab/chrombpnet"
CHROMBPNET_COMMIT = "ece97c93ccaa2d9ee5bc5687e62f4dbf8d055367"

# =============================================================================
# Logging Setup
# =============================================================================

def setup_logging(log_path=None):
    """Configure logging to both file and console."""
    root = logging.getLogger()
    for handler in root.handlers[:]:
        root.removeHandler(handler)

    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    if log_path:
        file_handler = logging.FileHandler(log_path, mode='w')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(formatter)
        root.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(formatter)

    root.setLevel(logging.INFO)
    root.addHandler(console_handler)

# =============================================================================
# Metadata Parsing (Local)
# =============================================================================

def parse_metadata_tsv(metadata_path):
    """
    Parse metadata.tsv file for local mode.
    
    Expected format (tab-separated):
    File accession    Biosample term name    Assay    bam_path
    rep1              TestCell               DNase-seq    /path/to/rep1.bam
    rep2              TestCell               DNase-seq    /path/to/rep2.bam
    
    Returns: dict mapping cell_line -> list of {"replicate_id": str, "bam_path": str}
    """
    cell_lines = {}
    
    with open(metadata_path, 'r') as f:
        header = f.readline().strip().split('\t')
        
        # For replicate/file accession
        if 'File accession' in header:
            replicate_idx = header.index('File accession')
        elif 'replicate_id' in header:
            replicate_idx = header.index('replicate_id')
        else:
            replicate_idx = 0
        
        # For cell line / biosample
        if 'Biosample term name' in header:
            cell_line_idx = header.index('Biosample term name')
        elif 'cell_line' in header:
            cell_line_idx = header.index('cell_line')
        else:
            cell_line_idx = 1
        
        # For BAM path
        if 'bam_path' in header:
            bam_idx = header.index('bam_path')
        else:
            bam_idx = len(header) - 1  # Assume last column
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < max(replicate_idx, cell_line_idx, bam_idx) + 1:
                continue
            
            replicate_id = fields[replicate_idx]
            cell_line = fields[cell_line_idx]
            bam_path = fields[bam_idx]
            
            if cell_line not in cell_lines:
                cell_lines[cell_line] = []
            
            cell_lines[cell_line].append({
                "replicate_id": replicate_id,
                "bam_path": bam_path
            })
    
    return cell_lines

# =============================================================================
# ENCODE Query
# =============================================================================

def get_json(url):
    """Get JSON from ENCODE API with rate limiting."""
    response = requests.get(url, headers={"Accept": "application/json"})
    time.sleep(0.1)
    return response.json()

def get_spot1_score(file_detail):
    """Get spot1_score from HotspotQualityMetric."""
    for qm in file_detail.get("quality_metrics", []):
        qm_type = qm.get("@type", [])
        if "HotspotQualityMetric" in qm_type:
            return qm.get("spot1_score")
    return None

def get_file_pipeline(file_detail):
    """Get pipeline name from file's analysis_step_version."""
    asv = file_detail.get("analysis_step_version")
    if isinstance(asv, dict):
        analysis_step = asv.get("analysis_step", {})
        if isinstance(analysis_step, dict):
            return analysis_step.get("name", "")
    return ""

def get_experiment_info(exp_acc):
    """Get experiment metadata including FASTQ read lengths and run type."""
    exp_data = get_json(f"https://www.encodeproject.org/experiments/{exp_acc}/?format=json")
    
    fastq_read_lengths = set()
    run_types = set()
    
    for f in exp_data.get("files", []):
        if isinstance(f, str):
            file_data = get_json(f"https://www.encodeproject.org{f}?format=json")
        else:
            file_data = f
        
        if file_data.get("file_format") == "fastq":
            rl = file_data.get("read_length")
            if rl is not None:
                fastq_read_lengths.add(rl)
            rt = file_data.get("run_type")
            if rt:
                run_types.add(rt)
    
    return {
        "fastq_read_lengths": fastq_read_lengths,
        "is_paired": "paired-ended" in run_types,
        "date_released": exp_data.get("date_released", "")
    }

def query_encode_for_cell_line(cell_line):
    """Query ENCODE for DNase BAM files for one cell line."""
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

    experiments = {}
    
    for f in files:
        file_acc = f.get("accession")
        file_detail = get_json(f"https://www.encodeproject.org/files/{file_acc}/?format=json")

        dataset = file_detail.get("dataset", "")
        exp_acc = dataset.split("/")[-2] if dataset else None
        if not exp_acc:
            continue

        spot1_score = get_spot1_score(file_detail)
        if spot1_score is None or spot1_score <= MIN_SPOT1_SCORE:
            continue

        pipeline = get_file_pipeline(file_detail)
        bio_reps = file_detail.get("biological_replicates", [])
        bio_rep = bio_reps[0] if bio_reps else None

        href = file_detail.get("href")
        download_url = f"https://www.encodeproject.org{href}" if href else \
                       f"https://www.encodeproject.org/files/{file_acc}/@@download/{file_acc}.bam"

        file_info = {
            "file_accession": file_acc,
            "experiment_accession": exp_acc,
            "bam_url": download_url,
            "spot1_score": spot1_score,
            "pipeline": pipeline,
            "biological_replicate": bio_rep,
            "date_created": file_detail.get("date_created", "")
        }

        if exp_acc not in experiments:
            experiments[exp_acc] = []
        experiments[exp_acc].append(file_info)

    logging.info(f"  Found {len(experiments)} experiments with files passing spot1 filter")

    valid_experiments = []
    
    for exp_acc, exp_files in experiments.items():
        exp_info = get_experiment_info(exp_acc)
        fastq_rls = exp_info["fastq_read_lengths"]
        
        if not fastq_rls:
            logging.info(f"  Excluding {exp_acc}: no FASTQ read length info")
            continue
        
        min_rl = min(fastq_rls)
        if min_rl < MIN_FASTQ_READ_LENGTH:
            logging.info(f"  Excluding {exp_acc}: min FASTQ read_length={min_rl} < {MIN_FASTQ_READ_LENGTH}")
            continue
        
        valid_experiments.append({
            "exp_acc": exp_acc,
            "files": exp_files,
            "is_paired": exp_info["is_paired"],
            "date_released": exp_info["date_released"],
            "fastq_read_lengths": fastq_rls
        })
        logging.info(f"  Valid: {exp_acc} (FASTQ read_lengths={fastq_rls}, paired={exp_info['is_paired']})")

    if not valid_experiments:
        logging.error("No experiments passed FASTQ read length filter")
        return []

    valid_experiments.sort(key=lambda x: (x["is_paired"], x["date_released"]), reverse=True)
    
    selected = valid_experiments[0]
    logging.info(f"  Selected experiment: {selected['exp_acc']} (paired={selected['is_paired']}, date={selected['date_released']})")

    files_by_bio_rep = {}
    for f in selected["files"]:
        bio_rep = f["biological_replicate"]
        if bio_rep is None:
            continue
        if bio_rep not in files_by_bio_rep:
            files_by_bio_rep[bio_rep] = []
        files_by_bio_rep[bio_rep].append(f)

    selected_files = []
    for bio_rep in sorted(files_by_bio_rep.keys()):
        rep_files = files_by_bio_rep[bio_rep]
        pipeline_files = [f for f in rep_files if TARGET_PIPELINE in f["pipeline"]]
        if pipeline_files:
            chosen = sorted(pipeline_files, key=lambda x: x["date_created"], reverse=True)[0]
        else:
            chosen = sorted(rep_files, key=lambda x: x["date_created"], reverse=True)[0]

        selected_files.append(chosen)
        logging.info(f"    Bio rep {bio_rep}: {chosen['file_accession']} (pipeline={chosen['pipeline']})")

    logging.info(f"  Selected {len(selected_files)} files (one per biological replicate)")
    return selected_files

# =============================================================================
# Download BAM Files
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
# Load Chromosome Sizes
# =============================================================================

def load_chrom_sizes_from_file(chrom_sizes_path):
    """Load chromosome sizes from file."""
    chrom_sizes = {}
    with open(chrom_sizes_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 2:
                chrom_sizes[fields[0]] = int(fields[1])
    return chrom_sizes

def load_valid_chroms_from_chrom_sizes(chrom_sizes_path):
    """Load valid chromosome names from chrom.sizes file."""
    return set(load_chrom_sizes_from_file(chrom_sizes_path).keys())

def load_valid_chroms_from_fasta(genome_fa):
    """Load valid chromosome names from reference FASTA."""
    import pyfaidx
    with pyfaidx.Fasta(genome_fa) as g:
        return set(g.keys())

# =============================================================================
# BAM to BigWig Conversion (ChromBPNet Methodology)
# =============================================================================

def bam_to_bigwig(bam_path, output_prefix, valid_chroms, chrom_sizes):
    """Convert BAM to base-resolution BigWig for DNase-seq."""
    output_bw = f"{output_prefix}_unstranded.bw"

    if os.path.exists(output_bw):
        logging.info(f"    Already exists: {os.path.basename(output_bw)}")
        return output_bw

    with tempfile.NamedTemporaryFile(mode='w', suffix='.bedGraph', delete=False) as tmp_bg:
        bedgraph_path = tmp_bg.name
    with tempfile.NamedTemporaryFile(mode='w', suffix='.bed', delete=False) as tmp_bed:
        tmp_bed_path = tmp_bed.name
    
    # Create temp chrom.sizes file for bedtools
    with tempfile.NamedTemporaryFile(mode='w', suffix='.chrom.sizes', delete=False) as tmp_cs:
        tmp_chrom_sizes_path = tmp_cs.name
        for chrom in sorted(chrom_sizes.keys()):
            tmp_cs.write(f"{chrom}\t{chrom_sizes[chrom]}\n")

    try:
        logging.info("    Converting BAM to shifted BED...")

        p1 = subprocess.Popen(
            ["bedtools", "bamtobed", "-i", bam_path],
            stdout=subprocess.PIPE
        )

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

                if strand == "+":
                    new_start = start + PLUS_SHIFT_DELTA
                    new_end = new_start + 1
                else:
                    new_start = end - 1 + MINUS_SHIFT_DELTA
                    new_end = new_start + 1

                out_f.write(f"{chrom}\t{new_start}\t{new_end}\t{name}\t{score}\t{strand}\n")

        p1.wait()

        logging.info("    Generating coverage bedGraph...")

        sorted_bed_path = tmp_bed_path + ".sorted"
        subprocess.run(
            f"sort -k1,1 -k2,2n {tmp_bed_path} > {sorted_bed_path}",
            shell=True, check=True
        )

        subprocess.run(
            f"bedtools genomecov -i {sorted_bed_path} -g {tmp_chrom_sizes_path} -bg > {bedgraph_path}",
            shell=True, check=True
        )

        logging.info("    Converting bedGraph to BigWig...")

        intervals_by_chrom = {}

        with open(bedgraph_path, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                value = float(fields[3])

                if chrom not in intervals_by_chrom:
                    intervals_by_chrom[chrom] = []
                intervals_by_chrom[chrom].append((start, end, value))

        write_bigwig_direct(intervals_by_chrom, chrom_sizes, output_bw)

        logging.info(f"    Created: {os.path.basename(output_bw)}")
        return output_bw

    finally:
        for tmp_file in [tmp_bed_path, tmp_bed_path + ".sorted", bedgraph_path, tmp_chrom_sizes_path]:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)

def write_bigwig_direct(intervals_by_chrom, chrom_sizes, output_path):
    """Write BigWig directly using pyBigWig."""
    sorted_chroms = sorted(intervals_by_chrom.keys())
    header = [(chrom, chrom_sizes[chrom]) for chrom in sorted_chroms if chrom in chrom_sizes]

    bw = pyBigWig.open(output_path, "w")
    bw.addHeader(header)

    for chrom in sorted_chroms:
        intervals = intervals_by_chrom[chrom]
        if not intervals:
            continue

        intervals = sorted(intervals, key=lambda x: x[0])

        starts = [i[0] for i in intervals]
        ends = [i[1] for i in intervals]
        values = [float(i[2]) for i in intervals]

        bw.addEntries([chrom] * len(starts), starts, ends=ends, values=values)

    bw.close()

def array_to_intervals(arr):
    """Convert numpy array to list of (start, end, value) intervals."""
    intervals = []
    n = len(arr)
    if n == 0:
        return intervals

    i = 0
    while i < n:
        val = arr[i]
        if val == 0:
            i += 1
            continue

        start = i
        while i < n and arr[i] == val:
            i += 1
        end = i

        intervals.append((start, end, val))

    return intervals

# =============================================================================
# Aggregation & Normalization
# =============================================================================

def average_bigwigs(bw_paths, output_path, chrom_sizes):
    """Average multiple BigWig files."""
    if len(bw_paths) == 1:
        import shutil
        shutil.copy(bw_paths[0], output_path)
        return {"n_replicates": 1}

    n_replicates = len(bw_paths)
    logging.info(f"  Averaging {n_replicates} replicates...")

    bws = [pyBigWig.open(p) for p in bw_paths]
    intervals_by_chrom = {}

    for chrom in sorted(chrom_sizes.keys()):
        size = chrom_sizes[chrom]

        if not all(chrom in bw.chroms() for bw in bws):
            continue

        all_vals = []
        for bw in bws:
            vals = bw.values(chrom, 0, size)
            vals = np.array(vals, dtype=np.float64)
            vals = np.nan_to_num(vals, nan=0.0)
            all_vals.append(vals)

        avg_vals = np.mean(all_vals, axis=0)
        intervals_by_chrom[chrom] = array_to_intervals(avg_vals)

    for bw in bws:
        bw.close()

    write_bigwig_direct(intervals_by_chrom, chrom_sizes, output_path)

    return {"n_replicates": n_replicates}

def normalize_bigwig(input_path, output_path, chrom_sizes, target_sum=TARGET_SUM):
    """Normalize BigWig so total signal = target_sum (100M)."""
    logging.info(f"  Normalizing to {target_sum:.0e} total counts...")

    bw_in = pyBigWig.open(input_path)

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

    scaling_factor = target_sum / pre_total
    logging.info(f"    Pre-total: {pre_total:.2e}, scaling factor: {scaling_factor:.6f}")

    intervals_by_chrom = {}
    for chrom in sorted(values_by_chrom.keys()):
        scaled_vals = values_by_chrom[chrom] * scaling_factor
        intervals_by_chrom[chrom] = array_to_intervals(scaled_vals)

    write_bigwig_direct(intervals_by_chrom, chrom_sizes, output_path)

    return {
        "pre_total": float(pre_total),
        "post_total": float(target_sum),
        "scaling_factor": float(scaling_factor)
    }

def compute_non_zero_average_mean(bw_path):
    """
    Compute non_zero_average_mean from a BigWig file.
    
    Computes the unweighted mean of non-zero interval values:
    np.mean([value for each non-zero interval])
    """
    bw = pyBigWig.open(bw_path)

    non_zero_values = []

    for chrom in bw.chroms():
        intervals = bw.intervals(chrom)
        if intervals:
            for start, end, value in intervals:
                if value != 0:
                    non_zero_values.append(value)

    bw.close()

    return float(np.mean(non_zero_values)) if non_zero_values else 0.0

# =============================================================================
# Validation
# =============================================================================

def validate_bigwig(bw_path, chrom_sizes, expected_sum=TARGET_SUM, tolerance=0.001):
    """Validate the final BigWig file."""
    logging.info(f"  Validating {os.path.basename(bw_path)}...")

    if not os.path.exists(bw_path):
        raise ValueError(f"File does not exist: {bw_path}")

    bw = pyBigWig.open(bw_path)
    if bw is None:
        raise ValueError(f"BigWig is None: {bw_path}")

    bw_chroms = set(bw.chroms().keys())
    expected_chroms = set(chrom_sizes.keys())
    unexpected = bw_chroms - expected_chroms
    if unexpected:
        bw.close()
        raise ValueError(f"Unexpected chromosomes: {unexpected}")

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

    logging.info(f"    Valid (total={total:.6e})")
    return total

# =============================================================================
# Manifest Generation
# =============================================================================

def write_manifest(output_path, cell_line, file_accessions, experiment_accession,
                   chrombpnet_commit, normalization_stats, aggregation_stats):
    """Write manifest.json with all required fields."""
    manifest = {
        "cell_line": cell_line,
        "experiment": experiment_accession,
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
# Main Pipeline (Local)
# =============================================================================

def process_cell_line_local(cell_line, bam_infos, chrom_sizes_path, output_dir, chrombpnet_commit=None):
    """Process one cell line in local mode (from metadata.tsv)."""
    out_dir = os.path.join(output_dir, cell_line)
    os.makedirs(out_dir, exist_ok=True)

    log_path = os.path.join(out_dir, "pipeline.log")
    setup_logging(log_path)

    logging.info(f"{'='*60}")
    logging.info(f"Processing {cell_line} (local mode)")
    logging.info(f"{'='*60}")
    logging.info(f"Started at: {datetime.now().isoformat()}")

    try:
        chrom_sizes = load_chrom_sizes_from_file(chrom_sizes_path)
        valid_chroms = set(chrom_sizes.keys())
        logging.info(f"  Loaded {len(valid_chroms)} chromosomes from {chrom_sizes_path}")

        # Create temp directory for intermediate files
        bigwig_dir = tempfile.mkdtemp(prefix="bigwigs_")
        
        replicate_bws = []
        replicate_ids = []

        for bam_info in bam_infos:
            rep_id = bam_info["replicate_id"]
            bam_path = bam_info["bam_path"]

            if not os.path.exists(bam_path):
                logging.warning(f"  BAM not found: {bam_path}, skipping")
                continue

            logging.info(f"  Processing {rep_id} ({bam_path})...")
            bw_prefix = os.path.join(bigwig_dir, rep_id)
            bw_path = bam_to_bigwig(bam_path, bw_prefix, valid_chroms, chrom_sizes)

            replicate_bws.append(bw_path)
            replicate_ids.append(rep_id)

        if not replicate_bws:
            raise ValueError("No BigWigs generated")

        replicate_bws = sorted(replicate_bws)
        replicate_ids = sorted(replicate_ids)

        # Average replicates
        logging.info("Averaging replicates...")
        avg_bw_path = os.path.join(out_dir, "dnase_avg.bw")
        aggregation_stats = average_bigwigs(replicate_bws, avg_bw_path, chrom_sizes)

        # Normalize to 100M
        logging.info("Normalizing to 100M counts...")
        final_bw_path = os.path.join(out_dir, "dnase_avg_norm100M.bw")
        normalization_stats = normalize_bigwig(avg_bw_path, final_bw_path, chrom_sizes)

        os.remove(avg_bw_path)

        # Validate
        logging.info("Validating final BigWig...")
        validate_bigwig(final_bw_path, chrom_sizes)

        # Compute non_zero_average_mean
        logging.info("Computing aggregation statistics...")
        non_zero_avg_mean = compute_non_zero_average_mean(final_bw_path)
        aggregation_stats["non_zero_average_mean"] = non_zero_avg_mean
        logging.info(f"  non_zero_average_mean: {non_zero_avg_mean:.6f}")

        # Write manifest
        logging.info("Writing manifest.json...")
        manifest_path = os.path.join(out_dir, "manifest.json")
        write_manifest(
            manifest_path,
            cell_line,
            replicate_ids,
            "local",  # No experiment accession in local mode
            chrombpnet_commit,
            normalization_stats,
            aggregation_stats
        )

        # Cleanup temp bigwigs
        import shutil
        shutil.rmtree(bigwig_dir, ignore_errors=True)

        logging.info(f"{'='*60}")
        logging.info(f"SUCCESS: {cell_line} completed")
        logging.info(f"Output: {final_bw_path}")
        logging.info(f"{'='*60}")

        return True

    except Exception as e:
        logging.error(f"FAILED: {e}")
        logging.exception("Traceback:")
        return False

# =============================================================================
# Main Pipeline (ENCODE)
# =============================================================================

def process_cell_line_encode(cell_line, skip_download=False, chrombpnet_commit=None):
    """Process one cell line in ENCODE mode (query ENCODE API)."""
    out_dir = f"out/{cell_line}"
    os.makedirs(out_dir, exist_ok=True)

    log_path = os.path.join(out_dir, "pipeline.log")
    setup_logging(log_path)

    logging.info(f"{'='*60}")
    logging.info(f"Processing {cell_line} (ENCODE mode)")
    logging.info(f"{'='*60}")
    logging.info(f"Started at: {datetime.now().isoformat()}")

    try:
        # Query ENCODE for BAM files
        file_infos = query_encode_for_cell_line(cell_line)

        if not file_infos:
            raise ValueError(f"No valid BAM files found for {cell_line}")

        experiment_accession = file_infos[0]["experiment_accession"]

        # Download BAMs
        bam_dir = "bams"
        os.makedirs(bam_dir, exist_ok=True)

        if not skip_download:
            logging.info("Downloading BAM files...")
            for file_info in file_infos:
                download_bam(file_info, bam_dir)
        else:
            logging.info("Skipping download (--skip-download)")

        # Convert BAMs to BigWigs
        logging.info("Converting BAMs to BigWigs...")
        valid_chroms = load_valid_chroms_from_fasta(GENOME_FA)
        chrom_sizes = load_chrom_sizes_from_file(CHROM_SIZES)
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
            bw_path = bam_to_bigwig(bam_path, bw_prefix, valid_chroms, chrom_sizes)

            replicate_bws.append(bw_path)
            file_accessions.append(file_acc)

        if not replicate_bws:
            raise ValueError("No BigWigs generated")

        replicate_bws = sorted(replicate_bws)
        file_accessions = sorted(file_accessions)

        # Average replicates
        logging.info("Averaging replicates...")
        avg_bw_path = os.path.join(out_dir, "dnase_avg.bw")
        aggregation_stats = average_bigwigs(replicate_bws, avg_bw_path, chrom_sizes)

        # Normalize to 100M
        logging.info("Normalizing to 100M counts...")
        final_bw_path = os.path.join(out_dir, "dnase_avg_norm100M.bw")
        normalization_stats = normalize_bigwig(avg_bw_path, final_bw_path, chrom_sizes)

        os.remove(avg_bw_path)

        # Validate
        logging.info("Validating final BigWig...")
        validate_bigwig(final_bw_path, chrom_sizes)

        # Compute non_zero_average_mean
        logging.info("Computing aggregation statistics...")
        non_zero_avg_mean = compute_non_zero_average_mean(final_bw_path)
        aggregation_stats["non_zero_average_mean"] = non_zero_avg_mean
        logging.info(f"  non_zero_average_mean: {non_zero_avg_mean:.6f}")

        # Write manifest
        logging.info("Writing manifest.json...")
        manifest_path = os.path.join(out_dir, "manifest.json")
        write_manifest(
            manifest_path,
            cell_line,
            file_accessions,
            experiment_accession,
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

# =============================================================================
# Main Entry Point
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

def main():
    parser = argparse.ArgumentParser(
        description="DNase preprocessing pipeline following AlphaGenome methodology"
    )
    
    # Local mode arguments
    parser.add_argument(
        "--metadata",
        type=str,
        default=None,
        help="Path to metadata.tsv file (enables local mode)"
    )
    parser.add_argument(
        "--chrom-sizes",
        type=str,
        default=None,
        help="Path to chromosome sizes file (required for local mode)"
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="out",
        help="Output directory (default: out)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads (currently unused, for compatibility)"
    )
    
    # ENCODE mode arguments
    parser.add_argument(
        "--cell-lines",
        nargs="+",
        default=CELL_LINES,
        help=f"Cell lines to process in ENCODE mode (default: {CELL_LINES})"
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip BAM download in ENCODE mode (use existing files in bams/)"
    )
    parser.add_argument(
        "--chrombpnet-commit",
        type=str,
        default=None,
        help="ChromBPNet commit SHA"
    )
    
    args = parser.parse_args()

    chrombpnet_commit = args.chrombpnet_commit or get_chrombpnet_commit() or CHROMBPNET_COMMIT
    print(f"ChromBPNet commit: {chrombpnet_commit}")

    # Determine mode based on arguments
    if args.metadata:
        # Local mode
        if not args.chrom_sizes:
            print("ERROR: --chrom-sizes is required when using --metadata")
            sys.exit(1)
        
        if not os.path.exists(args.metadata):
            print(f"ERROR: Metadata file not found: {args.metadata}")
            sys.exit(1)
        
        if not os.path.exists(args.chrom_sizes):
            print(f"ERROR: Chrom sizes file not found: {args.chrom_sizes}")
            sys.exit(1)
        
        print("Running in LOCAL mode")
        print(f"  Metadata: {args.metadata}")
        print(f"  Chrom sizes: {args.chrom_sizes}")
        print(f"  Output dir: {args.outdir}")
        
        # Parse metadata
        cell_lines_data = parse_metadata_tsv(args.metadata)
        
        success = True
        for cell_line, bam_infos in cell_lines_data.items():
            result = process_cell_line_local(
                cell_line,
                bam_infos,
                args.chrom_sizes,
                args.outdir,
                chrombpnet_commit=chrombpnet_commit
            )
            if not result:
                success = False
        
    else:
        # ENCODE mode
        print("Running in ENCODE mode")
        
        if not os.path.exists(GENOME_FA):
            print(f"ERROR: Reference genome not found: {GENOME_FA}")
            sys.exit(1)

        if not os.path.exists(CHROM_SIZES):
            print(f"ERROR: Chrom sizes not found: {CHROM_SIZES}")
            sys.exit(1)

        success = True
        for cell_line in args.cell_lines:
            result = process_cell_line_encode(
                cell_line,
                skip_download=args.skip_download,
                chrombpnet_commit=chrombpnet_commit
            )
            if not result:
                success = False

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
