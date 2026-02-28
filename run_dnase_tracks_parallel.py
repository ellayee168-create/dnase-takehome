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

import shutil

import requests
import numpy as np
import pyBigWig
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

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

def setup_logging(log_path=None, cell_line=None):
    """Configure logging to both file and console.
    
    Args:
        log_path: Path to log file (optional)
        cell_line: Cell line name to prefix log messages (optional).
                   When running cell lines in parallel, this makes it clear
                   which cell line each message belongs to.
    """
    root = logging.getLogger()
    for handler in root.handlers[:]:
        root.removeHandler(handler)

    # Include cell line prefix in format if provided
    if cell_line:
        fmt = f'%(asctime)s - %(levelname)s - [{cell_line}] %(message)s'
    else:
        fmt = '%(asctime)s - %(levelname)s - %(message)s'
    
    formatter = logging.Formatter(fmt, datefmt='%Y-%m-%d %H:%M:%S')

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

def get_json(url, retries=3):
    """Get JSON from ENCODE API with rate limiting and retry on transient errors.

    Retries up to `retries` times with exponential backoff on any
    requests.exceptions.RequestException (network errors, timeouts, HTTP
    error status codes).
    """
    for attempt in range(retries):
        try:
            response = requests.get(url, headers={"Accept": "application/json"})
            response.raise_for_status()
            time.sleep(0.1)
            return response.json()
        except requests.exceptions.RequestException:
            if attempt == retries - 1:
                raise
            wait = 2 ** attempt  # 1s, 2s, 4s
            logging.warning(f"ENCODE request failed (attempt {attempt + 1}/{retries}), retrying in {wait}s: {url}")
            time.sleep(wait)

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

def _fetch_file_detail(acc):
    """Picklable module-level helper: fetch one file's detail JSON from ENCODE."""
    return acc, get_json(f"https://www.encodeproject.org/files/{acc}/?format=json")


def _fetch_exp_info(acc):
    """Picklable module-level helper: fetch one experiment's info from ENCODE."""
    return acc, get_experiment_info(acc)


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

def query_encode_for_cell_line(cell_line, max_api_workers=4):
    """Query ENCODE for DNase BAM files for one cell line.

    PARALLELISM: file-detail fetches and experiment-info fetches are done
    concurrently with ThreadPoolExecutor (network I/O-bound). All downstream
    selection/sorting logic is unchanged from the original serial version.

    max_api_workers is passed in from --threads so the user has one consistent
    knob for all parallelism in the pipeline.
    """
    logging.info("Querying ENCODE for DNase-seq BAM files...")

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

    # Initial search uses the same retry pattern as get_json: up to 3 attempts
    # with exponential backoff on any RequestException (network error, HTTP error).
    for attempt in range(3):
        try:
            response = requests.get(search_url, params=params, headers={"Accept": "application/json"})
            response.raise_for_status()
            break
        except requests.exceptions.RequestException:
            if attempt == 2:
                raise
            wait = 2 ** attempt
            logging.warning(f"ENCODE search request failed (attempt {attempt + 1}/3), retrying in {wait}s")
            time.sleep(wait)
    files = response.json().get("@graph", [])
    logging.info(f"Found {len(files)} BAM files in ENCODE")

    # Preserve API-returned order so the experiments dict is built consistently.
    file_accs = [f.get("accession") for f in files if f.get("accession")]

    # --- PARALLEL: fetch all file details concurrently (network I/O-bound) ---
    file_details = {}
    with ThreadPoolExecutor(max_workers=max_api_workers) as pool:
        for acc, detail in pool.map(_fetch_file_detail, file_accs):
            file_details[acc] = detail

    experiments = {}

    # Iterate in original API order for stable ordering within each experiment's list.
    for file_acc in file_accs:
        file_detail = file_details[file_acc]

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

    logging.info(f"Found {len(experiments)} experiments with files passing spot1 filter")

    # --- PARALLEL: fetch experiment info for all candidates concurrently ---
    exp_accs = list(experiments.keys())

    exp_infos = {}
    with ThreadPoolExecutor(max_workers=max_api_workers) as pool:
        for exp_acc, info in pool.map(_fetch_exp_info, exp_accs):
            exp_infos[exp_acc] = info

    valid_experiments = []

    for exp_acc, exp_files in experiments.items():
        exp_info = exp_infos[exp_acc]
        fastq_rls = exp_info["fastq_read_lengths"]

        if not fastq_rls:
            logging.info(f"Excluding {exp_acc}: no FASTQ read length info")
            continue

        min_rl = min(fastq_rls)
        if min_rl < MIN_FASTQ_READ_LENGTH:
            logging.info(f"Excluding {exp_acc}: min FASTQ read_length={min_rl} < {MIN_FASTQ_READ_LENGTH}")
            continue

        valid_experiments.append({
            "exp_acc": exp_acc,
            "files": exp_files,
            "is_paired": exp_info["is_paired"],
            "date_released": exp_info["date_released"],
            "fastq_read_lengths": fastq_rls
        })
        logging.info(f"Valid: {exp_acc} (FASTQ read_lengths={fastq_rls}, paired={exp_info['is_paired']})")

    if not valid_experiments:
        logging.error("No experiments passed FASTQ read length filter")
        return []

    valid_experiments.sort(key=lambda x: (x["is_paired"], x["date_released"]), reverse=True)
    
    selected = valid_experiments[0]
    logging.info(f"Selected experiment: {selected['exp_acc']} (paired={selected['is_paired']}, date={selected['date_released']})")

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
        logging.info(f"Bio rep {bio_rep}: {chosen['file_accession']} (pipeline={chosen['pipeline']})")

    logging.info(f"Selected {len(selected_files)} files (one per biological replicate)")
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
        logging.info(f"Already downloaded: {file_acc}.bam")
        return output_path

    logging.info(f"Downloading {file_acc}.bam...")
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
        logging.info(f"Already exists: {os.path.basename(output_bw)}")
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
        logging.info(f"Converting BAM to shifted BED: {os.path.basename(bam_path)}")

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

        logging.info(f"Generating coverage bedGraph: {os.path.basename(bam_path)}")

        sorted_bed_path = tmp_bed_path + ".sorted"
        subprocess.run(
            f"sort -k1,1 -k2,2n {tmp_bed_path} > {sorted_bed_path}",
            shell=True, check=True
        )

        subprocess.run(
            f"bedtools genomecov -i {sorted_bed_path} -g {tmp_chrom_sizes_path} -bg > {bedgraph_path}",
            shell=True, check=True
        )

        logging.info(f"Converting bedGraph to BigWig: {os.path.basename(bam_path)}")

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

        logging.info(f"Created: {os.path.basename(output_bw)}")
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
        shutil.copy(bw_paths[0], output_path)
        return {"n_replicates": 1}

    n_replicates = len(bw_paths)
    logging.info(f"Averaging {n_replicates} replicates...")

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
    """Normalize BigWig so total signal = target_sum (100M).

    Uses a two-pass approach to avoid holding the entire genome in memory:
      Pass 1: read each chromosome to accumulate pre_total, then discard the array.
      Pass 2: read each chromosome again, apply the scaling factor, write, discard.
    """
    logging.info(f"Normalizing to {target_sum:.0e} total counts...")

    # --- Pass 1: compute pre_total only; discard arrays immediately ---
    bw_in = pyBigWig.open(input_path)

    pre_total = 0.0
    chroms_present = []  # ordered list of chroms actually in the BigWig

    for chrom in sorted(chrom_sizes.keys()):
        if chrom not in bw_in.chroms():
            continue
        size = chrom_sizes[chrom]
        vals = bw_in.values(chrom, 0, size)
        vals = np.array(vals, dtype=np.float64)
        vals = np.nan_to_num(vals, nan=0.0)
        pre_total += np.sum(vals)
        chroms_present.append(chrom)

    bw_in.close()

    if pre_total == 0.0:
        raise ValueError(
            f"Total signal in {input_path} is zero; cannot normalize. "
            "Check that your BAM file has reads mapping to valid chromosomes."
        )

    scaling_factor = target_sum / pre_total
    logging.info(f"Pre-total: {pre_total:.2e}, scaling factor: {scaling_factor:.6f}")

    # --- Pass 2: read, scale, convert to intervals, write; one chrom at a time ---
    bw_in = pyBigWig.open(input_path)

    intervals_by_chrom = {}
    for chrom in chroms_present:
        size = chrom_sizes[chrom]
        vals = bw_in.values(chrom, 0, size)
        vals = np.array(vals, dtype=np.float64)
        vals = np.nan_to_num(vals, nan=0.0)
        scaled_vals = vals * scaling_factor
        intervals_by_chrom[chrom] = array_to_intervals(scaled_vals)

    bw_in.close()

    write_bigwig_direct(intervals_by_chrom, chrom_sizes, output_path)

    return {
        "pre_total": float(pre_total),
        "post_total": float(target_sum),
        "scaling_factor": float(scaling_factor)
    }

def compute_non_zero_average_mean(bw_path):
    """Computes the weighted mean of a BigWig file."""
    bw = pyBigWig.open(bw_path)
    weighted_sum = 0.0
    total_length = 0
    
    for chrom in bw.chroms():
        intervals = bw.intervals(chrom)
        if intervals:
            for start, end, value in intervals:
                if value != 0:
                    length = end - start
                    weighted_sum += value * length
                    total_length += length
    
    bw.close()
    return weighted_sum / total_length if total_length > 0 else 0.0

# =============================================================================
# Validation
# =============================================================================

def validate_bigwig(bw_path, chrom_sizes, expected_sum=TARGET_SUM, tolerance=0.001):
    """Validate the final BigWig file."""
    logging.info(f"Validating {os.path.basename(bw_path)}...")

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

    logging.info(f"Valid (total={total:.6e})")
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

    logging.info(f"Wrote {output_path}")

# =============================================================================
# Top-level picklable worker wrappers for ProcessPoolExecutor
# (lambdas and nested functions cannot be pickled for inter-process use)
# =============================================================================

def _download_bam_worker(args):
    """Unpacking wrapper for download_bam used with pool.map."""
    file_info, output_dir = args
    return download_bam(file_info, output_dir)


def _bam_to_bigwig_worker(args):
    """Picklable wrapper for bam_to_bigwig."""
    bam_path, bw_prefix, valid_chroms, chrom_sizes = args
    return bam_to_bigwig(bam_path, bw_prefix, valid_chroms, chrom_sizes)


def _process_cell_line_local_worker(args):
    """Picklable wrapper for process_cell_line_local."""
    cell_line, bam_infos, chrom_sizes_path, output_dir, chrombpnet_commit, max_workers = args
    return process_cell_line_local(
        cell_line, bam_infos, chrom_sizes_path, output_dir,
        chrombpnet_commit=chrombpnet_commit, max_workers=max_workers,
    )


def _process_cell_line_encode_worker(args):
    """Picklable wrapper for process_cell_line_encode."""
    cell_line, skip_download, chrombpnet_commit, max_workers = args
    return process_cell_line_encode(
        cell_line, skip_download=skip_download,
        chrombpnet_commit=chrombpnet_commit, max_workers=max_workers,
    )


# =============================================================================
# Shared finalization: average -> normalize -> validate -> manifest
# =============================================================================

def _finalize_cell_line(replicate_bws, replicate_ids, chrom_sizes, out_dir,
                        cell_line, experiment_accession, chrombpnet_commit):
    """Average replicates, normalize, validate, compute stats, and write manifest.

    Called by both process_cell_line_local and process_cell_line_encode after
    each mode has assembled its list of per-replicate BigWig paths.

    Args:
        replicate_bws: list of per-replicate BigWig paths (already sorted)
        replicate_ids: list of replicate identifiers (already sorted)
        chrom_sizes: dict of {chrom: size}
        out_dir: output directory for this cell line
        cell_line: cell line name (for manifest)
        experiment_accession: experiment accession string (for manifest)
        chrombpnet_commit: ChromBPNet commit SHA (for manifest)
    """
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
    logging.info(f"non_zero_average_mean: {non_zero_avg_mean:.6f}")

    # Write manifest
    logging.info("Writing manifest.json...")
    manifest_path = os.path.join(out_dir, "manifest.json")
    write_manifest(
        manifest_path,
        cell_line,
        replicate_ids,
        experiment_accession,
        chrombpnet_commit,
        normalization_stats,
        aggregation_stats
    )

    return final_bw_path


# =============================================================================
# Main Pipeline (Local)
# =============================================================================

def process_cell_line_local(cell_line, bam_infos, chrom_sizes_path, output_dir, chrombpnet_commit=None, max_workers=4):
    """Process one cell line in local mode (from metadata.tsv)."""
    out_dir = os.path.join(output_dir, cell_line)
    os.makedirs(out_dir, exist_ok=True)

    log_path = os.path.join(out_dir, "pipeline.log")
    # Pass cell_line to setup_logging for prefixed messages
    setup_logging(log_path, cell_line=cell_line)

    logging.info("=" * 60)
    logging.info("Processing (local mode)")
    logging.info("=" * 60)
    logging.info(f"Started at: {datetime.now().isoformat()}")

    try:
        chrom_sizes = load_chrom_sizes_from_file(chrom_sizes_path)
        valid_chroms = set(chrom_sizes.keys())
        logging.info(f"Loaded {len(valid_chroms)} chromosomes from {chrom_sizes_path}")

        # Wrapped in try/finally so it is always cleaned up, even on failure.
        bigwig_dir = tempfile.mkdtemp(prefix="bigwigs_")
        try:
            replicate_bws = []
            replicate_ids = []

            # Pre-validate BAM existence before forking so warnings go to the log cleanly.
            valid_bam_infos = []
            for bam_info in bam_infos:
                if not os.path.exists(bam_info["bam_path"]):
                    logging.warning(f"BAM not found: {bam_info['bam_path']}, skipping")
                else:
                    logging.info(f"Queuing {bam_info['replicate_id']} ({bam_info['bam_path']})...")
                    valid_bam_infos.append(bam_info)

            # --- PARALLEL: convert each replicate BAM to BigWig concurrently ---
            worker_args = [
                (
                    bam_info["bam_path"],
                    os.path.join(bigwig_dir, bam_info["replicate_id"]),
                    valid_chroms,
                    chrom_sizes,
                )
                for bam_info in valid_bam_infos
            ]
            with ProcessPoolExecutor(max_workers=max_workers) as pool:
                bw_paths = list(pool.map(_bam_to_bigwig_worker, worker_args))

            for bam_info, bw_path in zip(valid_bam_infos, bw_paths):
                replicate_bws.append(bw_path)
                replicate_ids.append(bam_info["replicate_id"])

            if not replicate_bws:
                raise ValueError("No BigWigs generated")

            replicate_bws = sorted(replicate_bws)
            replicate_ids = sorted(replicate_ids)

            final_bw_path = _finalize_cell_line(
                replicate_bws, replicate_ids, chrom_sizes, out_dir,
                cell_line, "local", chrombpnet_commit
            )

        finally:
            # Always clean up temp bigwig dir, whether success or fail.
            shutil.rmtree(bigwig_dir, ignore_errors=True)

        logging.info("=" * 60)
        logging.info(f"SUCCESS: completed")
        logging.info(f"Output: {final_bw_path}")
        logging.info("=" * 60)

        return True

    except Exception as e:
        logging.error(f"FAILED: {e}")
        logging.exception("Traceback:")
        return False

# =============================================================================
# Main Pipeline (ENCODE)
# =============================================================================

def process_cell_line_encode(cell_line, skip_download=False, chrombpnet_commit=None, max_workers=4):
    """Process one cell line in ENCODE mode (query ENCODE API)."""
    out_dir = f"out/{cell_line}"
    os.makedirs(out_dir, exist_ok=True)

    log_path = os.path.join(out_dir, "pipeline.log")
    # Pass cell_line to setup_logging for prefixed messages
    setup_logging(log_path, cell_line=cell_line)

    logging.info("=" * 60)
    logging.info("Processing (ENCODE mode)")
    logging.info("=" * 60)
    logging.info(f"Started at: {datetime.now().isoformat()}")

    try:
        # Query ENCODE for BAM files
        file_infos = query_encode_for_cell_line(cell_line, max_api_workers=max_workers)

        if not file_infos:
            raise ValueError("No valid BAM files found")

        experiment_accession = file_infos[0]["experiment_accession"]
        bam_dir = os.path.join("bams", cell_line)
        os.makedirs(bam_dir, exist_ok=True)

        if not skip_download:
            logging.info("Downloading BAM files...")
            # --- PARALLEL: download BAM files concurrently (network I/O-bound) ---
            download_args = [(fi, bam_dir) for fi in file_infos]
            with ThreadPoolExecutor(max_workers=max_workers) as pool:
                list(pool.map(_download_bam_worker, download_args))
        else:
            logging.info("Skipping download (--skip-download)")

        # Convert BAMs to BigWigs
        logging.info("Converting BAMs to BigWigs...")
        valid_chroms = load_valid_chroms_from_fasta(GENOME_FA)
        chrom_sizes = load_chrom_sizes_from_file(CHROM_SIZES)
        logging.info(f"Loaded {len(valid_chroms)} chromosomes from reference")

        bigwig_dir = os.path.join("bigwigs", cell_line)
        os.makedirs(bigwig_dir, exist_ok=True)

        replicate_bws = []
        file_accessions = []

        # Pre-validate BAM existence before forking so warnings go to the log cleanly.
        valid_file_infos = []
        for file_info in file_infos:
            file_acc = file_info["file_accession"]
            bam_path = os.path.join(bam_dir, f"{file_acc}.bam")
            if not os.path.exists(bam_path):
                logging.warning(f"BAM not found: {bam_path}, skipping")
            else:
                logging.info(f"Queuing {file_acc}...")
                valid_file_infos.append((file_info, bam_path))

        # --- PARALLEL: convert each replicate BAM to BigWig concurrently ---
        worker_args = [
            (
                bam_path,
                os.path.join(bigwig_dir, fi["file_accession"]),
                valid_chroms,
                chrom_sizes,
            )
            for fi, bam_path in valid_file_infos
        ]
        with ProcessPoolExecutor(max_workers=max_workers) as pool:
            bw_paths = list(pool.map(_bam_to_bigwig_worker, worker_args))

        for (file_info, _), bw_path in zip(valid_file_infos, bw_paths):
            replicate_bws.append(bw_path)
            file_accessions.append(file_info["file_accession"])

        if not replicate_bws:
            raise ValueError("No BigWigs generated")

        replicate_bws = sorted(replicate_bws)
        file_accessions = sorted(file_accessions)

        final_bw_path = _finalize_cell_line(
            replicate_bws, file_accessions, chrom_sizes, out_dir,
            cell_line, experiment_accession, chrombpnet_commit
        )

        logging.info("=" * 60)
        logging.info("SUCCESS: completed")
        logging.info("Output: {final_bw_path}")
        logging.info("=" * 60)

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
        default=4,
        help="Workers per cell line: controls parallel BAM->BigWig conversion and "
             "concurrent ENCODE API fetches (default: 4)"
    )
    parser.add_argument(
        "--cell-line-workers",
        type=int,
        default=3,
        help="Number of cell lines to process in parallel (default: 3)"
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

    # Check required external tools are on PATH before doing any work.
    # bedtools is needed in all modes (BAM->BigWig conversion).
    # wget is only needed in ENCODE mode when not using --skip-download.
    tools_to_check = ["bedtools"]
    if not args.metadata and not args.skip_download:
        tools_to_check.append("wget")
    missing = [tool for tool in tools_to_check if shutil.which(tool) is None]
    if missing:
        print(f"ERROR: Required external tool(s) not found on PATH: {', '.join(missing)}")
        print("Please install them and ensure they are accessible before running.")
        sys.exit(1)

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

        # --- PARALLEL: each cell line writes to its own out/<cell_line>/ and
        worker_args = [
            (cell_line, bam_infos, args.chrom_sizes, args.outdir,
             chrombpnet_commit, args.threads)
            for cell_line, bam_infos in cell_lines_data.items()
        ]
        with ProcessPoolExecutor(max_workers=args.cell_line_workers) as pool:
            results = list(pool.map(_process_cell_line_local_worker, worker_args))
        success = all(results)
        
    else:
        # ENCODE mode
        print("Running in ENCODE mode")
        
        if not os.path.exists(GENOME_FA):
            print(f"ERROR: Reference genome not found: {GENOME_FA}")
            sys.exit(1)

        if not os.path.exists(CHROM_SIZES):
            print(f"ERROR: Chrom sizes not found: {CHROM_SIZES}")
            sys.exit(1)

        # --- PARALLEL: each cell line writes to out/<cell_line>/, bams/<cell_line>/,
        # and bigwigs/<cell_line>/.
        worker_args = [
            (cell_line, args.skip_download, chrombpnet_commit, args.threads)
            for cell_line in args.cell_lines
        ]
        with ProcessPoolExecutor(max_workers=args.cell_line_workers) as pool:
            results = list(pool.map(_process_cell_line_encode_worker, worker_args))
        success = all(results)

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
