#!/usr/bin/env python
"""
Query ENCODE for DNase-seq BAM files for the 3 cell lines.
Applies QC filters per AlphaGenome methodology.

Output: metadata.tsv
"""

import requests
import json
import time

CELL_LINES = ["GM12878", "HeLa-S3", "SK-N-SH"]

def get_json(url):
    """Get JSON from ENCODE API with rate limiting."""
    response = requests.get(url, headers={"Accept": "application/json"})
    time.sleep(0.1)  # Be nice to the API
    return response.json()

def trace_read_length(file_acc, depth=0):
    """
    Recursively trace derived_from to find read_length from FASTQ files.
    Returns list of read_lengths found.
    """
    if depth > 5:
        return []
    
    detail = get_json(f"https://www.encodeproject.org/files/{file_acc}/?format=json")
    file_format = detail.get("file_format")
    read_length = detail.get("read_length")
    
    # If this is a FASTQ with read_length, return it
    if file_format == "fastq" and read_length:
        return [read_length]
    
    # Otherwise, trace derived_from
    derived_from = detail.get("derived_from", [])
    read_lengths = []
    
    for df in derived_from:
        if isinstance(df, str):
            df_acc = df.split("/")[-2]
        else:
            continue
        read_lengths.extend(trace_read_length(df_acc, depth + 1))
    
    return read_lengths

def get_spot1_score(file_detail):
    """
    Get spot1_score from HotspotQualityMetric.
    spot1_score is the DNase-seq equivalent of FRiP.
    """
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

def query_cell_line(cell_line):
    """Query ENCODE for DNase BAM files for one cell line."""
    print(f"\n{'='*60}")
    print(f"Querying {cell_line}...")
    print(f"{'='*60}")
    
    # Search directly for BAM files
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
    print(f"Found {len(files)} BAM files")
    
    valid_files = []
    
    for f in files:
        file_acc = f.get("accession")
        print(f"\nChecking {file_acc}...")
        
        # Get full file details
        file_detail = get_json(f"https://www.encodeproject.org/files/{file_acc}/?format=json")
        
        # Get experiment accession
        dataset = file_detail.get("dataset", "")
        exp_acc = dataset.split("/")[-2] if dataset else None
        
        if not exp_acc:
            print(f"  SKIP: No experiment accession")
            continue
        
        # FILTER 1: No ERROR audits at experiment level
        if has_error_audit(exp_acc):
            print(f"  SKIP: Experiment {exp_acc} has ERROR audit")
            continue
        
        # FILTER 2: Get read_length from derived FASTQs (must be >= 36)
        read_lengths = trace_read_length(file_acc)
        if not read_lengths:
            print(f"  SKIP: Could not determine read_length")
            continue
        
        min_read_length = min(read_lengths)
        if min_read_length < 36:
            print(f"  SKIP: read_length={min_read_length} < 36")
            continue
        
        # FILTER 3: spot1_score > 0.10 (DNase equivalent of FRiP > 10%)
        spot1_score = get_spot1_score(file_detail)
        if spot1_score is None:
            print(f"  SKIP: No spot1_score found")
            continue
        
        if spot1_score <= 0.10:
            print(f"  SKIP: spot1_score={spot1_score:.3f} <= 0.10")
            continue
        
        # Get download URL
        href = file_detail.get("href")
        if href:
            download_url = f"https://www.encodeproject.org{href}"
        else:
            download_url = f"https://www.encodeproject.org/files/{file_acc}/@@download/{file_acc}.bam"
        
        # Get biological replicate info
        bio_reps = file_detail.get("biological_replicates", [])
        
        valid_files.append({
            "cell_line": cell_line,
            "experiment_accession": exp_acc,
            "file_accession": file_acc,
            "bam_url": download_url,
            "spot1_score": spot1_score,
            "read_length": min_read_length,
            "biological_replicates": bio_reps
        })
        
        print(f"  VALID: spot1_score={spot1_score:.3f}, read_length={min_read_length}")
    
    return valid_files

def main():
    all_files = []
    
    for cell_line in CELL_LINES:
        files = query_cell_line(cell_line)
        all_files.extend(files)
        print(f"\n{cell_line}: {len(files)} valid files")
    
    # Write to TSV
    with open("metadata.tsv", "w") as f:
        f.write("cell_line\texperiment_accession\tfile_accession\tbam_url\tspot1_score\tread_length\tbiological_replicates\n")
        for file_info in all_files:
            bio_reps = ",".join(map(str, file_info["biological_replicates"]))
            f.write(f"{file_info['cell_line']}\t{file_info['experiment_accession']}\t{file_info['file_accession']}\t{file_info['bam_url']}\t{file_info['spot1_score']}\t{file_info['read_length']}\t{bio_reps}\n")
    
    print(f"\n{'='*60}")
    print(f"Wrote {len(all_files)} files to metadata.tsv")
    
    # Summary by cell line
    for cell_line in CELL_LINES:
        count = len([f for f in all_files if f["cell_line"] == cell_line])
        print(f"  {cell_line}: {count} files")
    
    # Also save as JSON
    with open("metadata.json", "w") as f:
        json.dump(all_files, f, indent=2)

if __name__ == "__main__":
    main()