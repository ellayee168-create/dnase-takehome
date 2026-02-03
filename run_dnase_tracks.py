#!/usr/bin/env python3
"""
DNase Preprocessing Pipeline - Take-Home Assignment

Implement the DNase preprocessing pipeline following the Alphagenome methodology.

Your task:
1. Download DNase-seq BAM files from ENCODE for the specified cell lines
2. Convert BAM -> BigWig using ChromBPNet's reads_to_bigwig.py
3. Aggregate replicates by averaging per cell line
4. Normalize to 10^8 total counts
5. Output a manifest.json with provenance information

Target cell lines:
- GM12878
- HeLa-S3
- SK-N-SH

CLI usage:
    python run_dnase_tracks.py \
        --metadata metadata.tsv \
        --chrom-sizes reference.chrom.sizes \
        --outdir out/ \
        --threads 4

Expected output structure:
    out/
    ├── GM12878/
    │   ├── dnase_avg_norm100M.bw
    │   ├── manifest.json
    │   └── pipeline.log
    ├── HeLa-S3/
    │   └── ...
    └── SK-N-SH/
        └── ...

manifest.json must contain:
- cell_line: Name of the cell line
- replicates: List of {replicate_id, bam_path}
- chrombpnet: {repo_url, commit_sha}
- commands: List of commands run
- normalization: {target_total, pre_total, post_total, scaling_factor, tolerance}
- aggregation: {method, non_zero_average_mean, n_replicates}

Author: [Your Name]
"""

# TODO: Implement the pipeline
