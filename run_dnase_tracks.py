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

Target cell lines and their ENCODE experiment accessions:
- MCF-7: ENCSR000EJD
- K562: ENCSR000EMT
- SK-N-SH: ENCSR000ENH

Find BAM file URLs at: https://www.encodeproject.org/experiments/{ACCESSION}/
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

Author: [Your Name]
"""

# TODO: Implement the pipeline
