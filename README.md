# DNase Preprocessing Take-Home Assignment

## Overview

Implement the DNase preprocessing pipeline following the Alphagenome methodology to process DNase-seq data from raw BAM files to normalized BigWig tracks.

## Objective

Process DNase-seq data for 3 cell lines (MCF-7, K562, SK-N-SH) through the following pipeline:
1. Convert BAM → BigWig using ChromBPNet's `reads_to_bigwig.py`
2. Aggregate replicates by cell line (mean)
3. Normalize to 10^8 total counts

## Setup

```bash
# Create conda environment
conda env create -f environment.yml
conda activate dnase-pipeline
```

## Usage

```bash
python run_dnase_tracks.py \
    --metadata metadata.tsv \
    --chrom-sizes reference.chrom.sizes \
    --outdir out/ \
    --threads 4
```

## Output Structure

```
out/
├── MCF-7/
│   ├── dnase_avg_norm100M.bw
│   ├── manifest.json
│   └── pipeline.log
├── K562/
│   └── ...
└── SK-N-SH/
    └── ...
```

## Testing

```bash
# Run public tests
pytest tests/public/ -v

# Generate test data
python scripts/generate_test_data.py --output-dir test_data --small
```

## Requirements

- Python 3.10+
- ChromBPNet
- pysam, pyBigWig
- samtools, bedtools

See `environment.yml` for complete dependencies.
