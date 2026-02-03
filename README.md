# DNase Preprocessing Take-Home Assignment

## Overview

Implement the DNase preprocessing pipeline following the Alphagenome methodology to process DNase-seq data from ENCODE into normalized BigWig tracks suitable for machine learning.

## Objective

Process DNase-seq data for 3 cell lines through the following pipeline:

1. **Download** BAM files from ENCODE
2. **Convert** BAM → BigWig using ChromBPNet's `reads_to_bigwig.py`
3. **Aggregate** replicates by cell line (mean)
4. **Normalize** to 10^8 total counts
5. **Validate** output and generate provenance manifest

## Target Cell Lines

Download DNase-seq aligned BAM files (hg38) from ENCODE for:

| Cell Line | ENCODE Experiment |
|-----------|-------------------|
| MCF-7     | [ENCSR000EJD](https://www.encodeproject.org/experiments/ENCSR000EJD/) |
| K562      | [ENCSR000EMT](https://www.encodeproject.org/experiments/ENCSR000EMT/) |
| SK-N-SH   | [ENCSR000ENH](https://www.encodeproject.org/experiments/ENCSR000ENH/) |

For each experiment, find the "alignments" files in BAM format with **hg38** assembly. Download all biological replicates.

## Setup

```bash
# Create conda environment
conda env create -f environment.yml
conda activate dnase-pipeline

# Download hg38 chromosome sizes
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
```

## Your Task

### 1. Implement `run_dnase_tracks.py`

Complete the pipeline script to:

- Download BAM files from ENCODE (or accept local paths for testing)
- Validate inputs (BAM exists, has index, chromosome compatibility)
- Convert BAM → BigWig using ChromBPNet
- Aggregate replicates by averaging
- Normalize to 10^8 total counts
- Write manifest.json with full provenance

### 2. Implement `quality_checks.py`

Complete the quality checking script to validate:

- BigWig file integrity
- Normalization accuracy (within 1% tolerance)
- Manifest completeness
- Cross-cell-line consistency

## Expected Output

```
out/
├── MCF-7/
│   ├── dnase_avg_norm100M.bw    # Final normalized BigWig
│   ├── manifest.json             # Provenance and metadata
│   └── pipeline.log              # Processing log
├── K562/
│   └── ...
└── SK-N-SH/
    └── ...
```

### Manifest Format

```json
{
  "cell_line": "MCF-7",
  "replicates": [
    {"replicate_id": "ENCFF...", "bam_path": "..."}
  ],
  "chrombpnet": {
    "repo_url": "https://github.com/kundajelab/chrombpnet",
    "commit_sha": "..."
  },
  "commands": ["python -m chrombpnet.helpers..."],
  "normalization": {
    "target_total": 100000000,
    "pre_total": 12345678.0,
    "post_total": 100000000.0,
    "scaling_factor": 8.1,
    "tolerance": 0.01
  },
  "aggregation": {
    "method": "mean",
    "n_replicates": 2
  }
}
```

## Testing

```bash
# Generate synthetic test data
python scripts/generate_test_data.py --output-dir test_data --small

# Run public tests (uses synthetic data)
pytest tests/public/ -v

# Run your quality checks
python quality_checks.py --output-dir out/
```

## Evaluation Criteria

1. **Correctness**: Output matches expected format and normalization
2. **Reproducibility**: Deterministic output (same input → same output)
3. **Error Handling**: Graceful handling of missing/corrupted files
4. **Code Quality**: Clean, readable, well-documented code
5. **Provenance**: Complete manifest for reproducibility

## Requirements

- Python 3.10+
- ChromBPNet
- pysam, pyBigWig
- samtools (for BAM indexing if needed)

See `environment.yml` for complete dependencies.

## Hints

- Use `chrombpnet.helpers.preprocessing.reads_to_bigwig` for BAM → BigWig conversion
- Process chromosomes in sorted order for deterministic output
- The ENCODE portal provides direct download URLs for BAM files
- Test with synthetic data first before downloading large ENCODE files
