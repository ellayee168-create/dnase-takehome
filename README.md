# DNase Preprocessing Take-Home Assignment

## Overview

Implement the DNase preprocessing pipeline following the Alphagenome methodology to process DNase-seq data from ENCODE into normalized BigWig tracks suitable for machine learning.

## Objective

You will implement the DNase preprocessing step described in the Alphagenome paper for 3 cell lines (GM12878, HeLa-S3, and SK-N-SH). At a high level, you should:

1. **Download** BAM files from ENCODE
2. **Convert** BAM → BigWig using ChromBPNet's `reads_to_bigwig.py`
3. **Aggregate** replicates by cell line (mean)
4. **Normalize** to 10^8 total counts
5. **Validate** output and generate provenance manifest

Ensure that you apply your quality control similarly to how Alphagenome did.

## Target Cell Lines

Download DNase-seq aligned BAM files (hg38) from ENCODE for:

| Cell Line | Description |
|-----------|-------------|
| GM12878   | Lymphoblastoid cell line |
| HeLa-S3   | Cervical cancer cell line |
| SK-N-SH   | Neuroblastoma cell line |

For each experiment, find the "alignments" files in BAM format with **hg38** assembly. Download all biological replicates.

## Setup

```bash
# Create conda environment
conda env create -f environment.yml
conda activate dnase-pipeline

# Download hg38 chromosome sizes
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes
```

## CLI Requirements

One command that runs end-to-end:

```bash
python run_dnase_tracks.py \
    --metadata metadata.tsv \
    --chrom-sizes reference.chrom.sizes \
    --outdir out/ \
    --threads 4
```

Requirements:
- Exits with code 0 on success, non-zero on failure
- Writes logs to `out/<cell_line>/pipeline.log`

## Expected Output

```
out/
├── GM12878/
│   ├── dnase_avg_norm100M.bw    # Final normalized BigWig
│   ├── manifest.json             # Provenance and metadata
│   └── pipeline.log              # Processing log
├── HeLa-S3/
│   └── ...
└── SK-N-SH/
    └── ...
```

### Manifest Format

```json
{
  "cell_line": "GM12878",
  "replicates": [
    {"replicate_id": "ENCFF..."}
  ],
  "chrombpnet": {
    "repo_url": "https://github.com/kundajelab/chrombpnet",
    "commit_sha": "..."
  },
  "normalization": {
    "pre_total": 12345678.0,
    "post_total": 100000000.0,
    "scaling_factor": 8.1
  },
  "aggregation": {
    "non_zero_average_mean": 0.123,
    "n_replicates": 2
  }
}
```

## Requirements

### Determinism
- We will run your pipeline at least twice on the same inputs
- The final file must be identical between runs (same checksum)
- If you use any nondeterministic tools (multithreading, temp files, etc.), you are responsible for ensuring that your final results are consistent

### Methods
- Your pipeline must follow exactly what the paper describes for their DNase/ATAC preprocessing
- You may use any ChromBPNet version, but you must record which commit you used in manifest.json

### Validation (done by your pipeline)
- BigWig validity (file exists and is readable)
- Chromosome names match
- Total signal sum in the final BigWig is 10^8 AUC (within reasonable tolerance)
- Replicate aggregation: Record which sample IDs are aggregated for each cell line

### Failure Behavior
On failure:
- Exit non-zero
- Write a helpful error message to the appropriate pipeline.log file

## Testing

```bash
# Generate synthetic test data
python scripts/generate_test_data.py --output-dir test_data --small

# Run public tests (uses synthetic data)
pytest tests/public/ -v
```

## Grading Criteria

We grade primarily on operational quality and your ability to explain your choices:
- Should be able to rerun your code and get the same outputs with no debugging required

Should have:
- Deterministic reruns
- Conda env reproducibility
- Clear CLI + logs
- Strict and consistent output naming with the given keys
- Valid BigWigs that are correctly normalized

You will have visible test cases that you can run to verify parts of your output:
- All visible test cases will use your manifest.json file

We will manually verify parts such as the validation steps and the failure behavior to ensure they are implemented correctly.

## Submission

- Edit the given GitHub repo containing:
  - `run_dnase_tracks.py`
  - `environment.yml`
  - Any other relevant code
- Zip file containing your final DNase tracks

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
