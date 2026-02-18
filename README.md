# DNase-seq Preprocessing Pipeline

Replicates the DNase preprocessing methodology from the [AlphaGenome paper](https://doi.org/10.1038/s41586-025-10014-0) (Nature 2025) for generating normalized chromatin accessibility tracks.

## Features

- **Two modes**: ENCODE mode (query API) or Local mode (process pre-downloaded BAMs)
- **AlphaGenome QC filters**: FRiP > 10%, read length >= 36 nt
- **ChromBPNet methodology**: DNase shifts (0/+1), base-resolution BigWigs
- **100M normalization**: Rescales total signal to 10^8 counts

## Quick Start

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate dnase-pipeline

# 2. Run pipeline (ENCODE mode - downloads data automatically)
python run_dnase_tracks.py --cell-lines GM12878

# Or run with local BAM files
python run_dnase_tracks.py --metadata metadata.tsv --chrom-sizes ref.chrom.sizes --outdir out
```

## Installation

```bash
# Create conda environment
conda env create -f environment.yml
conda activate dnase-pipeline

# For ENCODE mode only: download reference files
mkdir -p reference
wget -O reference/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip reference/hg38.fa.gz
wget -O reference/hg38.chrom.sizes https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# Clone ChromBPNet for commit tracking
git clone https://github.com/kundajelab/chrombpnet.git
```

## Usage

### ENCODE Mode (Query ENCODE API)

Downloads and processes DNase-seq data directly from ENCODE:

```bash
# Process default cell lines (GM12878, HeLa-S3, SK-N-SH)
python run_dnase_tracks.py

# Process specific cell lines
python run_dnase_tracks.py --cell-lines GM12878 HeLa-S3

# Skip download (use existing BAMs in bams/)
python run_dnase_tracks.py --skip-download
```

### Local Mode (Pre-downloaded BAMs)

Process local BAM files specified in a metadata TSV:

```bash
python run_dnase_tracks.py \
    --metadata metadata.tsv \
    --chrom-sizes reference.chrom.sizes \
    --outdir output_dir \
    --threads 4
```

**Metadata TSV format:**
```
File accession	Biosample term name	Assay	bam_path
ENCFF001ABC	GM12878	DNase-seq	/path/to/rep1.bam
ENCFF002DEF	GM12878	DNase-seq	/path/to/rep2.bam
```

| Column | Description |
|--------|-------------|
| `File accession` | Replicate identifier |
| `Biosample term name` | Cell line name (used for grouping) |
| `Assay` | Assay type (informational) |
| `bam_path` | Absolute path to BAM file |

### CLI Reference

| Flag | Description | Default |
|------|-------------|---------|
| `--metadata` | Path to metadata TSV (enables local mode) | None |
| `--chrom-sizes` | Chromosome sizes file (required for local mode) | None |
| `--outdir` | Output directory | `out` |
| `--threads` | Number of threads | 1 |
| `--cell-lines` | Cell lines to process (ENCODE mode) | GM12878, HeLa-S3, SK-N-SH |
| `--skip-download` | Skip BAM download (ENCODE mode) | False |
| `--chrombpnet-commit` | ChromBPNet commit SHA | Auto-detected |

## Output

```
out/
├── GM12878/
│   ├── dnase_avg_norm100M.bw   # Final normalized BigWig
│   ├── manifest.json           # Metadata and statistics
│   └── pipeline.log            # Execution log
├── HeLa-S3/
│   └── ...
└── SK-N-SH/
    └── ...
```

### manifest.json

```json
{
  "cell_line": "GM12878",
  "experiment": "ENCSR000EMT",
  "replicates": ["ENCFF001ABC", "ENCFF002DEF"],
  "chrombpnet": {
    "repo_url": "https://github.com/kundajelab/chrombpnet",
    "commit_sha": "ece97c93..."
  },
  "normalization": {
    "pre_total": 1.5e8,
    "post_total": 1e8,
    "scaling_factor": 0.667
  },
  "aggregation": {
    "non_zero_average_mean": 1.97,
    "n_replicates": 2
  }
}
```

## Methodology

### Pipeline Steps

1. **Get BAM files** — Query ENCODE API or read from metadata TSV
2. **QC Filtering** — Apply AlphaGenome quality filters (ENCODE mode only)
3. **Experiment Selection** — Select one experiment per cell line
4. **BAM to BigWig** — Convert using ChromBPNet methodology with DNase shifts
5. **Average Replicates** — Combine biological replicates
6. **Normalize** — Scale to 100M total counts
7. **Validate** — Verify BigWig integrity and signal sum
8. **Generate Manifest** — Record metadata and statistics

### Quality Control Filters (ENCODE Mode)

From the [AlphaGenome supplementary methods](https://doi.org/10.1038/s41586-025-10014-0):

| Filter | Threshold | ENCODE Field | Rationale |
|--------|-----------|--------------|-----------|
| FRiP | > 10% | `spot1_score` | Signal-to-noise quality |
| Read length | >= 36 nt | FASTQ `read_length` | All FASTQs in experiment must pass |

### Experiment Selection (ENCODE Mode)

When multiple experiments pass QC:
1. Prefer paired-end over single-end
2. Prefer newest by `date_released` if tied

### DNase Shift Parameters

From ChromBPNet methodology:
- Plus strand: +0 shift
- Minus strand: +1 shift

### Normalization

All tracks normalized to 10^8 total counts (100M AUC).

## Dependencies

See `environment.yml` for full list. Key packages:

| Package | Purpose |
|---------|---------|
| `pyBigWig` | BigWig I/O |
| `pyfaidx` | FASTA indexing |
| `pysam` | BAM processing |
| `bedtools` | BAM/BED operations |
| `samtools` | BAM indexing |
| `requests` | ENCODE API queries |
| `numpy` | Numerical operations |

## References

- [AlphaGenome paper](https://doi.org/10.1038/s41586-025-10014-0) — DeepMind, Nature 2025
- [ChromBPNet](https://github.com/kundajelab/chrombpnet) — Kundaje Lab
- [ENCODE Portal](https://www.encodeproject.org/) — Data source
