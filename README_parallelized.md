# DNase-seq Preprocessing Pipeline (Parallelized)

Replicates the DNase preprocessing methodology from the [AlphaGenome paper](https://doi.org/10.1038/s41586-025-10014-0) (Nature 2025) for generating normalized chromatin accessibility tracks.

This parallelized version implements multi-threaded downloads, concurrent API fetches, and multi-process BAM→BigWig conversion to increase efficiency.

## Features

- **Two modes**: ENCODE (query API) or Local (process pre-downloaded BAMs)
- **AlphaGenome QC filters**: FRiP > 10%, read length >= 36 nt
- **ChromBPNet methodology**: DNase shifts (0/+1), base-resolution BigWigs
- **100M normalization**: Rescales total signal to 10^8 counts
- **Parallelization**: Multi-threaded I/O, multi-process computation

## Parallelization Overview

| Operation | Parallelization | Executor | Reasoning |
|-----------|-----------------|----------|-----|
| ENCODE API fetches | Concurrent | `ThreadPoolExecutor` | Network I/O-bound |
| BAM downloads | Concurrent | `ThreadPoolExecutor` | Network I/O-bound |
| BAM → BigWig conversion | Parallel | `ProcessPoolExecutor` | CPU-bound |
| Cell line processing | Parallel | `ProcessPoolExecutor` | Independent workloads |

## Quick Start

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate dnase-pipeline

# 2. Run pipeline (ENCODE mode - preferred for evaluation)
python run_dnase_tracks_parallel.py --cell-lines GM12878

# Or run with local BAM files (local mode)
python run_dnase_tracks_parallel.py --metadata metadata.tsv --chrom-sizes ref.chrom.sizes --outdir out
```

## Installation

```bash
# Create conda environment
conda env create -f environment.yml
conda activate dnase-pipeline

# For ENCODE mode: download reference files
mkdir -p reference
wget -O reference/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip reference/hg38.fa.gz
wget -O reference/hg38.chrom.sizes https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# Clone ChromBPNet for commit tracking
git clone https://github.com/kundajelab/chrombpnet.git
```

## Usage

### ENCODE Mode (Query ENCODE API, preferred for evaluation)

Downloads and processes DNase-seq data directly from ENCODE:

```bash
# Process default cell lines (GM12878, HeLa-S3, SK-N-SH) with 3 cell lines in parallel
python run_dnase_tracks_parallel.py --cell-lines GM12878 HeLa-S3 SK-N-SH --cell-line-workers 3

# Process specific cell lines with custom parallelism
python run_dnase_tracks_parallel.py --cell-lines GM12878 HeLa-S3 --threads 8 --cell-line-workers 2

# Skip download (use existing BAMs in bams/<cell_line>/)
python run_dnase_tracks_parallel.py --skip-download
```

### Local Mode (Pre-downloaded BAMs)

Process local BAM files specified in a metadata TSV:

```bash
python run_dnase_tracks_parallel.py \
    --metadata metadata.tsv \
    --chrom-sizes reference.chrom.sizes \
    --outdir output_dir \
    --threads 4 \
    --cell-line-workers 2
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
| `--threads` | Workers per cell line (BAM→BigWig, API fetches, downloads) | 4 |
| `--cell-line-workers` | Number of cell lines to process in parallel | 3 |
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

**Note**: This parallelized version stores intermediate files in cell-line-specific directories to prevent race conditions:
- `bams/<cell_line>/` — Downloaded BAM files
- `bigwigs/<cell_line>/` — Intermediate BigWig files

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

### Overview

1. **Get BAM files** — Query ENCODE API or read from metadata TSV
2. **QC Filtering** — Apply AlphaGenome quality filters (ENCODE mode only)
3. **Experiment Selection** — Select one experiment per cell line
4. **BAM to BigWig** — Convert using ChromBPNet methodology with DNase shifts
5. **Average Replicates** — Combine biological replicates
6. **Normalize** — Scale to 100M total counts
7. **Validate** — Verify BigWig integrity and signal sum
8. **Generate Manifest** — Record metadata and statistics

### Quality Control Filters (ENCODE Mode)

Based on [AlphaGenome supplementary methods](https://doi.org/10.1038/s41586-025-10014-0):

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

### Determinism

The parallelized version produces identical results to the serial version:

1. **API fetches**: Results stored in dict keyed by accession, then iterated in original API order
2. **Downloads**: Order irrelevant; files validated by existence and sorted before processing
3. **BAM→BigWig**: `pool.map()` preserves input order; results sorted before averaging
4. **Cell lines**: Each writes to isolated directories; no shared mutable state

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
