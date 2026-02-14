# DNase-seq Preprocessing Pipeline

Implements the DNase preprocessing methodology described in the [AlphaGenome paper](https://doi.org/10.1038/s41586-025-10014-0) for 3 cell lines: GM12878, HeLa-S3, and SK-N-SH.

## Quick Start

```bash
# 1. Create environment
conda env create -f environment.yml
conda activate dnase-pipeline

# 2. Clone ChromBPNet (for commit tracking)
git clone https://github.com/kundajelab/chrombpnet.git

# 3. Download reference files
mkdir -p reference
wget -O reference/hg38.fa.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip reference/hg38.fa.gz
wget -O reference/hg38.chrom.sizes https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# 4. Run pipeline
python run_dnase_tracks.py
```

## Usage

```bash
# Full pipeline (downloads BAMs from ENCODE)
python run_dnase_tracks.py

# Skip download (use existing BAMs in bams/)
python run_dnase_tracks.py --skip-download

# Process specific cell lines
python run_dnase_tracks.py --cell-lines GM12878 HeLa-S3

# Specify ChromBPNet commit manually
python run_dnase_tracks.py --chrombpnet-commit abc123
```

### CLI Flags

| Flag | Description |
|------|-------------|
| `--cell-lines` | Cell lines to process (default: GM12878, HeLa-S3, SK-N-SH) |
| `--skip-download` | Skip BAM download, use existing files in `bams/` |
| `--chrombpnet-commit` | ChromBPNet commit SHA (auto-detected from `chrombpnet/` if not specified) |

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

## Methodology

### Pipeline Steps

1. **Query ENCODE** - Find DNase-seq BAM files for each cell line
2. **QC Filtering** - Apply AlphaGenome quality filters (see below)
3. **BAM → BigWig** - Convert using ChromBPNet methodology with DNase shifts
4. **Average Replicates** - Combine all replicates per cell line
5. **Normalize** - Scale to 100M total counts
6. **Validate** - Verify BigWig integrity and signal sum
7. **Generate Manifest** - Record all metadata and statistics

### Quality Control Filters (AlphaGenome)

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| FRiP (spot1_score) | > 10% | Signal-to-noise quality |
| ENCODE audit | No ERROR | Data quality |

### DNase Shift Parameters

From AlphaGenome Paper: "Shifts of 0/+1 for DNase-seq were applied"
- Plus strand: +0 shift
- Minus strand: +1 shift

### Normalization

All tracks normalized to **10^8 total counts** (100M AUC) for cross-sample comparisons.

## Implementation Details

### Signal Preservation

The pipeline uses `pyBigWig.addEntries()` directly instead of the external `bedGraphToBigWig` tool to avoid signal loss from float-to-int compression during bedGraph conversion. ~12% signal loss was observed when using bedGraphToBigWig directly.

### `non_zero_average_mean` Calculation

Computed from the **final normalized BigWig** as the mean of all non-zero interval values:
```python
for chrom in bw.chroms():
    for start, end, value in bw.intervals(chrom):
        if value != 0:
            non_zero_values.append(value)
mean = np.mean(non_zero_values)
```

## Validation

The pipeline validates each output BigWig:

1. **File readable** - Opens successfully with pyBigWig
2. **Chromosomes valid** - All chromosome names in reference
3. **Signal sum** - Total = 10^8 ± 0.1%

## Failure Behavior

On failure:
- Exit code: non-zero
- Error message written to `out/<cell_line>/pipeline.log`
- Traceback included for debugging

## Dependencies

See `environment.yml`. Key packages:
- `pyBigWig` - BigWig I/O
- `pyfaidx` - FASTA indexing
- `bedtools` - BAM/BED operations
- `samtools` - BAM processing
- `requests` - ENCODE API queries

## References

- [AlphaGenome paper](https://doi.org/10.1038/s41586-025-10014-0) - DeepMind, Nature 2025
- [ChromBPNet](https://github.com/kundajelab/chrombpnet) - Kundaje Lab
- [ENCODE Portal](https://www.encodeproject.org/) - Data source
