#!/usr/bin/env python3
"""
Private Test Suite - Advanced Functionality Tests

These tests verify advanced functionality, edge cases, and error handling.
These are NOT provided to students and are used for grading.
"""

import hashlib
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest
import pysam
import pyBigWig

sys.path.insert(0, str(Path(__file__).parent.parent.parent))


# ============================================================================
# Helper Functions
# ============================================================================

def compute_file_hash(filepath: str) -> str:
    """Compute SHA256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(filepath, 'rb') as f:
        for block in iter(lambda: f.read(65536), b''):
            sha256.update(block)
    return sha256.hexdigest()


def compute_bigwig_total(bw_path: str) -> float:
    """Compute total signal from BigWig."""
    total = 0.0
    with pyBigWig.open(bw_path) as bw:
        for chrom in bw.chroms():
            intervals = bw.intervals(chrom)
            if intervals:
                for start, end, value in intervals:
                    total += value * (end - start)
    return total


def create_test_bam(output_path: str, chrom_sizes: dict, n_reads: int, seed: int):
    """Create a test BAM file."""
    import random
    random.seed(seed)
    np.random.seed(seed)

    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'LN': size, 'SN': chrom} for chrom, size in sorted(chrom_sizes.items())]
    }

    with pysam.AlignmentFile(output_path, "wb", header=header) as outf:
        read_id = 0
        for chrom, size in sorted(chrom_sizes.items()):
            for _ in range(n_reads // len(chrom_sizes)):
                pos = random.randint(0, max(0, size - 51))
                a = pysam.AlignedSegment()
                a.query_name = f"read_{read_id}"
                a.query_sequence = "A" * 50
                a.flag = 0
                a.reference_id = outf.get_tid(chrom)
                a.reference_start = pos
                a.mapping_quality = 60
                a.cigar = [(0, 50)]
                a.query_qualities = pysam.qualitystring_to_array("I" * 50)
                outf.write(a)
                read_id += 1

    sorted_path = output_path.replace('.bam', '.sorted.bam')
    pysam.sort("-o", sorted_path, output_path)
    os.rename(sorted_path, output_path)
    pysam.index(output_path)


# ============================================================================
# Test: Determinism
# ============================================================================

class TestDeterminism:
    """Tests for deterministic output."""

    @pytest.fixture(scope="class")
    def determinism_test_data(self, tmp_path_factory):
        """Create test data for determinism tests."""
        test_dir = tmp_path_factory.mktemp("determinism")

        chrom_sizes = {"chr1": 50000, "chr2": 40000}
        chrom_sizes_path = test_dir / "reference.chrom.sizes"
        with open(chrom_sizes_path, 'w') as f:
            for chrom, size in chrom_sizes.items():
                f.write(f"{chrom}\t{size}\n")

        bam_dir = test_dir / "bams"
        bam_dir.mkdir()

        create_test_bam(str(bam_dir / "rep1.bam"), chrom_sizes, 500, seed=42)
        create_test_bam(str(bam_dir / "rep2.bam"), chrom_sizes, 500, seed=43)

        metadata_path = test_dir / "metadata.tsv"
        with open(metadata_path, 'w') as f:
            f.write("File accession\tBiosample term name\tAssay\tbam_path\n")
            f.write(f"rep1\tDetTest\tDNase-seq\t{bam_dir}/rep1.bam\n")
            f.write(f"rep2\tDetTest\tDNase-seq\t{bam_dir}/rep2.bam\n")

        return test_dir

    def test_deterministic_reruns(self, determinism_test_data, tmp_path_factory):
        """Pipeline should produce identical output on reruns."""
        test_dir = determinism_test_data

        hashes = []
        for run in range(2):
            output_dir = tmp_path_factory.mktemp(f"output_run{run}")

            cmd = [
                "python", "run_dnase_tracks.py",
                "--metadata", str(test_dir / "metadata.tsv"),
                "--chrom-sizes", str(test_dir / "reference.chrom.sizes"),
                "--outdir", str(output_dir),
                "--threads", "1",
            ]

            result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
            assert result.returncode == 0, f"Run {run} failed: {result.stderr}"

            bw_path = output_dir / "DetTest" / "dnase_avg_norm100M.bw"
            file_hash = compute_file_hash(str(bw_path))
            hashes.append(file_hash)

        assert hashes[0] == hashes[1], "BigWig files differ between runs (non-deterministic)"


# ============================================================================
# Test: Error Handling
# ============================================================================

class TestErrorHandling:
    """Tests for proper error handling."""

    def test_missing_bam_file(self, tmp_path):
        """Pipeline should fail gracefully with missing BAM."""
        test_dir = tmp_path / "test"
        test_dir.mkdir()

        chrom_sizes_path = test_dir / "reference.chrom.sizes"
        with open(chrom_sizes_path, 'w') as f:
            f.write("chr1\t50000\n")

        metadata_path = test_dir / "metadata.tsv"
        with open(metadata_path, 'w') as f:
            f.write("File accession\tBiosample term name\tAssay\tbam_path\n")
            f.write("rep1\tTest\tDNase-seq\t/nonexistent/file.bam\n")

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(metadata_path),
            "--chrom-sizes", str(chrom_sizes_path),
            "--outdir", str(tmp_path / "output"),
        ]

        result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
        assert result.returncode != 0, "Pipeline should fail with missing BAM"

    def test_missing_bam_index(self, tmp_path):
        """Pipeline should fail gracefully with missing BAM index."""
        test_dir = tmp_path / "test"
        test_dir.mkdir()

        chrom_sizes = {"chr1": 50000}
        chrom_sizes_path = test_dir / "reference.chrom.sizes"
        with open(chrom_sizes_path, 'w') as f:
            f.write("chr1\t50000\n")

        bam_dir = test_dir / "bams"
        bam_dir.mkdir()
        bam_path = bam_dir / "rep1.bam"
        create_test_bam(str(bam_path), chrom_sizes, 100, seed=42)

        # Remove index
        index_path = Path(str(bam_path) + ".bai")
        if index_path.exists():
            index_path.unlink()

        metadata_path = test_dir / "metadata.tsv"
        with open(metadata_path, 'w') as f:
            f.write("File accession\tBiosample term name\tAssay\tbam_path\n")
            f.write(f"rep1\tTest\tDNase-seq\t{bam_path}\n")

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(metadata_path),
            "--chrom-sizes", str(chrom_sizes_path),
            "--outdir", str(tmp_path / "output"),
        ]

        result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
        assert result.returncode != 0, "Pipeline should fail with missing BAM index"

    def test_chromosome_mismatch(self, tmp_path):
        """Pipeline should fail with chromosome mismatch."""
        test_dir = tmp_path / "test"
        test_dir.mkdir()

        # Create BAM with different chromosomes
        bam_chrom_sizes = {"chrX": 50000, "chrY": 40000}
        ref_chrom_sizes = {"chr1": 50000, "chr2": 40000}

        chrom_sizes_path = test_dir / "reference.chrom.sizes"
        with open(chrom_sizes_path, 'w') as f:
            for chrom, size in ref_chrom_sizes.items():
                f.write(f"{chrom}\t{size}\n")

        bam_dir = test_dir / "bams"
        bam_dir.mkdir()
        bam_path = bam_dir / "rep1.bam"
        create_test_bam(str(bam_path), bam_chrom_sizes, 100, seed=42)

        metadata_path = test_dir / "metadata.tsv"
        with open(metadata_path, 'w') as f:
            f.write("File accession\tBiosample term name\tAssay\tbam_path\n")
            f.write(f"rep1\tTest\tDNase-seq\t{bam_path}\n")

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(metadata_path),
            "--chrom-sizes", str(chrom_sizes_path),
            "--outdir", str(tmp_path / "output"),
        ]

        result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
        assert result.returncode != 0, "Pipeline should fail with chromosome mismatch"

    def test_error_logged_to_file(self, tmp_path):
        """Errors should be written to pipeline.log."""
        test_dir = tmp_path / "test"
        test_dir.mkdir()

        chrom_sizes_path = test_dir / "reference.chrom.sizes"
        with open(chrom_sizes_path, 'w') as f:
            f.write("chr1\t50000\n")

        metadata_path = test_dir / "metadata.tsv"
        with open(metadata_path, 'w') as f:
            f.write("File accession\tBiosample term name\tAssay\tbam_path\n")
            f.write("rep1\tErrorTest\tDNase-seq\t/nonexistent/file.bam\n")

        output_dir = tmp_path / "output"

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(metadata_path),
            "--chrom-sizes", str(chrom_sizes_path),
            "--outdir", str(output_dir),
        ]

        subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)

        # Check for log file
        log_path = output_dir / "ErrorTest" / "pipeline.log"
        if log_path.exists():
            log_content = log_path.read_text()
            assert "ERROR" in log_content or "error" in log_content.lower()


# ============================================================================
# Test: Multiple Cell Lines
# ============================================================================

class TestMultipleCellLines:
    """Tests for processing multiple cell lines."""

    @pytest.fixture(scope="class")
    def multi_cell_data(self, tmp_path_factory):
        """Create test data with multiple cell lines."""
        test_dir = tmp_path_factory.mktemp("multi_cell")

        chrom_sizes = {"chr1": 50000, "chr2": 40000}
        chrom_sizes_path = test_dir / "reference.chrom.sizes"
        with open(chrom_sizes_path, 'w') as f:
            for chrom, size in chrom_sizes.items():
                f.write(f"{chrom}\t{size}\n")

        bam_dir = test_dir / "bams"
        bam_dir.mkdir()

        cell_lines = {
            "CellA": ["a_rep1", "a_rep2"],
            "CellB": ["b_rep1"],
            "CellC": ["c_rep1", "c_rep2", "c_rep3"],
        }

        seed = 42
        for cell_line, replicates in cell_lines.items():
            for rep_id in replicates:
                create_test_bam(str(bam_dir / f"{rep_id}.bam"), chrom_sizes, 500, seed=seed)
                seed += 1

        metadata_path = test_dir / "metadata.tsv"
        with open(metadata_path, 'w') as f:
            f.write("File accession\tBiosample term name\tAssay\tbam_path\n")
            for cell_line, replicates in cell_lines.items():
                for rep_id in replicates:
                    f.write(f"{rep_id}\t{cell_line}\tDNase-seq\t{bam_dir}/{rep_id}.bam\n")

        return test_dir, cell_lines

    def test_all_cell_lines_processed(self, multi_cell_data, tmp_path):
        """All cell lines in metadata should be processed."""
        test_dir, cell_lines = multi_cell_data
        output_dir = tmp_path / "output"

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(test_dir / "metadata.tsv"),
            "--chrom-sizes", str(test_dir / "reference.chrom.sizes"),
            "--outdir", str(output_dir),
            "--threads", "1",
        ]

        result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
        assert result.returncode == 0

        for cell_line in cell_lines.keys():
            bw_path = output_dir / cell_line / "dnase_avg_norm100M.bw"
            assert bw_path.exists(), f"Missing output for {cell_line}"

    def test_correct_replicate_counts(self, multi_cell_data, tmp_path):
        """Manifest should report correct replicate counts."""
        test_dir, cell_lines = multi_cell_data
        output_dir = tmp_path / "output"

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(test_dir / "metadata.tsv"),
            "--chrom-sizes", str(test_dir / "reference.chrom.sizes"),
            "--outdir", str(output_dir),
            "--threads", "1",
        ]

        subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)

        for cell_line, replicates in cell_lines.items():
            manifest_path = output_dir / cell_line / "manifest.json"
            with open(manifest_path) as f:
                manifest = json.load(f)

            assert manifest["aggregation"]["n_replicates"] == len(replicates)
            assert len(manifest["replicates"]) == len(replicates)


# ============================================================================
# Test: Normalization Precision
# ============================================================================

class TestNormalizationPrecision:
    """Tests for normalization accuracy."""

    TARGET_TOTAL = 1e8
    TOLERANCE = 0.01

    @pytest.fixture(scope="class")
    def norm_test_data(self, tmp_path_factory):
        """Create test data for normalization tests."""
        test_dir = tmp_path_factory.mktemp("norm_test")

        chrom_sizes = {"chr1": 100000, "chr2": 80000, "chr3": 60000}
        chrom_sizes_path = test_dir / "reference.chrom.sizes"
        with open(chrom_sizes_path, 'w') as f:
            for chrom, size in chrom_sizes.items():
                f.write(f"{chrom}\t{size}\n")

        bam_dir = test_dir / "bams"
        bam_dir.mkdir()

        # Create BAMs with different read counts to test scaling
        create_test_bam(str(bam_dir / "low_rep.bam"), chrom_sizes, 100, seed=42)
        create_test_bam(str(bam_dir / "high_rep.bam"), chrom_sizes, 2000, seed=43)

        metadata_path = test_dir / "metadata.tsv"
        with open(metadata_path, 'w') as f:
            f.write("File accession\tBiosample term name\tAssay\tbam_path\n")
            f.write(f"low_rep\tNormTest\tDNase-seq\t{bam_dir}/low_rep.bam\n")
            f.write(f"high_rep\tNormTest\tDNase-seq\t{bam_dir}/high_rep.bam\n")

        return test_dir

    def test_normalization_within_tolerance(self, norm_test_data, tmp_path):
        """Final BigWig should be normalized within tolerance."""
        output_dir = tmp_path / "output"

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(norm_test_data / "metadata.tsv"),
            "--chrom-sizes", str(norm_test_data / "reference.chrom.sizes"),
            "--outdir", str(output_dir),
            "--threads", "1",
        ]

        result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
        assert result.returncode == 0

        bw_path = output_dir / "NormTest" / "dnase_avg_norm100M.bw"
        total = compute_bigwig_total(str(bw_path))

        relative_error = abs(total - self.TARGET_TOTAL) / self.TARGET_TOTAL
        assert relative_error <= self.TOLERANCE, (
            f"Total {total:.2e} differs from target {self.TARGET_TOTAL:.2e} by {relative_error:.4f}"
        )

    def test_manifest_normalization_consistent(self, norm_test_data, tmp_path):
        """Manifest normalization values should be consistent."""
        output_dir = tmp_path / "output"

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(norm_test_data / "metadata.tsv"),
            "--chrom-sizes", str(norm_test_data / "reference.chrom.sizes"),
            "--outdir", str(output_dir),
            "--threads", "1",
        ]

        subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)

        manifest_path = output_dir / "NormTest" / "manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        norm = manifest["normalization"]

        # Check scaling factor consistency
        expected_factor = norm["target_total"] / norm["pre_total"]
        assert abs(norm["scaling_factor"] - expected_factor) < 1e-6

        # Check post_total is close to target
        relative_error = abs(norm["post_total"] - norm["target_total"]) / norm["target_total"]
        assert relative_error <= norm["tolerance"]


# ============================================================================
# Test: Single Replicate
# ============================================================================

class TestSingleReplicate:
    """Tests for cell lines with single replicate."""

    def test_single_replicate_processed(self, tmp_path):
        """Cell line with single replicate should be processed correctly."""
        test_dir = tmp_path / "test"
        test_dir.mkdir()

        chrom_sizes = {"chr1": 50000}
        chrom_sizes_path = test_dir / "reference.chrom.sizes"
        with open(chrom_sizes_path, 'w') as f:
            f.write("chr1\t50000\n")

        bam_dir = test_dir / "bams"
        bam_dir.mkdir()
        create_test_bam(str(bam_dir / "single.bam"), chrom_sizes, 500, seed=42)

        metadata_path = test_dir / "metadata.tsv"
        with open(metadata_path, 'w') as f:
            f.write("File accession\tBiosample term name\tAssay\tbam_path\n")
            f.write(f"single\tSingleRep\tDNase-seq\t{bam_dir}/single.bam\n")

        output_dir = tmp_path / "output"

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(metadata_path),
            "--chrom-sizes", str(chrom_sizes_path),
            "--outdir", str(output_dir),
            "--threads", "1",
        ]

        result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
        assert result.returncode == 0

        bw_path = output_dir / "SingleRep" / "dnase_avg_norm100M.bw"
        assert bw_path.exists()

        manifest_path = output_dir / "SingleRep" / "manifest.json"
        with open(manifest_path) as f:
            manifest = json.load(f)

        assert manifest["aggregation"]["n_replicates"] == 1


# ============================================================================
# Test: Output File Naming
# ============================================================================

class TestOutputNaming:
    """Tests for strict output file naming conventions."""

    def test_exact_bigwig_name(self, tmp_path):
        """BigWig must be named exactly 'dnase_avg_norm100M.bw'."""
        test_dir = tmp_path / "test"
        test_dir.mkdir()

        chrom_sizes = {"chr1": 50000}
        chrom_sizes_path = test_dir / "reference.chrom.sizes"
        with open(chrom_sizes_path, 'w') as f:
            f.write("chr1\t50000\n")

        bam_dir = test_dir / "bams"
        bam_dir.mkdir()
        create_test_bam(str(bam_dir / "rep1.bam"), chrom_sizes, 500, seed=42)

        metadata_path = test_dir / "metadata.tsv"
        with open(metadata_path, 'w') as f:
            f.write("File accession\tBiosample term name\tAssay\tbam_path\n")
            f.write(f"rep1\tNameTest\tDNase-seq\t{bam_dir}/rep1.bam\n")

        output_dir = tmp_path / "output"

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(metadata_path),
            "--chrom-sizes", str(chrom_sizes_path),
            "--outdir", str(output_dir),
            "--threads", "1",
        ]

        subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)

        # Check exact filename
        expected_path = output_dir / "NameTest" / "dnase_avg_norm100M.bw"
        assert expected_path.exists(), f"Expected file not found: {expected_path}"

        # Ensure no other .bw files in directory
        cell_dir = output_dir / "NameTest"
        bw_files = list(cell_dir.glob("*.bw"))
        assert len(bw_files) == 1, f"Expected 1 .bw file, found {len(bw_files)}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
