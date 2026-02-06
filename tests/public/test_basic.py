#!/usr/bin/env python3
"""
Public Test Suite - Basic Functionality Tests

These tests verify the core functionality of the DNase preprocessing pipeline.
Students can use these tests to validate their implementation.
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
import pyBigWig

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture(scope="module")
def test_data_dir(tmp_path_factory):
    """Generate test data in a temporary directory."""
    test_dir = tmp_path_factory.mktemp("test_data")

    # Generate test data
    from scripts.generate_test_data import (
        create_chrom_sizes,
        create_synthetic_bam,
        create_metadata_tsv,
    )

    chrom_sizes = {
        "chr1": 50000,
        "chr2": 40000,
        "chr3": 30000,
    }

    chrom_sizes_path = test_dir / "reference.chrom.sizes"
    create_chrom_sizes(str(chrom_sizes_path), chrom_sizes)

    bam_dir = test_dir / "bams"
    bam_dir.mkdir()

    cell_lines = {
        "TestCell": ["rep1", "rep2"],
    }

    for cell_line, replicates in cell_lines.items():
        for i, rep_id in enumerate(replicates):
            bam_path = bam_dir / f"{rep_id}.bam"
            create_synthetic_bam(str(bam_path), chrom_sizes, n_reads=1000, seed=42 + i)

    metadata_path = test_dir / "metadata.tsv"
    create_metadata_tsv(str(metadata_path), cell_lines, str(bam_dir))

    return test_dir


@pytest.fixture(scope="module")
def pipeline_output(test_data_dir, tmp_path_factory):
    """Run the pipeline and return output directory."""
    output_dir = tmp_path_factory.mktemp("output")

    cmd = [
        "python", "run_dnase_tracks.py",
        "--metadata", str(test_data_dir / "metadata.tsv"),
        "--chrom-sizes", str(test_data_dir / "reference.chrom.sizes"),
        "--outdir", str(output_dir),
        "--threads", "1",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path(__file__).parent.parent.parent)


    output_dir_path = Path(output_dir)
    (output_dir_path / "_stdout.txt").write_text(result.stdout)
    (output_dir_path / "_stderr.txt").write_text(result.stderr)
    (output_dir_path / "_returncode.txt").write_text(str(result.returncode))

    return output_dir


# ============================================================================
# Test: BigWig Validity
# ============================================================================

class TestBigWigValidity:
    """Tests for BigWig file validity."""

    def test_bigwig_readable(self, pipeline_output):
        """BigWig should be readable with pyBigWig."""
        bw_path = Path(pipeline_output) / "TestCell" / "dnase_avg_norm100M.bw"

        with pyBigWig.open(str(bw_path)) as bw:
            assert bw is not None
            assert bw.isBigWig()


# ============================================================================
# Test: Manifest Accuracy
# ============================================================================

class TestManifestAccuracy:
    """Tests for manifest.json accuracy."""

    def test_manifest_non_zero_average_mean_matches_bigwig(self, pipeline_output):
        """Manifest non_zero_average_mean should match computed value from BigWig."""
        manifest_path = Path(pipeline_output) / "TestCell" / "manifest.json"
        bw_path = Path(pipeline_output) / "TestCell" / "dnase_avg_norm100M.bw"

        with open(manifest_path) as f:
            manifest = json.load(f)

        non_zero_values = []
        with pyBigWig.open(str(bw_path)) as bw:
            for chrom in bw.chroms():
                intervals = bw.intervals(chrom)
                if intervals:
                    for start, end, value in intervals:
                        if value != 0:
                            non_zero_values.append(value)

        computed_mean = np.mean(non_zero_values) if non_zero_values else 0.0
        manifest_mean = manifest["aggregation"]["non_zero_average_mean"]

        assert abs(computed_mean - manifest_mean) < 1e-6, (
            f"non_zero_average_mean mismatch: manifest={manifest_mean}, computed={computed_mean}"
        )


# ============================================================================
# Test: Expected Values Validation (hash-based, values hidden)
# ============================================================================

class TestExpectedValues:
    """Tests that validate manifest values against reference data.

    Expected values are stored as hashes so candidates can verify correctness
    without seeing the actual expected values.
    """


    EXPECTED_HASHES = {
        "GM12878": [
            "2bb5a2194f062fe01073e1a714faf9e34032cf2afb6c224010aaa1b1090b510e", 
            "c558fa9c4a4a6ad4962c10fa6b220b5a9a7129ac62c2f17ab96ff96ff7c29024",
            "4348672d1a85fab51648c276be79cc8ee80608359b01705cd9a7f30dba2cb9dc",
        ],
        "HeLa-S3": [
            "90e683111c8a35ea92f86870185d5350d940315840e636019baa36cd37bb1dc7",
            "c85e09eb8176db9e00c43c8bdc564dd91330c9d14f2a52ad6bc7872931e4c7ef",
            "6669743c77ff5bb876bc7aa9cd5d54b997bee0bd5b6e1f86a1516feb574b0ff6",
        ],
        "SK-N-SH": [
            "18db98e3760ea46503aa9f79ad850cc408d0192446a5666220ca3b99346e3e8c",
            "67194cd2e7ca284251eab9f18c295f7f031fb464267547240d7350e2707895a9",
            "c68fa0c5dd8c9a2cf6555708f938a5d592d20b0d81f0d40ac6221d07500ae145",
        ],
    }

    @staticmethod
    def _hash_value(value: float) -> str:
        """Hash a float value rounded to 2 decimal places."""
        import hashlib
        rounded = f"{value:.2f}"
        return hashlib.sha256(rounded.encode()).hexdigest()

    @pytest.fixture(scope="class")
    def real_output_dir(self):
        """Return the output directory for real ENCODE data runs."""
        return Path(__file__).parent.parent.parent / "out"

    def test_gm12878_non_zero_mean(self, real_output_dir):
        """GM12878 non_zero_average_mean should match expected value (±0.01 tolerance)."""
        manifest_path = real_output_dir / "GM12878" / "manifest.json"
        if not manifest_path.exists():
            pytest.skip("GM12878 output not found - run pipeline first")

        with open(manifest_path) as f:
            manifest = json.load(f)

        actual = manifest["aggregation"]["non_zero_average_mean"]
        actual_hash = self._hash_value(actual)

        assert actual_hash in self.EXPECTED_HASHES["GM12878"], (
            f"GM12878 non_zero_average_mean incorrect (got {actual:.2f})"
        )

    def test_hela_s3_non_zero_mean(self, real_output_dir):
        """HeLa-S3 non_zero_average_mean should match expected value (±0.01 tolerance)."""
        manifest_path = real_output_dir / "HeLa-S3" / "manifest.json"
        if not manifest_path.exists():
            pytest.skip("HeLa-S3 output not found - run pipeline first")

        with open(manifest_path) as f:
            manifest = json.load(f)

        actual = manifest["aggregation"]["non_zero_average_mean"]
        actual_hash = self._hash_value(actual)

        assert actual_hash in self.EXPECTED_HASHES["HeLa-S3"], (
            f"HeLa-S3 non_zero_average_mean incorrect (got {actual:.2f})"
        )

    def test_sk_n_sh_non_zero_mean(self, real_output_dir):
        """SK-N-SH non_zero_average_mean should match expected value (±0.01 tolerance)."""
        manifest_path = real_output_dir / "SK-N-SH" / "manifest.json"
        if not manifest_path.exists():
            pytest.skip("SK-N-SH output not found - run pipeline first")

        with open(manifest_path) as f:
            manifest = json.load(f)

        actual = manifest["aggregation"]["non_zero_average_mean"]
        actual_hash = self._hash_value(actual)

        assert actual_hash in self.EXPECTED_HASHES["SK-N-SH"], (
            f"SK-N-SH non_zero_average_mean incorrect (got {actual:.2f})"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
