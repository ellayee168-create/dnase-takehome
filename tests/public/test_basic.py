#!/usr/bin/env python3
"""
Public Test Suite - Basic Functionality Tests

These tests verify the core functionality of the DNase preprocessing pipeline.
Students can use these tests to validate their implementation.
"""

import json
import os
import subprocess
import sys
import tempfile
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

    # Small chromosome sizes for fast testing
    chrom_sizes = {
        "chr1": 50000,
        "chr2": 40000,
        "chr3": 30000,
    }

    # Create chrom sizes
    chrom_sizes_path = test_dir / "reference.chrom.sizes"
    create_chrom_sizes(str(chrom_sizes_path), chrom_sizes)

    # Create BAMs
    bam_dir = test_dir / "bams"
    bam_dir.mkdir()

    cell_lines = {
        "TestCell": ["rep1", "rep2"],
    }

    for cell_line, replicates in cell_lines.items():
        for i, rep_id in enumerate(replicates):
            bam_path = bam_dir / f"{rep_id}.bam"
            create_synthetic_bam(str(bam_path), chrom_sizes, n_reads=1000, seed=42 + i)

    # Create metadata
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

    # Store result for debugging
    output_dir_path = Path(output_dir)
    (output_dir_path / "_stdout.txt").write_text(result.stdout)
    (output_dir_path / "_stderr.txt").write_text(result.stderr)
    (output_dir_path / "_returncode.txt").write_text(str(result.returncode))

    return output_dir


# ============================================================================
# Test: Output Structure
# ============================================================================

class TestOutputStructure:
    """Tests for correct output file structure."""

    def test_output_directory_exists(self, pipeline_output):
        """Output directory should be created."""
        assert Path(pipeline_output).exists()

    def test_cell_line_directory_exists(self, pipeline_output):
        """Cell line subdirectory should be created."""
        cell_dir = Path(pipeline_output) / "TestCell"
        assert cell_dir.exists(), f"Cell line directory not found: {cell_dir}"

    def test_bigwig_exists(self, pipeline_output):
        """Final BigWig file should exist."""
        bw_path = Path(pipeline_output) / "TestCell" / "dnase_avg_norm100M.bw"
        assert bw_path.exists(), f"BigWig not found: {bw_path}"

    def test_manifest_exists(self, pipeline_output):
        """Manifest JSON should exist."""
        manifest_path = Path(pipeline_output) / "TestCell" / "manifest.json"
        assert manifest_path.exists(), f"Manifest not found: {manifest_path}"

    def test_pipeline_log_exists(self, pipeline_output):
        """Pipeline log should exist."""
        log_path = Path(pipeline_output) / "TestCell" / "pipeline.log"
        assert log_path.exists(), f"Log not found: {log_path}"


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

    def test_bigwig_has_chromosomes(self, pipeline_output):
        """BigWig should contain chromosome data."""
        bw_path = Path(pipeline_output) / "TestCell" / "dnase_avg_norm100M.bw"

        with pyBigWig.open(str(bw_path)) as bw:
            chroms = bw.chroms()
            assert len(chroms) > 0, "BigWig has no chromosomes"

    def test_bigwig_has_signal(self, pipeline_output):
        """BigWig should contain signal values."""
        bw_path = Path(pipeline_output) / "TestCell" / "dnase_avg_norm100M.bw"

        with pyBigWig.open(str(bw_path)) as bw:
            total_signal = 0
            for chrom in bw.chroms():
                intervals = bw.intervals(chrom)
                if intervals:
                    for start, end, value in intervals:
                        total_signal += value * (end - start)

            assert total_signal > 0, "BigWig has zero total signal"


# ============================================================================
# Test: Normalization
# ============================================================================

class TestNormalization:
    """Tests for signal normalization."""

    TARGET_TOTAL = 1e8
    TOLERANCE = 0.01

    def test_normalized_to_target(self, pipeline_output):
        """BigWig total signal should be normalized to 10^8."""
        bw_path = Path(pipeline_output) / "TestCell" / "dnase_avg_norm100M.bw"

        with pyBigWig.open(str(bw_path)) as bw:
            total_signal = 0
            for chrom in bw.chroms():
                intervals = bw.intervals(chrom)
                if intervals:
                    for start, end, value in intervals:
                        total_signal += value * (end - start)

        relative_error = abs(total_signal - self.TARGET_TOTAL) / self.TARGET_TOTAL
        assert relative_error <= self.TOLERANCE, (
            f"Normalization error too high: {relative_error:.4f} > {self.TOLERANCE}"
        )


# ============================================================================
# Test: Manifest Contents
# ============================================================================

class TestManifest:
    """Tests for manifest.json contents."""

    REQUIRED_KEYS = [
        "cell_line",
        "replicates",
        "chrombpnet",
        "commands",
        "normalization",
        "aggregation",
    ]

    def test_manifest_valid_json(self, pipeline_output):
        """Manifest should be valid JSON."""
        manifest_path = Path(pipeline_output) / "TestCell" / "manifest.json"

        with open(manifest_path) as f:
            manifest = json.load(f)

        assert isinstance(manifest, dict)

    def test_manifest_has_required_keys(self, pipeline_output):
        """Manifest should contain all required keys."""
        manifest_path = Path(pipeline_output) / "TestCell" / "manifest.json"

        with open(manifest_path) as f:
            manifest = json.load(f)

        for key in self.REQUIRED_KEYS:
            assert key in manifest, f"Missing key in manifest: {key}"

    def test_manifest_cell_line(self, pipeline_output):
        """Manifest cell_line should match directory."""
        manifest_path = Path(pipeline_output) / "TestCell" / "manifest.json"

        with open(manifest_path) as f:
            manifest = json.load(f)

        assert manifest["cell_line"] == "TestCell"

    def test_manifest_replicates_structure(self, pipeline_output):
        """Manifest replicates should have correct structure."""
        manifest_path = Path(pipeline_output) / "TestCell" / "manifest.json"

        with open(manifest_path) as f:
            manifest = json.load(f)

        replicates = manifest["replicates"]
        assert isinstance(replicates, list)
        assert len(replicates) > 0

        for rep in replicates:
            assert "replicate_id" in rep
            assert "bam_path" in rep

    def test_manifest_chrombpnet_structure(self, pipeline_output):
        """Manifest chrombpnet should have repo_url and commit_sha."""
        manifest_path = Path(pipeline_output) / "TestCell" / "manifest.json"

        with open(manifest_path) as f:
            manifest = json.load(f)

        chrombpnet = manifest["chrombpnet"]
        assert "repo_url" in chrombpnet
        assert "commit_sha" in chrombpnet

    def test_manifest_normalization_structure(self, pipeline_output):
        """Manifest normalization should have required fields."""
        manifest_path = Path(pipeline_output) / "TestCell" / "manifest.json"

        with open(manifest_path) as f:
            manifest = json.load(f)

        norm = manifest["normalization"]
        required = ["target_total", "pre_total", "post_total", "scaling_factor", "tolerance"]

        for key in required:
            assert key in norm, f"Missing normalization key: {key}"

        assert norm["target_total"] == 1e8

    def test_manifest_aggregation_structure(self, pipeline_output):
        """Manifest aggregation should have method and n_replicates."""
        manifest_path = Path(pipeline_output) / "TestCell" / "manifest.json"

        with open(manifest_path) as f:
            manifest = json.load(f)

        agg = manifest["aggregation"]
        assert agg["method"] == "mean"
        assert agg["n_replicates"] == 2


# ============================================================================
# Test: Exit Codes
# ============================================================================

class TestExitCodes:
    """Tests for correct exit codes."""

    def test_success_exit_code(self, test_data_dir, tmp_path):
        """Pipeline should exit with code 0 on success."""
        output_dir = tmp_path / "output"

        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(test_data_dir / "metadata.tsv"),
            "--chrom-sizes", str(test_data_dir / "reference.chrom.sizes"),
            "--outdir", str(output_dir),
            "--threads", "1",
        ]

        result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
        assert result.returncode == 0, f"Expected exit code 0, got {result.returncode}"


# ============================================================================
# Test: CLI Arguments
# ============================================================================

class TestCLI:
    """Tests for CLI argument handling."""

    def test_missing_metadata_fails(self, test_data_dir, tmp_path):
        """Pipeline should fail with missing metadata."""
        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", "/nonexistent/metadata.tsv",
            "--chrom-sizes", str(test_data_dir / "reference.chrom.sizes"),
            "--outdir", str(tmp_path / "output"),
        ]

        result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
        assert result.returncode != 0

    def test_missing_chrom_sizes_fails(self, test_data_dir, tmp_path):
        """Pipeline should fail with missing chrom sizes."""
        cmd = [
            "python", "run_dnase_tracks.py",
            "--metadata", str(test_data_dir / "metadata.tsv"),
            "--chrom-sizes", "/nonexistent/chrom.sizes",
            "--outdir", str(tmp_path / "output"),
        ]

        result = subprocess.run(cmd, capture_output=True, cwd=Path(__file__).parent.parent.parent)
        assert result.returncode != 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
