"""
Pytest configuration and shared fixtures.
"""

import os
import sys
from pathlib import Path

import pytest


# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


def pytest_configure(config):
    """Configure pytest."""
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )


@pytest.fixture(scope="session")
def project_root():
    """Return project root directory."""
    return PROJECT_ROOT
