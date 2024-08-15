"""Conftest file containing pytest fixture for unit test of otterseq."""

import os

import pytest

from otterseq.qc import OtterQC
from otterseq.snp import OtterSNP


@pytest.fixture()
def otter_snp() -> OtterSNP:
    """OtterSNP mock class."""
    return OtterSNP()


@pytest.fixture()
def otter_qc() -> OtterQC:
    """OtterQC mock class."""
    return OtterQC()


@pytest.fixture()
def common_snps() -> list[str]:
    """Common SNPs expected output."""  # noqa: D401
    common_snps = [
        "rs0003",
        "rs0005",
        "rs0001",
    ]
    common_snps.sort()
    return common_snps


@pytest.fixture()
def filepath() -> str:
    """Path to test files."""
    return "tests/data/"


@pytest.fixture()
def filename(filepath: str) -> str:
    """Path fo PLINK1.9 file."""
    return os.path.join(filepath, "toy")


@pytest.fixture()
def no_filename(filepath: str) -> str:
    """Path to non existent PLINK1.9 file."""
    return os.path.join(filepath, "no_file")


@pytest.fixture()
def no_files_directory() -> str:
    """Directory containing no bim files."""
    return "tests"


@pytest.fixture()
def multiallelic_rsids(common_snps: list[str]):
    """Multiallelic (duplicated rsIDs) common snps."""
    multiallelic_rsids = [*common_snps, common_snps[0]]
    return multiallelic_rsids
