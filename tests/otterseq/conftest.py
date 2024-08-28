"""Conftest file containing pytest fixture for unit test of otterseq."""

import os

import pytest

from otterseq.pca import OtterPCA
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
def otter_pca() -> OtterPCA:
    """OtterPCA mock class."""
    return OtterPCA()


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
def all_snps() -> list[str]:
    """All unique rsIDs in toy files."""
    all_snps = [
        "rs0001",
        "rs0002",
        "rs0003",
        "rs0004",
        "rs0005",
        "rs0006",
        "rs0007",
    ]
    all_snps.sort()
    return all_snps


@pytest.fixture()
def filepath() -> str:
    """Path to test files."""
    return "tests/data/"


@pytest.fixture()
def filepath_qc() -> str:
    """Path to test files."""
    return "tests/data_qc/"


@pytest.fixture()
def filepath_pca() -> str:
    """Path to test PCA files."""
    return "tests/data_pca/"


@pytest.fixture()
def filename(filepath: str) -> str:
    """Path fo PLINK1.9 file."""
    return os.path.join(filepath, "toy")


@pytest.fixture()
def filename_dup(filepath_qc: str) -> str:
    """Path fo PLINK1.9 file."""
    return os.path.join(filepath_qc, "toy")


@pytest.fixture()
def filename_qc(filepath_qc: str) -> str:
    """Path fo PLINK1.9 file."""
    return os.path.join(filepath_qc, "toy_qc")


@pytest.fixture()
def filename_pca(filepath_pca: str) -> str:
    """Path fo PLINK1.9 file."""
    return os.path.join(filepath_pca, "toy")


@pytest.fixture()
def no_filename(filepath: str) -> str:
    """Path to non existent PLINK1.9 file."""
    return os.path.join(filepath, "no_file")


@pytest.fixture()
def no_files_directory() -> str:
    """Directory containing no bim files."""
    return "tests"


@pytest.fixture()
def multiallelic_rsids(common_snps: list[str]) -> list[str]:
    """Multiallelic (duplicated rsIDs) common snps."""
    multiallelic_rsids = [*common_snps, common_snps[0]]
    return multiallelic_rsids


@pytest.fixture()
def duplicated_variants() -> list[str]:
    """Duplicated variants based on coordinates and allele codes."""
    return ["rs0001_dup1", "rs0001_dup2", "rs0002_dup1"]


@pytest.fixture()
def duplicated_rsid() -> list[str]:
    """Duplicated variants based on coordinates and allele codes."""
    return ["rs0001", "rs0002"]
