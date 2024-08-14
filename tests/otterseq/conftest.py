"""Conftest file containing pytest fixture for unit test of otterseq."""

import pytest

from otterseq.snp import OtterSNP


@pytest.fixture()
def otter_snp() -> OtterSNP:
    """OtterSNP mock class."""
    return OtterSNP()


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
def no_bim_directory() -> str:
    """Directory containing no bim files."""
    return "tests"


@pytest.fixture()
def multiallelic_rsids(common_snps: list[str]):
    """Multiallelic (duplicated rsIDs) common snps."""
    multiallelic_rsids = [*common_snps, common_snps[0]]
    return multiallelic_rsids
