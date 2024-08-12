from otterseq.snp import OtterSNP
import pytest


@pytest.fixture()
def otter_snp():
    return OtterSNP()
