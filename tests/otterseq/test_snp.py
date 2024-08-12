import pytest
from otterseq.snp import OtterSNP
import os


@pytest.mark.parametrize(
    argnames=["filepath", "outpath"], argvalues=[(1000, "outpath"), ("filepath", 1000)]
)
def test_fail_binarize(otter_snp: OtterSNP, filepath: str, outpath: str):
    with pytest.raises(TypeError, match="must be of type str"):
        otter_snp.binarize_files(filepath, outpath)


@pytest.mark.parametrize(
    argnames=["filepath", "outpath"],
    argvalues=[
        ("tests/data", "tests/output_test"),
        ("tests/data/toy", "tests/output_test"),
        ("tests/data/toy.ped", "tests/output_test"),
    ],
)
def test_binarize(otter_snp: OtterSNP, filepath: str, outpath: str) -> None:
    otter_snp.binarize_files(filepath, outpath)

    filename = os.path.join(outpath, "toy")
    suffixes = [".bed", ".fam", ".bim", ".log"]
    for suffix in suffixes:
        outfile_path = filename + suffix
        assert os.path.isfile(outfile_path), f"{outfile_path} was not correctly created"
        os.remove(outfile_path)
