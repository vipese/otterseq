"""otterseq unit tests."""

import os

import pytest

from otterseq.errors import MultiAllelicError
from otterseq.snp import OtterSNP


@pytest.mark.parametrize(
    argnames=["filepath", "outpath"],
    argvalues=[(1000, "outpath"), ("filepath", 1000)],
)
def test_fail_binarize(otter_snp: OtterSNP, filepath: str, outpath: str):
    """Test correct raise of errors in OtterSNP.

    Args:
        otter_snp (OtterSNP): OtterSNP class
        filepath (str): Path to file or folder with files.
        outpath (str): Output path to write output.
    """
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
    """Test binarization of files.

    Use cases:
        - Folder with multiple PLINK files inside.
        - File without suffix.
        - File with .ped suffix.

    Args:
        otter_snp (OtterSNP): OtterSNP class
        filepath (str): Path to file or folder with files.
        outpath (str): Output path to write output.
    """
    os.makedirs(outpath, exist_ok=True)
    otter_snp.binarize_files(filepath, outpath)

    filenames = ["toy", "toy_2"] if filepath == "tests/data" else ["toy"]
    filenames = [os.path.join(outpath, fname) for fname in filenames]
    suffixes = [".pgen", ".psam", ".pvar", ".log"]
    for suffix in suffixes:
        for filename in filenames:
            outfile_path = filename + suffix
            assert os.path.isfile(
                outfile_path
            ), f"{outfile_path} was not correctly created"
            os.remove(outfile_path)
    os.removedirs(outpath)


@pytest.mark.parametrize(
    argnames=["filepath", "write", "outpath"],
    argvalues=[
        (1000, False, None),
        ("filepath", "write", None),
        ("filepath", False, 1000),
    ],
)
def test_type_get_common_snp(
    otter_snp: OtterSNP, filepath: str, write: bool, outpath: str
):
    """Test fails of get_common_snp."""
    with pytest.raises(TypeError):
        otter_snp.get_common_snp(filepath, write, outpath)


@pytest.mark.parametrize(
    argnames=["filepath", "write", "outpath"],
    argvalues=[
        ("filepath", True, None),
    ],
)
def test_value_get_common_snp(
    otter_snp: OtterSNP, filepath: str, write: bool, outpath: str
):
    """Test fails of get_common_snp."""
    with pytest.raises(
        ValueError, match="Write was set to True but no outpath was provided"
    ):
        otter_snp.get_common_snp(filepath, write, outpath)


def test_get_common_snps_no_bim(otter_snp: OtterSNP, no_bim_directory: str):
    """Test get_common_snp passing directory with no .bim files."""
    with pytest.raises(FileNotFoundError):
        otter_snp.get_common_snp(filepath=no_bim_directory)


def test_check_multi_allelic(  # noqa: D103
    otter_snp: OtterSNP, multiallelic_rsids: list[str]
) -> None:
    with pytest.raises(MultiAllelicError):
        otter_snp.check_multi_allelic(multiallelic_rsids)


@pytest.mark.parametrize(
    argnames=["filepath", "write", "outpath"],
    argvalues=[
        ("tests/data", False, None),
        ("tests/data", True, "tests/output_test"),
    ],
)
def test_get_common_snp(
    otter_snp: OtterSNP,
    filepath: str,
    write: bool,
    outpath: str | None,
    common_snps: list[str],
):
    """Test get_common snp.

    Use cases:
        - File path, no write.
        - File path and write.
    """
    if outpath is not None:
        os.makedirs(outpath, exist_ok=True)
    snps = otter_snp.get_common_snp(filepath, write, outpath)
    assert set(snps) == set(common_snps), "common_snps don't match expected"
    if outpath is not None:
        out_file = os.path.join(outpath, "common_snps.txt")
        assert os.path.isfile(out_file), "Common SNPs file not created"
        with open(out_file) as in_file:
            file_snps = in_file.read().splitlines()
        assert set(file_snps) == set(
            common_snps
        ), "common snps in file don't match expected"
    if outpath is not None:
        os.remove(out_file)
        os.removedirs(outpath)
