"""OtterSNP unit tests."""

import os
import tempfile
from typing import Any

import pytest

from otterseq.snp import OtterSNP


def test_read_snp_id_bim_success(otter_snp: OtterSNP) -> None:
    """Test successful reading of SNP IDs from .bim file."""
    # Create a temporary .bim file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bim", delete=False
    ) as f:
        f.write("1\trs0001\t0\t100100\tG\tA\n")
        f.write("1\trs0002\t0\t100200\tA\tT\n")
        f.write("1\trs0003\t0\t100300\tC\tG\n")
        temp_file = f.name

    try:
        snp_ids = otter_snp._read_snp_id_bim(temp_file)

        expected_ids = ["rs0001", "rs0002", "rs0003"]
        assert snp_ids == expected_ids

    finally:
        os.unlink(temp_file)


def test_read_snp_id_bim_empty_file(otter_snp: OtterSNP) -> None:
    """Test reading empty .bim file."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".bim", delete=False
    ) as f:
        temp_file = f.name

    try:
        snp_ids = otter_snp._read_snp_id_bim(temp_file)
        assert snp_ids == []

    finally:
        os.unlink(temp_file)


def test_read_snp_id_bim_file_not_found(otter_snp: OtterSNP) -> None:
    """Test that FileNotFoundError is raised for non-existent .bim file."""
    with pytest.raises(FileNotFoundError):
        otter_snp._read_snp_id_bim("non_existent_file.bim")


def test_process_bim_files_success(otter_snp: OtterSNP) -> None:
    """Test successful processing of .bim files."""
    # Create a temporary directory with .bim files
    with tempfile.TemporaryDirectory() as temp_dir:
        # First file with SNPs rs0001, rs0002, rs0003
        bim1 = os.path.join(temp_dir, "file1.bim")
        with open(bim1, "w") as f:
            f.write("1\trs0001\t0\t100100\tG\tA\n")
            f.write("1\trs0002\t0\t100200\tA\tT\n")
            f.write("1\trs0003\t0\t100300\tC\tG\n")

        # Second file with SNPs rs0002, rs0003, rs0004
        bim2 = os.path.join(temp_dir, "file2.bim")
        with open(bim2, "w") as f:
            f.write("1\trs0002\t0\t100200\tA\tT\n")
            f.write("1\trs0003\t0\t100300\tC\tG\n")
            f.write("1\trs0004\t0\t100400\tT\tC\n")

        snp_sets = otter_snp._process_bim_files(temp_dir)

        expected_sets = [
            {"rs0001", "rs0002", "rs0003"},
            {"rs0002", "rs0003", "rs0004"},
        ]
        assert len(snp_sets) == 2
        # Since os.listdir() doesn't guarantee order, we need to check that both expected sets are present
        # Convert to frozensets for comparison
        snp_sets_frozen = {frozenset(s) for s in snp_sets}
        expected_sets_frozen = {frozenset(s) for s in expected_sets}
        assert snp_sets_frozen == expected_sets_frozen


def test_process_bim_files_empty_file(otter_snp: OtterSNP) -> None:
    """Test processing .bim files with an empty file."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create an empty .bim file
        empty_bim = os.path.join(temp_dir, "empty.bim")
        with open(empty_bim, "w") as f:
            pass  # Create empty file

        # Create a normal .bim file
        normal_bim = os.path.join(temp_dir, "normal.bim")
        with open(normal_bim, "w") as f:
            f.write("1\trs0001\t0\t100100\tG\tA\n")
            f.write("1\trs0002\t0\t100200\tA\tT\n")

        snp_sets = otter_snp._process_bim_files(temp_dir)

        # Should only have one set (from the normal file)
        assert len(snp_sets) == 1
        assert snp_sets[0] == {"rs0001", "rs0002"}


def test_process_bim_files_no_bim_files(otter_snp: OtterSNP) -> None:
    """Test processing directory with no .bim files."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a file that's not .bim
        with open(os.path.join(temp_dir, "not_bim.txt"), "w") as f:
            f.write("some content")

        with pytest.raises(FileNotFoundError, match="No .bim files found"):
            otter_snp._process_bim_files(temp_dir)


def test_process_bim_files_file_error(otter_snp: OtterSNP) -> None:
    """Test processing .bim files with corrupted file."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a corrupted .bim file
        corrupted_bim = os.path.join(temp_dir, "corrupted.bim")
        with open(corrupted_bim, "w") as f:
            f.write("invalid\tformat\n")

        with pytest.raises(Exception):  # noqa: B017, PT011
            otter_snp._process_bim_files(temp_dir)


def test_find_common_snps_success(otter_snp: OtterSNP) -> None:
    """Test finding common SNPs across multiple sets."""
    snp_sets = [
        {"rs0001", "rs0002", "rs0003"},
        {"rs0002", "rs0003", "rs0004"},
        {"rs0003", "rs0004", "rs0005"},
    ]

    common_snps = otter_snp._find_common_snps(snp_sets)

    # Only rs0003 should be common across all three sets
    assert common_snps == ["rs0003"]


def test_find_common_snps_no_common(otter_snp: OtterSNP) -> None:
    """Test finding common SNPs when there are none."""
    snp_sets = [
        {"rs0001", "rs0002"},
        {"rs0003", "rs0004"},
        {"rs0005", "rs0006"},
    ]

    common_snps = otter_snp._find_common_snps(snp_sets)

    assert common_snps == []


def test_find_common_snps_empty_sets(otter_snp: OtterSNP) -> None:
    """Test finding common SNPs with empty sets."""
    snp_sets: list[set[str]] = []

    common_snps = otter_snp._find_common_snps(snp_sets)

    assert common_snps == []


def test_find_common_snps_single_set(otter_snp: OtterSNP) -> None:
    """Test finding common SNPs with only one set."""
    snp_sets = [{"rs0001", "rs0002", "rs0003"}]

    common_snps = otter_snp._find_common_snps(snp_sets)

    # Should return all SNPs from the single set, sorted
    assert common_snps == ["rs0001", "rs0002", "rs0003"]


def test_find_common_snps_all_common(otter_snp: OtterSNP) -> None:
    """Test finding common SNPs when all sets have the same SNPs."""
    snp_sets = [
        {"rs0001", "rs0002", "rs0003"},
        {"rs0001", "rs0002", "rs0003"},
        {"rs0001", "rs0002", "rs0003"},
    ]

    common_snps = otter_snp._find_common_snps(snp_sets)

    assert common_snps == ["rs0001", "rs0002", "rs0003"]


@pytest.mark.parametrize(
    argnames=("filepath", "outpath"),
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
    suffixes = [".bed", ".bim", ".fam", ".log"]
    for suffix in suffixes:
        for filename in filenames:
            outfile_path = filename + suffix
            assert os.path.isfile(
                outfile_path
            ), f"{outfile_path} was not correctly created"
            os.remove(outfile_path)
    os.removedirs(outpath)


def test_fnf__error_get_common_snps(
    otter_snp: OtterSNP, no_files_directory: str
) -> None:
    """Test get_common_snp passing directory with no .bim files."""
    with pytest.raises(FileNotFoundError):
        otter_snp.get_common_snp(filepath=no_files_directory)


@pytest.mark.parametrize(
    argnames=("filepath", "write", "outpath"),
    argvalues=[
        ("filepath", True, None),
    ],
)
def test_value_error_get_common_snp(
    otter_snp: OtterSNP, filepath: str, write: bool, outpath: str
) -> None:
    """Test fails of get_common_snp."""
    with pytest.raises(
        ValueError, match="Write was set to True but no outpath was provided"
    ):
        otter_snp.get_common_snp(filepath, write, outpath)


@pytest.mark.parametrize(
    argnames=("filepath", "write", "outpath"),
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
) -> None:
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


def test_get_common_snp_with_empty_file(otter_snp: OtterSNP) -> None:
    """Test get_common_snp with a directory containing an empty .bim file."""
    # Create a temporary directory with an empty .bim file
    with tempfile.TemporaryDirectory() as temp_dir:
        empty_bim_file = os.path.join(temp_dir, "empty.bim")
        with open(empty_bim_file, "w"):
            pass  # Create empty file

        # Should return empty list when one file has no SNPs
        snps = otter_snp.get_common_snp(temp_dir)
        assert snps == []


def test_get_common_snp_early_termination(otter_snp: OtterSNP) -> None:
    """Test that get_common_snp terminates early when no common SNPs remain."""
    # Create temporary directory with .bim files that have no common SNPs
    with tempfile.TemporaryDirectory() as temp_dir:
        # First file with SNPs rs0001, rs0002
        bim1 = os.path.join(temp_dir, "file1.bim")
        with open(bim1, "w") as f:
            f.write("1\trs0001\t0\t100100\tG\tA\n")
            f.write("1\trs0002\t0\t100200\tA\tT\n")

        # Second file with SNPs rs0003, rs0004 (no overlap)
        bim2 = os.path.join(temp_dir, "file2.bim")
        with open(bim2, "w") as f:
            f.write("1\trs0003\t0\t100300\tC\tG\n")
            f.write("1\trs0004\t0\t100400\tT\tC\n")

        # Should return empty list when no common SNPs
        snps = otter_snp.get_common_snp(temp_dir)
        assert snps == []


def test_type_merge_files_no_pgen(
    otter_snp: OtterSNP, no_files_directory: str
) -> None:
    """Test get_common_snp passing directory with no .bim files."""
    with pytest.raises(FileNotFoundError):
        otter_snp.merge_files(filepath=no_files_directory)


@pytest.mark.parametrize(
    argnames=("args", "expected_snps"),
    argvalues=[
        (
            {
                "filepath": "tests/data",
            },
            "common_snps",
        ),
        (
            {"filepath": "tests/data", "outpath": "tests/output_test"},
            "common_snps",
        ),
        (
            {"filepath": "tests/data", "prefix": "merged_pgen"},
            "common_snps",
        ),
        (
            {
                "filepath": "tests/data",
                "prefix": "merged_pgen",
                "only_common": False,
            },
            "all_snps",
        ),
        (
            {
                "filepath": "tests/data",
                "prefix": "merged_pgen",
                "only_common": True,
            },
            "common_snps",
        ),
    ],
)
def test_merge(
    otter_snp: OtterSNP,
    args: dict[str, Any],
    expected_snps: str,
    request: pytest.FixtureRequest,
) -> None:
    """Test merge_files.

    Use cases:
        - File path only.
        - File path and outpath.
        - File path and prefix.
        - File path, prefix and only_common=False.
        - File path, prefix and only_common=True.
    """
    expected_snps_list = request.getfixturevalue(expected_snps)
    outpath = args.get("outpath", args["filepath"])
    prefix = args.get("prefix", "merged_snps")
    os.makedirs(outpath, exist_ok=True)

    otter_snp.merge_files(**args)

    # Check that files were created
    suffixes = [".bed", ".bim", ".fam", ".log"]
    for suffix in suffixes:
        outfile_path = os.path.join(outpath, prefix + suffix)
        assert os.path.isfile(
            outfile_path
        ), f"{outfile_path} was not correctly created"
        os.remove(outfile_path)

    # Check that common_snps.txt was created if only_common=True
    if args.get("only_common", True):
        common_snps_path = os.path.join(outpath, "common_snps.txt")
        assert os.path.isfile(
            common_snps_path
        ), "common_snps.txt was not created"
        with open(common_snps_path) as in_file:
            file_snps = in_file.read().splitlines()
        assert set(file_snps) == set(
            expected_snps_list
        ), "common snps in file don't match expected"
        os.remove(common_snps_path)

    # Check that merge_list.txt was created
    merge_list_path = os.path.join(outpath, "merge_list.txt")
    assert os.path.isfile(merge_list_path), "merge_list.txt was not created"
    os.remove(merge_list_path)

    # Only remove the output directory if it's different from the input directory
    if outpath != args["filepath"]:
        os.removedirs(outpath)
