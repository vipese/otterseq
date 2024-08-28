"""OtterQC unit tests."""

import os
from typing import Any

import pandas as pd
import pytest
from beartype.roar import BeartypeCallHintParamViolation

from otterseq.qc import OtterQC


@pytest.mark.parametrize(
    argnames=["args", "error"],
    argvalues=[
        ({"filename": 1000}, TypeError),
        ({"filename": "tests", "threshold": "threshold"}, TypeError),
        ({"filename": "tests/toy"}, FileNotFoundError),
        ({"filename": "tests/data/toy", "threshold": -10}, ValueError),
        ({"filename": "tests/data/toy", "threshold": 1}, ValueError),
        ({"filename": "tests/data/toy", "threshold": 10}, ValueError),
    ],
)
def test_ibd_errors(  # noqa: D103
    otter_qc: OtterQC, args: dict[str, Any], error: Exception
) -> None:
    with pytest.raises(error):
        otter_qc.ibd(**args)


@pytest.mark.parametrize(
    argnames=["filename", "threshold", "expected_iid"],
    argvalues=[
        ("tests/data/toy", 0.25, ["7"]),
        ("tests/data/toy", 0.2, ["1", "3", "7", "10"]),
    ],
)
def test_ibd(
    otter_qc: OtterQC, filename: str, threshold: float, expected_iid: list[str]
) -> None:
    """Test Inheritance By Descent (IBD)."""
    indv_out = otter_qc.ibd(filename, threshold)

    assert set(indv_out.IID) == set(
        expected_iid
    ), f"wrong IID for threshold {threshold}. \
        Expected {expected_iid}, got {indv_out.IID}"

    king_file = filename + ".kin0"
    cut_in_file = filename + ".king.cutoff.in.id"
    cut_out_file = filename + ".king.cutoff.out.id"

    assert os.path.isfile(king_file), f"File {king_file} was not created"
    assert os.path.isfile(cut_in_file), f"File {cut_in_file} was not created"
    assert os.path.isfile(cut_out_file), f"File {cut_out_file} was not created"

    os.remove(king_file)
    os.remove(cut_in_file)
    os.remove(cut_out_file)


@pytest.mark.parametrize(
    argnames=["filename", "error"],
    argvalues=[
        (1000, BeartypeCallHintParamViolation),
        ("tests/toy_qc", FileNotFoundError),
    ],
)
def test_duplicate_vars_errors(  # noqa: D103
    otter_qc: OtterQC, filename: int | str, error: Exception
) -> None:
    with pytest.raises(error):
        otter_qc.get_duplicate_vars(filename)


def test_duplicate_vars(
    otter_qc: OtterQC, filename_dup: str, duplicated_variants: list[str]
) -> None:
    """Test `OtterQC.get_duplicate_vars`."""
    dup_vars = otter_qc.get_duplicate_vars(filename_dup)

    basename = os.path.splitext(filename_dup)[0]
    path_dup_list = basename + ".dupvar"
    assert os.path.isfile(
        path_dup_list
    ), "Duplicated list of variants not created"
    assert set(dup_vars) == set(
        duplicated_variants
    ), "Duplicated variants do not match expected"

    os.remove(path_dup_list)
    os.remove(basename + ".log")


@pytest.mark.parametrize(
    argnames=["filename", "error"],
    argvalues=[
        (1000, BeartypeCallHintParamViolation),
        ("tests/toy_qc", FileNotFoundError),
    ],
)
def test_duplicate_rsid_errors(  # noqa: D103
    otter_qc: OtterQC, filename: int | str, error: Exception
) -> None:
    with pytest.raises(error):
        otter_qc.get_duplicate_rsids(filename)


def test_duplicate_rsid(
    otter_qc: OtterQC, filename_dup: str, duplicated_rsid: list[str]
) -> None:
    """Test `OtterQC.get_duplicate_rsids`."""
    dup_rsids = otter_qc.get_duplicate_rsids(filename_dup)

    basename = os.path.splitext(filename_dup)[0]
    path_dup_list = basename + ".rmdup.mismatch"
    path_rm_dup_list = basename + ".rmdup.list"
    path_log = basename + ".log"
    assert os.path.isfile(
        path_dup_list
    ), "Duplicated list of rsids not created"
    assert set(dup_rsids) == set(
        duplicated_rsid
    ), "Duplicated variants do not match expected"

    os.remove(path_dup_list)
    os.remove(path_rm_dup_list)
    os.remove(path_log)


@pytest.mark.parametrize(
    argnames=["args", "error"],
    argvalues=[
        ({"filename": 1000}, BeartypeCallHintParamViolation),
        (
            {"filename": "filename", "outpath": 1000},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "exclude_vars": [1000]},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "exclude_indvs": [1000]},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "maf": "0.05"},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "geno_miss": "0.05"},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filename": "filename", "indv_miss": "0.05"},
            BeartypeCallHintParamViolation,
        ),
        ({"filename": "tests/data/toy", "maf": 10}, ValueError),
        ({"filename": "tests/data/toy", "maf": 1}, ValueError),
        ({"filename": "tests/data/toy", "maf": 0}, ValueError),
        ({"filename": "tests/data/toy", "maf": -10}, ValueError),
        ({"filename": "tests/data/toy", "geno_miss": 10}, ValueError),
        ({"filename": "tests/data/toy", "geno_miss": -10}, ValueError),
        ({"filename": "tests/data/data", "indv_miss": 10}, ValueError),
        ({"filename": "tests/data/toy", "indv_miss": -10}, ValueError),
        ({"filename": "filename"}, FileNotFoundError),
    ],
)
def test_qc_error(  # noqa: D103
    otter_qc: OtterQC, args: dict[str, Any], error: Exception
) -> None:
    with pytest.raises(error):
        otter_qc.qc(**args)


@pytest.mark.parametrize(
    argnames=[
        "outpath",
        "exclude_vars",
        "exclude_indvs",
        "maf",
        "geno_miss",
        "indv_miss",
        "rm_vars",
        "rm_indv",
    ],
    argvalues=[
        ("test", None, None, None, None, None, [], []),
        ("test", None, None, 0.1, None, None, ["rs0002"], []),
        ("test", None, None, None, 0.1, None, ["rs0007"], []),
        ("test", None, None, None, None, 0.1, [], ["1", "2", "7"]),
        ("test", ["rs0002"], None, None, None, None, ["rs0002"], []),
        (
            "test",
            None,
            pd.DataFrame({"FID": ["FAM001"], "IID": ["1"]}),
            None,
            None,
            None,
            [],
            ["1"],
        ),
        ("test", None, None, 0.1, 0.1, None, ["rs0002", "rs0007"], []),
        ("test", None, None, 0.1, None, 0.1, ["rs0002"], ["1", "2", "7"]),
        ("test", None, None, None, 0.1, 0.1, [], ["1", "2", "7"]),
        ("test", None, None, 0.1, 0.1, 0.1, [], ["1", "2", "7"]),
    ],
)
def test_qc(  # noqa: D103
    otter_qc: OtterQC,
    filename_qc: str,
    outpath: str,
    exclude_vars: list[str] | None,
    exclude_indvs: pd.DataFrame | None,
    maf: float,
    geno_miss: float,
    indv_miss: float,
    rm_indv: list[str],
    rm_vars: list[str],
) -> None:
    otter_qc.qc(
        filename_qc,
        outpath,
        maf=maf,
        geno_miss=geno_miss,
        indv_miss=indv_miss,
        exclude_indvs=exclude_indvs,
        exclude_vars=exclude_vars,
    )

    bim_file, fam_file = outpath + ".bim", outpath + ".fam"
    bim_data = pd.read_csv(bim_file, sep=r"\s+", names=["rsid"], usecols=[1])
    fam_data = pd.read_csv(
        fam_file, sep=r"\s+", names=["fid", "iid"], usecols=[0, 1], dtype=str
    )

    assert not any(
        bim_data.rsid.isin(rm_vars)
    ), "Variants not removed as expected"
    assert not any(
        fam_data.iid.isin(rm_indv)
    ), "Individuals not removed as expected"

    os.remove("test.bim")
    os.remove("test.fam")
    os.remove("test.bed")
    os.remove("test.log")
    if exclude_indvs is not None:
        rm_indv_file = outpath + ".rmindv"
        os.remove(rm_indv_file)
    if os.path.isfile("test.irem"):
        os.remove("test.irem")
    if exclude_vars is not None:
        rm_indv_file = outpath + ".rmvars"
        os.remove(rm_indv_file)
