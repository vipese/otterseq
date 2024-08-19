"""OtterQC unit tests."""

import os

import pytest

from otterseq.qc import OtterQC


@pytest.mark.parametrize(
    argnames=["filename", "threshold"],
    argvalues=[(1000, 0.25), ("filename", "threshold")],
)
def test_ibd_type_error(  # noqa: D103
    otter_qc: OtterQC, filename: int | str, threshold: float | str
) -> None:
    with pytest.raises(TypeError):
        otter_qc.ibd(filename, threshold)


def test_ibd_no_bed_file(  # noqa: D103
    otter_qc: OtterQC, no_filename: str
) -> None:
    with pytest.raises(FileNotFoundError):
        otter_qc.ibd(filename=no_filename)


@pytest.mark.parametrize(argnames="threshold", argvalues=[-10, 1, 10])
def test_ibd_threshold_error(  # noqa: D103
    otter_qc: OtterQC, filename: str, threshold: float
):
    with pytest.raises(ValueError, match="threshold not in range"):
        otter_qc.ibd(filename, threshold)


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
