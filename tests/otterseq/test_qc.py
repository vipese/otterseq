"""OtterQC unit tests."""

import os

import pytest

from otterseq.qc import OtterQC


@pytest.mark.parametrize(argnames="filename", argvalues=[1000])
def test_ibd_type_error(  # noqa: D103
    otter_qc: OtterQC, filename: str | int
) -> None:
    with pytest.raises(TypeError):
        otter_qc.ibd(filename)


def test_ibd_no_bed_file(  # noqa: D103
    otter_qc: OtterQC, no_filename: str
) -> None:
    with pytest.raises(FileNotFoundError):
        otter_qc.ibd(filename=no_filename)


def test_ibd(otter_qc: OtterQC, filename: str) -> None:
    """Test Inheritance By Descent (IBD)."""
    otter_qc.ibd(filename)
    outfile_path = filename + ".kin0"
    assert os.path.isfile(outfile_path), f"File {outfile_path} was not created"
    os.remove(outfile_path)
