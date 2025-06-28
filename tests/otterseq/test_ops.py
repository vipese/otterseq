"""Tests for the OtterOps class in the otterseq.operations module.

The tests cover the following functionalities:
- Running logistic regression on genetic data using PLINK.
- Handling missing files when running logistic regression.
- Running logistic regression with a custom output path.

Fixtures:
- otter_ops: Provides an instance of the OtterOps class.
- mock_data_dir: Creates a temporary directory for mock data files.

Test Functions:
- test_run_logistic_regression: Tests the run_logistic_regression method with mock data files.
- test_run_logistic_regression_file_not_found: Tests the run_logistic_regression method when required files are missing.
- test_run_logistic_regression_custom_outpath: Tests the run_logistic_regression method with a custom output path.
"""


import os
from unittest.mock import patch

import pandas as pd
import pytest

from otterseq.operations import OtterOps


def test_run_logistic_regression_file_not_found(
    otter_ops: OtterOps, no_files_directory: str
) -> None:
    """Test the run_logistic_regression method when required files are missing."""
    base_filename = no_files_directory + "nonexistent"
    with pytest.raises(FileNotFoundError):
        otter_ops.run_logistic_regression(base_filename)


@pytest.mark.parametrize(
    argnames="outpath", argvalues=[None, "data_pca/custom_output"]
)
def test_run_logistic_regression(
    otter_ops: OtterOps, filename_pca: str, outpath: str
) -> None:
    """Test the run_logistic_regression method with mock data files."""
    # mock the read functions and subprocess.run
    with patch("otterseq.operations.read_pca_file") as mock_read_pca, patch(
        "otterseq.operations.read_fam_file"
    ) as mock_read_fam, patch("subprocess.run") as mock_run, patch(
        "os.remove"
    ) as mock_remove:

        # set up mock return values
        mock_read_pca.return_value.columns = [
            "fid",
            "iid",
            "pc1",
            "pc2",
            "pc3",
            "pc4",
            "pc5",
        ]
        mock_read_fam.return_value.columns = [
            "fid",
            "iid",
            "father",
            "mother",
            "sex",
            "pheno",
        ]

        # call the function
        otter_ops.run_logistic_regression(filename_pca, outpath=outpath)

        # assertions
        mock_read_pca.assert_called_once_with(filename_pca)
        mock_read_fam.assert_called_once_with(filename_pca)
        mock_run.assert_called_once()
        assert mock_remove.call_count == 2


@pytest.mark.parametrize(
    argnames="outpath", argvalues=[None, "tests/data_pca/custom_output"]
)
def test_run_logistic_regression_result(
    otter_ops: OtterOps, filename_pca: str, outpath: str
) -> None:
    """Test the result of the run_logistic_regression method."""
    otter_ops.run_logistic_regression(filename_pca, outpath=outpath)

    outpath = outpath or filename_pca
    assert os.path.exists(outpath + ".assoc.logistic")

    log_df = pd.read_csv(outpath + ".assoc.logistic", sep=r"\s+", header=0)
    assert (
        not log_df.empty
    ), "The logistic regression results DataFrame should not be empty"

    required_columns = [
        "CHR",
        "SNP",
        "BP",
        "A1",
        "TEST",
        "NMISS",
        "OR",
        "STAT",
        "P",
    ]
    for col in required_columns:
        assert (
            col in log_df.columns
        ), f"Error: Column '{col}' is missing from the file."

    os.remove(outpath + ".assoc.logistic")


def test_plot_manhattan(otter_ops: OtterOps, filename_log: str) -> None:
    """Test plot manhattan."""
    otter_ops.plot_manhattan(filename_log)
