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

from unittest.mock import patch

import pandas as pd
import polars as pl
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
    with (
        patch("otterseq.operations.read_pca_file") as mock_read_pca,
        patch("otterseq.operations.read_fam_file") as mock_read_fam,
        patch("subprocess.run") as mock_run,
        patch("os.remove") as mock_remove,
    ):

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
    # Mock subprocess.run to avoid running actual PLINK commands
    with patch("subprocess.run") as mock_run:
        # Mock the read functions to return test data
        mock_pca_df = pl.DataFrame(
            {
                "fid": ["FAM001", "FAM001"],
                "iid": ["1", "2"],
                "pc1": [0.1, 0.2],
                "pc2": [0.3, 0.4],
            }
        )
        mock_fam_df = pl.DataFrame(
            {
                "fid": ["FAM001", "FAM001"],
                "iid": ["1", "2"],
                "father_id": [0, 0],
                "mother_id": [0, 0],
                "sex": [1, 2],
                "pheno": [1, 2],
            }
        )

        with (
            patch(
                "otterseq.operations.read_pca_file", return_value=mock_pca_df
            ) as mock_read_pca,
            patch(
                "otterseq.operations.read_fam_file", return_value=mock_fam_df
            ) as mock_read_fam,
            patch("os.remove") as mock_remove,
        ):
            # Call the function
            otter_ops.run_logistic_regression(filename_pca, outpath=outpath)

            # Verify subprocess was called
            mock_run.assert_called_once()

            # Verify the read functions were called
            mock_read_pca.assert_called_once_with(filename_pca)
            mock_read_fam.assert_called_once_with(filename_pca)

            # Verify temporary files were cleaned up
            assert mock_remove.call_count == 2

            # Verify the command structure
            call_args = mock_run.call_args
            command = call_args[0][0]
            assert command[0] == "bash"
            assert "--bfile" in command
            assert "--outpath" in command
            assert "--pheno" in command
            assert "--covar" in command


def test_plot_manhattan(otter_ops: OtterOps, filename_log: str) -> None:
    """Test plot manhattan."""
    # Mock pandas.read_csv to return test data
    mock_df = pd.DataFrame(
        {
            "CHR": [1, 1, 2, 2],
            "SNP": ["rs0001", "rs0002", "rs0003", "rs0004"],
            "BP": [100100, 100200, 200100, 200200],
            "TEST": ["ADD", "ADD", "ADD", "ADD"],
            "P": [0.001, 0.01, 0.05, 0.1],
        }
    )

    with (
        patch("pandas.read_csv", return_value=mock_df) as mock_read_csv,
        patch("plotly.express.scatter") as mock_scatter,
    ):
        mock_fig = mock_scatter.return_value
        mock_fig.show = lambda: None  # Mock the show method

        # Call the function
        otter_ops.plot_manhattan(filename_log)

        # Verify the file was read
        mock_read_csv.assert_called_once()

        # Verify the plot was created
        mock_scatter.assert_called_once()

        # Verify the plot data structure
        called_args = mock_scatter.call_args
        plot_df = called_args[0][0]  # First positional argument
        assert "chromosome" in plot_df.columns
        assert "position" in plot_df.columns
        assert "pvalue" in plot_df.columns
        assert "-log10(pvalue)" in plot_df.columns
