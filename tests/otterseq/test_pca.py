"""Unit test for `OtterPCA`."""
import os
from typing import Any

import pandas as pd
import plotly.graph_objects as go
import pytest
from beartype.roar import BeartypeCallHintParamViolation

from otterseq.pca import OtterPCA


@pytest.mark.parametrize(
    argnames=["args", "error"],
    argvalues=[
        ({"filepath": 100}, BeartypeCallHintParamViolation),
        (
            {"filepath": "filepath", "outpath": 100},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filepath": "filepath", "exclude_hla": "exlude_hla"},
            BeartypeCallHintParamViolation,
        ),
        (
            {"filepath": "filepath", "n_pcs": "n_pcs"},
            BeartypeCallHintParamViolation,
        ),
        ({"filepath": "filepath"}, FileNotFoundError),
        ({"filepath": "tests/data_qc/toy", "n_pcs": -10}, ValueError),
    ],
)
def test_pca_errors(  # noqa: D103
    otter_pca: OtterPCA, args: dict[str, Any], error: Exception
) -> None:
    with pytest.raises(error):
        otter_pca.pca(**args)


@pytest.mark.parametrize(
    argnames=["n_pcs", "expected_n_pcs"],
    argvalues=[
        (2, 2),
        (5, 5),
        (10, 5),
    ],
)
def test_pca(
    otter_pca: OtterPCA, filename: str, n_pcs: int, expected_n_pcs: int
) -> None:
    """Test `OtterPCA.pca`."""
    otter_pca.pca(filepath=filename, exclude_hla=False, n_pcs=n_pcs)

    eigenval_path = filename + ".eigenval"
    eigenvec_path = filename + ".eigenvec"
    eigenval = pd.read_csv(
        eigenval_path, names=["eigenval"], header=None, sep=r"\s+"
    )
    eigenvec = pd.read_csv(
        eigenvec_path,
        names=["FID"] + [f"PC{i}" for i in range(expected_n_pcs)],
        header=None,
        sep=r"\s+",
    )

    assert (
        len(eigenval) == expected_n_pcs
    ), f"The number of eigenvalues does not match expected ({expected_n_pcs}). Got {len(eigenval)}"
    assert (
        len(eigenvec) == 10
    ), f"The number of samples does not match expected (10). Got {len(eigenvec)}"
    assert (
        len(eigenvec.columns) == expected_n_pcs + 1
    ), f"Number of PCs does not match expected ({expected_n_pcs + 1}). Got {len(eigenvec.columns)}"

    os.remove(eigenval_path)
    os.remove(eigenvec_path)


def test_pca_plot(otter_pca: OtterPCA, filename_pca: str) -> None:
    """."""
    fig = otter_pca.plot_pca(filename=filename_pca, plot=False)

    assert isinstance(
        fig, go.Figure
    ), f"Output not of type go.Figure. Got{type(fig)}"
