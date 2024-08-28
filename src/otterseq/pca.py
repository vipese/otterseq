"""."""
import logging
import os
import subprocess
from typing import ClassVar

import plotly.express as px
import plotly.graph_objects as go
import polars as pl
from beartype import beartype

import ottersh


class OtterPCA:
    """Class used to compute PCA and plot results."""

    _OTTER_SH_PATH = ottersh.__path__[0]
    _PCA_SCRIPT = os.path.join(_OTTER_SH_PATH, "pca.sh")
    _SUFFIXES: ClassVar[list[str]] = [".bim", ".bed", ".fam"]

    def __init__(self) -> None:  # noqa: D107
        pass

    @beartype
    def pca(
        self,
        filepath: str,
        outpath: str | None = None,
        exclude_hla: bool = True,
        n_pcs: int = 20,
    ) -> None:
        """Compute Principal Component Analysis (PCA) on a PLINK binary file.

        Args:
            filepath (str): Path to the binary file, without suffix (e.g. data/toy)
            outpath (str | None, optional): Path to the output file. If None,
                is the same as `filepath`. Defaults to None.
            exclude_hla (bool, optional): True to exclude HLA region from the PCA.
                Defaults to True.
            n_pcs (int, optional): Number of Principal Components. Defaults to 20.

        Raises:
            FileNotFoundError: If any of the PLINK binary files are missing.
            ValueError: If the number of PCs `n_pcs` is lesser than 0.
        """
        for suf in self._SUFFIXES:
            if not os.path.isfile(filepath + suf):
                raise FileNotFoundError(f"{filepath} not found.")
        if n_pcs <= 0:
            raise ValueError(
                f"Number of PCs cannot be lower or equal to 0. Got {n_pcs}"
            )

        outpath = outpath or filepath

        command = [
            "bash",
            self._PCA_SCRIPT,
            "--bfile",
            filepath,
            "--outpath",
            outpath,
            "--exclude-hla",
            str(exclude_hla),
            "--pcs",
            str(n_pcs),
        ]
        subprocess.run(
            command, text=True, capture_output=True, check=False  # noqa: S603
        )

    @beartype
    def plot_pca(self, filename: str, plot: bool = True) -> go.Figure:
        """Plot PCA results.

        Args:
            filename (str | os.PathLike[str]): Path to the folder with files,
                without suffix (e.g, "data/toy")
            plot (bool, optional): If True, plots the PCA. Defaults to True.

        Raises:
            FileNotFoundError: If `filename`.eigenvec does not exist.

        Returns:
            go.Figure: plotly Figure which can be plotted with `.show()`
        """
        pca_filepath = filename + ".eigenvec"
        if not os.path.isfile(pca_filepath):
            raise FileNotFoundError(f"{pca_filepath} was not found.")

        fam_filepath = filename + ".fam"
        if not os.path.isfile(fam_filepath):
            logging.warning(
                FileNotFoundError(
                    f"{fam_filepath} file not found. Phenotypes will not be plotted"
                )
            )
            fam_df = None
        else:
            fam_df = pl.read_csv(
                fam_filepath,
                has_header=False,
                separator=" ",
                new_columns=[
                    "fid",
                    "iid",
                    "father_id",
                    "mother_id",
                    "sex",
                    "pheno",
                ],
            )
            fam_df = fam_df.with_columns(
                pl.when(pl.col("sex") == 2)
                .then(pl.lit("female"))
                .when(pl.col("sex") == 1)
                .then(pl.lit("male"))
                .otherwise(pl.lit("unknown"))
                .alias("sex_str")
            )
            fam_df = fam_df.with_columns(
                pl.when(pl.col("pheno") == 2)
                .then(pl.lit("control"))
                .when(pl.col("pheno") == 1)
                .then(pl.lit("case"))
                .otherwise(pl.lit("unknown"))
                .alias("pheno_str")
            )

        pca_df = pl.read_csv(pca_filepath, has_header=False, separator=" ")
        pca_df.columns = ["fid", "iid"] + [
            f"pc{i}" for i in range(len(pca_df.columns) - 2)
        ]

        if fam_df is not None:
            pca_df = pca_df.join(
                other=fam_df[["fid", "iid", "sex_str", "pheno_str"]],
                on=["fid", "iid"],
                how="left",
            )
            fig = px.scatter(
                data_frame=pca_df,
                x="pc1",
                y="pc2",
                color="pheno_str",
                symbol="sex_str",
                hover_data=["fid", "iid"],
                title=f"Principal Component Analysis of {pca_filepath}",
                labels={
                    "pc1": "PC1",
                    "pc2": "PC2",
                    "fid": "FID",
                    "iid": "IID",
                    "pheno_str": "Phenotype",
                    "sex_str": "Sex",
                },
            )
        else:
            fig = px.scatter(
                data_frame=pca_df,
                x="pc1",
                y="pc2",
                hover_data=["fid", "iid"],
                title=f"Principal Component Analysis of {pca_filepath}",
                labels={
                    "pc1": "PC1",
                    "pc2": "PC2",
                    "fid": "FID",
                    "iid": "IID",
                },
            )
        if plot:
            fig.show()
        return fig
