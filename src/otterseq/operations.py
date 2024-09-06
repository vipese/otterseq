"""Operations module for performing various genetic data analyses.

This module contains the OtterOps class, which provides methods for
running logistic regression and other genetic data analysis operations.
It utilizes external tools like PLINK and interfaces with other
otterseq modules to process genetic data files and perform statistical
analyses on genomic datasets.
"""

import math
import os
import subprocess
from typing import ClassVar, Literal

import pandas as pd
import plotly.express as px
from beartype import beartype

import ottersh
from otterseq.utils import read_fam_file, read_pca_file


class OtterOps:
    """Class used to perform various genetic data analysis."""

    _OTTER_SH_PATH = ottersh.__path__[0]
    _LOGISTIC_SCRIPT = os.path.join(_OTTER_SH_PATH, "logistic_regression.sh")
    _SUFFIXES: ClassVar[list[str]] = [".bim", ".bed", ".fam"]
    _PCA_SUFFIX: Literal[".eigenvec"] = ".eigenvec"
    _LOGISTIC_SUFFIX: Literal[".assoc.logistic"] = ".assoc.logistic"

    def __init__(self) -> None:
        """Initialize the OtterOps class."""

    @beartype
    def run_logistic_regression(
        self, filename: str, outpath: str | None = None
    ) -> None:
        """Run logistic regression on genetic data using PLINK.

        This function performs logistic regression on the provided genetic data
        using PLINK. It reads the necessary PCA and phenotype files, writes
        temporary covariates and phenotype files, runs the logistic regression
        script, and then removes the temporary files.

        Args:
            filename (str): Path to the base filename (without suffix) containing
                the .bim, .bed, .fam, and .eigenvec files.
            outpath (str | None, optional): Path to the output file. If None,
                the output path is the same as `filename`. Defaults to None.

        Raises:
            FileNotFoundError: If any of the required files (.bim, .bed, .fam, .eigenvec)
                are not found.

        """
        for suffix in self._SUFFIXES:
            if not os.path.exists(filename + suffix):
                raise FileNotFoundError(f"File {filename + suffix} not found")
        if not os.path.exists(filename + self._PCA_SUFFIX):
            raise FileNotFoundError(
                f"File {filename + self._PCA_SUFFIX} not found"
            )

        outpath = filename if outpath is None else outpath

        pca_df = read_pca_file(filename)
        fam_df = read_fam_file(filename)

        # Write covariates to temporary file
        covars_df = pca_df.select(
            ["fid", "iid"]
            + [col for col in pca_df.columns if col.startswith("pc")]
        )
        temp_covars_path = f"{outpath}.covars.tsv"
        covars_df.write_csv(
            temp_covars_path, separator="\t", include_header=False
        )

        # Write phenotype file
        pheno_df = fam_df.select(["fid", "iid", "pheno"])
        temp_pheno_path = f"{outpath}.pheno.tsv"
        pheno_df.write_csv(
            temp_pheno_path, separator="\t", include_header=False
        )

        # Run logistic regression
        cmd = [
            "bash",
            self._LOGISTIC_SCRIPT,
            "--bfile",
            filename,
            "--outpath",
            outpath,
            "--pheno",
            temp_pheno_path,
            "--covar",
            temp_covars_path,
        ]
        subprocess.run(
            cmd, check=False, text=True, capture_output=True  # noqa: S603
        )

        # Remove temporary files
        os.remove(temp_covars_path)
        os.remove(temp_pheno_path)

    @beartype
    def plot_manhattan(self, filename: str) -> None:
        """Plot a Manhattan plot for GWAS results.

        Args:
            filename (str): Path to the input file containing GWAS results.
                            The file should have a '.assoc.logistic' extension.

        Raises:
            FileNotFoundError: If the input file with '.assoc.logistic' extension is not found.
        """
        if not filename.endswith(self._LOGISTIC_SUFFIX):
            filename += self._LOGISTIC_SUFFIX

        if not os.path.exists(filename):
            raise FileNotFoundError(f"File {filename} not found")

        # Read and transform data
        log_df = pd.read_csv(filename, delimiter=r"\s+", na_values=["NA"])
        log_df = log_df[(log_df.TEST == "ADD") & (~log_df.P.isna())]
        log_df = log_df[["CHR", "SNP", "BP", "P"]]
        log_df = log_df.rename(
            columns={
                "CHR": "chromosome",
                "SNP": "snp",
                "BP": "position",
                "P": "pvalue",
            }
        )

        # Compute log of p-value for plot
        log_df = log_df.sort_values(["chromosome", "position"])
        log_df["-log10(pvalue)"] = log_df.pvalue.apply(
            lambda x: -math.log10(x)
        )

        # Define alternating colors for chromosomes
        colors = ["#1f77b4", "#ff7f0e"]  # Blue and orange
        color_map = {str(i): colors[i % 2] for i in range(1, 23)}
        color_map.update(
            {"X": "#2ca02c", "Y": "#d62728"}
        )  # Green for X, Red for Y

        # Plot Manhattan plot
        fig = px.scatter(
            log_df,
            x="position",
            y="-log10(pvalue)",
            color="chromosome",
            title="Manhattan Plot",
            labels={"pvalue": "p-value"},
        )
        fig.show()
