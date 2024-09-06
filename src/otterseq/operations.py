"""."""

import os
import subprocess
from typing import ClassVar, Literal

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
