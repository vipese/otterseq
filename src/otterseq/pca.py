"""."""
import os
import subprocess
from typing import ClassVar

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
