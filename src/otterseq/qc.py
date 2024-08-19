"""otterseq.qc is a module used to compute Quality Control on PLINK files.

Supported functionalities are:
- Inheritance By Descent (IBD)
"""
import os
import subprocess

import pandas as pd

import ottersh


class OtterQC:
    """Class used to conduct Quality Control operations on PLINk files."""

    OTTER_SH_PATH = ottersh.__path__[0]
    IBD_SCRIPT = os.path.join(OTTER_SH_PATH, "ibd.sh")

    def __init__(self):  # noqa: D107
        pass

    def ibd(
        self, filename: str, threshold: float | int = 0.25
    ) -> pd.DataFrame:
        """Compute Inheritance By Descent (IBD).

        This function leverages PLINK2.0 to compute IBD on a provided filename
        path (path to file with prefix and no suffix, such as "data/toy"). The output
        is written in the same directory as the input file ("data/toy.kin0").

        Args:
            filename (str): Path to the file with prefix (e.g. `data/toy`,
                where the "data" folder contains a "toy.map", "toy.bed", and
                "toy.bim" file).
            threshold (float | int): Kinship threshold used to filter out
                individuals.

        Raises:
            TypeError: If `filename` is not of type str.
            TypeError: If `threshold` is not of type int | float.
            ValueError: If `threshold` is not in between 0 and 1 (included).
        """
        if not isinstance(filename, str):
            raise TypeError(f"filepath not of type str. Got {type(filename)}")
        if not isinstance(threshold, float | int):
            raise TypeError(
                f"threshold not of type int or float. Got {type(threshold)}"
            )
        if not os.path.isfile(filename + ".bed"):
            raise FileNotFoundError(f"{filename}.bed file not found.")
        if not 0 <= threshold < 1:
            raise ValueError(f"threshold not in range (0,1). Got {threshold}")

        # Compute IBD
        command = [
            "bash",
            self.IBD_SCRIPT,
            "--bfile",
            filename,
            "--threshold",
            str(threshold),
        ]
        subprocess.run(
            command, capture_output=True, text=True, check=False  # noqa: S603
        )

        # Return filtered out patients
        indv_out = pd.read_csv(
            filename + ".king.cutoff.out.id",
            sep=r"\s+",
            header=0,
            names=["FID", "IID"],
            dtype={"FID": str, "IID": str},
        )

        return indv_out
