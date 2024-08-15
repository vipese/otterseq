"""otterseq.qc is a module used to compute Quality Control on PLINK files.

Supported functionalities are:
- Inheritance By Descent (IBD)
"""
import os
import subprocess

import ottersh


class OtterQC:
    """Class used to conduct Quality Control operations on PLINk files."""

    OTTER_SH_PATH = ottersh.__path__[0]
    IBD_SCRIPT = os.path.join(OTTER_SH_PATH, "ibd.sh")

    def __init__(self):  # noqa: D107
        pass

    def ibd(self, filename: str):
        """Compute Inheritance By Descent (IBD).

        Args:
            filename (str): Path to the file with prefix (e.g. `data/toy`,
                where the "data" folder contains a "toy.map", "toy.bed", and
                "toy.bim" file).

        Raises:
            TypeError: If `filename` is not of type str.
            FileNotFoundError: If the filename provided does not point to a `.bed`
                file.
        """
        if not isinstance(filename, str):
            raise TypeError(f"filepath not of type str. Got {type(filename)}")
        if not os.path.isfile(filename + ".bed"):
            raise FileNotFoundError(f"{filename}.bed file not found.")

        command = ["bash", self.IBD_SCRIPT, "--bfile", filename]
        subprocess.run(
            command, capture_output=True, text=True, check=False  # noqa: S603
        )
