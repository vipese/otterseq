"""."""

import subprocess
import ottersh
import logging
import os


class OtterSNP:
    def __init__(self) -> None:
        """."""

    def binarize_files(self, filepath: str, outpath: str) -> None:
        """Binarize PLINK raw files to binary PLINK files.

        Args:
            filepath (str): Path to files. Accepts a path to a file without suffix,
                or a directory to multiple files.
            outpath (str): Output path, without suffix.
        """

        if not isinstance(filepath, str):
            raise TypeError(f"filepath must be of type str. Got {type(filepath)}")
        if not isinstance(outpath, str):
            raise TypeError(f"filepath must be of type str. Got {type(filepath)}")

        # Handle file path to directory or file (including with suffix .ped or .map)
        ped_files: list[str]
        if os.path.isdir(filepath):
            ped_files = [
                os.path.join(filepath, f)
                for f in os.listdir(filepath)
                if f.endswith(".ped")
            ]
        else:
            if filepath.endswith(".ped"):
                ped_files = [filepath]
            else:
                ped_files = [filepath + ".ped"]
            if not os.path.isfile(ped_files[0]):
                raise FileNotFoundError(
                    f"Provided path did not point to an existing file. Got {ped_files[0]}"
                )
        file_list = [f.split(".ped")[-2] for f in ped_files]

        # Handle outfile
        if not os.path.isdir(outpath):
            os.makedirs(outpath)

        # Verbose
        logging.info("Binarizing files: \n")
        logging.info(f"{f} \n" for f in file_list)

        # Call bash script to binarize files in Data/GWAS
        script_path = os.path.join(ottersh.__path__[0], "binarize.sh")
        command = ["bash", script_path, "--files"] + file_list + ["--outpath", outpath]
        subprocess.run(command, capture_output=True, text=True)
