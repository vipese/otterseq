"""."""

import logging
import os
import subprocess

import ottersh
from otterseq.errors import MultiAllelicError


class OtterSNP:
    """Class used to handle PLINK files.

    OtterSNP provides all the necessary to manipulate, and convert PLINK files.
    """

    def __init__(self) -> None:  # noqa: D107
        pass

    def _read_snp_id_bim(self, filepath: str) -> list[str]:
        """Read the rsIDs from a bim file.

        Args:
            filepath (str): Path to the .bim file.

        Returns:
            list[str]: List of rsIDs.
        """
        with open(filepath) as in_file:
            snps: list[str] = [row.split("\t")[1] for row in in_file]
        return snps

    def check_multi_allelic(self, rs_ids: list[str]) -> None:
        """Check if there are multi-allelic (duplicated rsIDs) variants.

        Args:
            rs_ids (list[str]): List of rsIDs.

        Raises:
            MultiAllelicError: If there are duplicated rsIDs.
        """
        if len(rs_ids) > len(set(rs_ids)):
            raise MultiAllelicError("Multi-allelic variants found in file.")

    def binarize_files(self, filepath: str, outpath: str) -> None:
        """Binarize PLINK raw files to binary PLINK files.

        Args:
            filepath (str): Path to files. Accepts a path to a file without suffix,
                or a directory to multiple files.
            outpath (str): Output path, without suffix.

        Raises:
            TypeError: If filepath is not of type str.
            TypeError: If outpath if not of type str.
            FileNotFoundError: If the .ped file if not found.
        """
        if not isinstance(filepath, str):
            raise TypeError(
                f"filepath must be of type str. Got {type(filepath)}"
            )
        if not isinstance(outpath, str):
            raise TypeError(
                f"filepath must be of type str. Got {type(filepath)}"
            )

        # Handle file path to directory or file (including with suffix .ped)
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
        for file in file_list:
            command = [
                "bash",
                script_path,
                "--file",
                file,
                "--outpath",
                outpath,
            ]
            subprocess.run(
                command,  # noqa: S603
                capture_output=True,
                text=True,
                check=False,
            )

    def get_common_snp(
        self, filepath: str, write: bool = False, outpath: str | None = None
    ) -> list[str]:
        """Parse common SNPs of .bim files.

        Given a path to a directory, `get_common_snp` reads all the .bim files located
        inside the directory and returns the common SNPs.

        **Note**: If any of the .bim files contains duplicated rsIDs, the function
        raises a `MultiAllelicError` that stops the execution.

        Args:
            filepath (str): Path to directory containing .bim files
            write (bool, optional): True to write the common SNPs into a file.
                Defaults to False.
            outpath (str | None, optional): Path to the directory where the output
                should be written. Requires write = True. Default

        Raises:
            TypeError: If filepath is not of type str.
            TypeError: If write is not of type bool.
            TypeError: If outpath is not of type str
            ValueError: If write is False and outpath is not None.
            FileNotFoundError: If there are no .bim files in filepath.

        Returns:
            list[str]: _description_
        """
        # Enforce types
        if not isinstance(filepath, str):
            raise TypeError(
                f"filepath must be of type str. Got {type(filepath)}"
            )
        if not isinstance(write, bool):
            raise TypeError(f"write must be of type bool. Got {type(write)}")
        if outpath and not isinstance(outpath, str):
            raise TypeError(
                f"outpath must be of type str. Got {type(outpath)}"
            )
        if write is True and outpath is None:
            raise ValueError(
                "Write was set to True but no outpath was provided"
            )

        # Parse .bim files
        bim_files = [
            os.path.join(filepath, f)
            for f in os.listdir(filepath)
            if f.endswith(".bim")
        ]
        if not bim_files:
            raise FileNotFoundError(
                "No .bim files found in the provided filepath."
            )

        # Get first file to compute intersection iteratively
        total_snps = self._read_snp_id_bim(bim_files[0])
        self.check_multi_allelic(total_snps)
        total_snps = set(total_snps)

        # For the rest of the files, intersect with first file
        n_files = len(bim_files)
        for i in range(1, n_files):
            file_snps = self._read_snp_id_bim(bim_files[i])
            self.check_multi_allelic(file_snps)
            total_snps.intersection_update(set(file_snps))
        total_snps = list(total_snps)
        total_snps.sort()

        # Write to file
        if write is True and outpath is not None:
            with open(
                os.path.join(outpath, "common_snps.txt"), "w"
            ) as out_file:
                for snp in total_snps:
                    out_file.write(f"{snp}\n")

        return total_snps
