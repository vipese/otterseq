"""otterseq.snp is a library used to manipulate binary PLINK files.

Supported functionalities are:
- Binarizing PLINK text files.
- Compute common SNPs across .bim files.
- Merging binary PLINK files.

"""

import logging
import os
import subprocess

import ottersh
from otterseq.errors import MultiAllelicError


class OtterSNP:
    """Class used to handle PLINK files.

    OtterSNP provides all the necessary to manipulate, and convert PLINK files.
    """

    OTTER_SH_PATH = ottersh.__path__[0]
    BINARIZE_SCRIPT = os.path.join(OTTER_SH_PATH, "binarize.sh")
    MERGE_SCRIPT = os.path.join(OTTER_SH_PATH, "merge_files.sh")

    def __init__(self) -> None:  # noqa: D107
        pass

    def _read_snp_id_bim(self, filepath: str) -> list[str]:
        """Read the rsIDs from a .pvar file.

        Args:
            filepath (str): Path to the .pvar file.

        Returns:
            list[str]: List of rsIDs.
        """
        with open(filepath) as in_file:
            snps: list[str] = [row.split("\t")[1] for row in in_file]
        snps = snps[1:]  # Skip header
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
        script_path = self.BINARIZE_SCRIPT
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
            list[str]: List of common SNPs.
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
        var_files = [
            os.path.join(filepath, f)
            for f in os.listdir(filepath)
            if f.endswith(".bim")
        ]
        if not var_files:
            raise FileNotFoundError(
                "No .bim files found in the provided filepath."
            )

        # Get first file to compute intersection iteratively
        total_snps = self._read_snp_id_bim(var_files[0])
        self.check_multi_allelic(total_snps)
        total_snps = set(total_snps)

        # For the rest of the files, intersect with first file
        n_files = len(var_files)
        for i in range(1, n_files):
            file_snps = self._read_snp_id_bim(var_files[i])
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

    def merge_files(
        self,
        filepath: str,
        outpath: str | None = None,
        prefix: str | None = None,
    ) -> None:
        """Merge binary PLINK files.

        Given the path to a directory with binary PLINK1.9 files, compute
        the list of files to merge, and pass it to plink to merge files into
        a single binary PLINK file.

        Args:
            filepath (str): Path to directory with files to merge.
            outpath (str | None, optional): Path to directory to store the
                outputs. If None, uses `filepath`. Defaults to None.
            prefix (str | None, optional): Prefix used to save the binary PLINK
                files after merging. If None, defaults to "merged_snps".
                Defaults to None.

        Raises:
            TypeError: If `filepath` is not of type str.
            TypeError: If `outpath` is not None and not of type str.
            TypeError: If `prefix` is not None and not of type str.
            FileNotFoundError: If no `.bed` files were found in `filepath`.
        """
        if not isinstance(filepath, str):
            raise TypeError(f"filepath not of type str. Got {type(filepath)}")
        if outpath is not None and not isinstance(outpath, str):
            raise TypeError(f"outpath not of type str. Got {type(outpath)}")
        if prefix is not None and not isinstance(prefix, str):
            raise TypeError(f"prefix not of type str. Got {type(prefix)}")

        # Initialize paths
        outpath = outpath or filepath
        if not os.path.isdir(outpath):
            os.makedirs(outpath)
        merge_list_path = os.path.join(outpath, "merge_list.txt")

        # Parse pgen files
        pgen_files = [f for f in os.listdir(filepath) if f.endswith(".bed")]
        if not pgen_files:
            raise FileNotFoundError("No .bed files found in filepath")
        files_prefix = [
            os.path.join(filepath, f.split(".bed")[-2]) for f in pgen_files
        ]

        # Write merge list file
        with open(merge_list_path, "w") as out_file:
            for f in files_prefix:
                out_file.write(f"{f}\n")

        # Run script
        script_path = self.MERGE_SCRIPT
        command = [
            "bash",
            script_path,
            "--merge-file",
            merge_list_path,
            "--outpath",
            outpath,
        ]
        command = (
            [*command, "--prefix", prefix] if prefix is not None else command
        )
        subprocess.run(
            command, capture_output=True, text=True, check=False  # noqa: S603
        )
