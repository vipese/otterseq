"""otterseq.snp is a module used to manipulate binary PLINK files.

Supported functionalities are:
- Binarizing PLINK text files.
- Compute common SNPs across .bim files.
- Merging binary PLINK files.

"""

import logging
import os
import subprocess

import polars as pl
from beartype import beartype

import ottersh
from otterseq.utils import check_multi_allelic


class OtterSNP:
    """Class used to handle PLINK files.

    OtterSNP provides all the necessary to manipulate, and convert PLINK files.
    """

    OTTER_SH_PATH = ottersh.__path__[0]
    BINARIZE_SCRIPT = os.path.join(OTTER_SH_PATH, "binarize.sh")
    MERGE_SCRIPT = os.path.join(OTTER_SH_PATH, "merge_files.sh")

    @beartype
    def _read_snp_id_bim(self, filepath: str) -> list[str]:
        """Read the rsIDs from a .bim file.

        Args:
            filepath (str): Path to the .bim file.

        Returns:
            list[str]: List of rsIDs.
        """
        try:
            # Check if file is empty
            if os.path.getsize(filepath) == 0:
                return []

            # Read .bim file with polars (tab-separated, no header)
            df = pl.read_csv(
                filepath,
                has_header=False,
                separator="\t",
                new_columns=["chrom", "rsid", "cm", "pos", "a1", "a2"],
            )
            return df["rsid"].to_list()
        except Exception as e:
            logging.error(f"Error reading {filepath}: {e}")
            raise

    @beartype
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
            subprocess.run(  # noqa: S603
                command,
                capture_output=True,
                text=True,
                check=False,
            )

    @beartype
    def _validate_get_common_snp_params(
        self, write: bool, outpath: str | None
    ) -> None:
        """Validate parameters for get_common_snp method.

        Args:
            write (bool): Whether to write output to file.
            outpath (str | None): Output path for writing.

        Raises:
            ValueError: If write is True but outpath is None.
        """
        if write is True and outpath is None:
            raise ValueError(
                "Write was set to True but no outpath was provided"
            )

    @beartype
    def _process_bim_files(self, filepath: str) -> list[set[str]]:
        """Process .bim files and extract SNP sets.

        Args:
            filepath (str): Path to directory containing .bim files.

        Returns:
            list[set[str]]: List of SNP ID sets from each file.

        Raises:
            FileNotFoundError: If no .bim files found in directory.
            Exception: If error reading any file.
        """
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

        # Read all .bim files and extract SNP IDs
        snp_sets: list[set[str]] = []
        for bim_file in bim_files:
            try:
                if os.path.getsize(bim_file) == 0:
                    logging.warning(f"Empty .bim file: {bim_file}.")
                    snp_ids = []
                else:
                    df = pl.read_csv(
                        bim_file,
                        has_header=False,
                        separator="\t",
                        new_columns=["chrom", "rsid", "cm", "pos", "a1", "a2"],
                    )
                    snp_ids = df["rsid"].to_list()

                if not snp_ids:
                    logging.warning(f"No SNPs found in {bim_file}.")
                    continue

                # Check for multi-allelic variants
                check_multi_allelic(snp_ids)

                # Add SNP IDs to set
                snp_sets.append(set(snp_ids))

            except Exception as e:
                logging.error(f"Error reading {bim_file}: {e}")
                raise

        return snp_sets

    @beartype
    def _find_common_snps(self, snp_sets: list[set[str]]) -> list[str]:
        """Find common SNPs across all SNP sets.

        Args:
            snp_sets (list[set[str]]): List of SNP ID sets.

        Returns:
            list[str]: Sorted list of common SNPs.
        """
        if not snp_sets:
            return []

        # Start with the first set and intersect with all others
        common_snps = snp_sets[0]
        for snp_set in snp_sets[1:]:
            common_snps = common_snps.intersection(snp_set)

            # Early termination if no common SNPs remain
            if not common_snps:
                logging.warning("No common SNPs found across all files")
                break

        return sorted(common_snps)

    @beartype
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
                should be written. Requires write = True. Default to None.

        Raises:
            FileNotFoundError: If there are no .bim files in filepath.

        Returns:
            list[str]: List of common SNPs.
        """
        # Enforce logic
        if write is True and outpath is None:
            raise ValueError(
                "Write was set to True but no outpath was provided"
            )

        # Process .bim files and extract SNP sets
        snp_sets = self._process_bim_files(filepath)

        # Find common SNPs
        common_snps_list = self._find_common_snps(snp_sets)

        # Write to file if requested
        if write is True and outpath is not None:
            with open(
                os.path.join(outpath, "common_snps.txt"), "w"
            ) as out_file:
                for snp in common_snps_list:
                    out_file.write(f"{snp}\n")

        return common_snps_list

    @beartype
    def merge_files(
        self,
        filepath: str,
        outpath: str | None = None,
        prefix: str = "merged_snps",
        only_common: bool = True,
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
            only_common (bool): True to extract only common variants across files.
                Defaults to True.

        Raises:
            FileNotFoundError: If no `.bed` files were found in `filepath`.
        """
        # Handle outpath
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

        # Write common snps to file
        if only_common is True:
            _ = self.get_common_snp(
                filepath=filepath, outpath=outpath, write=True
            )
        common_snps_path = os.path.join(outpath, "common_snps.txt")

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
        command = (
            [*command, "--extract", common_snps_path]
            if only_common is True
            else command
        )
        subprocess.run(  # noqa: S603
            command, capture_output=True, text=True, check=False
        )
