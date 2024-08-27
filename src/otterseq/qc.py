"""otterseq.qc is a module used to compute Quality Control on PLINK files.

Supported functionalities are:
- Inheritance By Descent (IBD)
- Computation of duplicated variants by coordinates and allele codes.
- Computation of duplicated variant rsIDs.
"""
import os
import subprocess
from typing import ClassVar

import pandas as pd
from beartype import beartype

import ottersh


class OtterQC:
    """Class used to conduct Quality Control operations on PLINk files."""

    _OTTER_SH_PATH = ottersh.__path__[0]
    _IBD_SCRIPT = os.path.join(_OTTER_SH_PATH, "ibd.sh")
    _DUP_VARS = os.path.join(_OTTER_SH_PATH, "duplicate_vars.sh")
    _DUP_RSID = os.path.join(_OTTER_SH_PATH, "duplicate_rsid.sh")
    _QC_SCRIPT = os.path.join(_OTTER_SH_PATH, "qc.sh")
    _SUFFIXES: ClassVar[list[str]] = [".bim", ".bed", ".fam"]

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
            FileNotFoundError: If `.bed` are found in the provided `filename`.
            ValueError: If `threshold` is not in between 0 and 1 (included).

        Returns:
            pd.DataFrame: DataFrame containing the FID and IID of the excluded
                individuals.
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
            self._IBD_SCRIPT,
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

    @beartype
    def get_duplicate_vars(self, filename: str) -> list[str]:
        """Get the duplicated variants based on coordinates and allele codes.

        For instance, given a `.bim` file as follows:

        ```text
        1	rs0001	0	100100	G	A
        1	rs0001_dup1	0	100100	G	A
        1	rs0001_dup2	0	100100	G	A
        1	rs0002	0	100200	A	T
        1	rs0002_dup1	0	100200	A	T
        1	rs0003	0	100300	C	G
        1	rs0001	0	100300	T	G
        1	rs0004	0	100400	T	C
        1	rs0005	0	100500	T	A
        1	rs0002	0	100500	T	G
        ```

        This function creates the following `.dupvar` file:

        ```text
        CHR	POS	ALLELES	IDS
        1	100100	A,G	rs0001_dup1 rs0001_dup2
        1	100200	A,T	rs0002_dup1
        ```

        Indicating that `rs0001_dup1` and `rs0002_dup1` are duplicates. Then,
        it returns ["rs0001_dup1", "rs0002_dup1"]

        Args:
            filename (str): Path to the .bim file. Accepts the file without
                suffix (e.g, "toy"), or with suffix ("toy.bim")

        Returns:
            list[str]: List with duplicated variants rsIDs.
        """
        if not filename.endswith(".bim"):
            filename += ".bim"
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} not found.")

        basename = os.path.splitext(filename)[0]
        command = ["bash", self._DUP_VARS, "--bfile", basename]
        subprocess.run(
            command, capture_output=True, text=True, check=False  # noqa: S603
        )

        dup_vars = pd.read_csv(
            basename + ".dupvar",
            usecols=[3],
            names=["rsid"],
            sep="\t",
            header=0,
        )
        dup_vars = dup_vars.rsid.str.split().explode().to_list()
        return dup_vars

    @beartype
    def get_duplicate_rsids(self, filename: str) -> list[str]:
        """Get duplicated variants rsIDs.

        This function returns the variants that have duplicated rsIDs. Given the
        following `.bim` file:

        ```text
        1	rs0001	0	100100	G	A
        1	rs0001_dup1	0	100100	G	A
        1	rs0002	0	100200	A	T
        1	rs0002_dup1	0	100200	A	T
        1	rs0003	0	100300	C	G
        1	rs0001	0	100300	T	G
        1	rs0004	0	100400	T	C
        1	rs0005	0	100500	T	A
        1	rs0002	0	100500	T	G
        ```

        The function creates a `.rmdup.list` file with duplicated rsIDs:

        ```text
        rs0001
        rs0002
        ```

        Then, it returns a list ["rs0001", "rs0002"].

        Args:
            filename (str): Path to the .bim file Accepts path to the file with
                suffix ("toy.bim") or without suffix ("toy")

        Returns:
            list[str]: List of duplicated rsIDs.
        """
        if not filename.endswith(".bim"):
            filename += ".bim"
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} not found.")

        basename = os.path.splitext(filename)[0]
        command = ["bash", self._DUP_RSID, "--bfile", basename]
        subprocess.run(
            command, capture_output=True, text=True, check=False  # noqa: S603
        )

        dup_rsid = pd.read_csv(
            basename + ".rmdup.mismatch", usecols=[0], names=["rsid"]
        )
        return dup_rsid.rsid.to_list()

    @beartype
    def extract_duplicate_individuals(self, filename: str) -> pd.DataFrame:
        """Extract FID and IIDs from duplicated individuals from a .fam file.

        Args:
            filename (str): Path to the .fam file

        Returns:
            pd.DataFrame: DataFrame with the FID and IID of duplicated individuals
        """
        if not os.path.isfile(filename):
            raise FileNotFoundError(f"{filename} was not found,")

        fam_data = pd.read_csv(
            filename,
            sep=r"\s+",
            header=None,
            usecols=[0, 1],
            names=["FID", "IID"],
        )
        dup_data = fam_data[fam_data.duplicated()]
        return dup_data

    @beartype
    def qc(  # noqa: D102
        self,
        filename: str,
        outpath: str | None = None,
        exclude_vars: list[str] | None = None,
        exclude_indvs: pd.DataFrame | None = None,
        maf: float | int | None = None,
        geno_miss: float | int | None = None,
        indv_miss: float | int | None = None,
    ) -> None:

        if maf is not None and not 0 < maf < 0.5:
            raise ValueError(f"maf not in range (0,1). Got {maf}")
        if geno_miss is not None and not 0 <= geno_miss <= 1:
            raise ValueError(f"geno_miss not in range (0,1). Got {geno_miss}")
        if indv_miss is not None and not 0 <= indv_miss <= 1:
            raise ValueError(f"pheno_miss not in range (0,1). Got {indv_miss}")

        for suf in self._SUFFIXES:
            if not os.path.isfile(filename + suf):
                raise FileNotFoundError(f"{filename}{suf} not found")

        outpath = outpath or filename

        if exclude_vars is not None:
            rm_vars_path = outpath + ".rmvars"
            with open(rm_vars_path, "w") as f:
                for var in exclude_vars:
                    f.write(f"{var}\n")

        if exclude_indvs is not None:
            rm_indv_path = outpath + ".rmindv"
            exclude_indvs.to_csv(
                rm_indv_path, sep="\t", header=False, index=False
            )

        command = [
            "bash",
            self._QC_SCRIPT,
            "--bfile",
            filename,
            "--outpath",
            outpath,
        ]
        command = (
            [*command, "--indv-miss", str(indv_miss)] if indv_miss else command
        )
        command = (
            [*command, "--geno-miss", str(geno_miss)] if geno_miss else command
        )
        command = [*command, "--maf", str(maf)] if maf else command
        command = (
            [*command, "--rm-vars", rm_vars_path]
            if exclude_vars is not None
            else command
        )
        command = (
            [*command, "--rm-indv", rm_indv_path]
            if exclude_indvs is not None
            else command
        )

        subprocess.run(
            command, capture_output=True, text=True, check=False  # noqa: S603
        )
