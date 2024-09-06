"""Util functions sub-module used in otterseq."""
import os

import polars as pl
from beartype import beartype

from otterseq.errors import MultiAllelicError


@beartype
def check_multi_allelic(rs_ids: list[str]) -> None:  # noqa: D417
    """Check if there are multi-allelic (duplicated rsIDs) variants.

    Args:
        rs_ids (list[str]): List of rsIDs.

    Raises:
        MultiAllelicError: If there are duplicated rsIDs.
    """
    if len(rs_ids) > len(set(rs_ids)):
        raise MultiAllelicError("Multi-allelic variants found in file.")


@beartype
def read_pca_file(filename: str) -> pl.DataFrame:
    """Read eigenvalues file from PCA output.

    Args:
        filename (str): Path to the filename without suffix (e.g., "data/toy"),
            or with suffix (e.g., "data/toy.eigenvec").

    Raises:
        FileNotFoundError: If the .eigenvec file does not exist.

    Returns:
        pl.DataFrame: Polars DataFrame with .eigenvec content.
    """
    if not filename.endswith(".eigenvec"):
        filename += ".eigenvec"

    if not os.path.isfile(filename):
        raise FileNotFoundError(f"{filename} not found.")

    pca_df = pl.read_csv(filename, has_header=False, separator=" ")
    pca_df.columns = ["fid", "iid"] + [
        f"pc{i}" for i in range(1, len(pca_df.columns) - 1)
    ]

    return pca_df


@beartype
def read_fam_file(filename: str, vars_to_string: bool = True) -> pl.DataFrame:
    """Read .fam file in Polars format.

    Args:
        filename (str): Path to the file without suffix (e.g. "data/toy"),
            or with suffix (e.g. "data/toy.fam").
        vars_to_string (bool, optional): True to convert the `sex` and
            `pheno` variables into string versions. Defaults to True.

    Raises:
        FileNotFoundError: If the .fam file was not found.

    Returns:
        pl.DataFrame: Polars DataFrame with the .fam file content.
    """
    if not filename.endswith(".fam"):
        filename += ".fam"

    if not os.path.isfile(filename):
        raise FileNotFoundError(f"{filename} not found.")

    fam_df = pl.read_csv(
        filename,
        has_header=False,
        separator=" ",
        new_columns=[
            "fid",
            "iid",
            "father_id",
            "mother_id",
            "sex",
            "pheno",
        ],
    )

    if vars_to_string:
        fam_df = fam_df.with_columns(
            pl.when(pl.col("sex") == 2)
            .then(pl.lit("female"))
            .when(pl.col("sex") == 1)
            .then(pl.lit("male"))
            .otherwise(pl.lit("unknown"))
            .alias("sex_str")
        )
        fam_df = fam_df.with_columns(
            pl.when(pl.col("pheno") == 2)
            .then(pl.lit("control"))
            .when(pl.col("pheno") == 1)
            .then(pl.lit("case"))
            .otherwise(pl.lit("unknown"))
            .alias("pheno_str")
        )

    return fam_df
