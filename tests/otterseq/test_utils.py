"""Unit test of `otterseq.utils`."""

import os
import tempfile

import polars as pl
import pytest

from otterseq.errors import MultiAllelicError
from otterseq.utils import check_multi_allelic, read_fam_file, read_pca_file


def test_check_multi_allelic(  # noqa: D103
    multiallelic_rsids: list[str],
) -> None:
    """Test that multi-allelic variants raise the correct error."""
    with pytest.raises(MultiAllelicError):
        check_multi_allelic(multiallelic_rsids)


def test_check_multi_allelic_valid() -> None:
    """Test that valid (non-duplicated) rsIDs don't raise an error."""
    valid_rsids = ["rs0001", "rs0002", "rs0003"]
    # Should not raise any exception
    check_multi_allelic(valid_rsids)


def test_read_pca_file_success() -> None:
    """Test successful reading of PCA eigenvec file."""
    # Create a temporary PCA file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".eigenvec", delete=False
    ) as f:
        f.write("FAM001 1 0.1 0.2 0.3\n")
        f.write("FAM001 2 0.4 0.5 0.6\n")
        f.write("FAM001 3 0.7 0.8 0.9\n")
        temp_file = f.name

    try:
        df = read_pca_file(temp_file)

        # Check structure
        assert isinstance(df, pl.DataFrame)
        assert list(df.columns) == ["fid", "iid", "pc1", "pc2", "pc3"]
        assert df.shape == (3, 5)

        # Check data
        assert df["fid"].to_list() == ["FAM001", "FAM001", "FAM001"]
        assert df["iid"].to_list() == ["1", "2", "3"]
        assert df["pc1"].to_list() == [0.1, 0.4, 0.7]

    finally:
        os.unlink(temp_file)


def test_read_pca_file_without_suffix() -> None:
    """Test reading PCA file without .eigenvec suffix."""
    # Create a temporary PCA file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".eigenvec", delete=False
    ) as f:
        f.write("FAM001 1 0.1 0.2\n")
        f.write("FAM001 2 0.3 0.4\n")
        temp_file = f.name

    try:
        # Test without suffix
        base_name = temp_file.replace(".eigenvec", "")
        df = read_pca_file(base_name)

        assert isinstance(df, pl.DataFrame)
        assert list(df.columns) == ["fid", "iid", "pc1", "pc2"]
        assert df.shape == (2, 4)

    finally:
        os.unlink(temp_file)


def test_read_pca_file_not_found() -> None:
    """Test that FileNotFoundError is raised for non-existent file."""
    with pytest.raises(FileNotFoundError):
        read_pca_file("non_existent_file")


def test_read_fam_file_success() -> None:
    """Test successful reading of .fam file."""
    # Create a temporary .fam file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fam", delete=False
    ) as f:
        f.write("FAM001 1 0 0 1 2\n")
        f.write("FAM001 2 0 0 2 1\n")
        f.write("FAM001 3 0 0 1 1\n")
        f.write("FAM001 4 0 0 0 2\n")
        temp_file = f.name

    try:
        df = read_fam_file(temp_file)

        # Check structure
        assert isinstance(df, pl.DataFrame)
        expected_columns = [
            "fid",
            "iid",
            "father_id",
            "mother_id",
            "sex",
            "pheno",
            "sex_str",
            "pheno_str",
        ]
        assert list(df.columns) == expected_columns
        assert df.shape == (4, 8)

        # Check data
        assert df["fid"].to_list() == ["FAM001", "FAM001", "FAM001", "FAM001"]
        assert df["iid"].to_list() == ["1", "2", "3", "4"]
        assert df["sex_str"].to_list() == ["male", "female", "male", "unknown"]
        assert df["pheno_str"].to_list() == [
            "control",
            "case",
            "case",
            "control",
        ]

    finally:
        os.unlink(temp_file)


def test_read_fam_file_without_suffix() -> None:
    """Test reading .fam file without suffix."""
    # Create a temporary .fam file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fam", delete=False
    ) as f:
        f.write("FAM001 1 0 0 1 2\n")
        f.write("FAM001 2 0 0 2 1\n")
        temp_file = f.name

    try:
        # Test without suffix
        base_name = temp_file.replace(".fam", "")
        df = read_fam_file(base_name)

        assert isinstance(df, pl.DataFrame)
        assert df.shape == (2, 8)  # 6 original + 2 string columns

    finally:
        os.unlink(temp_file)


def test_read_fam_file_vars_to_string_false() -> None:
    """Test reading .fam file without string conversion."""
    # Create a temporary .fam file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fam", delete=False
    ) as f:
        f.write("FAM001 1 0 0 1 2\n")
        f.write("FAM001 2 0 0 2 1\n")
        temp_file = f.name

    try:
        df = read_fam_file(temp_file, vars_to_string=False)

        # Check structure - should not have string columns
        expected_columns = [
            "fid",
            "iid",
            "father_id",
            "mother_id",
            "sex",
            "pheno",
        ]
        assert list(df.columns) == expected_columns
        assert df.shape == (2, 6)

        # Check that sex and pheno are still numeric
        assert df["sex"].to_list() == [1, 2]
        assert df["pheno"].to_list() == [2, 1]

    finally:
        os.unlink(temp_file)


def test_read_fam_file_not_found() -> None:
    """Test that FileNotFoundError is raised for non-existent .fam file."""
    with pytest.raises(FileNotFoundError):
        read_fam_file("non_existent_file")


def test_read_fam_file_edge_cases() -> None:
    """Test .fam file with edge cases (unknown sex/pheno values)."""
    # Create a temporary .fam file with edge cases
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fam", delete=False
    ) as f:
        f.write("FAM001 1 0 0 1 2\n")  # male, control
        f.write("FAM001 2 0 0 2 1\n")  # female, case
        f.write("FAM001 3 0 0 0 0\n")  # unknown sex, unknown pheno
        f.write("FAM001 4 0 0 3 3\n")  # invalid sex, invalid pheno
        temp_file = f.name

    try:
        df = read_fam_file(temp_file)

        # Check string mappings
        assert df["sex_str"].to_list() == [
            "male",
            "female",
            "unknown",
            "unknown",
        ]
        assert df["pheno_str"].to_list() == [
            "control",
            "case",
            "unknown",
            "unknown",
        ]

    finally:
        os.unlink(temp_file)
