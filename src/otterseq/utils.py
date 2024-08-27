"""Util functions sub-module used in otterseq."""
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
