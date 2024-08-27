"""Unit test of `otterseq.utils`."""
import pytest

from otterseq.errors import MultiAllelicError
from otterseq.utils import check_multi_allelic


def test_check_multi_allelic(  # noqa: D103
    multiallelic_rsids: list[str],
) -> None:
    with pytest.raises(MultiAllelicError):
        check_multi_allelic(multiallelic_rsids)
