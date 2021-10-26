"""Test util.py functions."""
import numpy as np
import pytest

import mewarpx
from mewarpx.utils_store import util as mwxutil


def test_utils_check_version(caplog):
    """Note the version passed in here is the "script" version."""
    with pytest.raises(
        ValueError, match="This version of mewarpx .* is older"
    ):
        mwxutil.check_version((np.inf, 0, 0), -1)

    with pytest.raises(
        ValueError, match="This physics version of mewarpx .* is older"
    ):
        mwxutil.check_version((0, 0, 9), np.inf)

    # Check warnings printed
    caplog.clear()
    mwxutil.check_version(
        (-1, 0, 0), mewarpx.__physics_version__)
    assert len(caplog.record_tuples) == 1
    assert "is a newer API version" in caplog.text

    caplog.clear()
    mwxutil.check_version(
        mewarpx.__version_info__, -1
    )
    assert len(caplog.record_tuples) == 1
    assert "is newer than the version" in caplog.text

    # Check an OK setup
    caplog.clear()
    mwxutil.check_version(
        (mewarpx.__version_info__[0],
         mewarpx.__version_info__[1] - 1,
         mewarpx.__version_info__[2] - 1,
         ),
        mewarpx.__physics_version__
    )
    assert len(caplog.record_tuples) == 0
