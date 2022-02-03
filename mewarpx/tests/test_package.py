"""Test basic package aspects"""
# Native python imports
from builtins import next
import csv
import os

# Local imports
import mewarpx
from mewarpx.utils_store import util

# 3rd-party library imports


def test_version():
    with open(os.path.join(util.mewarpx_dir, '../changelog.csv'), 'r') as f:
        reader = csv.reader(f)
        next(reader)
        row = next(reader)
        assert row[0] == mewarpx.__version__
        assert row[1].strip() == str(mewarpx.__physics_version__)
