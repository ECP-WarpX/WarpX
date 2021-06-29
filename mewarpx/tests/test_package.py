"""Test basic package aspects"""
# Native python imports
from builtins import next
import csv
import os

# 3rd-party library imports

# Local imports
import mewarpx
from mewarpx import util


def test_version():
    with open(os.path.join(util.mewarpx_dir, '../changelog.csv'), 'r') as f:
        reader = csv.reader(f)
        next(reader)
        row = next(reader)
        assert row[0] == mewarpx.__version__
