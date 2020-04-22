# Copyright 2017-2020 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from .Bucket import Bucket

diagnostics = Bucket('diagnostics', ndiags=0, diags_names=[])
diagnostics_list = []

