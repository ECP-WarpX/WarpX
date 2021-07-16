#! /usr/bin/env python

# Copyright 2021 Modern Electron

# This script simply checks that the PICMI_inputs_2d.py run output
# diagnostics, which confirms that the PICMI MCC interface works otherwise
# the run would've crashed.

import glob

files = sorted(glob.glob('Python_background_mcc_plt*'))[1:]
assert len(files) > 0
