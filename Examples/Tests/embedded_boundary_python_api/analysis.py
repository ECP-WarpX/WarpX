#!/usr/bin/env python3

# This script just checks that the PICMI file executed successfully.
# If it did there will be a plotfile for the final step.

import sys

step = int(sys.argv[1][-5:])

assert step == 2
