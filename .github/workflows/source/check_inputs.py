#! /usr/bin/env python3

import os
import re
import sys

# mandatory prefixes for test names
testname_prefix = ["test_1d_", "test_2d_", "test_3d_", "test_rz_"]

# collect all test names and test input filenames from CMakeLists.txt files
testnames = []
testinputs = []
# walk through all files under Examples/, including subdirectories
for dirpath, dirnames, filenames in os.walk(top="./Examples"):
    # loop over CMakeLists.txt files
    for name in [filename for filename in filenames if filename == "CMakeLists.txt"]:
        filepath = os.path.join(dirpath, name)
        # open CMakeLists.txt file
        with open(filepath) as f:
            # loop over lines of CMakeLists.txt
            for line in f:
                # strip leading whitespaces
                line = line.lstrip()
                # find lines where 'add_warpx_test' is called
                if re.match("add_warpx_test", line):
                    # strip leading whitespaces, remove end-of-line comments
                    testname = next(f).lstrip().split(" ")[0]
                    testnames.append(testname)
                    # skip lines related to other function arguments
                    for _ in range(4):
                        next(f)
                    # strip leading whitespaces, remove end-of-line comments
                    testinput = next(f).lstrip().split(" ")[0]
                    # some Python input scripts are quoted
                    # to account for command-line arguments:
                    # strip initial quotation mark from string
                    if testinput.startswith('"'):
                        testinput = re.sub('"', "", testinput)
                    # extract filename from path
                    testinput = os.path.split(testinput)[1]
                    testinputs.append(testinput)

# check consistency of test names and test input filenames
print("\nCheck that test names and input names are correct...")
wrong_testname = False
wrong_testinput = False
for testname, testinput in zip(testnames, testinputs):
    if not testname.startswith(tuple(testname_prefix)):
        print(f"Wrong test  name: {testname}")
        wrong_testname = True
    if "restart" in testname:
        if not testname.endswith("_restart") and not testname.endswith("_restart.py"):
            print(f"Wrong test  name: {testname}")
            wrong_testname = True
    if not testinput == f"inputs_{testname}" and not testinput.endswith("_picmi.py"):
        print(f"Wrong input name: {testinput}")
        wrong_testinput = True

if wrong_testname:
    print(f"NOTE: Test names must start with one of {testname_prefix}.")
    print("      Test names must end with '_restart' for restart tests")
    print("      (with the extension '.py' in the case of PICMI input scripts).")
if wrong_testinput:
    print("NOTE: Test input names must start with 'inputs_' followed by the test name")
    print("      (with the extension '.py' in the case of PICMI input scripts).")

# check that all input files in Examples/ are tested
print("\nCheck that all test input files are tested...")
missing_input = False
# walk through all files under Examples/, including subdirectories
for dirpath, dirnames, filenames in os.walk(top="./Examples"):
    # loop over files starting with "inputs_test_"
    for name in [
        filename for filename in filenames if filename.startswith("inputs_test_")
    ]:
        if name not in testinputs:
            print(f"Input not tested: {os.path.join(dirpath, name)}")
            missing_input = True

if missing_input:
    print("NOTE: All test input files must be tested.\n")
else:
    print()

if wrong_testname or wrong_testinput or missing_input:
    sys.exit("FAILED\n")
