# Copyright 2023 Neil Zaim
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import json
import re
import sys

"""
This Python script updates the Azure benchmarks automatically using a raw Azure output textfile
that is given as the first and only argument of the script.

In the Azure output, we read the lines contained between
"New file for Test_Name:"
and the next occurrence of
"'----------------'"
And use these lines to update the benchmarks
"""

azure_output_filename = sys.argv[1]

pattern_test_name = "New file for (?P<testname>[\w\-]*)"
closing_string = "----------------"
benchmark_path = "../../Regression/Checksum/benchmarks_json/"
benchmark_suffix = ".json"

first_line_read = False
current_test = ""

with open(azure_output_filename, "r") as f:
    for line in f:
        if current_test == "":
            # Here we search lines that read, for example,
            # "New file for LaserAcceleration_BTD"
            # and we set current_test = "LaserAcceleration_BTD"
            match_test_name = re.search(pattern_test_name, line)
            if match_test_name:
                current_test = match_test_name.group("testname")
                new_file_string = ""

        else:
            # We add each line to the new file string until we find the line containing
            # "----------------"
            # which indicates that we have read the new file entirely

            if closing_string not in line:
                if not first_line_read:
                    # Raw Azure output comes with a prefix at the beginning of each line that we do
                    # not need here. The first line that we will read is the prefix followed by the
                    # "{" character, so we determine how long the prefix is by finding the last
                    # occurrence of the "{" character in this line.
                    azure_indent = line.rfind("{")
                    first_line_read = True
                new_file_string += line[azure_indent:]

            else:
                # We have read the new file entirely. Dump it in the json file.
                new_file_json = json.loads(new_file_string)
                json_filepath = benchmark_path + current_test + benchmark_suffix
                with open(json_filepath, "w") as f_json:
                    json.dump(new_file_json, f_json, indent=2)
                current_test = ""
