import re
import sys
import json

azure_output_filename = sys.argv[1]

first_char = 29 # first non prefix character, found empirically
test_name_str = "working on test:"
test_name_first_char = 50 # start of test name, found empirically
test_name_last_char = -5 # end of test name, found empirically
plotfile_str = "Plotfile : "
benchmark_path = "benchmarks_json/"
benchmark_suffix = ".json"

current_test = ""

with open(azure_output_filename, "r") as f:
    for line in f:
        # Here we parse lines that read, for example,
        # "working on test: bilinear_filter"
        # and we set current_test = "bilinear_filter"
        if (line.find(test_name_str)!= -1):
            current_test = line[test_name_first_char:test_name_last_char]
        line_noprefix = line[first_char:]
        # Here we parse lines that read, for example,
        # "Plotfile : [lev=0,jx] 1.012345678901234e+16"
        # and we set key1 = "lev=0", key2 = "jx", value = 1.012345678901234e+16
        if (line_noprefix.startswith(plotfile_str)):
            key1 = re.search('\[(.*),', line_noprefix).group(1)
            key2 = re.search(',(.*)\]', line_noprefix).group(1)
            value = float(re.search('\](.*)', line_noprefix).group(1)[1:])
            # Update json file with the new value
            json_filepath = benchmark_path+current_test+benchmark_suffix
            with open(json_filepath, "r") as f_json:
                json_data = json.load(f_json)
            json_data[key1][key2] = value
            with open(json_filepath, "w") as f_json:
                json.dump(json_data, f_json, indent=2)
