import json
import re
import sys

azure_output_filename = sys.argv[1]

pattern_test_name = 'working on test: (?P<testname>[\w\-]*)'
pattern_keys_value = 'Plotfile : \[(?P<key1>.*),(?P<key2>.*)\] (?P<value>\w.*)'
benchmark_path = "../../Regression/Checksum/benchmarks_json/"
benchmark_suffix = ".json"

json_file_updated = False
current_test = ""

def dump_json_data(json_filepath, json_data):
    with open(json_filepath, "w") as f_json:
        json.dump(json_data, f_json, indent=2)

with open(azure_output_filename, "r") as f:
    for line in f:
        # Here we search lines that read, for example,
        # "working on test: bilinear_filter"
        # and we set current_test = "bilinear_filter"
        match_test_name = re.search(pattern_test_name, line)
        if match_test_name:
            # If needed, update json file
            if json_file_updated:
                dump_json_data(json_filepath, json_data)
            current_test = match_test_name.group('testname')
            json_file_updated = False

        # Here we search lines that read, for example,
        # "Plotfile : [lev=0,jx] 1.012345678901234e+16"
        # and we update the corresponding json data with the new value
        match_keys_value = re.search(pattern_keys_value, line)
        if match_keys_value:
            # If not already done, read the json file and store it in json_data
            if not json_file_updated:
                json_filepath = benchmark_path+current_test+benchmark_suffix
                with open(json_filepath, "r") as f_json:
                    json_data = json.load(f_json)
            key1 = match_keys_value.group('key1') # "lev=0" in the example
            key2 = match_keys_value.group('key2') # "jx" in the example
            value = float(match_keys_value.group('value')) # "1.012345678901234e+16" in the example
            json_data[key1][key2] = value
            json_file_updated = True

# If needed, update json file for the last test
if json_file_updated:
    dump_json_data(json_filepath, json_data)
