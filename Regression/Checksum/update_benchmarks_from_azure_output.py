import re
import sys

azure_output_filename = sys.argv[1]

first_char = 29 # first non prefix character
test_name_str = "working on test:"
test_name_first_char = 50 # start of test name
test_name_last_char = -5 # end of test name
plotfile_str = "Plotfile : "
benchmark_path = "benchmarks_json/"
benchmark_suffix = ".json"

current_test = ""

f = open(azure_output_filename, "r")

for line in f:
    if (line.find(test_name_str)!= -1):
        current_test = line[test_name_first_char:test_name_last_char]
    line_noprefix = line[first_char:]
    if (line_noprefix.startswith(plotfile_str)):
        key1 = re.search('\[(.*),', line_noprefix).group(1)
        key2 = re.search(',(.*)\]', line_noprefix).group(1)
        value = re.search('\](.*)', line_noprefix).group(1)[1:]
        f_json = open(benchmark_path+current_test+benchmark_suffix, "r")
        json_data = f_json.readlines()
        f_json.close()
        key1_found = False
        for i, line2 in enumerate(json_data):
            if (line2.find(key1) != -1):
                key1_found = True
            if (key1_found and (line2.find(key2) != -1)):
                prefix = re.search('(.*):', line2).group(1)
                add_comma = True
                if (json_data[i+1].find("}") != -1):
                    add_comma = False
                if (add_comma):
                    json_data[i] = prefix+": "+value+",\n"
                else:
                    json_data[i] = prefix+": "+value+"\n"
                break
        f_json = open(benchmark_path+current_test+benchmark_suffix, "w")
        f_json.writelines(json_data)
        f_json.close()

f.close()

