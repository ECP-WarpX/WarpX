#! /usr/bin/env python3

"""
Copyright 2020

This file is part of WarpX.

License: BSD-3-Clause-LBNL
"""

from checksum import Checksum
from benchmark import Benchmark
import argparse
import glob
import sys
import os

'''
API for WarpX checksum tests. It can be used in two ways:

- Directly use functions below to make a checksum test from a python script.
  Example: add these lines to a WarpX CI Python analysis script to run the
           checksum test.
    > import checksumAPI
    > checksumAPI.evaluate_checksum(test_name, file_name)

- As a script, to evaluate or to reset a benchmark:
  * Evaluate a benchmark. From a bash terminal:
    $ ./checksumAPI.py --evaluate --file_name <path/to/file_name> \
                       --test-name <test name>
  * Reset a benchmark. From a bash terminal:
    $ ./checksumAPI.py --reset-benchmark --file_name <path/to/file_name> \
                       --test-name <test name>
'''


def evaluate_checksum(test_name, file_name, rtol=1.e-9, atol=1.e-40,
                      do_fields=True, do_particles=True):
    '''Compare IO file checksum with benchmark.

    Read checksum from input file_name, read benchmark
    corresponding to test_name, and assert their equality.

    @param test_name Name of test, as found between [] in .ini file.
    @param file_name IO file from which the checksum is computed.
    @param rtol Relative tolerance for the comparison.
    @param atol Absolute tolerance for the comparison.
    @param do_fields Whether to compare fields in the checksum.
    @param do_particles Whether to compare particles in the checksum.
    '''
    test_checksum = Checksum(test_name, file_name, do_fields=do_fields,
                             do_particles=do_particles)
    test_checksum.evaluate(rtol=rtol, atol=atol)


def reset_benchmark(test_name, file_name, do_fields=True, do_particles=True):
    '''Update the benchmark (overwrites reference json file).

    Overwrite value of benchmark corresponding to
    test_name with checksum read from input IO file.

    @param test_name Name of test, as found between [] in .ini file.
    @param file_name IO file from which the checksum is computed.
    @param do_fields Whether to write field checksums in the benchmark.
    @param do_particles Whether to write particles checksums in the benchmark.
    '''
    ref_checksum = Checksum(test_name, file_name, do_fields=do_fields,
                            do_particles=do_particles)
    ref_benchmark = Benchmark(test_name, ref_checksum.data)
    ref_benchmark.reset()


def reset_all_benchmarks(path_to_all_file_names):
    '''Update all benchmarks (overwrites reference json files)
    found in path_to_all_file_names

    @param path_to_all_file_names Path to all IO files for which the benchmarks
    are to be reset. The IO files should be named <test_name>_plt, which is
    what regression_testing.regtests.py does, provided we're careful enough.
    '''

    # Get list of IO files in path_to_all_file_names
    file_name_list = glob.glob(path_to_all_file_names + '*_plt?????',
                              recursive=True)
    file_name_list.sort()

    # Loop over IO files and reset the corresponding benchmark
    for file_name in file_name_list:
        test_name = os.path.split(file_name)[1][:-9]
        reset_benchmark(test_name, file_name)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    # Options relevant to evaluate a checksum or reset a benchmark
    parser.add_argument('--evaluate', dest='evaluate', action='store_true',
                        default=False, help='Evaluate a checksum.')
    parser.add_argument('--reset-benchmark', dest='reset_benchmark',
                        default=False,
                        action='store_true', help='Reset a benchmark.')
    parser.add_argument('--test-name', dest='test_name', type=str, default='',
                        required='--evaluate' in sys.argv or
                        '--reset-benchmark' in sys.argv,
                        help='Name of the test (as in WarpX-tests.ini)')
    parser.add_argument('--file_name', dest='file_name', type=str, default='',
                        required='--evaluate' in sys.argv or
                        '--reset-benchmark' in sys.argv,
                        help='Name of IO file')

    parser.add_argument('--skip-fields', dest='do_fields',
                        default=True, action='store_false',
                        help='If used, do not read/write field checksums')
    parser.add_argument('--skip-particles', dest='do_particles',
                        default=True, action='store_false',
                        help='If used, do not read/write particle checksums')

    # Fields and/or particles are read from IO file/written to benchmark?
    parser.add_argument('--rtol', dest='rtol',
                        type=float, default=1.e-9,
                        help='relative tolerance for comparison')
    parser.add_argument('--atol', dest='atol',
                        type=float, default=1.e-40,
                        help='absolute tolerance for comparison')

    # Option to reset all benchmarks present in a folder.
    parser.add_argument('--reset-all-benchmarks', dest='reset_all_benchmarks',
                        action='store_true', default=False,
                        help='Reset all benchmarks.')
    parser.add_argument('--path-to-all-IO-files',
                        dest='path_to_all_IO-files', type=str, default='',
                        required='--reset-all-benchmarks' in sys.argv,
                        help='Directory containing all benchmark IO files, \
                        typically WarpX-benchmarks generated by \
                        regression_testing/regtest.py')

    args = parser.parse_args()

    if args.reset_benchmark:
        reset_benchmark(args.test_name, args.file_name,
                        do_fields=args.do_fields,
                        do_particles=args.do_particles)

    if args.evaluate:
        evaluate_checksum(args.test_name, args.file_name, rtol=args.rtol,
                          atol=args.atol, do_fields=args.do_fields,
                          do_particles=args.do_particles)

    if args.reset_all_benchmarks:
        # WARNING: this mode does not support skip-fields/particles
        # and tolerances
        reset_all_benchmarks(args.path_to_all_file_names)
