#! /usr/bin/env python3

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
    > checksumAPI.evaluate_checksum(test_name, plotfile)

- As a script, to evaluate or to reset a benchmark:
  * Evaluate a benchmark. From a bash terminal:
    $ ./checksumAPI.py --evaluate --plotfile <path/to/plotfile> \
                       --test-name <test name>
  * Reset a benchmark. From a bash terminal:
    $ ./checksumAPI.py --reset-benchmark --plotfile <path/to/plotfile> \
                       --test-name <test name>
'''


def evaluate_checksum(test_name, plotfile, rtol=1.e-9, atol=1.e-40,
                      do_fields=True, do_particles=True):
    '''Compare plotfile checksum with benchmark.

    Read checksum from input plotfile, read benchmark
    corresponding to test_name, and assert their equality.

    @param test_name Name of test, as found between [] in .ini file.
    @param plotfile Plotfile from which the checksum is computed.
    @param rtol Relative tolerance for the comparison.
    @param atol Absolute tolerance for the comparison.
    @param do_fields Whether to compare fields in the checksum.
    @param do_particles Whether to compare particles in the checksum.
    '''
    test_checksum = Checksum(test_name, plotfile, do_fields=do_fields,
                             do_particles=do_particles)
    test_checksum.evaluate(rtol=rtol, atol=atol)


def reset_benchmark(test_name, plotfile, do_fields=True, do_particles=True):
    '''Update the benchmark (overwrites reference json file).

    Overwrite value of benchmark corresponding to
    test_name with checksum read from input plotfile.

    @param test_name Name of test, as found between [] in .ini file.
    @param plotfile Plotfile from which the checksum is computed.
    @param do_fields Whether to write field checksums in the benchmark.
    @param do_particles Whether to write particles checksums in the benchmark.
    '''
    ref_checksum = Checksum(test_name, plotfile, do_fields=do_fields,
                            do_particles=do_particles)
    ref_benchmark = Benchmark(test_name, ref_checksum.data)
    ref_benchmark.reset()


def reset_all_benchmarks(path_to_all_plotfiles):
    '''Update all benchmarks (overwrites reference json files)
    found in path_to_all_plotfiles

    @param path_to_all_plotfiles Path to all plotfiles for which the benchmarks
    are to be reset. The plotfiles should be named <test_name>_plt, which is
    what regression_testing.regtests.py does, provided we're careful enough.
    '''

    # Get list of plotfiles in path_to_all_plotfiles
    plotfile_list = glob.glob(path_to_all_plotfiles + '*_plt?????',
                              recursive=True)
    plotfile_list.sort()

    # Loop over plotfiles and reset the corresponding benchmark
    for plotfile in plotfile_list:
        test_name = os.path.split(plotfile)[1][:-9]
        reset_benchmark(test_name, plotfile)


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
    parser.add_argument('--plotfile', dest='plotfile', type=str, default='',
                        required='--evaluate' in sys.argv or
                        '--reset-benchmark' in sys.argv,
                        help='Name of WarpX plotfile')

    parser.add_argument('--skip-fields', dest='do_fields',
                        default=True, action='store_false',
                        help='If used, do not read/write field checksums')
    parser.add_argument('--skip-particles', dest='do_particles',
                        default=True, action='store_false',
                        help='If used, do not read/write particle checksums')

    # Fields and/or particles are read from plotfile/written to benchmark?
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
    parser.add_argument('--path-to-all-plotfiles',
                        dest='path_to_all_plotfiles', type=str, default='',
                        required='--reset-all-benchmarks' in sys.argv,
                        help='Directory containing all benchmark plotfiles, \
                        typically WarpX-benchmarks generated by \
                        regression_testing/regtest.py')

    args = parser.parse_args()

    if args.reset_benchmark:
        reset_benchmark(args.test_name, args.plotfile,
                        do_fields=args.do_fields,
                        do_particles=args.do_particles)

    if args.evaluate:
        evaluate_checksum(args.test_name, args.plotfile, rtol=args.rtol,
                          atol=args.atol, do_fields=args.do_fields,
                          do_particles=args.do_particles)

    if args.reset_all_benchmarks:
        # WARNING: this mode does not support skip-fields/particles
        # and tolerances
        reset_all_benchmarks(args.path_to_all_plotfiles)
