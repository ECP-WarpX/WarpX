import .checksum
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('evaluate', dest='evaluate', action='store_true', type=bool, default=False,
                    help='Evaluate a checksum. Necessitates arguments --test-name and --plotfile')
parser.add_argument('reset-benchmark', dest='reset_benchmark', action='store_true', type=bool, default=False,
                    help='Reset a benchmark. Necessitates arguments --test-name and --plotfile')
parser.add_argument('--test-name', dest='test_name', type=str, default='', required=True,
                    help='Name of the test (between [] in WarpX-tests.ini)')
parser.add_argument('--plotfile', dest='plotfile', type=str, default='', required=True,
                    help='Name of WarpX plotfile')
args = parser.parse_args()

def evaluate_checksum(test_name, plotfile):
    test_checksum = Checksum(test_name, plotfile)
    test_checksum.evaluate()

def reset_benchmark(test_name, plotfile):
    ref_checksum = Checksum(test_name, plotfile)
    benchmark = Benchmark(test_name, ref_checksum.data)
    benchmark.reset()

if __name__ == '__main__':

    if args.reset_benchmark:
        reset_benchmark(args.test_name, args.plotfile)

    if args.evaluate:
        evaluate_checksum(args.test_name, args.plotfile)
