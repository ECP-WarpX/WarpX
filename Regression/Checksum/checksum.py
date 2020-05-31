import .benchmark
import .utils

class Checksum:

    def __init__(self, test_name, plotfile):
        self.test_name = test_name
        self.plotfile = plotfile
        self.data = utils.read_plotfile(self.plotfile)

    def evaluate():
        benchmark = Benchmark(test_name)
        self.data == benchmark.data
