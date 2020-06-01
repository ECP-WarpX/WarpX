from benchmark import Benchmark
import yt
import sys
import numpy as np

yt.funcs.mylog.setLevel(50)

tolerance = 1.e-9


class Checksum:
    '''Class for checksum comparison of one test.
    '''

    def __init__(self, test_name, plotfile):
        '''Constructor

        Store test_name and plotfile name, and compute checksum
        from plotfile and store it in self.data.

        @param self The object pointer.
        @param test_name Name of test, as found between [] in .ini file.
        @param plotfile Plotfile from which the checksum is computed.
        '''

        self.test_name = test_name
        self.plotfile = plotfile
        self.data = self.read_plotfile()

    def read_plotfile(self):
        '''Get checksum from plotfile.

        Read an AMReX plotfile with yt, compute 1 checksum per field and return
        all checksums in a dictionary.
        The checksum of quantity Q is max(abs(Q)).

        @param self The object pointer.
        '''

        ds = yt.load(self.plotfile)
        grid_fields = [item for item in ds.field_list if item[0] == 'boxlib']
        species_list = set([item[0] for item in ds.field_list if
                            item[1][:9] == 'particle_' and item[0] != 'all'])
        data = {}

        # Compute checksum for field quantities
        for lev in range(ds.max_level+1):
            data_lev = {}
            lev_grids = [grid for grid in ds.index.grids if grid.Level == lev]
            # Warning: For now, we assume all levels are rectangular
            LeftEdge = np.min(
                np.array([grid.LeftEdge.v for grid in lev_grids]), axis=0)
            RightEdge = np.max(
                np.array([grid.RightEdge.v for grid in lev_grids]), axis=0)
            all_data_level = ds.covering_grid(
                level=lev, left_edge=LeftEdge, dims=ds.domain_dimensions)
            for field in grid_fields:
                Q = all_data_level[field].v.squeeze()
                data_lev[field[1]] = np.sum(np.abs(Q))
            data['lev=' + str(lev)] = data_lev
        ad = ds.all_data()

        # Compute checksum for particle quantities
        for species in species_list:
            part_fields = [item[1] for item in ds.field_list
                           if item[0] == species]
            data_species = {}
            for field in part_fields:
                Q = ad[(species, field)].v
                data_species[field] = np.sum(np.abs(Q))
            data[species] = data_species

        return data

    def evaluate(self):
        '''Compare plotfile checksum with benchmark.

        Read checksum from input plotfile, read benchmark
        corresponding to test_name, and assert that they are equal.
        Almost all the body of this functions is for
        user-readable print statements.

        @param self The object pointer.
        @param test_name Name of test, as found between [] in .ini file.
        @param plotfile Plotfile from which the checksum is computed.
        '''

        ref_benchmark = Benchmark(self.test_name)

        # Dictionaries have same outer keys (levels, species)?
        if (self.data.keys() != ref_benchmark.data.keys()):
            print("ERROR: Benchmark and plotfile checksum "
                  "have different outer keys:")
            print("Benchmark: %s" % ref_benchmark.data.keys())
            print("Plotfile : %s" % self.data.keys())
            sys.exit(1)

        # Dictionaries have same inner keys (field and particle quantities)?
        for key1 in ref_benchmark.data.keys():
            if (self.data[key1].keys() != ref_benchmark.data[key1].keys()):
                print("ERROR: Benchmark and plotfile checksum have "
                      "different inner keys:")
                print("Common outer keys: %s" % ref_benchmark.data.keys())
                print("Benchmark inner keys in %s: %s"
                      % (key1, ref_benchmark.data[key1].keys()))
                print("Plotfile  inner keys in %s: %s"
                      % (key1, self.data[key1].keys()))
                sys.exit(1)

        # Dictionaries have same values?
        for key1 in ref_benchmark.data.keys():
            for key2 in ref_benchmark.data[key1].keys():
                if (abs(self.data[key1][key2] - ref_benchmark.data[key1][key2])
                   > tolerance*abs(ref_benchmark.data[key1][key2])):
                    print("ERROR: Benchmark and plotfile checksum have "
                          "different value for key [%s,%s]" % (key1, key2))
                    print("Benchmark: [%s,%s] %.15f"
                          % (key1, key2, ref_benchmark.data[key1][key2]))
                    print("Plotfile : [%s,%s] %.15f"
                          % (key1, key2, self.data[key1][key2]))
                    sys.exit(1)
