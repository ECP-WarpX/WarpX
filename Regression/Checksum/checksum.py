"""
 Copyright 2020

 This file is part of WarpX.

 License: BSD-3-Clause-LBNL
 """

from benchmark import Benchmark
import yt
import re
import sys
import numpy as np

yt.funcs.mylog.setLevel(50)


class Checksum:
    '''Class for checksum comparison of one test.
    '''

    def __init__(self, test_name, plotfile, do_fields=True, do_particles=True):
        '''Constructor

        Store test_name and plotfile name, and compute checksum
        from plotfile and store it in self.data.

        @param self The object pointer.
        @param test_name Name of test, as found between [] in .ini file.
        @param plotfile Plotfile from which the checksum is computed.
        @param do_fields Whether to compare fields in the checksum.
        @param do_particles Whether to compare particles in the checksum.
        '''

        self.test_name = test_name
        self.plotfile = plotfile
        self.data = self.read_plotfile(do_fields=do_fields,
                                       do_particles=do_particles)

    def read_plotfile(self, do_fields=True, do_particles=True):
        '''Get checksum from plotfile.

        Read an AMReX plotfile with yt, compute 1 checksum per field and return
        all checksums in a dictionary.
        The checksum of quantity Q is max(abs(Q)).

        @param self The object pointer.
        @param do_fields Whether to read fields from the plotfile.
        @param do_particles Whether to read particles from the plotfile.
        '''

        ds = yt.load(self.plotfile)
        grid_fields = [item for item in ds.field_list if item[0] == 'boxlib']
        species_list = set([item[0] for item in ds.field_list if
                            item[1][:9] == 'particle_' and item[0] != 'all'])

        data = {}

        # Compute checksum for field quantities
        if do_fields:
            # Boolean variable parsed from test name needed for workaround below
            galilean = True if re.search('galilean', self.test_name) else False
            for lev in range(ds.max_level+1):
                data_lev = {}
                lev_grids = [grid for grid in ds.index.grids
                             if grid.Level == lev]
                if not (galilean):
                    # Warning: For now, we assume all levels are rectangular
                    LeftEdge = np.min(
                        np.array([grid.LeftEdge.v for grid in lev_grids]), axis=0)
                    all_data_level = ds.covering_grid(
                        level=lev, left_edge=LeftEdge, dims=ds.domain_dimensions)
                    for field in grid_fields:
                        Q = all_data_level[field].v.squeeze()
                        data_lev[field[1]] = np.sum(np.abs(Q))
                # Workaround for Galilean tests: the standard procedure above
                # does not seem to read 2D fields data correctly
                elif (galilean):
                    for field in grid_fields:
                        Q = ds.index.grids[lev][field].v.squeeze()
                        data_lev[field[1]] = np.sum(np.abs(Q))
                data['lev=' + str(lev)] = data_lev

        # Compute checksum for particle quantities
        if do_particles:
            ad = ds.all_data()
            for species in species_list:
                part_fields = [item[1] for item in ds.field_list
                               if item[0] == species]
                data_species = {}
                for field in part_fields:
                    Q = ad[(species, field)].v
                    data_species[field] = np.sum(np.abs(Q))
                data[species] = data_species

        return data

    def evaluate(self, rtol=1.e-9, atol=1.e-40):
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
        checksums_differ = False
        for key1 in ref_benchmark.data.keys():
            for key2 in ref_benchmark.data[key1].keys():
                passed = np.isclose(self.data[key1][key2],
                                    ref_benchmark.data[key1][key2],
                                    rtol=rtol, atol=atol)
                if not passed:
                    print("ERROR: Benchmark and plotfile checksum have "
                          "different value for key [%s,%s]" % (key1, key2))
                    print("Benchmark: [%s,%s] %.40f"
                          % (key1, key2, ref_benchmark.data[key1][key2]))
                    print("Plotfile : [%s,%s] %.40f"
                          % (key1, key2, self.data[key1][key2]))
                    checksums_differ = True
        if checksums_differ:
            sys.exit(1)
