"""
Copyright 2021

This file is part of WarpX

License: BSD-3-Clause-LBNL
"""

import re
import numpy as np
import yt
from yt.frontends.boxlib.data_structures import AMReXDataset

yt.funcs.mylog.setLevel(50)

class Backend:
    ''' Use AMReXDataset as the backend reader to read AMReX plotfiles
    '''

    def __init__(self, filename):
        ''' Constructor: store the dataset object
        '''

        self.dataset = AMReXDataset(filename)

    def fields_list(self):
        ''' Return the list of fields defined on the grid
        '''

        return [item for item in self.dataset.field_list if item[0] == 'boxlib']

    def species_list(self):
        ''' Return the list of species in the dataset
        '''

        return set([item[0] for item in self.dataset.field_list if
                            item[1][:9] == 'particle_' and item[0] != 'all'])

    def n_levels(self):
        ''' Return the number of MR levels in the dataset
        '''

        return self.dataset.max_level+1;

    def get_field_checksum(self, lev, field, test_name):
        ''' Calculate the checksum for a given field at a given level in the dataset
        '''

        lev_grids = [grid for grid in self.dataset.index.grids if grid.Level == lev]
        # Warning: For now, we assume all levels are rectangular
        LeftEdge = np.min(
            np.array([grid.LeftEdge.v for grid in lev_grids]), axis=0)
        all_data_level = self.dataset.covering_grid(
            level=lev, left_edge=LeftEdge, dims=self.dataset.domain_dimensions)
        for field in grid_fields:
            Q = all_data_level[field].v.squeeze()
            data_lev[field[1]] = np.sum(np.abs(Q))
        data['lev=' + str(lev)] = data_lev

    def get_species_attributes(self, species):
        ''' Return the list of attributes for a given species in the dataset
        '''

        return [item[1] for item in self.dataset.field_list if item[0] == species]

    def get_species_checksum(self, species, attribute):
        ''' Calculate the checksum for a given attribute of a given species in the dataset
        '''

        ad = self.dataset.all_data()
        Q = ad[(species, attribute)].v
        return np.sum(np.abs(Q))

