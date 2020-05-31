import yt
import numpy as np

def read_plotfile(plotfile):
    '''
    Read an AMReX plotfile with yt, compute 1 checksum per field and return all checksums in a
    dictionary. The checksum of quantity Q is max(abs(Q)).
    '''

    ds = yt.load( plotfile )
    grid_fields = [item for item in ds.field_list if item[0] == 'boxlib']
    part_fields  = [item for item in ds.field_list if item[1][:9] == 'particle_' and item[0] != 'all']
    data = {}

    # Compute checksum for field quantities
    for lev in range(ds.max_level+1):
        data_lev = {}
        lev_grids = [grid for grid in ds.index.grids if grid.Level == lev]
        # Warning: For now, we assume all levels are rectangular
        LeftEdge = np.min(np.array([grid.LeftEdge.v for grid in lev_grids]),axis=0)
        RightEdge = np.max(np.array([grid.RightEdge.v for grid in lev_grids]),axis=0)
        all_data_level = ds.covering_grid(level=lev,left_edge=LeftEdge, dims=ds.domain_dimensions)
        for field in grid_fields:
            Q = all_data_level[field].v.squeeze()
            data_lev[field[1]] = np.sum(np.abs(Q))
        data['lev=' + str(lev)] = data_lev
    ad = ds.all_data()

    # Compute checksum for particle quantities
    for species in [field[0] for field in part_fields]:
        data_species = {}
        for field in [field[1] for field in part_fields]:
            Q = ad[(species, field)].v
            data_species[field] = np.sum(np.abs(Q))
        data[species] = data_species

    return data
