#! /usr/bin/env python

# Copyright 2019 Luca Fedeli, Maxence Thevenet
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# -*- coding: utf-8 -*-


import sys
import yt

import analysis_core as ac

# This script is a frontend for the analysis routines
# in analysis_core.py (please refer to this file for
# a full description). It reads output files in yt
# format and convert and extracts the data needed for
# the analysis routines. 

def main():
    print("Opening yt ")
    filename_end = sys.argv[1]
    data_set_end = yt.load(filename_end)

    # get simulation time
    sim_time = data_set_end.current_time.to_value()
    
    # get particle data
    all_data_end = data_set_end.all_data()
    particle_data = {}

    for spec_name_type in zip(ac.get_all_species_names_and_types()):
        spec_name = spec_name_type[0]
        is_photon = spec_name_type[1]
        data = {}
        data["px"] = all_data_end[spec_name,"particle_momentum_x"].v
        data["py"] = all_data_end[spec_name,"particle_momentum_y"].v
        data["pz"] = all_data_end[spec_name,"particle_momentum_z"].v
        data["w"] = all_data_end[spec_name,"particle_weighting"].v

        if is_photon :
            try:
                data["opt"] =  all_data_end[spec_name,"particle_optical_depth_BW"].v
            except:
                pass
        else:
            try:
                data["opt"] = all_data_end[spec_name,"particle_optical_depth_QSR"].v
            except:
                pass
        particle_data[spec_name] = data
    
    ac.check(sim_time, particle_data)

if __name__ == "__main__":
    main()
