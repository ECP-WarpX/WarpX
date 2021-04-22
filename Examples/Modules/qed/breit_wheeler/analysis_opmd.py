#! /usr/bin/env python
# Copyright 2019 Luca Fedeli
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# -*- coding: utf-8 -*-


import sys
import openpmd_api as io
#sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
#import checksumAPI

import analysis_core as ac

# This script is a frontend for the analysis routines
# in analysis_core.py (please refer to this file for
# a full description). It reads output files in openPMD
# format and extracts the data needed for
# the analysis routines.

def main():
    print("Opening openPMD output")
    prefix = sys.argv[1]
    series = io.Series(prefix+"/openpmd_%T.h5", io.Access.read_only)
    data_set_end = series.iterations[1]

    # get simulation time
    sim_time = data_set_end.time

    # get particle data
    particle_data = {}

    names, types = ac.get_all_species_names_and_types()
    for spec_name_type in zip(names, types):
        spec_name = spec_name_type[0]
        is_photon = spec_name_type[1]
        data = {}

        px = data_set_end.particles[spec_name]["momentum"]["x"][:]
        py = data_set_end.particles[spec_name]["momentum"]["y"][:]
        pz = data_set_end.particles[spec_name]["momentum"]["z"][:]
        w = data_set_end.particles[spec_name]["weighting"][io.Mesh_Record_Component.SCALAR][:]

        if is_photon :
            opt = data_set_end.particles[spec_name]["opticalDepthBW"][io.Mesh_Record_Component.SCALAR][:]
        else:
            opt = data_set_end.particles[spec_name]["opticalDepthQSR"][io.Mesh_Record_Component.SCALAR][:]

        series.flush()

        data["px"] = px
        data["py"] = py
        data["pz"] = pz
        data["w"] = w
        data["opt"] = opt

        particle_data[spec_name] = data

    ac.check(sim_time, particle_data)

    #test_name = filename_end[:-9] # Could also be os.path.split(os.getcwd())[1]
    #checksumAPI.evaluate_checksum(test_name, filename_end)

if __name__ == "__main__":
    main()
