#! /usr/bin/env python

# Copyright 2021 Modern Electron

# Currently this script just checks that the simulation successfully
# finished. In the future it should be modified to check that the particle
# scraper buffer matches what was directly retrieved during the simulation
# using the Python wrappers.

from pathlib import Path

assert Path('Python_particle_scraper_plt00005').is_dir()

