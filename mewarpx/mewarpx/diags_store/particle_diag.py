"""Diagnostic code that wraps the picmi.ParticleDiagonstics class"""

import logging
import time
from mewarpx.diags_store.diag_base import WarpXDiagnostic

import numpy as np

from pywarpx import callbacks, picmi

from mewarpx.mwxrun import mwxrun

logger = logging.getLogger(__name__)

class WarpXParticleDiag(WarpXDiagnostic):
    """Output particle diagnostics every diagnostic period and produce
        plots from the diagnostics

    Contains:
        add_particle_diag (function): Function to add the particle diagnostic to
        the mwxrun.simulation

        post_processing: plot particle diagnostic data, not yet implemented
    """

    def __init__(self, diag_steps, name=None, species=None,
                data_list=None, **kwargs):

        """
        Arguments:
            Initializes the picmi.ParticleDiagnostic and adds the diagnostic to
            the simulation

            diag_steps (int): Number of steps between each diagnostic output

            name (str): name of the diag output folder, defaults to diag

            species (mepicmi.Species): species in the simulation, if None then
            uses all particles in the simulation

            data_list (str list): list of attributes to be outputted by he
            diagnostic, default uses ["position", "momentum", "weighting"]

        """
        self.name = name
        self.species = species
        self.data_list = data_list

        if self.name is None:
            self.name = "diag"
        if self.species is None:
            self.species = mwxrun.simulation.species
        if self.data_list is None:
            self.data_list = ["position", "momentum", "weighting"]

        super(WarpXParticleDiag, self).__init__(diag_steps=diag_steps, **kwargs)
        self.add_particle_diag()

    def add_particle_diag(self):
        particle_diag = picmi.ParticleDiagnostic(
            name=self.name,
            period=self.diag_steps,
            species=self.species,
            data_list=self.data_list
        )

        mwxrun.simulation.add_diagnostic(particle_diag)

    def post_processing():
        raise NotImplementedError
