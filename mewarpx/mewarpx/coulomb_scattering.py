"""
Module for Coulomb scattering in WarpX.
"""
import logging

import numpy as np
from pywarpx import callbacks

from mewarpx.mwxrun import mwxrun
import mewarpx.utils_store.mwxconstants as constants
import mewarpx.utils_store.util as mwxutil

# Get module-level logger
logger = logging.getLogger(__name__)


class LangevinElectronIonScattering(object):

    """
    Author: R. E. Groenewald

    Created: 12/30/2021

    Class to specifically handle electron-ion scattering using a simplified
    Langevin collision operator as derived by Manheimer et al. (1997)
    https://doi.org/10.1006/jcph.1997.5834. We assume infinitely massive ions
    with a velocity distribution function f(v) = delta(v) where delta is the
    kronecker delta. This allows massive simplification of the calculation of
    the diffusion coefficients. The friction coefficient is not calculated but
    instead energy conservation is used to determine scattering in the parallel
    direction (of a particle's velocity vector) given the diffusion scattering
    in the transverse direction."""

    def __init__(self, electron_species, ion_species, log_lambda=None,
                 subcycling_steps=1):
        """
        Arguments:
            electron_species (:class:`mewarpx.mespecies.Species`): Electron
                species that will be scattered off the ion species.
            ion_species (:class:`mewarpx.mespecies.Species`): Ion species from
                which electrons will be scattered.
            log_lambda (float): If specified, a fixed value for the Coulomb
                logarithm. If not specified it will be calculated according to
                the NRL formulary.
            subcycling_steps (int): Number of steps between updating the grid
                quantities. Default 1.
        """
        self.collider = electron_species
        self.field = ion_species

        self.log_lambda = log_lambda
        self.subcycling_steps = subcycling_steps

        self.nu_coef = (
            self.collider.sq**2 * self.field.sq**2
            / (4.0 * np.pi * constants.epsilon_0**2 * self.collider.sm**2)
        )

        if mwxrun.geom_str not in ['Z', 'XZ']:
            raise NotImplementedError(
                "Currently LangevinElectronIonScattering is only implemented "
                "for Z and XZ geometries."
            )

        print_str = (
            "Initialized electron-ion Coulomb scattering for species "
            f"{self.collider.name} and {self.field.name}."
        )
        if self.log_lambda is None:
            print_str += 'The Coulomb logarithm will be calculated.'
        else:
            print_str += (
                'The Coulomb logarithm is fixed at %.2f.' % self.log_lambda
            )
        logger.info(print_str)

        callbacks.installafterstep(self.run_scattering_method)

    def get_grid_quantities(self):
        """Function to return ion and electron density at the electron
        positions. Will also return the electron temperature if needed to
        calculate the Coulomb logarithm."""

        self.ion_density_grid = mwxrun.get_gathered_rho_grid(self.field.name)
        if mwxrun.geom_str == 'Z':
            self.ion_density_grid = self.ion_density_grid[:,0] / self.field.sq
        elif mwxrun.geom_str == 'XZ':
            self.ion_density_grid = self.ion_density_grid[:,:,0] / self.field.sq

        if self.log_lambda is None:
            raise NotImplementedError(
                "Calculation of the Coulomb logarithm is not yet supported."
            )
            # Instantiate a particle processor for the electron species used
            # to get the electron density and temperature on grid if the
            # Coulomb logarithm needs to be calculated. Note that the
            # temperature is given in eV.
            '''
            mypp = particles.WarpParticleProcessor(js=self.collider.jslist[0])
            mypp.collect_velocities()
            self.elec_temp_grid = (mypp.Tx + mypp.Ty + mypp.Tz) / 3.0
            mypp.collect_density()
            self.elec_density_grid = mypp.density
            '''

    def get_coulomb_log(self, coords):
        """Function to return the Coulomb logarithm for electron-ion
        collisions in the case appropriate to thermionic converters, see p 34
        of 2019 NRL Plasma Formulary."""
        if self.log_lambda is not None:
            return self.log_lambda

        elec_density = mwxutil.interpolate_from_grid(
            coords, self.elec_density_grid)
        elec_temperature = mwxutil.interpolate_from_grid(
            coords, self.elec_temp_grid)

        log_lambda = (
            23.0 - np.log(
                np.sqrt(elec_density * 1e-6 / elec_temperature**3)
                * self.field.sq / minutil.e
            )
        )
        log_lambda[np.logical_not(log_lambda > 0)] = 0.0
        return log_lambda

    def run_scattering_method(self):
        """Function to execute the Langevin Coulomb collision operator
        scattering method."""
        if (mwxrun.get_it() - 1) % self.subcycling_steps == 0:
            self.get_grid_quantities()

        # Short-circuit if no particles present. This is both more efficient
        # and avoids crashes empty particle lists can cause.
        if mwxrun.sim_ext.get_particle_count(self.collider.name) == 0:
            return

        # collect electron particle positions (structs-of-arrays)
        structs = mwxrun.sim_ext.get_particle_structs(self.collider.name, 0)

        # collect electron particle velocities (array-of-structs)
        ux_arrays = mwxrun.sim_ext.get_particle_ux(self.collider.name)
        uy_arrays = mwxrun.sim_ext.get_particle_uy(self.collider.name)
        uz_arrays = mwxrun.sim_ext.get_particle_uz(self.collider.name)

        # loop over tiles and scatter the electrons appropriately
        for ii in range(len(structs)):

            ux = ux_arrays[ii]
            uy = uy_arrays[ii]
            uz = uz_arrays[ii]

            # create a new array of the velocity components for convenience
            v = np.array([ux, uy, uz])
            v_mag = np.sqrt(np.sum(v**2, axis=0))
            v_perp = np.sqrt(v[0]**2 + v[1]**2)

            # interpolate ion density to electron positions
            coords = np.zeros((mwxrun.dim, len(ux)))
            if mwxrun.geom_str == 'Z':
                coords[0] = structs[ii]['x']
            elif mwxrun.geom_str == 'XZ':
                coords[0] = structs[ii]['x']
                coords[1] = structs[ii]['y']

            density = mwxutil.interpolate_from_grid(
                coords, self.ion_density_grid
            )

            # calculate diffusion coefficient assuming infinitely massive ions
            coulomb_log = self.get_coulomb_log(coords)
            d11 = self.nu_coef * density * coulomb_log / v_mag

            # generate diffusion scattering vectors in the perpendicular plane
            sigma = np.sqrt(mwxrun.get_dt()*d11)
            Q1 = np.random.normal(0, sigma, len(v_mag))
            Q2 = np.random.normal(0, sigma, len(v_mag))

            # calculate rotation angles to parallel coordinates frame
            cos_theta = v[2] / v_mag
            sin_theta = v_perp / v_mag
            cos_phi = v[0] / v_perp
            sin_phi = v[1] / v_perp

            # enforce energy conservation
            dif = v_mag**2 - Q1**2 - Q2**2
            # pick out unphysical points - for these we use isotropic scattering
            idx = np.where(dif <= 0)[0]
            isotropized_vels = mwxutil.get_vel_vector(v_mag[idx])
            Q1[idx] = isotropized_vels[:, 0]
            Q2[idx] = isotropized_vels[:, 1]
            dif[idx] = isotropized_vels[:, 2]**2
            Q3 = np.sqrt(dif) - v_mag

            # add the dynamical friction component
            # F = self.nu_coef * density * coulomb_log / v_mag**2
            # Q3 -= F * warp.top.dt

            # transform Q from the parallel coordinates to the lab frame
            Q = np.array([
                Q1 * cos_theta * cos_phi - Q2 * sin_phi + Q3 * sin_theta * cos_phi,
                Q1 * cos_theta * sin_phi + Q2 * cos_phi + Q3 * sin_theta * sin_phi,
                -Q1 * sin_theta + Q3 * cos_theta
            ])

            # compare initial and final kinetic energy for sanity checking
            # E_init = 0.5 * constants.m_e / constants.eV_SI * (
            #    ux**2 + uy**2 + uz**2
            # )

            ux[:] += Q[0]
            uy[:] += Q[1]
            uz[:] += Q[2]

            # E_final = 0.5 * constants.m_e / constants.eV_SI * (
            #    ux**2 + uy**2 + uz**2
            # )
            # print(f"delta energy = {np.sum(E_init - E_final):.3e} eV")
            # assert np.allclose(E_init, E_final)
