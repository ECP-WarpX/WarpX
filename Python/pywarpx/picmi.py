# Copyright 2018-2022 Andrew Myers, David Grote, Ligia Diana Amorim
# Maxence Thevenet, Remi Lehe, Revathi Jambunathan, Lorenzo Giacomel
#
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Classes following the PICMI standard
"""
import os
import re

import numpy as np
import periodictable

import picmistandard
import pywarpx

codename = 'warpx'
picmistandard.register_codename(codename)

# dictionary to map field boundary conditions from picmistandard to WarpX
BC_map = {
    'open':'pml', 'dirichlet':'pec', 'periodic':'periodic', 'damped':'damped',
    'absorbing_silver_mueller':'absorbing_silver_mueller',
    'neumann':'neumann', 'none':'none', None:'none'
}

class constants:
    # --- Put the constants in their own namespace
    # --- Values from WarpXConst.H
    c = 299792458.
    ep0 = 8.8541878128e-12
    mu0 = 1.25663706212e-06
    q_e = 1.602176634e-19
    m_e = 9.1093837015e-31
    m_p = 1.67262192369e-27
    hbar = 1.054571817e-34
    kb = 1.380649e-23

picmistandard.register_constants(constants)


class Species(picmistandard.PICMI_Species):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ----------
    warpx_boost_adjust_transverse_positions: bool, default=False
        Whether to adjust transverse positions when apply the boost
        to the simulation frame

    warpx_self_fields_required_precision: float, default=1.e-11
        Relative precision on the electrostatic solver
        (when using the relativistic solver)

    warpx_self_fields_absolute_tolerance: float, default=0.
        Absolute precision on the electrostatic solver
        (when using the relativistic solver)

    warpx_self_fields_max_iters: integer, default=200
        Maximum number of iterations for the electrostatic
        solver for the species

    warpx_self_fields_verbosity: integer, default=2
        Level of verbosity for the electrostatic solver

    warpx_save_previous_position: bool, default=False
        Whether to save the old particle positions

    warpx_do_not_deposit: bool, default=False
        Whether or not to deposit the charge and current density for
        for this species

    warpx_do_not_push: bool, default=False
        Whether or not to push this species

    warpx_do_not_gather: bool, default=False
        Whether or not to gahter the fields from grids for this species

    warpx_random_theta: bool, default=True
        Whether or not to add random angle to the particles in theta
        when in RZ mode.

    warpx_reflection_model_xlo: string, default='0.'
        Expression (in terms of the velocity "v") specifying the probability
        that the particle will reflect on the lower x boundary

    warpx_reflection_model_xhi: string, default='0.'
        Expression (in terms of the velocity "v") specifying the probability
        that the particle will reflect on the upper x boundary

    warpx_reflection_model_ylo: string, default='0.'
        Expression (in terms of the velocity "v") specifying the probability
        that the particle will reflect on the lower y boundary

    warpx_reflection_model_yhi: string, default='0.'
        Expression (in terms of the velocity "v") specifying the probability
        that the particle will reflect on the upper y boundary

    warpx_reflection_model_zlo: string, default='0.'
        Expression (in terms of the velocity "v") specifying the probability
        that the particle will reflect on the lower z boundary

    warpx_reflection_model_zhi: string, default='0.'
        Expression (in terms of the velocity "v") specifying the probability
        that the particle will reflect on the upper z boundary

    warpx_save_particles_at_xlo: bool, default=False
        Whether to save particles lost at the lower x boundary

    warpx_save_particles_at_xhi: bool, default=False
        Whether to save particles lost at the upper x boundary

    warpx_save_particles_at_ylo: bool, default=False
        Whether to save particles lost at the lower y boundary

    warpx_save_particles_at_yhi: bool, default=False
        Whether to save particles lost at the upper y boundary

    warpx_save_particles_at_zlo: bool, default=False
        Whether to save particles lost at the lower z boundary

    warpx_save_particles_at_zhi: bool, default=False
        Whether to save particles lost at the upper z boundary

    warpx_save_particles_at_eb: bool, default=False
        Whether to save particles lost at the embedded boundary
    """
    def init(self, kw):

        if self.particle_type == 'electron':
            if self.charge is None: self.charge = '-q_e'
            if self.mass is None: self.mass = 'm_e'
        elif self.particle_type == 'positron':
            if self.charge is None: self.charge = 'q_e'
            if self.mass is None: self.mass = 'm_e'
        elif self.particle_type == 'proton':
            if self.charge is None: self.charge = 'q_e'
            if self.mass is None: self.mass = 'm_p'
        elif self.particle_type == 'anti-proton':
            if self.charge is None: self.charge = '-q_e'
            if self.mass is None: self.mass = 'm_p'
        else:
            if self.charge is None and self.charge_state is not None:
                if self.charge_state == +1.:
                    self.charge = 'q_e'
                elif self.charge_state == -1.:
                    self.charge = '-q_e'
                else:
                    self.charge = self.charge_state*constants.q_e
            if self.particle_type is not None:
                # Match a string of the format '#nXx', with the '#n' optional isotope number.
                m = re.match(r'(?P<iso>#[\d+])*(?P<sym>[A-Za-z]+)', self.particle_type)
                if m is not None:
                    element = periodictable.elements.symbol(m['sym'])
                    if m['iso'] is not None:
                        element = element[m['iso'][1:]]
                    if self.charge_state is not None:
                        assert self.charge_state <= element.number, Exception('%s charge state not valid'%self.particle_type)
                        try:
                            element = element.ion[self.charge_state]
                        except ValueError:
                            # Note that not all valid charge states are defined in elements,
                            # so this value error can be ignored.
                            pass
                    self.element = element
                    if self.mass is None:
                        self.mass = element.mass*periodictable.constants.atomic_mass_constant
                else:
                    raise Exception('The species "particle_type" is not known')

        self.boost_adjust_transverse_positions = kw.pop('warpx_boost_adjust_transverse_positions', None)

        # For the relativistic electrostatic solver
        self.self_fields_required_precision = kw.pop('warpx_self_fields_required_precision', None)
        self.self_fields_absolute_tolerance = kw.pop('warpx_self_fields_absolute_tolerance', None)
        self.self_fields_max_iters = kw.pop('warpx_self_fields_max_iters', None)
        self.self_fields_verbosity = kw.pop('warpx_self_fields_verbosity', None)
        self.save_previous_position = kw.pop('warpx_save_previous_position', None)
        self.do_not_deposit = kw.pop('warpx_do_not_deposit', None)
        self.do_not_push = kw.pop('warpx_do_not_push', None)
        self.do_not_gather = kw.pop('warpx_do_not_gather', None)
        self.random_theta = kw.pop('warpx_random_theta', None)

        # For particle reflection
        self.reflection_model_xlo = kw.pop('warpx_reflection_model_xlo', None)
        self.reflection_model_xhi = kw.pop('warpx_reflection_model_xhi', None)
        self.reflection_model_ylo = kw.pop('warpx_reflection_model_ylo', None)
        self.reflection_model_yhi = kw.pop('warpx_reflection_model_yhi', None)
        self.reflection_model_zlo = kw.pop('warpx_reflection_model_zlo', None)
        self.reflection_model_zhi = kw.pop('warpx_reflection_model_zhi', None)
        # self.reflection_model_eb = kw.pop('warpx_reflection_model_eb', None)

        # For the scraper buffer
        self.save_particles_at_xlo = kw.pop('warpx_save_particles_at_xlo', None)
        self.save_particles_at_xhi = kw.pop('warpx_save_particles_at_xhi', None)
        self.save_particles_at_ylo = kw.pop('warpx_save_particles_at_ylo', None)
        self.save_particles_at_yhi = kw.pop('warpx_save_particles_at_yhi', None)
        self.save_particles_at_zlo = kw.pop('warpx_save_particles_at_zlo', None)
        self.save_particles_at_zhi = kw.pop('warpx_save_particles_at_zhi', None)
        self.save_particles_at_eb = kw.pop('warpx_save_particles_at_eb', None)

    def initialize_inputs(self, layout,
                          initialize_self_fields = False,
                          injection_plane_position = None,
                          injection_plane_normal_vector = None):
        self.species_number = len(pywarpx.particles.species_names)

        if self.name is None:
            self.name = 'species{}'.format(self.species_number)

        pywarpx.particles.species_names.append(self.name)

        if initialize_self_fields is None:
            initialize_self_fields = False

        self.species = pywarpx.Bucket.Bucket(self.name,
                                             mass = self.mass,
                                             charge = self.charge,
                                             injection_style = None,
                                             initialize_self_fields = int(initialize_self_fields),
                                             boost_adjust_transverse_positions = self.boost_adjust_transverse_positions,
                                             self_fields_required_precision = self.self_fields_required_precision,
                                             self_fields_absolute_tolerance = self.self_fields_absolute_tolerance,
                                             self_fields_max_iters = self.self_fields_max_iters,
                                             self_fields_verbosity = self.self_fields_verbosity,
                                             save_particles_at_xlo = self.save_particles_at_xlo,
                                             save_particles_at_xhi = self.save_particles_at_xhi,
                                             save_particles_at_ylo = self.save_particles_at_ylo,
                                             save_particles_at_yhi = self.save_particles_at_yhi,
                                             save_particles_at_zlo = self.save_particles_at_zlo,
                                             save_particles_at_zhi = self.save_particles_at_zhi,
                                             save_particles_at_eb = self.save_particles_at_eb,
                                             save_previous_position = self.save_previous_position,
                                             do_not_deposit = self.do_not_deposit,
                                             do_not_push = self.do_not_push,
                                             do_not_gather = self.do_not_gather,
                                             random_theta = self.random_theta)

        # add reflection models
        self.species.add_new_attr("reflection_model_xlo(E)", self.reflection_model_xlo)
        self.species.add_new_attr("reflection_model_xhi(E)", self.reflection_model_xhi)
        self.species.add_new_attr("reflection_model_ylo(E)", self.reflection_model_ylo)
        self.species.add_new_attr("reflection_model_yhi(E)", self.reflection_model_yhi)
        self.species.add_new_attr("reflection_model_zlo(E)", self.reflection_model_zlo)
        self.species.add_new_attr("reflection_model_zhi(E)", self.reflection_model_zhi)
        # self.species.add_new_attr("reflection_model_eb(E)", self.reflection_model_eb)

        pywarpx.Particles.particles_list.append(self.species)

        if self.initial_distribution is not None:
            self.initial_distribution.initialize_inputs(self.species_number, layout, self.species, self.density_scale)

        if injection_plane_position is not None:
            if injection_plane_normal_vector is not None:
                assert injection_plane_normal_vector[0] == 0. and injection_plane_normal_vector[1] == 0.,\
                    Exception('Rigid injection can only be done along z')
            pywarpx.particles.rigid_injected_species.append(self.name)
            self.species.rigid_advance = 1
            self.species.zinject_plane = injection_plane_position


picmistandard.PICMI_MultiSpecies.Species_class = Species
class MultiSpecies(picmistandard.PICMI_MultiSpecies):
    def initialize_inputs(self, layout,
                          initialize_self_fields = False,
                          injection_plane_position = None,
                          injection_plane_normal_vector = None):
        for species in self.species_instances_list:
            species.initialize_inputs(layout,
                                      initialize_self_fields,
                                      injection_plane_position,
                                      injection_plane_normal_vector)


class GaussianBunchDistribution(picmistandard.PICMI_GaussianBunchDistribution):
    def initialize_inputs(self, species_number, layout, species, density_scale):
        species.injection_style = "gaussian_beam"
        species.x_m = self.centroid_position[0]
        species.y_m = self.centroid_position[1]
        species.z_m = self.centroid_position[2]
        species.x_rms = self.rms_bunch_size[0]
        species.y_rms = self.rms_bunch_size[1]
        species.z_rms = self.rms_bunch_size[2]

        # --- Only PseudoRandomLayout is supported
        species.npart = layout.n_macroparticles

        # --- Calculate the total charge. Note that charge might be a string instead of a number.
        charge = species.charge
        if charge == 'q_e' or charge == '+q_e':
            charge = constants.q_e
        elif charge == '-q_e':
            charge = -constants.q_e
        species.q_tot = self.n_physical_particles*charge
        if density_scale is not None:
            species.q_tot *= density_scale

        # --- The PICMI standard doesn't yet have a way of specifying these values.
        # --- They should default to the size of the domain. They are not typically
        # --- necessary though since any particles outside the domain are rejected.
        #species.xmin
        #species.xmax
        #species.ymin
        #species.ymax
        #species.zmin
        #species.zmax

        # --- Note that WarpX takes gamma*beta as input
        if np.any(np.not_equal(self.velocity_divergence, 0.)):
            species.momentum_distribution_type = "radial_expansion"
            species.u_over_r = self.velocity_divergence[0]/constants.c
            #species.u_over_y = self.velocity_divergence[1]/constants.c
            #species.u_over_z = self.velocity_divergence[2]/constants.c
        elif np.any(np.not_equal(self.rms_velocity, 0.)):
            species.momentum_distribution_type = "gaussian"
            species.ux_m = self.centroid_velocity[0]/constants.c
            species.uy_m = self.centroid_velocity[1]/constants.c
            species.uz_m = self.centroid_velocity[2]/constants.c
            species.ux_th = self.rms_velocity[0]/constants.c
            species.uy_th = self.rms_velocity[1]/constants.c
            species.uz_th = self.rms_velocity[2]/constants.c
        else:
            species.momentum_distribution_type = "constant"
            species.ux = self.centroid_velocity[0]/constants.c
            species.uy = self.centroid_velocity[1]/constants.c
            species.uz = self.centroid_velocity[2]/constants.c


class DensityDistributionBase(object):
    """This is a base class for several predefined density distributions. It
    captures universal initialization logic."""

    def set_mangle_dict(self):
        if not hasattr(self, 'mangle_dict'):
            self.mangle_dict = None

        if hasattr(self, "user_defined_kw") and self.mangle_dict is None:
            # Only do this once so that the same variables can be used multiple
            # times
            self.mangle_dict = pywarpx.my_constants.add_keywords(self.user_defined_kw)

    def set_species_attributes(self, species, layout):
        if isinstance(layout, GriddedLayout):
            # --- Note that the grid attribute of GriddedLayout is ignored
            species.injection_style = "nuniformpercell"
            species.num_particles_per_cell_each_dim = layout.n_macroparticle_per_cell
        elif isinstance(layout, PseudoRandomLayout):
            assert (layout.n_macroparticles_per_cell is not None), Exception('WarpX only supports n_macroparticles_per_cell for the PseudoRandomLayout with this distribution')
            species.injection_style = "nrandompercell"
            species.num_particles_per_cell = layout.n_macroparticles_per_cell
        else:
            raise Exception('WarpX does not support the specified layout for this distribution')

        species.xmin = self.lower_bound[0]
        species.xmax = self.upper_bound[0]
        species.ymin = self.lower_bound[1]
        species.ymax = self.upper_bound[1]
        species.zmin = self.lower_bound[2]
        species.zmax = self.upper_bound[2]

        if self.fill_in:
            species.do_continuous_injection = 1

        # --- Note that WarpX takes gamma*beta as input
        if (hasattr(self, "momentum_expressions")
            and np.any(np.not_equal(self.momentum_expressions, None))
        ):
            species.momentum_distribution_type = 'parse_momentum_function'
            self.setup_parse_momentum_functions(species)
        elif np.any(np.not_equal(self.rms_velocity, 0.)):
            species.momentum_distribution_type = "gaussian"
            species.ux_m = self.directed_velocity[0]/constants.c
            species.uy_m = self.directed_velocity[1]/constants.c
            species.uz_m = self.directed_velocity[2]/constants.c
            species.ux_th = self.rms_velocity[0]/constants.c
            species.uy_th = self.rms_velocity[1]/constants.c
            species.uz_th = self.rms_velocity[2]/constants.c
        else:
            species.momentum_distribution_type = "constant"
            species.ux = self.directed_velocity[0]/constants.c
            species.uy = self.directed_velocity[1]/constants.c
            species.uz = self.directed_velocity[2]/constants.c

    def setup_parse_momentum_functions(self, species):
        for sdir, idir in zip(['x', 'y', 'z'], [0, 1, 2]):
            if self.momentum_expressions[idir] is not None:
                expression = pywarpx.my_constants.mangle_expression(self.momentum_expressions[idir], self.mangle_dict)
            else:
                expression = f'{self.directed_velocity[idir]}'
            species.__setattr__(f'momentum_function_u{sdir}(x,y,z)', f'({expression})/{constants.c}')


class UniformFluxDistribution(picmistandard.PICMI_UniformFluxDistribution, DensityDistributionBase):
    def initialize_inputs(self, species_number, layout, species, density_scale):

        self.fill_in = False
        self.set_mangle_dict()
        self.set_species_attributes(species, layout)

        species.profile = "constant"
        species.density = self.flux
        if density_scale is not None:
            species.density *= density_scale
        species.flux_normal_axis = self.flux_normal_axis
        species.surface_flux_pos = self.surface_flux_position
        species.flux_direction = self.flux_direction
        species.flux_tmin = self.flux_tmin
        species.flux_tmax = self.flux_tmax

        # --- Use specific attributes for flux injection
        species.injection_style = "nfluxpercell"
        assert (isinstance(layout, PseudoRandomLayout)), Exception('UniformFluxDistribution only supports the PseudoRandomLayout in WarpX')
        if species.momentum_distribution_type == "gaussian":
            species.momentum_distribution_type = "gaussianflux"


class UniformDistribution(picmistandard.PICMI_UniformDistribution, DensityDistributionBase):
    def initialize_inputs(self, species_number, layout, species, density_scale):

        self.set_mangle_dict()
        self.set_species_attributes(species, layout)

        # --- Only constant density is supported by this class
        species.profile = "constant"
        species.density = self.density
        if density_scale is not None:
            species.density *= density_scale


class AnalyticDistribution(picmistandard.PICMI_AnalyticDistribution, DensityDistributionBase):
    def initialize_inputs(self, species_number, layout, species, density_scale):

        self.set_mangle_dict()
        self.set_species_attributes(species, layout)

        species.profile = "parse_density_function"
        expression = pywarpx.my_constants.mangle_expression(self.density_expression, self.mangle_dict)
        if density_scale is None:
            species.__setattr__('density_function(x,y,z)', expression)
        else:
            species.__setattr__('density_function(x,y,z)', "{}*({})".format(density_scale, expression))


class ParticleListDistribution(picmistandard.PICMI_ParticleListDistribution):
    def init(self, kw):
        pass

    def initialize_inputs(self, species_number, layout, species, density_scale):

        species.injection_style = "multipleparticles"
        species.multiple_particles_pos_x = self.x
        species.multiple_particles_pos_y = self.y
        species.multiple_particles_pos_z = self.z
        species.multiple_particles_ux = np.array(self.ux)/constants.c
        species.multiple_particles_uy = np.array(self.uy)/constants.c
        species.multiple_particles_uz = np.array(self.uz)/constants.c
        species.multiple_particles_weight = self.weight
        if density_scale is not None:
            species.multiple_particles_weight = self.weight*density_scale


class ParticleDistributionPlanarInjector(picmistandard.PICMI_ParticleDistributionPlanarInjector):
    pass


class GriddedLayout(picmistandard.PICMI_GriddedLayout):
    pass


class PseudoRandomLayout(picmistandard.PICMI_PseudoRandomLayout):
    def init(self, kw):
        if self.seed is not None:
            print('Warning: WarpX does not support specifying the random number seed in PseudoRandomLayout')


class BinomialSmoother(picmistandard.PICMI_BinomialSmoother):

    def initialize_inputs(self, solver):
        pywarpx.warpx.use_filter = 1
        pywarpx.warpx.use_filter_compensation = bool(np.all(self.compensation))
        if self.n_pass is None:
            # If not specified, do at least one pass in each direction.
            self.n_pass = 1
        try:
            # Check if n_pass is a vector
            len(self.n_pass)
        except TypeError:
            # If not, make it a vector
            self.n_pass = solver.grid.number_of_dimensions*[self.n_pass]
        pywarpx.warpx.filter_npass_each_dir = self.n_pass


class CylindricalGrid(picmistandard.PICMI_CylindricalGrid):
    """
    This assumes that WarpX was compiled with USE_RZ = TRUE

    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ---------
    warpx_max_grid_size: integer, default=32
       Maximum block size in either direction

    warpx_max_grid_size_x: integer, optional
       Maximum block size in radial direction

    warpx_max_grid_size_y: integer, optional
       Maximum block size in longitudinal direction

    warpx_blocking_factor: integer, optional
       Blocking factor (which controls the block size)

    warpx_blocking_factor_x: integer, optional
       Blocking factor (which controls the block size) in the radial direction

    warpx_blocking_factor_y: integer, optional
       Blocking factor (which controls the block size) in the longitudinal direction

    warpx_potential_lo_r: float, default=0.
       Electrostatic potential on the lower radial boundary

    warpx_potential_hi_r: float, default=0.
       Electrostatic potential on the upper radial boundary

    warpx_potential_lo_z: float, default=0.
       Electrostatic potential on the lower longitudinal boundary

    warpx_potential_hi_z: float, default=0.
       Electrostatic potential on the upper longitudinal boundary

    warpx_reflect_all_velocities: bool default=False
        Whether the sign of all of the particle velocities are changed upon
        reflection on a boundary, or only the velocity normal to the surface
    """
    def init(self, kw):
        self.max_grid_size = kw.pop('warpx_max_grid_size', 32)
        self.max_grid_size_x = kw.pop('warpx_max_grid_size_x', None)
        self.max_grid_size_y = kw.pop('warpx_max_grid_size_y', None)
        self.blocking_factor = kw.pop('warpx_blocking_factor', None)
        self.blocking_factor_x = kw.pop('warpx_blocking_factor_x', None)
        self.blocking_factor_y = kw.pop('warpx_blocking_factor_y', None)

        self.potential_xmin = kw.pop('warpx_potential_lo_r', None)
        self.potential_xmax = kw.pop('warpx_potential_hi_r', None)
        self.potential_ymin = None
        self.potential_ymax = None
        self.potential_zmin = kw.pop('warpx_potential_lo_z', None)
        self.potential_zmax = kw.pop('warpx_potential_hi_z', None)
        self.reflect_all_velocities = kw.pop('warpx_reflect_all_velocities', None)

        # Geometry
        # Set these as soon as the information is available
        # (since these are needed to determine which shared object to load)
        pywarpx.geometry.dims = 'RZ'
        pywarpx.geometry.prob_lo = self.lower_bound  # physical domain
        pywarpx.geometry.prob_hi = self.upper_bound

    def initialize_inputs(self):
        pywarpx.amr.n_cell = self.number_of_cells

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = self.max_grid_size
        pywarpx.amr.max_grid_size_x = self.max_grid_size_x
        pywarpx.amr.max_grid_size_y = self.max_grid_size_y
        pywarpx.amr.blocking_factor = self.blocking_factor
        pywarpx.amr.blocking_factor_x = self.blocking_factor_x
        pywarpx.amr.blocking_factor_y = self.blocking_factor_y

        assert self.lower_bound[0] >= 0., Exception('Lower radial boundary must be >= 0.')
        assert self.lower_boundary_conditions[0] != 'periodic' and self.upper_boundary_conditions[0] != 'periodic', Exception('Radial boundaries can not be periodic')

        pywarpx.warpx.n_rz_azimuthal_modes = self.n_azimuthal_modes

        # Boundary conditions
        pywarpx.boundary.field_lo = [BC_map[bc] for bc in self.lower_boundary_conditions]
        pywarpx.boundary.field_hi = [BC_map[bc] for bc in self.upper_boundary_conditions]
        pywarpx.boundary.particle_lo = self.lower_boundary_conditions_particles
        pywarpx.boundary.particle_hi = self.upper_boundary_conditions_particles
        pywarpx.boundary.reflect_all_velocities = self.reflect_all_velocities

        if self.moving_window_velocity is not None and np.any(np.not_equal(self.moving_window_velocity, 0.)):
            pywarpx.warpx.do_moving_window = 1
            if self.moving_window_velocity[0] != 0.:
                pywarpx.warpx.moving_window_dir = 'r'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[0]/constants.c  # in units of the speed of light
            if self.moving_window_velocity[1] != 0.:
                pywarpx.warpx.moving_window_dir = 'z'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[1]/constants.c  # in units of the speed of light

        if self.refined_regions:
            assert len(self.refined_regions) == 1, Exception('WarpX only supports one refined region.')
            assert self.refined_regions[0][0] == 1, Exception('The one refined region can only be level 1')
            pywarpx.amr.max_level = 1
            pywarpx.warpx.fine_tag_lo = self.refined_regions[0][1]
            pywarpx.warpx.fine_tag_hi = self.refined_regions[0][2]
            # The refinement_factor is ignored (assumed to be [2,2])
        else:
            pywarpx.amr.max_level = 0


class Cartesian1DGrid(picmistandard.PICMI_Cartesian1DGrid):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ---------
    warpx_max_grid_size: integer, default=32
       Maximum block size in either direction

    warpx_max_grid_size_x: integer, optional
       Maximum block size in longitudinal direction

    warpx_blocking_factor: integer, optional
       Blocking factor (which controls the block size)

    warpx_blocking_factor_x: integer, optional
       Blocking factor (which controls the block size) in the longitudinal direction

    warpx_potential_lo_z: float, default=0.
       Electrostatic potential on the lower longitudinal boundary

    warpx_potential_hi_z: float, default=0.
       Electrostatic potential on the upper longitudinal boundary
    """
    def init(self, kw):
        self.max_grid_size = kw.pop('warpx_max_grid_size', 32)
        self.max_grid_size_x = kw.pop('warpx_max_grid_size_x', None)
        self.blocking_factor = kw.pop('warpx_blocking_factor', None)
        self.blocking_factor_x = kw.pop('warpx_blocking_factor_x', None)

        self.potential_xmin = None
        self.potential_xmax = None
        self.potential_ymin = None
        self.potential_ymax = None
        self.potential_zmin = kw.pop('warpx_potential_lo_z', None)
        self.potential_zmax = kw.pop('warpx_potential_hi_z', None)

        # Geometry
        # Set these as soon as the information is available
        # (since these are needed to determine which shared object to load)
        pywarpx.geometry.dims = '1'
        pywarpx.geometry.prob_lo = self.lower_bound  # physical domain
        pywarpx.geometry.prob_hi = self.upper_bound

    def initialize_inputs(self):
        pywarpx.amr.n_cell = self.number_of_cells

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = self.max_grid_size
        pywarpx.amr.max_grid_size_x = self.max_grid_size_x
        pywarpx.amr.blocking_factor = self.blocking_factor
        pywarpx.amr.blocking_factor_x = self.blocking_factor_x

        # Boundary conditions
        pywarpx.boundary.field_lo = [BC_map[bc] for bc in self.lower_boundary_conditions]
        pywarpx.boundary.field_hi = [BC_map[bc] for bc in self.upper_boundary_conditions]
        pywarpx.boundary.particle_lo = self.lower_boundary_conditions_particles
        pywarpx.boundary.particle_hi = self.upper_boundary_conditions_particles

        if self.moving_window_velocity is not None and np.any(np.not_equal(self.moving_window_velocity, 0.)):
            pywarpx.warpx.do_moving_window = 1
            if self.moving_window_velocity[0] != 0.:
                pywarpx.warpx.moving_window_dir = 'z'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[0]/constants.c  # in units of the speed of light

        if self.refined_regions:
            assert len(self.refined_regions) == 1, Exception('WarpX only supports one refined region.')
            assert self.refined_regions[0][0] == 1, Exception('The one refined region can only be level 1')
            pywarpx.amr.max_level = 1
            pywarpx.warpx.fine_tag_lo = self.refined_regions[0][1]
            pywarpx.warpx.fine_tag_hi = self.refined_regions[0][2]
            # The refinement_factor is ignored (assumed to be [2,2])
        else:
            pywarpx.amr.max_level = 0

class Cartesian2DGrid(picmistandard.PICMI_Cartesian2DGrid):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ---------
    warpx_max_grid_size: integer, default=32
       Maximum block size in either direction

    warpx_max_grid_size_x: integer, optional
       Maximum block size in x direction

    warpx_max_grid_size_y: integer, optional
       Maximum block size in z direction

    warpx_blocking_factor: integer, optional
       Blocking factor (which controls the block size)

    warpx_blocking_factor_x: integer, optional
       Blocking factor (which controls the block size) in the x direction

    warpx_blocking_factor_y: integer, optional
       Blocking factor (which controls the block size) in the z direction

    warpx_potential_lo_x: float, default=0.
       Electrostatic potential on the lower x boundary

    warpx_potential_hi_x: float, default=0.
       Electrostatic potential on the upper x boundary

    warpx_potential_lo_z: float, default=0.
       Electrostatic potential on the lower z boundary

    warpx_potential_hi_z: float, default=0.
       Electrostatic potential on the upper z boundary
    """
    def init(self, kw):
        self.max_grid_size = kw.pop('warpx_max_grid_size', 32)
        self.max_grid_size_x = kw.pop('warpx_max_grid_size_x', None)
        self.max_grid_size_y = kw.pop('warpx_max_grid_size_y', None)
        self.blocking_factor = kw.pop('warpx_blocking_factor', None)
        self.blocking_factor_x = kw.pop('warpx_blocking_factor_x', None)
        self.blocking_factor_y = kw.pop('warpx_blocking_factor_y', None)

        self.potential_xmin = kw.pop('warpx_potential_lo_x', None)
        self.potential_xmax = kw.pop('warpx_potential_hi_x', None)
        self.potential_ymin = None
        self.potential_ymax = None
        self.potential_zmin = kw.pop('warpx_potential_lo_z', None)
        self.potential_zmax = kw.pop('warpx_potential_hi_z', None)

        # Geometry
        # Set these as soon as the information is available
        # (since these are needed to determine which shared object to load)
        pywarpx.geometry.dims = '2'
        pywarpx.geometry.prob_lo = self.lower_bound  # physical domain
        pywarpx.geometry.prob_hi = self.upper_bound

    def initialize_inputs(self):
        pywarpx.amr.n_cell = self.number_of_cells

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = self.max_grid_size
        pywarpx.amr.max_grid_size_x = self.max_grid_size_x
        pywarpx.amr.max_grid_size_y = self.max_grid_size_y
        pywarpx.amr.blocking_factor = self.blocking_factor
        pywarpx.amr.blocking_factor_x = self.blocking_factor_x
        pywarpx.amr.blocking_factor_y = self.blocking_factor_y

        # Boundary conditions
        pywarpx.boundary.field_lo = [BC_map[bc] for bc in self.lower_boundary_conditions]
        pywarpx.boundary.field_hi = [BC_map[bc] for bc in self.upper_boundary_conditions]
        pywarpx.boundary.particle_lo = self.lower_boundary_conditions_particles
        pywarpx.boundary.particle_hi = self.upper_boundary_conditions_particles

        if self.moving_window_velocity is not None and np.any(np.not_equal(self.moving_window_velocity, 0.)):
            pywarpx.warpx.do_moving_window = 1
            if self.moving_window_velocity[0] != 0.:
                pywarpx.warpx.moving_window_dir = 'x'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[0]/constants.c  # in units of the speed of light
            if self.moving_window_velocity[1] != 0.:
                pywarpx.warpx.moving_window_dir = 'z'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[1]/constants.c  # in units of the speed of light

        if self.refined_regions:
            assert len(self.refined_regions) == 1, Exception('WarpX only supports one refined region.')
            assert self.refined_regions[0][0] == 1, Exception('The one refined region can only be level 1')
            pywarpx.amr.max_level = 1
            pywarpx.warpx.fine_tag_lo = self.refined_regions[0][1]
            pywarpx.warpx.fine_tag_hi = self.refined_regions[0][2]
            # The refinement_factor is ignored (assumed to be [2,2])
        else:
            pywarpx.amr.max_level = 0


class Cartesian3DGrid(picmistandard.PICMI_Cartesian3DGrid):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ---------
    warpx_max_grid_size: integer, default=32
       Maximum block size in either direction

    warpx_max_grid_size_x: integer, optional
       Maximum block size in x direction

    warpx_max_grid_size_y: integer, optional
       Maximum block size in z direction

    warpx_max_grid_size_z: integer, optional
       Maximum block size in z direction

    warpx_blocking_factor: integer, optional
       Blocking factor (which controls the block size)

    warpx_blocking_factor_x: integer, optional
       Blocking factor (which controls the block size) in the x direction

    warpx_blocking_factor_y: integer, optional
       Blocking factor (which controls the block size) in the z direction

    warpx_blocking_factor_z: integer, optional
       Blocking factor (which controls the block size) in the z direction

    warpx_potential_lo_x: float, default=0.
       Electrostatic potential on the lower x boundary

    warpx_potential_hi_x: float, default=0.
       Electrostatic potential on the upper x boundary

    warpx_potential_lo_y: float, default=0.
       Electrostatic potential on the lower z boundary

    warpx_potential_hi_y: float, default=0.
       Electrostatic potential on the upper z boundary

    warpx_potential_lo_z: float, default=0.
       Electrostatic potential on the lower z boundary

    warpx_potential_hi_z: float, default=0.
       Electrostatic potential on the upper z boundary
    """
    def init(self, kw):
        self.max_grid_size = kw.pop('warpx_max_grid_size', 32)
        self.max_grid_size_x = kw.pop('warpx_max_grid_size_x', None)
        self.max_grid_size_y = kw.pop('warpx_max_grid_size_y', None)
        self.max_grid_size_z = kw.pop('warpx_max_grid_size_z', None)
        self.blocking_factor = kw.pop('warpx_blocking_factor', None)
        self.blocking_factor_x = kw.pop('warpx_blocking_factor_x', None)
        self.blocking_factor_y = kw.pop('warpx_blocking_factor_y', None)
        self.blocking_factor_z = kw.pop('warpx_blocking_factor_z', None)

        self.potential_xmin = kw.pop('warpx_potential_lo_x', None)
        self.potential_xmax = kw.pop('warpx_potential_hi_x', None)
        self.potential_ymin = kw.pop('warpx_potential_lo_y', None)
        self.potential_ymax = kw.pop('warpx_potential_hi_y', None)
        self.potential_zmin = kw.pop('warpx_potential_lo_z', None)
        self.potential_zmax = kw.pop('warpx_potential_hi_z', None)

        # Geometry
        # Set these as soon as the information is available
        # (since these are needed to determine which shared object to load)
        pywarpx.geometry.dims = '3'
        pywarpx.geometry.prob_lo = self.lower_bound  # physical domain
        pywarpx.geometry.prob_hi = self.upper_bound

    def initialize_inputs(self):
        pywarpx.amr.n_cell = self.number_of_cells

        # Maximum allowable size of each subdomain in the problem domain;
        #    this is used to decompose the domain for parallel calculations.
        pywarpx.amr.max_grid_size = self.max_grid_size
        pywarpx.amr.max_grid_size_x = self.max_grid_size_x
        pywarpx.amr.max_grid_size_y = self.max_grid_size_y
        pywarpx.amr.max_grid_size_z = self.max_grid_size_z
        pywarpx.amr.blocking_factor = self.blocking_factor
        pywarpx.amr.blocking_factor_x = self.blocking_factor_x
        pywarpx.amr.blocking_factor_y = self.blocking_factor_y
        pywarpx.amr.blocking_factor_z = self.blocking_factor_z

        # Boundary conditions
        pywarpx.boundary.field_lo = [BC_map[bc] for bc in self.lower_boundary_conditions]
        pywarpx.boundary.field_hi = [BC_map[bc] for bc in self.upper_boundary_conditions]
        pywarpx.boundary.particle_lo = self.lower_boundary_conditions_particles
        pywarpx.boundary.particle_hi = self.upper_boundary_conditions_particles

        if self.moving_window_velocity is not None and np.any(np.not_equal(self.moving_window_velocity, 0.)):
            pywarpx.warpx.do_moving_window = 1
            if self.moving_window_velocity[0] != 0.:
                pywarpx.warpx.moving_window_dir = 'x'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[0]/constants.c  # in units of the speed of light
            if self.moving_window_velocity[1] != 0.:
                pywarpx.warpx.moving_window_dir = 'y'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[1]/constants.c  # in units of the speed of light
            if self.moving_window_velocity[2] != 0.:
                pywarpx.warpx.moving_window_dir = 'z'
                pywarpx.warpx.moving_window_v = self.moving_window_velocity[2]/constants.c  # in units of the speed of light

        if self.refined_regions:
            assert len(self.refined_regions) == 1, Exception('WarpX only supports one refined region.')
            assert self.refined_regions[0][0] == 1, Exception('The one refined region can only be level 1')
            pywarpx.amr.max_level = 1
            pywarpx.warpx.fine_tag_lo = self.refined_regions[0][1]
            pywarpx.warpx.fine_tag_hi = self.refined_regions[0][2]
            # The refinement_factor is ignored (assumed to be [2,2,2])
        else:
            pywarpx.amr.max_level = 0

class ElectromagneticSolver(picmistandard.PICMI_ElectromagneticSolver):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ----------
    warpx_pml_ncell: integer, optional
        The depth of the PML, in number of cells

    warpx_periodic_single_box_fft: bool, default=False
        Whether to do the spectral solver FFTs assuming a single
        simulation block

    warpx_current_correction: bool, default=True
        Whether to do the current correction for the spectral solver.
        See documentation for exceptions to the default value.

    warpx_psatd_update_with_rho: bool, optional
        Whether to update with the actual rho for the spectral solver

    warpx_psatd_do_time_averaging: bool, optional
        Whether to do the time averaging for the spectral solver

    warpx_psatd_J_in_time: {'constant', 'linear'}, default='constant'
        This determines whether the current density is assumed to be constant
        or linear in time, within the time step over which the electromagnetic
        fields are evolved.

    warpx_psatd_rho_in_time: {'linear'}, default='linear'
        This determines whether the charge density is assumed to be linear
        in time, within the time step over which the electromagnetic fields are evolved.

    warpx_do_pml_in_domain: bool, default=False
        Whether to do the PML boundaries within the domain (versus
        in the guard cells)

    warpx_pml_has_particles: bool, default=False
        Whether to allow particles in the PML region

    warpx_do_pml_j_damping: bool, default=False
        Whether to do damping of J in the PML
    """
    def init(self, kw):
        assert self.method is None or self.method in ['Yee', 'CKC', 'PSATD', 'ECT'], Exception("Only 'Yee', 'CKC', 'PSATD', and 'ECT' are supported")

        self.pml_ncell = kw.pop('warpx_pml_ncell', None)

        if self.method == 'PSATD':
            self.psatd_periodic_single_box_fft = kw.pop('warpx_periodic_single_box_fft', None)
            self.psatd_current_correction = kw.pop('warpx_current_correction', None)
            self.psatd_update_with_rho = kw.pop('warpx_psatd_update_with_rho', None)
            self.psatd_do_time_averaging = kw.pop('warpx_psatd_do_time_averaging', None)
            self.psatd_J_in_time = kw.pop('warpx_psatd_J_in_time', None)
            self.psatd_rho_in_time = kw.pop('warpx_psatd_rho_in_time', None)

        self.do_pml_in_domain = kw.pop('warpx_do_pml_in_domain', None)
        self.pml_has_particles = kw.pop('warpx_pml_has_particles', None)
        self.do_pml_j_damping = kw.pop('warpx_do_pml_j_damping', None)

    def initialize_inputs(self):

        self.grid.initialize_inputs()

        pywarpx.warpx.pml_ncell = self.pml_ncell

        if self.method == 'PSATD':
            pywarpx.psatd.periodic_single_box_fft = self.psatd_periodic_single_box_fft
            pywarpx.psatd.current_correction = self.psatd_current_correction
            pywarpx.psatd.update_with_rho = self.psatd_update_with_rho
            pywarpx.psatd.do_time_averaging = self.psatd_do_time_averaging
            pywarpx.psatd.J_in_time = self.psatd_J_in_time
            pywarpx.psatd.rho_in_time = self.psatd_rho_in_time

            if self.grid.guard_cells is not None:
                pywarpx.psatd.nx_guard = self.grid.guard_cells[0]
                if self.grid.number_of_dimensions == 3:
                    pywarpx.psatd.ny_guard = self.grid.guard_cells[1]
                pywarpx.psatd.nz_guard = self.grid.guard_cells[-1]

            if self.stencil_order is not None:
                pywarpx.psatd.nox = self.stencil_order[0]
                if self.grid.number_of_dimensions == 3:
                    pywarpx.psatd.noy = self.stencil_order[1]
                pywarpx.psatd.noz = self.stencil_order[-1]

            if self.galilean_velocity is not None:
                if self.grid.number_of_dimensions == 2:
                    self.galilean_velocity = [self.galilean_velocity[0], 0., self.galilean_velocity[1]]
                pywarpx.psatd.v_galilean = np.array(self.galilean_velocity)/constants.c

        # --- Same method names are used, though mapped to lower case.
        pywarpx.algo.maxwell_solver = self.method

        if self.cfl is not None:
            pywarpx.warpx.cfl = self.cfl

        if self.source_smoother is not None:
            self.source_smoother.initialize_inputs(self)

        pywarpx.warpx.do_dive_cleaning = self.divE_cleaning
        pywarpx.warpx.do_divb_cleaning = self.divB_cleaning

        pywarpx.warpx.do_pml_dive_cleaning = self.pml_divE_cleaning
        pywarpx.warpx.do_pml_divb_cleaning = self.pml_divB_cleaning

        pywarpx.warpx.do_pml_in_domain = self.do_pml_in_domain
        pywarpx.warpx.pml_has_particles = self.pml_has_particles
        pywarpx.warpx.do_pml_j_damping = self.do_pml_j_damping

class ElectrostaticSolver(picmistandard.PICMI_ElectrostaticSolver):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ----------
    warpx_relativistic: bool, default=False
        Whether to use the relativistic solver or lab frame solver

    warpx_absolute_tolerance: float, default=0.
        Absolute tolerance on the lab fram solver

    warpx_self_fields_verbosity: integer, default=2
        Level of verbosity for the lab frame solver
    """
    def init(self, kw):
        self.relativistic = kw.pop('warpx_relativistic', False)
        self.absolute_tolerance = kw.pop('warpx_absolute_tolerance', None)
        self.self_fields_verbosity = kw.pop('warpx_self_fields_verbosity', None)
        self.magnetostatic = kw.pop('warpx_magnetostatic', False)

    def initialize_inputs(self):

        self.grid.initialize_inputs()

        if self.relativistic:
            pywarpx.warpx.do_electrostatic = 'relativistic'
        else:
            if self.magnetostatic:
                pywarpx.warpx.do_electrostatic = 'labframe-electromagnetostatic'
            else:
                pywarpx.warpx.do_electrostatic = 'labframe'
            pywarpx.warpx.self_fields_required_precision = self.required_precision
            pywarpx.warpx.self_fields_absolute_tolerance = self.absolute_tolerance
            pywarpx.warpx.self_fields_max_iters = self.maximum_iterations
            pywarpx.warpx.self_fields_verbosity = self.self_fields_verbosity
            pywarpx.boundary.potential_lo_x = self.grid.potential_xmin
            pywarpx.boundary.potential_lo_y = self.grid.potential_ymin
            pywarpx.boundary.potential_lo_z = self.grid.potential_zmin
            pywarpx.boundary.potential_hi_x = self.grid.potential_xmax
            pywarpx.boundary.potential_hi_y = self.grid.potential_ymax
            pywarpx.boundary.potential_hi_z = self.grid.potential_zmax


class GaussianLaser(picmistandard.PICMI_GaussianLaser):
    def initialize_inputs(self):
        self.laser_number = len(pywarpx.lasers.names) + 1
        if self.name is None:
            self.name = 'laser{}'.format(self.laser_number)

        self.laser = pywarpx.Lasers.newlaser(self.name)

        self.laser.profile = "Gaussian"
        self.laser.wavelength = self.wavelength  # The wavelength of the laser (in meters)
        self.laser.e_max = self.E0  # Maximum amplitude of the laser field (in V/m)
        self.laser.polarization = self.polarization_direction  # The main polarization vector
        self.laser.profile_waist = self.waist  # The waist of the laser (in meters)
        self.laser.profile_duration = self.duration  # The duration of the laser (in seconds)
        self.laser.direction = self.propagation_direction
        self.laser.zeta = self.zeta
        self.laser.beta = self.beta
        self.laser.phi2 = self.phi2
        self.laser.phi0 = self.phi0

        self.laser.do_continuous_injection = self.fill_in

class AnalyticLaser(picmistandard.PICMI_AnalyticLaser):
    def init(self, kw):
        self.mangle_dict = None

    def initialize_inputs(self):
        self.laser_number = len(pywarpx.lasers.names) + 1
        if self.name is None:
            self.name = 'laser{}'.format(self.laser_number)

        self.laser = pywarpx.Lasers.newlaser(self.name)

        self.laser.profile = "parse_field_function"
        self.laser.wavelength = self.wavelength  # The wavelength of the laser (in meters)
        self.laser.e_max = self.Emax  # Maximum amplitude of the laser field (in V/m)
        self.laser.polarization = self.polarization_direction  # The main polarization vector
        self.laser.direction = self.propagation_direction
        self.laser.do_continuous_injection = self.fill_in

        if self.mangle_dict is None:
            # Only do this once so that the same variables are used in this distribution
            # is used multiple times
            self.mangle_dict = pywarpx.my_constants.add_keywords(self.user_defined_kw)
        expression = pywarpx.my_constants.mangle_expression(self.field_expression, self.mangle_dict)
        self.laser.__setattr__('field_function(X,Y,t)', expression)


class LaserAntenna(picmistandard.PICMI_LaserAntenna):
    def initialize_inputs(self, laser):
        laser.laser.position = self.position  # This point is on the laser plane
        if (
            self.normal_vector is not None
            and not np.allclose(laser.laser.direction, self.normal_vector)
        ):
            raise AttributeError(
                'The specified laser direction does not match the '
                'specified antenna normal.'
            )
        self.normal_vector = laser.laser.direction # The plane normal direction
        if isinstance(laser, GaussianLaser):
            # Focal distance from the antenna (in meters)
            laser.laser.profile_focal_distance = np.sqrt(
                (laser.focal_position[0] - self.position[0])**2 +
                (laser.focal_position[1] - self.position[1])**2 +
                (laser.focal_position[2] - self.position[2])**2
            )
            # The time at which the laser reaches its peak (in seconds)
            laser.laser.profile_t_peak = np.sqrt(
                (self.position[0] - laser.centroid_position[0])**2 +
                (self.position[1] - laser.centroid_position[1])**2 +
                (self.position[2] - laser.centroid_position[2])**2
            ) / constants.c


class AnalyticInitialField(picmistandard.PICMI_AnalyticAppliedField):
    def init(self, kw):
        self.mangle_dict = None

    def initialize_inputs(self):
        # Note that lower and upper_bound are not used by WarpX

        if self.mangle_dict is None:
            # Only do this once so that the same variables are used in this distribution
            # is used multiple times
            self.mangle_dict = pywarpx.my_constants.add_keywords(self.user_defined_kw)

        if (self.Ex_expression is not None or
            self.Ey_expression is not None or
            self.Ez_expression is not None):
            pywarpx.warpx.E_ext_grid_init_style = 'parse_e_ext_grid_function'
            for sdir, expression in zip(['x', 'y', 'z'], [self.Ex_expression, self.Ey_expression, self.Ez_expression]):
                expression = pywarpx.my_constants.mangle_expression(expression, self.mangle_dict)
                pywarpx.warpx.__setattr__(f'E{sdir}_external_grid_function(x,y,z)', expression)

        if (self.Bx_expression is not None or
            self.By_expression is not None or
            self.Bz_expression is not None):
            pywarpx.warpx.B_ext_grid_init_style = 'parse_b_ext_grid_function'
            for sdir, expression in zip(['x', 'y', 'z'], [self.Bx_expression, self.By_expression, self.Bz_expression]):
                expression = pywarpx.my_constants.mangle_expression(expression, self.mangle_dict)
                pywarpx.warpx.__setattr__(f'B{sdir}_external_grid_function(x,y,z)', expression)


class ConstantAppliedField(picmistandard.PICMI_ConstantAppliedField):
    def initialize_inputs(self):
        # Note that lower and upper_bound are not used by WarpX

        if (self.Ex is not None or
            self.Ey is not None or
            self.Ez is not None):
            pywarpx.particles.E_ext_particle_init_style = 'constant'
            pywarpx.particles.E_external_particle = [self.Ex or 0., self.Ey or 0., self.Ez or 0.]

        if (self.Bx is not None or
            self.By is not None or
            self.Bz is not None):
            pywarpx.particles.B_ext_particle_init_style = 'constant'
            pywarpx.particles.B_external_particle = [self.Bx or 0., self.By or 0., self.Bz or 0.]


class AnalyticAppliedField(picmistandard.PICMI_AnalyticAppliedField):
    def init(self, kw):
        self.mangle_dict = None

    def initialize_inputs(self):
        # Note that lower and upper_bound are not used by WarpX

        if self.mangle_dict is None:
            # Only do this once so that the same variables are used in this distribution
            # is used multiple times
            self.mangle_dict = pywarpx.my_constants.add_keywords(self.user_defined_kw)

        if (self.Ex_expression is not None or
            self.Ey_expression is not None or
            self.Ez_expression is not None):
            pywarpx.particles.E_ext_particle_init_style = 'parse_e_ext_particle_function'
            for sdir, expression in zip(['x', 'y', 'z'], [self.Ex_expression, self.Ey_expression, self.Ez_expression]):
                expression = pywarpx.my_constants.mangle_expression(expression, self.mangle_dict)
                pywarpx.particles.__setattr__(f'E{sdir}_external_particle_function(x,y,z,t)', expression)

        if (self.Bx_expression is not None or
            self.By_expression is not None or
            self.Bz_expression is not None):
            pywarpx.particles.B_ext_particle_init_style = 'parse_b_ext_particle_function'
            for sdir, expression in zip(['x', 'y', 'z'], [self.Bx_expression, self.By_expression, self.Bz_expression]):
                expression = pywarpx.my_constants.mangle_expression(expression, self.mangle_dict)
                pywarpx.particles.__setattr__(f'B{sdir}_external_particle_function(x,y,z,t)', expression)


class Mirror(picmistandard.PICMI_Mirror):
    def initialize_inputs(self):
        try:
            pywarpx.warpx.num_mirrors
        except AttributeError:
            pywarpx.warpx.num_mirrors = 0
            pywarpx.warpx.mirror_z = []
            pywarpx.warpx.mirror_z_width = []
            pywarpx.warpx.mirror_z_npoints = []

        pywarpx.warpx.num_mirrors += 1
        pywarpx.warpx.mirror_z.append(self.z_front_location)
        pywarpx.warpx.mirror_z_width.append(self.depth)
        pywarpx.warpx.mirror_z_npoints.append(self.number_of_cells)


class FieldIonization(picmistandard.PICMI_FieldIonization):
    """
    WarpX only has ADK ionization model implemented.
    """
    def initialize_inputs(self):
        assert self.model == 'ADK', 'WarpX only has ADK ionization model implemented'
        self.ionized_species.species.do_field_ionization = 1
        self.ionized_species.species.physical_element = self.ionized_species.particle_type
        self.ionized_species.species.ionization_product_species = self.product_species.name
        self.ionized_species.species.ionization_initial_level = self.ionized_species.charge_state
        self.ionized_species.species.charge = 'q_e'


class CoulombCollisions(picmistandard.base._ClassWithInit):
    """
    Custom class to handle setup of binary Coulmb collisions in WarpX. If
    collision initialization is added to picmistandard this can be changed to
    inherit that functionality.

    Parameters
    ----------
    name: string
        Name of instance (used in the inputs file)

    species: list of species instances
        The species involved in the collision. Must be of length 2.

    CoulombLog: float, optional
        Value of the Coulomb log to use in the collision cross section.
        If not supplied, it is calculated from the local conditions.

    ndt: integer, optional
        The collisions will be applied every "ndt" steps. Must be 1 or larger.
    """
    def __init__(self, name, species, CoulombLog=None, ndt=None, **kw):
        self.name = name
        self.species = species
        self.CoulombLog = CoulombLog
        self.ndt = ndt

        self.handle_init(kw)

    def initialize_inputs(self):
        collision = pywarpx.Collisions.newcollision(self.name)
        collision.type = 'pairwisecoulomb'
        collision.species = [species.name for species in self.species]
        collision.CoulombLog = self.CoulombLog
        collision.ndt = self.ndt


class MCCCollisions(picmistandard.base._ClassWithInit):
    """
    Custom class to handle setup of MCC collisions in WarpX. If collision
    initialization is added to picmistandard this can be changed to inherit
    that functionality.

    Parameters
    ----------
    name: string
        Name of instance (used in the inputs file)

    species: species instance
        The species involved in the collision

    background_density: float
        The density of the background

    background_temperature: float
        The temperature of the background

    scattering_processes: dictionary
        The scattering process to use and any needed information

    background_mass: float, optional
        The mass of the background particle. If not supplied, the default depends
        on the type of scattering process.

    ndt: integer, optional
        The collisions will be applied every "ndt" steps. Must be 1 or larger.
    """

    def __init__(self, name, species, background_density,
                 background_temperature, scattering_processes,
                 background_mass=None, ndt=None, **kw):
        self.name = name
        self.species = species
        self.background_density = background_density
        self.background_temperature = background_temperature
        self.background_mass = background_mass
        self.scattering_processes = scattering_processes
        self.ndt = ndt

        self.handle_init(kw)

    def initialize_inputs(self):
        collision = pywarpx.Collisions.newcollision(self.name)
        collision.type = 'background_mcc'
        collision.species = self.species.name
        collision.background_density = self.background_density
        collision.background_temperature = self.background_temperature
        collision.background_mass = self.background_mass
        collision.ndt = self.ndt

        collision.scattering_processes = self.scattering_processes.keys()
        for process, kw in self.scattering_processes.items():
            for key, val in kw.items():
                if key == 'species':
                    val = val.name
                collision.add_new_attr(process+'_'+key, val)


class EmbeddedBoundary(picmistandard.base._ClassWithInit):
    """
    Custom class to handle set up of embedded boundaries specific to WarpX.
    If embedded boundary initialization is added to picmistandard this can be
    changed to inherit that functionality. The geometry can be specified either as
    an implicit function or as an STL file (ASCII or binary). In the latter case the
    geometry specified in the STL file can be scaled, translated and inverted.

    Parameters
    ----------
    implicit_function: string
        Analytic expression describing the embedded boundary

    stl_file: string
        STL file path (string), file contains the embedded boundary geometry

    stl_scale: float
        Factor by which the STL geometry is scaled

    stl_center: vector of floats
        Vector by which the STL geometry is translated (in meters)

    stl_reverse_normal: bool
        If True inverts the orientation of the STL geometry

    potential: string, default=0.
        Analytic expression defining the potential. Can only be specified
        when the solver is electrostatic.

    cover_multiple_cuts: bool, default=None
        Whether to cover cells with multiple cuts.
        (If False, this will raise an error if some cells have multiple cuts)

    Parameters used in the analytic expressions should be given as additional keyword arguments.

    """
    def __init__(self, implicit_function=None, stl_file=None, stl_scale=None, stl_center=None, stl_reverse_normal=False,
                 potential=None, cover_multiple_cuts=None, **kw):

        assert stl_file is None or implicit_function is None, Exception('Only one between implicit_function and '
                                                                            'stl_file can be specified')

        self.implicit_function = implicit_function
        self.stl_file = stl_file

        if stl_file is None:
            assert stl_scale is None, Exception('EB can only be scaled only when using an stl file')
            assert stl_center is None, Exception('EB can only be translated only when using an stl file')
            assert stl_reverse_normal is False, Exception('EB can only be reversed only when using an stl file')

        self.stl_scale = stl_scale
        self.stl_center = stl_center
        self.stl_reverse_normal = stl_reverse_normal

        self.potential = potential

        self.cover_multiple_cuts = cover_multiple_cuts

        # Handle keyword arguments used in expressions
        self.user_defined_kw = {}
        for k in list(kw.keys()):
            if (implicit_function is not None and re.search(r'\b%s\b'%k, implicit_function) or
               (potential is not None and re.search(r'\b%s\b'%k, potential))):
                self.user_defined_kw[k] = kw[k]
                del kw[k]

        self.handle_init(kw)

    def initialize_inputs(self, solver):

        # Add the user defined keywords to my_constants
        # The keywords are mangled if there is a conflicting variable already
        # defined in my_constants with the same name but different value.
        self.mangle_dict = pywarpx.my_constants.add_keywords(self.user_defined_kw)

        if self.implicit_function is not None:
            expression = pywarpx.my_constants.mangle_expression(self.implicit_function, self.mangle_dict)
            pywarpx.warpx.eb_implicit_function = expression

        if self.stl_file is not None:
            pywarpx.eb2.geom_type = "stl"
            pywarpx.eb2.stl_file = self.stl_file
            pywarpx.eb2.stl_scale = self.stl_scale
            pywarpx.eb2.stl_center = self.stl_center
            pywarpx.eb2.stl_reverse_normal = self.stl_reverse_normal

        pywarpx.eb2.cover_multiple_cuts = self.cover_multiple_cuts

        if self.potential is not None:
            assert isinstance(solver, ElectrostaticSolver), Exception('The potential is only supported with the ElectrostaticSolver')
            expression = pywarpx.my_constants.mangle_expression(self.potential, self.mangle_dict)
            pywarpx.warpx.__setattr__('eb_potential(x,y,z,t)', expression)


class PlasmaLens(picmistandard.base._ClassWithInit):
    """
    Custom class to setup a plasma lens lattice.
    The applied fields are dependent only on the transverse position.

    Parameters
    ----------
    period: float
        Periodicity of the lattice (in lab frame, in meters)

    starts: list of floats
        The start of each lens relative to the periodic repeat

    lengths: list of floats
        The length of each lens

    strengths_E=None: list of floats, default = 0.
        The electric field strength of each lens

    strengths_B=None: list of floats, default = 0.
        The magnetic field strength of each lens


    The field that is applied depends on the transverse position of the particle, (x,y)

    - Ex = x*stengths_E

    - Ey = y*stengths_E

    - Bx = +y*stengths_B

    - By = -x*stengths_B

    """
    def __init__(self, period, starts, lengths, strengths_E=None, strengths_B=None, **kw):
        self.period = period
        self.starts = starts
        self.lengths = lengths
        self.strengths_E = strengths_E
        self.strengths_B = strengths_B

        assert (self.strengths_E is not None) or (self.strengths_B is not None),\
               Exception('One of strengths_E or strengths_B must be supplied')

        self.handle_init(kw)

    def initialize_inputs(self):

        pywarpx.particles.E_ext_particle_init_style = 'repeated_plasma_lens'
        pywarpx.particles.B_ext_particle_init_style = 'repeated_plasma_lens'
        pywarpx.particles.repeated_plasma_lens_period = self.period
        pywarpx.particles.repeated_plasma_lens_starts = self.starts
        pywarpx.particles.repeated_plasma_lens_lengths = self.lengths
        pywarpx.particles.repeated_plasma_lens_strengths_E = self.strengths_E
        pywarpx.particles.repeated_plasma_lens_strengths_B = self.strengths_B


class Simulation(picmistandard.PICMI_Simulation):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ----------
    warpx_current_deposition_algo: {'direct', 'esirkepov', and 'vay'}, optional
        Current deposition algorithm. The default depends on conditions.

    warpx_charge_deposition_algo: {'standard'}, optional
        Charge deposition algorithm.

    warpx_field_gathering_algo: {'energy-conserving', 'momentum-conserving'}, optional
        Field gathering algorithm. The default depends on conditions.

    warpx_particle_pusher_algo: {'boris', 'vay', 'higuera'}, default='boris'
        Particle pushing algorithm.

    warpx_use_filter: bool, optional
        Whether to use filtering. The default depends on the conditions.

    warpx_do_multi_J: bool, default=0
        Whether to use the multi-J algorithm, where current deposition and
        field update are performed multiple times within each time step.

    warpx_do_multi_J_n_depositions: integer
        Number of sub-steps to use with the multi-J algorithm, when ``warpx_do_multi_J=1``.
        Note that this input parameter is not optional and must always be set in all
        input files where ``warpx.do_multi_J=1``. No default value is provided automatically.

    warpx_grid_type: {'collocated', 'staggered', 'hybrid'}, default='staggered'
        Whether to use a collocated grid (all fields defined at the cell nodes),
        a staggered grid (fields defined on a Yee grid), or a hybrid grid
        (fields and currents are interpolated back and forth between a staggered grid
        and a collocated grid, must be used with momentum-conserving field gathering algorithm).

    warpx_do_current_centering: bool, optional
        If true, the current is deposited on a nodal grid and then centered
        to a staggered grid (Yee grid), using finite-order interpolation.
        Default: warpx.do_current_centering=0 with collocated or staggered grids,
        warpx.do_current_centering=1 with hybrid grids.

    warpx_field_centering_nox/noy/noz: integer, optional
        The order of interpolation used with staggered or hybrid grids (``warpx_grid_type=staggered``
        or ``warpx_grid_type=hybrid``) and momentum-conserving field gathering
        (``warpx_field_gathering_algo=momentum-conserving``) to interpolate the
        electric and magnetic fields from the cell centers to the cell nodes,
        before gathering the fields from the cell nodes to the particle positions.
        Default: ``warpx_field_centering_no<x,y,z>=2`` with staggered grids,
        ``warpx_field_centering_no<x,y,z>=8`` with hybrid grids (typically necessary
        to ensure stability in boosted-frame simulations of relativistic plasmas and beams).

    warpx_current_centering_nox/noy/noz: integer, optional
        The order of interpolation used with hybrid grids (``warpx_grid_type=hybrid``)
        to interpolate the currents from the cell nodes to the cell centers when
        ``warpx_do_current_centering=1``, before pushing the Maxwell fields on staggered grids.
        Default: ``warpx_current_centering_no<x,y,z>=8`` with hybrid grids (typically necessary
        to ensure stability in boosted-frame simulations of relativistic plasmas and beams).

    warpx_serialize_initial_conditions: bool, default=False
        Controls the random numbers used for initialization.
        This parameter should only be used for testing and continuous integration.

    warpx_random_seed: string or int, optional
        (See documentation)

    warpx_do_dynamic_scheduling: bool, default=True
        Whether to do dynamic scheduling with OpenMP

    warpx_load_balance_intervals: string, default='0'
        The intervals for doing load balancing

    warpx_load_balance_efficiency_ratio_threshold: float, default=1.1
        (See documentation)

    warpx_load_balance_with_sfc: bool, default=0
        (See documentation)

    warpx_load_balance_knapsack_factor: float, default=1.24
        (See documentation)

    warpx_load_balance_costs_update: {'heuristic' or 'timers' or 'gpuclock'}, optional
        (See documentation)

    warpx_costs_heuristic_particles_wt: float, optional
        (See documentation)

    warpx_costs_heuristic_cells_wt: float, optional
        (See documentation)

    warpx_use_fdtd_nci_corr: bool, optional
        Whether to use the NCI correction when using the FDTD solver

    warpx_amr_check_input: bool, optional
        Whether AMReX should perform checks on the input
        (primarily related to the max grid size and blocking factors)

    warpx_amr_restart: string, optional
        The name of the restart to use

    warpx_amrex_the_arena_is_managed: bool, optional
        Whether to use managed memory in the AMReX Arena

    warpx_amrex_the_arena_init_size: long int, optional
        The amount of memory in bytes to allocate in the Arena.

    warpx_zmax_plasma_to_compute_max_step: float, optional
        Sets the simulation run time based on the maximum z value

    warpx_compute_max_step_from_btd: bool, default=0
        If specified, automatically calculates the number of iterations
        required in the boosted frame for all back-transformed diagnostics
        to be completed.

    warpx_collisions: collision instance, optional
        The collision instance specifying the particle collisions

    warpx_embedded_boundary: embedded boundary instance, optional

    warpx_break_signals: list of strings
        Signals on which to break

    warpx_checkpoint_signals: list of strings
        Signals on which to write out a checkpoint
    """

    # Set the C++ WarpX interface (see _libwarpx.LibWarpX) as an extension to
    # Simulation objects. In the future, LibWarpX objects may actually be owned
    # by Simulation objects to permit multiple WarpX runs simultaneously.
    extension = pywarpx.libwarpx

    def init(self, kw):

        self.current_deposition_algo = kw.pop('warpx_current_deposition_algo', None)
        self.charge_deposition_algo = kw.pop('warpx_charge_deposition_algo', None)
        self.field_gathering_algo = kw.pop('warpx_field_gathering_algo', None)
        self.particle_pusher_algo = kw.pop('warpx_particle_pusher_algo', None)
        self.use_filter = kw.pop('warpx_use_filter', None)
        self.do_multi_J = kw.pop('warpx_do_multi_J', None)
        self.do_multi_J_n_depositions = kw.pop('warpx_do_multi_J_n_depositions', None)
        self.grid_type = kw.pop('warpx_grid_type', None)
        self.do_current_centering = kw.pop('warpx_do_current_centering', None)
        self.field_centering_order = kw.pop('warpx_field_centering_order', None)
        self.current_centering_order = kw.pop('warpx_current_centering_order', None)
        self.serialize_initial_conditions = kw.pop('warpx_serialize_initial_conditions', None)
        self.random_seed = kw.pop('warpx_random_seed', None)
        self.do_dynamic_scheduling = kw.pop('warpx_do_dynamic_scheduling', None)
        self.load_balance_intervals = kw.pop('warpx_load_balance_intervals', None)
        self.load_balance_efficiency_ratio_threshold = kw.pop('warpx_load_balance_efficiency_ratio_threshold', None)
        self.load_balance_with_sfc = kw.pop('warpx_load_balance_with_sfc', None)
        self.load_balance_knapsack_factor = kw.pop('warpx_load_balance_knapsack_factor', None)
        self.load_balance_costs_update = kw.pop('warpx_load_balance_costs_update', None)
        self.costs_heuristic_particles_wt = kw.pop('warpx_costs_heuristic_particles_wt', None)
        self.costs_heuristic_cells_wt = kw.pop('warpx_costs_heuristic_cells_wt', None)
        self.use_fdtd_nci_corr = kw.pop('warpx_use_fdtd_nci_corr', None)
        self.amr_check_input = kw.pop('warpx_amr_check_input', None)
        self.amr_restart = kw.pop('warpx_amr_restart', None)
        self.amrex_the_arena_is_managed = kw.pop('warpx_amrex_the_arena_is_managed', None)
        self.amrex_the_arena_init_size = kw.pop('warpx_amrex_the_arena_init_size', None)
        self.zmax_plasma_to_compute_max_step = kw.pop('warpx_zmax_plasma_to_compute_max_step', None)
        self.compute_max_step_from_btd = kw.pop('warpx_compute_max_step_from_btd', None)

        self.collisions = kw.pop('warpx_collisions', None)
        self.embedded_boundary = kw.pop('warpx_embedded_boundary', None)

        self.break_signals = kw.pop('warpx_break_signals', None)
        self.checkpoint_signals = kw.pop('warpx_checkpoint_signals', None)

        self.inputs_initialized = False
        self.warpx_initialized = False

    def initialize_inputs(self):
        if self.inputs_initialized:
            return

        self.inputs_initialized = True

        pywarpx.warpx.verbose = self.verbose
        if self.time_step_size is not None:
            pywarpx.warpx.const_dt = self.time_step_size

        if self.gamma_boost is not None:
            pywarpx.warpx.gamma_boost = self.gamma_boost
            pywarpx.warpx.boost_direction = 'z'

        pywarpx.warpx.zmax_plasma_to_compute_max_step = self.zmax_plasma_to_compute_max_step
        pywarpx.warpx.compute_max_step_from_btd = self.compute_max_step_from_btd

        pywarpx.algo.current_deposition = self.current_deposition_algo
        pywarpx.algo.charge_deposition = self.charge_deposition_algo
        pywarpx.algo.field_gathering = self.field_gathering_algo
        pywarpx.algo.particle_pusher = self.particle_pusher_algo
        pywarpx.algo.load_balance_intervals = self.load_balance_intervals
        pywarpx.algo.load_balance_efficiency_ratio_threshold = self.load_balance_efficiency_ratio_threshold
        pywarpx.algo.load_balance_with_sfc = self.load_balance_with_sfc
        pywarpx.algo.load_balance_knapsack_factor = self.load_balance_knapsack_factor
        pywarpx.algo.load_balance_costs_update = self.load_balance_costs_update
        pywarpx.algo.costs_heuristic_particles_wt = self.costs_heuristic_particles_wt
        pywarpx.algo.costs_heuristic_cells_wt = self.costs_heuristic_cells_wt

        pywarpx.warpx.grid_type = self.grid_type
        pywarpx.warpx.do_current_centering = self.do_current_centering
        pywarpx.warpx.use_filter = self.use_filter
        pywarpx.warpx.do_multi_J = self.do_multi_J
        pywarpx.warpx.do_multi_J_n_depositions = self.do_multi_J_n_depositions
        pywarpx.warpx.serialize_initial_conditions = self.serialize_initial_conditions
        pywarpx.warpx.random_seed = self.random_seed

        pywarpx.warpx.do_dynamic_scheduling = self.do_dynamic_scheduling

        pywarpx.particles.use_fdtd_nci_corr = self.use_fdtd_nci_corr

        pywarpx.amr.check_input = self.amr_check_input

        pywarpx.warpx.break_signals = self.break_signals
        pywarpx.warpx.checkpoint_signals = self.checkpoint_signals

        particle_shape = self.particle_shape
        for s in self.species:
            if s.particle_shape is not None:
                assert particle_shape is None or particle_shape == s.particle_shape, Exception('WarpX only supports one particle shape for all species')
                # --- If this was set for any species, use that value.
                particle_shape = s.particle_shape

        if particle_shape is not None and (len(self.species) > 0 or len(self.lasers) > 0):
            if isinstance(particle_shape, str):
                interpolation_order = {'NGP':0, 'linear':1, 'quadratic':2, 'cubic':3}[particle_shape]
            else:
                interpolation_order = particle_shape
            pywarpx.algo.particle_shape = interpolation_order

        self.solver.initialize_inputs()

        # Initialize warpx.field_centering_no<x,y,z> and warpx.current_centering_no<x,y,z>
        # if set by the user in the input (need to access grid info from solver attribute)
        # warpx.field_centering_no<x,y,z>
        if self.field_centering_order is not None:
            pywarpx.warpx.field_centering_nox = self.field_centering_order[0]
            if self.solver.grid.number_of_dimensions == 3:
                pywarpx.warpx.field_centering_noy = self.field_centering_order[1]
            pywarpx.warpx.field_centering_noz = self.field_centering_order[-1]
        # warpx.current_centering_no<x,y,z>
        if self.current_centering_order is not None:
            pywarpx.warpx.current_centering_nox = self.current_centering_order[0]
            if self.solver.grid.number_of_dimensions == 3:
                pywarpx.warpx.current_centering_noy = self.current_centering_order[1]
            pywarpx.warpx.current_centering_noz = self.current_centering_order[-1]

        for i in range(len(self.species)):
            self.species[i].initialize_inputs(self.layouts[i],
                                              self.initialize_self_fields[i],
                                              self.injection_plane_positions[i],
                                              self.injection_plane_normal_vectors[i])

        for interaction in self.interactions:
            assert(isinstance(interaction, FieldIonization))
            interaction.initialize_inputs()

        if self.collisions is not None:
            pywarpx.collisions.collision_names = []
            for collision in self.collisions:
                pywarpx.collisions.collision_names.append(collision.name)
                collision.initialize_inputs()

        if self.embedded_boundary is not None:
            self.embedded_boundary.initialize_inputs(self.solver)

        for i in range(len(self.lasers)):
            self.lasers[i].initialize_inputs()
            self.laser_injection_methods[i].initialize_inputs(self.lasers[i])

        for applied_field in self.applied_fields:
            applied_field.initialize_inputs()

        for diagnostic in self.diagnostics:
            diagnostic.initialize_inputs()

        if self.amr_restart:
            pywarpx.amr.restart = self.amr_restart

        if self.amrex_the_arena_is_managed is not None:
            pywarpx.amrex.the_arena_is_managed = self.amrex_the_arena_is_managed

        if self.amrex_the_arena_init_size is not None:
            pywarpx.amrex.the_arena_init_size = self.amrex_the_arena_init_size

    def initialize_warpx(self, mpi_comm=None):
        if self.warpx_initialized:
            return

        self.warpx_initialized = True
        pywarpx.warpx.init(mpi_comm, max_step=self.max_steps, stop_time=self.max_time)

    def write_input_file(self, file_name='inputs'):
        self.initialize_inputs()
        pywarpx.warpx.write_inputs(file_name, max_step=self.max_steps, stop_time=self.max_time)

    def step(self, nsteps=None, mpi_comm=None):
        self.initialize_inputs()
        self.initialize_warpx(mpi_comm)
        if nsteps is None:
            if self.max_steps is not None:
                nsteps = self.max_steps
            else:
                nsteps = -1
        pywarpx.warpx.evolve(nsteps)

    def finalize(self):
        if self.warpx_initialized:
            self.warpx_initialized = False
            pywarpx.warpx.finalize()


# ----------------------------
# Simulation frame diagnostics
# ----------------------------

class WarpXDiagnosticBase(object):
    """
    Base class for all WarpX diagnostic containing functionality shared by
    all WarpX diagnostic installations.
    """
    def add_diagnostic(self):
        # reduced diagnostics go in a different bucket than regular diagnostics
        if isinstance(self, ReducedDiagnostic):
            bucket = pywarpx.reduced_diagnostics
            name_template = 'reduced_diag'
        else:
            bucket = pywarpx.diagnostics
            name_template = 'diag'

        name = getattr(self, 'name', None)
        if name is None:
            diagnostics_number = (len(bucket._diagnostics_dict) + 1)
            self.name = f'{name_template}{diagnostics_number}'

        try:
            self.diagnostic = bucket._diagnostics_dict[self.name]
        except KeyError:
            self.diagnostic = pywarpx.Diagnostics.Diagnostic(
                self.name, _species_dict={}
            )
            bucket._diagnostics_dict[self.name] = self.diagnostic

    def set_write_dir(self):
        if self.write_dir is not None or self.file_prefix is not None:
            write_dir = (self.write_dir or 'diags')
            file_prefix = (self.file_prefix or self.name)
            self.diagnostic.file_prefix = os.path.join(write_dir, file_prefix)


class FieldDiagnostic(picmistandard.PICMI_FieldDiagnostic, WarpXDiagnosticBase):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ----------
    warpx_plot_raw_fields: bool, optional
        Flag whether to dump the raw fields

    warpx_plot_raw_fields_guards: bool, optional
        Flag whether the raw fields should include the guard cells

    warpx_write_species: bool, optional
        Flag whether to output particle data with the diagnostic

    warpx_format: {plotfile, checkpoint, openpmd, ascent, sensei}, optional
        Diagnostic file format

    warpx_openpmd_backend: {bp, h5, json}, optional
        Openpmd backend file format

    warpx_file_prefix: string, optional
        Prefix on the diagnostic file name

    warpx_file_min_digits: integer, optional
        Minimum number of digits for the time step number in the file name

    warpx_dump_rz_modes: bool, optional
        Flag whether to dump the data for all RZ modes
    """
    def init(self, kw):

        self.plot_raw_fields = kw.pop('warpx_plot_raw_fields', None)
        self.plot_raw_fields_guards = kw.pop('warpx_plot_raw_fields_guards', None)
        self.plot_finepatch = kw.pop('warpx_plot_finepatch', None)
        self.plot_crsepatch = kw.pop('warpx_plot_crsepatch', None)
        self.write_species = kw.pop('warpx_write_species', None)
        self.format = kw.pop('warpx_format', 'plotfile')
        self.openpmd_backend = kw.pop('warpx_openpmd_backend', None)
        self.file_prefix = kw.pop('warpx_file_prefix', None)
        self.file_min_digits = kw.pop('warpx_file_min_digits', None)
        self.dump_rz_modes = kw.pop('warpx_dump_rz_modes', None)

    def initialize_inputs(self):

        self.add_diagnostic()

        self.diagnostic.diag_type = 'Full'
        self.diagnostic.format = self.format
        self.diagnostic.openpmd_backend = self.openpmd_backend
        self.diagnostic.file_min_digits = self.file_min_digits
        self.diagnostic.dump_rz_modes = self.dump_rz_modes
        self.diagnostic.intervals = self.period
        self.diagnostic.diag_lo = self.lower_bound
        self.diagnostic.diag_hi = self.upper_bound
        if self.number_of_cells is not None:
            self.diagnostic.coarsening_ratio = (np.array(self.grid.number_of_cells)/np.array(self.number_of_cells)).astype(int)

        # --- Use a set to ensure that fields don't get repeated.
        fields_to_plot = set()

        if pywarpx.geometry.dims == 'RZ':
            E_fields_list = ['Er', 'Et', 'Ez']
            B_fields_list = ['Br', 'Bt', 'Bz']
            J_fields_list = ['Jr', 'Jt', 'Jz']
            A_fields_list = ['Ar', 'At', 'Az']
        else:
            E_fields_list = ['Ex', 'Ey', 'Ez']
            B_fields_list = ['Bx', 'By', 'Bz']
            J_fields_list = ['Jx', 'Jy', 'Jz']
            A_fields_list = ['Ax', 'Ay', 'Az']
        if self.data_list is not None:
            for dataname in self.data_list:
                if dataname == 'E':
                    for field_name in E_fields_list:
                        fields_to_plot.add(field_name)
                elif dataname == 'B':
                    for field_name in B_fields_list:
                        fields_to_plot.add(field_name)
                elif dataname == 'J':
                    for field_name in J_fields_list:
                        fields_to_plot.add(field_name.lower())
                elif dataname == 'A':
                    for field_name in A_fields_list:
                        fields_to_plot.add(field_name)
                elif dataname in E_fields_list:
                    fields_to_plot.add(dataname)
                elif dataname in B_fields_list:
                    fields_to_plot.add(dataname)
                elif dataname in A_fields_list:
                    fields_to_plot.add(dataname)
                elif dataname in ['rho', 'phi', 'F', 'proc_number', 'part_per_cell']:
                    fields_to_plot.add(dataname)
                elif dataname in J_fields_list:
                    fields_to_plot.add(dataname.lower())
                elif dataname.startswith('rho_'):
                    # Adds rho_species diagnostic
                    fields_to_plot.add(dataname)
                elif dataname == 'dive':
                    fields_to_plot.add('divE')
                elif dataname == 'divb':
                    fields_to_plot.add('divB')
                elif dataname == 'raw_fields':
                    self.plot_raw_fields = 1
                elif dataname == 'raw_fields_guards':
                    self.plot_raw_fields_guards = 1
                elif dataname == 'finepatch':
                    self.plot_finepatch = 1
                elif dataname == 'crsepatch':
                    self.plot_crsepatch = 1
                elif dataname == 'none':
                    fields_to_plot = set(('none',))

            # --- Convert the set to a sorted list so that the order
            # --- is the same on all processors.
            fields_to_plot = list(fields_to_plot)
            fields_to_plot.sort()
            self.diagnostic.fields_to_plot = fields_to_plot

        self.diagnostic.plot_raw_fields = self.plot_raw_fields
        self.diagnostic.plot_raw_fields_guards = self.plot_raw_fields_guards
        self.diagnostic.plot_finepatch = self.plot_finepatch
        self.diagnostic.plot_crsepatch = self.plot_crsepatch
        self.diagnostic.write_species = self.write_species

        self.set_write_dir()


ElectrostaticFieldDiagnostic = FieldDiagnostic


class Checkpoint(picmistandard.base._ClassWithInit, WarpXDiagnosticBase):
    """
    Sets up checkpointing of the simulation, allowing for later restarts

    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ----------
    warpx_file_prefix: string
        The prefix to the checkpoint directory names

    warpx_file_min_digits: integer
        Minimum number of digits for the time step number in the checkpoint
        directory name.
    """

    def __init__(self, period = 1, write_dir = None, name = None, **kw):

        self.period = period
        self.write_dir = write_dir
        self.file_prefix = kw.pop('warpx_file_prefix', None)
        self.file_min_digits = kw.pop('warpx_file_min_digits', None)
        self.name = name

        if self.name is None:
            self.name = 'chkpoint'

        self.handle_init(kw)

    def initialize_inputs(self):

        self.add_diagnostic()

        self.diagnostic.intervals = self.period
        self.diagnostic.diag_type = 'Full'
        self.diagnostic.format = 'checkpoint'
        self.diagnostic.file_min_digits = self.file_min_digits

        self.set_write_dir()


class ParticleDiagnostic(picmistandard.PICMI_ParticleDiagnostic, WarpXDiagnosticBase):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ----------
    warpx_format: {plotfile, checkpoint, openpmd, ascent, sensei}, optional
        Diagnostic file format

    warpx_openpmd_backend: {bp, h5, json}, optional
        Openpmd backend file format

    warpx_file_prefix: string, optional
        Prefix on the diagnostic file name

    warpx_file_min_digits: integer, optional
        Minimum number of digits for the time step number in the file name

    warpx_random_fraction: float, optional
        Random fraction of particles to include in the diagnostic

    warpx_uniform_stride: integer, optional
        Stride to down select to the particles to include in the diagnostic

    warpx_plot_filter_function: string, optional
        Analytic expression to down select the particles to in the diagnostic
    """
    def init(self, kw):

        self.format = kw.pop('warpx_format', 'plotfile')
        self.openpmd_backend = kw.pop('warpx_openpmd_backend', None)
        self.file_prefix = kw.pop('warpx_file_prefix', None)
        self.file_min_digits = kw.pop('warpx_file_min_digits', None)
        self.random_fraction = kw.pop('warpx_random_fraction', None)
        self.uniform_stride = kw.pop('warpx_uniform_stride', None)
        self.plot_filter_function = kw.pop('warpx_plot_filter_function', None)

        self.user_defined_kw = {}
        if self.plot_filter_function is not None:
            # This allows variables to be used in the plot_filter_function, but
            # in order not to break other codes, the variables must begin with "warpx_"
            for k in list(kw.keys()):
                if k.startswith('warpx_') and re.search(r'\b%s\b'%k, self.plot_filter_function):
                    self.user_defined_kw[k] = kw[k]
                    del kw[k]

        self.mangle_dict = None

    def initialize_inputs(self):

        self.add_diagnostic()

        self.diagnostic.diag_type = 'Full'
        self.diagnostic.format = self.format
        self.diagnostic.openpmd_backend = self.openpmd_backend
        self.diagnostic.file_min_digits = self.file_min_digits
        self.diagnostic.intervals = self.period

        self.set_write_dir()

        # --- Use a set to ensure that fields don't get repeated.
        variables = set()

        if self.data_list is not None:
            for dataname in self.data_list:
                if dataname == 'position':
                    # --- The positions are alway written out anyway
                    pass
                elif dataname == 'momentum':
                    variables.add('ux')
                    variables.add('uy')
                    variables.add('uz')
                elif dataname == 'weighting':
                    variables.add('w')
                elif dataname == 'fields':
                    variables.add('Ex')
                    variables.add('Ey')
                    variables.add('Ez')
                    variables.add('Bx')
                    variables.add('By')
                    variables.add('Bz')
                elif dataname in ['ux', 'uy', 'uz', 'Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz']:
                    variables.add(dataname)

            # --- Convert the set to a sorted list so that the order
            # --- is the same on all processors.
            variables = list(variables)
            variables.sort()

        if np.iterable(self.species):
            species_list = self.species
        else:
            species_list = [self.species]

        if self.mangle_dict is None:
            # Only do this once so that the same variables are used in this distribution
            # is used multiple times
            self.mangle_dict = pywarpx.my_constants.add_keywords(self.user_defined_kw)

        for specie in species_list:
            diag = pywarpx.Bucket.Bucket(self.name + '.' + specie.name,
                                         variables = variables,
                                         random_fraction = self.random_fraction,
                                         uniform_stride = self.uniform_stride)
            expression = pywarpx.my_constants.mangle_expression(self.plot_filter_function, self.mangle_dict)
            diag.__setattr__('plot_filter_function(t,x,y,z,ux,uy,uz)', expression)
            self.diagnostic._species_dict[specie.name] = diag

# ----------------------------
# Lab frame diagnostics
# ----------------------------


class LabFrameFieldDiagnostic(picmistandard.PICMI_LabFrameFieldDiagnostic,
                              WarpXDiagnosticBase):
    """
    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html>`_ for more information.

    Parameters
    ----------
    warpx_format: string, optional
        Passed to <diagnostic name>.format

    warpx_openpmd_backend: string, optional
        Passed to <diagnostic name>.openpmd_backend

    warpx_file_prefix: string, optional
        Passed to <diagnostic name>.file_prefix

    warpx_file_min_digits: integer, optional
        Passed to <diagnostic name>.file_min_digits

    warpx_buffer_size: integer, optional
        Passed to <diagnostic name>.buffer_size

    warpx_lower_bound: vector of floats, optional
        Passed to <diagnostic name>.lower_bound

    warpx_upper_bound: vector of floats, optional
        Passed to <diagnostic name>.upper_bound
    """
    def init(self, kw):
        # The user is using the new BTD
        self.format = kw.pop('warpx_format', None)
        self.openpmd_backend = kw.pop('warpx_openpmd_backend', None)
        self.file_prefix = kw.pop('warpx_file_prefix', None)
        self.file_min_digits = kw.pop('warpx_file_min_digits', None)
        self.buffer_size = kw.pop('warpx_buffer_size', None)
        self.lower_bound = kw.pop('warpx_lower_bound', None)
        self.upper_bound = kw.pop('warpx_upper_bound', None)

    def initialize_inputs(self):

        self.add_diagnostic()

        self.diagnostic.diag_type = 'BackTransformed'
        self.diagnostic.format = self.format
        self.diagnostic.openpmd_backend = self.openpmd_backend
        self.diagnostic.file_min_digits = self.file_min_digits
        self.diagnostic.diag_lo = self.lower_bound
        self.diagnostic.diag_hi = self.upper_bound

        self.diagnostic.do_back_transformed_fields = 1
        self.diagnostic.num_snapshots_lab = self.num_snapshots
        self.diagnostic.dt_snapshots_lab = self.dt_snapshots
        self.diagnostic.buffer_size = self.buffer_size

        # --- Use a set to ensure that fields don't get repeated.
        fields_to_plot = set()

        if pywarpx.geometry.dims == 'RZ':
            E_fields_list = ['Er', 'Et', 'Ez']
            B_fields_list = ['Br', 'Bt', 'Bz']
            J_fields_list = ['Jr', 'Jt', 'Jz']
        else:
            E_fields_list = ['Ex', 'Ey', 'Ez']
            B_fields_list = ['Bx', 'By', 'Bz']
            J_fields_list = ['Jx', 'Jy', 'Jz']
        if self.data_list is not None:
            for dataname in self.data_list:
                if dataname == 'E':
                    for field_name in E_fields_list:
                        fields_to_plot.add(field_name)
                elif dataname == 'B':
                    for field_name in B_fields_list:
                        fields_to_plot.add(field_name)
                elif dataname == 'J':
                    for field_name in J_fields_list:
                        fields_to_plot.add(field_name.lower())
                elif dataname in E_fields_list:
                    fields_to_plot.add(dataname)
                elif dataname in B_fields_list:
                    fields_to_plot.add(dataname)
                elif dataname in J_fields_list:
                    fields_to_plot.add(dataname.lower())
                elif dataname.startswith('rho_'):
                    # Adds rho_species diagnostic
                    fields_to_plot.add(dataname)

            # --- Convert the set to a sorted list so that the order
            # --- is the same on all processors.
            fields_to_plot = list(fields_to_plot)
            fields_to_plot.sort()
            self.diagnostic.fields_to_plot = fields_to_plot

        self.set_write_dir()

class ReducedDiagnostic(picmistandard.base._ClassWithInit, WarpXDiagnosticBase):
    """
    Sets up a reduced diagnostic in the simulation.

    See `Input Parameters <https://warpx.readthedocs.io/en/latest/usage/parameters.html#reduced-diagnostics>`_
    for more information.

    Parameters
    ----------
    diag_type: string
        The type of reduced diagnostic. See the link above for all the different
        types of reduced diagnostics available.

    name: string
        The name of this diagnostic which will also be the name of the data
        file written to disk.

    period: integer
        The simulation step interval at which to output this diagnostic.

    path: string
        The file path in which the diagnostic file should be written.

    extension: string
        The file extension used for the diagnostic output.

    separator: string
        The separator between row values in the output file.

    species: species instance
        The name of the species for which to calculate the diagnostic, required for
        diagnostic types 'BeamRelevant', 'ParticleHistogram', and 'ParticleExtrema'

    bin_number: integer
        For diagnostic type 'ParticleHistogram', the number of bins used for the histogram

    bin_max: float
        For diagnostic type 'ParticleHistogram', the maximum value of the bins

    bin_min: float
        For diagnostic type 'ParticleHistogram', the minimum value of the bins

    normalization: {'unity_particle_weight', 'max_to_unity', 'area_to_unity'}, optional
        For diagnostic type 'ParticleHistogram', normalization method of the histogram.

    histogram_function: string
        For diagnostic type 'ParticleHistogram', the function evaluated to produce the histogram data

    filter_function: string, optional
        For diagnostic type 'ParticleHistogram', the function to filter whether particles are included in the histogram

    reduced_function: string
        For diagnostic type 'FieldReduction', the function of the fields to evaluate

    weighting_function: string, optional
        For diagnostic type 'ChargeOnEB', the function to weight contributions to the total charge

    reduction_type: {'Maximum', 'Minimum', or 'Integral'}
        For diagnostic type 'FieldReduction', the type of reduction

    probe_geometry: {'Point', 'Line', 'Plane'}, defaut='Point'
        For diagnostic type 'FieldProbe', the geometry of the probe

    integrate: bool, default=false
        For diagnostic type 'FieldProbe', whether the field is integrated

    do_moving_window_FP: bool, default=False
        For diagnostic type 'FieldProbe', whether the moving window is followed

    x_probe, y_probe, z_probe: floats
        For diagnostic type 'FieldProbe', a probe location. For 'Point', the location of the point. For 'Line', the start of the
        line. For 'Plane', the center of the square detector.

    interp_order: integer
        For diagnostic type 'FieldProbe', the interpolation order for 'Line' and 'Plane'

    resolution: integer
        For diagnostic type 'FieldProbe', the number of points along the 'Line' or along each edge of the square 'Plane'

    x1_probe, y1_probe, z1_probe: floats
        For diagnostic type 'FieldProbe', the end point for 'Line'

    detector_radius: float
        For diagnostic type 'FieldProbe', the detector "radius" (half edge length) of the 'Plane'

    target_normal_x, target_normal_y, target_normal_z: floats
        For diagnostic type 'FieldProbe', the normal vector to the 'Plane'. Only applicable in 3D

    target_up_x, target_up_y, target_up_z: floats
        For diagnostic type 'FieldProbe', the vector specifying up in the 'Plane'
    """

    def __init__(self, diag_type, name=None, period=1, path=None,
                 extension=None, separator=None, **kw):

        self.name = name
        self.type = diag_type
        self.intervals = period
        self.path = path
        self.extension = extension
        self.separator = separator

        self.user_defined_kw = {}

        # Now we need to handle all the specific inputs required for the
        # different reduced diagnostic types.

        # The simple diagnostics do not require any additional arguments
        self._simple_reduced_diagnostics = [
            'ParticleEnergy', 'ParticleMomentum', 'FieldEnergy',
            'FieldMomentum', 'FieldMaximum', 'RhoMaximum', 'ParticleNumber',
            'LoadBalanceCosts', 'LoadBalanceEfficiency'
        ]
        # The species diagnostics require a species to be provided
        self._species_reduced_diagnostics = [
            'BeamRelevant', 'ParticleHistogram', 'ParticleExtrema'
        ]

        if self.type in self._simple_reduced_diagnostics:
            pass
        elif self.type in self._species_reduced_diagnostics:
            species = kw.pop('species')
            self.species = species.name
            if self.type == 'ParticleHistogram':
                kw = self._handle_particle_histogram(**kw)
        elif self.type == "FieldProbe":
            kw = self._handle_field_probe(**kw)
        elif self.type == "FieldReduction":
            kw = self._handle_field_reduction(**kw)
        elif self.type == "ChargeOnEB":
            kw = self._handle_charge_on_eb(**kw)
        else:
            raise RuntimeError(
                f"{self.type} reduced diagnostic is not yet supported "
                "in pywarpx."
            )

        self.handle_init(kw)

    def _handle_field_probe(self, **kw):
        """Utility function to grab required inputs for a field probe from kw"""
        self.probe_geometry = kw.pop("probe_geometry")
        self.x_probe = kw.pop("x_probe", None)
        self.y_probe = kw.pop("y_probe", None)
        self.z_probe = kw.pop("z_probe")

        self.interp_order = kw.pop("interp_order", None)
        self.integrate = kw.pop("integrate", None)
        self.do_moving_window_FP = kw.pop("do_moving_window_FP", None)

        if self.probe_geometry.lower() != 'point':
            self.resolution = kw.pop("resolution")

        if self.probe_geometry.lower() == 'line':
            self.x1_probe = kw.pop("x1_probe", None)
            self.y1_probe = kw.pop("y1_probe", None)
            self.z1_probe = kw.pop("z1_probe")

        if self.probe_geometry.lower() == 'plane':
            self.detector_radius = kw.pop("detector_radius")

            self.target_normal_x = kw.pop("target_normal_x", None)
            self.target_normal_y = kw.pop("target_normal_y", None)
            self.target_normal_z = kw.pop("target_normal_z", None)

            self.target_up_x = kw.pop("target_up_x", None)
            self.target_up_y = kw.pop("target_up_y", None)
            self.target_up_z = kw.pop("target_up_z", None)

        return kw

    def _handle_particle_histogram(self, **kw):
        self.bin_number = kw.pop("bin_number")
        self.bin_max = kw.pop("bin_max")
        self.bin_min = kw.pop("bin_min")
        self.normalization = kw.pop("normalization", None)
        if self.normalization not in [None, "unity_particle_weight", "max_to_unity", "area_to_unity"]:
            raise AttributeError(
               "The ParticleHistogram normalization must be one of 'unity_particle_weight', 'max_to_unity', or 'area_to_unity'")

        histogram_function = kw.pop("histogram_function")
        filter_function = kw.pop("filter_function", None)

        self.__setattr__("histogram_function(t,x,y,z,ux,uy,uz)", histogram_function)
        self.__setattr__("filter_function(t,x,y,z,ux,uy,uz)", filter_function)

        # Check the reduced function expressions for constants
        for k in list(kw.keys()):
            if (re.search(r'\b%s\b'%k, histogram_function) or
                (filter_function is not None and re.search(r'\b%s\b'%k, filter_function))):
                self.user_defined_kw[k] = kw[k]
                del kw[k]

        return kw

    def _handle_field_reduction(self, **kw):
        self.reduction_type = kw.pop("reduction_type")
        reduced_function = kw.pop("reduced_function")

        self.__setattr__("reduced_function(x,y,z,Ex,Ey,Ez,Bx,By,Bz)", reduced_function)

        # Check the reduced function expression for constants
        for k in list(kw.keys()):
            if re.search(r'\b%s\b'%k, reduced_function):
                self.user_defined_kw[k] = kw[k]
                del kw[k]

        return kw

    def _handle_charge_on_eb(self, **kw):
        weighting_function = kw.pop("weighting_function", None)

        self.__setattr__("weighting_function(x,y,z)", weighting_function)

        # Check the reduced function expression for constants
        for k in list(kw.keys()):
            if re.search(r'\b%s\b'%k, weighting_function):
                self.user_defined_kw[k] = kw[k]
                del kw[k]

        return kw

    def initialize_inputs(self):

        self.add_diagnostic()

        self.mangle_dict = pywarpx.my_constants.add_keywords(self.user_defined_kw)

        for key, value in self.__dict__.items():
            if not key.startswith('_') and key not in ['name', 'diagnostic']:
                if key.endswith(")"):
                    # Analytic expressions require processing to deal with constants
                    expression = pywarpx.my_constants.mangle_expression(value, self.mangle_dict)
                    self.diagnostic.__setattr__(key, expression)
                else:
                    self.diagnostic.__setattr__(key, value)
