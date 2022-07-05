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
    'open':'pml', 'dirichlet':'pec', 'periodic':'periodic', 'damped':'damped', 'none':'none', None:'none'
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

        self.boost_adjust_transverse_positions = kw.pop('warpx_boost_adjust_transverse_positions', None)

        # For the relativistic electrostatic solver
        self.self_fields_required_precision = kw.pop('warpx_self_fields_required_precision', None)
        self.self_fields_absolute_tolerance = kw.pop('warpx_self_fields_absolute_tolerance', None)
        self.self_fields_max_iters = kw.pop('warpx_self_fields_max_iters', None)
        self.self_fields_verbosity = kw.pop('warpx_self_fields_verbosity', None)
        self.save_previous_position = kw.pop('warpx_save_previous_position', None)
        self.do_not_deposit = kw.pop('warpx_do_not_deposit', None)

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
                                             do_not_deposit = self.do_not_deposit)

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

        for interaction in self.interactions:
            assert interaction[0] == 'ionization'
            assert interaction[1] == 'ADK', 'WarpX only has ADK ionization model implemented'
            self.species.do_field_ionization=1
            self.species.physical_element=self.particle_type
            self.species.ionization_product_species = interaction[2].name
            self.species.ionization_initial_level = self.charge_state
            self.species.charge = 'q_e'

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


class UniformDistribution(picmistandard.PICMI_UniformDistribution):
    def initialize_inputs(self, species_number, layout, species, density_scale):

        if isinstance(layout, GriddedLayout):
            # --- Note that the grid attribute of GriddedLayout is ignored
            species.injection_style = "nuniformpercell"
            species.num_particles_per_cell_each_dim = layout.n_macroparticle_per_cell
        elif isinstance(layout, PseudoRandomLayout):
            assert (layout.n_macroparticles_per_cell is not None), Exception('WarpX only supports n_macroparticles_per_cell for the PseudoRandomLayout with UniformDistribution')
            species.injection_style = "nrandompercell"
            species.num_particles_per_cell = layout.n_macroparticles_per_cell
        else:
            raise Exception('WarpX does not support the specified layout for UniformDistribution')

        species.xmin = self.lower_bound[0]
        species.xmax = self.upper_bound[0]
        species.ymin = self.lower_bound[1]
        species.ymax = self.upper_bound[1]
        species.zmin = self.lower_bound[2]
        species.zmax = self.upper_bound[2]

        # --- Only constant density is supported at this time
        species.profile = "constant"
        species.density = self.density
        if density_scale is not None:
            species.density *= density_scale

        # --- Note that WarpX takes gamma*beta as input
        if np.any(np.not_equal(self.rms_velocity, 0.)):
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

        if self.fill_in:
            species.do_continuous_injection = 1


class AnalyticDistribution(picmistandard.PICMI_AnalyticDistribution):
    def init(self, kw):
        self.mangle_dict = None

    def initialize_inputs(self, species_number, layout, species, density_scale):

        if isinstance(layout, GriddedLayout):
            # --- Note that the grid attribute of GriddedLayout is ignored
            species.injection_style = "nuniformpercell"
            species.num_particles_per_cell_each_dim = layout.n_macroparticle_per_cell
        elif isinstance(layout, PseudoRandomLayout):
            assert (layout.n_macroparticles_per_cell is not None), Exception('WarpX only supports n_macroparticles_per_cell for the PseudoRandomLayout with UniformDistribution')
            species.injection_style = "nrandompercell"
            species.num_particles_per_cell = layout.n_macroparticles_per_cell
        else:
            raise Exception('WarpX does not support the specified layout for UniformDistribution')

        species.xmin = self.lower_bound[0]
        species.xmax = self.upper_bound[0]
        species.ymin = self.lower_bound[1]
        species.ymax = self.upper_bound[1]
        species.zmin = self.lower_bound[2]
        species.zmax = self.upper_bound[2]

        if self.mangle_dict is None:
            # Only do this once so that the same variables are used in this distribution
            # is used multiple times
            self.mangle_dict = pywarpx.my_constants.add_keywords(self.user_defined_kw)
        expression = pywarpx.my_constants.mangle_expression(self.density_expression, self.mangle_dict)

        species.profile = "parse_density_function"
        if density_scale is None:
            species.__setattr__('density_function(x,y,z)', expression)
        else:
            species.__setattr__('density_function(x,y,z)', "{}*({})".format(density_scale, expression))

        # --- Note that WarpX takes gamma*beta as input
        if np.any(np.not_equal(self.momentum_expressions, None)):
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

        if self.fill_in:
            species.do_continuous_injection = 1

    def setup_parse_momentum_functions(self, species):
        for sdir, idir in zip(['x', 'y', 'z'], [0, 1, 2]):
            if self.momentum_expressions[idir] is not None:
                expression = pywarpx.my_constants.mangle_expression(self.momentum_expressions[idir], self.mangle_dict)
            else:
                expression = f'{self.directed_velocity[idir]}'
            species.__setattr__(f'momentum_function_u{sdir}(x,y,z)', f'({expression})/{constants.c}')

class ParticleListDistribution(picmistandard.PICMI_ParticleListDistribution):
    def init(self, kw):
        pass

    def initialize_inputs(self, species_number, layout, species, density_scale):

        species.injection_style = "multipleparticles"
        species.multiple_particles_pos_x = self.x
        species.multiple_particles_pos_y = self.y
        species.multiple_particles_pos_z = self.z
        species.multiple_particles_vel_x = np.array(self.ux)/constants.c
        species.multiple_particles_vel_y = np.array(self.uy)/constants.c
        species.multiple_particles_vel_z = np.array(self.uz)/constants.c
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
    """This assumes that WarpX was compiled with USE_RZ = TRUE
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
    def init(self, kw):
        assert self.method is None or self.method in ['Yee', 'CKC', 'PSATD', 'ECT'], Exception("Only 'Yee', 'CKC', 'PSATD', and 'ECT' are supported")

        self.pml_ncell = kw.pop('warpx_pml_ncell', None)

        if self.method == 'PSATD':
            self.psatd_periodic_single_box_fft = kw.pop('warpx_periodic_single_box_fft', None)
            self.psatd_current_correction = kw.pop('warpx_current_correction', None)
            self.psatd_update_with_rho = kw.pop('warpx_psatd_update_with_rho', None)
            self.psatd_do_time_averaging = kw.pop('warpx_psatd_do_time_averaging', None)

        self.do_pml_in_domain = kw.pop('warpx_do_pml_in_domain', None)
        self.pml_has_particles = kw.pop('warpx_pml_has_particles', None)
        self.do_pml_j_damping = kw.pop('warpx_do_pml_j_damping', None)

    def initialize_inputs(self):

        self.grid.initialize_inputs()

        pywarpx.warpx.pml_ncell = self.pml_ncell
        pywarpx.warpx.do_nodal = self.l_nodal

        if self.method == 'PSATD':
            pywarpx.psatd.periodic_single_box_fft = self.psatd_periodic_single_box_fft
            pywarpx.psatd.current_correction = self.psatd_current_correction
            pywarpx.psatd.update_with_rho = self.psatd_update_with_rho
            pywarpx.psatd.do_time_averaging = self.psatd_do_time_averaging

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
    def init(self, kw):
        self.relativistic = kw.pop('warpx_relativistic', False)
        self.absolute_tolerance = kw.pop('warpx_absolute_tolerance', None)
        self.self_fields_verbosity = kw.pop('warpx_self_fields_verbosity', None)

    def initialize_inputs(self):

        self.grid.initialize_inputs()

        if self.relativistic:
            pywarpx.warpx.do_electrostatic = 'relativistic'
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
        laser.laser.direction = self.normal_vector  # The plane normal direction
        if isinstance(laser, GaussianLaser):
            laser.laser.profile_focal_distance = laser.focal_position[2] - self.position[2]  # Focal distance from the antenna (in meters)
            laser.laser.profile_t_peak = (self.position[2] - laser.centroid_position[2])/constants.c  # The time at which the laser reaches its peak (in seconds)


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


class CoulombCollisions(picmistandard.base._ClassWithInit):
    """Custom class to handle setup of binary Coulmb collisions in WarpX. If
    collision initialization is added to picmistandard this can be changed to
    inherit that functionality."""

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
    """Custom class to handle setup of MCC collisions in WarpX. If collision
    initialization is added to picmistandard this can be changed to inherit
    that functionality."""

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
    - implicit_function: Analytic expression describing the embedded boundary
    - stl_file: STL file path (string), file contains the embedded boundary geometry
    - stl_scale: factor by which the STL geometry is scaled (pure number)
    - stl_center: vector by which the STL geometry is translated (in meters)
    - stl_reverse_normal: if True inverts the orientation of the STL geometry
    - potential: Analytic expression defining the potential. Can only be specified
                 when the solver is electrostatic. Optional, defaults to 0.
     Parameters used in the expressions should be given as additional keyword arguments.
    """
    def __init__(self, implicit_function=None, stl_file=None, stl_scale=None, stl_center=None, stl_reverse_normal=False,
                 potential=None, **kw):

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

        if self.potential is not None:
            assert isinstance(solver, ElectrostaticSolver), Exception('The potential is only supported with the ElectrostaticSolver')
            expression = pywarpx.my_constants.mangle_expression(self.potential, self.mangle_dict)
            pywarpx.warpx.__setattr__('eb_potential(x,y,z,t)', expression)


class PlasmaLens(picmistandard.base._ClassWithInit):
    """
    Custom class to setup a plasma lens lattice.
    The applied fields are dependent on the transverse position
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
        self.serialize_initial_conditions = kw.pop('warpx_serialize_initial_conditions', None)
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

        pywarpx.warpx.use_filter = self.use_filter
        pywarpx.warpx.serialize_initial_conditions = self.serialize_initial_conditions

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

        for i in range(len(self.species)):
            self.species[i].initialize_inputs(self.layouts[i],
                                              self.initialize_self_fields[i],
                                              self.injection_plane_positions[i],
                                              self.injection_plane_normal_vectors[i])

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

    def initialize_warpx(self, mpi_comm=None):
        if self.warpx_initialized:
            return

        self.warpx_initialized = True
        pywarpx.warpx.init(mpi_comm)

    def write_input_file(self, file_name='inputs'):
        self.initialize_inputs()
        kw = {}
        if self.max_steps is not None:
            kw['max_step'] = self.max_steps
        if self.max_time is not None:
            kw['stop_time'] = self.max_time
        pywarpx.warpx.write_inputs(file_name, **kw)

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


class FieldDiagnostic(picmistandard.PICMI_FieldDiagnostic):
    def init(self, kw):

        self.plot_raw_fields = kw.pop('warpx_plot_raw_fields', None)
        self.plot_raw_fields_guards = kw.pop('warpx_plot_raw_fields_guards', None)
        self.plot_finepatch = kw.pop('warpx_plot_finepatch', None)
        self.plot_crsepatch = kw.pop('warpx_plot_crsepatch', None)
        self.format = kw.pop('warpx_format', 'plotfile')
        self.openpmd_backend = kw.pop('warpx_openpmd_backend', None)
        self.file_prefix = kw.pop('warpx_file_prefix', None)
        self.file_min_digits = kw.pop('warpx_file_min_digits', None)
        self.dump_rz_modes = kw.pop('warpx_dump_rz_modes', None)

    def initialize_inputs(self):

        name = getattr(self, 'name', None)
        if name is None:
            diagnostics_number = len(pywarpx.diagnostics._diagnostics_dict) + 1
            self.name = 'diag{}'.format(diagnostics_number)

        try:
            self.diagnostic = pywarpx.diagnostics._diagnostics_dict[self.name]
        except KeyError:
            self.diagnostic = pywarpx.Diagnostics.Diagnostic(self.name, _species_dict={})
            pywarpx.diagnostics._diagnostics_dict[self.name] = self.diagnostic

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

        if self.data_list is not None:
            for dataname in self.data_list:
                if dataname == 'E':
                    fields_to_plot.add('Ex')
                    fields_to_plot.add('Ey')
                    fields_to_plot.add('Ez')
                elif dataname == 'B':
                    fields_to_plot.add('Bx')
                    fields_to_plot.add('By')
                    fields_to_plot.add('Bz')
                elif dataname == 'J':
                    fields_to_plot.add('jx')
                    fields_to_plot.add('jy')
                    fields_to_plot.add('jz')
                elif dataname in ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz', 'rho', 'phi', 'F', 'proc_number', 'part_per_cell']:
                    fields_to_plot.add(dataname)
                elif dataname in ['Jx', 'Jy', 'Jz']:
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

        if self.write_dir is not None or self.file_prefix is not None:
            write_dir = (self.write_dir or 'diags')
            file_prefix = (self.file_prefix or self.name)
            self.diagnostic.file_prefix = os.path.join(write_dir, file_prefix)


ElectrostaticFieldDiagnostic = FieldDiagnostic


class Checkpoint(picmistandard.base._ClassWithInit):

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

        try:
            self.diagnostic = pywarpx.diagnostics._diagnostics_dict[self.name]
        except KeyError:
            self.diagnostic = pywarpx.Diagnostics.Diagnostic(self.name, _species_dict={})
            pywarpx.diagnostics._diagnostics_dict[self.name] = self.diagnostic

        self.diagnostic.intervals = self.period
        self.diagnostic.diag_type = 'Full'
        self.diagnostic.format = 'checkpoint'
        self.diagnostic.file_min_digits = self.file_min_digits

        if self.write_dir is not None or self.file_prefix is not None:
            write_dir = (self.write_dir or 'diags')
            file_prefix = (self.file_prefix or self.name)
            self.diagnostic.file_prefix = os.path.join(write_dir, file_prefix)

class ParticleDiagnostic(picmistandard.PICMI_ParticleDiagnostic):
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

        name = getattr(self, 'name', None)
        if name is None:
            diagnostics_number = len(pywarpx.diagnostics._diagnostics_dict) + 1
            self.name = 'diag{}'.format(diagnostics_number)

        try:
            self.diagnostic = pywarpx.diagnostics._diagnostics_dict[self.name]
        except KeyError:
            self.diagnostic = pywarpx.Diagnostics.Diagnostic(self.name, _species_dict={})
            pywarpx.diagnostics._diagnostics_dict[self.name] = self.diagnostic

        self.diagnostic.diag_type = 'Full'
        self.diagnostic.format = self.format
        self.diagnostic.openpmd_backend = self.openpmd_backend
        self.diagnostic.file_min_digits = self.file_min_digits
        self.diagnostic.intervals = self.period

        if self.write_dir is not None or self.file_prefix is not None:
            write_dir = (self.write_dir or 'diags')
            file_prefix = (self.file_prefix or self.name)
            self.diagnostic.file_prefix = os.path.join(write_dir, file_prefix)

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


class LabFrameFieldDiagnostic(picmistandard.PICMI_LabFrameFieldDiagnostic):
    """
    Warp specific arguments:
      - warpx_new_BTD: Use the new BTD diagnostics
      - warpx_format: Passed to <diagnostic name>.format
      - warpx_openpmd_backend: Passed to <diagnostic name>.openpmd_backend
      - warpx_file_prefix: Passed to <diagnostic name>.file_prefix
      - warpx_file_min_digits: Passed to <diagnostic name>.file_min_digits
      - warpx_buffer_size: Passed to <diagnostic name>.buffer_size
      - warpx_lower_bound: Passed to <diagnostic name>.lower_bound
      - warpx_upper_bound: Passed to <diagnostic name>.upper_bound
    """
    __doc__ = picmistandard.PICMI_LabFrameFieldDiagnostic.__doc__ + __doc__
    def init(self, kw):
        self.use_new_BTD = kw.pop('warpx_new_BTD', False)
        if self.use_new_BTD:
            # The user is using the new BTD
            self.format = kw.pop('warpx_format', None)
            self.openpmd_backend = kw.pop('warpx_openpmd_backend', None)
            self.file_prefix = kw.pop('warpx_file_prefix', None)
            self.file_min_digits = kw.pop('warpx_file_min_digits', None)
            self.buffer_size = kw.pop('warpx_buffer_size', None)
            self.lower_bound = kw.pop('warpx_lower_bound', None)
            self.upper_bound = kw.pop('warpx_upper_bound', None)

    def initialize_inputs(self):
        if self.use_new_BTD:
            self.initialize_inputs_new()
        else:
            self.initialize_inputs_old()

    def initialize_inputs_old(self):

        pywarpx.warpx.check_consistency('num_snapshots_lab', self.num_snapshots, 'The number of snapshots must be the same in all lab frame diagnostics')
        pywarpx.warpx.check_consistency('dt_snapshots_lab', self.dt_snapshots, 'The time between snapshots must be the same in all lab frame diagnostics')
        pywarpx.warpx.check_consistency('lab_data_directory', self.write_dir, 'The write directory must be the same in all lab frame diagnostics')

        pywarpx.warpx.do_back_transformed_diagnostics = 1
        pywarpx.warpx.num_snapshots_lab = self.num_snapshots
        pywarpx.warpx.dt_snapshots_lab = self.dt_snapshots
        pywarpx.warpx.do_back_transformed_fields = 1
        pywarpx.warpx.lab_data_directory = self.write_dir

    def initialize_inputs_new(self):

        name = getattr(self, 'name', None)
        if name is None:
            diagnostics_number = len(pywarpx.diagnostics._diagnostics_dict) + 1
            self.name = 'diag{}'.format(diagnostics_number)

        try:
            self.diagnostic = pywarpx.diagnostics._diagnostics_dict[self.name]
        except KeyError:
            self.diagnostic = pywarpx.Diagnostics.Diagnostic(self.name, _species_dict={})
            pywarpx.diagnostics._diagnostics_dict[self.name] = self.diagnostic

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

        if self.data_list is not None:
            for dataname in self.data_list:
                if dataname == 'E':
                    fields_to_plot.add('Ex')
                    fields_to_plot.add('Ey')
                    fields_to_plot.add('Ez')
                elif dataname == 'B':
                    fields_to_plot.add('Bx')
                    fields_to_plot.add('By')
                    fields_to_plot.add('Bz')
                elif dataname == 'J':
                    fields_to_plot.add('jx')
                    fields_to_plot.add('jy')
                    fields_to_plot.add('jz')
                elif dataname in ['Ex', 'Ey', 'Ez', 'Bx', 'By', 'Bz', 'rho']:
                    fields_to_plot.add(dataname)
                elif dataname in ['Jx', 'Jy', 'Jz']:
                    fields_to_plot.add(dataname.lower())
                elif dataname.startswith('rho_'):
                    # Adds rho_species diagnostic
                    fields_to_plot.add(dataname)

            # --- Convert the set to a sorted list so that the order
            # --- is the same on all processors.
            fields_to_plot = list(fields_to_plot)
            fields_to_plot.sort()
            self.diagnostic.fields_to_plot = fields_to_plot

        if self.write_dir is not None or self.file_prefix is not None:
            write_dir = (self.write_dir or 'diags')
            file_prefix = (self.file_prefix or self.name)
            self.diagnostic.file_prefix = os.path.join(write_dir, file_prefix)


class LabFrameParticleDiagnostic(picmistandard.PICMI_LabFrameParticleDiagnostic):
    def initialize_inputs(self):

        pywarpx.warpx.check_consistency('num_snapshots_lab', self.num_snapshots, 'The number of snapshots must be the same in all lab frame diagnostics')
        pywarpx.warpx.check_consistency('dt_snapshots_lab', self.dt_snapshots, 'The time between snapshots must be the same in all lab frame diagnostics')
        pywarpx.warpx.check_consistency('lab_data_directory', self.write_dir, 'The write directory must be the same in all lab frame diagnostics')

        pywarpx.warpx.do_back_transformed_diagnostics = 1

        if isinstance(self.species, Species):
            self.species.do_back_transformed_diagnostics = 1
        else:
            try:
                for specie in self.species:
                    specie.do_back_transformed_diagnostics = 1
            except TypeError:
                pass

        pywarpx.warpx.num_snapshots_lab = self.num_snapshots
        pywarpx.warpx.dt_snapshots_lab = self.dt_snapshots
        pywarpx.warpx.lab_data_directory = self.write_dir
