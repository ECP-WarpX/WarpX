/* Copyright 2024 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FIELDS_H_
#define WARPX_FIELDS_H_

#include <ablastr/fields/MultiFabRegister.H>

#include <AMReX_Enum.H>

#include <algorithm>
#include <iterator>


namespace warpx::fields
{
    AMREX_ENUM(FieldType,
        None,
        Efield_aux,
        Bfield_aux,
        Efield_fp,
        Bfield_fp,
        Efield_fp_external,
        Bfield_fp_external,
        current_fp,
        current_fp_nodal,
        current_fp_vay,
        current_buf,
        current_store,
        rho_buf,
        rho_fp,
        F_fp,
        G_fp,
        phi_fp,
        vector_potential_fp,
        vector_potential_fp_nodal,
        vector_potential_grad_buf_e_stag,
        vector_potential_grad_buf_b_stag,
        hybrid_electron_pressure_fp,
        hybrid_rho_fp_temp,
        hybrid_current_fp_temp,
        hybrid_current_fp_ampere,
        hybrid_current_fp_external,
        Efield_cp,
        Bfield_cp,
        current_cp,
        rho_cp,
        F_cp,
        G_cp,
        Efield_cax,
        Bfield_cax,
        E_external_particle_field,
        B_external_particle_field,
        distance_to_eb,
        edge_lengths,
        face_areas,
        area_mod,
        pml_E_fp,
        pml_B_fp,
        pml_j_fp,
        pml_F_fp,
        pml_G_fp,
        pml_E_cp,
        pml_B_cp,
        pml_j_cp,
        pml_F_cp,
        pml_G_cp,
        pml_edge_lengths,
        Efield_avg_fp,
        Bfield_avg_fp,
        Efield_avg_cp,
        Bfield_avg_cp,
        Bold,
        ECTRhofield,
        Venl
    );

    /** these are vector fields */
    constexpr FieldType ArrayFieldTypes[] = {
        FieldType::Efield_aux,
        FieldType::Bfield_aux,
        FieldType::Efield_fp,
        FieldType::Bfield_fp,
        FieldType::current_fp,
        FieldType::current_fp_nodal,
        FieldType::current_fp_vay,
        FieldType::current_buf,
        FieldType::current_store,
        FieldType::vector_potential_fp,
        FieldType::vector_potential_fp_nodal,
        FieldType::vector_potential_grad_buf_e_stag,
        FieldType::vector_potential_grad_buf_b_stag,
        FieldType::hybrid_current_fp_temp,
        FieldType::hybrid_current_fp_ampere,
        FieldType::hybrid_current_fp_external,
        FieldType::Efield_cp,
        FieldType::Bfield_cp,
        FieldType::current_cp,
        FieldType::Efield_cax,
        FieldType::Bfield_cax,
        FieldType::E_external_particle_field,
        FieldType::B_external_particle_field,
        FieldType::pml_E_fp,
        FieldType::pml_B_fp,
        FieldType::pml_j_fp,
        FieldType::pml_E_cp,
        FieldType::pml_B_cp,
        FieldType::pml_j_cp,
        FieldType::Efield_avg_fp,
        FieldType::Bfield_avg_fp,
        FieldType::Efield_avg_cp,
        FieldType::Bfield_avg_cp,
        FieldType::Bold,
        FieldType::ECTRhofield,
        FieldType::Venl
    };

    /** Returns true if a FieldType represents a vector field */
    inline bool
    isFieldArray (const FieldType field_type)
    {
        return std::any_of( std::begin(ArrayFieldTypes), std::end(ArrayFieldTypes),
            [field_type](const FieldType& f) { return f == field_type; });
    }

}

#endif //WARPX_FIELDS_H_
