/* Copyright 2024 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FIELDS_H_
#define WARPX_FIELDS_H_

namespace warpx::fields
{
    enum struct FieldType : int
    {
        None,
        Efield_aux,
        Bfield_aux,
        Efield_fp,
        Bfield_fp,
        Efield_fp_external,
        Bfield_fp_external,
        current_fp,
        current_fp_nodal,
        rho_fp,
        F_fp,
        G_fp,
        phi_fp,
        vector_potential_fp,
        Efield_cp,
        Bfield_cp,
        current_cp,
        rho_cp,
        F_cp,
        G_cp,
        edge_lengths,
        face_areas,
        Efield_avg_fp,
        Bfield_avg_fp,
        Efield_avg_cp,
        Bfield_avg_cp
    };

    inline bool
    isFieldArray (const FieldType field_type)
    {
        if (field_type == FieldType::Efield_aux || field_type == FieldType::Bfield_aux ||
            field_type == FieldType::Efield_fp  || field_type == FieldType::Bfield_fp ||
            field_type == FieldType::current_fp || field_type == FieldType::current_fp_nodal ||
            field_type == FieldType::vector_potential_fp ||
            field_type == FieldType::Efield_cp  || field_type == FieldType::Bfield_cp ||
            field_type == FieldType::current_cp ||
            field_type == FieldType::Efield_avg_fp || field_type == FieldType::Bfield_avg_fp ||
            field_type == FieldType::Efield_avg_cp || field_type == FieldType::Bfield_avg_cp) { return true; }
        else { return false; }
    }
}

#endif //WARPX_FIELDS_H_
