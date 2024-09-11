/* Copyright 2024 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_FIELDS_H_
#define WARPX_FIELDS_H_

#include <algorithm>
#include <iterator>

namespace warpx::fields
{
    enum struct FieldType : int
    {
        None,
        Efield_fp,
        Bfield_fp,
        current_fp,
        current_fp_nodal,
        phi_fp,
        vector_potential_fp,
        Efield_cp,
        Bfield_cp,
        current_cp,
        edge_lengths,
        face_areas,
        Efield_avg_fp,
        Bfield_avg_fp,
        Efield_avg_cp,
        Bfield_avg_cp
    };

    constexpr FieldType ArrayFieldTypes[] = {
            FieldType::Efield_fp, FieldType::Bfield_fp,
            FieldType::current_fp, FieldType::current_fp_nodal, FieldType::vector_potential_fp,
            FieldType::Efield_cp, FieldType::Bfield_cp, FieldType::current_cp,
            FieldType::Efield_avg_fp, FieldType::Bfield_avg_fp, FieldType::Efield_avg_cp, FieldType::Bfield_avg_cp};

    inline bool
    isFieldArray (const FieldType field_type)
    {
        return std::any_of( std::begin(ArrayFieldTypes), std::end(ArrayFieldTypes),
            [field_type](const FieldType& f) { return f == field_type; });
    }

}

#endif //WARPX_FIELDS_H_
