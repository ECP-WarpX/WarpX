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
        vector_potential_fp,
        phi_fp
    };

    constexpr FieldType ArrayFieldTypes[] = {
            FieldType::Efield_fp,
            FieldType::Bfield_fp,
            FieldType::vector_potential_fp };

    inline bool
    isFieldArray (const FieldType field_type)
    {
        return std::any_of( std::begin(ArrayFieldTypes), std::end(ArrayFieldTypes),
            [field_type](const FieldType& f) { return f == field_type; });
    }

}

#endif //WARPX_FIELDS_H_
