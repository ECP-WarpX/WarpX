/* Copyright 2024 Axel Huebl
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef ABLASTR_UTILS_ENUMS_H_
#define ABLASTR_UTILS_ENUMS_H_

namespace ablastr::utils::enums
{
    /** Type of grids used in a simulation:
     *
     * Collocated at the same location (AMReX: all "NODAL"), staggered (Yee-style), or hybrid.
     */
    enum struct GridType {
        Collocated = 0,
        Staggered = 1,
        Hybrid = 2
    };

    /** Mesh-refinement patch
     *
     * The fine or coarse patch (in terms of spatial resolution) on the same MR level.
     * https://warpx.readthedocs.io/en/latest/theory/amr.html
     */
    enum struct PatchType {
        fine,
        coarse
    };

} // namespace ablastr::utils::enums

#endif // ABLASTR_UTILS_ENUMS_H_
