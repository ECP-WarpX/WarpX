/* Copyright 2019-2022 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_INDEX_HANDLING_H
#define ABLASTR_INDEX_HANDLING_H

#include <cstdint>


namespace ablastr::particles {

    /** A helper function to derive a globally unique particle ID
     *
     * @param[in] id  AMReX particle ID (on local cpu/rank), AoS .id
     * @param[in] cpu AMReX particle CPU (rank) at creation of the particle, AoS .cpu
     * @return global particle ID that is unique and permanent in the whole simulation
     */
    constexpr uint64_t
    localIDtoGlobal (int const id, int const cpu)
    {
        static_assert(sizeof(int) * 2u <= sizeof(uint64_t),
                      "int size might cause collisions in global IDs");
        // implementation:
        //   - we cast both 32-bit (or smaller) ints to a 64bit unsigned int
        //   - this will leave half of the "upper range" bits in the 64bit unsigned int zeroed out
        //     because the corresponding (extended) value range was not part of the value range in
        //     the int representation
        //   - we bit-shift the cpu into the upper half of zero bits in the 64 bit unsigned int
        //     (imagine this step as "adding a per-cpu/rank offset to the local integers")
        //   - then we add this offset
        //     note: the add is expressed as bitwise OR (|) since this saves us from writing
        //           brackets for operator precedence between + and <<
        return uint64_t(id) | uint64_t(cpu) << 32u;
    }

} // namespace ablastr::particles

#endif // ABLASTR_INDEX_HANDLING_H
