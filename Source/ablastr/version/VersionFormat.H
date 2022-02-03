/* Copyright 2022 Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_VERSION_FORMAT
#define ABLASTR_VERSION_FORMAT


namespace ablastr::version {

    /** Collect software dependencies of ABLASTR/WarpX
     *
     * amrex::Print() the versions of dependent software in ABLASTR/WarpX
     * for output and reproducibility.
     */
    void
    printSoftwareVersions();

} // namespace ablastr::version

#endif // ABLASTR_VERSION_FORMAT
