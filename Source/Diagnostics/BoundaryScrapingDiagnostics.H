/* Copyright 2022 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_BOUNDARYSCRAPINGDIAGNOSTICS_H_
#define WARPX_BOUNDARYSCRAPINGDIAGNOSTICS_H_

#include "Diagnostics.H"

#include <string>

class
BoundaryScrapingDiagnostics final : public Diagnostics
{
public:

    BoundaryScrapingDiagnostics (int i, std::string name);

private:
    /** Read relevant parameters for BoundaryScraping */
    void ReadParameters ();
    /** \brief Flush particles to file. */
    void Flush (int i_buffer) override;
    /** whether to write output files at this time step
     *  The data is flushed when the buffer is full and/or
     *  when the simulation ends or when forced.
     * \param[in] step current time step
     * \param[in] i_buffer index of a back-transformed snapshot
     * \param[in] force_flush if true, return true for any step
     *  The return bool is used to determine when to write buffer data out
     */
    bool DoDump (int step, int i_buffer, bool force_flush=false) override;
    /** TODO
     * \param[in] step current time step, return true for any step >= 0
     * \param[in] force_flush if true, return true for any step
     */
    bool DoComputeAndPack (int step, bool force_flush=false) override;
    /** Initialize functors that store pointers to the fields requested by the user.
     *  Additionally, the cell-center functors that stores pointers to all fields,
     *  namely, Ex, Ey, Ez, Bx, By, Bz, jx, jy, jz, and rho is also initialized.
     * \param[in] lev level on which the vector of unique_ptrs to field functors
     *                is initialized.
     */
    void InitializeFieldFunctors (int lev) override;
    /** This function allocates and initializes particle buffers for all the snapshots.
     * This is currently an empty function:
     * The particle containers required for this must be added to populate this function.
     */
    void InitializeParticleBuffer () override;
    /** Initialize m_mf_output vectors and data required to construct the buffers
     * \param[in] i_buffer index of buffer for which the output MultiFab is defined.
     * \param[in] lev level on which the output MultiFab is defined
     */
     void InitializeBufferData (int i_buffer, int lev) override;

};
#endif // WARPX_BOUNDARYSCRAPINGDIAGNOSTICS_H_
