/* Copyright 2021 Lorenzo Giacomel, Neil Zaim, Yinjian Zhao
 *                Tiberius Rheaume, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_

#include "ReducedDiags.H"
#include "FieldProbeParticleContainer.H"

#include <AMReX.H>
#include <AMReX_Vector.H>

#include <unordered_map>
#include <string>
#include <vector>

/**
 * This enumeration is used for assigning structural geometry levels (point vs line vs plane)
 */
enum struct DetectorGeometry
{
    Point = 0,
    Line,
    Plane
};

/**
 *  This class mainly contains a function that computes the value of each component
 * of the EM field at a given point
 */
class FieldProbe : public ReducedDiags
{
public:

    /**
     * constructor
     * @param[in] rd_name reduced diags names
     */
    FieldProbe (std::string rd_name);

    /**
     * This function assins test/data particles to constructed environemnt
     */
    void InitData () override final;

    /** Redistribute parallel data structures during load balance
     */
    void LoadBalance () override final;

    /**
     * This function computes the value of Ex, Ey, Ez, Bx, By, Bz and at a given point
     *
     * @param[in] step current time step
     */
    void ComputeDiags (int step) override final;

    /*
     * Define constants used throughout FieldProbe
     */

    //! noutputs is 10 (x, y, z, Ex, Ey, Ez, Bx, By, Bz, S)
    static constexpr int noutputs = FieldProbePIdx::nattribs + 3;

private:
    amrex::Real x_probe, y_probe, z_probe;
    amrex::Real x1_probe, y1_probe, z1_probe;
    amrex::Real target_normal_x, target_normal_y, target_normal_z;
    amrex::Real target_up_x, target_up_y, target_up_z;
    amrex::Real detector_radius;

    //! counts number of particles for all MPI ranks
    long m_valid_particles {0};

    //! determines geometry of detector point distribution
    DetectorGeometry m_probe_geometry = DetectorGeometry::Point;

    //! determines number of particles places for non-point geometries
    int m_resolution = 0;

    //! Empty vector for to which data is pushed
    amrex::Vector<amrex::Real> m_data;

    //! Empty array to be used by IOProcessor node to store and output data
    amrex::Vector<amrex::Real> m_data_out;

    //! this is the particle container in which probe particles are stored
    FieldProbeParticleContainer m_probe;

    //! if true, integrate values over time instead of probing instantaneous values
    bool m_field_probe_integrate = false;

    //! particle shape used for field gather
    int interp_order = 1;

    //! Judges whether to gather raw fields or interpolated data
    bool raw_fields = false;

    /**
     * Built-in function in ReducedDiags to write out test data
     */
    virtual void WriteToFile (int step) const override;

    /** Check if the probe is in the simulation domain boundary
     */
    bool ProbeInDomain () const;
};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_FIELDPROBE_H_
