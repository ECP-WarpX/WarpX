/* Copyright 2019-2020 Andrew Myers, Burlen Loring, Luca Fedeli
 * Maxence Thevenet, Remi Lehe, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "WarpXAlgorithmSelection.H"
#include "WarpXConst.H"
#include "WarpXUtil.H"

#include <AMReX_ParmParse.H>

#include <cmath>
#include <fstream>
#include <set>
#include <string>


using namespace amrex;

void ParseGeometryInput()
{
    ParmParse pp_geometry("geometry");

    Vector<Real> prob_lo(AMREX_SPACEDIM);
    Vector<Real> prob_hi(AMREX_SPACEDIM);

    getArrWithParser(pp_geometry, "prob_lo", prob_lo, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_lo.size() == AMREX_SPACEDIM);
    getArrWithParser(pp_geometry, "prob_hi", prob_hi, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

    pp_geometry.addarr("prob_lo", prob_lo);
    pp_geometry.addarr("prob_hi", prob_hi);
}

void ReadBoostedFrameParameters(Real& gamma_boost, Real& beta_boost,
                                Vector<int>& boost_direction)
{
    ParmParse pp_warpx("warpx");
    queryWithParser(pp_warpx, "gamma_boost", gamma_boost);
    if( gamma_boost > 1. ) {
        beta_boost = std::sqrt(1.-1./pow(gamma_boost,2));
        std::string s;
        pp_warpx.get("boost_direction", s);
        if (s == "x" || s == "X") {
            boost_direction[0] = 1;
        }
#if (AMREX_SPACEDIM == 3)
        else if (s == "y" || s == "Y") {
            boost_direction[1] = 1;
        }
#endif
        else if (s == "z" || s == "Z") {
            boost_direction[2] = 1;
        }
        else {
            const std::string msg = "Unknown boost_dir: "+s;
            Abort(msg.c_str());
        }

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( s == "z" || s == "Z" ,
                                          "The boost must be in the z direction.");
    }
}

void ConvertLabParamsToBoost()
{
    Real gamma_boost = 1., beta_boost = 0.;
    int max_level = 0;
    Vector<int> boost_direction {0,0,0};

    ReadBoostedFrameParameters(gamma_boost, beta_boost, boost_direction);

    if (gamma_boost <= 1.) return;

    Vector<Real> prob_lo(AMREX_SPACEDIM);
    Vector<Real> prob_hi(AMREX_SPACEDIM);
    Vector<Real> fine_tag_lo(AMREX_SPACEDIM);
    Vector<Real> fine_tag_hi(AMREX_SPACEDIM);
    Vector<Real> slice_lo(AMREX_SPACEDIM);
    Vector<Real> slice_hi(AMREX_SPACEDIM);

    ParmParse pp_geometry("geometry");
    ParmParse pp_warpx("warpx");
    ParmParse pp_amr("amr");
    ParmParse pp_slice("slice");

    getArrWithParser(pp_geometry, "prob_lo", prob_lo, 0, AMREX_SPACEDIM);
    getArrWithParser(pp_geometry, "prob_hi", prob_hi, 0, AMREX_SPACEDIM);

    queryArrWithParser(pp_slice, "dom_lo", slice_lo, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(slice_lo.size() == AMREX_SPACEDIM);
    queryArrWithParser(pp_slice, "dom_hi", slice_hi, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(slice_hi.size() == AMREX_SPACEDIM);


    pp_amr.query("max_level", max_level);
    if (max_level > 0){
      getArrWithParser(pp_warpx, "fine_tag_lo", fine_tag_lo);
      getArrWithParser(pp_warpx, "fine_tag_hi", fine_tag_hi);
    }


#if (AMREX_SPACEDIM == 3)
    Vector<int> dim_map {0, 1, 2};
#else
    Vector<int> dim_map {0, 2};
#endif

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (boost_direction[dim_map[idim]]) {
            amrex::Real convert_factor;
            // Assume that the window travels with speed +c
            convert_factor = 1./( gamma_boost * ( 1 - beta_boost ) );
            prob_lo[idim] *= convert_factor;
            prob_hi[idim] *= convert_factor;
            if (max_level > 0){
              fine_tag_lo[idim] *= convert_factor;
              fine_tag_hi[idim] *= convert_factor;
            }
            slice_lo[idim] *= convert_factor;
            slice_hi[idim] *= convert_factor;
            break;
        }
    }

    pp_geometry.addarr("prob_lo", prob_lo);
    pp_geometry.addarr("prob_hi", prob_hi);
    if (max_level > 0){
      pp_warpx.addarr("fine_tag_lo", fine_tag_lo);
      pp_warpx.addarr("fine_tag_hi", fine_tag_hi);
    }

    pp_slice.addarr("dom_lo",slice_lo);
    pp_slice.addarr("dom_hi",slice_hi);

}

/* \brief Function that sets the value of MultiFab MF to zero for z between
 * zmin and zmax.
 */
void NullifyMF(amrex::MultiFab& mf, int lev, amrex::Real zmin, amrex::Real zmax){
    WARPX_PROFILE("WarpXUtil::NullifyMF()");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for(amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi){
        const amrex::Box& bx = mfi.tilebox();
        // Get box lower and upper physical z bound, and dz
#if (AMREX_SPACEDIM == 3)
            amrex::Array<amrex::Real,3> galilean_shift = { 0., 0., 0., };
#elif (AMREX_SPACEDIM == 2)
            amrex::Array<amrex::Real,3> galilean_shift = { 0., std::numeric_limits<Real>::quiet_NaN(),  0., } ;
#endif
        const amrex::Real zmin_box = WarpX::LowerCorner(bx, galilean_shift, lev)[2];
        const amrex::Real zmax_box = WarpX::UpperCorner(bx, lev)[2];
        amrex::Real dz  = WarpX::CellSize(lev)[2];
        // Get box lower index in the z direction
#if (AMREX_SPACEDIM==3)
        const int lo_ind = bx.loVect()[2];
#else
        const int lo_ind = bx.loVect()[1];
#endif
        // Check if box intersect with [zmin, zmax]
        if ( (zmax>zmin_box && zmin<=zmax_box) ){
            Array4<Real> arr = mf[mfi].array();
            // Set field to 0 between zmin and zmax
            ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept{
#if (AMREX_SPACEDIM==3)
                    const Real z_gridpoint = zmin_box+(k-lo_ind)*dz;
#else
                    const Real z_gridpoint = zmin_box+(j-lo_ind)*dz;
#endif
                    if ( (z_gridpoint >= zmin) && (z_gridpoint < zmax) ) {
                        arr(i,j,k) = 0.;
                    }
                }
            );
        }
    }
}

namespace WarpXUtilIO{
    bool WriteBinaryDataOnFile(std::string filename, const amrex::Vector<char>& data)
    {
        std::ofstream of{filename, std::ios::binary};
        of.write(data.data(), data.size());
        of.close();
        return  of.good();
    }
}

void Store_parserString(const amrex::ParmParse& pp, std::string query_string,
                        std::string& stored_string)
{
    std::vector<std::string> f;
    pp.getarr(query_string.c_str(), f);
    stored_string.clear();
    for (auto const& s : f) {
        stored_string += s;
    }
    f.clear();
}

WarpXParser makeParser (std::string const& parse_function, std::vector<std::string> const& varnames)
{
    // Since queryWithParser recursively calls this routine, keep track of symbols
    // in case an infinite recursion is found (a symbol's value depending on itself).
    static std::set<std::string> recursive_symbols;

    WarpXParser parser(parse_function);
    parser.registerVariables(varnames);
    ParmParse pp_my_constants("my_constants");
    std::set<std::string> symbols = parser.symbols();
    for (auto const& v : varnames) symbols.erase(v.c_str());
    for (auto it = symbols.begin(); it != symbols.end(); ) {
        Real v;

        WarpXUtilMsg::AlwaysAssert(recursive_symbols.count(*it)==0, "Expressions contains recursive symbol "+*it);
        recursive_symbols.insert(*it);
        const bool is_input = queryWithParser(pp_my_constants, it->c_str(), v);
        recursive_symbols.erase(*it);

        if (is_input) {
            parser.setConstant(*it, v);
            it = symbols.erase(it);
        } else if (std::strcmp(it->c_str(), "q_e") == 0) {
            parser.setConstant(*it, PhysConst::q_e);
            it = symbols.erase(it);
        } else if (std::strcmp(it->c_str(), "m_e") == 0) {
            parser.setConstant(*it, PhysConst::m_e);
            it = symbols.erase(it);
        } else if (std::strcmp(it->c_str(), "m_p") == 0) {
            parser.setConstant(*it, PhysConst::m_p);
            it = symbols.erase(it);
        } else if (std::strcmp(it->c_str(), "epsilon0") == 0) {
            parser.setConstant(*it, PhysConst::ep0);
            it = symbols.erase(it);
        }  else if (std::strcmp(it->c_str(), "mu0") == 0) {
            parser.setConstant(*it, PhysConst::mu0);
            it = symbols.erase(it);
        } else if (std::strcmp(it->c_str(), "clight") == 0) {
            parser.setConstant(*it, PhysConst::c);
            it = symbols.erase(it);
        } else if (std::strcmp(it->c_str(), "pi") == 0) {
            parser.setConstant(*it, MathConst::pi);
            it = symbols.erase(it);
        } else {
            ++it;
        }
    }
    for (auto const& s : symbols) {
        amrex::Abort("makeParser::Unknown symbol "+s);
    }
    return parser;
}

int
queryWithParser (const amrex::ParmParse& a_pp, char const * const str, amrex::Real& val)
{
    // call amrex::ParmParse::query, check if the user specified str.
    std::string tmp_str;
    int is_specified = a_pp.query(str, tmp_str);
    if (is_specified)
    {
        // If so, create a parser object and apply it to the value provided by the user.
        std::string str_val;
        Store_parserString(a_pp, str, str_val);

        auto parser = makeParser(str_val, {});
        val = parser.eval();
    }
    // return the same output as amrex::ParmParse::query
    return is_specified;
}

void
getWithParser (const amrex::ParmParse& a_pp, char const * const str, amrex::Real& val)
{
    // If so, create a parser object and apply it to the value provided by the user.
    std::string str_val;
    Store_parserString(a_pp, str, str_val);

    auto parser = makeParser(str_val, {});
    val = parser.eval();
}

int
queryArrWithParser (const amrex::ParmParse& a_pp, char const * const str, std::vector<amrex::Real>& val,
                    const int start_ix, const int num_val)
{
    // call amrex::ParmParse::query, check if the user specified str.
    std::vector<std::string> tmp_str_arr;
    int is_specified = a_pp.queryarr(str, tmp_str_arr, start_ix, num_val);
    if (is_specified)
    {
        // If so, create parser objects and apply them to the values provided by the user.
        int const n = static_cast<int>(tmp_str_arr.size());
        val.resize(n);
        for (int i=0 ; i < n ; i++) {
            auto parser = makeParser(tmp_str_arr[i], {});
            val[i] = parser.eval();
        }
    }
    // return the same output as amrex::ParmParse::query
    return is_specified;
}

void
getArrWithParser (const amrex::ParmParse& a_pp, char const * const str, std::vector<amrex::Real>& val)
{
    // Create parser objects and apply them to the values provided by the user.
    std::vector<std::string> tmp_str_arr;
    a_pp.getarr(str, tmp_str_arr);

    int const n = static_cast<int>(tmp_str_arr.size());
    val.resize(n);
    for (int i=0 ; i < n ; i++) {
        auto parser = makeParser(tmp_str_arr[i], {});
        val[i] = parser.eval();
    }
}

void
getArrWithParser (const amrex::ParmParse& a_pp, char const * const str, std::vector<amrex::Real>& val,
                    const int start_ix, const int num_val)
{
    // Create parser objects and apply them to the values provided by the user.
    std::vector<std::string> tmp_str_arr;
    a_pp.getarr(str, tmp_str_arr, start_ix, num_val);

    int const n = static_cast<int>(tmp_str_arr.size());
    val.resize(n);
    for (int i=0 ; i < n ; i++) {
        auto parser = makeParser(tmp_str_arr[i], {});
        val[i] = parser.eval();
    }
}

/**
 * \brief Ensures that the blocks are setup correctly for the RZ spectral solver
 * When using the RZ spectral solver, the Hankel transform cannot be
 * divided among multiple blocks. Each block must extend over the
 * entire radial extent.
 * The grid can be divided up along z, but the number of blocks
 * must be >= the number of processors.
 */
void CheckGriddingForRZSpectral ()
{
#ifndef WARPX_DIM_RZ
    amrex::Abort("CheckGriddingForRZSpectral: WarpX was not built with RZ geometry.");
#else

    // Ensure that geometry.coord_sys is set properly.
    ParmParse pp_geometry("geometry");
    int coord_sys = 1;
    pp_geometry.query("coord_sys", coord_sys);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(coord_sys == 1, "geometry.coord_sys needs to be 1 when using cylindrical geometry");
    pp_geometry.add("coord_sys", coord_sys);

    ParmParse pp_algo("algo");
    int maxwell_solver_id = GetAlgorithmInteger(pp_algo, "maxwell_solver");

    // only check for PSATD in RZ
    if (maxwell_solver_id != MaxwellSolverAlgo::PSATD)
        return;

    int max_level;
    Vector<int> n_cell(AMREX_SPACEDIM, -1);

    ParmParse pp_amr("amr");

    pp_amr.get("max_level",max_level);
    pp_amr.getarr("n_cell",n_cell,0,AMREX_SPACEDIM);

    Vector<int> blocking_factor_x(max_level+1);
    Vector<int> max_grid_size_x(max_level+1);

    // Set the radial block size to be equal to the radial grid size.
    blocking_factor_x[0] = n_cell[0];
    max_grid_size_x[0] = n_cell[0];

    for (int lev=1 ; lev <= max_level ; lev++) {
        // For this to be correct, this needs to read in any user specified refinement ratios.
        // But since that is messy and unlikely to be needed anytime soon, the ratio is
        // fixed to 2 which will be the most likely value.
        blocking_factor_x[lev] = blocking_factor_x[lev-1]*2; // refRatio(lev-1);
        max_grid_size_x[lev] = max_grid_size_x[lev-1]*2; // refRatio(lev-1);
    }

    // Note that any user input values for these parameters are discarded.
    pp_amr.addarr("blocking_factor_x", blocking_factor_x);
    pp_amr.addarr("max_grid_size_x", max_grid_size_x);

    // Adjust the longitudinal block sizes, making sure that there are
    // more blocks than processors.
    // The factor of 8 is there to make some room for higher order
    // shape factors and filtering.
    int nprocs = ParallelDescriptor::NProcs();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_cell[1] >= 8*nprocs,
                                     "With RZ spectral, there must be at least eight z-cells per processor so that there can be at least one block per processor.");

    // Get the longitudinal blocking factor in case it was set by the user.
    // If not set, use the default value of 8.
    Vector<int> bf;
    pp_amr.queryarr("blocking_factor",bf);
    pp_amr.queryarr("blocking_factor_y",bf);
    bf.resize(std::max(static_cast<int>(bf.size()),1),8);

    // Modify the default or any user input, making sure that the blocking factor
    // is small enough so that there will be at least as many blocks as there are
    // processors. Because of the ASSERT above, bf will never be less than 8.
    while (n_cell[1] < nprocs*bf[0]) {
        bf[0] /= 2;
    }
    pp_amr.addarr("blocking_factor_y", bf);

    // Get the longitudinal max grid size in case it was set by the user.
    // If not set, use the default value of 128.
    Vector<int> mg;
    pp_amr.queryarr("max_grid_size",mg);
    pp_amr.queryarr("max_grid_size_y",mg);
    mg.resize(std::max(static_cast<int>(mg.size()),1),128);

    // Modify the default or any user input, making sure that the max grid size
    // (of the coarsest level) is small enough so that there will be at least
    // as many blocks as there are processors.
    while (n_cell[1] < nprocs*mg[0]) {
        mg[0] /= 2;
    }
    pp_amr.addarr("max_grid_size_y", mg);
#endif
}


void ReadBCParams ()
{

    amrex::Vector<std::string> field_BC_lo(AMREX_SPACEDIM,"default");
    amrex::Vector<std::string> field_BC_hi(AMREX_SPACEDIM,"default");
    amrex::Vector<std::string> particle_BC_lo(AMREX_SPACEDIM,"default");
    amrex::Vector<std::string> particle_BC_hi(AMREX_SPACEDIM,"default");
    amrex::Vector<int> geom_periodicity(AMREX_SPACEDIM,0);
    ParmParse pp_geometry("geometry");
    ParmParse pp_warpx("warpx");
    ParmParse pp_algo("algo");
    int maxwell_solver_id = GetAlgorithmInteger(pp_algo, "maxwell_solver");
    if (pp_geometry.queryarr("is_periodic", geom_periodicity)) {
        // set default field and particle boundary appropriately
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (geom_periodicity[idim] == 1) {
                // set boundary to periodic based on user-defined periodicity
                WarpX::field_boundary_lo[idim] = FieldBoundaryType::Periodic;
                WarpX::field_boundary_hi[idim] = FieldBoundaryType::Periodic;
                WarpX::particle_boundary_lo[idim] = ParticleBoundaryType::Periodic;
                WarpX::particle_boundary_hi[idim] = ParticleBoundaryType::Periodic;
            } else {
                // if non-periodic and do_pml=0, then set default boundary to PEC
                int pml_input = 1;
                int silverMueller_input = 0;
                pp_warpx.query("do_pml", pml_input);
                pp_warpx.query("do_silver_mueller", silverMueller_input);
                if (pml_input == 0 and silverMueller_input == 0) {
                    if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
                        WarpX::field_boundary_lo[idim] = FieldBoundaryType::None;
                        WarpX::field_boundary_hi[idim] = FieldBoundaryType::None;
                    } else {
                        WarpX::field_boundary_lo[idim] = FieldBoundaryType::PEC;
                        WarpX::field_boundary_hi[idim] = FieldBoundaryType::PEC;
                    }
#ifdef WARPX_DIM_RZ
                    if (idim == 0) WarpX::field_boundary_lo[idim] = FieldBoundaryType::None;
#endif
                }
            }
        }
        // Temporarily setting default boundary to Damped until new boundary interface is introduced
        if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
            ParmParse pp_psatd("psatd");
            int do_moving_window = 0;
            pp_warpx.query("do_moving_window", do_moving_window);
            if (do_moving_window == 1) {
                std::string s;
                pp_warpx.get("moving_window_dir", s);
                int zdir;
                if (s == "z" || s == "Z") {
                    zdir = AMREX_SPACEDIM-1;
                    WarpX::field_boundary_lo[zdir] = FieldBoundaryType::Damped;
                    WarpX::field_boundary_hi[zdir] = FieldBoundaryType::Damped;
                }
            }
        }
        return;
        // When all boundary conditions are supported, the abort statement below will be introduced
        //amrex::Abort("geometry.is_periodic is not supported. Please use `boundary.field_lo`, `boundary.field_hi` to specifiy field boundary conditions and 'boundary.particle_lo', 'boundary.particle_hi'  to specify particle boundary conditions.");
    }
    // particle boundary may not be explicitly specified for some applications
    bool particle_boundary_specified = false;
    ParmParse pp_boundary("boundary");
    pp_boundary.queryarr("field_lo", field_BC_lo, 0, AMREX_SPACEDIM);
    pp_boundary.queryarr("field_hi", field_BC_hi, 0, AMREX_SPACEDIM);
    if (pp_boundary.queryarr("particle_lo", particle_BC_lo, 0, AMREX_SPACEDIM))
        particle_boundary_specified = true;
    if (pp_boundary.queryarr("particle_hi", particle_BC_hi, 0, AMREX_SPACEDIM))
        particle_boundary_specified = true;
    AMREX_ALWAYS_ASSERT(field_BC_lo.size() == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(field_BC_hi.size() == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(particle_BC_lo.size() == AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(particle_BC_hi.size() == AMREX_SPACEDIM);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        // Get field boundary type
        WarpX::field_boundary_lo[idim] = GetBCTypeInteger(field_BC_lo[idim], true);
        WarpX::field_boundary_hi[idim] = GetBCTypeInteger(field_BC_hi[idim], true);
        // Get particle boundary type
        WarpX::particle_boundary_lo[idim] = GetBCTypeInteger(particle_BC_lo[idim], false);
        WarpX::particle_boundary_hi[idim] = GetBCTypeInteger(particle_BC_hi[idim], false);

        if (WarpX::field_boundary_lo[idim] == FieldBoundaryType::Periodic ||
            WarpX::field_boundary_hi[idim] == FieldBoundaryType::Periodic ||
            WarpX::particle_boundary_lo[idim] == ParticleBoundaryType::Periodic ||
            WarpX::particle_boundary_hi[idim] == ParticleBoundaryType::Periodic ) {
            geom_periodicity[idim] = 1;
            // to ensure both lo and hi are set to periodic consistently for both field and particles.
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                (WarpX::field_boundary_lo[idim]  == FieldBoundaryType::Periodic) &&
                (WarpX::field_boundary_hi[idim]  == FieldBoundaryType::Periodic),
            "field boundary must be consistenly periodic in both lo and hi");
            if (particle_boundary_specified) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                    (WarpX::particle_boundary_lo[idim] == ParticleBoundaryType::Periodic) &&
                    (WarpX::particle_boundary_hi[idim] == ParticleBoundaryType::Periodic),
               "field and particle boundary must be periodic in both lo and hi");
            } else {
                // set particle boundary to periodic
                WarpX::particle_boundary_lo[idim] = ParticleBoundaryType::Periodic;
                WarpX::particle_boundary_hi[idim] = ParticleBoundaryType::Periodic;
            }
        }
        if (maxwell_solver_id == MaxwellSolverAlgo::PSATD) {
            if (WarpX::field_boundary_lo[idim] == FieldBoundaryType::PEC ||
                WarpX::field_boundary_hi[idim] == FieldBoundaryType::PEC) {
                amrex::Abort(" PEC boundary not implemented for PSATD, yet!");
            }
        }
    }
#ifdef WARPX_DIM_RZ
    // Ensure code aborts if PEC is specified at r=0 for RZ
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE( WarpX::field_boundary_lo[0] == FieldBoundaryType::None,
        "Error : Field boundary at r=0 must be ``none``. \n");
#endif

    pp_geometry.addarr("is_periodic", geom_periodicity);
}

namespace WarpXUtilMsg{

void AlwaysAssert(bool is_expression_true, const std::string& msg = "ERROR!")
{
    if(is_expression_true) return;

    amrex::Abort(msg);
}

}

namespace WarpXUtilStr
{
    bool is_in(const std::vector<std::string>& vect,
               const std::string& elem)
    {
        bool value = false;
        if (std::find(vect.begin(), vect.end(), elem) != vect.end()){
            value = true;
        }
        return value;
    }

    bool is_in(const std::vector<std::string>& vect,
               const std::vector<std::string>& elems)
    {
        bool value = false;
        for (auto elem : elems){
            if (is_in(vect, elem)) value = true;
        }
        return value;
    }

}
