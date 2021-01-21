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


using namespace amrex;

void ReadBoostedFrameParameters(Real& gamma_boost, Real& beta_boost,
                                Vector<int>& boost_direction)
{
    ParmParse pp("warpx");
    queryWithParser(pp, "gamma_boost", gamma_boost);
    if( gamma_boost > 1. ) {
        beta_boost = std::sqrt(1.-1./pow(gamma_boost,2));
        std::string s;
        pp.get("boost_direction", s);
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

    ParmParse pp_geom("geometry");
    ParmParse pp_wpx("warpx");
    ParmParse pp_amr("amr");
    ParmParse pp_slice("slice");

    pp_geom.getarr("prob_lo",prob_lo,0,AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_lo.size() == AMREX_SPACEDIM);
    pp_geom.getarr("prob_hi",prob_hi,0,AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

    pp_slice.queryarr("dom_lo",slice_lo,0,AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(slice_lo.size() == AMREX_SPACEDIM);
    pp_slice.queryarr("dom_hi",slice_hi,0,AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(slice_hi.size() == AMREX_SPACEDIM);


    pp_amr.query("max_level", max_level);
    if (max_level > 0){
      pp_wpx.getarr("fine_tag_lo", fine_tag_lo);
      pp_wpx.getarr("fine_tag_hi", fine_tag_hi);
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

    pp_geom.addarr("prob_lo", prob_lo);
    pp_geom.addarr("prob_hi", prob_hi);
    if (max_level > 0){
      pp_wpx.addarr("fine_tag_lo", fine_tag_lo);
      pp_wpx.addarr("fine_tag_hi", fine_tag_hi);
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
    WarpXParser parser(parse_function);
    parser.registerVariables(varnames);
    ParmParse pp("my_constants");
    std::set<std::string> symbols = parser.symbols();
    for (auto const& v : varnames) symbols.erase(v.c_str());
    for (auto it = symbols.begin(); it != symbols.end(); ) {
        Real v;
        if (pp.query(it->c_str(), v)) {
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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "CheckGriddingForRZSpectral: WarpX was not built with RZ geometry.");
#endif

    ParmParse pp("algo");
    int maxwell_solver_id = GetAlgorithmInteger(pp, "maxwell_solver");

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
