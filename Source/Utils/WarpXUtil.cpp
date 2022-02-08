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
#include "WarpXProfilerWrapper.H"
#include "WarpXUtil.H"

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Config.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <set>
#include <string>
#include <limits>

using namespace amrex;

void PreparseAMReXInputIntArray(amrex::ParmParse& a_pp, char const * const input_str, const bool replace)
{
    const int cnt = a_pp.countval(input_str);
    if (cnt > 0) {
        Vector<int> input_array;
        getArrWithParser(a_pp, input_str, input_array);
        if (replace) {
            a_pp.remove(input_str);
        }
        a_pp.addarr(input_str, input_array);
    }
}

void ParseGeometryInput()
{
    // Ensure that geometry.dims is set properly.
    CheckDims();

    // Parse prob_lo and hi, evaluating any expressions since geometry does not
    // parse its input
    ParmParse pp_geometry("geometry");

    Vector<Real> prob_lo(AMREX_SPACEDIM);
    Vector<Real> prob_hi(AMREX_SPACEDIM);

    getArrWithParser(pp_geometry, "prob_lo", prob_lo, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_lo.size() == AMREX_SPACEDIM);
    getArrWithParser(pp_geometry, "prob_hi", prob_hi, 0, AMREX_SPACEDIM);
    AMREX_ALWAYS_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

#ifdef WARPX_DIM_RZ
    ParmParse pp_algo("algo");
    int maxwell_solver_id = GetAlgorithmInteger(pp_algo, "maxwell_solver");
    if (maxwell_solver_id == MaxwellSolverAlgo::PSATD)
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(prob_lo[0] == 0.,
            "Lower bound of radial coordinate (prob_lo[0]) with RZ PSATD solver must be zero");
    }
    else
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(prob_lo[0] >= 0.,
            "Lower bound of radial coordinate (prob_lo[0]) with RZ FDTD solver must be non-negative");
    }
#endif

    pp_geometry.addarr("prob_lo", prob_lo);
    pp_geometry.addarr("prob_hi", prob_hi);

    // Parse amr input, evaluating any expressions since amr does not parse its input
    ParmParse pp_amr("amr");

    // Note that n_cell is replaced so that only the parsed version is written out to the
    // warpx_job_info file. This must be done since yt expects to be able to parse
    // the value of n_cell from that file. For the rest, this doesn't matter.
    PreparseAMReXInputIntArray(pp_amr, "n_cell", true);
    PreparseAMReXInputIntArray(pp_amr, "max_grid_size", false);
    PreparseAMReXInputIntArray(pp_amr, "max_grid_size_x", false);
    PreparseAMReXInputIntArray(pp_amr, "max_grid_size_y", false);
    PreparseAMReXInputIntArray(pp_amr, "max_grid_size_z", false);
    PreparseAMReXInputIntArray(pp_amr, "blocking_factor", false);
    PreparseAMReXInputIntArray(pp_amr, "blocking_factor_x", false);
    PreparseAMReXInputIntArray(pp_amr, "blocking_factor_y", false);
    PreparseAMReXInputIntArray(pp_amr, "blocking_factor_z", false);
}

void ReadBoostedFrameParameters(Real& gamma_boost, Real& beta_boost,
                                Vector<int>& boost_direction)
{
    ParmParse pp_warpx("warpx");
    queryWithParser(pp_warpx, "gamma_boost", gamma_boost);
    if( gamma_boost > 1. ) {
        beta_boost = std::sqrt(1._rt-1._rt/std::pow(gamma_boost,2._rt));
        std::string s;
        pp_warpx.get("boost_direction", s);
        if (s == "x" || s == "X") {
            boost_direction[0] = 1;
        }
#if defined(WARPX_DIM_3D)
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


#if defined(WARPX_DIM_3D)
    Vector<int> dim_map {0, 1, 2};
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    Vector<int> dim_map {0, 2};
#else
    Vector<int> dim_map {2};
#endif

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (boost_direction[dim_map[idim]]) {
            amrex::Real convert_factor;
            // Assume that the window travels with speed +c
            convert_factor = 1._rt/( gamma_boost * ( 1 - beta_boost ) );
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
#if defined(WARPX_DIM_3D)
            amrex::Array<amrex::Real,3> galilean_shift = { 0._rt, 0._rt, 0._rt, };
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            amrex::Array<amrex::Real,3> galilean_shift = { 0._rt, std::numeric_limits<amrex::Real>::quiet_NaN(),  0._rt, } ;
#elif defined(WARPX_DIM_1D_Z)
            amrex::Array<amrex::Real,3> galilean_shift = {std::numeric_limits<amrex::Real>::quiet_NaN(), std::numeric_limits<amrex::Real>::quiet_NaN(),  0._rt, } ;
#endif
        const amrex::Real zmin_box = WarpX::LowerCorner(bx, galilean_shift, lev)[2];
        const amrex::Real zmax_box = WarpX::UpperCorner(bx, lev)[2];
        amrex::Real dz  = WarpX::CellSize(lev)[2];
        // Get box lower index in the z direction
#if defined(WARPX_DIM_3D)
        const int lo_ind = bx.loVect()[2];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        const int lo_ind = bx.loVect()[1];
#else
        const int lo_ind = bx.loVect()[0];
#endif
        // Check if box intersect with [zmin, zmax]
        if ( (zmax>zmin_box && zmin<=zmax_box) ){
            Array4<Real> arr = mf[mfi].array();
            // Set field to 0 between zmin and zmax
            ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept{
#if defined(WARPX_DIM_3D)
                    const Real z_gridpoint = zmin_box+(k-lo_ind)*dz;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                    const Real z_gridpoint = zmin_box+(j-lo_ind)*dz;
#else
                    const Real z_gridpoint = zmin_box+(i-lo_ind)*dz;
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

int safeCastToInt(const amrex::Real x, const std::string& real_name) {
    int result = 0;
    bool error_detected = false;
    std::string assert_msg;
    // (2.0*(numeric_limits<int>::max()/2+1)) converts numeric_limits<int>::max()+1 to a real ensuring accuracy to all digits
    // This accepts x = 2**31-1 but rejects 2**31.
    using namespace amrex::literals;
    constexpr amrex::Real max_range = (2.0_rt*static_cast<amrex::Real>(std::numeric_limits<int>::max()/2+1));
    if (x < max_range) {
        if (std::ceil(x) >= std::numeric_limits<int>::min()) {
            result = static_cast<int>(x);
        } else {
            error_detected = true;
            assert_msg = "Error: Negative overflow detected when casting " + real_name + " = " + std::to_string(x) + " to int";
        }
    } else if (x > 0) {
        error_detected = true;
        assert_msg =  "Error: Overflow detected when casting " + real_name + " = " + std::to_string(x) + " to int";
    } else {
        error_detected = true;
        assert_msg =  "Error: NaN detected when casting " + real_name + " to int";
    }
    WarpXUtilMsg::AlwaysAssert(!error_detected, assert_msg);
    return result;
}

Parser makeParser (std::string const& parse_function, amrex::Vector<std::string> const& varnames)
{
    // Since queryWithParser recursively calls this routine, keep track of symbols
    // in case an infinite recursion is found (a symbol's value depending on itself).
    static std::set<std::string> recursive_symbols;

    Parser parser(parse_function);
    parser.registerVariables(varnames);

    std::set<std::string> symbols = parser.symbols();
    for (auto const& v : varnames) symbols.erase(v.c_str());

    // User can provide inputs under this name, through which expressions
    // can be provided for arbitrary variables. PICMI inputs are aware of
    // this convention and use the same prefix as well. This potentially
    // includes variable names that match physical or mathematical
    // constants, in case the user wishes to enforce a different
    // system of units or some form of quasi-physical behavior in the
    // simulation. Thus, this needs to override any built-in
    // constants.
    ParmParse pp_my_constants("my_constants");

    // Physical / Numerical Constants available to parsed expressions
    static std::map<std::string, amrex::Real> warpx_constants =
      {
       {"clight", PhysConst::c},
       {"epsilon0", PhysConst::ep0},
       {"mu0", PhysConst::mu0},
       {"q_e", PhysConst::q_e},
       {"m_e", PhysConst::m_e},
       {"m_p", PhysConst::m_p},
       {"m_u", PhysConst::m_u},
       {"kb", PhysConst::kb},
       {"pi", MathConst::pi},
      };

    for (auto it = symbols.begin(); it != symbols.end(); ) {
        // Always parsing in double precision avoids potential overflows that may occur when parsing
        // user's expressions because of the limited range of exponentials in single precision
        double v;

        WarpXUtilMsg::AlwaysAssert(recursive_symbols.count(*it)==0, "Expressions contains recursive symbol "+*it);
        recursive_symbols.insert(*it);
        const bool is_input = queryWithParser(pp_my_constants, it->c_str(), v);
        recursive_symbols.erase(*it);

        if (is_input) {
            parser.setConstant(*it, v);
            it = symbols.erase(it);
            continue;
        }

        auto constant = warpx_constants.find(*it);
        if (constant != warpx_constants.end()) {
          parser.setConstant(*it, constant->second);
          it = symbols.erase(it);
          continue;
        }

        ++it;
    }
    for (auto const& s : symbols) {
        amrex::Abort("makeParser::Unknown symbol "+s);
    }
    return parser;
}

double
parseStringtoReal(std::string str)
{
    auto parser = makeParser(str, {});
    auto exe = parser.compileHost<0>();
    double result = exe();
    return result;
}

int
parseStringtoInt(std::string str, std::string name)
{
    auto const rval = static_cast<amrex::Real>(parseStringtoReal(str));
    int ival = safeCastToInt(std::round(rval), name);
    return ival;
}

// Overloads for float/double instead of amrex::Real to allow makeParser() to query for
// my_constants as double even in single precision mode
// Always parsing in double precision avoids potential overflows that may occur when parsing user's
// expressions because of the limited range of exponentials in single precision
int
queryWithParser (const amrex::ParmParse& a_pp, char const * const str, float& val)
{
    // call amrex::ParmParse::query, check if the user specified str.
    std::string tmp_str;
    int is_specified = a_pp.query(str, tmp_str);
    if (is_specified)
    {
        // If so, create a parser object and apply it to the value provided by the user.
        std::string str_val;
        Store_parserString(a_pp, str, str_val);
        val = static_cast<float>(parseStringtoReal(str_val));
    }
    // return the same output as amrex::ParmParse::query
    return is_specified;
}

void
getWithParser (const amrex::ParmParse& a_pp, char const * const str, float& val)
{
    // If so, create a parser object and apply it to the value provided by the user.
    std::string str_val;
    Store_parserString(a_pp, str, str_val);
    val = static_cast<float>(parseStringtoReal(str_val));
}

int
queryWithParser (const amrex::ParmParse& a_pp, char const * const str, double& val)
{
    // call amrex::ParmParse::query, check if the user specified str.
    std::string tmp_str;
    int is_specified = a_pp.query(str, tmp_str);
    if (is_specified)
    {
        // If so, create a parser object and apply it to the value provided by the user.
        std::string str_val;
        Store_parserString(a_pp, str, str_val);
        val = parseStringtoReal(str_val);
    }
    // return the same output as amrex::ParmParse::query
    return is_specified;
}

void
getWithParser (const amrex::ParmParse& a_pp, char const * const str, double& val)
{
    // If so, create a parser object and apply it to the value provided by the user.
    std::string str_val;
    Store_parserString(a_pp, str, str_val);
    val = parseStringtoReal(str_val);
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
            val[i] = static_cast<amrex::Real>(parseStringtoReal(tmp_str_arr[i]));
        }
    }
    // return the same output as amrex::ParmParse::query
    return is_specified;
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
        val[i] = static_cast<amrex::Real>(parseStringtoReal(tmp_str_arr[i]));
    }
}

int queryWithParser (const amrex::ParmParse& a_pp, char const * const str, int& val) {
    amrex::Real rval;
    const int result = queryWithParser(a_pp, str, rval);
    if (result) {
        val = safeCastToInt(std::round(rval), str);
    }
    return result;
}

void getWithParser (const amrex::ParmParse& a_pp, char const * const str, int& val) {
    amrex::Real rval;
    getWithParser(a_pp, str, rval);
    val = safeCastToInt(std::round(rval), str);
}

int queryArrWithParser (const amrex::ParmParse& a_pp, char const * const str, std::vector<int>& val,
                        const int start_ix, const int num_val) {
    std::vector<amrex::Real> rval;
    const int result = queryArrWithParser(a_pp, str, rval, start_ix, num_val);
    if (result) {
        val.resize(rval.size());
        for (unsigned long i = 0 ; i < val.size() ; i++) {
            val[i] = safeCastToInt(std::round(rval[i]), str);
        }
    }
    return result;
}

void getArrWithParser (const amrex::ParmParse& a_pp, char const * const str, std::vector<int>& val,
                       const int start_ix, const int num_val) {
    std::vector<amrex::Real> rval;
    getArrWithParser(a_pp, str, rval, start_ix, num_val);
    val.resize(rval.size());
    for (unsigned long i = 0 ; i < val.size() ; i++) {
        val[i] = safeCastToInt(std::round(rval[i]), str);
    }
}

void CheckDims ()
{
    // Ensure that geometry.dims is set properly.
#if defined(WARPX_DIM_3D)
    std::string const dims_compiled = "3";
#elif defined(WARPX_DIM_XZ)
    std::string const dims_compiled = "2";
#elif defined(WARPX_DIM_1D_Z)
    std::string const dims_compiled = "1";
#elif defined(WARPX_DIM_RZ)
    std::string const dims_compiled = "RZ";
#endif
    ParmParse pp_geometry("geometry");
    std::string dims;
    pp_geometry.get("dims", dims);
    std::string dims_error = "ERROR: The selected WarpX executable was built as '";
    dims_error.append(dims_compiled).append("'-dimensional, but the ");
    dims_error.append("inputs file declares 'geometry.dims = ").append(dims).append("'.\n");
    dims_error.append("Please re-compile with a different WarpX_DIMS option or select the right executable name.");
    WarpXUtilMsg::AlwaysAssert(dims == dims_compiled, dims_error);
}

void CheckGriddingForRZSpectral ()
{
#ifdef WARPX_DIM_RZ
    // Ensure that geometry.dims is set properly.
    CheckDims();

    ParmParse pp_algo("algo");
    int maxwell_solver_id = GetAlgorithmInteger(pp_algo, "maxwell_solver");

    // only check for PSATD in RZ
    if (maxwell_solver_id != MaxwellSolverAlgo::PSATD)
        return;

    int max_level;
    Vector<int> n_cell(AMREX_SPACEDIM, -1);

    ParmParse pp_amr("amr");

    pp_amr.get("max_level",max_level);
    pp_amr.getarr("n_cell", n_cell, 0, AMREX_SPACEDIM);

    Vector<int> blocking_factor_x(max_level+1);
    Vector<int> max_grid_size_x(max_level+1);

    // Set the radial block size to be the power of 2 greater than or equal to
    // the number of grid cells. The blocking_factor must be a power of 2
    // and the max_grid_size should be a multiple of the blocking_factor.
    int k = 1;
    while (k < n_cell[0]) {
        k *= 2;
    }
    blocking_factor_x[0] = k;
    max_grid_size_x[0] = k;

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
    pp_amr.queryarr("blocking_factor", bf);
    pp_amr.queryarr("blocking_factor_y", bf);
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
    pp_amr.queryarr("max_grid_size", mg);
    pp_amr.queryarr("max_grid_size_y", mg);
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
        amrex::Abort("geometry.is_periodic is not supported. Please use `boundary.field_lo`, `boundary.field_hi` to specifiy field boundary conditions and 'boundary.particle_lo', 'boundary.particle_hi'  to specify particle boundary conditions.");
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
        WarpX::field_boundary_lo[idim] = GetFieldBCTypeInteger(field_BC_lo[idim]);
        WarpX::field_boundary_hi[idim] = GetFieldBCTypeInteger(field_BC_hi[idim]);
        // Get particle boundary type
        WarpX::particle_boundary_lo[idim] = GetParticleBCTypeInteger(particle_BC_lo[idim]);
        WarpX::particle_boundary_hi[idim] = GetParticleBCTypeInteger(particle_BC_hi[idim]);

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

    // Appending periodicity information to input so that it can be used by amrex
    // to set parameters necessary to define geometry and perform communication
    // such as FillBoundary. The periodicity is 1 if user-define boundary condition is
    // periodic else it is set to 0.
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
        return (std::find(vect.begin(), vect.end(), elem) != vect.end());
    }

    bool is_in(const std::vector<std::string>& vect,
               const std::vector<std::string>& elems)
    {
        return std::any_of(elems.begin(), elems.end(),
            [&](const auto elem){return is_in(vect, elem);});
    }

}

namespace WarpXUtilLoadBalance
{
    bool doCosts (const amrex::LayoutData<amrex::Real>* costs, const amrex::BoxArray ba,
                  const amrex::DistributionMapping& dm)
    {
        bool consistent = costs && (dm == costs->DistributionMap()) &&
            (ba.CellEqual(costs->boxArray())) &&
            (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers);
        return consistent;
    }
}
