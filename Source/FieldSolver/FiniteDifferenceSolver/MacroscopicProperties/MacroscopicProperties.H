#ifndef WARPX_MACROSCOPICPROPERTIES_H_
#define WARPX_MACROSCOPICPROPERTIES_H_


#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>

#include <AMReX_MultiFab.H>
#include <algorithm>
#include <array>

#include "Parser/WarpXParser.H"
#include "Parser/WarpXParserWrapper.H"


/**
 * \brief This class contains the macroscopic parameters of the medium needed to
 * evaluate macroscopic Maxwell equation.
 */
class
MacroscopicProperties
{
public:
     MacroscopicProperties (); // constructor
     /** \brief Read user-defined macroscopic properties. Called in constructor. */
     void ReadParameters ();
     /** Initialize multifabs storing macroscopic multifabs */
     void InitData ();

//     /** return Real, sigma (conductivity) of the medium. */
//     amrex::Real sigma () const noexcept {return m_sigma;}
//     /** return Real, epsilon (permittivity) of the medium. */
//     amrex::Real epsilon () const noexcept {return m_epsilon;}
//     /** return Real, mu (permeability) of the medium. */
//     amrex::Real mu () const noexcept {return m_mu;}
     /** */
     amrex::MultiFab& getsigma_mf  () {return (*m_sigma_mf);}
     amrex::MultiFab& getepsilon_mf  () {return (*m_eps_mf);}
     amrex::MultiFab& getmu_mf  () {return (*m_mu_mf);}

private:
     /** Conductivity, sigma, of the medium */
     amrex::Real m_sigma;
     /** Permittivity, epsilon, of the medium */
     amrex::Real m_epsilon;
     /** Permeability, mu, of the medium */
     amrex::Real m_mu;
     /** Multifab for m_sigma */
     std::unique_ptr<amrex::MultiFab> m_sigma_mf;
     /** Multifab for m_epsilon */
     std::unique_ptr<amrex::MultiFab> m_eps_mf;
     /** Multifab for m_mu */
     std::unique_ptr<amrex::MultiFab> m_mu_mf;
     std::string m_sigma_s;
     std::string m_epsilon_s;
     std::string m_mu_s;

     void InitializeMacroMultiFabUsingParser (amrex::MultiFab *macro_mf,
                                  ParserWrapper<3> *macro_parser, int lev);

     /** string for storing parser function */
     std::string m_str_sigma_function;
     std::string m_str_epsilon_function;
     std::string m_str_mu_function;
     /** Parser Wrappers */
     std::unique_ptr<ParserWrapper<3> > m_sigma_parser;
     std::unique_ptr<ParserWrapper<3> > m_epsilon_parser;
     std::unique_ptr<ParserWrapper<3> > m_mu_parser;
};

struct LaxWendroffAlgo {

     // iv(1,0,0) for cell-centered macroparameters to Ex edge
     // iv(0,1,0) for cell-centered macroparameters to Ey edge
     // iv(0,0,1) for cell-centered macroparameters to Ez edge
     // iv(0,1,1) for nodal macroparameters to Ex edge
     // iv(1,0,1) for nodal macroparameters to Ey edge
     // iv(1,1,0) for nodal macroparameters to Ez edge
     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real macro_avg_to_edge (int i, int j, int k, amrex::IntVect iv,
                                           amrex::Array4<amrex::Real> const& macro_prop) {
         using namespace amrex;
         return ( 0.125_rt * ( macro_prop(i,         j,         k        )
                             + macro_prop(i+iv[0]-1, j,         k        )
                             + macro_prop(i,         j+iv[1]-1, k        )
                             + macro_prop(i,         j,         k+iv[2]-1)
                             + macro_prop(i+iv[0]-1, j+iv[1]-1, k        )
                             + macro_prop(i+iv[0]-1, j        , k+iv[2]-1)
                             + macro_prop(i,         j+iv[1]-1, k+iv[2]-1)
                             + macro_prop(i+iv[0]-1, j+iv[1]-1, k+iv[2]-1) ) );
     }

     /** */
     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real alpha (amrex::Array4<amrex::Real> const& sigma,
                               amrex::Array4<amrex::Real> const& epsilon,
                               amrex::Real dt, int i, int j, int k,
                               amrex::IntVect iv_cc) {
         using namespace amrex;
         amrex::Real fac1 = 0.5_rt * macro_avg_to_edge(i, j, k, iv_cc, sigma) * dt
                                   / macro_avg_to_edge(i, j, k, iv_cc, epsilon);
         amrex::Real alpha = (1._rt - fac1)/(1._rt + fac1);
         return alpha;
     };

     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real beta (amrex::Array4<amrex::Real> const& sigma,
                              amrex::Array4<amrex::Real> const& epsilon,
                              amrex::Real dt, int i, int j, int k,
                              amrex::IntVect iv_cc) {
         using namespace amrex;
         amrex::Real fac1 = 0.5_rt * macro_avg_to_edge(i, j, k, iv_cc, sigma) * dt
                                   / macro_avg_to_edge(i, j, k, iv_cc, epsilon);
         amrex::Real beta = dt / ( macro_avg_to_edge(i, j, k, iv_cc, epsilon) * (1._rt + fac1) );
         return beta;
     };
};

struct BackwardEulerAlgo {

     // iv(1,0,0) for cell-centered macroparameters to Ex edge
     // iv(0,1,0) for cell-centered macroparameters to Ey edge
     // iv(0,0,1) for cell-centered macroparameters to Ez edge
     // iv(2,1,1) for nodal macroparameters to Ex edge
     // iv(1,2,1) for nodal macroparameters to Ey edge
     // iv(1,1,2) for nodal macroparameters to Ez edge
     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real macro_avg_to_edge (int i, int j, int k, amrex::IntVect iv,
                                           amrex::Array4<amrex::Real> const& macro_prop) {
         using namespace amrex;
         return ( 0.125_rt * ( macro_prop(i,         j,         k        )
                             + macro_prop(i+iv[0]-1, j,         k        )
                             + macro_prop(i,         j+iv[1]-1, k        )
                             + macro_prop(i,         j,         k+iv[2]-1)
                             + macro_prop(i+iv[0]-1, j+iv[1]-1, k        )
                             + macro_prop(i+iv[0]-1, j        , k+iv[2]-1)
                             + macro_prop(i,         j+iv[1]-1, k+iv[2]-1)
                             + macro_prop(i+iv[0]-1, j+iv[1]-1, k+iv[2]-1) ) );
     }

     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real alpha (amrex::Array4<amrex::Real> const& sigma,
                               amrex::Array4<amrex::Real> const& epsilon,
                               amrex::Real dt, int i, int j, int k,
                               amrex::IntVect iv_cc) {
         using namespace amrex;
         amrex::Real fac1 = macro_avg_to_edge(i, j, k, iv_cc, sigma) * dt
                            / macro_avg_to_edge(i, j, k, iv_cc, epsilon);
         amrex::Real alpha = 1._rt/ (1._rt + fac1);
         return alpha;
     };

     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     static amrex::Real beta (amrex::Array4<amrex::Real> const& sigma,
                              amrex::Array4<amrex::Real> const& epsilon,
                              amrex::Real dt, int i, int j, int k,
                              amrex::IntVect iv_cc) {
         using namespace amrex;
         amrex::Real fac1 = macro_avg_to_edge(i, j, k, iv_cc, sigma) * dt
                            / macro_avg_to_edge(i, j, k, iv_cc, epsilon);
         amrex::Real beta = dt / ( macro_avg_to_edge(i, j, k, iv_cc, epsilon) * (1._rt + fac1) );
         return beta;
     };
};

#endif // WARPX_MACROSCOPIC_PROPERTIES_H_
