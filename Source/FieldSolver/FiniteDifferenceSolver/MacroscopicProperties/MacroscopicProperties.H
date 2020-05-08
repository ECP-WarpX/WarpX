#ifndef WARPX_MACROSCOPICPROPERTIES_H_
#define WARPX_MACROSCOPICPROPERTIES_H_


#include <AMReX_MultiFab.H>


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
     /** return Real, sigma (conductivity) of the medium. */
     amrex::Real sigma () const noexcept {return m_sigma;}
     /** return Real, epsilon (permittivity) of the medium. */
     amrex::Real epsilon () const noexcept {return m_epsilon;}
     /** return Real, mu (permeability) of the medium. */
     amrex::Real mu () const noexcept {return m_mu;}

private:
     /** Conductivity, sigma, of the medium */
     amrex::Real m_sigma;
     /** Permittivity, epsilon, of the medium */
     amrex::Real m_epsilon;
     /** Permeability, mu, of the medium */
     amrex::Real m_mu;
};

#endif // WARPX_MACROSCOPIC_PROPERTIES_H_
