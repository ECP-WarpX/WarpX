#ifndef WARPX_MACROSCOPICPROPERTIES_H_
#define WARPX_MACROSCOPICPROPERTIES_H_


#include <AMReX_MultiFab.H>

class
MacroscopicProperties
{
public :
     MacroscopicProperties ();
     /** \brief Read input parameters. This function is called in class constructor. */
     void ReadParameters ();
     amrex::Real sigma () {return m_sigma;}
     amrex::Real epsilon () {return m_epsilon;}
     amrex::Real mu () {return m_mu;}
private :
     amrex::Real m_sigma; // conductivity (constant for now)
     amrex::Real m_epsilon; // permittivity (constant for now)
     amrex::Real m_mu;  // permeability (constant for now)
};

#endif // WARPX_MACROSCOPIC_PROPERTIES_H_
