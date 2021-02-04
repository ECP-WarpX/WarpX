#ifndef BTD_PLOTFILE_HEADER_IMPL_H
#define BTD_PLOTFILE_HEADER_IMPL_H

#include <string>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>


class BTDPlotfileHeaderImpl
{
public:
    BTDPlotfileHeaderImpl (std::string const& Headerfile_path);
    ~BTDPlotfileHeaderImpl () = default;

    std::string fileVersion () const noexcept {return m_file_version; }
    int ncomp () const noexcept {return m_nComp; }
    const amrex::Vector<std::string>& varnames () const noexcept {return m_varnames; }

    int spaceDim () const noexcept {return m_spacedim; }

    amrex::Real time () const noexcept {return m_time; }
  
    int finestLevel () const noexcept { return m_finest_level; }

    amrex::Real problo (int dim) const noexcept {return m_prob_lo[dim]; }
    amrex::Real probhi (int dim) const noexcept {return m_prob_hi[dim]; }


    amrex::Array<amrex::Real, AMREX_SPACEDIM> probLo () const noexcept { return m_prob_lo; } 
    amrex::Array<amrex::Real, AMREX_SPACEDIM> probHo () const noexcept { return m_prob_hi; } 
    amrex::Array<amrex::Real, AMREX_SPACEDIM> cellsize () const noexcept {return m_cell_size; }
    amrex::Box probDomain () const noexcept {return m_prob_domain; }
    int timestep () const noexcept {return m_timestep; }
    int coordsys () const noexcept {return m_coordsys; }
    int bwidth () const noexcept {return m_bwidth; }
    int cur_level () const noexcept {return m_cur_level; }
    int numFabs () const noexcept {return m_numFabs; }
    
    const amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM> >& glo () const noexcept {return m_glo; }
    const amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM> >& ghi () const noexcept {return m_ghi; }

    amrex::Array<amrex::Real, AMREX_SPACEDIM> FabLo (int iFab) const noexcept {return m_glo[iFab]; }
    amrex::Array<amrex::Real, AMREX_SPACEDIM> FabHi (int iFab) const noexcept {return m_ghi[iFab]; }
    std::string CellPath () const noexcept {return m_CellPath; }

    void ReadHeaderData ();
    void WriteHeader ();
 
    void set_time ( amrex::Real new_time) { m_time = new_time;}
    void set_timestep (int new_timestep) { m_timestep = new_timestep; }
    void set_problo (int dim, amrex::Real lo) { m_prob_lo[dim] = lo;}
    void set_probhi (int dim, amrex::Real hi) { m_prob_hi[dim] = hi;}
    void set_probDomain (amrex::Box domainBox) {m_prob_domain = domainBox; }
    void IncrementNumFabs () { m_numFabs++;}
    void ResizeFabLo () { m_glo.resize(m_numFabs); }
    void ResizeFabHi () { m_ghi.resize(m_numFabs); }
    void AppendNewFabLo (amrex::Array<amrex::Real, AMREX_SPACEDIM> newFabLo);
    void AppendNewFabHi (amrex::Array<amrex::Real, AMREX_SPACEDIM> newFabHi);
private:
    std::string m_Header_path;
    std::string m_file_version;
    int m_nComp;
    amrex::Vector<std::string> m_varnames;
    int m_spacedim;
    amrex::Real m_time;
    int m_finest_level, m_nlevel;
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_prob_lo {{AMREX_D_DECL(0.,0.,0.)}};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_prob_hi {{AMREX_D_DECL(1.,1.,1.)}};
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_cell_size {{AMREX_D_DECL(1.,1.,1.)}};
    amrex::Box m_prob_domain;
    int m_timestep;
    int m_coordsys;
    int m_bwidth;
    int m_cur_level;
    int m_numFabs;
    amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM> > m_glo {{AMREX_D_DECL(0.,0.,0.)}};
    amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM> > m_ghi {{AMREX_D_DECL(1.,1.,1.)}};
    std::string m_CellPath;
        
};

class BTDMultiFabHeaderImpl
{
    public:
    BTDMultiFabHeaderImpl (std::string const& Headerfile_path);
    ~BTDMultiFabHeaderImpl () = default;
    void ReadMultiFabHeader ();
    void WriteMultiFabHeader ();

    int ba_size () {return m_ba_size;}
    amrex::Box ba_box (int ibox) {return m_ba[ibox]; }
    std::string fodPrefix (int ifab) {return m_FabOnDiskPrefix[ifab]; }
    std::string FabName (int ifab) {return m_fabname[ifab]; }
    int FabHead (int ifab) {return m_fabhead[ifab]; }
    amrex::Vector<amrex::Real> minval(int ifab) { return m_minval[ifab];}
    amrex::Vector<amrex::Real> maxval(int ifab) { return m_maxval[ifab];}
    void ResizeFabData ();
    void IncreaseMultiFabSize (int add_mf_size) {m_ba_size += add_mf_size;}
    void SetBox (int ibox, amrex::Box ba_box) { m_ba.set(ibox, ba_box); }
    void SetFabName (int ifab, std::string fodPrefix, std::string FabName,
                     int FabHead);
    void SetMinVal (int ifab, amrex::Vector<amrex::Real> minval);
    void SetMaxVal (int ifab, amrex::Vector<amrex::Real> maxval);
    private:
    std::string m_Header_path;
    int m_vers; 
    int m_how;
    int m_ncomp;
    int m_ngrow; 
    int m_ba_size;
    amrex::BoxArray m_ba;
    amrex::Vector<std::string> m_FabOnDiskPrefix;
    amrex::Vector<std::string> m_fabname;
    amrex::Vector<int> m_fabhead;
    /** The min of each component of each FAB in the BoxArray, m_ba.
     *  To access the min value of the ith fab and jth component [ifab][jcomp]*/
    amrex::Vector<amrex::Vector< amrex::Real> > m_minval;
    /** The max of each component of each FAB in the BoxArray, m_ba.
     *  To access the max value of the ith fab and jth component [ifab][jcomp]*/
    amrex::Vector<amrex::Vector< amrex::Real> > m_maxval;
    void CopyVec (amrex::Vector<amrex::Real>& dst,
                  amrex::Vector<amrex::Real> src);
};

#endif

