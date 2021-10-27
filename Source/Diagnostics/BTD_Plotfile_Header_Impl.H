#ifndef BTD_PLOTFILE_HEADER_IMPL_H
#define BTD_PLOTFILE_HEADER_IMPL_H

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

#include <string>

/**
 * \brief Class to read, modify, and write plotfile header when
 * back-transformed diag format is selected as plotfile.
 * This class enables multiple BTD buffers to be interweaved and stitched into
 * a single plotfile with a single Header.
 */

class BTDPlotfileHeaderImpl
{
public:
    /** Constructor.
     * \param[in] Headerfile_path string containing path of Headerfile
     */
    BTDPlotfileHeaderImpl (std::string const& Headerfile_path);
    /** Destructor */
    ~BTDPlotfileHeaderImpl () = default;

    /** Returns the Header file version for plotfile */
    std::string fileVersion () const noexcept {return m_file_version; }
    /** Returns the number of components written in the Headerfile */
    int ncomp () const noexcept {return m_nComp; }
    /** Returns the names of components in the Headerfile */
    const amrex::Vector<std::string>& varnames () const noexcept {return m_varnames; }
    /** Returns the number of dimensions in the Headerfile */
    int spaceDim () const noexcept {return m_spacedim; }
    /** Returns the physical time in the simulation in the boosted-frame */
    amrex::Real time () const noexcept {return m_time; }
    /** Returns finest level output in the Headerfile */
    int finestLevel () const noexcept { return m_finest_level; }
    /** Returns the physical co-ordinates of the lower corner in dimension, idim,
     *  that corresponds the to the respective plotfile data
     */
    amrex::Real problo (int dim) const noexcept {return m_prob_lo[dim]; }
    /** Returns the physical co-ordinates of the upper corner in dimension, idim,
     *  that corresponds the to the respective plotfile data.
     */
    amrex::Real probhi (int dim) const noexcept {return m_prob_hi[dim]; }

    /** Returns the bounding box of the domain spanned in the plotfile */
    amrex::Box probDomain () const noexcept {return m_prob_domain; }
    /** Returns timestep at which the plotfile was written */
    int timestep () const noexcept {return m_timestep; }
    int numFabs () const noexcept {return m_numFabs; }

    /** Returns physical co-ordinates of the lower-corner for the i-th Fab.
     *  \param[in] iFab id of the ith Fab in the list of Multifabs
     *  \returns   Array of lower-corner physical co-ordinates corresponding to the ith Fab
     */
    amrex::Array<amrex::Real, AMREX_SPACEDIM> FabLo (int iFab) const noexcept {return m_glo[iFab]; }
    /** Returns physical co-ordinates of the lower-corner for the i-th Fab.
     *  \param[in] iFab id of the ith Fab in the list of Multifabs
     *  \returns   Array of lower-corner physical co-ordinates corresponding to the ith Fab
     */
    amrex::Array<amrex::Real, AMREX_SPACEDIM> FabHi (int iFab) const noexcept {return m_ghi[iFab]; }
    /** Returns path to location of multifabs */
    std::string CellPath () const noexcept {return m_CellPath; }

    /** Reads the Header file data for BTD */
    void ReadHeaderData ();
    /** Writes Header file data for BTD */
    void WriteHeader ();

    /** Sets the physical simulation time, m_time, in the Header file to a new_time.
     *  \param[in] new_time current time in the boosted-frame read from Header plotfile
     **/
    void set_time ( amrex::Real new_time) { m_time = new_time;}
    /** Sets the timestep, m_timestep, in the Header file to a new_timestep
     *  \param[in] new_timestep current timestep in the boosted-frame (read from the Header file of the newly flushed buffer)
     **/
    void set_timestep (int new_timestep) { m_timestep = new_timestep; }
    /** Set the ith-dimension physical co-ordinate of the lower corner.
     *  \param[in] dim dimension modified.
     *  \param[in] lo  lower-corner physical coordinate to be stored in dimension, dim.
     **/
    void set_problo (int dim, amrex::Real lo) { m_prob_lo[dim] = lo;}
    /** Set the ith-dimension physical co-ordinate of the upper corner.
     *  \param[in] dim dimension modified.
     *  \param[in] hi  upper-corner physical coordinate to be stored in dimension, dim.
     **/
    void set_probhi (int dim, amrex::Real hi) { m_prob_hi[dim] = hi;}
    /** Set the Domain box spanning the plotfile data in the modified plotfile and Header.
     * \param[in] domainBox spanning the domain corresponding to the plotfile.
     */
    void set_probDomain (amrex::Box domainBox) {m_prob_domain = domainBox; }
    /** Increments the number of fabs stored in the plotfile by one. */
    void IncrementNumFabs () { m_numFabs++;}
    void ResizeFabLo () { m_glo.resize(m_numFabs); }
    void ResizeFabHi () { m_ghi.resize(m_numFabs); }
    /** Append array of lower-corner physical coordinates corresponding to a new fab to
     *  the existing list of coordinates, m_glo.
     *  \param[in] newFabLo contains physical coordinates
     *                      of the newly appended fab-data to the plotfile.
     */
    void AppendNewFabLo (amrex::Array<amrex::Real, AMREX_SPACEDIM> newFabLo);
    /** Append array of upper-corner physical coordinates corresponding to a new fab to
     *  the existing list of coordinates, m_ghi.
     *  \param[in] newFabHi contains physical coordinates
     *                      of the newly appended fab-data to the plotfile.
     */
    void AppendNewFabHi (amrex::Array<amrex::Real, AMREX_SPACEDIM> newFabHi);
private:
    /** string containing path of the Header file. */
    std::string m_Header_path;
    /** String containing file version of the plotfile. */
    std::string m_file_version;
    /** Number of components in the plotfile. */
    int m_nComp;
    /** Names of components stored in the plotfile. */
    amrex::Vector<std::string> m_varnames;
    /** Number of dimensions stored in the plotfile, should be same as AMREX_SPACEDIM */
    int m_spacedim;
    /** Physical time */
    amrex::Real m_time;
    int m_finest_level, m_nlevel;
    /** Lower cornder physical coordinates of the domain spanned in the plotfile */
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_prob_lo {{AMREX_D_DECL(0.,0.,0.)}};
    /** Upper corner physical coordinates of the domain spanned in the plotfile */
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_prob_hi {{AMREX_D_DECL(1.,1.,1.)}};
    /** Cell size */
    amrex::Array<amrex::Real, AMREX_SPACEDIM> m_cell_size {{AMREX_D_DECL(1.,1.,1.)}};
    /** Box covering the span of the physical domain in the plotfile */
    amrex::Box m_prob_domain;
    int m_timestep;
    int m_coordsys;
    int m_bwidth;
    int m_cur_level;
    /** Number of Fabs in the plotfile. */
    int m_numFabs;
    /** Lower corner physical coordinates of each fab in the plotfile. */
    amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM> > m_glo {{AMREX_D_DECL(0.,0.,0.)}};
    /** Upper corner physical coordinates of each fab in the plotfile. */
    amrex::Vector<amrex::Array<amrex::Real, AMREX_SPACEDIM> > m_ghi {{AMREX_D_DECL(1.,1.,1.)}};
    std::string m_CellPath;

};

/**
 * \brief Class to read, modify, and write MultiFab header in Level_0/Cell_H when
 * back-transformed diag format is selected as plotfile.
 * This class enables multiple fabs to be interweaved and stitched into
 * a single plotfile with a single Header, Cell_H.
 */
class BTDMultiFabHeaderImpl
{
    public:
    /** Constructor.
     * \param[in] Headerfile_path string containing path of Headerfile
     */
    BTDMultiFabHeaderImpl (std::string const& Headerfile_path);
    ~BTDMultiFabHeaderImpl () = default;
    /** Reads the Multifab Header file and stores its data. */
    void ReadMultiFabHeader ();
    /** Writes the meta-data of the Multifab in a header file, with path, m_Header_path. */
    void WriteMultiFabHeader ();

    /** Returns size, m_ba_size, of the Box Array, m_ba.*/
    int ba_size () {return m_ba_size;}
    /** Returns box corresponding to the ith box in the BoxArray, m_ba.
     *  \param[in] ibox index of the box in the BoxArray.
     */
    amrex::Box ba_box (int ibox) {return m_ba[ibox]; }
    /** Returns prefix of the ith-fab on disk, i.e., ith fab of the MultiFab data.
     *  \param[in] ifab index of the ith fab in the MultiFab data.
     */
    std::string fodPrefix (int ifab) {return m_FabOnDiskPrefix[ifab]; }
    /** Returns name of the ith fab stored in the MultiFab. */
    std::string FabName (int ifab) {return m_fabname[ifab]; }
    /** Returns the starting byte of the ith fab data */
    int FabHead (int ifab) {return m_fabhead[ifab]; }
    /** Returns minimum value of all the components stored in the ith fab. */
    amrex::Vector<amrex::Real> minval(int ifab) { return m_minval[ifab];}
    /** Returns maximum value of all the components stored in the ith fab. */
    amrex::Vector<amrex::Real> maxval(int ifab) { return m_maxval[ifab];}
    void ResizeFabData ();
    /** Increments MultiFab size, m_ba_size, by add_mf_size.
     *  \param[in] add_mf_size number of new multifabs to be appended to the existing
     *             Box Array.
     */
    void IncreaseMultiFabSize (int add_mf_size) {m_ba_size += add_mf_size;}
    /** Set Box indices of the ith-box in Box Array, m_ba, to the new Box, ba_box.
     *  \param[in] ibox   index of the ith box in BoxArray, m_ba.
     *  \param[in] ba_box dimensions corresponding to the ith Fab.
     */
    void SetBox (int ibox, amrex::Box ba_box) { m_ba.set(ibox, ba_box); }
    /** Set Fab name of the ith fab to be written in the multifab Header file.*/
    void SetFabName (int ifab, std::string fodPrefix, std::string FabName,
                     int FabHead);
    /** Set minimum value of all the components for the ith fab. */
    void SetMinVal (int ifab, amrex::Vector<amrex::Real> minval);
    /** Set maximum value of all the components for the ith fab. */
    void SetMaxVal (int ifab, amrex::Vector<amrex::Real> maxval);
    private:
    /** Header file path */
    std::string m_Header_path;
    int m_vers;
    int m_how;
    /** number of components stored in the multifab. */
    int m_ncomp;
    /** number of guard cells in the multifab. */
    int m_ngrow;
    /** Size of the BoxArray, m_ba.*/
    int m_ba_size;
    /** BoxArray corresponding to the multifab stored in the plotfile.*/
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
    /** Copy values of the vector from the src vector, src, to dst vector. */
    void CopyVec (amrex::Vector<amrex::Real>& dst,
                  amrex::Vector<amrex::Real> src);
};

#endif

