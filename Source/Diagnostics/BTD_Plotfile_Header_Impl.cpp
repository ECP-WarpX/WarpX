#include "BTD_Plotfile_Header_Impl.H"
#include "WarpX.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FileSystem.H>
#include <AMReX_INT.H>
#include <memory>

using namespace amrex::literals;


BTDPlotfileHeaderImpl::BTDPlotfileHeaderImpl (std::string const & Headerfile_path)
    : m_Header_path(Headerfile_path)
{

}

void
BTDPlotfileHeaderImpl::ReadHeaderData ()
{

    // Read existing snapshot Header first
    amrex::Vector<char> HeaderCharPtr;
    amrex::Long fileLength(0), fileLengthPadded(0);
    std::ifstream iss;

    iss.open(m_Header_path.c_str(), std::ios::in);
    iss.seekg(0, std::ios::end);
    fileLength = static_cast<std::streamoff>(iss.tellg());
    iss.seekg(0, std::ios::beg);

    fileLengthPadded = fileLength + 1;
    HeaderCharPtr.resize(fileLengthPadded);
    iss.read(HeaderCharPtr.dataPtr(), fileLength);
    iss.close();
    HeaderCharPtr[fileLength] = '\0';

    std::istringstream is(HeaderCharPtr.dataPtr(), std::istringstream::in);

    is >> m_file_version;

    is >> m_nComp;
    m_varnames.resize(m_nComp);
    for (int icomp=0; icomp < m_nComp; ++icomp) {
        is >> m_varnames[icomp];
    }
    is >> m_spacedim ;
    is >> m_time;
    is >> m_finest_level;
    m_nlevel = m_finest_level + 1;
    for (int idim = 0; idim < m_spacedim; ++idim) {
        is >> m_prob_lo[idim];
    }
    for (int idim = 0; idim < m_spacedim; ++idim) {
        is >> m_prob_hi[idim];
    }
    WarpX::GotoNextLine(is);

    is >> m_prob_domain;

    is >> m_timestep;
    for (int idim = 0; idim < m_spacedim; ++idim) {
        is >> m_cell_size[idim];
    }
    is >> m_coordsys;
    is >> m_bwidth;
    is >> m_cur_level >> m_numFabs >> m_time;
    is >> m_timestep;
    m_glo.resize(m_numFabs);
    m_ghi.resize(m_numFabs);
    for (int igrid = 0; igrid < m_numFabs; ++igrid) {
        for (int idim = 0; idim < m_spacedim; ++idim ) {
            is >> m_glo[igrid][idim] >> m_ghi[igrid][idim];
        }
    }
    is >> m_CellPath;
}

void
BTDPlotfileHeaderImpl::AppendNewFabLo (amrex::Array<amrex::Real, AMREX_SPACEDIM> newFabLo)
{
    ResizeFabLo();
    m_glo[m_numFabs-1] = newFabLo;
}

void
BTDPlotfileHeaderImpl::AppendNewFabHi (amrex::Array<amrex::Real, AMREX_SPACEDIM> newFabHi)
{
    ResizeFabHi();
    m_ghi[m_numFabs-1] = newFabHi;
}

void
BTDPlotfileHeaderImpl::WriteHeader ()
{
    if ( amrex::FileExists(m_Header_path) ) {
        amrex::Print() << " removing this file : " << m_Header_path << "\n";
        amrex::FileSystem::Remove(m_Header_path);
    }
    std::ofstream HeaderFile;
    HeaderFile.open(m_Header_path.c_str(), std::ofstream::out |
                                           std::ofstream::trunc |
                                           std::ofstream::binary);
    if ( !HeaderFile.good()) amrex::FileOpenFailed(m_Header_path);

    HeaderFile.precision(17);

    // Generic Plotfile type name
    HeaderFile << m_file_version << '\n';
    // number of components
    HeaderFile << m_varnames.size() << '\n';
    // write the component string
    for (int icomp = 0; icomp < m_varnames.size(); ++icomp ) {
        HeaderFile << m_varnames[icomp] << '\n';
    }
    // space dim
    HeaderFile << m_spacedim << '\n';
    // time
    HeaderFile << m_time << '\n';
    // finest level
    HeaderFile << m_finest_level << '\n';
    // Physical coordinate of the lower corner
    for (int idim = 0; idim < m_spacedim; ++idim) {
        HeaderFile << m_prob_lo[idim] << ' ';
    }
    HeaderFile << '\n';
    // Physical cooridnate of the upper corner
    for (int idim = 0; idim < m_spacedim; ++idim) {
        HeaderFile << m_prob_hi[idim] << ' ';
    }
    HeaderFile << '\n';
    // since nlevels=0, not writing ref_ratio as seen in a typical MR Header
    HeaderFile << '\n';
    // Indices of the box covering the entire domain
    HeaderFile << m_prob_domain << '\n';
    // timestep in boosted frame
    HeaderFile << m_timestep << '\n';
    // cell-size in the back-transformed lab-frame
    for (int idim = 0; idim < m_spacedim; ++idim) {
        HeaderFile << m_cell_size[idim] << ' ';
    }
    HeaderFile << '\n';
    // coordinate system (Cartesian)
    HeaderFile << m_coordsys << '\n';
    //
    HeaderFile << m_bwidth << '\n';
    // current level, number of Fabs, current time -- for a single level (m_level = 0)
    HeaderFile << m_cur_level << ' ' << m_numFabs << ' ' << m_time << '\n';
    // timestep for level=0
    HeaderFile << m_timestep << '\n';
    // Physical (lo,hi) in each dimension for all the Fabs in the snapshot
    for (int iFab = 0; iFab < m_numFabs; ++iFab) {
        for (int idim = 0; idim < m_spacedim; ++idim) {
            HeaderFile << m_glo[iFab][idim] << ' ' << m_ghi[iFab][idim] << '\n';
        }
    }
    // MultiFabHeaderPath
    HeaderFile << m_CellPath << '\n';

}


BTDMultiFabHeaderImpl::BTDMultiFabHeaderImpl (std::string const & Headerfile_path)
    : m_Header_path(Headerfile_path)
{

}

void
BTDMultiFabHeaderImpl::ReadMultiFabHeader ()
{
    // Read existing fab Header first
    amrex::Vector<char> HeaderCharPtr;
    amrex::Long fileLength(0), fileLengthPadded(0);
    std::ifstream iss;

    iss.open(m_Header_path.c_str(), std::ios::in);
    iss.seekg(0, std::ios::end);
    fileLength = static_cast<std::streamoff>(iss.tellg());
    iss.seekg(0, std::ios::beg);

    fileLengthPadded = fileLength + 1;
    HeaderCharPtr.resize(fileLengthPadded);
    iss.read(HeaderCharPtr.dataPtr(), fileLength);
    iss.close();
    HeaderCharPtr[fileLength] = '\0';

    std::istringstream is(HeaderCharPtr.dataPtr(), std::istringstream::in);

    is >> m_vers;
    is >> m_how;
    is >> m_ncomp;
    is >> m_ngrow;
    // can also call readBoxArray(m_ba, is, True);
    int in_hash;
    int bl_ignore_max = 100000;
    is.ignore(bl_ignore_max,'(') >> m_ba_size >> in_hash;
    m_ba.resize(m_ba_size);
    for (int ibox = 0; ibox < m_ba.size(); ++ibox) {
        amrex::Box bx;
        is >> bx;
        m_ba.set(ibox, bx);
    }
    is.ignore(bl_ignore_max, ')');

    is >> in_hash; // repeat of reading ba_size
    m_FabOnDiskPrefix.resize(m_ba.size());
    m_fabname.resize(m_ba.size());
    m_fabhead.resize(m_ba.size());
    for (int ifab = 0; ifab < m_ba.size(); ++ifab) {
        is >> m_FabOnDiskPrefix[ifab] >> m_fabname[ifab] >> m_fabhead[ifab];
    }
    WarpX::GotoNextLine(is);
    char ch;
    is >> in_hash >> ch >> in_hash;
    m_minval.resize(m_ba.size());
    for (int ifab = 0; ifab < m_ba.size(); ++ifab) {
        m_minval[ifab].resize(m_ncomp);
        for (int icomp = 0; icomp < m_ncomp; ++icomp) {
            is >> m_minval[ifab][icomp] >> ch;
            if( ch != ',' ) amrex::Error("Expected a ',' got something else");
        }
    }
    WarpX::GotoNextLine(is);
    is >> in_hash >> ch >> in_hash;
    m_maxval.resize(m_ba.size());
    for (int ifab = 0; ifab < m_ba.size(); ++ifab) {
        m_maxval[ifab].resize(m_ncomp);
        for (int icomp = 0; icomp < m_ncomp; ++icomp) {
            is >> m_maxval[ifab][icomp] >> ch;
            if( ch != ',' ) amrex::Error("Expected a ',' got something else");
        }
    }


}

void
BTDMultiFabHeaderImpl::WriteMultiFabHeader ()
{
    if ( amrex::FileExists(m_Header_path) ) {
        amrex::Print() << " removing this file : " << m_Header_path << "\n";
        amrex::FileSystem::Remove(m_Header_path);
    }
    std::ofstream FabHeaderFile;
    FabHeaderFile.open(m_Header_path.c_str(), std::ofstream::out |
                                           std::ofstream::trunc |
                                           std::ofstream::binary);
    if ( !FabHeaderFile.good()) amrex::FileOpenFailed(m_Header_path);

    FabHeaderFile.precision(17);

    // multifab header version
    FabHeaderFile << m_vers << '\n';
    // VisMF :: how
    FabHeaderFile << m_how << '\n';
    // number of components
    FabHeaderFile << m_ncomp << '\n';
    // number of guard cells
    FabHeaderFile << m_ngrow << '\n';
    // WriteBoxArray
    m_ba.writeOn(FabHeaderFile);
    FabHeaderFile << '\n';
    // Write ba size
    FabHeaderFile << m_ba.size() << '\n';
    for (int ifab = 0; ifab < m_ba.size(); ++ifab) {
        FabHeaderFile << m_FabOnDiskPrefix[ifab] << ' ' << m_fabname[ifab] << ' ' << m_fabhead[ifab];
        FabHeaderFile << '\n';
    }
    FabHeaderFile << '\n';
    // Write minvalue of all the components for each fab in the multifab
    FabHeaderFile << m_ba.size() << ',' << m_ncomp << '\n';
    for (int ifab = 0; ifab < m_ba.size(); ++ifab) {
        for (int icomp = 0; icomp < m_ncomp; ++icomp) {
            FabHeaderFile << m_minval[ifab][icomp] << ',';
        }
        FabHeaderFile << '\n';
    }
    FabHeaderFile << '\n';
    // Write minvalue of all the components for each fab in the multifab
    FabHeaderFile << m_ba.size() << ',' << m_ncomp << '\n';
    for (int ifab = 0; ifab < m_ba.size(); ++ifab) {
        for (int icomp = 0; icomp < m_ncomp; ++icomp) {
            FabHeaderFile << m_maxval[ifab][icomp] << ',';
        }
        FabHeaderFile << '\n';
    }
}

void
BTDMultiFabHeaderImpl::ResizeFabData ()
{
    m_ba.resize(m_ba_size);
    m_FabOnDiskPrefix.resize(m_ba_size);
    m_fabname.resize(m_ba_size);
    m_fabhead.resize(m_ba_size);
    m_minval.resize(m_ba_size);
    m_maxval.resize(m_ba_size);
}

void
BTDMultiFabHeaderImpl::SetFabName (int ifab, std::string fodPrefix, std::string FabName,
                                   int FabHead)
{
    m_FabOnDiskPrefix[ifab] = fodPrefix;
    m_fabname[ifab] = FabName;
    m_fabhead[ifab] = FabHead;

}

void
BTDMultiFabHeaderImpl::SetMinVal (int ifab, amrex::Vector<amrex::Real> minval)
{
    CopyVec(m_minval[ifab], minval);
}

void
BTDMultiFabHeaderImpl::SetMaxVal (int ifab, amrex::Vector<amrex::Real> maxval)
{
    CopyVec(m_maxval[ifab], maxval);
}

void
BTDMultiFabHeaderImpl::CopyVec(amrex::Vector<amrex::Real>& dst,
                               amrex::Vector<amrex::Real> src)
{
    dst.resize(src.size());
    for (int i = 0; i < src.size(); ++i) {
        dst[i] = src[i];
    }
}


