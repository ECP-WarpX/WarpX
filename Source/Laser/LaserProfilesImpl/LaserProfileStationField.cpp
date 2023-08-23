#include "Laser/LaserProfiles.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"

#include <AMReX_VisMF.H>

namespace WarpXLaserProfiles {

void StationFieldLaserProfile::init (const amrex::ParmParse& ppl,
                                     CommonLaserParameters params)
{
    m_common_params = std::move(params);
    if (m_common_params.p_X[0] == amrex::Real(1.0)) {
        m_xy = 0;
    } else {
        m_xy = 1;
    }

    ppl.get("station_file_name", m_station_file);

    m_ibuffer = -1;

    amrex::Vector<char> headerfile;
    amrex::ParallelDescriptor::ReadAndBcastFile(m_station_file+"/StationHeader",
                                                headerfile);
    std::istringstream is(std::string(headerfile.data()));
    amrex::Real tmp;
    is >> tmp;
    std::string tt;
    std::getline(is,tt); // eat \n
    while (std::getline(is, tt)) {
        std::istringstream istt(tt);
        amrex::Real t0, t1;
        istt >> t0;
        istt >> t1;
        m_times.emplace_back(t0,t1);
    }
}

void StationFieldLaserProfile::update (amrex::Real t)
{
    const auto nbuffers = int(m_times.size());

    bool data_is_ready = (m_ibuffer >= 0 && m_ibuffer < nbuffers &&
                          t >= m_times[m_ibuffer].first &&
                          t <  m_times[m_ibuffer].second);

    if (! data_is_ready) {
        for (m_ibuffer = 0; m_ibuffer < nbuffers; ++m_ibuffer) {
            if (t >= m_times[m_ibuffer].first &&
                t <  m_times[m_ibuffer].second) {
                break;
            }
        }

        if (m_ibuffer >= 0 && m_ibuffer < nbuffers) {
            amrex::MultiFab mf;
            amrex::VisMF::Read(mf, m_station_file+"/Level_0/buffer-"
                               +std::to_string(m_ibuffer));

            amrex::Print() << "xxxxx LaserProfileStationField: reading "
                           << m_station_file+"/Level_0/buffer-"+std::to_string(m_ibuffer)
                           << "\n";

            m_station_mf = std::make_unique<amrex::MultiFab>
                (mf.boxArray(), mf.DistributionMap(), 1, 0);

            auto const& sa = mf.const_arrays();
            auto const& da = m_station_mf->arrays();
            if (m_xy == 0) {
                amrex::ParallelFor(*m_station_mf,
                [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                {
                    // 0.5*(Ex+c*By)
                    da[bno](i,j,k) = amrex::Real(0.5) *
                        (sa[bno](i,j,k,0) + PhysConst::c*sa[bno](i,j,k,4));
                });
            } else {
                amrex::ParallelFor(*m_station_mf,
                [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
                {
                    // 0.5*(Ey-c*Bx)
                    da[bno](i,j,k) = amrex::Real(0.5) *
                        (sa[bno](i,j,k,1) - PhysConst::c*sa[bno](i,j,k,3));
                });
            }
        } else {
            m_station_mf.reset();
            m_slice_fab.reset();

            amrex::Print() << "xxxxx LaserProfileStationField: no more data\n";

            return;
        }
    }

    m_station_domain = m_station_mf->boxArray().minimalBox();
    auto const nz = m_station_domain.length(AMREX_SPACEDIM-1);
    AMREX_ASSERT(nz > 1);

    if (m_slice_fab == nullptr) {
        auto slicebox = m_station_domain;
        slicebox.setSmall(AMREX_SPACEDIM-1,0);
        if (amrex::ParallelDescriptor::MyProc() == 0) {
            slicebox.setBig(AMREX_SPACEDIM-1,1);
        } else {
            slicebox.setBig(AMREX_SPACEDIM-1,0);
        }
        m_slice_fab = std::make_unique<amrex::FArrayBox>(slicebox,1);
    }

    auto const dt = (m_times[m_ibuffer].second - m_times[m_ibuffer].first)
        / amrex::Real(nz-1);
    auto islice = static_cast<int>(std::floor((t-m_times[m_ibuffer].first)/dt));
    islice = std::min(islice, nz-2);

    {
        auto slicebox = m_station_domain;
        slicebox.setSmall(AMREX_SPACEDIM-1, islice);
        slicebox.setBig(AMREX_SPACEDIM-1, islice+1);
        amrex::BoxArray sliceba(slicebox);
        amrex::DistributionMapping slicedm(amrex::Vector<int>({0}));
        amrex::MultiFab slicemf(sliceba, slicedm, 1, 0, amrex::MFInfo().SetAlloc(false));
        if (amrex::ParallelDescriptor::MyProc() == 0) {
            slicemf.setFab(0, amrex::FArrayBox(slicebox, 1, m_slice_fab->dataPtr()));
        }
        slicemf.ParallelCopy(*m_station_mf);
    }

    if (amrex::ParallelDescriptor::MyProc() == 0) {
        amrex::Real w = (t - m_times[m_ibuffer].first)/dt - amrex::Real(islice);
        amrex::Box b = m_slice_fab->box();
        b.setBig(AMREX_SPACEDIM-1,0);
        auto const& a = m_slice_fab->array();
        amrex::ParallelFor(b, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
#if (AMREX_SPACEDIM == 2)
            a(i,j,0) = (a(i,0,0)*(amrex::Real(1.0)-w) + a(i,1,0)*w);
#elif (AMREX_SPACEDIM == 3)
            a(i,j,0) = (a(i,j,0)*(amrex::Real(1.0)-w) + a(i,j,1)*w);
#endif
        });
    }

    amrex::Box b = m_slice_fab->box();
    b.setBig(AMREX_SPACEDIM-1,0);
    amrex::ParallelDescriptor::Bcast(m_slice_fab->dataPtr(), b.numPts());
}

void StationFieldLaserProfile::fill_amplitude (
    const int np,
    amrex::Real const * AMREX_RESTRICT const Xp,
    amrex::Real const * AMREX_RESTRICT const Yp,
    amrex::Real t,
    amrex::Real * AMREX_RESTRICT const amplitude) const
{
    using namespace amrex::literals;
    if (m_ibuffer >= 0 && m_ibuffer < int(m_times.size())) {
        auto const& geom = WarpX::GetInstance().Geom(0);
        auto const plo = geom.ProbLoArray();
        auto const phi = geom.ProbHiArray();
        auto const dxi = amrex::Real(m_station_domain.length(0)-1) / (phi[0]-plo[0]); // -1: because box is nodal
#if (AMREX_SPACEDIM == 3)
        auto const dyi = amrex::Real(m_station_domain.length(1)-1) / (phi[1]-plo[1]);
#endif
        auto const& a = m_slice_fab->const_array();
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int ip)
        {
            auto const xi = (Xp[ip]-plo[0])*dxi;
            auto const i = static_cast<int>(std::floor(xi));
            amrex::Real const wx = xi - amrex::Real(i);
#if (AMREX_SPACEDIM == 2)
            amplitude[ip] = (1._rt-wx)*a(i  ,0,0)
                +                  wx *a(i+1,0,0);
#elif (AMREX_SPACEDIM == 3)
            auto const yj = (Yp[ip]-plo[1])*dyi;
            auto const j = static_cast<int>(std::floor(yj));
            amrex::Real const wy = yi - amrex::Real(j);
            amplitude[ip] = (1._rt-wx)*(1._rt-wy)*a(i  ,j  ,0)
                +                  wx *(1._rt-wy)*a(i+1,j  ,0)
                +           (1._rt-wx)*       wy *a(i  ,j+1,0)
                +                  wx *       wy *a(i+1,j+1,0);
#endif
        });
    } else {
        amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int ip)
        {
            amplitude[ip] = 0.0;
        });
    }
}

}
