#ifndef WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXPOLYGONIF_H
#define WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXPOLYGONIF_H

#include "WarpX.H"

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_EB2_IF_Base.H>

#include <cmath>
#include <algorithm>
#include <regex>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>

typedef boost::geometry::model::d2::point_xy<amrex::Real> point_t;
typedef boost::geometry::model::polygon<point_t> polygon_t;
typedef boost::geometry::model::multi_polygon<polygon_t> mpolygon_t;

namespace
{
    // For all implicit functions, >0: body; =0: boundary; <0: fluid

    inline void add_jump(std::vector<amrex::Real> &r, std::vector<amrex::Real> &z, std::vector<size_t> &jumps)
    {
        if (jumps.size() > 0 && jumps[jumps.size() - 1] == r.size() - 1)
        {
            r.pop_back();
            z.pop_back();
        }
        else
        {
            jumps.push_back(r.size());
        }
    } 
    
    template <typename T>
    inline void expand_ring(const T &polygon, std::vector<amrex::Real> &r, std::vector<amrex::Real> &z, std::vector<size_t> &jumps)
    {
        for (auto p : polygon)
        {
            if (p.x() == 0 && r[r.size() - 1] == 0)
                add_jump(r, z, jumps);
            r.push_back(p.x());
            z.push_back(p.y());
        }
        add_jump(r, z, jumps);
    }
    
    void expand_multipolygon(const mpolygon_t &multipolygon, std::vector<amrex::Real> &r, std::vector<amrex::Real> &z, std::vector<size_t> &jumps)
    {
        for (auto polygon : multipolygon)
        {
            expand_ring(polygon.outer(), r, z, jumps);
            for (auto ring : polygon.inners())
                expand_ring(ring, r, z, jumps);
        }
    }

    void parse_multipolygon(const std::string &wkt, std::vector<amrex::Real> &r, std::vector<amrex::Real> &z, std::vector<size_t> &jumps) {
        mpolygon_t multipolygon;

        boost::geometry::read_wkt(wkt, multipolygon);
        
        expand_multipolygon(multipolygon, r, z, jumps);
    }

    class PolygonXYIF
        : public amrex::GPUable
    {
    public:
        PolygonXYIF(amrex::Real* r, amrex::Real *z, size_t n, size_t *jumps, size_t n_jumps)
            : m_r(r),
              m_z(z),
              m_n(n),
              m_jumps(jumps),
              m_n_jumps(n_jumps)
        {
        }

        ~PolygonXYIF() {}

        PolygonXYIF(const PolygonXYIF &rhs) noexcept = default;
        PolygonXYIF(PolygonXYIF &&rhs) noexcept = default;
        PolygonXYIF &operator=(const PolygonXYIF &rhs) = delete;
        PolygonXYIF &operator=(PolygonXYIF &&rhs) = delete;

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real
        operator()(AMREX_D_DECL(amrex::Real x, amrex::Real y, amrex::Real z))
            const noexcept
        {
            amrex::Real dist = 1;

            size_t i_jump = 0;

            for (size_t i = 0; i < m_n - 1; i++) {
                if(i_jump < m_n_jumps && m_jumps[i_jump] == i+1) {
                    i_jump++;
                    continue;
                }
                dist *= intersects_outwards(x, y, i);
            }

            return dist;
        }

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::Real
        operator()(const amrex::RealArray &p) const noexcept
        {
            return this->operator()(AMREX_D_DECL(p[0], p[1], p[2]));
        }

    protected:
        amrex::Real *m_r;
        amrex::Real *m_z;
        size_t m_n;

        size_t *m_jumps;
        size_t m_n_jumps;

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
            amrex::Real
            intersects_outwards(amrex::Real r, amrex::Real z, int i) const noexcept
        {
            // Determine if the point is radially inside the line segment from curve[lineseg_i] to curve[lineseg_i+1].

            if (m_r[i] <= r && m_r[i + 1] <= r)                                // point is completely outside the linesegment
                return 1;                                                      // not inside
            if (m_z[i] < z && m_z[i + 1] < z)                                  // point is completely above the linesegment
                return 1;                                                      // not inside
            if (m_z[i] >= z && m_z[i + 1] >= z)                                // point is completely below the linesegment
                return 1;                                                      // not inside
            if (m_r[i] > r && m_r[i + 1] > r && m_z[i] < z && z <= m_z[i + 1]) // point is completely inside the linesegment
                return -1;                                                     // inside

            // no clear bounding box decision can be made. Check for real if the point lies to the left or right of the line segment

            auto p_minus_a_r = r - m_r[i];
            auto p_minus_a_z = z - m_z[i];
            auto b_minus_a_r = m_r[i + 1] - m_r[i];
            auto b_minus_a_z = m_z[i + 1] - m_z[i];

            // check the sign of the third component of the cross product of (p-a) x (b-a)

            auto cross = (p_minus_a_r * b_minus_a_z) - (p_minus_a_z * b_minus_a_r);
            if (m_z[i] < m_z[i + 1])
                return cross;
            else
                return -cross;
        }
    };
}

#endif