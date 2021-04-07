
#include "WarpX.H"
#include <boost/math/special_functions/bessel.hpp>

amrex::RealArray
WarpX::AnalyticSolSphere (amrex::Real x, amrex::Real y, amrex::Real z, amrex::Real t, amrex::Real r_sphere) {
  amrex::Real k=2.7437/r_sphere;
  amrex::Real r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  amrex::Real theta = atan2(sqrt(pow(x, 2) + pow(y, 2)), z);
  amrex::Real phi = atan2(y, x);
  amrex::Real H_r = 0;
  amrex::Real H_theta = 0;
  amrex::Real mu_r = PhysConst::mu0;

  amrex::Real H_phi = k / (r_sphere * mu_r) * boost::math::sph_bessel(1, k*r) * sin(theta) * cos(k * PhysConst::c * t);

  amrex::Real H_x = H_r * sin(theta) * cos(phi) + H_theta * cos(theta) * cos(phi) - H_phi * sin(phi);
  amrex::Real H_y = H_r * cos(theta) * cos(phi) + H_theta * cos(theta) * sin(phi) + H_phi * cos(phi);
  amrex::Real H_z = H_r * cos(theta) - H_theta * sin(theta);

  return {PhysConst::mu0*H_x/1e6, PhysConst::mu0*H_y/1e6, PhysConst::mu0*H_z/1e6};
}

void
WarpX::InitEMFieldSphere (){
  	amrex::Real r_sphere;
  	amrex::ParmParse pp_eb2("eb2");
  	amrex::ParmParse pp_geometry("geometry");
  	pp_eb2.get("sphere_radius", r_sphere);
  	amrex::RealArray prob_lo;
  	pp_geometry.get("prob_lo", prob_lo);
  	auto &Bfield = Bfield_fp[maxLevel()];
  	auto &Efield = Efield_fp[maxLevel()];

  	auto &face_areas_0 = m_face_areas[maxLevel()];
  	for (amrex::MFIter mfi(*Bfield[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    	amrex::Array4<amrex::Real> const &Bx = Bfield[0]->array(mfi);
		amrex::Array4<amrex::Real> const &By = Bfield[1]->array(mfi);
		amrex::Array4<amrex::Real> const &Bz = Bfield[2]->array(mfi);
		amrex::Array4<amrex::Real> const &Sx = face_areas_0[0]->array(mfi);
		amrex::Array4<amrex::Real> const &Sy = face_areas_0[1]->array(mfi);
		amrex::Array4<amrex::Real> const &Sz = face_areas_0[2]->array(mfi);
		amrex::Box const &box = mfi.validbox();
		const auto hi = amrex::ubound(box);
		amrex::Real x, y, z;
		for (int i = 0; i < hi.x; i++) {
	  		for (int j = 0; j < hi.y; j++) {
				for (int k = 0; k < hi.z; k++) {
		  			if (Sx(i, j, k) > 0) {
						x = (i + 0.5) * CellSize(maxLevel())[0] + prob_lo[0];
						y = (j + 0.5) * CellSize(maxLevel())[1] + prob_lo[1];
						z = k * CellSize(maxLevel())[2] + prob_lo[2];

						Bx(i, j, k) = AnalyticSolSphere(x, y, z, 0., r_sphere)[0];
		  			}

		  			if (Sy(i, j, k) > 0) {
						x = (i + 0.5) * CellSize(maxLevel())[0] + prob_lo[0];
						y = j * CellSize(maxLevel())[1] + prob_lo[1];
						z = (k + 0.5) * CellSize(maxLevel())[2] + prob_lo[2];

						By(i, j, k) = AnalyticSolSphere(x, y, z, 0., r_sphere)[1];
		  			}

		  			if (Sz(i, j, k) > 0){
						x = i * CellSize(maxLevel())[0] + prob_lo[0];
						y = (j + 0.5) * CellSize(maxLevel())[1] + prob_lo[1];
						z = (k + 0.5) * CellSize(maxLevel())[2] + prob_lo[2];

						Bz(i, j, k) = AnalyticSolSphere(x, y, z, 0., r_sphere)[2];
		  			}
				}
	  		}
		}
  	}
}