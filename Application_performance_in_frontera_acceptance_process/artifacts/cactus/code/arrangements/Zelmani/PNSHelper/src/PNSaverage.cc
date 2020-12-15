#include <stdio.h>
#include <cassert>
#include <cmath>
#include <mpi.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>
#include <cctk_Functions.h>
#include <cctk_Faces.h>
#include <Symmetry.h>
#include <loopcontrol.h>
#include <time.h>
#include <PNSinterp.hh>


using namespace PNSinterp;

extern "C" {
  void PNSHelper_AverageProfiles(CCTK_ARGUMENTS);
}


void PNSHelper_AverageProfiles(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(!*have_interp_data) return;
  
  double theta1 = theta[0];
  double theta2 = theta[ntheta-1];
  double phi1 = phi[0];
  double phi2 = phi[nphi-1];
   
  // average density
  pns_ang_ave(theta1,
	      theta2,
	      ntheta,
	      theta,
	      phi1,
	      phi2,
	      nphi,
	      phi,
	      nrad,
	      nrad_outer,
	      rad,
	      pns_rho,
	      pns_av_rho)  ;

  // average pressure
  pns_ang_ave(theta1,
	      theta2,
	      ntheta,
	      theta,
	      phi1,
	      phi2,
	      nphi,
	      phi,
	      nrad,
	      nrad_outer,
	      rad,
	      pns_press,
	      pns_av_press)  ;

  // average eps
  pns_ang_ave(theta1,
	      theta2,
	      ntheta,
	      theta,
	      phi1,
	      phi2,
	      nphi,
	      phi,
	      nrad,
	      nrad_outer,
	      rad,
	      pns_eps,
	      pns_av_eps)  ;

  // average Y_e
  pns_ang_ave(theta1,
	      theta2,
	      ntheta,
	      theta,
	      phi1,
	      phi2,
	      nphi,
	      phi,
	      nrad,
	      nrad_outer,
	      rad,
	      pns_ye,
	      pns_av_ye)  ;

  // average temperature
  pns_ang_ave(theta1,
	      theta2,
	      ntheta,
	      theta,
	      phi1,
	      phi2,
	      nphi,
	      phi,
	      nrad,
	      nrad_outer,
	      rad,
	      pns_temp,
	      pns_av_temp)  ;

  // average w_lorentz
  pns_ang_ave(theta1,
	      theta2,
	      ntheta,
	      theta,
	      phi1,
	      phi2,
	      nphi,
	      phi,
	      nrad,
	      nrad_outer,
	      rad,
	      pns_w_lorentz,
	      pns_av_w_lorentz)  ;

  // average old lapse
  pns_ang_ave(theta1,
	      theta2,
	      ntheta,
	      theta,
	      phi1,
	      phi2,
	      nphi,
	      phi,
	      nrad,
	      nrad_outer,
	      rad,
	      pns_old_alp,
	      pns_av_alp)  ;

  // average gxx
  pns_ang_ave(theta1,
	      theta2,
	      ntheta,
	      theta,
	      phi1,
	      phi2,
	      nphi,
	      phi,
	      nrad,
	      nrad_outer,
	      rad,
	      pns_old_gxx,
	      pns_av_gxx)  ;
 

}
