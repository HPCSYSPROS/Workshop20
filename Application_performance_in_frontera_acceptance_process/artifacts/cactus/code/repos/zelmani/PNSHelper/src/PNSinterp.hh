#ifndef PNS_HH
#define PNS_HH

#define PRESS_GF 1.7982953469278e-39
#define RHO_GF 1.61620075314614e-18
#define EPS_GF 1.11265006e-21
#define LENGTH_GF 6.77140812e-06
#define INV_RHO_GF 6.18735016707159e17
#define INV_PRESS_GF 5.5608218177753538
#define INV_LENGTH_GF 147679.77092481e0

#include <cmath>
#include <cctk_Parameters.h>

namespace PNSinterp {
  using namespace std;
 
  inline
  void
  spher2cart (CCTK_REAL rad, CCTK_REAL theta, CCTK_REAL phi,
              CCTK_REAL* x, CCTK_REAL* y, CCTK_REAL* z,
	      CCTK_REAL* x0, CCTK_REAL* y0, CCTK_REAL* z0)
    
  {
    DECLARE_CCTK_PARAMETERS;
    using namespace std;
    
    *x = *x0 + rad*sin(theta)*cos(phi);
    *y = *y0 + rad*sin(theta)*sin(phi);
    *z = *z0 + rad*cos(theta);
  }
  
  
  
  inline
  void
  cart2spher (CCTK_REAL x, CCTK_REAL y, CCTK_REAL z,
              CCTK_REAL& rad, CCTK_REAL& theta, CCTK_REAL& phi,
	      CCTK_REAL& x0, CCTK_REAL& y0, CCTK_REAL& z0)
  {
    DECLARE_CCTK_PARAMETERS;
    
    x -= x0;
    y -= y0;
    z -= z0;
    
    rad   = sqrt (x*x + y*y + z*z);
    theta = acos (z / rad);
    phi   = atan2 (y, x);
  }
  
  inline
  void pns_ang_ave(double theta1,
		     double theta2,
		     int ntheta,
		     double *theta,
		     double phi1,
		     double phi2,
		     int nphi,
		     double *phi,
		     int nrad,
		     int nrad_outer,
		     double *rad,
		     double *data,
		     double *average)

  {

    for(int irad=0;irad<nrad+nrad_outer;irad++) {
      int nntheta = 0;  
      double tav = 0.0;
      for(int itheta=0;itheta<ntheta;itheta++) {
	if(theta[itheta]>=theta1 && theta[itheta]<=theta2) {
	  nntheta++;
	  double pav = 0.0;
	  int nnphi = 0;
	  for(int iphi=0;iphi<nphi;iphi++) {
	    if(phi[iphi]>=phi1 && phi[iphi]<=phi2) {
	      nnphi++;
	      int const iind3d = irad+(nrad+nrad_outer)*(itheta+ntheta*iphi);
	      pav += data[iind3d];
	    }
	  }
	  tav += pav/(nnphi+1.0e-10);
	}
      }
      tav = tav / (nntheta+1.0e-10);
      average[irad] = tav;
    }

  }


} // namespace PNSinterp

#endif // PNS_HH
