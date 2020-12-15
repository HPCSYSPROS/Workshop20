#ifndef ZST_HH
#define ZST_HH

#define PRESS_GF 1.7982953469278e-39
#define RHO_GF 1.61620075314614e-18
#define EPS_GF 1.11265006e-21
#define LENGTH_GF 6.77140812e-06
#define INV_RHO_GF 6.18735016707159e17
#define INV_PRESS_GF 5.5608218177753538
#define INV_LENGTH_GF 147679.77092481e0

#include <cmath>
#include <cctk_Parameters.h>

namespace ZST {
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
  
} // namespace ZST

#endif // ZST_HH
