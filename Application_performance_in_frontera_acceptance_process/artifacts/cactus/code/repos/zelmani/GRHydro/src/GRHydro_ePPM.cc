#include <cmath>
#include <algorithm>
#include <vector>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MIN3(a,b,c) (MIN(a, MIN(b, c)))
#define MIN4(a,b,c,d) (MIN(a, MIN(b, MIN(c, d))))

#define MAX5(a,b,c,d,e) (MAX(a, MAX(b, MAX(c, MAX(d, e)))))

using namespace std;


/**
   ePPM Reconstruction (see Reisswig et al. 2013) based on McCorquodale & Colella (2011).
   Note: This routine differs from the Fortran routine (and the algorithm described in Reisswig et al. 2013)
         in the following way: We do not treat the specific internal energy "epsilon" in any special way.
         Rather, to avoid reconstruction to eps<0 for gamma-laws, we rely on
         'eps<0'-error checks after reconstruction (a similar treatment is also necessary for other reconstruction methods such as WENO5 and MP5).
*/


/*
  Cases that must be considered:
  * basic hydro
  * hydro + temperature + ye
  * hydro + ye
  * basic mhd
  * mhd + temperature + ye 
  * mhd + ye 
  * mppm (not supported right now)
  * not supporting trivial_rp
  * with or without divergence cleaning
 */


static inline double approx_at_cell_interface(const double* const restrict a, const int i) {
  return 7.0/12.0*(a[i]+a[i+1]) - 1.0/12.0*(a[i-1]+a[i+2]);
}

static inline double limit(const double* const restrict a, const double ah, const double C, const int i) 
{
  if ((MIN(a[i],a[i+1]) <= ah) && (ah <= MAX(a[i],a[i+1]))) {
     return ah;
  } else {
     const double D2a  = 3.0 * ((a[i]   + a[i+1]) - 2.0*ah    );
     const double D2aL =       ((a[i-1] + a[i+1]) - 2.0*a[i]  );
     const double D2aR =       ((a[i]   + a[i+2]) - 2.0*a[i+1]);
     const double D2aLim = copysign(1.0, D2a)*MIN3(C*fabs(D2aL), C*fabs(D2aR), fabs(D2a)) * 1.0/3.0;
     if (D2a*D2aR >= 0 && D2a*D2aL >= 0)
        return 0.5*(a[i]+a[i+1]) - D2aLim;
     else
        return 0.5*(a[i]+a[i+1]);
  }
  return 0;
}


/// Monotonicity of PPM 2011 (McCorquodale & Colella 2011), Sec. 2.4.1, Eq. 23-34.
/// This does not use the check for deviations from a cubic. Thus, it gets away with only 3 stencil points!
static inline void monotonize(double* const restrict aminus,
			      const double* const restrict a,
			      double* const restrict aplus,
                              const double C,
			      const int i) 
{
   double D2aLim = 0; 
   double rhi = 0;
   const double daplus  = aplus[i]-a[i];
   const double daminus = a[i]-aminus[i];
   if (daplus*daminus <= 0 || (a[i-2]-a[i])*(a[i]-a[i+2]) <= 0) {
      const double D2a  = - (12.0*a[i] - 6.0*(aminus[i]+aplus[i]));
      const double D2aC = (a[i-1] + a[i+1]) - 2.0*a[i];
      const double D2aL = (a[i-2] + a[i]  ) - 2.0*a[i-1];
      const double D2aR = (a[i]   + a[i+2]) - 2.0*a[i+1];
      if (copysign(1.0, D2a) == copysign(1.0, D2aC) && copysign(1.0, D2a) == copysign(1.0, D2aL) && copysign(1.0, D2a) == copysign(1.0, D2aR))
         D2aLim = copysign(1.0, D2a) * MIN4(C*fabs(D2aL), C*fabs(D2aR), C*fabs(D2aC), fabs(D2a));
      if (!(fabs(D2a) <= 1e-12*MAX5(fabs(a[i-2]), fabs(a[i-1]), fabs(a[i]), fabs(a[i+1]), fabs(a[i+2]))))
         rhi = D2aLim / D2a;
      if (! (rhi >= 1.0 - 1e-12)) {
            if (daplus*daminus < 0) {
               aplus[i]  = a[i] + daplus * rhi;
               aminus[i] = a[i] - daminus * rhi;
            } else if (fabs(daminus) >= 2.0*fabs(daplus)) {
               aminus[i]  = a[i] - (2.0*(1.0-rhi)*daplus + rhi*daminus);
            } else if (fabs(daplus) >= 2.0*fabs(daminus)) {
               aplus[i]  = a[i] + (2.0*(1.0-rhi)*daminus + rhi*daplus);
            }
      }
      //trivial_rp(i-1) = .false.
      //trivial_rp(i) = .false.
   } else {
         if (fabs(daplus) >= 2.0*fabs(daminus))
            aplus[i] = a[i] + 2.0*daminus;
         else if (fabs(daminus) >= 2.0*fabs(daplus))
            aminus[i] = a[i] - 2.0*daplus;
      //trivial_rp(i-1) = .false.
      //trivial_rp(i) = .false.
   }
   
   return;
}


template<bool do_temp, bool do_ye, bool do_mhd, 
	 bool dc_flag, bool do_ppm_detect>
void GRHydro_eppm1d_cxx(const int nx, 
		      const double dx, 
		      const double* const restrict rho, 
		      const double* const restrict velx,
		      const double* const restrict vely, 
		      const double* const restrict velz, 
		      const double* const restrict eps, 
		      const double* const restrict press,
		      const double* const restrict temp,
		      const double* const restrict ye,
		      const double* const restrict Bvcx,
		      const double* const restrict Bvcy,
		      const double* const restrict Bvcz,
		      const double* const restrict psidc,
		      double* const restrict rhominus, 
		      double* const restrict velxminus, 
		      double* const restrict velyminus,
		      double* const restrict velzminus, 
		      double* const restrict epsminus, 
		      double* const restrict tempminus, 
		      double* const restrict yeminus, 
		      double* const restrict Bvcxminus,
		      double* const restrict Bvcyminus,
		      double* const restrict Bvczminus,
		      double* const restrict psidcminus,
		      double* const restrict rhoplus, 
		      double* const restrict velxplus, 
		      double* const restrict velyplus,
		      double* const restrict velzplus, 
		      double* const restrict epsplus,
		      double* const restrict tempplus, 
		      double* const restrict yeplus,
		      double* const restrict Bvcxplus,
		      double* const restrict Bvcyplus,
		      double* const restrict Bvczplus,
		      double* const restrict psidcplus)
{
  DECLARE_CCTK_PARAMETERS;


  // We initialize "plus" \equiv a_j+1/2 with (16) via APPROX_AT_CELL_INTERFACE, 
  // then checking for (13) of Colella & Sekora 2008 and applying
  // (18) and (19) if (13) is not satisfied. This is done with LIMIT.
  for (int i=1; i < nx-2; ++i) {
         rhoplus[i] = approx_at_cell_interface(rho, i);
         rhoplus[i] = limit(rho, rhoplus[i], enhanced_ppm_C2, i);
         rhominus[i+1] = rhoplus[i];
  
         epsplus[i] = approx_at_cell_interface(eps, i);
         epsplus[i] = limit(eps, epsplus[i], enhanced_ppm_C2, i);
         epsminus[i+1] = epsplus[i];
         
         velxplus[i] = approx_at_cell_interface(velx, i);
         velxplus[i] = limit(velx, velxplus[i], enhanced_ppm_C2, i);
         velxminus[i+1] = velxplus[i];
         
         velyplus[i] = approx_at_cell_interface(vely, i);
         velyplus[i] = limit(vely, velyplus[i], enhanced_ppm_C2, i);
         velyminus[i+1] = velyplus[i];
         
         velzplus[i] = approx_at_cell_interface(velz, i);
         velzplus[i] = limit(velz, velzplus[i], enhanced_ppm_C2, i);
         velzminus[i+1] = velzplus[i];
         
         if(do_ye) {
            if(do_temp) {
              tempplus[i] = approx_at_cell_interface(temp, i);
              tempplus[i] = limit(temp, tempplus[i], enhanced_ppm_C2, i);
              tempminus[i+1] = tempplus[i];
            }
            yeplus[i] = approx_at_cell_interface(ye, i);
            yeplus[i] = limit(ye, yeplus[i], enhanced_ppm_C2, i);
            yeminus[i+1] = yeplus[i];
         }
         if(do_mhd) {
            Bvcxplus[i] = approx_at_cell_interface(Bvcx, i);
            Bvcxplus[i] = limit(Bvcx, Bvcxplus[i], enhanced_ppm_C2, i);
            Bvcxminus[i+1] = Bvcxplus[i];
            
            Bvcyplus[i] = approx_at_cell_interface(Bvcy, i);
            Bvcyplus[i] = limit(Bvcy, Bvcyplus[i], enhanced_ppm_C2, i);
            Bvcyminus[i+1] = Bvcyplus[i];
            
            Bvczplus[i] = approx_at_cell_interface(Bvcz, i);
            Bvczplus[i] = limit(Bvcz, Bvczplus[i], enhanced_ppm_C2, i);
            Bvczminus[i+1] = Bvczplus[i];
            
            if(dc_flag) {
               psidcplus[i] = approx_at_cell_interface(psidc, i);
               psidcplus[i] = limit(psidc, psidcplus[i], enhanced_ppm_C2, i);
               psidcminus[i+1] = psidcplus[i];
            }
         }
         
   }

  // Monotonize
  for(int i=GRHydro_stencil-1; i<nx-GRHydro_stencil+1; ++i) {
    monotonize(rhominus,rho,rhoplus,enhanced_ppm_C2,i);
    monotonize(epsminus,eps,epsplus,enhanced_ppm_C2,i);
    monotonize(velxminus,velx,velxplus,enhanced_ppm_C2,i);
    monotonize(velyminus,vely,velyplus,enhanced_ppm_C2,i);
    monotonize(velzminus,velz,velzplus,enhanced_ppm_C2,i);
    if(do_ye) {
      if(do_temp) {
	monotonize(tempminus,temp,tempplus,enhanced_ppm_C2,i);
      }
      monotonize(yeminus,ye,yeplus,enhanced_ppm_C2,i);
    }
    if(do_mhd) {
      monotonize(Bvcxminus,Bvcx,Bvcxplus,enhanced_ppm_C2,i);
      monotonize(Bvcyminus,Bvcy,Bvcyplus,enhanced_ppm_C2,i);
      monotonize(Bvczminus,Bvcz,Bvczplus,enhanced_ppm_C2,i);
      if(dc_flag) {
	monotonize(psidcminus,psidc,psidcplus,enhanced_ppm_C2,i);
      }
    }
  }

  
  // Finally, apply flattening
  vector<double> dpress(nx,0);
  for(int i=1; i<nx-1; ++i)
    dpress[i] = (press[i+1]-press[i-1]);
  
  vector<double> tilde_flatten(nx,0);
  for(int i=2; i<nx-2; ++i) {
    const double dpress2 = press[i+2] - press[i-2];
    const double dvel = velx[i+1] - velx[i-1];
    double w=0.0;
    if ( (fabs(dpress[i]) >  ppm_epsilon * MIN(press[i-1],press[i+1])) 
         && (dvel < 0.0) ) 
      {
        w = 1.0;
      } 
    if (fabs(dpress2) < ppm_small) 
      {
        tilde_flatten[i] = 1.0;
      } 
    else
      {
        tilde_flatten[i] = MAX(0.0, 1.0 - w * MAX(0.0, 
                           ppm_omega2 * (dpress[i] / dpress2 - ppm_omega1)));
      }
  } 
  
  for(int i=2; i<nx-2; ++i) {
    const double flatten = tilde_flatten[i];
    rhoplus[i] = flatten * rhoplus[i] + (1.0 - flatten) * rho[i];
    rhominus[i] = flatten * rhominus[i] + (1.0 - flatten) * rho[i];
    epsplus[i] = flatten * epsplus[i] + (1.0 - flatten) * eps[i];
    epsminus[i] = flatten * epsminus[i] + (1.0 - flatten) * eps[i];
    velxplus[i] = flatten * velxplus[i] + (1.0 - flatten) * velx[i];
    velxminus[i] = flatten * velxminus[i] + (1.0 - flatten) * velx[i];
    velyplus[i] = flatten * velyplus[i] + (1.0 - flatten) * vely[i];
    velyminus[i] = flatten * velyminus[i] + (1.0 - flatten) * vely[i];
    velzplus[i] = flatten * velzplus[i] + (1.0 - flatten) * velz[i];
    velzminus[i] = flatten * velzminus[i] + (1.0 - flatten) * velz[i];

    if(do_ye) {
      if(do_temp) {
        tempplus[i] = flatten * tempplus[i] + (1.0 - flatten) * temp[i];
        tempminus[i] = flatten * tempminus[i] + (1.0 - flatten) * temp[i];
      }
      yeplus[i] = flatten * yeplus[i] + (1.0 - flatten) * ye[i];
      yeminus[i] = flatten * yeminus[i] + (1.0 - flatten) * ye[i];
    }

    if(do_mhd) {
      Bvcxplus[i] = flatten * Bvcxplus[i] + (1.0 - flatten) * Bvcx[i];
      Bvcxminus[i] = flatten * Bvcxminus[i] + (1.0 - flatten) * Bvcx[i];
      Bvcyplus[i] = flatten * Bvcyplus[i] + (1.0 - flatten) * Bvcy[i];
      Bvcyminus[i] = flatten * Bvcyminus[i] + (1.0 - flatten) * Bvcy[i];
      Bvczplus[i] = flatten * Bvczplus[i] + (1.0 - flatten) * Bvcz[i];
      Bvczminus[i] = flatten * Bvczminus[i] + (1.0 - flatten) * Bvcz[i];
      if(dc_flag) {
        psidcplus[i] = flatten * psidcplus[i] + (1.0 - flatten) * psidc[i];
        psidcminus[i] = flatten * psidcminus[i] + (1.0 - flatten) * psidc[i];
      }
    }
  } // flattening
  
  return;
}

//template template<bool do_temp, bool do_ye, bool do_mhd, bool dc_flag, bool do_ppm_detect>
#define INSTANTIATE(do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect)			\
template void GRHydro_eppm1d_cxx<do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect>(const int nx, \
		      const double dx,					\
		      const double* restrict rho,			\
		      const double* restrict velx,			\
		      const double* restrict vely,			\
		      const double* restrict velz,			\
		      const double* restrict eps,			\
		      const double* restrict press,			\
		      const double* restrict temp,			\
		      const double* restrict ye,			\
		      const double* restrict Bvcx,			\
		      const double* restrict Bvcy,			\
		      const double* restrict Bvcz,			\
		      const double* restrict psidc,			\
		      double* restrict rhominus,			\
		      double* restrict velxminus,			\
		      double* restrict velyminus,			\
		      double* restrict velzminus,			\
		      double* restrict epsminus,			\
		      double* restrict tempminus,			\
		      double* restrict yeminus,				\
		      double* restrict Bvcxminus,			\
		      double* restrict Bvcyminus,			\
		      double* restrict Bvczminus,			\
		      double* restrict psidcminus,			\
		      double* restrict rhoplus,				\
		      double* restrict velxplus,			\
		      double* restrict velyplus,			\
		      double* restrict velzplus,			\
		      double* restrict epsplus,				\
		      double* restrict tempplus,			\
		      double* restrict yeplus,				\
		      double* restrict Bvcxplus,			\
		      double* restrict Bvcyplus,			\
		      double* restrict Bvczplus,			\
		      double* restrict psidcplus);

//template template<bool do_temp, bool do_ye, bool do_mhd, bool dc_flag, bool do_ppm_detect>


// first stuff without MHD
INSTANTIATE(false,false,false,false,false)
INSTANTIATE(false,false,false,false,true)
// do_temp can only be true of do_ye is also true
INSTANTIATE(true,true,false,false,false)
INSTANTIATE(true,true,false,false,true)
// but do_ye can be true and do_temp can be false
INSTANTIATE(false,true,false,false,false)
INSTANTIATE(false,true,false,false,true)

// with MHD, but dc_flag false
INSTANTIATE(false,false,true,false,false)
INSTANTIATE(false,false,true,false,true)
// do_temp can only be true of do_ye is also true
INSTANTIATE(true,true,true,false,false)
INSTANTIATE(true,true,true,false,true)
// but do_ye can be true and do_temp can be false
INSTANTIATE(false,true,true,false,false)
INSTANTIATE(false,true,true,false,true)

// with MHD, dc_flag true
INSTANTIATE(false,false,true,true,false)
INSTANTIATE(false,false,true,true,true)
// do_temp can only be true of do_ye is also true
INSTANTIATE(true,true,true,true,false)
INSTANTIATE(true,true,true,true,true)
// but do_ye can be true and do_temp can be false
INSTANTIATE(false,true,true,true,false)
INSTANTIATE(false,true,true,true,true)
