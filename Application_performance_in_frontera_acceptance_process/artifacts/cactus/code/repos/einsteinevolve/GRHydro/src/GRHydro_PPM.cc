#include <cmath>
#include <algorithm>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

using namespace std;

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


static inline void steep(const double * const restrict x, 
			 const double * const restrict dx, 
			 double * const restrict dmx, 
			 const int i) {
  if ( (x[i+1] - x[i]) * (x[i]-x[i-1]) > 0.0 ) {
    dmx[i] = copysign(1.0,dx[i]) * MIN(fabs(dx[i]),
					MIN(2.0*fabs(x[i]-x[i-1]),
					     2.0*fabs(x[i+1]-x[i])));
  } else {
    dmx[i] = 0.0;
  }
}


static inline void monotonize(double* const restrict xminus,
			      const double* const restrict x,
			      double* const restrict xplus,
			      const int i) {

  if (  !(xplus[i]==x[i] && x[i]==xminus[i]) 
	&& ( (xplus[i]-x[i])*(x[i]-xminus[i]) <= 0.0 ) ) 
    {
      xminus[i] = x[i];
      xplus[i] = x[i];
    }  else if( 6.0 * (xplus[i]-xminus[i]) * 
		(x[i]-0.5*(xplus[i]+xminus[i])) >
		(xplus[i]-xminus[i])*(xplus[i]-xminus[i]) )
    {
      xminus[i] = 3.0*x[i]-2.0*xplus[i]; 
    } else if( 6.0 * (xplus[i]-xminus[i]) * 
	       (x[i]-0.5*(xplus[i]+xminus[i])) <
	       -(xplus[i]-xminus[i])*(xplus[i]-xminus[i]) ) 
    {
      xplus[i] = 3.0*x[i]-2.0*xminus[i]; 
    }
  
  return;
}



template<bool do_temp, bool do_ye, bool do_mhd, 
	 bool dc_flag, bool do_ppm_detect>
void GRHydro_ppm1d_cxx(const int nx, 
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

  double drho[nx], dpress[nx], deps[nx], d2rho[nx]; 
  double dmrho[nx], dmeps[nx];
  double dmtemp[nx], dmye[nx];
  double dtemp[nx], dye[nx];
  double dvelx[nx], dvely[nx], dvelz[nx];
  double dmvelx[nx], dmvely[nx], dmvelz[nx];
  double dBvcx[nx], dBvcy[nx], dBvcz[nx], dpsidc[nx];
  double dmBvcx[nx], dmBvcy[nx], dmBvcz[nx], dmpsidc[nx];
  double tilde_flatten[nx];
  
  const double onesixth = 1.0/6.0;

  //  Average slopes delta_m(a). See (1.7) of Colella and Woodward, p.178
  //  This is the expression for an even grid.
  for(int i=1; i < nx-1; ++i) {
    drho[i]   = 0.5 * (rho[i+1]-rho[i-1]);
    d2rho[i]  = rho[i+1] - 2.0 * rho[i] + rho[i-1];
    dpress[i] = (press[i+1]-press[i-1]);
    deps[i]   = 0.5 * (eps[i+1]-eps[i-1]);
    dvelx[i]  = 0.5 * (velx[i+1]-velx[i-1]);
    dvely[i]  = 0.5 * (vely[i+1]-vely[i-1]);
    dvelz[i]  = 0.5 * (velz[i+1]-velz[i-1]);
    if(do_ye) {
      if(do_temp) {
	dtemp[i]  = 0.5 * (temp[i+1]-temp[i-1]);
      }
      dye[i]    = 0.5 * (ye[i+1]-ye[i-1]);
    }
    if(do_mhd) {
      dBvcx[i]  = 0.5 * (Bvcx[i+1]-Bvcx[i-1]);
      dBvcy[i]  = 0.5 * (Bvcy[i+1]-Bvcy[i-1]);
      dBvcz[i]  = 0.5 * (Bvcz[i+1]-Bvcz[i-1]);
      if(dc_flag) {
	dpsidc[i] = 0.5 * (psidc[i+1]-psidc[i-1]);
      }
    }
  }

  //  Steepened slope. See (1.8) of Colella and Woodward, p.178
  for(int i=1; i<nx-1; ++i) {
    steep(rho,drho,dmrho,i);
    steep(eps,deps,dmeps,i);
    steep(velx,dvelx,dmvelx,i);
    steep(vely,dvely,dmvely,i);
    steep(velz,dvelz,dmvelz,i);
    if(do_ye) {
      if(do_temp) {
	steep(temp,dtemp,dmtemp,i);
      }
      steep(ye,dye,dmye,i);
    }
    if(do_mhd) {
      steep(Bvcx,dBvcx,dmBvcx,i);
      steep(Bvcy,dBvcy,dmBvcy,i);
      steep(Bvcz,dBvcz,dmBvcz,i);
      if(dc_flag) {
	steep(psidc,dpsidc,dmpsidc,i);
      }
    }
  }

  // Initial boundary states. See (1.9) of Colella and Woodward, p.178
  for(int i=1; i<nx-2; ++i) {
    rhoplus[i] = 0.5 * (rho[i] + rho[i+1]) + 
      (dmrho[i] - dmrho[i+1]) * onesixth;
    epsplus[i] = 0.5 * (eps[i] + eps[i+1]) + 
      (dmeps[i] - dmeps[i+1]) * onesixth;
    velxplus[i] = 0.5 * (velx[i] + velx[i+1]) + 
      (dmvelx[i] - dmvelx[i+1]) * onesixth;
    velyplus[i] = 0.5 * (vely[i] + vely[i+1]) + 
      (dmvely[i] - dmvely[i+1]) * onesixth;
    velzplus[i] = 0.5 * (velz[i] + velz[i+1]) + 
      (dmvelz[i] - dmvelz[i+1]) * onesixth;
    if(do_ye) {
      if(do_temp) {
	tempplus[i] = 0.5 * (temp[i] + temp[i+1]) + 
	  (dmtemp[i] - dmtemp[i+1]) * onesixth;
      }
      yeplus[i] = 0.5 * (ye[i] + ye[i+1]) + 
	(dmye[i] - dmye[i+1]) * onesixth;
    }
    if(do_mhd) {
      Bvcxplus[i] = 0.5 * (Bvcx[i] + Bvcx[i+1]) + 
	(dmBvcx[i] - dmBvcx[i+1]) * onesixth;
      Bvcyplus[i] = 0.5 * (Bvcy[i] + Bvcy[i+1]) + 
	(dmBvcy[i] - dmBvcy[i+1]) * onesixth;
      Bvczplus[i] = 0.5 * (Bvcz[i] + Bvcz[i+1]) + 
	(dmBvcz[i] - dmBvcz[i+1]) * onesixth;
      if(dc_flag) {
	psidcplus[i] = 0.5 * (psidc[i] + psidc[i+1]) + 
	  (dmpsidc[i] - dmpsidc[i+1]) * onesixth;
      }
    }
  }

  // fill minus states
  for(int i=1; i<nx-2; ++i) {
    rhominus[i+1] = rhoplus[i];
    epsminus[i+1] = epsplus[i];
    velxminus[i+1] = velxplus[i];
    velyminus[i+1] = velyplus[i];
    velzminus[i+1] = velzplus[i];
    if(do_ye) {
      if(do_temp) {
	tempminus[i+1] = tempplus[i];
      }
      yeminus[i+1] = yeplus[i];
    }
    if(do_mhd) {
      Bvcxminus[i+1] = Bvcxplus[i];
      Bvcyminus[i+1] = Bvcyplus[i];
      Bvczminus[i+1] = Bvczplus[i];
      if(dc_flag) {
	psidcminus[i+1] = psidcplus[i];
      }
    }
  }

  /*
    Discontinuity steepening. See (1.14-17) of C&W.
    This is the detect routine which mat be 
    activated with the ppm_detect parameter
    Note that this part really also depends on the grid being even. 
    Note also that we don''t have access to the gas constant gamma.
    So this is just dropped from eq. (3.2) of C&W.
    We can get around this by just rescaling the constant k0 (ppm_k0 here).
  */

  if(do_ppm_detect) {
    for(int i=2; i<nx-2; ++i) {
      double etatilde = 0.0;
      if ( d2rho[i+1]*d2rho[i-1] < 0.0 
	   && ( fabs(rho[i+1]-rho[i-1]) - ppm_epsilon_shock 
		* MIN(fabs(rho[i+1]), 
		      fabs(rho[i-1])) > 0.0) )
	{
	  etatilde = (rho[i-2] - rho[i+2] + 4.0 * drho[i]) / (drho[i] * 12.0);
	}
      double eta = MAX(0.0, MIN(1.0, ppm_eta1 * (etatilde - ppm_eta2)));
      if (ppm_k0 * fabs(drho[i]) * MIN(press[i-1],press[i+1]) 
	  < fabs(dpress[i]) * MIN(rho[i-1], rho[i+1])) 
	{
	  eta = 0.0;
	}
      rhominus[i] = rhominus[i] * (1.0 - eta) + 
	(rho[i-1] + 0.5 * dmrho[i-1]) * eta;
      rhoplus[i] = rhoplus[i] * (1.0 - eta) + 
	(rho[i+1] - 0.5 * dmrho[i+1]) * eta;
    }
  } 
  
  // flattening
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

  for(int i=GRHydro_stencil-1; i<nx-GRHydro_stencil+1; ++i) {
    monotonize(rhominus,rho,rhoplus,i);
    monotonize(epsminus,eps,epsplus,i);
    monotonize(velxminus,velx,velxplus,i);
    monotonize(velyminus,vely,velyplus,i);
    monotonize(velzminus,velz,velzplus,i);
    if(do_ye) {
      if(do_temp) {
	monotonize(tempminus,temp,tempplus,i);
      }
      monotonize(yeminus,ye,yeplus,i);
    }
    if(do_mhd) {
      monotonize(Bvcxminus,Bvcx,Bvcxplus,i);
      monotonize(Bvcyminus,Bvcy,Bvcyplus,i);
      monotonize(Bvczminus,Bvcz,Bvczplus,i);
      if(dc_flag) {
	monotonize(psidcminus,psidc,psidcplus,i);
      }
    }
  }

  return;
}

//template template<bool do_temp, bool do_ye, bool do_mhd, bool dc_flag, bool do_ppm_detect>
#define INSTANTIATE(do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect)			\
template void GRHydro_ppm1d_cxx<do_temp,do_ye,do_mhd,dc_flag,do_ppm_detect>(const int nx, \
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
