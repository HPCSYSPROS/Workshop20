#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"


#include "util_Table.h"
#include <math.h>

#include "ZelmaniQuadWaveExtract.hh"


void ZelmaniQuadWaveLegacy(CCTK_REAL rho,
			  CCTK_REAL sdetg,
			  CCTK_REAL alp,
			  CCTK_REAL betax,
			  CCTK_REAL betay,
			  CCTK_REAL betaz,
			  CCTK_REAL w_lo,
			  CCTK_REAL press,
			  CCTK_REAL eps,
			  CCTK_REAL velx,
			  CCTK_REAL vely,
			  CCTK_REAL velz,
			  CCTK_REAL x,
			  CCTK_REAL y,
			  CCTK_REAL z,
			  CCTK_REAL *dIdotxx,
			  CCTK_REAL *dIdotxy,
			  CCTK_REAL *dIdotxz,
			  CCTK_REAL *dIdotyy,
			  CCTK_REAL *dIdotyz,
			  CCTK_REAL *dIdotzz)
{

  CCTK_REAL rhostar;

  CCTK_REAL localweight = 1.0e0;
  CCTK_REAL mariage_term;

  rhostar = localweight * rho * sdetg * w_lo;

  mariage_term = 2.0e0 / 3.0e0 * 
    ( velx*x + vely*y + velz*z );
	  
  
  *dIdotxx = rhostar * ( 2.0e0 * velx*x - mariage_term );
  
  *dIdotyy = rhostar * ( 2.0e0 * vely*y - mariage_term );
  
  *dIdotzz = rhostar * ( 2.0e0 * velz*z - mariage_term );

  *dIdotxy = rhostar * (velx*y + vely*x);
  
  *dIdotxz = rhostar * (velx*z + velz*x);	
  
  *dIdotyz = rhostar * (vely*z + velz*y);


  return;

}

void ZelmaniQuadWaveShibata(CCTK_REAL rho,
			  CCTK_REAL sdetg,
			  CCTK_REAL alp,
			  CCTK_REAL betax,
			  CCTK_REAL betay,
			  CCTK_REAL betaz,
			  CCTK_REAL w_lo,
			  CCTK_REAL press,
			  CCTK_REAL eps,
			  CCTK_REAL velx,
			  CCTK_REAL vely,
			  CCTK_REAL velz,
			  CCTK_REAL x,
			  CCTK_REAL y,
			  CCTK_REAL z,
			  CCTK_REAL *dIdotxx,
			  CCTK_REAL *dIdotxy,
			  CCTK_REAL *dIdotxz,
			  CCTK_REAL *dIdotyy,
			  CCTK_REAL *dIdotyz,
			  CCTK_REAL *dIdotzz)
{

  CCTK_REAL rhostar;

  CCTK_REAL localweight = 1.0e0;
  CCTK_REAL mariage_term;

  CCTK_REAL svx,svy,svz;

  svx = alp*velx-betax;
  svy = alp*vely-betay;
  svz = alp*velz-betaz;

  rhostar = localweight * rho * sdetg * w_lo;

  mariage_term = 2.0e0 / 3.0e0 * 
    ( svx*x + svy*y + svz*z );
	  
  
  *dIdotxx = rhostar * ( 2.0e0 * svx*x - mariage_term );
  
  *dIdotyy = rhostar * ( 2.0e0 * svy*y - mariage_term );
  
  *dIdotzz = rhostar * ( 2.0e0 * svz*z - mariage_term );

  *dIdotxy = rhostar * (svx*y + svy*x);
  
  *dIdotxz = rhostar * (svx*z + svz*x);	
  
  *dIdotyz = rhostar * (svy*z + svz*y);


  return;

}


void ZelmaniQuadWaveLegacyRadiusCriterion(CCTK_REAL rho,
					 CCTK_REAL sdetg,
					 CCTK_REAL alp,
					 CCTK_REAL betax,
					 CCTK_REAL betay,
					 CCTK_REAL betaz,
					 CCTK_REAL w_lo,
					 CCTK_REAL press,
					 CCTK_REAL eps,
					 CCTK_REAL velx,
					 CCTK_REAL vely,
					 CCTK_REAL velz,
					 CCTK_REAL x,
					 CCTK_REAL y,
					 CCTK_REAL z,
					 CCTK_REAL *dIdotxx,
					 CCTK_REAL *dIdotxy,
					 CCTK_REAL *dIdotxz,
					 CCTK_REAL *dIdotyy,
					 CCTK_REAL *dIdotyz,
					 CCTK_REAL *dIdotzz)
{

  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL rhostar;

  CCTK_REAL localweight = 1.0e0;
  CCTK_REAL mariage_term;
  CCTK_REAL radius;

  radius = sqrt(x*x + y*y + z*z);
  
  if(radius >= outer_radius) {
    localweight = 0.0e0;
  }

  rhostar = localweight * rho * sdetg * w_lo;

  mariage_term = 2.0e0 / 3.0e0 * 
    ( velx*x + vely*y + velz*z );
	  
  
  *dIdotxx = rhostar * ( 2.0e0 * velx*x - mariage_term );
  
  *dIdotyy = rhostar * ( 2.0e0 * vely*y - mariage_term );
  
  *dIdotzz = rhostar * ( 2.0e0 * velz*z - mariage_term );

  *dIdotxy = rhostar * (velx*y + vely*x);
  
  *dIdotxz = rhostar * (velx*z + velz*x);	
  
  *dIdotyz = rhostar * (vely*z + velz*y);


  return;

}

void ZelmaniQuadWaveEdensRadiusCriterion(CCTK_REAL rho,
					 CCTK_REAL sdetg,
					 CCTK_REAL alp,
					 CCTK_REAL betax,
					 CCTK_REAL betay,
					 CCTK_REAL betaz,
					 CCTK_REAL w_lo,
					 CCTK_REAL press,
					 CCTK_REAL eps,
					 CCTK_REAL velx,
					 CCTK_REAL vely,
					 CCTK_REAL velz,
					 CCTK_REAL x,
					 CCTK_REAL y,
					 CCTK_REAL z,
					 CCTK_REAL *dIdotxx,
					 CCTK_REAL *dIdotxy,
					 CCTK_REAL *dIdotxz,
					 CCTK_REAL *dIdotyy,
					 CCTK_REAL *dIdotyz,
					 CCTK_REAL *dIdotzz)
{

  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL rhostar;
  CCTK_REAL h;

  CCTK_REAL localweight = 1.0e0;
  CCTK_REAL mariage_term;
  CCTK_REAL radius;

  radius = sqrt(x*x + y*y + z*z);
  
  if(radius >= outer_radius) {
    localweight = 0.0e0;
  }

  h = 1.0e0 + eps + press/(rho+1.0e-50);

  rhostar = localweight * sdetg * (rho * h * w_lo * w_lo - press);

  mariage_term = 2.0e0 / 3.0e0 * 
    ( velx*x + vely*y + velz*z );
	  
  
  *dIdotxx = rhostar * ( 2.0e0 * velx*x - mariage_term );
  
  *dIdotyy = rhostar * ( 2.0e0 * vely*y - mariage_term );
  
  *dIdotzz = rhostar * ( 2.0e0 * velz*z - mariage_term );

  *dIdotxy = rhostar * (velx*y + vely*x);
  
  *dIdotxz = rhostar * (velx*z + velz*x);	
  
  *dIdotyz = rhostar * (vely*z + velz*y);


  return;

}

void ZelmaniQuadWaveLegacyDensityCriterion(CCTK_REAL rho,
					  CCTK_REAL sdetg,
					  CCTK_REAL alp,
					  CCTK_REAL betax,
					  CCTK_REAL betay,
					  CCTK_REAL betaz,
					  CCTK_REAL w_lo,
					  CCTK_REAL press,
					  CCTK_REAL eps,
					  CCTK_REAL velx,
					  CCTK_REAL vely,
					  CCTK_REAL velz,
					  CCTK_REAL x,
					  CCTK_REAL y,
					  CCTK_REAL z,
					  CCTK_REAL *dIdotxx,
					  CCTK_REAL *dIdotxy,
					  CCTK_REAL *dIdotxz,
					  CCTK_REAL *dIdotyy,
					  CCTK_REAL *dIdotyz,
					  CCTK_REAL *dIdotzz)
{

  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL rhostar;

  CCTK_REAL localweight = 1.0e0;
  CCTK_REAL mariage_term;

  if(rho <= rho_min*CGS_RHO_TO_CACTUS) {
    localweight = 0.0e0;
  }

  rhostar = localweight * rho * sdetg * w_lo;

  mariage_term = 2.0e0 / 3.0e0 * 
    ( velx*x + vely*y + velz*z );
	  
  
  *dIdotxx = rhostar * ( 2.0e0 * velx*x - mariage_term );
  
  *dIdotyy = rhostar * ( 2.0e0 * vely*y - mariage_term );
  
  *dIdotzz = rhostar * ( 2.0e0 * velz*z - mariage_term );

  *dIdotxy = rhostar * (velx*y + vely*x);
  
  *dIdotxz = rhostar * (velx*z + velz*x);	
  
  *dIdotyz = rhostar * (vely*z + velz*y);


  return;

}
