#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"


#include "util_Table.h"
#include <math.h>

#include "ZelmaniQuadWaveExtract.hh"


#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

extern "C" { void ZelmaniQuadWaveExtract_LocalStuff(CCTK_ARGUMENTS);
}

static
void SpatialDeterminant(CCTK_REAL gxx,
			CCTK_REAL gxy,
			CCTK_REAL gxz,
			CCTK_REAL gyy,
			CCTK_REAL gyz,
			CCTK_REAL gzz,
			CCTK_REAL *detg);

void ZelmaniQuadWaveExtract_LocalStuff(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  // here we do all the stuff that is required
  // locally for the later reduction operation


  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  if (*volume_form_state == 0)
  {
     CCTK_WARN(0, "Activate storage for the volume form in your 'Coordinates' implementation!");
  }


  if(!vphys) {
#pragma omp parallel for 
    for (int k=0;k<nz;k++)
      for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++) {
	  
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double detg = 0.0;
	  

	  
	  SpatialDeterminant(gxx[index],gxy[index],gxz[index],
			     gyy[index],gyz[index],gzz[index],
			     &detg);

	  dBaryMass[index] = w_lorentz[index]*rho[index]*sqrt(detg)*volume_form[index];
	  
	  // call function pointer that returns integrands
	  ZelmaniQuadWaveIntegrandPointer(rho[index],
					  sqrt(detg),
					  alp[index],
					  betax[index],
					  betay[index],
					  betaz[index],
					  w_lorentz[index],
					  press[index],
					  eps[index],
					  velx[index],
					  vely[index],
					  velz[index],
					  x[index],
					  y[index],
					  z[index],
					  &dIdotxx[index],
					  &dIdotxy[index],
					  &dIdotxz[index],
					  &dIdotyy[index],
					  &dIdotyz[index],
					  &dIdotzz[index]);	    

          dIdotxx[index] *= volume_form[index];
          dIdotxy[index] *= volume_form[index];
          dIdotxz[index] *= volume_form[index];
          dIdotyy[index] *= volume_form[index];
          dIdotyz[index] *= volume_form[index];
          dIdotzz[index] *= volume_form[index];

	}
  } else
  {
#pragma omp parallel for 
    for (int k=0;k<nz;k++)
      for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++) {
	  
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double detg = 0.0;

	  SpatialDeterminant(gxx[index],gxy[index],gxz[index],
			     gyy[index],gyz[index],gzz[index],
			     &detg);
	  
	  dBaryMass[index] = w_lorentz[index]*rho[index]*sqrt(detg)*volume_form[index];	  

	  double lvx = sqrt(gxx[index])*velx[index];
	  
	  double lvy = sqrt(gyy[index])*vely[index];
	  
	  double lvz = sqrt(gzz[index])*velz[index];

	  // call function pointer that returns integrands
	  ZelmaniQuadWaveIntegrandPointer(rho[index],
					  sqrt(detg),
					  alp[index],
					  betax[index],
					  betay[index],
					  betaz[index],
					  w_lorentz[index],
					  press[index],
					  eps[index],
					  lvx,
					  lvy,
					  lvz,
					  x[index],
					  y[index],
					  z[index],
					  &dIdotxx[index],
					  &dIdotxy[index],
					  &dIdotxz[index],
					  &dIdotyy[index],
					  &dIdotyz[index],
					  &dIdotzz[index]);	    
		
	   dIdotxx[index] *= volume_form[index];
           dIdotxy[index] *= volume_form[index];
           dIdotxz[index] *= volume_form[index];
           dIdotyy[index] *= volume_form[index];
           dIdotyz[index] *= volume_form[index];
           dIdotzz[index] *= volume_form[index];
	}
    
  }


  return;



}


static
void SpatialDeterminant(CCTK_REAL gxx,
			CCTK_REAL gxy,
			CCTK_REAL gxz,
			CCTK_REAL gyy,
			CCTK_REAL gyz,
			CCTK_REAL gzz,
			CCTK_REAL *detg)
{

  *detg = -gxz*gxz*gyy + 2.0*gxy*gxz*gyz 
    - gxx*gyz*gyz - gxy*gxy*gzz + gxx*gyy*gzz;

  return;
}
