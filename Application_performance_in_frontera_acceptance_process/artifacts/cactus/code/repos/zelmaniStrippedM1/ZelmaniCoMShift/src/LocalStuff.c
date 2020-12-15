#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include <math.h>
#include <stdio.h>

#define RHOCGSTOCACTUS (1.0e0/6.1755e17)

void ZelmaniCoMShift_Init(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];
  int i,j,k,index;

  *Mass = 0.0;
  *Mr = 0.0;
  *Mx = 0.0;
  *My = 0.0;
  *Mz = 0.0;
  *have_good_data = 0;

#pragma omp parallel for private(index,i,j,k)
  for (k=0;k<nz;k++)
    for (j=0;j<ny;j++)
      for (i=0;i<nx;i++) {
	
	index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	dMass[index] = 0.0e0;
	dMx[index] = 0.0e0;
	dMy[index] = 0.0e0;
	dMz[index] = 0.0e0;

	dMass_p[index] = 0.0e0;
	dMx_p[index] = 0.0e0;
	dMy_p[index] = 0.0e0;
	dMz_p[index] = 0.0e0;

	dMass_p_p[index] = 0.0e0;
	dMx_p_p[index] = 0.0e0;
	dMy_p_p[index] = 0.0e0;
	dMz_p_p[index] = 0.0e0;
	
      }


}

void ZelmaniCoMShift_LocalStuff(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  int i,j,k;
  int index;

  CCTK_REAL localweight = 1.0e0;

  CCTK_REAL tiny,rt,rt2,rt4,x2,y2;

  CCTK_REAL cutrho,cutradius;
  CCTK_REAL cutrhoCoM,cutradiusCoM;

  cutradius = 1.0e20;
  cutrhoCoM = 0.0e0;

  if(use_CoM_radius) {
    cutrhoCoM = 1.0e0*RHOCGSTOCACTUS;
    cutradiusCoM = CoM_radius;
  } else {
    cutradiusCoM = cutradius;
  }

  tiny = 1.e-50;

  // We don`t want to do anything until shortly before our start time...

  if(verbose_level > 1) {
    CCTK_INFO("Calculating local expressions");
  }

#pragma omp parallel for private(index,i,j,k,localweight)
  for (k=0;k<nz;k++)
    for (j=0;j<ny;j++)
      for (i=0;i<nx;i++) {
	
	index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	localweight = 1.0e0;

	if( (rho[index] <= cutrhoCoM) || (r[index] > cutradiusCoM) ) localweight = 0.0e0;

	dMass[index] = localweight*dens[index];
	dMx[index] = localweight*dens[index]*x[index];
	dMy[index] = localweight*dens[index]*y[index];
	dMz[index] = localweight*dens[index]*z[index];
	
      }




  return;
}

void ZelmaniCoMShift_LocalStuffRecover(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;


  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  int i,j,k;
  int index;

  CCTK_REAL localweight = 1.0e0;

  CCTK_REAL tiny,rt,rt2,rt4,x2,y2;

  CCTK_REAL cutrho,cutradius;
  CCTK_REAL cutrhoCoM,cutradiusCoM;

  cutradius = 1.0e20;
  cutrhoCoM = 0.0e0;

  if(use_CoM_radius) {
    cutrhoCoM = 1.0e0*RHOCGSTOCACTUS;
    cutradiusCoM = CoM_radius;
  } else {
    cutradiusCoM = cutradius;
  }

  tiny = 1.e-50;

  // We don`t want to do anything until shortly before our start time...

  if(verbose_level > 1) {
    CCTK_INFO("Calculating local expressions at recovery");
  }

  if(setup_all_timelevels_after_recovery) {

    *Mass = 0.0;
    *Mr = 0.0;
    *Mx = 0.0;
    *My = 0.0;
    *Mz = 0.0;
    *have_good_data = 0;

#pragma omp parallel for private(index,i,j,k,localweight)
    for (k=0;k<nz;k++)
      for (j=0;j<ny;j++)
	for (i=0;i<nx;i++) {
	  
	  index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  
	  localweight = 1.0e0;
	  
	  if( (rho[index] <= cutrhoCoM) || (r[index] > cutradiusCoM) ) localweight = 0.0e0;
	  
	  dMass[index] = localweight*dens[index];
	  dMx[index] = localweight*dens[index]*x[index];
	  dMy[index] = localweight*dens[index]*y[index];
	  dMz[index] = localweight*dens[index]*z[index];

	  dMass_p[index] = localweight*dens_p[index];
	  dMx_p[index] = localweight*dens_p[index]*x[index];
	  dMy_p[index] = localweight*dens_p[index]*y[index];
	  dMz_p[index] = localweight*dens_p[index]*z[index];

	  dMass_p_p[index] = localweight*dens_p_p[index];
	  dMx_p_p[index] = localweight*dens_p_p[index]*x[index];
	  dMy_p_p[index] = localweight*dens_p_p[index]*y[index];
	  dMz_p_p[index] = localweight*dens_p_p[index]*z[index];
	  
	}
  }
  else
    {
#pragma omp parallel for private(index,i,j,k,localweight)
    for (k=0;k<nz;k++)
      for (j=0;j<ny;j++)
	for (i=0;i<nx;i++) {
	  
	  index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  
	  localweight = 1.0e0;
	  
	  if( (rho[index] <= cutrhoCoM) || (r[index] > cutradiusCoM) ) localweight = 0.0e0;
	  
	  dMass[index] = localweight*dens[index];
	  dMx[index] = localweight*dens[index]*x[index];
	  dMy[index] = localweight*dens[index]*y[index];
	  dMz[index] = localweight*dens[index]*z[index];
      }
    }


  return;
}
