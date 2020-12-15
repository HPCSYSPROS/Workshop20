#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include <math.h>
#include <stdio.h>

#define RHOCGSTOCACTUS (1.0e0/6.1755e17)

void ZelmaniCoMShift_ChangeShift(CCTK_ARGUMENTS)
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

  // We don`t want to do anything until shortly before our start time...


  for (k=0;k<nz;k++)
    for (j=0;j<ny;j++)
      for (i=0;i<nx;i++) {
	
	index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	localweight = 1.0e0;

	betax[index] = betax[index] - *Mx/CoM_shift_factor;
	betay[index] = betay[index] - *My/CoM_shift_factor;
	betaz[index] = betaz[index] - *Mz/CoM_shift_factor;

      }




  return;
}
