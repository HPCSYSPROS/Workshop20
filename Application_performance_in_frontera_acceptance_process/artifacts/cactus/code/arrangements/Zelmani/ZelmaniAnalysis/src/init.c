#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "Symmetry.h"

#include "assert.h"
#include "math.h"



#include "util_Table.h"



void ZelmaniAnalysis_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int nx = cctkGH->cctk_lsh[0];
  int ny = cctkGH->cctk_lsh[1];
  int nz = cctkGH->cctk_lsh[2];
  int i,j,k,is,ijk,ijks;

  for(k=0;k<nz;k++) 
    for(j=0;j<ny;j++)
      for(i=0;i<nx;i++) {
	ijk = CCTK_GFINDEX3D(cctkGH,i,j,k);

	proper_mass_local[ijk] = 0.0e0;
	gravitational_mass_local[ijk] = 0.0e0;
	gravitational_mass_local_in1e12[ijk] = 0.0e0;
	baryonic_mass_local[ijk] = 0.0e0;
	baryonic_mass_local_in1e12[ijk] = 0.0e0;
	angular_momentum_local[ijk] = 0.0e0;
	kinetic_energy_local[ijk] = 0.0e0;
	total_kinetic_energy_local[ijk] = 0.0e0;

	if(do_ic_analysis) {
	  proper_mass_ic_local[ijk] = 0.0e0;
	  kinetic_energy_ic_local[ijk] = 0.0e0;
	  angular_momentum_ic_local[ijk] = 0.0e0;
	  gravitational_mass_ic_local[ijk] = 0.0e0;
	}
      }

  if(number_of_spheres>0) {
    for(int m=0;m<number_of_spheres;m++) {
      adm_mass_volume[m] = 0.0e0;
      baryonic_mass_interior[m] = 0.0e0;
      if(do_enu_interior) {
	enu_interior[m] = 0.0e0;
      }
      for(k=0;k<nz;k++) 
	for(j=0;j<ny;j++)
	  for(i=0;i<nx;i++) {
	    int i4D = CCTK_VECTGFINDEX3D(cctkGH,i,j,k,m);
	    adm_mass_volume_local[i4D] = 0.0e0;
	    baryonic_mass_interior_local[i4D] = 0.0e0;
	    if(do_enu_interior) {
	      enu_interior_local[i4D] = 0.0e0;
	    }
	  }
    }
  }

#if 0
  for(is=0;is < number_of_spheres;is++) {  
    for(i=0; i< nx; i++)
      for(j=0; j< ny; j++)
	for(k=0; k <nz; k++)
	  {
	    ijks = CCTK_GFINDEX4D(cctkGH,i,j,k,is);
	    //	    adm_mass_volume_local[ijks] = 0.0e0;
	  }
  }
#endif

}

