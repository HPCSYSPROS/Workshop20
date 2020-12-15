 /*@@
   @file      adm_mass_volume.c
   @date      
   @author    Frank Loeffler and Luca Baiotti and Christian D. Ott
   @desc 
          Computes the ADM mass as a volume integral.
	  Stolen from AEIDevelopment/ADMMass
   @enddesc 
 @@*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>


#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117068

void ZelmaniAnalysis_ADMMass_Volume(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT i,j,k, ijk, ghost, ti, tj, tk, tl;
    CCTK_REAL detg, idetg;
    CCTK_REAL u[3][3], dg[3][3][3];
    CCTK_REAL radius;

    int is,ijks;

#include "CactusEinstein/ADMMacros/src/macro/UPPERMET_declare.h"

    /* grid-function strides */
    const CCTK_INT di = 1;
    const CCTK_INT dj = cctk_lsh[0];
    const CCTK_INT dk = cctk_lsh[0]*cctk_lsh[1];
    
    /* deonminators for derivatives */
    const CCTK_REAL OneOverTwoDX = 1.0 / (2.0 * CCTK_DELTA_SPACE(0));
    const CCTK_REAL OneOverTwoDY = 1.0 / (2.0 * CCTK_DELTA_SPACE(1));
    const CCTK_REAL OneOverTwoDZ = 1.0 / (2.0 * CCTK_DELTA_SPACE(2));
    
    CCTK_REAL* ADMMass_VolumeMass_pot_x;
    CCTK_REAL* ADMMass_VolumeMass_pot_y;
    CCTK_REAL* ADMMass_VolumeMass_pot_z;
    
    ADMMass_VolumeMass_pot_x = (CCTK_REAL*) 
      malloc(sizeof(CCTK_REAL)*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]);
    ADMMass_VolumeMass_pot_y = (CCTK_REAL*) 
      malloc(sizeof(CCTK_REAL)*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]);
    ADMMass_VolumeMass_pot_z = (CCTK_REAL*) 
      malloc(sizeof(CCTK_REAL)*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]);
							 

    const void *ghost_ptr = CCTK_ParameterGet("ghost_size","Carpet",NULL);
    assert( ghost_ptr != NULL );
    ghost = *(const int *)ghost_ptr;

    for(i=1; i<cctk_lsh[0]-1; i++)
      for(j=1; j<cctk_lsh[1]-1; j++)
	for(k=1; k<cctk_lsh[2]-1; k++)   
	  {
	    ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
	    ADMMass_VolumeMass_pot_x[ijk] = 0.0;
	    ADMMass_VolumeMass_pot_y[ijk] = 0.0;
	    ADMMass_VolumeMass_pot_z[ijk] = 0.0;
	      
#include "CactusEinstein/ADMMacros/src/macro/UPPERMET_guts.h"
	    u[0][0] = UPPERMET_UXX;
	    u[0][1] = UPPERMET_UXY;
	    u[0][2] = UPPERMET_UXZ;
	    u[1][0] = UPPERMET_UXY;
	    u[1][1] = UPPERMET_UYY;
	    u[1][2] = UPPERMET_UYZ;
	    u[2][0] = UPPERMET_UXZ;
	    u[2][1] = UPPERMET_UYZ;
	    u[2][2] = UPPERMET_UZZ;
	    
	    dg[0][0][0] = ( gxx[di+ijk] - gxx[-di+ijk] ) * OneOverTwoDX;
	    dg[0][0][1] = ( gxx[dj+ijk] - gxx[-dj+ijk] ) * OneOverTwoDY;
	    dg[0][0][2] = ( gxx[dk+ijk] - gxx[-dk+ijk] ) * OneOverTwoDZ;
	    
	    dg[0][1][0] = ( gxy[di+ijk] - gxy[-di+ijk] ) * OneOverTwoDX;
	    dg[0][1][1] = ( gxy[dj+ijk] - gxy[-dj+ijk] ) * OneOverTwoDY;
	    dg[0][1][2] = ( gxy[dk+ijk] - gxy[-dk+ijk] ) * OneOverTwoDZ;
	      
	    dg[0][2][0] = ( gxz[di+ijk] - gxz[-di+ijk] ) * OneOverTwoDX;
	    dg[0][2][1] = ( gxz[dj+ijk] - gxz[-dj+ijk] ) * OneOverTwoDY;
	    dg[0][2][2] = ( gxz[dk+ijk] - gxz[-dk+ijk] ) * OneOverTwoDZ;
	      
	    dg[1][0][0] = dg[0][1][0];
	    dg[1][0][1] = dg[0][1][1];
	    dg[1][0][2] = dg[0][1][2];
	      
	    dg[1][1][0] = ( gyy[di+ijk] - gyy[-di+ijk] ) * OneOverTwoDX;
	    dg[1][1][1] = ( gyy[dj+ijk] - gyy[-dj+ijk] ) * OneOverTwoDY;
	    dg[1][1][2] = ( gyy[dk+ijk] - gyy[-dk+ijk] ) * OneOverTwoDZ;
	    
	    dg[1][2][0] = ( gyz[di+ijk] - gyz[-di+ijk] ) * OneOverTwoDX;
	    dg[1][2][1] = ( gyz[dj+ijk] - gyz[-dj+ijk] ) * OneOverTwoDY;
	    dg[1][2][2] = ( gyz[dk+ijk] - gyz[-dk+ijk] ) * OneOverTwoDZ;

	    dg[2][0][0] = dg[0][2][0];
	    dg[2][0][1] = dg[0][2][1];
	    dg[2][0][2] = dg[0][2][2];
	      
	    dg[2][1][0] = dg[1][2][0];
	    dg[2][1][1] = dg[1][2][1];
	    dg[2][1][2] = dg[1][2][2];
	    
	    dg[2][2][0] = ( gzz[di+ijk] - gzz[-di+ijk] ) * OneOverTwoDX;
	    dg[2][2][1] = ( gzz[dj+ijk] - gzz[-dj+ijk] ) * OneOverTwoDY;
	    dg[2][2][2] = ( gzz[dk+ijk] - gzz[-dk+ijk] ) * OneOverTwoDZ;
	    
	    for (ti = 0; ti < 3; ti++)
	      for (tj = 0; tj < 3; tj++)
		for (tk = 0; tk < 3; tk++)
		  {
		    ADMMass_VolumeMass_pot_x[ijk] +=
		      u[ti][tj] * u[tk][0] *
		      ( dg[ti][tk][tj] - dg[ti][tj][tk] );
		    ADMMass_VolumeMass_pot_y[ijk] +=
		      u[ti][tj] * u[tk][1] *
		      ( dg[ti][tk][tj] - dg[ti][tj][tk] );
		    ADMMass_VolumeMass_pot_z[ijk] +=
		      u[ti][tj] * u[tk][2] *
		      ( dg[ti][tk][tj] - dg[ti][tj][tk] );
		  }
	    ADMMass_VolumeMass_pot_x[ijk] *= alp[ijk] * sqrt(DETG_DETG);
	    ADMMass_VolumeMass_pot_y[ijk] *= alp[ijk] * sqrt(DETG_DETG);
	    ADMMass_VolumeMass_pot_z[ijk] *= alp[ijk] * sqrt(DETG_DETG);
	  }

      /* Do not compute in ghost zones */
    for(is=0;is < number_of_spheres;is++) {
      for(i=ghost; i<cctk_lsh[0]-ghost; i++)
	for(j=ghost; j<cctk_lsh[1]-ghost; j++)
	  for(k=ghost; k<cctk_lsh[2]-ghost; k++)
	    {
	      ijk = CCTK_GFINDEX3D(cctkGH, i, j, k);
	      ijks = CCTK_GFINDEX4D(cctkGH,i,j,k,is);
	      adm_mass_volume_local[ijks] = 0.0e0;
	      if(r[ijk] <= radii[is]) {
		adm_mass_volume_local[ijks] =
		  ((ADMMass_VolumeMass_pot_x[CCTK_GFINDEX3D(cctkGH,i+1,j,k)]-
		    ADMMass_VolumeMass_pot_x[CCTK_GFINDEX3D(cctkGH,i-1,j,k)])*OneOverTwoDX+
		   (ADMMass_VolumeMass_pot_y[CCTK_GFINDEX3D(cctkGH,i,j+1,k)]-
		    ADMMass_VolumeMass_pot_y[CCTK_GFINDEX3D(cctkGH,i,j-1,k)])*OneOverTwoDY+
		   (ADMMass_VolumeMass_pot_z[CCTK_GFINDEX3D(cctkGH,i,j,k+1)]-
		    ADMMass_VolumeMass_pot_z[CCTK_GFINDEX3D(cctkGH,i,j,k-1)])*OneOverTwoDZ);
	      }
	    }
    }
    free(ADMMass_VolumeMass_pot_x);
    free(ADMMass_VolumeMass_pot_y);
    free(ADMMass_VolumeMass_pot_z);
      
#include "CactusEinstein/ADMMacros/src/macro/UPPERMET_undefine.h"
}      

void ZelmaniAnalysis_ADMMass_Volume_Global(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS

    CCTK_INT reduction_handle;
    CCTK_REAL radius, grid_spacing_product;
    int is;
    char varname[100];



    reduction_handle = CCTK_ReductionHandle("sum");
    if (reduction_handle < 0)
        CCTK_WARN(0, "Unable to get reduction handle.");

    for(is=0;is<number_of_spheres;is++) {

      sprintf(varname,"%s[%d]","ZelmaniAnalysis::adm_mass_volume_local",is);
    
      if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
		      CCTK_VARIABLE_REAL,
		      (void*)&(adm_mass_volume[is]), 1,
		      CCTK_VarIndex(varname)))
        CCTK_WARN(0, "Error while reducing ZelmaniAnalysis::adm_mass_volume_local");
      
      grid_spacing_product = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];
      
      adm_mass_volume[is] *=
	grid_spacing_product / (16.0*PI);
      
    }

}
