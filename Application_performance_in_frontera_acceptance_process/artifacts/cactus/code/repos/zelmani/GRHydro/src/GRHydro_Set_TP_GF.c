 /*@@
   @file      Set_Trivial_Riemann_Problem_Grid_Function
   @date      Thu May 08
   @author    Frank Loeffler
   @desc 
   This routine sets the grid function for the trivial rieman problem
   bits. This is only done for debugging purposes

   @enddesc 
 @@*/
   
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "SpaceMask.h"

void Set_Trivial_Riemann_Problem_Grid_Function(CCTK_ARGUMENTS)
{
  
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int nx, ny, nz, i, j, k, point;
  int type, state;

  nx = cctk_lsh[0];
  ny = cctk_lsh[1];
  nz = cctk_lsh[2];

  if (flux_direction[0]==1)
  {
    type =SpaceMask_GetTypeBits("Hydro_RiemannProblemX");
    state=SpaceMask_GetStateBits("Hydro_RiemannProblemX", "trivial");
    for (k=0 ; k<nz ; k++)
      for (j=0 ; j<ny ; j++)
        for (i=0 ; i<nx ; i++)
        {
          point = CCTK_GFINDEX3D(cctkGH,i,j,k);
          if (SpaceMask_CheckStateBits(space_mask, point, type, state))
          {
            GRHydro_trivial_rp_gf_x[point]=1;
          }
          else
          {
            GRHydro_trivial_rp_gf_x[point]=0;
          }
        }
  }
  if (flux_direction[0]==2)
  {
    type =SpaceMask_GetTypeBits("Hydro_RiemannProblemY");
    state=SpaceMask_GetStateBits("Hydro_RiemannProblemY", "trivial");
    for (k=0 ; k<nz ; k++)
      for (j=0 ; j<ny ; j++)
        for (i=0 ; i<nx ; i++)
        {
          point = CCTK_GFINDEX3D(cctkGH,i,j,k);
          if (SpaceMask_CheckStateBits(space_mask, point, type, state))
          {
            GRHydro_trivial_rp_gf_y[point]=1;
          }
          else
          {
            GRHydro_trivial_rp_gf_y[point]=0;
          }
        }
  }
  if (flux_direction[0]==3)
  {
    type =SpaceMask_GetTypeBits("Hydro_RiemannProblemZ");
    state=SpaceMask_GetStateBits("Hydro_RiemannProblemZ", "trivial");
    for (k=0 ; k<nz ; k++)
      for (j=0 ; j<ny ; j++)
        for (i=0 ; i<nx ; i++)
        {
          point = CCTK_GFINDEX3D(cctkGH,i,j,k);
          if (SpaceMask_CheckStateBits(space_mask, point, type, state))
          {
            GRHydro_trivial_rp_gf_z[point]=1;
          }
          else
          {
            GRHydro_trivial_rp_gf_z[point]=0;
          }
        }
  }
  return;
}
