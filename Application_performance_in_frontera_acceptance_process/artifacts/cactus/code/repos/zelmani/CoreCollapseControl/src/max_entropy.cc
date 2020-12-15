#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include <math.h>
#include <stdio.h>

#include "util_Table.h"
#include <assert.h>


void CoreCollapseControl_LocalSMax(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int nx = cctk_lsh[0];
  const int ny = cctk_lsh[1];
  const int nz = cctk_lsh[2];

  *local_entropy_max = 0.;;

  if(!(*in_prebounce) && !(*force_check)) return;

#pragma omp parallel
  {
    CCTK_REAL entropy_max = 0.;
    for (int k=0;k<nz;k++)
      for (int j=0;j<ny;j++)
        for (int i=0;i<nx;i++) {

          int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

          if(r[index] <= 100.0e0) {
            if(entropy[index] > entropy_max) {
              entropy_max = entropy[index];
            }
          }
        }
#pragma omp critical
    if(entropy_max > *local_entropy_max) {
      *local_entropy_max = entropy_max;
    }
  }

  return;
}

void CoreCollapseControl_GlobalSMax(CCTK_ARGUMENTS) 
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if( ((((cctk_iteration) % rho_max_every) != 0) && !(*force_check) )|| (*bounce)) 
    return;

  if(!(*in_prebounce) && !(*force_check)) return;

// lets find the maximum entropy on the grid
  const int varindex = CCTK_VarIndex("CoreCollapseControl::local_entropy_max");
  assert(varindex >= 0);
  const int reduction_handle = CCTK_ReductionHandle("maximum");
  assert(reduction_handle >= 0);
  const int ierr = CCTK_Reduce(cctkGH,-1,reduction_handle,1,CCTK_VARIABLE_REAL,global_entropy_max,1,varindex);
  assert(!ierr);

  return;
}
