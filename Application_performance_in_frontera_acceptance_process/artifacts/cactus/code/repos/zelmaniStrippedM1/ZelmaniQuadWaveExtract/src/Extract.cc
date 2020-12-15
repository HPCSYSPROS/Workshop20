#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"


#include "util_Table.h"


extern "C" { void ZelmaniQuadWaveExtract_Extract(CCTK_ARGUMENTS);
}

void ZelmaniQuadWaveExtract_Extract(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  // hey. this is our global mode routine that
  // performs the reductions.

  // here we go: we need some handle on things!

  if(!(*dostuff) && !(cctk_iteration==0)) {
    return;
  }

  if(verbose) {
    CCTK_INFO("Performing reductions");
  }

  int varindex = -1;
  int ierr = 0;

  int cf = 1.0;

  int reduction_handle = CCTK_ReductionHandle("sum");

  varindex = CCTK_VarIndex("ZelmaniQuadWaveExtract::dBaryMass");
  assert(varindex>=0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)BaryMass, 1, varindex);
  assert(!ierr);

  varindex = CCTK_VarIndex("ZelmaniQuadWaveExtract::dIdotxx");
  assert(varindex>=0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)Idotxx, 1, varindex);
  assert(!ierr);


  varindex = CCTK_VarIndex("ZelmaniQuadWaveExtract::dIdotxy");

  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)Idotxy, 1, varindex);
  assert(!ierr);

  varindex = CCTK_VarIndex("ZelmaniQuadWaveExtract::dIdotxz");
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)Idotxz, 1, varindex);
  assert(!ierr);

  varindex = CCTK_VarIndex("ZelmaniQuadWaveExtract::dIdotyy");
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)Idotyy, 1, varindex);
  assert(!ierr);

  varindex = CCTK_VarIndex("ZelmaniQuadWaveExtract::dIdotyz");
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)Idotyz, 1, varindex);
  assert(!ierr);

  varindex = CCTK_VarIndex("ZelmaniQuadWaveExtract::dIdotzz");
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)Idotzz, 1, varindex);
  assert(!ierr);

  //  varindex = CCTK_VarIndex("grid::x");
  //ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
  //		     1, CCTK_VARIABLE_REAL, (void *)xsum, 1, varindex);
  //assert(!ierr);

  // let's count all the cells
  reduction_handle = CCTK_ReductionHandle("count");

  varindex = CCTK_VarIndex("ZelmaniQuadWaveExtract::dIdotxx");
  assert(varindex>=0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)volume, 1, varindex);
  assert(!ierr);

  // let's find the maximum density on the grid
  varindex = CCTK_VarIndex("hydrobase::rho");
  int vartype = CCTK_VarTypeI(varindex);
  reduction_handle = CCTK_ReductionHandle("maximum");
  ierr = CCTK_Reduce(cctkGH,-1,reduction_handle,1,vartype,rhomax,1,varindex);
  assert(!ierr);

  *rhomax = (*rhomax)*6.1755e+17;

  // now multiply by the coordinate volume element

  //CCTK_REAL d3x = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];

  int symfac = 1.0;
  if (CCTK_EQUALS(domain,"bitant")){
    symfac = 2.0e0;
  } else if (CCTK_EQUALS(domain,"octant")){
    symfac = 8.0e0;
  }

  //*volume = (*volume) * d3x;

  *BaryMass = cf * (*BaryMass);
  *Idotxx = cf * symfac * (*Idotxx);
  *Idotxy = cf * symfac * (*Idotxy);
  *Idotxz = cf * symfac * (*Idotxz);
  *Idotyy = cf * symfac * (*Idotyy);
  *Idotyz = cf * symfac * (*Idotyz);
  *Idotzz = cf * symfac * (*Idotzz);


  return;
}

