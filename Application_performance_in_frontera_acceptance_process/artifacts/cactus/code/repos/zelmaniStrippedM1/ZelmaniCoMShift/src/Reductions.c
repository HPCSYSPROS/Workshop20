#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void ZelmaniCoMShift_Reductions(CCTK_ARGUMENTS)
{


  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int varindex = -1;
  int ierr = 0;
  int cf = 1.0;
  int reduction_handle;
  int vartype;

  CCTK_REAL d3x,tiny;

  CCTK_REAL sym_factor1,sym_factor2,sym_factor3;

  tiny = 1.0e-20;

  if (!(cctk_time >= (start_time))){
    return;
  }

  if ( (cctk_iteration % do_every) != 0  && *have_good_data ) return;

  CCTK_REAL total_mass = 0.0e0;
  sym_factor1 = 1.0e0;
  sym_factor2 = 1.0e0;
  sym_factor3 = 1.0e0;


  d3x = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];

  if (CCTK_EQUALS(domain,"bitant")){
    sym_factor1 = 2.0e0;
    sym_factor2 = 2.0e0;
    sym_factor3 = 0.0e0;
  } else if (CCTK_EQUALS(domain,"octant")){
    sym_factor1 = 8.0e0;
    sym_factor2 = 0.0e0;
    sym_factor3 = 0.0e0;
  }

  reduction_handle = CCTK_ReductionHandle("sum");


  varindex = CCTK_VarIndex("ZelmaniCoMShift::dMass");
  assert(varindex>=0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)&total_mass, 1, varindex);
  assert(!ierr);
  
  total_mass = sym_factor1 * d3x * (total_mass);
  
  if(isnan(total_mass)) {
    CCTK_WARN(0,"Total Mass is NAN");
  }
  
  *Mass = total_mass;


  varindex = CCTK_VarIndex("ZelmaniCoMShift::dMx");
  assert(varindex>=0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)Mx, 1, varindex);
  assert(!ierr);
  
  varindex = CCTK_VarIndex("ZelmaniCoMShift::dMy");
  assert(varindex>=0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)My, 1, varindex);
  assert(!ierr);
  
  varindex = CCTK_VarIndex("ZelmaniCoMShift::dMz");
  assert(varindex>=0);
  ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		     1, CCTK_VARIABLE_REAL, (void *)Mz, 1, varindex);
  assert(!ierr);
  
  *Mx = sym_factor2*d3x*(*Mx)/(*Mass);
  *My = sym_factor2*d3x*(*My)/(*Mass);
  *Mz = sym_factor3*d3x*(*Mz)/(*Mass);
  
  *Mr = sqrt( (*Mx)*(*Mx)+(*My)*(*My)+(*Mz)*(*Mz) );
  
  if(verbose_level > 0) {
    CCTK_VInfo(CCTK_THORNSTRING,"CoM of r<=%15.6E at x=%15.6E y=%15.6E z=%15.6E r=%15.6E",
	       CoM_radius,*Mx,*My,*Mz,*Mr);
  }

  *have_good_data = 1;

  return;
}
