#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

void CCCCGlobalModes_Reductions(CCTK_ARGUMENTS)
{


  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int varindex = -1;
  int ierr = 0;
  int cf = 1.0;
  int reduction_handle;
  int vartype;

  CCTK_REAL tiny;

  CCTK_REAL sym_factor1,sym_factor2,sym_factor3;

  tiny = 1.0e-20;

  int do_stuff = ( (((cctk_iteration) % compute_every) == 0) || cctk_iteration == 0) &&
    (compute_every != -1) &&
    (cctk_time >= (start_time*2.03e2)) && (do_shibata || do_saijo || do_CoM || do_qlm);

  if(!do_stuff) return;


  sym_factor1 = 1.0e0;
  sym_factor2 = 1.0e0;
  sym_factor3 = 1.0e0;


  //d3x = cctk_delta_space[0]*cctk_delta_space[1]*cctk_delta_space[2];

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


  if(do_saijo||do_CoM||do_P) {

    varindex = CCTK_VarIndex("CCCCGlobalModes::dMass");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)total_mass, 1, varindex);
    assert(!ierr);

    *total_mass = sym_factor1 * (*total_mass);

    if(isnan(*total_mass)) {
      CCTK_WARN(NaN_Warnlevel,"Total Mass is NAN");
    }

  }


  if(do_saijo) {

    varindex = CCTK_VarIndex("CCCCGlobalModes::ddi_re");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)di_re, 1, varindex);
    assert(!ierr);


    varindex = CCTK_VarIndex("CCCCGlobalModes::ddi_im");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)di_im, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::dquad_re");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)quad_re, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::dquad_im");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)quad_im, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::dsextu_re");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)sextu_re, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::dsextu_im");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)sextu_im, 1, varindex);
    assert(!ierr);


    *di_re = sym_factor2 *(*di_re)/(*total_mass);
    *di_im = sym_factor2 *(*di_im)/(*total_mass);
    *quad_re = sym_factor2 *(*quad_re)/(*total_mass);
    *quad_im = sym_factor2 *(*quad_im)/(*total_mass);
    *sextu_re = sym_factor2 *(*sextu_re)/(*total_mass);
    *sextu_im = sym_factor2 *(*sextu_im)/(*total_mass);


  }


  if (do_shibata) {

    varindex = CCTK_VarIndex("CCCCGlobalModes::dIxx");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)Ixx, 1, varindex);
    assert(!ierr);


    varindex = CCTK_VarIndex("CCCCGlobalModes::dIyy");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)Iyy, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::dIxy");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)Ixy, 1, varindex);
    assert(!ierr);

    *Ixx = sym_factor2 * (*Ixx);
    *Iyy = sym_factor2 * (*Iyy);
    *Ixy = sym_factor2 * (*Ixy);
    
    *eta_plus = (*Ixx - *Iyy)/(*Ixx + *Iyy + tiny);
    *eta_cross = (2.0e0*(*Ixy))/(*Ixx + *Iyy + tiny);
    *eta = sqrt((*eta_plus)*(*eta_plus) + (*eta_cross)*(*eta_cross));


  }

  if (do_CoM) {


    varindex = CCTK_VarIndex("CCCCGlobalModes::dMx");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)Mx, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::dMy");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)My, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::dMz");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)Mz, 1, varindex);
    assert(!ierr);

    *Mx = sym_factor2*(*Mx)/(*total_mass);
    *My = sym_factor2*(*My)/(*total_mass);
    *Mz = sym_factor3*(*Mz)/(*total_mass);
    
    *Mr = sqrt( (*Mx)*(*Mx)+(*My)*(*My)+(*Mz)*(*Mz) );

    if(do_CoM_abort) {
      if(*Mr > do_CoM_abort_threshold) {
       CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
		  "CoM violation above termination threshold! %15.6E %15.6E", *Mr, do_CoM_abort_threshold);
      }
    }
    
    *Mass = *total_mass;
  }

  if (do_P) {


    varindex = CCTK_VarIndex("CCCCGlobalModes::dPx");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)Px, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::dPy");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)Py, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::dPz");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)Pz, 1, varindex);
    assert(!ierr);

    *Px = sym_factor2*(*Px)/(*total_mass);
    *Py = sym_factor2*(*Py)/(*total_mass);
    *Pz = sym_factor3*(*Pz)/(*total_mass);
    
    *Pr = sqrt( (*Px)*(*Px)+(*Py)*(*Py)+(*Pz)*(*Pz) );
  }

  if (do_qlm) {

    varindex = CCTK_VarIndex("CCCCGlobalModes::qlm_integrand_re");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)qlm_re, 1, varindex);
    assert(!ierr);

    varindex = CCTK_VarIndex("CCCCGlobalModes::qlm_integrand_im");
    assert(varindex>=0);
    ierr = CCTK_Reduce(cctkGH, -1, reduction_handle, 
		       1, CCTK_VARIABLE_REAL, (void *)qlm_im, 1, varindex);
    assert(!ierr);

    if (sym_factor1 != 1 || sym_factor2 != 1 || sym_factor3 != 1)
    {
       CCTK_WARN(1, "Computation of qlm in presence of symmetries currently not supported!");
    }
  }


  return;
}
