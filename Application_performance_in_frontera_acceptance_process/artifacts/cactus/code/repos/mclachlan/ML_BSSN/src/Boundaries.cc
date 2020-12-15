/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"


/* the boundary treatment is split into 3 steps:    */
/* 1. excision                                      */
/* 2. symmetries                                    */
/* 3. "other" boundary conditions, e.g. radiative */

/* to simplify scheduling and testing, the 3 steps  */
/* are currently applied in separate functions      */


extern "C" void ML_BSSN_CheckBoundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void ML_BSSN_SelectBoundConds(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(ML_log_confac_bound, "none"  ) ||
      CCTK_EQUALS(ML_log_confac_bound, "static") ||
      CCTK_EQUALS(ML_log_confac_bound, "flat"  ) ||
      CCTK_EQUALS(ML_log_confac_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::ML_log_confac", ML_log_confac_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_log_confac_bound BC for ML_BSSN::ML_log_confac!");
  }
  
  if (CCTK_EQUALS(ML_metric_bound, "none"  ) ||
      CCTK_EQUALS(ML_metric_bound, "static") ||
      CCTK_EQUALS(ML_metric_bound, "flat"  ) ||
      CCTK_EQUALS(ML_metric_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::ML_metric", ML_metric_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_metric_bound BC for ML_BSSN::ML_metric!");
  }
  
  if (CCTK_EQUALS(ML_Gamma_bound, "none"  ) ||
      CCTK_EQUALS(ML_Gamma_bound, "static") ||
      CCTK_EQUALS(ML_Gamma_bound, "flat"  ) ||
      CCTK_EQUALS(ML_Gamma_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::ML_Gamma", ML_Gamma_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_Gamma_bound BC for ML_BSSN::ML_Gamma!");
  }
  
  if (CCTK_EQUALS(ML_trace_curv_bound, "none"  ) ||
      CCTK_EQUALS(ML_trace_curv_bound, "static") ||
      CCTK_EQUALS(ML_trace_curv_bound, "flat"  ) ||
      CCTK_EQUALS(ML_trace_curv_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::ML_trace_curv", ML_trace_curv_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_trace_curv_bound BC for ML_BSSN::ML_trace_curv!");
  }
  
  if (CCTK_EQUALS(ML_curv_bound, "none"  ) ||
      CCTK_EQUALS(ML_curv_bound, "static") ||
      CCTK_EQUALS(ML_curv_bound, "flat"  ) ||
      CCTK_EQUALS(ML_curv_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::ML_curv", ML_curv_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_curv_bound BC for ML_BSSN::ML_curv!");
  }
  
  if (CCTK_EQUALS(ML_lapse_bound, "none"  ) ||
      CCTK_EQUALS(ML_lapse_bound, "static") ||
      CCTK_EQUALS(ML_lapse_bound, "flat"  ) ||
      CCTK_EQUALS(ML_lapse_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::ML_lapse", ML_lapse_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_lapse_bound BC for ML_BSSN::ML_lapse!");
  }
  
  if (CCTK_EQUALS(ML_dtlapse_bound, "none"  ) ||
      CCTK_EQUALS(ML_dtlapse_bound, "static") ||
      CCTK_EQUALS(ML_dtlapse_bound, "flat"  ) ||
      CCTK_EQUALS(ML_dtlapse_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::ML_dtlapse", ML_dtlapse_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_dtlapse_bound BC for ML_BSSN::ML_dtlapse!");
  }
  
  if (CCTK_EQUALS(ML_shift_bound, "none"  ) ||
      CCTK_EQUALS(ML_shift_bound, "static") ||
      CCTK_EQUALS(ML_shift_bound, "flat"  ) ||
      CCTK_EQUALS(ML_shift_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::ML_shift", ML_shift_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_shift_bound BC for ML_BSSN::ML_shift!");
  }
  
  if (CCTK_EQUALS(ML_dtshift_bound, "none"  ) ||
      CCTK_EQUALS(ML_dtshift_bound, "static") ||
      CCTK_EQUALS(ML_dtshift_bound, "flat"  ) ||
      CCTK_EQUALS(ML_dtshift_bound, "zero"  ) )
  {
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::ML_dtshift", ML_dtshift_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register ML_dtshift_bound BC for ML_BSSN::ML_dtshift!");
  }
  
  if (CCTK_EQUALS(phi_bound, "none"  ) ||
      CCTK_EQUALS(phi_bound, "static") ||
      CCTK_EQUALS(phi_bound, "flat"  ) ||
      CCTK_EQUALS(phi_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::phi", phi_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register phi_bound BC for ML_BSSN::phi!");
  }
  
  if (CCTK_EQUALS(gt11_bound, "none"  ) ||
      CCTK_EQUALS(gt11_bound, "static") ||
      CCTK_EQUALS(gt11_bound, "flat"  ) ||
      CCTK_EQUALS(gt11_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::gt11", gt11_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt11_bound BC for ML_BSSN::gt11!");
  }
  
  if (CCTK_EQUALS(gt12_bound, "none"  ) ||
      CCTK_EQUALS(gt12_bound, "static") ||
      CCTK_EQUALS(gt12_bound, "flat"  ) ||
      CCTK_EQUALS(gt12_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::gt12", gt12_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt12_bound BC for ML_BSSN::gt12!");
  }
  
  if (CCTK_EQUALS(gt13_bound, "none"  ) ||
      CCTK_EQUALS(gt13_bound, "static") ||
      CCTK_EQUALS(gt13_bound, "flat"  ) ||
      CCTK_EQUALS(gt13_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::gt13", gt13_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt13_bound BC for ML_BSSN::gt13!");
  }
  
  if (CCTK_EQUALS(gt22_bound, "none"  ) ||
      CCTK_EQUALS(gt22_bound, "static") ||
      CCTK_EQUALS(gt22_bound, "flat"  ) ||
      CCTK_EQUALS(gt22_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::gt22", gt22_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt22_bound BC for ML_BSSN::gt22!");
  }
  
  if (CCTK_EQUALS(gt23_bound, "none"  ) ||
      CCTK_EQUALS(gt23_bound, "static") ||
      CCTK_EQUALS(gt23_bound, "flat"  ) ||
      CCTK_EQUALS(gt23_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::gt23", gt23_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt23_bound BC for ML_BSSN::gt23!");
  }
  
  if (CCTK_EQUALS(gt33_bound, "none"  ) ||
      CCTK_EQUALS(gt33_bound, "static") ||
      CCTK_EQUALS(gt33_bound, "flat"  ) ||
      CCTK_EQUALS(gt33_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::gt33", gt33_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register gt33_bound BC for ML_BSSN::gt33!");
  }
  
  if (CCTK_EQUALS(Xt1_bound, "none"  ) ||
      CCTK_EQUALS(Xt1_bound, "static") ||
      CCTK_EQUALS(Xt1_bound, "flat"  ) ||
      CCTK_EQUALS(Xt1_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::Xt1", Xt1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Xt1_bound BC for ML_BSSN::Xt1!");
  }
  
  if (CCTK_EQUALS(Xt2_bound, "none"  ) ||
      CCTK_EQUALS(Xt2_bound, "static") ||
      CCTK_EQUALS(Xt2_bound, "flat"  ) ||
      CCTK_EQUALS(Xt2_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::Xt2", Xt2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Xt2_bound BC for ML_BSSN::Xt2!");
  }
  
  if (CCTK_EQUALS(Xt3_bound, "none"  ) ||
      CCTK_EQUALS(Xt3_bound, "static") ||
      CCTK_EQUALS(Xt3_bound, "flat"  ) ||
      CCTK_EQUALS(Xt3_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::Xt3", Xt3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register Xt3_bound BC for ML_BSSN::Xt3!");
  }
  
  if (CCTK_EQUALS(trK_bound, "none"  ) ||
      CCTK_EQUALS(trK_bound, "static") ||
      CCTK_EQUALS(trK_bound, "flat"  ) ||
      CCTK_EQUALS(trK_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::trK", trK_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register trK_bound BC for ML_BSSN::trK!");
  }
  
  if (CCTK_EQUALS(At11_bound, "none"  ) ||
      CCTK_EQUALS(At11_bound, "static") ||
      CCTK_EQUALS(At11_bound, "flat"  ) ||
      CCTK_EQUALS(At11_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::At11", At11_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At11_bound BC for ML_BSSN::At11!");
  }
  
  if (CCTK_EQUALS(At12_bound, "none"  ) ||
      CCTK_EQUALS(At12_bound, "static") ||
      CCTK_EQUALS(At12_bound, "flat"  ) ||
      CCTK_EQUALS(At12_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::At12", At12_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At12_bound BC for ML_BSSN::At12!");
  }
  
  if (CCTK_EQUALS(At13_bound, "none"  ) ||
      CCTK_EQUALS(At13_bound, "static") ||
      CCTK_EQUALS(At13_bound, "flat"  ) ||
      CCTK_EQUALS(At13_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::At13", At13_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At13_bound BC for ML_BSSN::At13!");
  }
  
  if (CCTK_EQUALS(At22_bound, "none"  ) ||
      CCTK_EQUALS(At22_bound, "static") ||
      CCTK_EQUALS(At22_bound, "flat"  ) ||
      CCTK_EQUALS(At22_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::At22", At22_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At22_bound BC for ML_BSSN::At22!");
  }
  
  if (CCTK_EQUALS(At23_bound, "none"  ) ||
      CCTK_EQUALS(At23_bound, "static") ||
      CCTK_EQUALS(At23_bound, "flat"  ) ||
      CCTK_EQUALS(At23_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::At23", At23_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At23_bound BC for ML_BSSN::At23!");
  }
  
  if (CCTK_EQUALS(At33_bound, "none"  ) ||
      CCTK_EQUALS(At33_bound, "static") ||
      CCTK_EQUALS(At33_bound, "flat"  ) ||
      CCTK_EQUALS(At33_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::At33", At33_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register At33_bound BC for ML_BSSN::At33!");
  }
  
  if (CCTK_EQUALS(alpha_bound, "none"  ) ||
      CCTK_EQUALS(alpha_bound, "static") ||
      CCTK_EQUALS(alpha_bound, "flat"  ) ||
      CCTK_EQUALS(alpha_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::alpha", alpha_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register alpha_bound BC for ML_BSSN::alpha!");
  }
  
  if (CCTK_EQUALS(A_bound, "none"  ) ||
      CCTK_EQUALS(A_bound, "static") ||
      CCTK_EQUALS(A_bound, "flat"  ) ||
      CCTK_EQUALS(A_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::A", A_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register A_bound BC for ML_BSSN::A!");
  }
  
  if (CCTK_EQUALS(beta1_bound, "none"  ) ||
      CCTK_EQUALS(beta1_bound, "static") ||
      CCTK_EQUALS(beta1_bound, "flat"  ) ||
      CCTK_EQUALS(beta1_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::beta1", beta1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta1_bound BC for ML_BSSN::beta1!");
  }
  
  if (CCTK_EQUALS(beta2_bound, "none"  ) ||
      CCTK_EQUALS(beta2_bound, "static") ||
      CCTK_EQUALS(beta2_bound, "flat"  ) ||
      CCTK_EQUALS(beta2_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::beta2", beta2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta2_bound BC for ML_BSSN::beta2!");
  }
  
  if (CCTK_EQUALS(beta3_bound, "none"  ) ||
      CCTK_EQUALS(beta3_bound, "static") ||
      CCTK_EQUALS(beta3_bound, "flat"  ) ||
      CCTK_EQUALS(beta3_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::beta3", beta3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register beta3_bound BC for ML_BSSN::beta3!");
  }
  
  if (CCTK_EQUALS(B1_bound, "none"  ) ||
      CCTK_EQUALS(B1_bound, "static") ||
      CCTK_EQUALS(B1_bound, "flat"  ) ||
      CCTK_EQUALS(B1_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::B1", B1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register B1_bound BC for ML_BSSN::B1!");
  }
  
  if (CCTK_EQUALS(B2_bound, "none"  ) ||
      CCTK_EQUALS(B2_bound, "static") ||
      CCTK_EQUALS(B2_bound, "flat"  ) ||
      CCTK_EQUALS(B2_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::B2", B2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register B2_bound BC for ML_BSSN::B2!");
  }
  
  if (CCTK_EQUALS(B3_bound, "none"  ) ||
      CCTK_EQUALS(B3_bound, "static") ||
      CCTK_EQUALS(B3_bound, "flat"  ) ||
      CCTK_EQUALS(B3_bound, "zero"  ) )
  {
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "ML_BSSN::B3", B3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register B3_bound BC for ML_BSSN::B3!");
  }
  
  if (CCTK_EQUALS(ML_log_confac_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_log_confac_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_log_confac_bound < 0) handle_ML_log_confac_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_log_confac_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_log_confac_bound , ML_log_confac_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_log_confac_bound ,ML_log_confac_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_log_confac_bound, 
                      "ML_BSSN::ML_log_confac", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::ML_log_confac!");
  
  }
  
  if (CCTK_EQUALS(ML_metric_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_metric_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_metric_bound < 0) handle_ML_metric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_metric_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_metric_bound , ML_metric_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_metric_bound ,ML_metric_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_metric_bound, 
                      "ML_BSSN::ML_metric", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::ML_metric!");
  
  }
  
  if (CCTK_EQUALS(ML_Gamma_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_Gamma_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_Gamma_bound < 0) handle_ML_Gamma_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_Gamma_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_Gamma_bound , ML_Gamma_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_Gamma_bound ,ML_Gamma_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_Gamma_bound, 
                      "ML_BSSN::ML_Gamma", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::ML_Gamma!");
  
  }
  
  if (CCTK_EQUALS(ML_trace_curv_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_trace_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_trace_curv_bound < 0) handle_ML_trace_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_trace_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_trace_curv_bound , ML_trace_curv_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_trace_curv_bound ,ML_trace_curv_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_trace_curv_bound, 
                      "ML_BSSN::ML_trace_curv", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::ML_trace_curv!");
  
  }
  
  if (CCTK_EQUALS(ML_curv_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_curv_bound < 0) handle_ML_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_curv_bound , ML_curv_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_curv_bound ,ML_curv_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_curv_bound, 
                      "ML_BSSN::ML_curv", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::ML_curv!");
  
  }
  
  if (CCTK_EQUALS(ML_lapse_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_lapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_lapse_bound < 0) handle_ML_lapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_lapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_lapse_bound , ML_lapse_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_lapse_bound ,ML_lapse_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_lapse_bound, 
                      "ML_BSSN::ML_lapse", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::ML_lapse!");
  
  }
  
  if (CCTK_EQUALS(ML_dtlapse_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_dtlapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_dtlapse_bound < 0) handle_ML_dtlapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_dtlapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_dtlapse_bound , ML_dtlapse_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_dtlapse_bound ,ML_dtlapse_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_dtlapse_bound, 
                      "ML_BSSN::ML_dtlapse", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::ML_dtlapse!");
  
  }
  
  if (CCTK_EQUALS(ML_shift_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_shift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_shift_bound < 0) handle_ML_shift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_shift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_shift_bound , ML_shift_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_shift_bound ,ML_shift_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_shift_bound, 
                      "ML_BSSN::ML_shift", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::ML_shift!");
  
  }
  
  if (CCTK_EQUALS(ML_dtshift_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_ML_dtshift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_dtshift_bound < 0) handle_ML_dtshift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_dtshift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_dtshift_bound , ML_dtshift_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_ML_dtshift_bound ,ML_dtshift_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_dtshift_bound, 
                      "ML_BSSN::ML_dtshift", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::ML_dtshift!");
  
  }
  
  if (CCTK_EQUALS(phi_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_phi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phi_bound < 0) handle_phi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phi_bound , phi_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_phi_bound ,phi_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phi_bound, 
                      "ML_BSSN::phi", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::phi!");
  
  }
  
  if (CCTK_EQUALS(gt11_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt11_bound < 0) handle_gt11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt11_bound , gt11_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt11_bound ,gt11_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt11_bound, 
                      "ML_BSSN::gt11", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::gt11!");
  
  }
  
  if (CCTK_EQUALS(gt12_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt12_bound < 0) handle_gt12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt12_bound , gt12_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt12_bound ,gt12_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt12_bound, 
                      "ML_BSSN::gt12", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::gt12!");
  
  }
  
  if (CCTK_EQUALS(gt13_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt13_bound < 0) handle_gt13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt13_bound , gt13_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt13_bound ,gt13_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt13_bound, 
                      "ML_BSSN::gt13", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::gt13!");
  
  }
  
  if (CCTK_EQUALS(gt22_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt22_bound < 0) handle_gt22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt22_bound , gt22_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt22_bound ,gt22_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt22_bound, 
                      "ML_BSSN::gt22", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::gt22!");
  
  }
  
  if (CCTK_EQUALS(gt23_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt23_bound < 0) handle_gt23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt23_bound , gt23_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt23_bound ,gt23_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt23_bound, 
                      "ML_BSSN::gt23", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::gt23!");
  
  }
  
  if (CCTK_EQUALS(gt33_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_gt33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt33_bound < 0) handle_gt33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt33_bound , gt33_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_gt33_bound ,gt33_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt33_bound, 
                      "ML_BSSN::gt33", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::gt33!");
  
  }
  
  if (CCTK_EQUALS(Xt1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Xt1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt1_bound < 0) handle_Xt1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt1_bound , Xt1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Xt1_bound ,Xt1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt1_bound, 
                      "ML_BSSN::Xt1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::Xt1!");
  
  }
  
  if (CCTK_EQUALS(Xt2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Xt2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt2_bound < 0) handle_Xt2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt2_bound , Xt2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Xt2_bound ,Xt2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt2_bound, 
                      "ML_BSSN::Xt2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::Xt2!");
  
  }
  
  if (CCTK_EQUALS(Xt3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_Xt3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt3_bound < 0) handle_Xt3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt3_bound , Xt3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_Xt3_bound ,Xt3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt3_bound, 
                      "ML_BSSN::Xt3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::Xt3!");
  
  }
  
  if (CCTK_EQUALS(trK_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_trK_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_trK_bound < 0) handle_trK_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_trK_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_trK_bound , trK_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_trK_bound ,trK_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_trK_bound, 
                      "ML_BSSN::trK", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::trK!");
  
  }
  
  if (CCTK_EQUALS(At11_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At11_bound < 0) handle_At11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At11_bound , At11_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At11_bound ,At11_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At11_bound, 
                      "ML_BSSN::At11", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::At11!");
  
  }
  
  if (CCTK_EQUALS(At12_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At12_bound < 0) handle_At12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At12_bound , At12_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At12_bound ,At12_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At12_bound, 
                      "ML_BSSN::At12", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::At12!");
  
  }
  
  if (CCTK_EQUALS(At13_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At13_bound < 0) handle_At13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At13_bound , At13_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At13_bound ,At13_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At13_bound, 
                      "ML_BSSN::At13", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::At13!");
  
  }
  
  if (CCTK_EQUALS(At22_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At22_bound < 0) handle_At22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At22_bound , At22_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At22_bound ,At22_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At22_bound, 
                      "ML_BSSN::At22", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::At22!");
  
  }
  
  if (CCTK_EQUALS(At23_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At23_bound < 0) handle_At23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At23_bound , At23_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At23_bound ,At23_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At23_bound, 
                      "ML_BSSN::At23", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::At23!");
  
  }
  
  if (CCTK_EQUALS(At33_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_At33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At33_bound < 0) handle_At33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At33_bound , At33_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_At33_bound ,At33_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At33_bound, 
                      "ML_BSSN::At33", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::At33!");
  
  }
  
  if (CCTK_EQUALS(alpha_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_alpha_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_alpha_bound < 0) handle_alpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_alpha_bound , alpha_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_alpha_bound ,alpha_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_bound, 
                      "ML_BSSN::alpha", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::alpha!");
  
  }
  
  if (CCTK_EQUALS(A_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_A_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A_bound < 0) handle_A_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A_bound , A_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_A_bound ,A_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A_bound, 
                      "ML_BSSN::A", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::A!");
  
  }
  
  if (CCTK_EQUALS(beta1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta1_bound < 0) handle_beta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta1_bound , beta1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta1_bound ,beta1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta1_bound, 
                      "ML_BSSN::beta1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::beta1!");
  
  }
  
  if (CCTK_EQUALS(beta2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta2_bound < 0) handle_beta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta2_bound , beta2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta2_bound ,beta2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta2_bound, 
                      "ML_BSSN::beta2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::beta2!");
  
  }
  
  if (CCTK_EQUALS(beta3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_beta3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta3_bound < 0) handle_beta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta3_bound , beta3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_beta3_bound ,beta3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta3_bound, 
                      "ML_BSSN::beta3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::beta3!");
  
  }
  
  if (CCTK_EQUALS(B1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B1_bound < 0) handle_B1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B1_bound , B1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B1_bound ,B1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B1_bound, 
                      "ML_BSSN::B1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::B1!");
  
  }
  
  if (CCTK_EQUALS(B2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B2_bound < 0) handle_B2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B2_bound , B2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B2_bound ,B2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B2_bound, 
                      "ML_BSSN::B2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::B2!");
  
  }
  
  if (CCTK_EQUALS(B3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_B3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B3_bound < 0) handle_B3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B3_bound , B3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_B3_bound ,B3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B3_bound, 
                      "ML_BSSN::B3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for ML_BSSN::B3!");
  
  }
  
  if (CCTK_EQUALS(ML_log_confac_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_log_confac_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_log_confac_bound < 0) handle_ML_log_confac_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_log_confac_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_log_confac_bound ,ML_log_confac_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_log_confac_bound, 
                      "ML_BSSN::ML_log_confac", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_BSSN::ML_log_confac!");
  
  }
  
  if (CCTK_EQUALS(ML_metric_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_metric_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_metric_bound < 0) handle_ML_metric_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_metric_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_metric_bound ,ML_metric_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_metric_bound, 
                      "ML_BSSN::ML_metric", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_BSSN::ML_metric!");
  
  }
  
  if (CCTK_EQUALS(ML_Gamma_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_Gamma_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_Gamma_bound < 0) handle_ML_Gamma_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_Gamma_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_Gamma_bound ,ML_Gamma_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_Gamma_bound, 
                      "ML_BSSN::ML_Gamma", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_BSSN::ML_Gamma!");
  
  }
  
  if (CCTK_EQUALS(ML_trace_curv_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_trace_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_trace_curv_bound < 0) handle_ML_trace_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_trace_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_trace_curv_bound ,ML_trace_curv_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_trace_curv_bound, 
                      "ML_BSSN::ML_trace_curv", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_BSSN::ML_trace_curv!");
  
  }
  
  if (CCTK_EQUALS(ML_curv_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_curv_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_curv_bound < 0) handle_ML_curv_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_curv_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_curv_bound ,ML_curv_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_curv_bound, 
                      "ML_BSSN::ML_curv", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_BSSN::ML_curv!");
  
  }
  
  if (CCTK_EQUALS(ML_lapse_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_lapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_lapse_bound < 0) handle_ML_lapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_lapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_lapse_bound ,ML_lapse_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_lapse_bound, 
                      "ML_BSSN::ML_lapse", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_BSSN::ML_lapse!");
  
  }
  
  if (CCTK_EQUALS(ML_dtlapse_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_dtlapse_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_dtlapse_bound < 0) handle_ML_dtlapse_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_dtlapse_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_dtlapse_bound ,ML_dtlapse_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_dtlapse_bound, 
                      "ML_BSSN::ML_dtlapse", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_BSSN::ML_dtlapse!");
  
  }
  
  if (CCTK_EQUALS(ML_shift_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_shift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_shift_bound < 0) handle_ML_shift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_shift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_shift_bound ,ML_shift_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_shift_bound, 
                      "ML_BSSN::ML_shift", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_BSSN::ML_shift!");
  
  }
  
  if (CCTK_EQUALS(ML_dtshift_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_ML_dtshift_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_ML_dtshift_bound < 0) handle_ML_dtshift_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_ML_dtshift_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_ML_dtshift_bound ,ML_dtshift_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_ML_dtshift_bound, 
                      "ML_BSSN::ML_dtshift", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for ML_BSSN::ML_dtshift!");
  
  }
  
  if (CCTK_EQUALS(phi_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_phi_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_phi_bound < 0) handle_phi_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_phi_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_phi_bound ,phi_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_phi_bound, 
                      "ML_BSSN::phi", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::phi!");
  
  }
  
  if (CCTK_EQUALS(gt11_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt11_bound < 0) handle_gt11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt11_bound ,gt11_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt11_bound, 
                      "ML_BSSN::gt11", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::gt11!");
  
  }
  
  if (CCTK_EQUALS(gt12_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt12_bound < 0) handle_gt12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt12_bound ,gt12_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt12_bound, 
                      "ML_BSSN::gt12", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::gt12!");
  
  }
  
  if (CCTK_EQUALS(gt13_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt13_bound < 0) handle_gt13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt13_bound ,gt13_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt13_bound, 
                      "ML_BSSN::gt13", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::gt13!");
  
  }
  
  if (CCTK_EQUALS(gt22_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt22_bound < 0) handle_gt22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt22_bound ,gt22_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt22_bound, 
                      "ML_BSSN::gt22", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::gt22!");
  
  }
  
  if (CCTK_EQUALS(gt23_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt23_bound < 0) handle_gt23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt23_bound ,gt23_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt23_bound, 
                      "ML_BSSN::gt23", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::gt23!");
  
  }
  
  if (CCTK_EQUALS(gt33_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_gt33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_gt33_bound < 0) handle_gt33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_gt33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_gt33_bound ,gt33_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_gt33_bound, 
                      "ML_BSSN::gt33", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::gt33!");
  
  }
  
  if (CCTK_EQUALS(Xt1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Xt1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt1_bound < 0) handle_Xt1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt1_bound ,Xt1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt1_bound, 
                      "ML_BSSN::Xt1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::Xt1!");
  
  }
  
  if (CCTK_EQUALS(Xt2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Xt2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt2_bound < 0) handle_Xt2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt2_bound ,Xt2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt2_bound, 
                      "ML_BSSN::Xt2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::Xt2!");
  
  }
  
  if (CCTK_EQUALS(Xt3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_Xt3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_Xt3_bound < 0) handle_Xt3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_Xt3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_Xt3_bound ,Xt3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_Xt3_bound, 
                      "ML_BSSN::Xt3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::Xt3!");
  
  }
  
  if (CCTK_EQUALS(trK_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_trK_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_trK_bound < 0) handle_trK_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_trK_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_trK_bound ,trK_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_trK_bound, 
                      "ML_BSSN::trK", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::trK!");
  
  }
  
  if (CCTK_EQUALS(At11_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At11_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At11_bound < 0) handle_At11_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At11_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At11_bound ,At11_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At11_bound, 
                      "ML_BSSN::At11", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::At11!");
  
  }
  
  if (CCTK_EQUALS(At12_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At12_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At12_bound < 0) handle_At12_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At12_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At12_bound ,At12_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At12_bound, 
                      "ML_BSSN::At12", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::At12!");
  
  }
  
  if (CCTK_EQUALS(At13_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At13_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At13_bound < 0) handle_At13_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At13_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At13_bound ,At13_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At13_bound, 
                      "ML_BSSN::At13", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::At13!");
  
  }
  
  if (CCTK_EQUALS(At22_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At22_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At22_bound < 0) handle_At22_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At22_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At22_bound ,At22_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At22_bound, 
                      "ML_BSSN::At22", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::At22!");
  
  }
  
  if (CCTK_EQUALS(At23_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At23_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At23_bound < 0) handle_At23_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At23_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At23_bound ,At23_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At23_bound, 
                      "ML_BSSN::At23", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::At23!");
  
  }
  
  if (CCTK_EQUALS(At33_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_At33_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_At33_bound < 0) handle_At33_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_At33_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_At33_bound ,At33_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_At33_bound, 
                      "ML_BSSN::At33", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::At33!");
  
  }
  
  if (CCTK_EQUALS(alpha_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_alpha_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_alpha_bound < 0) handle_alpha_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_alpha_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_alpha_bound ,alpha_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_alpha_bound, 
                      "ML_BSSN::alpha", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::alpha!");
  
  }
  
  if (CCTK_EQUALS(A_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_A_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_A_bound < 0) handle_A_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_A_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_A_bound ,A_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_A_bound, 
                      "ML_BSSN::A", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::A!");
  
  }
  
  if (CCTK_EQUALS(beta1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta1_bound < 0) handle_beta1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta1_bound ,beta1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta1_bound, 
                      "ML_BSSN::beta1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::beta1!");
  
  }
  
  if (CCTK_EQUALS(beta2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta2_bound < 0) handle_beta2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta2_bound ,beta2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta2_bound, 
                      "ML_BSSN::beta2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::beta2!");
  
  }
  
  if (CCTK_EQUALS(beta3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_beta3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_beta3_bound < 0) handle_beta3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_beta3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_beta3_bound ,beta3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_beta3_bound, 
                      "ML_BSSN::beta3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::beta3!");
  
  }
  
  if (CCTK_EQUALS(B1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B1_bound < 0) handle_B1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B1_bound ,B1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B1_bound, 
                      "ML_BSSN::B1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::B1!");
  
  }
  
  if (CCTK_EQUALS(B2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B2_bound < 0) handle_B2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B2_bound ,B2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B2_bound, 
                      "ML_BSSN::B2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::B2!");
  
  }
  
  if (CCTK_EQUALS(B3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_B3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_B3_bound < 0) handle_B3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_B3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_B3_bound ,B3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = Boundary_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_B3_bound, 
                      "ML_BSSN::B3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for ML_BSSN::B3!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#ML_BSSN::ML_log_confac_bound       = "skip"
#$bound$#ML_BSSN::ML_log_confac_bound_speed = 1.0
#$bound$#ML_BSSN::ML_log_confac_bound_limit = 0.0
#$bound$#ML_BSSN::ML_log_confac_bound_scalar = 0.0

#$bound$#ML_BSSN::ML_metric_bound       = "skip"
#$bound$#ML_BSSN::ML_metric_bound_speed = 1.0
#$bound$#ML_BSSN::ML_metric_bound_limit = 0.0
#$bound$#ML_BSSN::ML_metric_bound_scalar = 0.0

#$bound$#ML_BSSN::ML_Gamma_bound       = "skip"
#$bound$#ML_BSSN::ML_Gamma_bound_speed = 1.0
#$bound$#ML_BSSN::ML_Gamma_bound_limit = 0.0
#$bound$#ML_BSSN::ML_Gamma_bound_scalar = 0.0

#$bound$#ML_BSSN::ML_trace_curv_bound       = "skip"
#$bound$#ML_BSSN::ML_trace_curv_bound_speed = 1.0
#$bound$#ML_BSSN::ML_trace_curv_bound_limit = 0.0
#$bound$#ML_BSSN::ML_trace_curv_bound_scalar = 0.0

#$bound$#ML_BSSN::ML_curv_bound       = "skip"
#$bound$#ML_BSSN::ML_curv_bound_speed = 1.0
#$bound$#ML_BSSN::ML_curv_bound_limit = 0.0
#$bound$#ML_BSSN::ML_curv_bound_scalar = 0.0

#$bound$#ML_BSSN::ML_lapse_bound       = "skip"
#$bound$#ML_BSSN::ML_lapse_bound_speed = 1.0
#$bound$#ML_BSSN::ML_lapse_bound_limit = 0.0
#$bound$#ML_BSSN::ML_lapse_bound_scalar = 0.0

#$bound$#ML_BSSN::ML_dtlapse_bound       = "skip"
#$bound$#ML_BSSN::ML_dtlapse_bound_speed = 1.0
#$bound$#ML_BSSN::ML_dtlapse_bound_limit = 0.0
#$bound$#ML_BSSN::ML_dtlapse_bound_scalar = 0.0

#$bound$#ML_BSSN::ML_shift_bound       = "skip"
#$bound$#ML_BSSN::ML_shift_bound_speed = 1.0
#$bound$#ML_BSSN::ML_shift_bound_limit = 0.0
#$bound$#ML_BSSN::ML_shift_bound_scalar = 0.0

#$bound$#ML_BSSN::ML_dtshift_bound       = "skip"
#$bound$#ML_BSSN::ML_dtshift_bound_speed = 1.0
#$bound$#ML_BSSN::ML_dtshift_bound_limit = 0.0
#$bound$#ML_BSSN::ML_dtshift_bound_scalar = 0.0

#$bound$#ML_BSSN::phi_bound       = "skip"
#$bound$#ML_BSSN::phi_bound_speed = 1.0
#$bound$#ML_BSSN::phi_bound_limit = 0.0
#$bound$#ML_BSSN::phi_bound_scalar = 0.0

#$bound$#ML_BSSN::gt11_bound       = "skip"
#$bound$#ML_BSSN::gt11_bound_speed = 1.0
#$bound$#ML_BSSN::gt11_bound_limit = 0.0
#$bound$#ML_BSSN::gt11_bound_scalar = 0.0

#$bound$#ML_BSSN::gt12_bound       = "skip"
#$bound$#ML_BSSN::gt12_bound_speed = 1.0
#$bound$#ML_BSSN::gt12_bound_limit = 0.0
#$bound$#ML_BSSN::gt12_bound_scalar = 0.0

#$bound$#ML_BSSN::gt13_bound       = "skip"
#$bound$#ML_BSSN::gt13_bound_speed = 1.0
#$bound$#ML_BSSN::gt13_bound_limit = 0.0
#$bound$#ML_BSSN::gt13_bound_scalar = 0.0

#$bound$#ML_BSSN::gt22_bound       = "skip"
#$bound$#ML_BSSN::gt22_bound_speed = 1.0
#$bound$#ML_BSSN::gt22_bound_limit = 0.0
#$bound$#ML_BSSN::gt22_bound_scalar = 0.0

#$bound$#ML_BSSN::gt23_bound       = "skip"
#$bound$#ML_BSSN::gt23_bound_speed = 1.0
#$bound$#ML_BSSN::gt23_bound_limit = 0.0
#$bound$#ML_BSSN::gt23_bound_scalar = 0.0

#$bound$#ML_BSSN::gt33_bound       = "skip"
#$bound$#ML_BSSN::gt33_bound_speed = 1.0
#$bound$#ML_BSSN::gt33_bound_limit = 0.0
#$bound$#ML_BSSN::gt33_bound_scalar = 0.0

#$bound$#ML_BSSN::Xt1_bound       = "skip"
#$bound$#ML_BSSN::Xt1_bound_speed = 1.0
#$bound$#ML_BSSN::Xt1_bound_limit = 0.0
#$bound$#ML_BSSN::Xt1_bound_scalar = 0.0

#$bound$#ML_BSSN::Xt2_bound       = "skip"
#$bound$#ML_BSSN::Xt2_bound_speed = 1.0
#$bound$#ML_BSSN::Xt2_bound_limit = 0.0
#$bound$#ML_BSSN::Xt2_bound_scalar = 0.0

#$bound$#ML_BSSN::Xt3_bound       = "skip"
#$bound$#ML_BSSN::Xt3_bound_speed = 1.0
#$bound$#ML_BSSN::Xt3_bound_limit = 0.0
#$bound$#ML_BSSN::Xt3_bound_scalar = 0.0

#$bound$#ML_BSSN::trK_bound       = "skip"
#$bound$#ML_BSSN::trK_bound_speed = 1.0
#$bound$#ML_BSSN::trK_bound_limit = 0.0
#$bound$#ML_BSSN::trK_bound_scalar = 0.0

#$bound$#ML_BSSN::At11_bound       = "skip"
#$bound$#ML_BSSN::At11_bound_speed = 1.0
#$bound$#ML_BSSN::At11_bound_limit = 0.0
#$bound$#ML_BSSN::At11_bound_scalar = 0.0

#$bound$#ML_BSSN::At12_bound       = "skip"
#$bound$#ML_BSSN::At12_bound_speed = 1.0
#$bound$#ML_BSSN::At12_bound_limit = 0.0
#$bound$#ML_BSSN::At12_bound_scalar = 0.0

#$bound$#ML_BSSN::At13_bound       = "skip"
#$bound$#ML_BSSN::At13_bound_speed = 1.0
#$bound$#ML_BSSN::At13_bound_limit = 0.0
#$bound$#ML_BSSN::At13_bound_scalar = 0.0

#$bound$#ML_BSSN::At22_bound       = "skip"
#$bound$#ML_BSSN::At22_bound_speed = 1.0
#$bound$#ML_BSSN::At22_bound_limit = 0.0
#$bound$#ML_BSSN::At22_bound_scalar = 0.0

#$bound$#ML_BSSN::At23_bound       = "skip"
#$bound$#ML_BSSN::At23_bound_speed = 1.0
#$bound$#ML_BSSN::At23_bound_limit = 0.0
#$bound$#ML_BSSN::At23_bound_scalar = 0.0

#$bound$#ML_BSSN::At33_bound       = "skip"
#$bound$#ML_BSSN::At33_bound_speed = 1.0
#$bound$#ML_BSSN::At33_bound_limit = 0.0
#$bound$#ML_BSSN::At33_bound_scalar = 0.0

#$bound$#ML_BSSN::alpha_bound       = "skip"
#$bound$#ML_BSSN::alpha_bound_speed = 1.0
#$bound$#ML_BSSN::alpha_bound_limit = 0.0
#$bound$#ML_BSSN::alpha_bound_scalar = 0.0

#$bound$#ML_BSSN::A_bound       = "skip"
#$bound$#ML_BSSN::A_bound_speed = 1.0
#$bound$#ML_BSSN::A_bound_limit = 0.0
#$bound$#ML_BSSN::A_bound_scalar = 0.0

#$bound$#ML_BSSN::beta1_bound       = "skip"
#$bound$#ML_BSSN::beta1_bound_speed = 1.0
#$bound$#ML_BSSN::beta1_bound_limit = 0.0
#$bound$#ML_BSSN::beta1_bound_scalar = 0.0

#$bound$#ML_BSSN::beta2_bound       = "skip"
#$bound$#ML_BSSN::beta2_bound_speed = 1.0
#$bound$#ML_BSSN::beta2_bound_limit = 0.0
#$bound$#ML_BSSN::beta2_bound_scalar = 0.0

#$bound$#ML_BSSN::beta3_bound       = "skip"
#$bound$#ML_BSSN::beta3_bound_speed = 1.0
#$bound$#ML_BSSN::beta3_bound_limit = 0.0
#$bound$#ML_BSSN::beta3_bound_scalar = 0.0

#$bound$#ML_BSSN::B1_bound       = "skip"
#$bound$#ML_BSSN::B1_bound_speed = 1.0
#$bound$#ML_BSSN::B1_bound_limit = 0.0
#$bound$#ML_BSSN::B1_bound_scalar = 0.0

#$bound$#ML_BSSN::B2_bound       = "skip"
#$bound$#ML_BSSN::B2_bound_speed = 1.0
#$bound$#ML_BSSN::B2_bound_limit = 0.0
#$bound$#ML_BSSN::B2_bound_scalar = 0.0

#$bound$#ML_BSSN::B3_bound       = "skip"
#$bound$#ML_BSSN::B3_bound_speed = 1.0
#$bound$#ML_BSSN::B3_bound_limit = 0.0
#$bound$#ML_BSSN::B3_bound_scalar = 0.0

*/

