/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ML_BSSN_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At11"),  CCTK_VarIndex("ML_BSSN::At11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At12"),  CCTK_VarIndex("ML_BSSN::At12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At13"),  CCTK_VarIndex("ML_BSSN::At13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At22"),  CCTK_VarIndex("ML_BSSN::At22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At23"),  CCTK_VarIndex("ML_BSSN::At23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::At33"),  CCTK_VarIndex("ML_BSSN::At33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::A"),  CCTK_VarIndex("ML_BSSN::Arhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::B1"),  CCTK_VarIndex("ML_BSSN::B1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::B2"),  CCTK_VarIndex("ML_BSSN::B2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::B3"),  CCTK_VarIndex("ML_BSSN::B3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::Xt1"),  CCTK_VarIndex("ML_BSSN::Xt1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::Xt2"),  CCTK_VarIndex("ML_BSSN::Xt2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::Xt3"),  CCTK_VarIndex("ML_BSSN::Xt3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::alpha"),  CCTK_VarIndex("ML_BSSN::alpharhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::phi"),  CCTK_VarIndex("ML_BSSN::phirhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt11"),  CCTK_VarIndex("ML_BSSN::gt11rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt12"),  CCTK_VarIndex("ML_BSSN::gt12rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt13"),  CCTK_VarIndex("ML_BSSN::gt13rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt22"),  CCTK_VarIndex("ML_BSSN::gt22rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt23"),  CCTK_VarIndex("ML_BSSN::gt23rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::gt33"),  CCTK_VarIndex("ML_BSSN::gt33rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::beta1"),  CCTK_VarIndex("ML_BSSN::beta1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::beta2"),  CCTK_VarIndex("ML_BSSN::beta2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::beta3"),  CCTK_VarIndex("ML_BSSN::beta3rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("ML_BSSN::trK"),  CCTK_VarIndex("ML_BSSN::trKrhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
