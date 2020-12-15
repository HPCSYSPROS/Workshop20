#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>
#include <cctk_Functions.h>
#include <cctk_Faces.h>
#include <Symmetry.h>
#include <loopcontrol.h>
#include <time.h>

void PNSHelper_MoLRegister(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr=0;

  CCTK_INFO("Registering metric/shift/lapse variables with MoL.");

  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::betax"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::betay"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::betaz"));


  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gxx"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gxy"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gxz"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gyy"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gyz"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::gzz"));

  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kxx"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kxy"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kxz"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kyy"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kyz"));
  ierr += MoLRegisterConstrained(CCTK_VarIndex("ADMBase::kzz"));

  if (ierr) CCTK_WARN(0,"Problems registering variables with MoL");

}
