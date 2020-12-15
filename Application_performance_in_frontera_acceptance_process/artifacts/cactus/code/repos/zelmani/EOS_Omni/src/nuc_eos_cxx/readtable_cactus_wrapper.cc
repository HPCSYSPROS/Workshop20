#include <stdlib.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
//#include "cctk_Functions.h"
//#include "nuc_eos.hh"

extern "C"
void nuc_eos_C_ReadTable(const char *nuceos_table_name, const cGH* cctkGH );

extern "C"
void nuc_eos_readtable_cactus_wrapper(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  nuc_eos_C_ReadTable(nuceos_table_name, cctkGH);

}
