/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void WeylScal4_RegisterSymmetries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* array holding symmetry definitions */
  int sym[3];
  
  /* Register symmetries of grid functions */
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi0r");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi0i");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi1r");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi1i");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi2r");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi2i");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi3r");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi3i");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi4r");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::Psi4i");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::curvIr");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::curvIi");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::curvJr");
  
  sym[0] = -1;
  sym[1] = -1;
  sym[2] = -1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::curvJi");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::curvJ1");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::curvJ2");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::curvJ3");
  
  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;
  SetCartSymVN(cctkGH, sym, "WeylScal4::curvJ4");
  
}
