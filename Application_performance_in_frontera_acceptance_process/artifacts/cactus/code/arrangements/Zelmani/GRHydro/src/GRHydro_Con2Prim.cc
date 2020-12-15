 /*@@
   @file      GRHydro_Reconstruct.cc
   @date      April 2016
   @author    Christian D. Ott
   @desc
   Wrapper routine for conservative to primitive conversion
   @enddesc
 @@*/


#include <iostream>
#include <cassert>
#include <cstring>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "SpaceMask.h"

// prototype
void GRHydro_Con2Prim_simpleEOS_CXX(CCTK_ARGUMENTS);

static inline __attribute__((always_inline))
  void upper_metric(const double gxx, const double gxy,
		    const double gxz, const double gyy,
		    const double gyz, const double gzz,
		    const double sdetg,
		    double* uxx, double* uxy, double* uxz,
		    double* uyy, double* uyz, double* uzz){

  double invdet = 1.0 / (sdetg*sdetg);

  *uxx = (-gyz*gyz + gyy*gzz)*invdet;
  *uxy = ( gxz*gyz - gxy*gzz) *invdet;
  *uxz = (-gxz*gyy + gxy*gyz)*invdet;
  *uyy = (-gxz*gxz + gxx*gzz)*invdet;
  *uyz = ( gxy*gxz - gxx*gyz) *invdet;
  *uzz = (-gxy*gxy + gxx*gyy)*invdet;
}

extern "C" void GRHydro_Con2Prim_CXX(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(GRHydro_disable_hydro_update) {
    return;
  }

  if(*evolve_temper == 1) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
	       "Hot EOS not yet implemented for C++ C2P");
  } else {
    GRHydro_Con2Prim_simpleEOS_CXX(CCTK_PASS_CTOC);
  }

  if(evolve_tracer) { // this is a boolean parameter
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
	       "Simple EOS C2P does not handle tracers");
  }

  if(*evolve_Y_e != 0) {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
	       "Simple EOS C2P does not handle Y_e");
  }
}

// currently this helper function exists only in F90. This declaration has to
// be in this .cc file and cannot be moved into an .h file since the Cactus
// preprocessor does not parse .h files.
extern "C"
void CCTK_FNAME(GRHydro_DeclareEvolutionMask)(cGH* const & cctkGH,
                                              CCTK_REAL*& evolution_mask,
                                              CCTK_INT& evolution_mask_valid);
void GRHydro_Con2Prim_simpleEOS_CXX(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // helper vars
  CCTK_INT myproc = CCTK_MyProc(cctkGH);

  // NOTE: need to add MP support 
  
  // CarpetEvolutionMask support
  CCTK_REAL* evolution_mask = NULL;
  CCTK_INT check_evolution_mask;
  CCTK_FNAME(GRHydro_DeclareEvolutionMask)(cctkGH, evolution_mask, check_evolution_mask);
}
