 /*@@
   @file      RK4-RK2.c
   @date      2012-03-25
   @author    Christian Reisswig
   @desc 
   A routine to perform homegrown RK4RK2 evolution. Mostly copied from
   genericRK.c
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <stdio.h>

#include "ExternalVariables.h"

/* #define MOLDEBUG */

 /*@@
   @routine    MoL_RK4_RK2_Add
   @date       
   @author     
   @desc 
   Performs a single step of a RK4_RK2 type time
   integration.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 @@*/ 


void MoL_RK4_RK2_Add(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
    
  static int scratchspace_firstindex      = -1;
  static int scratchspace_firstindex_slow = -1;
  if (scratchspace_firstindex < 0) {
    scratchspace_firstindex      = CCTK_FirstVarIndex("MOL::SCRATCHSPACE");
    scratchspace_firstindex_slow = CCTK_FirstVarIndex("MOL::SCRATCHSPACESLOW");
  }
  
  int const step = MoL_Intermediate_Steps - *MoL_Intermediate_Step;
  
  int totalsize = 1;
  for (int d=0; d<cctk_dim; ++d) totalsize *= cctk_ash[d];
  
  CCTK_REAL const dt = *Original_Delta_Time / cctkGH->cctk_timefac;
  
  
  
  int const allvar1 = MoLNumEvolvedVariables + MoLNumEvolvedVariablesSlow;
#pragma omp parallel for schedule(dynamic)
  for (int var1=0; var1<allvar1; ++var1) {
    
    if (var1 < MoLNumEvolvedVariables) {
      /* a fast variable */
      int const var = var1;
      
      CCTK_REAL *restrict const UpdateVar =
        CCTK_VarDataPtrI(cctkGH, 0, EvolvedVariableIndex[var]);
      CCTK_REAL const *restrict const OldVar =
        CCTK_VarDataPtrI(cctkGH, 1, EvolvedVariableIndex[var]);
      CCTK_REAL const *restrict const RHSVar =
        CCTK_VarDataPtrI(cctkGH, 0, RHSVariableIndex[var]);
      
#define SCRATCHINDEX(step)                                              \
      (scratchspace_firstindex + (step))
      CCTK_REAL *restrict const ScratchVar =
        CCTK_VarDataPtrI(cctkGH, var, SCRATCHINDEX(0));
      
      switch (step) {
        
      case 0:
        for (int i=0; i<totalsize; ++i) {
          CCTK_REAL const scaled_rhs = dt * RHSVar[i];
          ScratchVar[i] = OldVar[i] + 1.0/6.0 * scaled_rhs;
          UpdateVar[i] = OldVar[i] + 0.5 * scaled_rhs;
        }
        break;
        
      case 1:
        for (int i=0; i<totalsize; ++i) {
          CCTK_REAL const scaled_rhs = dt * RHSVar[i];
          ScratchVar[i] += 1.0/3.0 * scaled_rhs;
          UpdateVar[i] = OldVar[i] + 0.5 * scaled_rhs;
        }
        break;
        
      case 2:
        for (int i=0; i<totalsize; ++i) {
          CCTK_REAL const scaled_rhs = dt * RHSVar[i];
          ScratchVar[i] += 1.0/3.0 * scaled_rhs;
          UpdateVar[i] = OldVar[i] + scaled_rhs;
        }
        break;
        
      case 3:
        for (int i=0; i<totalsize; ++i) {
          CCTK_REAL const scaled_rhs = dt * RHSVar[i];
          /* ScratchVar contains OldVar */
          UpdateVar[i] = ScratchVar[i] + 1.0/6.0 * scaled_rhs;
        }
        break;
        
      default:
        assert(0);
      }
#undef SCRATCHINDEX
      
    } else {
      /* a slow variable */
      int const var = var1 - MoLNumEvolvedVariables;
      
      CCTK_REAL *restrict const UpdateVar =
        CCTK_VarDataPtrI(cctkGH, 0, EvolvedVariableIndexSlow[var]);
      CCTK_REAL const *restrict const OldVar =
        CCTK_VarDataPtrI(cctkGH, 1, EvolvedVariableIndexSlow[var]);
      CCTK_REAL const *restrict const RHSVar =
        CCTK_VarDataPtrI(cctkGH, 0, RHSVariableIndexSlow[var]);
      
#define SCRATCHINDEX(step)                                              \
      (scratchspace_firstindex_slow + (step))
      CCTK_REAL *restrict const ScratchVar =
        CCTK_VarDataPtrI(cctkGH, var, SCRATCHINDEX(0));
      
      switch (step) {
        
      case 0:
        for (int i=0; i<totalsize; ++i) {
          CCTK_REAL const scaled_rhs = dt * RHSVar[i];
          CCTK_REAL const scratchval = OldVar[i] + scaled_rhs;
          ScratchVar[i] = scratchval;
          UpdateVar[i] = scratchval;
        }
        break;
        
      case 1:
      case 2:
        for (int i=0; i<totalsize; ++i) {
          /* This is the same value as for the previous MoL step.
             However, MoL_PostStep may have modified it (e.g. enforced
             a constraint), so we need to recreate the original value
             here for consistency. */
          /* ScratchVar contains OldVar */
          UpdateVar[i] = ScratchVar[i];
        }
        break;
        
      case 3:
        for (int i=0; i<totalsize; ++i) {
          CCTK_REAL const scaled_rhs = dt * RHSVar[i];
          /* ScratchVar contains OldVar */
          UpdateVar[i] =
            0.5 * OldVar[i] + 0.5 * ScratchVar[i] + 0.5 * scaled_rhs;
        }
        break;
        
      default:
        assert(0);
      }
#undef SCRATCHINDEX
      
    } /* if fast or slow */
  }   /* for var */
  
}
