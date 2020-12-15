 /*@@
   @file      RK45.c
   @date      Sun May 26 03:47:15 2002
   @author    Ian Hawke
   @desc 
   RK45 following Forsythe, Malcolm and Moler
   (Computer Methods for Mathematical Computations).
   @enddesc 
   @version   $Header$
 @@*/

#include "ExternalVariables.h"
#include "Operators.h"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <assert.h>

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusNumerical_MoL_RK45_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RK45Add(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    MoL_RK45Add
   @date       Sun May 26 03:50:44 2002
   @author     Ian Hawke
   @desc 
   Performs a single step of a Runge-Kutta 45 type time
   integration, storing the error estimate.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_RK45Add(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
    
  /* Fehlberg coefficients
     <https://en.wikipedia.org/wiki/Runge–Kutta–Fehlberg_method> */
  static const CCTK_REAL beta_array_F[5][5] = {
    { 1.0/4.0, 0.0, 0.0, 0.0, 0.0 },
    { 3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0 },
    { 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0 },
    { 439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0 },
    { -8.0/27.0, 2.0, -3544.0/2565.0, 1859.0/4104.0, -11.0/40.0 }
  };

  static const CCTK_REAL gamma_array_F[6] = 
    { 16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0 };

  static const CCTK_REAL gammastar_array_F[6] = 
    { 25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, -1.0/5.0, 0.0 };
  
  /* Cash-Karp coefficients */
  static const CCTK_REAL beta_array_CK[5][5] = {
    { 1.0/5.0,        0.0,         0.0,           0.0,              0.0,      },
    { 3.0/40.0,       9.0/40.0,    0.0,           0.0,              0.0,      },
    { 3.0/10.0,       -9.0/10.0,   6.0/5.0,       0.0,              0.0,      },
    { -11.0/54.0,     5.0/2.0,     -70.0/27.0,    35.0/27.0,        0.0,      },
    { 1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0,
      253.0/4096.0, },
  };

  static const CCTK_REAL gamma_array_CK[6] = {
    37.0/378.0,
    0.0,
    250.0/621.0,
    125.0/594.0,
    0.0,
    512.0/1771.0,
  };

  static const CCTK_REAL gammastar_array_CK[6] = {
    2825.0/27648.0,
    0.0,
    18575.0/48384.0,
    13525.0/55296.0,
    277.0/14336.0,
    1.0/4.0,
  };

  const CCTK_REAL (* restrict beta_array)[5];
  const CCTK_REAL * restrict gamma_array;
  const CCTK_REAL * restrict gammastar_array;

  if (CCTK_EQUALS(ODE_Method, "RK45"))
  {
    beta_array = beta_array_F;
    gamma_array = gamma_array_F;
    gammastar_array = gammastar_array_F;
  }
  else if (CCTK_EQUALS(ODE_Method, "RK45CK"))
  {
    beta_array = beta_array_CK;
    gamma_array = gamma_array_CK;
    gammastar_array = gammastar_array_CK;
  }
  else
  {
    CCTK_ERROR ("internal error");
    /* Avoid compiler warnings */
    beta_array = NULL; 
    gammastar_array = NULL;
    gamma_array = NULL;
  }

  CCTK_INT rl = 0;
  if (CCTK_IsFunctionAliased("GetRefinementLevel")) {
    rl = GetRefinementLevel(cctkGH);
  }

  const int scratchvarindex = CCTK_FirstVarIndex("MOL::SCRATCHSPACE");
  if (scratchvarindex < 0)
  {
    CCTK_ERROR("Internal error");
  }
  const int errorvarindex = CCTK_FirstVarIndex("MOL::ERRORESTIMATE");
  if (errorvarindex < 0)
  {
    CCTK_ERROR("Internal error");
  }

  const int mol_step = MoL_Intermediate_Steps - *MoL_Intermediate_Step;

  /* Real GFs */

  /* First store (dt times) the RHS in the scratch array. */
  for (int var = 0; var < MoLNumEvolvedVariables; var++)
  {
    const CCTK_INT nsrcs = 1;
    const CCTK_INT srcs[] = {RHSVariableIndex[var]};
    const CCTK_INT tls[] = {0};
    const CCTK_REAL facts[] = {*Original_Delta_Time / cctkGH->cctk_timefac};
    const CCTK_INT dst = scratchvarindex + mol_step;
    const CCTK_INT tl = var;
    MoL_LinearCombination(cctkGH, dst, rl, tl, 0.0, srcs, tls, facts, nsrcs);
  }

  if (mol_step < MoL_Intermediate_Steps - 1)
  { 
    for (int var = 0; var < MoLNumEvolvedVariables; var++)
    {
      const int num_scratchsteps = mol_step + 1;
      CCTK_INT srcs[num_scratchsteps + 1];
      CCTK_INT tls[num_scratchsteps + 1];
      CCTK_REAL facts[num_scratchsteps + 1];
      CCTK_INT nsrcs = 0;
      srcs[nsrcs] = EvolvedVariableIndex[var];
      tls[nsrcs] = 1;
      facts[nsrcs] = 1.0;
      ++ nsrcs;
      for (int scratchstep = 0; scratchstep < num_scratchsteps; ++ scratchstep)
      {
        const CCTK_REAL beta = beta_array[mol_step][scratchstep];
        if (beta != 0.0)
        {
          srcs[nsrcs] = scratchvarindex + scratchstep;
          tls[nsrcs] = var;
          facts[nsrcs] = beta;
          ++ nsrcs;
        }
      }
      const CCTK_INT dst = EvolvedVariableIndex[var];
      const CCTK_INT tl = 0;
      MoL_LinearCombination(cctkGH, dst, rl, tl, 0.0, srcs, tls, facts, nsrcs);
    }
  }
  else
  {
    assert(mol_step == 5);
    for (int var = 0; var < MoLNumEvolvedVariables; var++)
    {
      const int num_scratchsteps = 6;
      CCTK_INT srcs[num_scratchsteps + 1];
      CCTK_INT tls[num_scratchsteps + 1];
      CCTK_REAL facts[num_scratchsteps + 1];
      CCTK_INT nsrcs = 0;
      srcs[nsrcs] = EvolvedVariableIndex[var];
      tls[nsrcs] = 1;
      facts[nsrcs] = 1.0;
      ++ nsrcs;
      for (int scratchstep = 0; scratchstep < num_scratchsteps; ++ scratchstep)
      {
        const CCTK_REAL gamma = gamma_array[scratchstep];
        if (gamma != 0.0)
        {
          srcs[nsrcs] = scratchvarindex + scratchstep;
          tls[nsrcs] = var;
          facts[nsrcs] = gamma;
          ++ nsrcs;
        }
      }
      const CCTK_INT dst = EvolvedVariableIndex[var];
      const CCTK_INT tl = 0;
      MoL_LinearCombination(cctkGH, dst, rl, tl, 0.0, srcs, tls, facts, nsrcs);
    }
    for (int var = 0; var < MoLNumEvolvedVariables; var++)
    {
      const int num_scratchsteps = 6;
      CCTK_INT srcs[num_scratchsteps];
      CCTK_INT tls[num_scratchsteps];
      CCTK_REAL facts[num_scratchsteps];
      CCTK_INT nsrcs = 0;
      for (int scratchstep = 0; scratchstep < num_scratchsteps; ++ scratchstep)
      {
        const CCTK_REAL gamma_error =
          gamma_array[scratchstep] - gammastar_array[scratchstep];
        if (gamma_error != 0.0)
        {
          srcs[nsrcs] = scratchvarindex + scratchstep;
          tls[nsrcs] = var;
          facts[nsrcs] = gamma_error;
          ++ nsrcs;
        }
      }
      const CCTK_INT dst = errorvarindex;
      const CCTK_INT tl = var;
      MoL_LinearCombination(cctkGH, dst, rl, tl, 0.0, srcs, tls, facts, nsrcs);
    }
  }

  /* Real arrays */

  for (int var = 0; var < MoLNumEvolvedArrayVariables; var++)
  {
    CCTK_ERROR ("Ian has been too lazy to write the RK45 routine "
                "out for array variables. Better send him an email...");
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

