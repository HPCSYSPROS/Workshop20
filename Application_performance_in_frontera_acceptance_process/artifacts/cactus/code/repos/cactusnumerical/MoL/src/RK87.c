 /*@@
   @file      RK87.c
   @date      Sun May 26 03:47:15 2002
   @author    Peter Diener (based on RK45.c by Ian Hawke)
   @desc 
   RK87 following P. J. Prince and J. R. Dormand
   Journal of Computational and Applied Mathematics, volume 7, no 1, 1981
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"
#include "Operators.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_RK87_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RK87Add(CCTK_ARGUMENTS);

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
   @routine    MoL_RK87Add
   @date       Sun May 26 03:50:44 2002
   @author     Peter Diener (based on MoL_RK45Add by Ian Hawke)
   @desc 
   Performs a single step of a Runge-Kutta 87 type time
   integration, storing the error estimate.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_RK87Add(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
    
  CCTK_INT arraydim;

  static CCTK_INT scratchspace_firstindex = -99;
  static CCTK_INT error_firstindex = -99;
  CCTK_INT var, scratchstep;
  CCTK_INT totalsize;

  static const CCTK_REAL beta_array[12][12] = {
    { 1.0/18.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 1.0/48.0, 1.0/16.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 1.0/32.0, 0.0, 3.0/32.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 5.0/16.0, 0.0, -75.0/64.0, 75.0/64.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 3.0/80.0, 0.0, 0.0, 3.0/16.0, 3.0/20.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 29443841.0/614563906.0, 0.0, 0.0, 77736538.0/692538347.0, -28693883.0/1125000000.0, 23124283.0/1800000000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 16016141.0/946692911.0, 0.0, 0.0, 61564180.0/158732637.0, 22789713.0/633445777.0, 545815736.0/2771057229.0, -180193667.0/1043307555.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 39632708.0/573591083.0, 0.0, 0.0, -433636366.0/683701615.0, -421739975.0/2616292301.0, 100302831.0/723423059.0, 790204164.0/839813087.0, 800635310.0/3783071287.0, 0.0, 0.0, 0.0, 0.0 },
    { 246121993.0/1340847787.0, 0.0, 0.0, -37695042795.0/15268766246.0, -309121744.0/1061227803.0, -12992083.0/490766935.0, 6005943493.0/2108947869.0, 393006217.0/1396673457.0, 123872331.0/1001029789.0, 0.0, 0.0, 0.0 },
    { -1028468189.0/846180014.0, 0.0, 0.0, 8478235783.0/508512852.0, 1311729495.0/1432422823.0, -10304129995.0/1701304382.0, -48777925059.0/3047939560.0, 15336726248.0/1032824649.0, -45442868181.0/3398467696.0, 3065993473.0/597172653.0, 0.0, 0.0 },
    { 185892177.0/718116043.0, 0.0, 0.0, -3185094517.0/667107341.0, -477755414.0/1098053517.0, -703635378.0/230739211.0, 5731566787.0/1027545527.0, 5232866602.0/850066563.0, -4093664535.0/808688257.0, 3962137247.0/1805957418.0, 65686358.0/487910083.0, 0.0 },
    { 403863854.0/491063109.0, 0.0, 0.0, -5068492393.0/434740067.0, -411421997.0/543043805.0, 652783627.0/914296604.0, 11173962825.0/925320556.0, -13158990841.0/6184727034.0, 3936647629.0/1978049680.0, -160528059.0/685178525.0, 248638103.0/1413531060.0, 0.0 }
  };

  static const CCTK_REAL gamma_array[13] = 
    { 14005451.0/335480064.0,
      0.0,
      0.0,
      0.0,
      0.0,
      -59238493.0/1068277825.0,
      181606767.0/758867731.0,
      561292985.0/797845732.0,
      -1041891430.0/1371343529.0,
      760417239.0/1151165299.0,
      118820643.0/751138087.0,
      -528747749.0/2220607170.0,
      1.0/4.0
    };

  static const CCTK_REAL gammastar_array[13] = 
    { 13451932.0/455176623.0,
      0.0,
      0.0,
      0.0,
      0.0,
      -808719846.0/976000145.0,
      1757004468.0/5645159321.0,
      656045339.0/265891186.0,
      -3867574721.0/1518517206.0,
      465885868.0/322736535.0,
      53011238.0/667516719.0,
      2.0/45.0,
      0.0
    };

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_ash[arraydim];
  }  

  if (scratchspace_firstindex == -99)
  {
    scratchspace_firstindex = CCTK_FirstVarIndex("MOL::SCRATCHSPACE");
  }
  if (error_firstindex == -99)
  {
    error_firstindex = CCTK_FirstVarIndex("MOL::ERRORESTIMATE");
  }

  /* Real GFs */

  /* First store (dt times) the rhs in the scratch array. */

  CCTK_INT rl = 0;
  if (CCTK_IsFunctionAliased("GetRefinementLevel")) {
    rl = GetRefinementLevel(cctkGH);
  }
  CCTK_INT tl = 0;
  
  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    int const step = MoL_Intermediate_Steps - (*MoL_Intermediate_Step);
    int const scratchvarindex = scratchspace_firstindex + step;

    int const nsrcs = 1;
    CCTK_INT const srcs[] = {RHSVariableIndex[var]};
    CCTK_INT const tls[] = {0};
    CCTK_REAL const facts[] = {(*Original_Delta_Time) / cctkGH->cctk_timefac};
    MoL_LinearCombination(cctkGH,
                          scratchvarindex, rl, var, 0.0,
                          srcs, tls, facts, nsrcs);
  }

  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    int const step = MoL_Intermediate_Steps - (*MoL_Intermediate_Step);

    if (*MoL_Intermediate_Step - 1)
    {
      int const nsrcs = step + 2;
      CCTK_INT srcs[nsrcs];
      CCTK_INT tls[nsrcs];
      CCTK_REAL facts[nsrcs];
      srcs[0]  = EvolvedVariableIndex[var];
      tls[0]   = 1;
      facts[0] = 1.0;

      for (scratchstep = 0; scratchstep < step + 1; scratchstep++)
      {
        int const scratchvarindex = scratchspace_firstindex + scratchstep;
        srcs [scratchstep+1] = scratchvarindex;
        tls  [scratchstep+1] = var;
        facts[scratchstep+1] = beta_array[step][scratchstep];
      }

      MoL_LinearCombination(cctkGH,
                            EvolvedVariableIndex[var], rl, tl, 0.0,
                            srcs, tls, facts, nsrcs);
    }
    else
    {
      {
        int const nsrcs = 14;
        CCTK_INT srcs[nsrcs];
        CCTK_INT tls[nsrcs];
        CCTK_REAL facts[nsrcs];
        srcs[0]  = EvolvedVariableIndex[var];
        tls[0]   = 1;
        facts[0] = 1.0;
        
        for (scratchstep = 0; scratchstep < 13; scratchstep++)
        {
          int const scratchvarindex = scratchspace_firstindex + scratchstep;
          srcs [scratchstep+1] = scratchvarindex;
          tls  [scratchstep+1] = var;
          facts[scratchstep+1] = gamma_array[scratchstep];
        }

        MoL_LinearCombination(cctkGH,
                              EvolvedVariableIndex[var], rl, tl, 0.0,
                              srcs, tls, facts, nsrcs);
      }

      {
        int const nsrcs = 13;
        CCTK_INT srcs[nsrcs];
        CCTK_INT tls[nsrcs];
        CCTK_REAL facts[nsrcs];
        
        for (scratchstep = 0; scratchstep < 13; scratchstep++)
        {
          int const scratchvarindex = scratchspace_firstindex + scratchstep;
          srcs [scratchstep] = scratchvarindex;
          tls  [scratchstep] = var;
          facts[scratchstep] =
            gamma_array[scratchstep] - gammastar_array[scratchstep];
        }

        MoL_LinearCombination(cctkGH,
                              error_firstindex, rl, var, 0.0,
                              srcs, tls, facts, nsrcs);
      }

    }

  }

  /* Arrays */

  if (MoLNumEvolvedArrayVariables > 0)
  {
   CCTK_ERROR("Peter has been too lazy to write the RK87 routine "
              "out for array variables. Better send him an email...");
  }
        
  return;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

