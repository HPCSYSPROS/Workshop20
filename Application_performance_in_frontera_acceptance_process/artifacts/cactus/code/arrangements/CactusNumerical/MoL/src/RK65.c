 /*@@
   @file      RK65.c
   @date      Sun May 26 03:47:15 2002
   @author    Peter Diener (based on RK45.c by Ian Hawke)
   @desc 
   RK65 following P. J. Prince and J. R. Dormand
   Journal of Computational and Applied Mathematics, volume 7, no 1, 1981
   @enddesc 
   @version   $Header$
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "ExternalVariables.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_RK65_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_RK65Add(CCTK_ARGUMENTS);

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
   @routine    MoL_RK65Add
   @date       Sun May 26 03:50:44 2002
   @author     Peter Diener (based on MoL_RK45Add by Ian Hawke)
   @desc 
   Performs a single step of a Runge-Kutta 65 type time
   integration, storing the error estimate.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_RK65Add(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
    
  CCTK_INT arraydim;

  CCTK_INT index, var, scratchstep;
  CCTK_INT totalsize;

  CCTK_REAL       * restrict UpdateVar;
  CCTK_REAL const * restrict RHSVar;
  CCTK_REAL       * restrict ScratchVar;
  CCTK_REAL       * restrict ErrorVar;
  CCTK_REAL const * restrict OldVar;

  CCTK_REAL beta, gamma, gamma_error;

  static const CCTK_REAL beta_array[7][7] = {
    { 1.0/10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { -2.0/81.0, 20.0/81.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    { 615.0/1372.0, -270.0/343.0, 1053.0/1372.0, 0.0, 0.0, 0.0, 0.0 },
    { 3243.0/5500.0, -54.0/55.0, 50949.0/71500.0, 4998.0/17875.0, 0.0, 0.0, 0.0 },
    { -26492.0/37125.0, 72.0/55.0, 2808.0/23375.0, -24206.0/37125.0, 338.0/459.0, 0.0, 0.0 },
    { 5561.0/2376.0, -35.0/11.0, -24117.0/31603.0, 899983.0/200772.0, -5225.0/1836.0, 3925.0/4056.0, 0.0 },
    { 465467.0/266112.0, -2945.0/1232.0, -5610201.0/14158144.0, 10513573.0/3212352.0, -424325.0/205632.0, 376225.0/454272.0, 0.0 }
  };

  static const CCTK_REAL gamma_array[8] = 
    { 61.0/864.0, 
      0.0, 
      98415.0/321776.0, 
      16807.0/146016.0, 
      1375.0/7344.0, 
      1375.0/5408.0,
      -37.0/1120.0,
      1.0/10.0
    };

  static const CCTK_REAL gammastar_array[8] = 
    { 821.0/10800.0,
      0.0, 
      19683.0/71825.0, 
      175273.0/912600.0, 
      395.0/3672.0, 
      785.0/2704.0,
      3.0/50.0,
      0.0 
    };

  totalsize = 1;
  for (arraydim = 0; arraydim < cctk_dim; arraydim++)
  {
    totalsize *= cctk_ash[arraydim];
  }  

  /* Real GFs */

  /* First store (dt times) the rhs in the scratch array. */

  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    const CCTK_REAL tmp = (*Original_Delta_Time) / cctkGH->cctk_timefac;

    UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedVariableIndex[var]);
    RHSVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, 0, 
						 RHSVariableIndex[var]);
    ScratchVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, var, 
                                              CCTK_FirstVarIndex("MOL::SCRATCHSPACE")
                                              + (MoL_Intermediate_Steps - 
                                               (*MoL_Intermediate_Step)));
#pragma omp parallel for
    for (index = 0; index < totalsize; index++)
    {
      ScratchVar[index] = tmp * RHSVar[index];
    }
  }


  for (var = 0; var < MoLNumEvolvedVariables; var++)
  {
    OldVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, 1, 
						 EvolvedVariableIndex[var]);
    UpdateVar = (CCTK_REAL *)CCTK_VarDataPtrI(cctkGH, 0, 
                                              EvolvedVariableIndex[var]);
    RHSVar = (CCTK_REAL const *)CCTK_VarDataPtrI(cctkGH, var, 
						 CCTK_FirstVarIndex("MOL::SCRATCHSPACE")
						 + (MoL_Intermediate_Steps - 
						  (*MoL_Intermediate_Step)));
    ErrorVar =  (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, var, 
                                             CCTK_FirstVarIndex("MOL::ERRORESTIMATE"));

    if (*MoL_Intermediate_Step - 1)
    { 

#pragma omp parallel for
      for (index = 0; index < totalsize; index++)
      {
        UpdateVar[index] = OldVar[index];
      }

      for (scratchstep = 0; 
           scratchstep < MoL_Intermediate_Steps - (*MoL_Intermediate_Step) + 1;
           scratchstep++)
      {
        
        ScratchVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, var, 
                                                  CCTK_FirstVarIndex("MOL::SCRATCHSPACE")
                                                  + scratchstep);
        
        beta = beta_array[MoL_Intermediate_Steps - (*MoL_Intermediate_Step)][scratchstep];
        
        if ( (beta > MoL_Tiny)||(beta < -MoL_Tiny) )
        {
#pragma omp parallel for
          for (index = 0; index < totalsize; index++)
          {
            UpdateVar[index] += beta * ScratchVar[index];
          }
        }
        
      }
      
    }
    else
    {

#pragma omp parallel for
      for (index = 0; index < totalsize; index++)
      {
        UpdateVar[index] = OldVar[index];
        ErrorVar[index] = 0;
      }
    
      for (scratchstep = 0; scratchstep < 8; scratchstep++)
      {
        
        ScratchVar = (CCTK_REAL*)CCTK_VarDataPtrI(cctkGH, var, 
                                                  CCTK_FirstVarIndex("MOL::SCRATCHSPACE")
                                                  + scratchstep);
        
        gamma = gamma_array[scratchstep];
        gamma_error = gamma - gammastar_array[scratchstep];
        
        if ( (gamma > MoL_Tiny)||(gamma < -MoL_Tiny) )
        {
#pragma omp parallel for
          for (index = 0; index < totalsize; index++)
          {
            UpdateVar[index] += gamma * ScratchVar[index];
            ErrorVar[index] += gamma_error * ScratchVar[index];
          }
        }
      }
      
    }
    
  }

  /* Real arrays */

  for (var = 0; var < MoLNumEvolvedArrayVariables; var++)
  {

    CCTK_ERROR("Peter has been too lazy to write the RK65 routine "
               "out for array variables. Better send him an email...");

  }
        
  return;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

