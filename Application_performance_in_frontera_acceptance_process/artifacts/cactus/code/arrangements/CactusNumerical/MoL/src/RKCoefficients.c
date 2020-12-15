 /*@@
   @file      RKCoefficients.c
   @date      Tue May 21 02:46:51 2002
   @author    Ian Hawke
   @desc 
   The routine setting up the coefficients for the generic Runge-Kutta
   style integrator. At some point this should be extended so that
   these can be set from the parameter file.
   @version   $Header$
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_ErrorCodes.h"
#include "util_Table.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusBase_MoL_RKCoefficients_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MoL_SetupRKCoefficients(CCTK_ARGUMENTS);

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
   @routine    MoL_SetupRKCoefficients
   @date       Tue May 21 02:49:06 2002
   @author     Ian Hawke
   @desc 
   Sets up the coefficients of the RKAlpha and Beta arrays. These
   are currently set "by hand" for the Runge-Kutta and generic ICN
   methods. Should add the ability to set from a parameter file.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void MoL_SetupRKCoefficients(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT i, j;

  CCTK_INT ierr, options_table;

  if (CCTK_EQUALS(Generic_Type,"Classic RK3"))
  {
    if (MoL_Num_Scratch_Levels != 2)
    {
      CCTK_ERROR("For Classic RK3, MoL_Num_Scratch_Levels "
                 "should be at least 2");
    }
    if (MoL_Intermediate_Steps != 3)
    {
      CCTK_ERROR("For Classic RK3, MoL_Intermediate_Steps "
                 "should be at least 3");
    }
    for (i = 0; i < MoL_Intermediate_Steps; i++) 
    {
      for (j = 0; j < MoL_Num_Scratch_Levels + 1; j++) 
      {  
        RKAlphaCoefficients[i * MoL_Intermediate_Steps + j] = 0.0;
      }
      RKBetaCoefficients[i] = 0.0;
    }
    RKAlphaCoefficients[0] = 1.0;
    RKAlphaCoefficients[3] = 1.0;
    RKAlphaCoefficients[6] = 1.0 / 9.0;
    RKAlphaCoefficients[7] = 4.0 / 9.0;
    RKAlphaCoefficients[8] = 4.0 / 9.0;
    RKBetaCoefficients[0] = 0.5;
    RKBetaCoefficients[1] = 0.75;
    RKBetaCoefficients[2] = 4.0 / 9.0;    
  }
  else if (CCTK_EQUALS(Generic_Type,"ICN")) 
  {
    for (i = 0; i < MoL_Intermediate_Steps; i++) 
    {
      RKAlphaCoefficients[i * MoL_Intermediate_Steps] = 1.0;
      for (j = 1; j < MoL_Num_Scratch_Levels + 1; j++) 
      {
        RKAlphaCoefficients[i * MoL_Intermediate_Steps + j] = 0.0;
      }
      if (i == MoL_Intermediate_Steps-1)
      {
        RKBetaCoefficients[i] = 1.0;
      }
      else
      {
        RKBetaCoefficients[i] = 0.5;
      } 
    }
  }
  else if (CCTK_EQUALS(Generic_Type,"RK")) 
  {
    if (MoL_Num_Scratch_Levels < MoL_Intermediate_Steps - 1)
    {
      CCTK_ERROR("For generic RK methods, MoL_Num_Scratch_Levels "
                 "should be at least MoL_Intermediate_Steps - 1");
    }
    for (i = 0; i < MoL_Intermediate_Steps; i++) 
    {
      for (j = 0; j < MoL_Num_Scratch_Levels + 1; j++) 
      {  
        RKAlphaCoefficients[i * MoL_Intermediate_Steps + j] = 0.0;
      }
      RKBetaCoefficients[i] = 0.0;
    }
    if (MoL_Intermediate_Steps == 1)
    {
      RKAlphaCoefficients[0] = 1.0;
      RKBetaCoefficients[0] = 1.0;
    }
    else if (MoL_Intermediate_Steps == 2)
    {
      RKAlphaCoefficients[0] = 1.0;
      RKAlphaCoefficients[2] = 0.5;
      RKAlphaCoefficients[3] = 0.5;
      RKBetaCoefficients[0] = 1.0;
      RKBetaCoefficients[1] = 0.5;

    }
    else if (MoL_Intermediate_Steps == 3)
    {
      RKAlphaCoefficients[0] = 1.0;
      RKAlphaCoefficients[3] = 0.75;
      RKAlphaCoefficients[4] = 0.25;
      RKAlphaCoefficients[6] = 1.0 / 3.0;
      RKAlphaCoefficients[8] = 2.0 / 3.0;
      RKBetaCoefficients[0] = 1.0;
      RKBetaCoefficients[1] = 0.25;
      RKBetaCoefficients[2] = 2.0 / 3.0;
    }
    else if (MoL_Intermediate_Steps == 4)
    {
      RKAlphaCoefficients[0] = 1.0;
      RKAlphaCoefficients[4] = 1.0;
      RKAlphaCoefficients[8] = 1.0;
      RKAlphaCoefficients[12] = -1.0 / 3.0;
      RKAlphaCoefficients[13] = 1.0 / 3.0;
      RKAlphaCoefficients[14] = 2.0 / 3.0;
      RKAlphaCoefficients[15] = 1.0 / 3.0;
      RKBetaCoefficients[0] = 0.5;
      RKBetaCoefficients[1] = 0.5;
      RKBetaCoefficients[2] = 1.0;
      RKBetaCoefficients[3] = 1.0 / 6.0;
    }
    else 
    {
    CCTK_ERROR("RKCoefficients cannot do generic RK methods "
               "with MoL_Intermediate_Steps greater than 4");
    }
  }
  else if (CCTK_EQUALS(Generic_Type,"Table")) 
  {
    if (MoL_Num_Scratch_Levels < MoL_Intermediate_Steps - 1)
    {
      CCTK_ERROR("For generic methods, MoL_Num_Scratch_Levels "
                 "should be at least MoL_Intermediate_Steps - 1");
    }
    options_table =
      Util_TableCreateFromString(Generic_Method_Descriptor);
    if (options_table < 0)
    {
      CCTK_ERROR("Failed to create table from "
                 "Generic_Method_Descriptor!");
    }
    ierr = Util_TableGetRealArray(options_table,
                                  (MoL_Num_Scratch_Levels + 1) * 
                                  MoL_Intermediate_Steps,
                                  RKAlphaCoefficients,
                                  "GenericAlphaCoeffs");
    if (ierr < (MoL_Num_Scratch_Levels + 1) * MoL_Intermediate_Steps )
    {
      if (ierr >= 0)
      {
        CCTK_ERROR("Insufficient elements in the specified "
                   "GenericAlphaCoeffs array");
      }
      else if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY)
      {
        CCTK_ERROR("When using the generic table options you "
                   "must set \"GenericAlphaCoeffs\" in the options table");
      }
      else
      {
        CCTK_ERROR("Table error - check with Ian.");
      }
    }
    ierr = Util_TableGetRealArray(options_table,
                                  MoL_Intermediate_Steps,
                                  RKBetaCoefficients,
                                  "GenericBetaCoeffs");
    if (ierr < MoL_Intermediate_Steps)
    {
      if (ierr >= 0)
      {
        CCTK_ERROR("Insufficient elements in the specified "
                  "GenericBetaCoeffs array");
      }
      else if (ierr == UTIL_ERROR_TABLE_NO_SUCH_KEY)
      {
        CCTK_ERROR("When using the generic table options you "
                  "must set \"GenericBetaCoeffs\" in the options table");
      }
      else
      {
        CCTK_ERROR("Table error - check with Ian.");
      }
    }
    ierr = Util_TableDestroy(options_table);
  }
  else
  {
    CCTK_ERROR("RKCoefficients does not recognize the value "
               "of Generic_Type");
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
