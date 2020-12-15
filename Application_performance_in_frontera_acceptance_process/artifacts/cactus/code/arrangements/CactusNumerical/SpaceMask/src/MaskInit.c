 /*@@
   @file      MaskInit.c
   @date      Fri Apr 26 17:20:13 2002
   @author    Miguel Alcubierre
   @desc 
   Initialise the mask (I just copied LapseInits.c)
   @enddesc 
   @version $Header$
 @@*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "cctk.h"

#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"
#include "SpaceMask.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_Einstein_MaskInit_c);


/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void MaskOne(CCTK_ARGUMENTS);
void MaskSym(CCTK_ARGUMENTS);
void MaskSym_emask(CCTK_ARGUMENTS);
void MaskZero(CCTK_ARGUMENTS);

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

/*@@
   @routine    MaskSym
   @date       October 2002
   @author     Denis Pollney
   @desc 
   Scheduled routine to set symmetries for mask
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
@@*/
void MaskSym(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS

  int one;
  int sym[3];

  one = 1;
  
  sym[0] = one;
  sym[1] = one;
  sym[2] = one;

  SetCartSymVN(cctkGH, sym, "spacemask::space_mask");

  return;
}

/*@@
   @routine    MaskSym_emask
   @date       Fri 3 May 2002
   @author     Gabrielle Allen
   @desc 
   Scheduled routine to set symmetries for mask
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
@@*/
void MaskSym_emask(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS

  int one;
  int sym[3];

  one = 1;
  
  sym[0] = one;
  sym[1] = one;
  sym[2] = one;

  SetCartSymVN(cctkGH, sym, "spacemask::emask");

  return;
}

/*@@
   @routine    MaskOne
   @date       
   @author     Miguel Alcubierre
   @desc 
   Scheduled routine to initialise the mask to one.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
@@*/
void MaskOne(CCTK_ARGUMENTS)
{
  int i;
  DECLARE_CCTK_ARGUMENTS;

  for(i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
  {
    emask[i] = 1.0;
  }

  return;
}

/*@@
   @routine    CheckMask
   @date       
   @author     Erik Schnetter
   @desc 
   Ensure that all mask values are legal.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
@@*/
void CheckMask(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS
  int i;

  for(i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
  {
    if (   fabs(emask[i] - 1.0) > 1.0e-12
        && fabs(emask[i] - 0.5) > 1.0e-12 
        && fabs(emask[i] - 0.0) > 1.0e-12)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Illegal mask value %g detected at grid point %d",
                  (double)emask[i], i);
    }
  }
}

/*@@
   @routine    MaskZero
   @date       
   @author     Denis Pollney
   @desc 
               Initialise the mask to zero.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
 
@@*/
void MaskZero(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS

  int i;

  for(i = 0; i < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; i++)
    space_mask[i] = 0;

  return;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
