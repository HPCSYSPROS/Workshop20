 /*@@
   @file      FortranBindingsIO.c
   @date      Thu Feb  18 14:27:18 1999
   @author    Gabrielle Allen
   @desc 
              Fortran bindings for the IO functions
   @enddesc 
   @version   $Id$
 @@*/ 

#include <stdlib.h>
#include "cctk_Flesh.h"
#include "cctk_FortranString.h"
#include "cctk_IO.h"
#include "cctk_IOMethods.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(IO_FortranBindingsIO_c);

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

void CCTK_FCALL CCTK_FNAME (CCTK_OutputGH)
                           (int *istat, const cGH **GH);
void CCTK_FCALL CCTK_FNAME (CCTKi_RegisterIOMethod)
                           (int *handle, TWO_FORTSTRING_ARG);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTK_OutputGH
   @date       
   @author     
   @desc
   Fortran binding for CCTK_OutputGH
   @enddesc
   @var        istat
   @vdesc      error code
   @vtype      int
   @vio        out
   @endvar
   @var        GH
   @vdesc      grid hierachy
   @vtype      const GH *
   @vio        in
   @endvar

   @returntype int
   @returndesc

   @endreturndesc
@@*/
void CCTK_FCALL CCTK_FNAME (CCTK_OutputGH)
                           (int *istat, const cGH **GH)
{
  *istat = CCTK_OutputGH (*GH);
}

 /*@@
   @routine    CCTKi_RegisterIOMethod
   @date       
   @author     
   @desc

   @enddesc
@@*/

/* FIXME: This function is no longer needed, test can remove */

void CCTK_FCALL CCTK_FNAME (CCTKi_RegisterIOMethod)
                           (int *handle, TWO_FORTSTRING_ARG)
{
  TWO_FORTSTRING_CREATE (thorn, name);
  *handle = CCTKi_RegisterIOMethod (thorn, name);
  free (thorn);
  free (name);
}

