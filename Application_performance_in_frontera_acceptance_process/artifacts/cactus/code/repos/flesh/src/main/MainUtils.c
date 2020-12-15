 /*@@
   @file      MainUtils.c
   @date      Sep 22 1999
   @author    Thomas Radke, Gabrielle Allen
   @desc
              Utility Flesh routines
   @enddesc
   @version   $Id$
 @@*/

#include <string.h>

#include "cctk_Flesh.h"
#include "cctk_Parameter.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_MainUtils_c);

/********************************************************************
 ********************* External Routine Prototypes ******************
 ********************************************************************/
int CCTK_RunTitle(int len, char *title);


 /*@@
   @routine    CCTK_RunTitle
   @date       Sun Sep 17 2000
   @author     Gabrielle Allen
   @desc
               Returns the simulation description
   @enddesc

   @var        len
   @vdesc      The size (in characters) of the title buffer
   @vtype      int
   @vio        in
   @endvar
   @var        title
   @vdesc      The title buffer
   @vtype      char *
   @vio        out
   @endvar

   @returntype int
   @returndesc
               The length of the (copied) title, or<BR>
               -1 if parameter value of "Cactus::cctk_run_title" is unknown
   @endreturndesc
@@*/
int CCTK_RunTitle(int len, char *title)
{
  int retval;
  const char *cctk_title;

  retval = -1;

  cctk_title = *(const char *const *)
               CCTK_ParameterGet("cctk_run_title", "Cactus", NULL);
  if (cctk_title)
  {
    strncpy (title, *cctk_title ? cctk_title : "Cactus Simulation", len-1);
    title[len-1] = 0;
    retval = strlen(title);
  }
  return retval;
}
