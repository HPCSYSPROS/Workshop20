/*@@
  @file      Initialise.c
  @date      Thu Apr 25 19:53:02 2002
  @author    Tom Goodale
  @desc
  Code to initialise stuff for StaticConfomal
  @enddesc
  @version $Header$
@@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_StaticConformal_Initialise_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void StaticConformal_InitialiseState(CCTK_ARGUMENTS);

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
  @routine    StaticConformal_InitialiseState
  @date       Thu Apr 25 19:11:32 2002
  @author     Tom Goodale
  @desc
  Scheduled routine to set the conformal_state variable to 0
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

@@*/
void StaticConformal_InitialiseState(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *conformal_state = 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
