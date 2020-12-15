/*@@
  @file      ParamCheck.c
  @date      Thu Apr 25 19:02:51 2002
  @author    Tom Goodale
  @desc
  Parameter checking stuff for ADMBase
  @enddesc
  @version $Header$
@@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_ADMBase_ParamCheck_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

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
  @routine    ADMBase_ParamCheck
  @date       Thu Apr 25 19:04:06 2002
  @author     Tom Goodale
  @desc
  Scheduled routine to detect invalid parameter settings.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

@@*/
void ADMBase_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_data, "Cartesian Minkowski") &&
      !CCTK_EQUALS(metric_type, "physical")) {
    CCTK_PARAMWARN("ADMBase asked to setup Cartesian Minkowski initial data, "
                   "but metric_type is not \"physical\".  Perhaps you want "
                   "\"Conformal Minkowski\" data, provided by IDSimple?");
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
