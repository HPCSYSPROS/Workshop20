/*@@
  @file      Evolve.c
  @date      Thu Jul 29 2004
  @author    Erik Schnetter
  @desc
  Evolve the static conformal factor.
  @enddesc
  @version $Header$
@@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <string.h>

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_StaticConformal_Evolve_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void StaticConformal_Evolve(CCTK_ARGUMENTS);

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
  @routine    StaticConformal_Evolve
  @date       Thu Jul 29 2004
  @author     Erik Schnetter
  @desc
  Evolve the static conformal factor.
  @enddesc
  @calls
  @calledby

@@*/
void StaticConformal_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  size_t const npoints = cctk_lsh[0] * cctk_lsh[1] * cctk_lsh[2];

  if (*conformal_state >= 1) {
    memcpy(psi, psi_p, npoints * sizeof *psi);
  }

  if (*conformal_state >= 2) {
    memcpy(psix, psix_p, npoints * sizeof *psix);
    memcpy(psiy, psiy_p, npoints * sizeof *psiy);
    memcpy(psiz, psiz_p, npoints * sizeof *psiz);
  }

  if (*conformal_state >= 3) {
    memcpy(psixx, psixx_p, npoints * sizeof *psixx);
    memcpy(psixy, psixy_p, npoints * sizeof *psixy);
    memcpy(psixz, psixz_p, npoints * sizeof *psixz);
    memcpy(psiyy, psiyy_p, npoints * sizeof *psiyy);
    memcpy(psiyz, psiyz_p, npoints * sizeof *psiyz);
    memcpy(psizz, psizz_p, npoints * sizeof *psizz);
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
