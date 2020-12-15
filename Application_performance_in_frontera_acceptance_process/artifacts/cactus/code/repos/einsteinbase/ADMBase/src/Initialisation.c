/*@@
  @header    Initialisation.c
  @date      Thu Apr 25 22:58:03 2002
  @author    Tom Goodale
  @desc
  Do all the lapse, shift, metric and curvature initialisation known about by
  ADMBase.
  @enddesc
  @version $Header$
@@*/

#include "cctk.h"

#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_ADMBase_Initialisation_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

void ADMBase_CartesianMinkowski(CCTK_ARGUMENTS);
void ADMBase_LapseOne(CCTK_ARGUMENTS);
void ADMBase_ShiftZero(CCTK_ARGUMENTS);
void ADMBase_DtLapseZero(CCTK_ARGUMENTS);
void ADMBase_DtShiftZero(CCTK_ARGUMENTS);

void ADMBase_SetShiftStateOn(CCTK_ARGUMENTS);
void ADMBase_SetShiftStateOff(CCTK_ARGUMENTS);
void ADMBase_SetDtLapseStateOn(CCTK_ARGUMENTS);
void ADMBase_SetDtLapseStateOff(CCTK_ARGUMENTS);
void ADMBase_SetDtShiftStateOn(CCTK_ARGUMENTS);
void ADMBase_SetDtShiftStateOff(CCTK_ARGUMENTS);

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
  @routine    ADMBase_CartesianMinkowski
  @date       Thu Apr 25 23:12:18 2002
  @author     Tom Goodale
  @desc
  Scheduled routine to initialise the metric and extrinsic curvature to
  Minkowski space in cartesian coordinate values.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_CartesianMinkowski(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
    gxx[i] = 1.0;
    gyy[i] = 1.0;
    gzz[i] = 1.0;

    gxy[i] = 0.0;
    gxz[i] = 0.0;
    gyz[i] = 0.0;

    kxx[i] = 0.0;
    kyy[i] = 0.0;
    kzz[i] = 0.0;

    kxy[i] = 0.0;
    kxz[i] = 0.0;
    kyz[i] = 0.0;
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric") > 1) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      gxx_p[i] = 1.0;
      gyy_p[i] = 1.0;
      gzz_p[i] = 1.0;

      gxy_p[i] = 0.0;
      gxz_p[i] = 0.0;
      gyz_p[i] = 0.0;
    }
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::metric") > 2) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      gxx_p_p[i] = 1.0;
      gyy_p_p[i] = 1.0;
      gzz_p_p[i] = 1.0;

      gxy_p_p[i] = 0.0;
      gxz_p_p[i] = 0.0;
      gyz_p_p[i] = 0.0;
    }
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::curv") > 1) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      kxx_p[i] = 0.0;
      kyy_p[i] = 0.0;
      kzz_p[i] = 0.0;

      kxy_p[i] = 0.0;
      kxz_p[i] = 0.0;
      kyz_p[i] = 0.0;
    }
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::curv") > 2) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {

      kxx_p_p[i] = 0.0;
      kyy_p_p[i] = 0.0;
      kzz_p_p[i] = 0.0;

      kxy_p_p[i] = 0.0;
      kxz_p_p[i] = 0.0;
      kyz_p_p[i] = 0.0;
    }
  }
}

/*@@
  @routine    ADMBase_LapseOne
  @date       Thu Apr 25 23:12:18 2002
  @author     Tom Goodale
  @desc
  Scheduled routine to initialise the lapse to one.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_LapseOne(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
    alp[i] = 1.0;
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse") > 1) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      alp_p[i] = 1.0;
    }
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::lapse") > 2) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      alp_p_p[i] = 1.0;
    }
  }
}

/*@@
  @routine    ADMBase_ShiftZero
  @date       Thu Apr 25 23:12:18 2002
  @author     Tom Goodale
  @desc
  Scheduled routine to initialise the shift to zero.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_ShiftZero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
    betax[i] = 0.0;
    betay[i] = 0.0;
    betaz[i] = 0.0;
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift") > 1) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      betax_p[i] = 0.0;
      betay_p[i] = 0.0;
      betaz_p[i] = 0.0;
    }
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::shift") > 2) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      betax_p_p[i] = 0.0;
      betay_p_p[i] = 0.0;
      betaz_p_p[i] = 0.0;
    }
  }
}

/*@@
  @routine    ADMBase_DtLapseZero
  @date       Oct 18 2007
  @author     Erik Schnetter
  @desc
  Scheduled routine to initialise the dtlapse to zero.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_DtLapseZero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
    dtalp[i] = 0.0;
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::dtlapse") > 1) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      dtalp_p[i] = 0.0;
    }
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::dtlapse") > 2) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      dtalp_p_p[i] = 0.0;
    }
  }
}

/*@@
  @routine    ADMBase_DtShiftZero
  @date       Oct 18 2007
  @author     Erik Schnetter
  @desc
  Scheduled routine to initialise the dtshift to zero.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_DtShiftZero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
    dtbetax[i] = 0.0;
    dtbetay[i] = 0.0;
    dtbetaz[i] = 0.0;
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::dtshift") > 1) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      dtbetax_p[i] = 0.0;
      dtbetay_p[i] = 0.0;
      dtbetaz_p[i] = 0.0;
    }
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "ADMBase::dtshift") > 2) {
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
      dtbetax_p_p[i] = 0.0;
      dtbetay_p_p[i] = 0.0;
      dtbetaz_p_p[i] = 0.0;
    }
  }
}

/*@@
  @routine    ADMBase_SetShiftStateOn
  @date       Thu Apr 25 23:12:18 2002
  @author     Tom Goodale
  @desc
  Scheduled routine to set the value of the shift state to on.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_SetShiftStateOn(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *shift_state = 1;
}

/*@@
  @routine    ADMBase_SetShiftStateOff
  @date       Thu Apr 25 23:12:18 2002
  @author     Tom Goodale
  @desc
  Scheduled routine to set the value of the shift state to off.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_SetShiftStateOff(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *shift_state = 0;
}

/*@@
  @routine    ADMBase_SetDtLapseStateOn
  @date       Oct 18 2007
  @author     Erik Schnetter
  @desc
  Scheduled routine to set the value of the dtlapse state to on.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_SetDtLapseStateOn(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *dtlapse_state = 1;
}

/*@@
  @routine    ADMBase_SetDtLapseStateOff
  @date       Oct 18 2007
  @author     Erik Schnetter
  @desc
  Scheduled routine to set the value of the dtlapse state to off.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_SetDtLapseStateOff(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *dtlapse_state = 0;
}

/*@@
  @routine    ADMBase_SetDtShiftStateOn
  @date       Oct 18 2007
  @author     Erik Schnetter
  @desc
  Scheduled routine to set the value of the dtshift state to on.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_SetDtShiftStateOn(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *dtshift_state = 1;
}

/*@@
  @routine    ADMBase_SetDtShiftStateOff
  @date       Oct 18 2007
  @author     Erik Schnetter
  @desc
  Scheduled routine to set the value of the dtshift state to off.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory

 @@*/
void ADMBase_SetDtShiftStateOff(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  *dtshift_state = 0;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
