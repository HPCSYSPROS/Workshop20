/*@@
   @file      SetParams.c
   @date      Tue Jan 12 19:16:38 1999
   @author    Tom Goodale
   @desc
              Callback functions for the parameter file parser.
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <stdlib.h>

#include "cctk_Misc.h"
#include "cctk_Flesh.h"
#include "cctk_WarnLevel.h"
#include "cctk_CommandLine.h"
#include "cctk_Parameter.h"
#include "cctki_Parameter.h"
#include "cctk_ActiveThorns.h"
#include "cctki_ActiveThorns.h"


static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_SetParams_c);

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/
static int ReallySetParameter (const char *parameter, const char *value);


/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/
static int num_0errors = 0; /* number of level 0 errors in parameter file */
static int num_1errors = 0; /* number of level 1 errors in parameter file */


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTKi_SetParameter
   @date       Tue Jan 12 19:25:37 1999
   @author     Tom Goodale
   @desc
               Callback routine for @seeroutine ParseFile.<BR>
               It sets a single parameter to a given value (treating the
               "ActiveThorns" parameter as a special one).
               In case of an error a warning is printed, and the internal
               counter for the appropriate warning level is incremented.
   @enddesc
   @calls      CCTKi_ActivateThorns
               ReallySetParameter

   @var        parameter
   @vdesc      Name of a parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        value
   @vdesc      Value of the parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        lineno
   @vdesc      current line number in the parameter file
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               0 for success, or return code of @seeroutine ReallySetParameter
   @endreturndesc
@@*/
int CCTKi_SetParameter (const char *parameter, const char *value, int lineno)
{
  int retval, parameter_check;
  char parfile[256];


  parameter_check = CCTK_ParameterLevel ();
  CCTK_ParameterFilename (sizeof (parfile), parfile);

  if (CCTK_Equals (parameter, "ActiveThorns"))
  {
    if (CCTKi_ActivateThorns (value))
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                 "CCTKi_SetParameter: Error at line %d in parameter file %s "
                 "while activating thorns",lineno,parfile);
    }
    retval = 0;
  }
  else
  {
    retval = ReallySetParameter (parameter, value);
  }

  if (retval == -1)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "Major error in parameter file '%s' line %d: "
                "Range error setting parameter "
                "'%s' to '%s'", parfile, lineno, parameter, value);
    num_0errors++;
  }
  else if (retval == -2)
  {
    /* Parameter not defined in thorn */
    if (parameter_check == CCTK_PARAMETER_RELAXED)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Minor error in parameter file '%s' line %d: "
                  "Parameter '%s' not found",
                  parfile, lineno, parameter);
      num_1errors++;
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Major error in parameter file '%s' line %d: "
                  "Parameter '%s' not found",
                  parfile, lineno, parameter);
      num_0errors++;
    }
  }
  else if (retval == -4)
  {
    /* Setting parameter twice */
    if (parameter_check == CCTK_PARAMETER_RELAXED)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Minor error in parameter file '%s' line %d: "
                  "Parameter '%s' set in two "
                  "different thorns", parfile, lineno, parameter);
      num_1errors++;
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Major error in parameter file '%s' line %d: "
                  "Parameter '%s' set in two "
                  "different thorns", parfile, lineno, parameter);
      num_0errors++;
    }
  }
  else if (retval == -5)
  {
    /* Parameter not defined by any active thorn */
    if (parameter_check == CCTK_PARAMETER_STRICT)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Major error in parameter file '%s' line %d: "
                  "Parameter '%s' is not "
                  "associated with an active thorn",
                  parfile, lineno, parameter);
      num_0errors++;
    }
    else if (parameter_check == CCTK_PARAMETER_NORMAL)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Minor error in parameter file '%s' line %d: "
                  "Parameter '%s' is not "
                  "associated with an active thorn",
                  parfile, lineno, parameter);
      num_1errors++;
    }
  }
  else if (retval == -6)
  {
    /* Parameter value is not a valid number */
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "Major error in parameter file '%s' line %d: "
                "Error setting parameter '%s' "
                "'%s' is not a valid number",
                parfile, lineno, parameter, value);
    num_0errors++;
  }
  else if (retval == -7)
  {
    /* Tried to set an accumulator parameter directly */
    if (parameter_check == CCTK_PARAMETER_STRICT)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Major error in parameter file '%s' line %d: "
                  "Parameter '%s' is an "
                  "accumulator parameter; it cannot be set directly",
                  parfile, lineno, parameter);
      num_0errors++;
    }
    else if (parameter_check == CCTK_PARAMETER_NORMAL)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Minor error in parameter file '%s' line %d: "
                  "Parameter '%s' is an "
                  "accumulator parameter; it cannot be set directly",
                  parfile, lineno, parameter);
      num_1errors++;
    }
  }
  else if (retval == -8)
  {
    /* Tried to set an array base parameter */
    if (parameter_check == CCTK_PARAMETER_RELAXED)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Minor error in parameter file '%s' line %d: "
                  "Parameter '%s' is an "
                  "array base parameter; please use %s[<number>]",
                  parfile, lineno, parameter, parameter);
      num_1errors++;
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Major error in parameter file '%s' line %d: "
                  "Parameter '%s' is an "
                  "array base parameter; please use %s[<number>]",
                  parfile, lineno, parameter, parameter);
      num_0errors++;
    }
  }
  else if (retval == -10)
  {
    /* parameter has already been set to a different value in the parfile */
    if (parameter_check == CCTK_PARAMETER_RELAXED)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Minor error in parameter file '%s' line %d: "
                  "Parameter '%s' has already been set to a different value "
                  "before", parfile, lineno, parameter);
      num_1errors++;
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Major error in parameter file '%s' line %d: "
                  "Parameter '%s' has already been set to a different value "
                  "before", parfile, lineno, parameter);
      num_0errors++;
    }
  }
  else if (retval == -11)
  {
    /* parameter has already been set to the same value in the parfile */
    if (parameter_check != CCTK_PARAMETER_STRICT)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Minor error in parameter file '%s' line %d: "
                  "Parameter '%s' has already been set to the same value "
                  "before", parfile, lineno, parameter);
      num_1errors++;
    }
    else
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Major error in parameter file '%s' line %d: "
                  "Parameter '%s' has already been set to the same value "
                  "before", parfile, lineno, parameter);
      num_0errors++;
    }
  }

  return (retval);
}


 /*@@
   @routine    CCTKi_NumParameterFileErrors
   @date       Mon Dec 3 2001
   @author     Gabrielle Allen
   @desc
               Returns the number of errors in the parameter file with a
               given level
   @enddesc

   @var        level
   @vdesc      Level of errors to report
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               number of errors in parameter file with this level
   @endreturndesc
@@*/
int CCTKi_NumParameterFileErrors (int level)
{
  int retval;


  if (level == 0)
  {
    retval = num_0errors;
  }
  else if (level == 1)
  {
    retval = num_1errors;
  }
  else
  {
    retval = 0;
  }

  return (retval);
}


/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
 /*@@
   @routine    ReallySetParameter
   @date       Tue Jan 12 19:25:37 1999
   @author     Tom Goodale
   @desc
               Really sets the parameter value.
   @enddesc
   @calls      CCTK_IsImplementationActive
               CCTK_IsThornActive
               CCTK_ActivatingThorn
               CCTK_ParameterSet

   @var        parameter
   @vdesc      Name of a parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        value
   @vdesc      Value of the parameter
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
                0  = success,<BR>
               -10 = parameter is set more than once in the parfile<BR>
               -9  = final value of accumulator out of range<br>
               -6  = not a valid integer or float<br>
               -5  = thorn/imp not active,<BR>
               -4  = tried to set parameter in two different thorns,<BR>
               -3  = tried to steer nonsteerable parameter,<BR>
               -2  = parameter not defined in the active thorn,<BR>
               -1  = parameter out of range
   @endreturndesc
@@*/
static int ReallySetParameter(const char *parameter, const char *value)
{
  int retval;
  int retval_thorn=0;
  int found = 0;
  int retval_imp=0;
  const char *thorn;
  char *param, *imp;


  retval = 0;

  Util_SplitString (&imp, &param, parameter, "::");

  if (! param)
  {
    /* must be global parameter */
    retval = CCTK_ParameterSet (parameter, imp, value);
  }
  else
  {
    /* try and set parameter from implementation */
    if (CCTK_IsImplementationActive (imp))
    {
      thorn = CCTK_ActivatingThorn (imp);
      /* only do it if the thorn name is different to the imp */
      if (! CCTK_Equals (thorn, imp))
      {
        found++;
        retval_imp = CCTK_ParameterSet (param, thorn, value);
      }
    }

    /* try and set parameter from thorn */
    if (CCTK_IsThornActive (imp))
    {
      found++;
      retval_thorn = CCTK_ParameterSet (param, imp, value);
    }

    if (! found)
    {
      /* imp or thorn not found */
      retval = -5;
    }
    else if (found == 2 && retval_imp >= 0 && retval_thorn >= 0)
    {
      /* tried to set parameter for both imp and thorn */
      retval = -4;
    }
    else if (found == 2 && retval_imp < 0 && retval_thorn < 0)
    {
      /* failed to set parameter for both imp and thorn */
      /* FIXME: returning imp but loosing thorn info*/
      retval = retval_imp;
    }
    else if (found == 2)
    {
      /* Only one succeeded */
      retval = retval_imp >= 0 ? retval_imp : retval_thorn;
    }
    else if (found == 1)
    {
      retval = retval_imp + retval_thorn;
    }
  }

  /* Free any allocated memory. */
  free (imp);
  free (param);

  return (retval);
}
