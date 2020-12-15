/*@@
   @file      Parameters.c
   @date      Mon Jun 28 21:44:17 1999
   @author    Tom Goodale
   @desc
              Routines to deal with the parameters.
   @enddesc
   @version   $Id$
 @@*/

#include <sys/types.h>
#include <regex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <assert.h>

#include "cctk_Flesh.h"
#include "cctk_ActiveThorns.h"
#include "cctk_Constants.h"
#include "cctk_WarnLevel.h"
#include "cctk_Parameter.h"
#include "cctk_GNU.h"
#include "cctk_FortranString.h"

#include "StoreHandledData.h"
#include "util_String.h"
#include "util_Expression.h"

#include "cctki_Parameter.h"

#include "ParameterBindings.h"

static const char *rcsid="$Header$";

CCTK_FILEVERSION(main_Parameters_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/* what is a parameter:
 * - properties
 * - data
 * nothing else!
 * Props are supposed to be public to thorn programmers,
 * Data are supposed to be private to flesh!
 */
typedef struct PARAM
{
  cParamData *props;
  void       *data;

  int n_accumulator_sources;
  struct PARAM **accumulates_from;

  int n_accumulator_bases;
  struct PARAM **accumulator_bases;

  struct PARAM *array;
} t_param;

/* what is a list of parameters:
 * just a normal linked list...
 */
typedef struct PARAMLIST
{
  t_param          *param;
  struct PARAMLIST *last;
  struct PARAMLIST *next;
} t_paramlist;

/* what is a node of parameter tree
 * look at it: a list of linked lists... ;-)
 */
typedef struct PARAMTREENODE
{
  t_paramlist *paramlist;
} t_paramtreenode;

/* structure to describe a registered parameter set notify callback */
typedef struct NOTIFYCALLBACK
{
  cParameterSetNotifyCallbackFn fn;
  void *data;
  regex_t thorn_preg,
          param_preg;
  int has_thorn_regex,
      has_param_regex;
} t_notifyCallback;

/********************************************************************
 ************************* Static Variables *************************
 ********************************************************************/

/* mask for CCTK_ParameterSet() */
static int cctk_parameter_set_mask;

/* list of registered parameter set notify callbacks */
static cHandledData *notifyCallbackList = NULL;

/* number of registered callbacks */
static int num_notify_callbacks = 0;

/********************************************************************
 ********************* Fortran Wrapper Prototypes *******************
 ********************************************************************/

void CCTK_FCALL CCTK_FNAME (CCTK_ParameterSetNotifyRegister)
                           (int *status, cParameterSetNotifyCallbackFn callback,
                            void *data, THREE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_ParameterSetNotifyUnregister)
                           (int *status, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_ParameterValString)
                           (int *nchars, THREE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_ParameterSet)
                           (int *status, THREE_FORTSTRING_ARG);

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static t_param *ParameterFind (const char *name,
                               const char *thorn,
                               int scope);

static t_param *ParameterNew (const char *thorn,
                              const char *name,
                              const char *type,
                              const char *scope,
                              int        steerable,
                              const char *description,
                              const char *defval,
                              void       *data,
                              int         arraysize,
                              const char *accumulator_expression);

static cParamData *ParamDataNew(char *thorn,
                                char *name,
                                char *description,
                                char *defval,
                                int   scope,
                                int   type,
                                int   steerable,
                                int   array_size,
                                int   array_index,
                                char *accumulator_expression);

static const void *ParameterGetSimple (const t_param *param, int *type);

static int ParameterSet(t_param *param, const char *value);

static int ParameterSetSimple (t_param *param, const char *value);

static t_paramtreenode *ParameterPTreeNodeFind (t_sktree *tree,
                                                const char *name);

static int ParameterGetScope (const char *scope);
static int ParameterGetType (const char *type);

static int ParameterInsert (t_sktree **tree, t_param *newparam);

#ifdef TEST_PARAMETERS
static int ParameterPrintSimple (t_param *param,
                                 const char *format,
                                 FILE *file);
#endif

static int ParameterExtend (t_param *param,
                            const char *range_origin,
                            const char *range,
                            const char *range_description);

static int ParameterListAddParam (t_paramlist **paramlist,
                                  t_param *newparam);


static int ParameterSetAccumulator(t_param *param);

static int ParameterSetKeyword  (t_param *param, const char *value);
static int ParameterSetString   (t_param *param, const char *value);
static int ParameterSetSentence (t_param *param, const char *value);
static int ParameterSetInteger  (t_param *param, const char *value);
static int ParameterSetReal     (t_param *param, const char *value);
static int ParameterSetBoolean  (t_param *param, const char *value);
static int SetVarEvaluator(int nvars, const char * const *vars, uExpressionValue *vals, const void *data);

static void GetBaseName(const char *name, char **basename, int *array_index);
static char *ArrayParamName(const char *basename,int array_index);
static int AccVarEvaluator(int nvars, const char * const *vars, uExpressionValue *vals, const void *data);
static void AddAccumulators(t_param *base,
                            t_param *extra,
                            const char *thorn,
                            const char *parameter,
                            const char *imp,
                            const char *baseparam);
static int LinkAccumulators(t_param *base, t_param *extra);
static void ParameterActivate(t_param *param);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

extern void CCTKi_SetParameterSetMask (int mask);

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static t_sktree *paramtree = NULL;

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/


/*@@
   @routine    CCTKi_ParameterCreate
   @date       Tue Jun 29 10:08:56 1999
   @author     Tom Goodale
   @desc
               Creates a parameter originating from a thorn/implementation.
   @enddesc
   @calls      ParameterFind
               ParameterNew
               ParameterExtend
               ParameterSetSimple

   @var        name
   @vdesc      The name of the parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        thorn
   @vdesc      The originating thorn
   @vtype      const char *
   @vio        in
   @endvar
   @var        scope
   @vdesc      The scope of the parameter
   @vtype      const char *
   @vio        in
   @vcomment
               Private, Restricted, or Global
   @endvar
   @var        steerable
   @vdesc      Is the parameter steerable ?
   @vtype      int
   @vio        in
   @endvar
   @var        description
   @vdesc      A description of the parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        defval
   @vdesc      The default value
   @vtype      const char *
   @vio        in
   @endvar
   @var        data
   @vdesc      A pointer to the memory holding the parameter
   @vtype      void *
   @vio        in
   @endvar
   @var        array
   @vdesc      array size
   @vtype      int
   @vio        in
   @vcomment
     0 means no array
   @endvar
   @var        accumulator
   @vdesc      accumulator expression
   @vtype      const char *
   @vio        in
   @vcomment
     NULL means not an accumulator
   @endvar
   @var        n_ranges
   @vdesc      number of basic ranges
   @vtype      int
   @vio        in
   @endvar
   @var        ...
   @vdesc      range data
   @vtype      variable argument list of const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
                0 on success, or<br>
               -1 on failure to create the parameter<br>
               -2 if parameter inconsistent
   @endreturndesc
@@*/
int CCTKi_ParameterCreate (const char *name,
                           const char *thorn,
                           const char *type,
                           const char *scope,
                           int        steerable,
                           const char *description,
                           const char *defval,
                           void       *data,
                           int        array,
                           const char *accumulator,
                           int n_ranges,
                           ...)
{
  int i, retval;
  t_param *param;
  va_list ranges;
  const char *rangeval, *rangedesc;

  param = ParameterFind (name, thorn, ParameterGetScope (scope));

  if (! param)
  {
    param = ParameterNew (thorn, name, type, scope, steerable, description,
                          defval, data, array, accumulator);
    if (n_ranges)
    {
      va_start (ranges, n_ranges);

      for (i = 0; i < n_ranges; i++)
      {
        rangeval = (const char *) va_arg (ranges, const char *);
        rangedesc = (const char *) va_arg (ranges, const char *);

        ParameterExtend (param, thorn, rangeval, rangedesc);
      }
      va_end (ranges);
    }

    if(array)
    {
      retval = 0;
      for(i = 0; i < array; i++)
      {
        /* Setup ranges on the array element */
        param->array[i].props->range = param->props->range;
        /* Set its default value */
        retval = ParameterSetSimple (&(param->array[i]), defval);
      }
    }
    else
    {
      retval = ParameterSetSimple (param, defval);
    }
  }
  else
  {
    retval = -2;
  }

  return (retval);
}


/*@@
   @routine    CCTKi_ParameterAddRange
   @date       Tue Jun 29 10:15:53 1999
   @author     Tom Goodale
   @desc
               Adds a range.
               Only allowed to add a range if in appropriate scope.
   @enddesc
   @calls      CCTK_ImpThornList
               SKTreeFindFirst
               ParameterFind
               ParameterExtend

   @var        origin
   @vdesc      The originating implementation or thorn of the parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        name
   @vdesc      The name of the parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        range_origin
   @vdesc      The originating implementation or thorn of the range
   @vtype      const char *
   @vio        in
   @endvar
   @var        range
   @vdesc      The new range
   @vtype      const char *
   @vio        in
   @endvar
   @var        range_description
   @vdesc      A description of the new range
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
                0 on success, or<br>
               -1 if no such parameter exists
   @endreturndesc
@@*/
int CCTKi_ParameterAddRange (const char *implementation,
                             const char *name,
                             const char *range_origin,
                             const char *range,
                             const char *range_description)
{
  int retval;
  t_param *param;
  /* For the moment do this in the quick and dirty way 8-(  FIXME */
  t_sktree *thornlist, *node;


  retval = -1;

  thornlist = CCTK_ImpThornList (implementation);
  if (thornlist)
  {
    for (node = SKTreeFindFirst (thornlist); node; node = node->next)
    {
      param = ParameterFind (name, node->key, SCOPE_RESTRICTED);
      if (param)
      {
        retval = ParameterExtend (param, range_origin, range,range_description);
      }
      else
      {
        retval = -1;
      }
    }
  }

  return (retval);
}

 /*@@
   @routine    CCTKi_ParameterAccumulatorBase
   @date       Mon May 20 02:57:55 2002
   @author     Tom Goodale
   @desc
   Sets the accumulator base for a parameter.
   @enddesc

   @var     thorn
   @vdesc   Name of thorn providing an extra value to accumulator
   @vtype   const char *
   @vio     in
   @endvar
   @var     thorn
   @vdesc   Name of parameter providing an extra value to accumulator
   @vtype   const char *
   @vio     in
   @vcomment
   @endvar
   @var     importhorn
   @vdesc   Name of thorn or implementation providing accumulator
   @vtype   const char *
   @vio     in
   @endvar
   @var     thorn
   @vdesc   Name of parameter providing accumulator
   @vtype   const char *
   @vio     in
   @vcomment
   @endvar

 @@*/
void CCTKi_ParameterAccumulatorBase(const char *thorn,
                                    const char *parameter,
                                    const char *importhorn,
                                    const char *acc_base)
{
  t_param *param;
  t_param *accumulator_base;

  param = ParameterFind (parameter, thorn, SCOPE_ANY);


  if(param)
  {
    if(CCTK_Equals(thorn, importhorn))
    {
      accumulator_base = ParameterFind(acc_base, thorn, SCOPE_ANY);

      if(accumulator_base)
      {
        AddAccumulators(accumulator_base, param, thorn,parameter,thorn,acc_base);
      }
    }
    else
    {
      t_sktree *thornlist;
      t_sktree *node;

      thornlist = CCTK_ImpThornList (importhorn);

      for (node = SKTreeFindFirst (thornlist); node; node = node->next)
      {
        accumulator_base = ParameterFind(acc_base, node->key, SCOPE_RESTRICTED);

        if(accumulator_base)
        {
          AddAccumulators(accumulator_base, param, thorn,parameter,node->key,acc_base);
        }
      }
    }
  }
}

/*@@
   @routine    CCTK_ParameterSet
   @date       Tue Jun 29 10:22:22 1999
   @author     Tom Goodale
   @desc
               Sets the value (checks for steerable if not initialisation).
   @enddesc
   @calls      ParameterFind
               ParameterSet

   @var        name
   @vdesc      The name of the parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        thorn
   @vdesc      The originating thorn
   @vtype      const char *
   @vio        in
   @endvar
   @var        value
   @vdesc      The value of the parameter
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
                0  for success, or<br>
               return code of @seeroutine ParameterSet, or<br>
               -1  if parameter is out of range<br>
               -2  if parameter was not found<br>
               -3  if trying to steer a non-steerable parameter<br>
               -6  if not a valid integer or float<br>
               -7  if tried to set an accumulator parameter directly<br>
               -8  if tried to set an accumulator parameter directly<br>
               -9  if final value of accumulator out of range<br>
               -10 if parameter has already been set to a different value<br>
               -11 if parameter has already been set to the same value<br>
   @endreturndesc
@@*/
int CCTK_ParameterSet (const char *name, const char *thorn, const char *value)
{
  int handle, retval = 0;
  char *old_value, *new_value;
  t_param *param;
  const t_notifyCallback *cb;


  param = ParameterFind (name, thorn, SCOPE_ANY);

  if (param)
  {
    old_value = CCTK_ParameterValString (param->props->name,
                                         param->props->thorn);

    if(param->array)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "CCTK_ParameterSet: Cannot set base array parameter '%s::%s' "
                  , thorn, name);
      retval = -8;
    }
    else if(param->props->accumulator_expression)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "CCTK_ParameterSet: Cannot set accumulator parameter '%s::%s'"
                  " directly", thorn, name);
      retval = -7;
    }
    /* before parameter recovery (which is while parsing the parameter file)
       all parameters can be set */
    else if (cctk_parameter_set_mask == PARAMETER_RECOVERY_POST &&
             param->props->steerable != CCTK_STEERABLE_ALWAYS)
    {
      /* after parameter recovery only steerable parameters can be set */
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "CCTK_ParameterSet: Cannot set parameter '%s::%s' "
                  "(non-steerable)", thorn, name);
      retval = -3;
    }
    else if (cctk_parameter_set_mask == PARAMETER_RECOVERY_IN &&
             param->props->n_set > 0)
    {
      /* during parameter recovery STEERABLE_NEVER parameters which were set
       from the parameter file are overwritten by the checkpoint file */
      if (param->props->steerable == CCTK_STEERABLE_NEVER)
      {
        if (strcmp (old_value, value))
        {
          CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                      "CCTK_ParameterSet: Non-steerable parameter '%s::%s' "
                      "cannot be set from the parameter file but is recovered "
                      "from the checkpoint file",
                      thorn, name);
          /* we ignore a possible error code from ParameterSet, instead always
           * returning "this parameter has already been set" to our caller.
           * We expect the caller to terminate, but in case they do not, we use
           * the value of the parameter recovered from the checkpoint file.
           */
          ParameterSet (param, value);
          retval = -10;
        }
        else
        {
          /* unregister the previous set operation
             if the restored value is the same as the one previously set */
          param->props->n_set--;
        }
      }
      else
      {
        /* do not restore the original value
           but register another set operation */
        retval = 0;
        param->props->n_set++;
      }
    }
    else
    {
      /* call parameter set notify listeners */
      if (cctk_parameter_set_mask == PARAMETER_RECOVERY_POST)
      {
        for (handle = 0; handle < num_notify_callbacks; handle++)
        {
          cb = (const t_notifyCallback *)
               Util_GetHandledData (notifyCallbackList, handle);
          if (cb)
          {
            if (cb->has_thorn_regex && regexec (&cb->thorn_preg, thorn, 0,
                                                NULL, 0))
            {
              continue;
            }
            if (cb->has_param_regex && regexec (&cb->param_preg, name, 0,
                                                NULL, 0))
            {
              continue;
            }
            cb->fn (cb->data, thorn, name, value);
          }
        }
      }

      retval = ParameterSet (param, value);

      new_value = CCTK_ParameterValString (param->props->name,
                                           param->props->thorn);

      /* check if a parameter is set more than once in a parfile */
      if (cctk_parameter_set_mask == PARAMETER_RECOVERY_PRE &&
          param->props->n_set > 0)
      {
        if (retval == 0)
        {
          retval = strcmp (old_value, new_value) ? -10 : -11;
        }
      }

      /* register another set operation
         if the parameter has not simply been restored to the same value */
      if (! (cctk_parameter_set_mask == PARAMETER_RECOVERY_IN &&
             strcmp (old_value, new_value) == 0))
      {
        param->props->n_set++;
      }

      free (new_value);
    }

    free (old_value);
  }
  else
  {
    retval = -2;
  }

  return (retval);
}


void CCTK_FCALL CCTK_FNAME (CCTK_ParameterSet)
                           (int *status, THREE_FORTSTRING_ARG)
{
  THREE_FORTSTRING_CREATE (name, thorn, value)

  *status = CCTK_ParameterSet (name, thorn, value);

  free (name);
  free (thorn);
  free (value);
}


/*@@
   @routine    CCTK_ParameterGet
   @date       Tue Jun 29 10:28:20 1999
   @author     Tom Goodale
   @desc
               Gets the pointer to the parameter.
               Should be used for checkpointing and recovery.
   @enddesc
   @calls      ParameterFind
               ParameterGetSimple

   @var        name
   @vdesc      The name of the parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        thorn
   @vdesc      The originating thorn
   @vtype      const char *
   @vio        in
   @endvar
   @var        type
   @vdesc      If not NULL, a pointer to an integer which will hold the type
               of the parameter
   @vtype      int *
   @vio        out
   @endvar

   @returntype const void *
   @returndesc
               The address of the parameter, or<br>
               NULL if no parameter by that name was found
   @endreturndesc
@@*/
const void *CCTK_ParameterGet (const char *name, const char *thorn, int *type)
{
  const void *retval;
  t_param *param;


  param = ParameterFind (name, thorn, SCOPE_ANY);
  if (param)
  {
    retval = ParameterGetSimple (param, type);
  }
  else
  {
    retval = NULL;
  }

  return (retval);
}



/*@@
   @routine    CCTK_ParameterQueryTimesSet
   @date       Fri July 7 2000
   @author     Gabrielle Allen
   @desc
               Returns the number of times that a parameter has been set,
               that is if it returns 0 the parameter was not set in a
               parameter file.
   @enddesc
   @calls      CCTK_ParameterData

   @var        name
   @vdesc      The name of the parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        thorn
   @vdesc      The originating thorn
   @vtype      const char *
   @vio        in
   @endvar

   @returntype int
   @returndesc
               number of times set, or<br>
               -1 if no parameter by that name exists
   @endreturndesc
@@*/
int CCTK_ParameterQueryTimesSet (const char *name, const char *thorn)
{
  const cParamData *param;


  param = CCTK_ParameterData (name, thorn);

  return (param ? param->n_set : -1);
}


/*@@
   @routine    CCTK_ParameterValString
   @date       Thu Jan 21 2000
   @author     Thomas Radke
   @desc
               Gets the string representation of a parameter's value
               - can be used for checkpointing and recovery.
   @enddesc
   @calls      CCTK_ParameterGet

   @var        param_name
   @vdesc      The name of the parameter
   @vtype      const char *
   @vio        in
   @endvar
   @var        thorn
   @vdesc      The originating thorn
   @vtype      const char *
   @vio        in
   @endvar

   @returntype char *
   @returndesc
               the allocated string (should be freed after usage), or<br>
               NULL if no parameter by that name was found
   @endreturndesc
@@*/
char *CCTK_ParameterValString (const char *param_name, const char *thorn)
{
  int param_type;
  const void *param_data;
  char *retval;
  char buffer[80];


  retval = NULL;

  param_data = CCTK_ParameterGet (param_name, thorn, &param_type);
  if (param_data)
  {
    switch (param_type)
    {
      case PARAMETER_KEYWORD:
      case PARAMETER_STRING:
      case PARAMETER_SENTENCE:
        retval = strdup (*(const char *const *) param_data);
        break;

      case PARAMETER_BOOLEAN:
        retval = strdup ((int) (*(const CCTK_INT *) param_data) ? "yes" : "no");
        break;

      case PARAMETER_INT:
        sprintf (buffer, "%d", (int) (*(const CCTK_INT *) param_data));
        retval = strdup (buffer);
        break;

      case PARAMETER_REAL:
        sprintf (buffer, "%.20g", (double) (*(const CCTK_REAL *) param_data));
        retval = strdup (buffer);
        break;

      default:
        CCTK_VWarn (3, __LINE__, __FILE__, "Cactus",
                    "CCTK_ParameterValString: Unknown type %d for parameter "
                    "'%s::%s'", param_type, thorn, param_name);
        break;
    }
  }

  return (retval);
}

/*@@
   @routine    CCTK_PARAMETERVALSTRING
   @date       Thu Jan 21 2000
   @author     Thomas Radke
   @desc
               Copies the stringified parameter value into a fortran string.
   @enddesc
   @calls      CCTK_ParameterValString

   @var        nchars
   @vdesc      Number of characters in the stringified parameter value, or<BR>
               -1 if the parameter doesn't exist
   @vtype      CCTK_INT *
   @vio        out
   @vcomment
               It will copy only as many characters as fit into the fortran
               string. You should check for truncation by comparing 'nchars'
               against the length of your fortran string.
   @endvar
   @var        THREE_FORTSTRING_ARG
   @vdesc      three fortran string arguments denoting parameter name,
               thorn name, and the string buffer to copy the parameter value to
   @vtype      FORTRAN string macro
   @vio        inout
   @endvar
@@*/
void CCTK_FCALL CCTK_FNAME (CCTK_ParameterValString)
                           (int *nchars, THREE_FORTSTRING_ARG)
{
  size_t c_strlen;
  char *c_string;
  THREE_FORTSTRING_PTR (unused1, unused2, fortran_string)
  THREE_FORTSTRING_CREATE (param, thorn, value)


  /* get rid of compiler warnings about unused variables */
  unused1 = unused1;
  unused2 = unused2;

  c_string = CCTK_ParameterValString (param, thorn);
  if (c_string)
  {
    *nchars = c_strlen = strlen (c_string);
    if (c_strlen > (size_t) cctk_strlen3)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "CCTK_ParameterValString: fortran string buffer is too short "
                  "to hold value '%s' of parameter '%s::%s', string will be "
                  "truncated", c_string, thorn, param);
      c_strlen = (size_t) cctk_strlen3;
    }
    /* copy up to the size of the fortran string
       and pad remaining chars in the fortran string with spaces */
    memcpy (fortran_string, c_string, c_strlen);
    memset (fortran_string + c_strlen, ' ', (size_t) cctk_strlen3 - c_strlen);
  }
  else
  {
    /* no such parameter exists */
    *nchars = -1;
  }

  free (param);
  free (thorn);
  free (value);
}

/*@@
   @routine    CCTK_ParameterWalk
   @date       Tue Jun 29 10:30:21 1999
   @author     Tom Goodale
   @desc
               Gets parameters in order, restricted to ones from 'origin',
               or all if 'origin' is NULL.  Starts with the first parameter
               if 'first' is true, otherwise gets the next one.
               Can be used for generating full help file, or for walking
               the list and checkpointing.
   @enddesc
   @calls      CCTK_ThornImplementation
               SKTreeFindFirst

   @var        first
   @vdesc      Flag to indicate get first parameter or not
   @vtype      int
   @vio        in
   @endvar
   @var        origin
   @vdesc      The thorn to walk, or NULL if to walk all params
   @vtype      const char *
   @vio        in
   @endvar
   @var        pfullname
   @vdesc      address of pointer to string where the fullname of the parameter
               will be stored
   @vtype      char *
   @vio        out
   @vcomment
               If not NULL, the full name is written into an allocated string
               which can be refered to via *pfullname. Should be freed after
               usage.  Full name is "implementation::name" for global and
               restricted parameters, and "thorn::name" for private parameters.
   @endvar
   @var        pdata
   @vdesc      address of pointer to parameter data structure
   @vtype      cParamData **
   @vio        out
   @vcomment
               If not NULL, the pointer to the parameter data structure
               is stored in *pdata.
   @endvar

   @returntype int
   @returndesc
            - zero for success
            - positive if parameter was not found
            - negative if initial startpoint was not set
   @endreturndesc

@@*/
int CCTK_ParameterWalk (int first,
                        const char *origin,
                        char **pfullname,
                        const cParamData **pdata)
{
  int             return_found;
  const char     *prefix;
  t_sktree        *tnode;
  t_paramtreenode *node;
  t_paramlist     *paramlist;
  t_param         *startpoint;
  static t_param  *prev_startpoint_all = NULL;
  static t_param  *prev_startpoint_thorn = NULL;
  static int      next_index = 0;
  /* FIXME : This routine has become extremely ugly:
   *         It should only have one return in it.
   *         The malloc failure should be flagged.
   */

  /* determine the startpoint for search */
  if (! first)
  {
    startpoint = origin ? prev_startpoint_thorn : prev_startpoint_all;

    if (startpoint == NULL)
    {
      CCTK_Warn (2, __LINE__, __FILE__, "Cactus",
                 "CCTK_ParameterWalk: Cannot walk through parameter list "
                 "without setting a startpoint at first");
      return (-1);
    }
    if (startpoint->props->array_size > 0 &&
        startpoint->props->array_size > next_index)
    {
      /* return next (array parameter) match at this startpoint */
      return_found = 1;
    }
    else
    {
      /* return next match AFTER startpoint */
      next_index = 0;
      return_found = 0;
    }
  }
  else
  {
    /* return next match which also becomes startpoint for following searches */
    return_found = 1;
    startpoint = NULL;
  }

  /* iterate over nodes */
  for (tnode = SKTreeFindFirst (paramtree) ; tnode ; tnode = tnode->next)
  {
    /* get node data */
    node  = (t_paramtreenode *) tnode->data;

    /* iterate over parameters in list */
    for (paramlist = node->paramlist; paramlist; paramlist = paramlist->next)
    {
      /* if startpoint is still unassigned set it to next match in list */
      if (startpoint == NULL)
      {
        if (! origin || CCTK_Equals (origin, paramlist->param->props->thorn))
        {
          startpoint = paramlist->param;
        }
      }

      /* Hey, we've found the startpoint ! */
      if (startpoint == paramlist->param)
      {
        /* Do we have to return this one ?
           If not prepare finding the next matching param. */
        if (return_found)
        {
          /* save the last startpoint */
          prev_startpoint_all = prev_startpoint_thorn = startpoint;

          /* advance to the next if this is an array parameter */
          if (startpoint->props->array_size > 0)
          {
            startpoint = &startpoint->array[next_index++];
          }

          /* set the parameter's fullname if requested by the caller */
          if (pfullname)
          {
            prefix = startpoint->props->thorn;
            if (startpoint->props->scope != SCOPE_PRIVATE)
            {
              prefix = CCTK_ThornImplementation (prefix);
            }

            *pfullname = malloc (strlen (prefix) +
                                 strlen (startpoint->props->name) + 3);
            if(*pfullname)
            {
                sprintf (*pfullname, "%s::%s",
                         prefix, startpoint->props->name);
            }
          }

          if (pdata)
          {
            *pdata = startpoint->props;
          }

          return (0);
        }
        else
        {
          startpoint = NULL;
          return_found = 1;
        }
      }
    } /* end looping over parameter list */
  } /* end looping over all nodes */

  return (1);
}

 /*@@
   @routine    CCTKi_ParameterActivateThornParameters
   @date       Mon May 20 07:00:59 2002
   @author     Tom Goodale
   @desc
   Does any activations necessary for a thorn's parameters
   @enddesc

   @var     thorn
   @vdesc   The thorn who's parameters are to be activated
   @vtype   const char *
   @vio     in
   @endvar

 @@*/
void CCTKi_ParameterActivateThornParameters(const char *thorn)
{
  t_sktree        *tnode;
  t_paramtreenode *node;
  t_paramlist     *paramlist;
  t_param         *current;

  for (tnode = SKTreeFindFirst (paramtree) ; tnode ; tnode = tnode->next)
  {
    /* get node data */
    node  = (t_paramtreenode *) tnode->data;

    /* iterate over parameters in list */
    for (paramlist = node->paramlist; paramlist; paramlist = paramlist->next)
    {
      current = paramlist->param;

      if(CCTK_Equals(current->props->thorn, thorn))
      {
        ParameterActivate(current);
      }
    }
  }
}

/*@@
  @routine    CCTK_ParameterData
  @date       Tue Aug 31 18:10:46 MSZ 1999
  @author     Andre Merzky
  @desc
              For a given parameter, return description, type and range.
  @enddesc
  @calls      ParameterFind

  @var        name
  @vdesc      parameter name
  @vtype      const char *
  @vio        in
  @endvar
  @var        thorn
  @vdesc      thorn parameter belongs to
  @vtype      const char *
  @vio        in
  @endvar

  @returntype const cParamData *
  @returndesc
              pointer to parameter data structure on success,<br>
              NULL on failure.
  @endreturndesc
@@*/
const cParamData *CCTK_ParameterData (const char *name,
                                      const char *thorn)
{
  const t_param *param;


  param = ParameterFind (name, thorn, SCOPE_ANY);

  return (param ? param->props : NULL);
}


/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    ParameterFind
   @date       Sun Oct 17 16:20:48 1999
   @author     Tom Goodale
   @desc
               Finds a parameter given by its name, thorn, and scope.
   @enddesc
   @calls      ParameterPTreeNodeFind
@@*/
static t_param *ParameterFind (const char *name,
                               const char *thorn,
                               int scope)
{
  t_param *retval;
  t_paramtreenode *node;
  t_paramlist *list;

  char *basename;
  int   array_index;

  /* Split an array parameter into its name and index */
  GetBaseName(name, &basename, &array_index);

  /* Find the parameter */
  node = ParameterPTreeNodeFind (paramtree, basename);

  free(basename);

  list = NULL;

  if (node)
  {
    /* Parameter exists for some thorn;  check thorn */
    for (list = node->paramlist; list; list = list->next)
    {
      if (! thorn)
      {
        if (list->param->props->scope == SCOPE_GLOBAL)
        {
          break;
        }
      }
      else if (scope == SCOPE_ANY)
      {
        if (thorn && ! Util_StrCmpi (thorn, list->param->props->thorn))
        {
          break;
        }
      }
      else if (! Util_StrCmpi (thorn, list->param->props->thorn) &&
               list->param->props->scope == scope)
      {
        break;
      }
    }
  }

  /* Now get the correct parameter structure */
  if(list)
  {
    if(!list->param->array && array_index > -1)
    {
      /* Trying to treat a non-array parameter as an array */
      retval = NULL;
    }
    else if(list->param->array && array_index > -1)
    {
      /* It's an array parameter */
      if(array_index < list->param->props->array_size)
      {
        retval = &(list->param->array[array_index]);
      }
      else
      {
        /* It's out of bounds */
        retval = NULL;
      }
    }
    else
    {
      /* Just a normal parameter */
      retval = list->param;
    }
  }
  else
  {
    retval = NULL;
  }

  return retval;
}


/*@@
   @routine    ParameterNew
   @date       Mon Jul 26 10:59:42 1999
   @author     Tom Goodale
   @desc
               Creates a new parameter.
   @enddesc
   @calls      ParameterGetScope
               ParameterGetType
               ParameterInsert
   @history

   @endhistory
   @var     thorn
   @vdesc   Name of the thorn
   @vtype   const char *
   @vio     in
   @endvar
   @var     name
   @vdesc   Name of the parameter
   @vtype   const char *
   @vio     in
   @endvar
   @var     type
   @vdesc   Parameter type
   @vtype   int
   @vio     in
   @endvar
   @var     scope
   @vdesc   Parameter scope
   @vtype   int
   @vio     in
   @endvar
   @var     steerable
   @vdesc   Is the parameter steerable
   @vtype   int
   @vio     in
   @endvar
   @var     description
   @vdesc   Description of the parameter
   @vtype   const char *
   @vio     in
   @endvar
   @var     defval
   @vdesc   Default value of the parameter
   @vtype   const char *
   @vio     in
   @endvar
   @var     data
   @vdesc   Pointer to parameter data
   @vtype   void *
   @vio     inout
   @endvar
   @var     arraysize
   @vdesc   Size of parameter array, if any
   @vtype   int
   @vio     in
   @endvar
   @var     accumulator_expression
   @vdesc   expression to accumulate on
   @vtype   const char *
   @vio     in
   @endvar

   @returntype t_param *
   @returndesc
    The new parameter data structure.
   @endreturndesc
@@*/
static t_param *ParameterNew (const char *thorn,
                              const char *name,
                              const char *type,
                              const char *scope,
                              int         steerable,
                              const char *description,
                              const char *defval,
                              void       *data,
                              int         arraysize,
                              const char *accumulator_expression)
{
  t_param *newparam;
  int i;
  char  *newname;

  newparam = (t_param *) malloc (sizeof (t_param));
  if (newparam)
  {
    newparam->props = (cParamData *) malloc (sizeof (cParamData));
    if (newparam->props)
    {
      newparam->props->thorn       = Util_Strdup (thorn);
      newparam->props->name        = Util_Strdup (name);
      newparam->props->description = Util_Strdup (description);
      newparam->props->defval      = Util_Strdup (defval);
      newparam->props->scope       = ParameterGetScope(scope);
      newparam->props->type        = ParameterGetType(type);
      newparam->props->steerable   = steerable;
      newparam->props->range       = NULL;
      newparam->props->n_set       = 0;
      newparam->props->array_size  = arraysize;
      newparam->props->array_index = -1;

      if(accumulator_expression)
      {
        newparam->props->accumulator_expression = Util_Strdup(accumulator_expression);
      }
      else
      {
        newparam->props->accumulator_expression = NULL;
      }

      newparam->n_accumulator_sources = 0;
      newparam->accumulates_from = NULL;
      newparam->n_accumulator_bases = 0;
      newparam->accumulator_bases = NULL;
      newparam->data = data;

      ParameterInsert (&paramtree, newparam);

      if(arraysize == 0)
      {
        if (newparam->props->type == PARAMETER_STRING ||
            newparam->props->type == PARAMETER_SENTENCE ||
            newparam->props->type == PARAMETER_KEYWORD)
        {
          *(char **) data = NULL;
        }
        newparam->array = NULL;
      }
      else
      {

        /* It's an array parameter, so setup a new t_param for each array element. */
        newparam->array = (t_param *)malloc(arraysize*sizeof(t_param));

        if(newparam->array)
        {
          for(i = 0; i < arraysize; i++)
          {
            newname = ArrayParamName(newparam->props->name,i);

            newparam->array[i].props = ParamDataNew(newparam->props->thorn,
                                                    newname,
                                                    newparam->props->description,
                                                    newparam->props->defval,
                                                    newparam->props->scope,
                                                    newparam->props->type,
                                                    newparam->props->steerable,
                                                    newparam->props->array_size,
                                                    i,
                                                    newparam->props->accumulator_expression);

            newparam->array[i].n_accumulator_sources = 0;
            newparam->array[i].accumulates_from      = NULL;
            newparam->array[i].n_accumulator_bases   = 0;
            newparam->array[i].accumulator_bases     = NULL;

            newparam->array[i].array            = NULL;

            switch(newparam->props->type)
            {
              case PARAMETER_BOOLEAN : /*Fall through */
              case PARAMETER_INT  : newparam->array[i].data = &(((CCTK_INT *)data)[i]);
                                    break;
              case PARAMETER_REAL : newparam->array[i].data = &(((CCTK_REAL *)data)[i]);
                                    break;
              default :
                /* All remaining types are strings */
                newparam->array[i].data = &(((CCTK_CHAR **)data)[i]);
                *(char **) (newparam->array[i].data) = NULL;
            }
          }
        }
      }
    }
  }

  return (newparam);
}

 /*@@
   @routine    ParamDataNew
   @date       Mon May 20 07:15:40 2002
   @author     Tom Goodale
   @desc
   Creates a new cParamData structure, just setting
   the values to its inputs rather than duplicating them.
   @enddesc

   @var     thorn
   @vdesc   Name of the thorn
   @vtype   char *
   @vio     in
   @endvar
   @var     name
   @vdesc   Name of the parameter
   @vtype   char *
   @vio     in
   @endvar
   @var     description
   @vdesc   Description of the parameter
   @vtype   char *
   @vio     in
   @endvar
   @var     defval
   @vdesc   Default value of the parameter
   @vtype   char *
   @vio     in
   @endvar
   @var     scope
   @vdesc   Parameter scope
   @vtype   int
   @vio     in
   @endvar
   @var     type
   @vdesc   Parameter type
   @vtype   int
   @vio     in
   @endvar
   @var     steerable
   @vdesc   Is the parameter steerable
   @vtype   int
   @vio     in
   @endvar
   @var     array_size
   @vdesc   Size of parameter array, if any
   @vtype   int
   @vio     in
   @endvar
   @var     array_index
   @vdesc   Index into parameter array
   @vtype   int
   @vio     in
   @endvar
   @var     accumulator_expression
   @vdesc   expression to accumulate on
   @vtype   const char *
   @vio     in
   @endvar

   @returntype cParamData *
   @returndesc
    The new data structure.
   @endreturndesc
 @@*/
static cParamData *ParamDataNew(char *thorn,
                                char *name,
                                char *description,
                                char *defval,
                                int   scope,
                                int   type,
                                int   steerable,
                                int   array_size,
                                int   array_index,
                                char *accumulator_expression)
{
  cParamData *props;

  props = (cParamData *) malloc (sizeof (cParamData));

  if(props)
  {
    props->thorn       = thorn;
    props->name        = name;
    props->description = description;
    props->defval      = defval;
    props->scope       = scope;
    props->type        = type;
    props->steerable   = steerable;
    props->range       = NULL;
    props->n_set       = 0;
    props->array_size  = array_size;
    props->array_index = array_index;

    props->accumulator_expression = accumulator_expression;
  }

  return props;
}

static t_paramtreenode *ParameterPTreeNodeFind (t_sktree *tree,
                                                const char *name)
{
  const t_sktree *node;


  node = SKTreeFindNode (tree, name);

  return (node ? (t_paramtreenode *) node->data : NULL);
}


static int ParameterGetScope (const char *scope)
{
  int retval;


  if (! Util_StrCmpi (scope, "GLOBAL"))
  {
    retval = SCOPE_GLOBAL;
  }
  else if (! Util_StrCmpi(scope, "RESTRICTED"))
  {
    retval = SCOPE_RESTRICTED;
  }
  else if (! Util_StrCmpi(scope, "PRIVATE"))
  {
    retval = SCOPE_PRIVATE;
  }
  else
  {
    retval = -1;
  }

  return (retval);
}


static int ParameterGetType (const char *type)
{
  int retval;


  if (! Util_StrCmpi (type, "KEYWORD"))
  {
    retval = PARAMETER_KEYWORD;
  }
  else if (! Util_StrCmpi (type, "STRING"))
  {
    retval = PARAMETER_STRING;
  }
  else if (! Util_StrCmpi (type, "SENTENCE"))
  {
    retval = PARAMETER_SENTENCE;
  }
  else if (! Util_StrCmpi (type, "INT"))
  {
    retval = PARAMETER_INT;
  }
  else if (! Util_StrCmpi (type, "REAL"))
  {
    retval = PARAMETER_REAL;
  }
  else if (! Util_StrCmpi (type, "BOOLEAN"))
  {
    retval = PARAMETER_BOOLEAN;
  }
  else
  {
    retval = -1;
  }

  return (retval);
}


/*@@
   @routine    ParameterInsert
   @date       Fri Jul 16 10:08:25 1999
   @author     Tom Goodale
   @desc

   @enddesc
@@*/
static int ParameterInsert (t_sktree **tree, t_param *newparam)
{
  int retval;
  t_sktree *treenode;
  t_paramtreenode *node;
  t_paramlist *list;


  treenode = SKTreeFindNode (*tree, newparam->props->name);
  if (treenode)
  {
    retval = ParameterListAddParam (&((t_paramtreenode *) treenode->data)->paramlist, newparam);
  }
  else
  {
    node = (t_paramtreenode *) malloc (sizeof (t_paramtreenode));
    list = (t_paramlist *) malloc (sizeof (t_paramlist));

    if (node && list)
    {
      node->paramlist = list;
      list->param = newparam;
      list->last = list->next = NULL;
      treenode = SKTreeStoreData (*tree, *tree, newparam->props->name, node);
      if (! *tree)
      {
        *tree = treenode;
      }
      retval = 0;
    }
    else
    {
      retval = -1;
      free (list);
      free (node);
    }
  }

  return (retval);
}


/*@@
   @routine    ParameterGetSimple
   @date       Fri Jul 16 10:07:46 1999
   @author     Tom Goodale
   @desc
               Gets the value of a parameter.
   @enddesc
@@*/
static const void *ParameterGetSimple (const t_param *param, int *type)
{
  if (type)
  {
    *type = param->props->type;
  }

  return (param->data);
}


/*@@
   @routine    ParameterExtend
   @date       Thu Jul 15 12:55:06 1999
   @author     Tom Goodale
   @desc
               Adds a range to a parameter.
   @enddesc
@@*/
static int ParameterExtend (t_param *param,
                            const char *range_origin,
                            const char *range,
                            const char *range_description)
{
  int order, retcode;
  t_range *newrange, *rangenode, *lastnode;


  newrange = (t_range *) malloc (sizeof (t_range));
  if (newrange)
  {
    /* Fill out the data */
    newrange->range = strdup (range);
    newrange->origin = strdup (range_origin);
    newrange->description = strdup (range_description);
    newrange->last = newrange->next = NULL;
    newrange->active = 0;

    lastnode = NULL;

    /* Search the list for the right place to insert it. */
    for (rangenode = param->props->range; rangenode;rangenode = rangenode->next)
    {
      lastnode = rangenode;

      order = Util_StrCmpi (range_origin, rangenode->origin);

      if (order <= 0)
      {
        /* Insert before this node */
        newrange->next = rangenode;
        newrange->last = rangenode->last;
        rangenode->last = newrange;
        if (param->props->range == rangenode)
        {
          param->props->range=newrange;
        }
        if (newrange->last)
        {
          newrange->last->next = newrange;
        }
        break;
      }
    }

    if (! rangenode)
    {
      /* Insert at the end of the list */
      newrange->next = NULL;
      newrange->last = lastnode;
      if (param->props->range == NULL)
      {
        param->props->range = newrange;
      }
      if (newrange->last)
      {
        newrange->last->next = newrange;
      }
    }

    retcode = 0;
  }
  else
  {
    retcode = -1;
  }

  return (retcode);
}

 /*@@
   @routine    ParameterSet
   @date       Mon May 20 07:08:03 2002
   @author     Tom Goodale
   @desc
   Sets the value of a parameter
   @enddesc

   @var     param
   @vdesc   The parameter to be set
   @vtype   t_param *
   @vio     in
   @endvar
   @var     value
   @vdesc   The value to set the parameter to
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
    The return code of ParameterSetSimple or
    ParameterSetAccumulator
   @endreturndesc
 @@*/
static int ParameterSet(t_param *param, const char *value)
{
  int retval;

  if(! param->accumulator_bases)
  {
    retval = ParameterSetSimple(param, value);
  }
  else
  {
    char *oldval;

    /* Save old value */
    oldval = CCTK_ParameterValString(param->props->name,
                                     param->props->thorn);

    /* Now try to set this parameter */
    retval = ParameterSetSimple(param, value);

    if(! retval)
    {
      /* OK, that worked. Now try to set accumulator */

      int i;

      /* Loop over possible bases - only one will be active. */
      for(i = 0; i < param->n_accumulator_bases; i++)
      {
        if(CCTK_IsThornActive(param->accumulator_bases[i]->props->thorn))
        {
          retval = ParameterSetAccumulator(param->accumulator_bases[i]);
          break;
        }
      }

      if(retval)
      {
        /* That didn't work, so restore old val. */
        ParameterSetSimple(param, oldval);
      }
    }

    free(oldval);
  }

  return retval;
}


 /*@@
   @routine    ParameterSetAccumulator
   @date       Mon May 20 07:10:53 2002
   @author     Tom Goodale
   @desc
   Sets the value of an accumulator parameter.
   @enddesc

   @var     param
   @vdesc   The parameter to be set
   @vtype   t_param *
   @vio     in
   @endvar

   @returntype int
   @returndesc
     0 - success
    -7 - unsupported accumulator parameter type
    -8 - unsupported extra parameter type
    -9 - final value out of range
   @endreturndesc
 @@*/
static int ParameterSetAccumulator(t_param *param)
{
  int retval;
  uExpression parsed_expression;
  uExpressionValue value;
  uExpressionValue xy[2];

  parsed_expression = Util_ExpressionParse(param->props->accumulator_expression);

  retval = 0;
  switch(param->props->type)
  {
    case PARAMETER_INT :
      xy[0].type = ival;
      xy[0].value.ival = atoi(param->props->defval);
      break;

    case PARAMETER_REAL :
      xy[0].type = rval;
      xy[0].value.rval = atof(param->props->defval);
      break;

    default :
      retval = -7;
  }

  if(! retval)
  {
    int i;
    /* Assume no real or int is larger than 100 characters */
    char newval[100];
    for(i = 0; i < param->n_accumulator_sources; i++)
    {
      if(CCTK_IsThornActive(param->accumulates_from[i]->props->thorn))
      {
        switch(param->accumulates_from[i]->props->type)
        {
          case PARAMETER_INT :
            xy[1].type = ival;
            xy[1].value.ival = *((CCTK_INT *)param->accumulates_from[i]->data);
            break;

          case PARAMETER_REAL :
            xy[1].type = rval;
            xy[1].value.rval = *((CCTK_REAL *)param->accumulates_from[i]->data);;
            break;

          default :
            retval = -8;
        }

        retval = Util_ExpressionEvaluate(parsed_expression,
                                         &value,
                                         AccVarEvaluator,
                                         xy);

        xy[0] = value;
      }
    }

    switch(value.type)
    {
      case ival:
        sprintf(newval,"%ld",(long)value.value.ival);
        break;
      case rval:
        sprintf(newval,"%.20g",(double)value.value.rval);
        break;
      default :
        ;
    }

    retval = ParameterSetSimple (param, newval);

    if(retval)
    {
      retval = -9;
    }

  }

  if (parsed_expression)
  {
    Util_ExpressionFree (parsed_expression);
  }

  return retval;
}

 /*@@
   @routine    ParameterSetSimple
   @date       Mon May 20 07:08:03 2002
   @author     Tom Goodale
   @desc
   Sets the value of a parameter
   @enddesc

   @var     param
   @vdesc   The parameter to be set
   @vtype   t_param *
   @vio     in
   @endvar
   @var     value
   @vdesc   The value to set the parameter to
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
    The return code of the basic setting routine, or -2 if
    the type is unknown.
   @endreturndesc
 @@*/
static int ParameterSetSimple (t_param *param, const char *value)
{
  int retval;


  switch (param->props->type)
  {
    case PARAMETER_KEYWORD:
      retval = ParameterSetKeyword (param, value); break;
    case PARAMETER_STRING:
      retval = ParameterSetString (param, value); break;
    case PARAMETER_SENTENCE:
      retval = ParameterSetSentence (param, value); break;
    case PARAMETER_INT:
      retval = ParameterSetInteger (param, value); break;
    case PARAMETER_REAL:
      retval = ParameterSetReal (param, value); break;
    case PARAMETER_BOOLEAN:
      retval = ParameterSetBoolean (param, value); break;
    default:
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Unknown parameter type %d", param->props->type);
      retval = -2;
  }

  return (retval);
}


static int ParameterListAddParam (t_paramlist **paramlist, t_param *newparam)
{
  t_paramlist *node;


  node = (t_paramlist *) malloc (sizeof (t_paramlist));
  if (node)
  {
    node->param = newparam;

    /* Place at beginning of list for now. */
    node->next = *paramlist;
    node->last = NULL;
    (*paramlist)->last = node;

    *paramlist = node;
  }

  return (node != NULL);
}


static int ParameterSetKeyword (t_param *param, const char *value)
{
  int retval;
  const t_range *range;


  retval = -1;
  for (range = param->props->range; range; range = range->next)
  {
    if (CCTK_IsThornActive (range->origin) ||
        CCTK_Equals (param->props->thorn, range->origin))
    {
      if (!Util_StrCmpi(value, range->range))
      {
        retval = CCTK_SetString (param->data, value);
        break;
      }
    }
  }

  if (retval == -1)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                "ParameterSetKeyword: Unable to set keyword '%s::%s', "
                "new value '%s' is not in any active range",
                param->props->thorn, param->props->name, value);

    if (*(char **) param->data == NULL)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Since this was the default value, "
                  "setting anyway - please fix!");
      CCTK_SetString (param->data, value);
    }
  }

  return (retval);
}


static int ParameterSetString (t_param *param, const char *value)
{
  int retval;
  const t_range *range;


  retval = -1;
  for (range = param->props->range; range; range = range->next)
  {
    if (CCTK_IsThornActive (range->origin) ||
        CCTK_Equals (param->props->thorn, range->origin))
    {
#ifndef CCTK_PARAMUNCHECKED
      const int matched = CCTK_RegexMatch (value, range->range, 0, NULL);
      if (matched > 0)
      {
#endif
        retval = CCTK_SetString (param->data, value);
        break;
#ifndef CCTK_PARAMUNCHECKED
      }
      else if (matched < 0)
      {
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                   "Invalid regular expression '%s' used as range for string %s::%s",
                   range->range, param->props->thorn, param->props->name);
      }
#endif
    }
  }

  if (retval == -1)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                "ParameterSetString: Unable to set string '%s::%s', "
                "new value '%s' is not in any active range",
                param->props->thorn, param->props->name, value);

    if (*(char **) param->data == NULL)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Since this was the default value, "
                  "setting anyway - please fix!");
      CCTK_SetString (param->data, value);
    }
  }

  return (retval);
}


static int ParameterSetSentence (t_param *param, const char *value)
{
  int retval;
  const t_range *range;


  retval = -1;
  for (range = param->props->range; range; range = range->next)
  {
    if (CCTK_IsThornActive (range->origin) ||
        CCTK_Equals (param->props->thorn, range->origin))
    {
#ifndef CCTK_PARAMUNCHECKED
      const int matched = CCTK_RegexMatch (value, range->range, 0, NULL);
      if (matched > 0)
      {
#endif
        retval = CCTK_SetString (param->data, value);
        break;
#ifndef CCTK_PARAMUNCHECKED
      }
      else if (matched < 0)
      {
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                   "Invalid regular expression '%s' used as range for sequence %s::%s",
                   range->range, param->props->thorn, param->props->name);
      }
#endif
    }
  }

  if (retval == -1)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                "ParameterSetSentence: Unable to set sentence '%s::%s' = '%s' "
                "not in any active range",
                param->props->thorn, param->props->name, value);

    if (*(char **) param->data == NULL)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Since this was the default value, "
                  "setting anyway - please fix!");
      CCTK_SetString (param->data, value);
    }
  }

  return (retval);
}


static int ParameterSetInteger (t_param *param, const char *value)
{
  int inval, retval;
  const t_range *range;
  char *endptr;

  retval = 0;

  /* try parsing as number */
  if (strcmp(value, "") == 0) /* strtol returns 0 in this case */
    retval = -6;
  else
    inval = strtol (value, &endptr, 0);

  if (!retval && *endptr) /* if we cannot parse as a number, try expression */
  {
    int type = PARAMETER_INT;
    uExpressionValue val;
    uExpression *expr;

    expr = Util_ExpressionParse(value);
    if (expr != NULL)
    {
      retval = Util_ExpressionEvaluate(expr, &val, SetVarEvaluator, &type);
      Util_ExpressionFree(expr);
    }
    else
    {
      retval = -6;
    }

    if (retval == 0)
    {
      assert(val.type == ival || val.type == rval);

      if(val.type == ival)
      {
        inval =(int)val.value.ival;
      }
      else if(fabs(round(val.value.rval) - val.value.rval) < 1e-12) /* enforce integer result */
      {
        inval = (int)lrint(val.value.rval);
      }
      else
      {
        retval = -6;
      }
    }
    else
    {
      retval = -6;
    }
  }

  if (!retval)
  {
    retval = -1;
    for (range = param->props->range; range; range = range->next)
    {
      if (CCTK_IsThornActive (range->origin) ||
          CCTK_Equals (param->props->thorn, range->origin))
      {
#ifndef CCTK_PARAMUNCHECKED
        if (Util_IntInRange (inval, range->range))
        {
#endif
          *(CCTK_INT *) param->data = inval;
          retval = 0;
          break;
#ifndef CCTK_PARAMUNCHECKED
        }
#endif
      }
    }
  }

  if (retval == -1)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                "ParameterSetInteger: Unable to set integer '%s::%s' = '%s' "
                "not in any active range",
                param->props->thorn, param->props->name, value);
  }
  else if (retval == -6)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                "ParameterSetInteger: Unable to set integer '%s::%s' - '%s' "
                "does not evaluate to a valid integer",
                param->props->thorn, param->props->name, value);
  }

  return (retval);
}


static int ParameterSetReal (t_param *param, const char *value)
{
  int retval;
  const t_range *range;
  double inval;
  char *temp;
  char *endptr;

  retval = 0;

  /*
   * Canonicalize the string by converting exponent letters
   * (we allow [eEdD]) to 'e', since strtod(3) (which we will use
   * to do the actual conversion) only groks [eE].
   */
  temp = strdup (value);
  for (char *p = temp; *p; ++p)
  {
    if (*p == 'd' || *p == 'D')
    {
      *p = 'e';
      break;
    }
  }

  /* try parsing as number */
  if (strcmp (value, "") == 0) /* strtod returns 0. in this case */
    retval = -6;
  else
    inval = strtod (temp, &endptr);

  if (!retval && *endptr) /* if we cannot parse as a number, try expression */
  {
    int type = PARAMETER_REAL;
    uExpressionValue val;
    uExpression *expr;

    expr = Util_ExpressionParse(value);
    if (expr != NULL)
    {
      retval = Util_ExpressionEvaluate(expr, &val, SetVarEvaluator, &type);
      Util_ExpressionFree(expr);
    }
    else
    {
      retval = -6;
    }

    if (retval == 0)
    {
      assert(val.type == ival || val.type == rval);

      if (val.type == ival)
        inval = (int)val.value.ival;
      else
        inval = val.value.rval;
    }
    else
    {
      retval = -6;
    }
  }

  free(temp); /* must be after last access to *endptr */

  if (!retval)
  {
    retval = -1;
    for (range = param->props->range; range; range = range->next)
    {
      if (CCTK_IsThornActive (range->origin) ||
          CCTK_Equals (param->props->thorn, range->origin))
      {
#ifndef CCTK_PARAMUNCHECKED
        if(Util_DoubleInRange (inval, range->range))
        {
#endif
          *(CCTK_REAL *) param->data = inval;
          retval = 0;
          break;
#ifndef CCTK_PARAMUNCHECKED
        }
#endif
      }
    }
  }

  if (retval == -1)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                "ParameterSetReal: Unable to set real '%s::%s' = '%s' not in "
                "any active range",
                param->props->thorn, param->props->name, value);
  }
  else if (retval == -6)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                "ParameterSetInteger: Unable to set real '%s::%s' - '%s' "
                "does not evaluate to a valid floating point number",
                param->props->thorn, param->props->name, value);
  }

  return (retval);
}


static int ParameterSetBoolean (t_param *param, const char *value)
{
  int type = PARAMETER_BOOLEAN;
  int retval;
  CCTK_INT inval;
  uExpressionValue val;
  uExpression *expr;

  /* first try parsing as yes/no/true/false */
  retval = CCTK_SetBoolean (&inval, value);

  if (retval < 0) /* if we cannot parse as a boolean, try expression */ 
  {
    expr = Util_ExpressionParse(value);
    if(expr != NULL)
    {
      retval = Util_ExpressionEvaluate(expr, &val, SetVarEvaluator, &type);
      Util_ExpressionFree(expr);
    }
    else
    {
      retval = -6;
    }

    if (retval == 0)
    {
      assert(val.type == ival || val.type == rval);
      if (val.type == ival)
      {
        inval =(int)val.value.ival;
        retval = 0;
      }
      else if(fabs(round(val.value.rval) - val.value.rval) < 1e-12) /* enforce integer result */
      {
        inval = (int)lrint(val.value.rval);
        retval = 0;
      }
      else
      {
        retval = -1;
      }
    }
    else
    {
      retval = -1;
    }
  }

  if (!retval)
    *(CCTK_INT *)param->data = inval != 0;
  else
  {
    CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                "ParameterSetBoolean: Unable to set boolean '%s::%s' = '%s' "
                "not recognised as boolean",
                param->props->thorn, param->props->name, value);
  }

  return (retval);
}


void CCTKi_SetParameterSetMask (int mask)
{
  cctk_parameter_set_mask = mask;
}

 /*@@
   @routine    GetBaseName
   @date       Sun May 19 19:08:28 2002
   @author     Tom Goodale
   @desc
   Gets the basename of a (possibly) array parameter.
   @enddesc

   @var     name
   @vdesc   Parameter name
   @vtype   const char *
   @vio     in
   @endvar
   @var     basename
   @vdesc   Name of array base parameter
   @vtype   char **
   @vio     out
   @endvar
   @var     array_index
   @vdesc   Array index
   @vtype   int
   @vio     out
   @endvar
 @@*/
static void GetBaseName(const char *name, char **basename, int *array_index)
{
  int baselen;
  const char *pos;

  pos = strchr(name, '[');

  if(pos)
  {
    baselen = pos-name;
    *array_index  = atoi(pos+1);
  }
  else
  {
    baselen = strlen(name);
    *array_index = -1;
  }

  *basename = (char *)malloc(baselen+1);

  if(basename)
  {
    strncpy(*basename, name, baselen);
    (*basename)[baselen] = 0;
  }
}
 /*@@
   @routine    ArrayParamName
   @date       Sun May 19 22:03:44 2002
   @author     Tom Goodale
   @desc
   Takes the basename of an array parameter, and an array index,
   and returns the name in the form basename[index].
   @enddesc

   @var     basename
   @vdesc   Basename of an array parameter
   @vtype   const char *
   @vio     in
   @endvar
   @var     array_index
   @vdesc   index of parameter into array
   @vtype   int
   @vio     in
   @endvar

   @returntype char *
   @returndesc
    string contining the name of the array parameter.
   @endreturndesc
 @@*/

static char *ArrayParamName(const char *basename,int array_index)
{
  char *retval;

  /* Assume the string representation of an integer is no greater than 40 chars */
  retval = (char *)malloc(strlen(basename)+2+40+1);

  if(retval)
  {
    sprintf(retval, "%s[%d]",basename,array_index);
  }

  return retval;
}

 /*@@
   @routine    AccVarEvaluator
   @date       Mon May 20 07:04:03 2002
   @author     Tom Goodale
   @desc
   Routine called from the expression parser to evaluate
   the vars in an accumulator expression
   @enddesc

   @var     nvars
   @vdesc   Number of variables to evaluate
   @vtype   int
   @vio     in
   @vcomment

   @endvar
   @var     vars
   @vdesc   an array of variable names
   @vtype   const char * const *
   @vio     in
   @vcomment
   Should be just x and y.
   @endvar
   @var     vals
   @vdesc   Output array to hold values
   @vtype   uExpressionValue *
   @vio     out
   @vcomment

   @endvar
   @var     data
   @vdesc   Data passed from expression evaluator
   @vtype   const void *
   @vio     in
   @vcomment
     Should be an array of two uExpressionValues,
     one for x, one for y.
   @endvar

   @returntype int
   @returndesc
   0
   @endreturndesc
 @@*/
static int AccVarEvaluator(int nvars, const char * const *vars, uExpressionValue *vals, const void *data)
{
  int i;
  const uExpressionValue *exps = (const uExpressionValue *)data;

  for(i=0; i < nvars; i++)
  {
    if(strcmp(vars[i], "x"))
    {
      vals[i] = exps[0];
    }
    else if(strcmp(vars[i], "y"))
    {
      vals[i] = exps[1];
    }
    else
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "AccVarEvaluator: unrecognised '%s' in expression",
                  vars[i]);
    }
  }

  return 0;
}

 /*@@
   @routine    SetVarEvaluator
   @date       Wed Oct 26 16:25:49 PDT 2011
   @author     Roland Haas
   @desc
   Routine called from the expression parser to evaluate
   the vars in an parameter expression
   @enddesc

   @var     nvars
   @vdesc   Number of variables to evaluate
   @vtype   int
   @vio     in
   @vcomment

   @endvar
   @var     vars
   @vdesc   an array of variable names
   @vtype   const char * const *
   @vio     in
   @vcomment
   Can be anything that GetParamter will recognize. Must be of type CCTK_INT,
   CCTK_BOOLEAN or CCTK_REAL.
   @endvar
   @var     vals
   @vdesc   Output array to hold values
   @vtype   uExpressionValue *
   @vio     out
   @vcomment

   @endvar
   @var     data
   @vdesc   Data passed from expression evaluator
   @vtype   const void *
   @vio     in
   @vcomment
   Not used.
   @endvar

   @returntype int
   @returndesc
   0
   @endreturndesc
 @@*/
static int SetVarEvaluator(int nvars, const char * const *vars, uExpressionValue *vals, const void *data)
{
  const int restype = *(const int *)data;
  int retval = 0;

  for (int i=0; i < nvars; i++)
  {
    int ierr;

    if (strstr(vars[i], "::")) /* a variable [parameter] */
    {
      const void *paramval;
      int type;
      char *name, *thorn;

      ierr = Util_SplitString(&thorn, &name, vars[i], "::");
      if (!ierr)
      {
        paramval = CCTK_ParameterGet(name, thorn, &type);
        if (paramval)
        {
          switch(type)
          {
            case PARAMETER_REAL:
              vals[i].type = rval;
              vals[i].value.rval = *(const CCTK_REAL *)paramval;
              ierr = 0;
              break;
            case PARAMETER_INT:
              vals[i].type = ival;
              vals[i].value.ival = *(const CCTK_INT *)paramval;
              ierr = 0;
              break;
            case PARAMETER_BOOLEAN:
              vals[i].type = ival;
              vals[i].value.ival = *(const CCTK_INT *)paramval;
              ierr = 0;
              break;
            default:
              CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                          "SetVarEvaluator: cannot handle type %d for parameter '%s::%s'. Only REAL, INT and BOOLEAN are supported.",
                          type, thorn, name);
              ierr = -1;
              break;
          }
        }
        else
        {
          CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                      "SetVarEvaluator: could not find '%s::%s'",
                      thorn, name);
          ierr = -1;
        }

        free(thorn);
        free(name);
      }
      else
      {
        CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                    "SetVarEvaluator: cannot split '%s' into thorn::parameter: %d",
                    vars[i], ierr);
        ierr = -1;
      }
    }
    else /* a direct value */
    {
      char *endptr, *temp;
      if(restype == PARAMETER_BOOLEAN && CCTK_SetBoolean(&vals[i].value.ival, vars[i]) == 0)
      {
          vals[i].type = ival;
          ierr = 0;
      }
      else if (strpbrk(vars[i], "eDdD.") ||
               Util_StrCmpi(vars[i], "nan") == 0 ||
               Util_StrCmpi(vars[i], "inf") == 0)
      {
        /*
         * Canonicalize the string by converting all exponent letters
         * (we allow [eEdD]) to 'e', since strtod(3) (which we will use
         * to do the actual conversion) only groks [eE].
         */
        temp = strdup (vars[i]);
        for (unsigned int p = 0; p < strlen (temp); p++)
        {
          if (temp[p] == 'E' || temp[p] == 'd' || temp[p] == 'D')
          {
            temp[p] = 'e';
            break;
          }
        }

        vals[i].type = rval;
        vals[i].value.rval = (CCTK_REAL)strtod (temp, &endptr);
        if (!*endptr)
          ierr = 0;
        else
          ierr = -1;

        free(temp);
      }
      else
      {
          vals[i].type = ival;
          vals[i].value.ival = (CCTK_INT)strtol (vars[i], &endptr, 0);
          if (!*endptr)
            ierr = 0;
          else
            ierr = -1;
      }
    }

    if (ierr)
    {
      CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                  "SetVarEvaluator: Unable to set value - '%s' "
                  "does not evaluate to a real or integer or boolean",
                  vars[i]);
      vals[i].type = rval;
      vals[i].value.rval = (double)atof("nan");
      retval = ierr;
    }
  }

  return retval;
}

 /*@@
   @routine    AddAccumulators
   @date       Mon May 20 07:27:09 2002
   @author     Tom Goodale
   @desc
   Adds accumulator data to parameters
   @enddesc

   @var     base
   @vdesc   accumulator base parameter
   @vtype   t_param *
   @vio     in
   @endvar
   @var     extra
   @vdesc   accumulator source parameter
   @vtype   t_param *
   @vio     in
   @endvar
   @var     thorn
   @vdesc   thorn providing source parameter
   @vtype   const char *
   @vio     in
   @endvar
   @var     parameter
   @vdesc   source parameter
   @vtype   const char *
   @vio     in
   @endvar
   @var     imp
   @vdesc   implementationm or thorn providing base parameter
   @vtype   const char *
   @vio     in
   @endvar
   @var     baseparam
   @vdesc   base parameter
   @vtype   const char *
   @vio     in
   @endvar

 @@*/
static void AddAccumulators(t_param *base,
                            t_param *extra,
                            const char *thorn,
                            const char *parameter,
                            const char *imp,
                            const char *baseparam)
{
  int errcode;

  errcode = LinkAccumulators(base, extra);

  if(! errcode)
  {
    /* If the base parameter is an array parameter, copy data to its members. */
    if(base->array && extra->array &&
       base->props->array_size == extra->props->array_size)
    {
      int i;

      for(i=0; i < base->props->array_size; i++)
      {
        if(LinkAccumulators(&(base->array[i]),&(extra->array[i])))
        {
          CCTK_Warn (0, __LINE__, __FILE__, "Cactus",
                     "CCTKi_ParameterAccumulatorBase: error, probably out of memory");
        }
      }
    }
    else if(base->array || extra->array)
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "Accumulator base parameter %s::%s and parameter %s::%s have different array sizes",
                  imp,baseparam,thorn,parameter);
    }
  }
  else
  {
    CCTK_Warn (0, __LINE__, __FILE__, "Cactus",
               "CCTKi_ParameterAccumulatorBase: error, probably out of memory");
  }
}

 /*@@
   @routine    LinkAccumulators
   @date       Mon May 20 07:29:48 2002
   @author     Tom Goodale
   @desc
   Links the accumulator data on two parameters.
   @enddesc

   @var     base
   @vdesc   accumulator base parameter
   @vtype   t_param *
   @vio     in
   @endvar
   @var     extra
   @vdesc   accumulator source parameter
   @vtype   t_param *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0  - success
   -1 - out of memory
   -2 - out of memory
   @endreturndesc
 @@*/
static int LinkAccumulators(t_param *base, t_param *extra)
{
  int retcode;

  t_param **temp;

  temp = (t_param **)realloc(base->accumulates_from,
                             (base->n_accumulator_sources+1)*sizeof(t_param *));

  if(temp)
  {
    /* Update the base parameter */
    base->accumulates_from = temp;
    base->accumulates_from[base->n_accumulator_sources++] = extra;

    /* Update the extra parameter */
    temp = (t_param **)realloc(extra->accumulator_bases,
                               (extra->n_accumulator_bases+1)*sizeof(t_param *));


    if(temp)
    {
      /* Update the base parameter */
      extra->accumulator_bases = temp;
      extra->accumulator_bases[extra->n_accumulator_bases++] = base;

      retcode = 0;
    }
    else
    {
      retcode = -2;
    }
  }
  else
  {
    retcode = -1;
  }

  return retcode;
}


 /*@@
   @routine    ParameterActivate
   @date       Mon May 20 06:59:51 2002
   @author     Tom Goodale
   @desc
   Does any necessary activations on a parameter
   @enddesc

   @var     param
   @vdesc   parameter to be activated
   @vtype   t_param *
   @vio     in
   @endvar

 @@*/
static void ParameterActivate(t_param *param)
{
  if(param->array)
  {
    int i;

    for(i=0; i < param->props->array_size; i++)
    {
      ParameterActivate(&(param->array[i]));
    }

  }
  else if(param->accumulator_bases)
  {
    ParameterSet(param, param->props->defval);
  }
}


 /*@@
   @routine CCTK_ParameterSetNotifyRegister
   @date    19 September 2006
   @author  Thomas Radke
   @desc
            Registers a parameter set notify callback
   @enddesc 
   
   @var     callback
   @vdesc   user-supplied callback function
   @vtype   cParameterSetNotifyCallbackFn
   @vio     in
   @endvar
   @var     data
   @vdesc   pointer to an optional user-supplied data structure
   @vtype   void *
   @vio     in
   @endvar
   @var     name
   @vdesc   unique name of the notify callback to register
   @vtype   const char *
   @vio     in
   @endvar
   @var     thorn_regex
   @vdesc   optional regular expression for the thorn name to match
   @vtype   const char *
   @vio     in
   @endvar
   @var     param_regex
   @vdesc   optional regular expression for the parameter name to match
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
                0 - success
               -1 - another callback has already been registered under this name
               -2 - out of memory
               -3 - invalid regular expression given for thorn_regex/param_regex
   @endreturndesc
 @@*/
int CCTK_ParameterSetNotifyRegister (cParameterSetNotifyCallbackFn callback,
                                     void *data,
                                     const char *name,
                                     const char *thorn_regex,
                                     const char *param_regex)
{
  int handle;
  t_notifyCallback *cb;


  handle = Util_GetHandle (notifyCallbackList, name, NULL);
  if (handle >= 0)
  {
    return (-1);
  }

  cb = (t_notifyCallback *) malloc (sizeof (t_notifyCallback));
  if (cb == NULL)
  {
    return (-2);
  }

  cb->fn = callback;
  cb->data = data;
  cb->has_thorn_regex = thorn_regex && *thorn_regex;
  cb->has_param_regex = param_regex && *param_regex;
  if ((cb->has_thorn_regex && regcomp (&cb->thorn_preg, thorn_regex,
                                       REG_EXTENDED | REG_ICASE | REG_NOSUB)) ||
      (cb->has_param_regex && regcomp (&cb->param_preg, param_regex, 
                                       REG_EXTENDED | REG_ICASE | REG_NOSUB)))
  {
    return (-3);
  }

  handle = Util_NewHandle (&notifyCallbackList, name, cb);
  num_notify_callbacks++;

  return (0);
}

void CCTK_FCALL CCTK_FNAME (CCTK_ParameterSetNotifyRegister)
                           (int *status, cParameterSetNotifyCallbackFn callback,
                            void *data, THREE_FORTSTRING_ARG)
{
  THREE_FORTSTRING_CREATE (name, thorn_regex, param_regex)
  *status = CCTK_ParameterSetNotifyRegister (callback, data, name,
                                             thorn_regex, param_regex);
  free (param_regex);
  free (thorn_regex);
  free (name);
}


 /*@@
   @routine CCTK_ParameterSetNotifyUnregister
   @date    19 September 2006
   @author  Thomas Radke
   @desc
            Unregisters a parameter set notify callback
   @enddesc

   @var     name
   @vdesc   unique name of the notify callback to unregister
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
                0 - success
               -1 - no callback was registered under this name
   @endreturndesc
 @@*/
int CCTK_ParameterSetNotifyUnregister (const char *name)
{
  int handle;
  t_notifyCallback *cb;


  handle = Util_GetHandle (notifyCallbackList, name, (void **) &cb);
  if (handle >= 0)
  {
    Util_DeleteHandle (notifyCallbackList, handle);
    if (cb->has_thorn_regex)
    {
      regfree (&cb->thorn_preg);
    }
    if (cb->has_param_regex)
    {
      regfree (&cb->param_preg);
    }
    free (cb);
  }

  return (handle >= 0 ? 0 : -1);
}

void CCTK_FCALL CCTK_FNAME (CCTK_ParameterSetNotifyUnregister)
                           (int *status, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (name)
  *status = CCTK_ParameterSetNotifyUnregister (name);
  free (name);
}


/*****************************************************************************/

/*#define TEST_PARAMETERS*/
#ifdef TEST_PARAMETERS

/*@@
   @routine    ParameterPrintDescription
   @date       Tue Jun 29 10:24:49 1999
   @author     Tom Goodale
   @desc
   Prints out a description on the given file descriptor with the given
   format.  Should include all data - i.e. ranges, range descriptions,
   range origins, and default value.

   @enddesc

   @var     name
   @vdesc   The name of the parameter
   @vtype   const char *
   @vio     in
   @vcomment
   name
   @endvar
   @var     thorn
   @vdesc   The originating thorn
   @vtype   const char *
   @vio     in
   @vcomment

   @endvar
   @var     implementation
   @vdesc   The originating implementation
   @vtype   const char *
   @vio     in
   @vcomment

   @endvar
   @var     format
   @vdesc   The printf format string
   @vtype   const char *
   @vio     in
   @vcomment
   This is for each line of the test, and should have one %s in it.
   @endvar
   @var     file
   @vdesc   File descriptor
   @vtype   FILE *
   @vio     in
   @vcomment
   The file to pront out on.
   @endvar

   @returntype int
   @returndesc
   1 on success, 0 on failure.
   @endreturndesc

@@*/
int ParameterPrintDescription(const char *name,
                              const char *thorn,
                              const char *format,
                              FILE *file)
{
  int retval;
  t_param *param;

  param = NULL;

  param = ParameterFind(name, thorn, SCOPE_ANY);

  if(param)
  {
    retval = ParameterPrintSimple(param, format, file);
  }
  else
  {
    retval = -1;
  }

  return retval;
}


static int ParameterPrintSimple(t_param *param,
                                const char *format,
                                FILE *file)
{
  t_range *range;

  fprintf(file, format, "Parameter", param->props->name);
  fprintf(file, format, "Thorn", param->props->thorn);
  fprintf(file, format, "Desc", param->props->description);
  fprintf(file, format, "Def", param->props->defval);

  for(range=param->props->range; range; range=range->next)
  {
    fprintf(file, format, "Range:", range->range);
  }

  return 0;
}
struct
{
  int a;
  char *foo;
  double beta;
} params;

int main(void)
{
  CCTKi_ParameterCreate("a", "thorn1", "imp1", "int", "global", 0,
                        "The a param", "2", &(params.a));

  CCTKi_ParameterCreate("foo", "thorn2", "imp2", "keyword", "private", 0,
                        "The foo param", "bingo", &(params.foo));

  printf("Testing thorn,null\n");

  ParameterPrintDescription("a",
                            "thorn1", /*const char *thorn,*/
                            NULL, /* const char *implementation,*/
                            "..%s..%s\n",/*  const char *format,*/
                            stdout);

  printf("Testing null,imp\n");

  ParameterPrintDescription("a",
                            NULL, /*const char *thorn,*/
                            "imp1", /* const char *implementation,*/
                            "..%s..%s\n",/*  const char *format,*/
                            stdout);

  printf("Testing thorn,thorn\n");

  ParameterPrintDescription("a",
                            "thorn1", /*const char *thorn,*/
                            "thorn1", /* const char *implementation,*/
                            "..%s..%s\n",/*  const char *format,*/
                            stdout);

  printf("Testing imp,imp\n");

  ParameterPrintDescription("a",
                            "imp1", /*const char *thorn,*/
                            "imp1", /* const char *implementation,*/
                            "..%s..%s\n",/*  const char *format,*/
                            stdout);

  printf("Testing imp,null\n");
  ParameterPrintDescription("a",
                            "imp1", /*const char *thorn,*/
                            NULL, /* const char *implementation,*/
                            "..%s..%s\n",/*  const char *format,*/
                            stdout);

  printf("Adding a range to a\n");
  ParameterAddRange("imp1",
                    "a",
                    "imp1",
                    "1:7:0",
                    "A nice range for a");

  printf("Adding another range to a\n");
  ParameterAddRange("imp1",
                    "a",
                    "imp1",
                    "1:7:1",
                    "Another nice range for a");

  printf("a is now\n");

  ParameterPrintDescription("a",
                            "thorn1", /*const char *thorn,*/
                            NULL, /* const char *implementation,*/
                            "..%s..%s\n",/*  const char *format,*/
                            stdout);

  return 0;
}

#endif
