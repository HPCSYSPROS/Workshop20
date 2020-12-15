/*@@
   @file      Groups.c
   @date      Mon Feb  1 12:16:28 1999
   @author    Tom Goodale
   @desc
              Routines to deal with groups.
   @enddesc
   @version   $Id$
 @@*/

#include <assert.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "cctk_Constants.h"
#include "cctk_WarnLevel.h"
#include "cctk_Flesh.h"
#include "cctk_FortranString.h"
#include "cctk_Groups.h"
#include "cctk_Parameter.h"
#include "cctk_Types.h"
#include "cctk_ActiveThorns.h"

#include "cctki_Groups.h"

#include "util_Expression.h"
#include "util_Table.h"

#include "util_String.h"

/*#define DEBUG_GROUPS*/

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_Groups_c);


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/
/* prototypes for external C routines are declared in header cctk_Groups.h
   here only follow the fortran wrapper prototypes */
void CCTK_FCALL CCTK_FNAME (CCTK_FirstVarIndex)
                           (int *first, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_FirstVarIndexI)
                           (int *first, const int *group);
void CCTK_FCALL CCTK_FNAME (CCTK_FullName)
                           (int *nchars, int *var, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupDimI)
                           (int *dim, const int *group);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupDimFromVarI)
                           (int *dim, const int *vi);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupDistribNumber)
                           (int *number, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupIndex)
                           (int *vindex, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupIndexFromVar)
                           (int *vindex, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupIndexFromVarI)
                           (int *gindex, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupName)
                           (int *nchars, int *var, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupNameFromVarI)
                           (int *nchars, int *var, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupScopeNumber)
                           (int *number, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupTagsTable)
                           (int *table, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupTagsTableI)
                           (int *table, int *group);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupTypeFromVarI)
                           (int *type, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupTypeI)
                           (int *type, const int *group);
void CCTK_FCALL CCTK_FNAME (CCTK_GroupTypeNumber)
                           (int *number, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_ImpFromVarI)
                           (int *nchars, int *var, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_MaxDim)
                           (int *dim);
void CCTK_FCALL CCTK_FNAME (CCTK_MaxGFDim)
                           (int *dim);
void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevels)
                           (int *num, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevelsGI)
                           (int *num, const int *group);
void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevelsGN)
                           (int *num, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevelsVI)
                           (int *num, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevelsVN)
                           (int *num, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_NumGroups)
                           (int *num_groups);
void CCTK_FCALL CCTK_FNAME (CCTK_NumVars)
                           (int *num_vars);
void CCTK_FCALL CCTK_FNAME (CCTK_NumVarsInGroup)
                           (int *num, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_NumVarsInGroupI)
                           (int *numvars, int *group);
void CCTK_FCALL CCTK_FNAME (CCTK_PrintGroup)
                           (const int *group);
void CCTK_FCALL CCTK_FNAME (CCTK_PrintVar)
                           (const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_VarIndex)
                           (int *vindex, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_VarName)
                           (int *nchars, int *var, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_VarTypeI)
                           (int *type, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_VarTypeName)
                           (int *nchars, int *var, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_VarTypeNumber)
                           (int *number, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_VarTypeSize)
                           (int *size, int *vtype);

/* DEPRECATED IN BETA 13 */
void CCTK_FCALL CCTK_FNAME (CCTK_NumTimeLevelsFromVarI)
                           (int *num, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_NumTimeLevelsFromVar)
                           (int *num, ONE_FORTSTRING_ARG);
void CCTK_FCALL CCTK_FNAME (CCTK_NumTimeLevelsI)
                           (int *num, const int *var);
void CCTK_FCALL CCTK_FNAME (CCTK_NumTimeLevels)
                           (int *num, ONE_FORTSTRING_ARG);

/* prototype for CCTKi_VarDataPtr() doesn't appear in a header file because it
   is only used in the variable bindings (see grdoc for CCTKi_VarDataPtr()) */
void *CCTKi_VarDataPtr(const cGH *GH, int timelevel,
                       const char *implementation, const char *varname);

/********************************************************************
 ********************    Internal Typedefs   ************************
 ********************************************************************/
typedef struct
{
  char *name;
  int number;

  /* dimensional_comm_array[dim] */
  char *dimensional_comm_array;
} cVariableDefinition;

typedef struct
{
  /* The various names of the thing. */
  char *thorn,
       *implementation,
       *name;

  /* The group number. */
  int number;

  /* The types. */
  int gtype,
      vtype,
      dtype;

  int gscope;

  int dim;

  int n_timelevels;

  int n_variables;

  /* 1 for scalar groups, number of vector elements for vector groups */
  int vectorlength;

  /* *size[dim]  - pointers to parameter data*/
  CCTK_INT **size;

  /* *ghostsize[dim]  - pointers to parameter data*/
  CCTK_INT **ghostsize;

  /* variables[n_variables] */
  cVariableDefinition *variables;

  /* Variable array size parameter */
  const char *vararraysize;

  /* Stuff for group tags.  Really want table on a per-cGH basis,
   * but that means we need a cGH to get data which would rule
   * out Startup until 4.1
   */
  const char *tags_string;
  int tags_table;

} cGroupDefinition;


/********************************************************************
 ********************    Static Variables   *************************
 ********************************************************************/
/* Static variables needed to hold group and variable data. */

static int n_groups = 0;
static cGroupDefinition *groups = NULL;

static int total_variables = 0;

static int *group_of_variable = NULL;

static int maxdim = 0;
static int gfdim = 0;


/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/
static cGroupDefinition *CCTKi_SetupGroup (const char *implementation,
                                           const char *name,
                                           int n_variables,
                                           int vectorlength);
static CCTK_INT **CCTKi_ExtractSize (int dimension,
                                     const char *thorn,
                                     const char *sizestring,
                                     const char *gname);

static int CCTKi_ParamExpressionToInt(const char *expression, const char *thorn);


 /*@@
   @routine    CCTK_GroupIndex
   @date       Fri Jan 29 08:43:48 1999
   @author     Tom Goodale
   @desc
               Gets the global index number for the specified group.
   @enddesc

   @returntype int
   @returndesc
               the index of the given group, or negative for failure:
               -1 if no group of such name exists, or
                  error return code of @seeroutine CCTK_DecomposeName
   @endreturndesc
@@*/
int CCTK_GroupIndex (const char *fullgroupname)
{
  int group;
  int retval;
  char *impname,
       *groupname;


  impname = groupname = NULL;
  retval = CCTK_DecomposeName (fullgroupname, &impname, &groupname);
  if (! retval)
  {
    retval = -1;
    for (group = 0; group < n_groups; group++)
    {
      if (CCTK_Equals (impname, groups[group].implementation) &&
          CCTK_Equals (groupname, groups[group].name))
      {
        retval = group;
        break;
      }
    }

    if (retval < 0)
    {
      CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                  "CCTK_GroupIndex: No group named '%s' found",
                  fullgroupname);
    }
  }

  /* Free memory from CCTK_DecomposeName */
  free (impname);
  free (groupname);

  return (retval);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupIndex)
                           (int *vindex, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (name)
  *vindex = CCTK_GroupIndex (name);
  free (name);
}


 /*@@
   @routine    CCTK_VarIndex
   @date       Mon Feb  8 12:03:22 1999
   @author     Tom Goodale
   @desc
               For a given variable name, return its associated global index.
   @enddesc

   @returntype int
   @returndesc
               the index of the given variable, or negative for failure:
               -1 if no variable of such name exists, or
                  error return code of @seeroutine CCTK_DecomposeName
   @endreturndesc
@@*/
int CCTK_VarIndex (const char *fullvarname)
{
  int retval,
      group,
      variable;
  char *impname = NULL;
  char *varname = NULL;


  impname = varname = NULL;
  retval = CCTK_DecomposeName (fullvarname, &impname, &varname);
  if (! retval)
  {
    retval = -1;
    for (group = 0; group < n_groups && retval < 0; group++)
    {
      if (CCTK_Equals (impname, groups[group].implementation))
      {
        for (variable = 0; variable < groups[group].n_variables; variable++)
        {
          if (CCTK_Equals (varname, groups[group].variables[variable].name))
          {
            retval = groups[group].variables[variable].number;
            break;
          }
        }
      }
    }

    if (retval < 0)
    {
      CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                  "CCTK_VarIndex: No variable named '%s' found",
                  fullvarname);
    }
  }

#ifdef DEBUG_GROUPS
  printf (" In VarIndex\n", " ------------\n");
  printf ("   impname -%s-\n", impname);
  printf ("   varname -%s-\n", varname);
#endif

  /* free strings allocated in CCTK_DecomposeName() */
  free (impname);
  free (varname);

  return (retval);
}

void CCTK_FCALL CCTK_FNAME (CCTK_VarIndex)
                           (int *vindex, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (name)
  *vindex = CCTK_VarIndex (name);
  free (name);
}


 /*@@
   @routine    CCTKi_VarDataPtr
   @date       Wed 19 June 2002
   @author     Thomas Radke
   @desc
               For a grid variable given by its timelevel, implementation,
               and name, return its data pointer on the cGH.
               This internal function is only called by the variable bindings
               in the DECLARE_CCTK_ARGUMENTS macro. It does not print a warning
               if the variable doesn't exist (eg. because the providing thorn
               wasn't activated) - it passes back a NULL pointer in this case.
               <P>
               Note that for performance reasons the variable's full name
               is passed in as separate implementation and varname strings
               (to save a call to CCTK_Decompose()), and we also don't use
               CCTK_Equals() for the string comparisons but make use of the
               implicit assumption that the CST-generated macros pass in the
               names in the same notation as they are stored in the database.
   @enddesc

   @returntype void *
   @returndesc
               data pointer on the cGH, or NULL if variable doesn't exist
   @endreturndesc
@@*/
void *CCTKi_VarDataPtr(const cGH *GH, int timelevel,
                       const char *implementation, const char *varname)
{
  int group, var;


  for (group = 0; group < n_groups; group++)
  {
    if (strcmp (implementation, groups[group].implementation) == 0)
    {
      for (var = 0; var < groups[group].n_variables; var++)
      {
        if (strcmp (varname, groups[group].variables[var].name) == 0)
        {
          return (GH->data[groups[group].variables[var].number][timelevel]);
        }
      }
    }
  }

  return (NULL);
}


 /*@@
   @routine    CCTK_MaxDim
   @date       Mon Feb  8 12:04:01 1999
   @author     Tom Goodale
   @desc
               Gets the maximum dimension of all groups.
   @enddesc

   @returntype int
   @returndesc
               the maximum dimension of all groups
   @endreturndesc
@@*/
int CCTK_MaxDim (void)
{
  return (maxdim);
}

void CCTK_FCALL CCTK_FNAME (CCTK_MaxDim)
                           (int *dim)
{
  *dim = CCTK_MaxDim ();
}


 /*@@
   @routine    CCTK_MaxGFDim
   @date       2008-03-19
   @author     Erik Schnetter
   @desc
               Gets the maximum dimension of all grid function groups.
   @enddesc

   @returntype int
   @returndesc
               the maximum dimension of all grid function groups
   @endreturndesc
@@*/
int CCTK_MaxGFDim (void)
{
  return (gfdim);
}

void CCTK_FCALL CCTK_FNAME (CCTK_MaxGFDim)
                           (int *dim)
{
  *dim = CCTK_MaxGFDim ();
}


 /*@@
   @routine    CCTK_NumVars
   @date       Mon Feb  8 12:04:50 1999
   @author     Tom Goodale
   @desc
               Gets the total number of variables.
   @enddesc

   @returntype int
   @returndesc
               total number of variables created so far
   @endreturndesc
@@*/
int CCTK_NumVars (void)
{
  return (total_variables);
}

void CCTK_FCALL CCTK_FNAME (CCTK_NumVars)
                           (int *num_vars)
{
  *num_vars = CCTK_NumVars ();
}


 /*@@
   @routine    CCTK_NumGroups
   @date       Mon Feb  8 12:04:50 1999
   @author     Tom Goodale
   @desc
               Gets the total number of groups.
   @enddesc

   @returntype int
   @returndesc
               total number of groups created so far
   @endreturndesc
@@*/
int CCTK_NumGroups (void)
{
  return (n_groups);
}

void CCTK_FCALL CCTK_FNAME (CCTK_NumGroups)
                           (int *num_groups)
{
  *num_groups = CCTK_NumGroups ();
}


 /*@@
   @routine    CCTK_GroupNameFromVarI
   @date       Mon Feb 22
   @author     Gabrielle Allen
   @desc
               Given a variable index return a group name.
   @enddesc

   @returntype char *
   @returndesc
               the full group name of the given variable (which should be freed
               if not needed anymore), or
               NULL if the given variable index is invalid, or out of memory
   @endreturndesc
@@*/
char *CCTK_GroupNameFromVarI (int var)
{
  int group;
  char *fullname;


  if (0 <= var && var < total_variables)
  {
    group = group_of_variable[var];
    fullname = malloc (strlen (groups[group].name) +
                       strlen (groups[group].implementation) + 3);
    if (fullname)
    {
      sprintf (fullname, "%s::%s",
               groups[group].implementation, groups[group].name);
    }
  }
  else
  {
    fullname = NULL;
  }

  return (fullname);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupNameFromVarI) (int *nchars, int *var, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_PTR (groupnameptr)
  ONE_FORTSTRING_LEN (groupnamelen)
  char *groupname;


  groupname = CCTK_GroupNameFromVarI (*var);
  *nchars = CCTK_FortranString (groupname ? groupname : "", groupnameptr,
                                groupnamelen);
  free (groupname);
}


 /*@@
   @routine    CCTK_GroupIndexFromVarI
   @date       Mon Feb 22
   @author     Gabrielle Allen
   @desc
               Given a variable index return a group index.
   @enddesc

   @returntype int
   @returndesc
               the group index of the given variable, or
               -1 if the given variable index is invalid
   @endreturndesc
@@*/
int CCTK_GroupIndexFromVarI (int var)
{
  return ((0 <= var && var < total_variables) ? group_of_variable[var] : -1);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupIndexFromVarI)
                           (int *gindex, const int *var)
{
  *gindex = CCTK_GroupIndexFromVarI (*var);
}



 /*@@
   @routine    CCTK_GroupIndexFromVar
   @date       Mon Feb 22
   @author     Gabrielle Allen
   @desc
               Given a variable name returns a group index.
   @enddesc

   @returntype int
   @returndesc
               return code of @seeroutine CCTK_GroupIndexFromVarI
               -1 if the given variable index is invalid
   @endreturndesc
@@*/
int CCTK_GroupIndexFromVar (const char *var)
{
  return CCTK_GroupIndexFromVarI (CCTK_VarIndex (var));
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupIndexFromVar)
                           (int *vindex, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (var)
  *vindex = CCTK_GroupIndexFromVar (var);
  free (var);
}


 /*@@
   @routine    CCTK_ImpFromVarI
   @date       Mon Feb 22
   @author     Gabrielle Allen
   @desc
               Given a variable index return the implementation name.
   @enddesc

   @returntype const char *
   @returndesc
               the implementation name of the given variable, or
               NULL if given variable index is invalid
   @endreturndesc
@@*/
const char *CCTK_ImpFromVarI (int var)
{
  return ((0 <= var && var < total_variables) ?
          groups[group_of_variable[var]].implementation : NULL);
}

void CCTK_FCALL CCTK_FNAME (CCTK_ImpFromVarI) (int *nchars, int *var, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_PTR (impptr)
  ONE_FORTSTRING_LEN (implen)
  const char *imp;


  imp = CCTK_ImpFromVarI (*var);
  if (!imp)
  {
    *nchars = CCTK_FortranString ("", impptr, implen);
  }
  else
  {
    *nchars = CCTK_FortranString (imp, impptr, implen);
  }
}


 /*@@
   @routine    CCTK_FullName
   @date       Mon Feb 22
   @author     Gabrielle Allen
   @desc
               Given a variable index return the variable's full name,
               ie. <implementation name>::<variable name>.
   @enddesc

   @returntype char *
   @returndesc
               the full name of the given variable (which should be freed
               if not needed anymore), or
               NULL if given variable index is invalid, or out out memory
   @endreturndesc
@@*/
char *CCTK_FullName (int var)
{
  const char *impname,
             *varname;
  char *fullname;


  varname = CCTK_VarName (var);
  if (varname)
  {
    impname = groups[group_of_variable[var]].implementation;
    fullname = malloc (strlen (varname) + strlen (impname) + 3);
    if (fullname)
    {
      sprintf (fullname, "%s::%s", impname, varname);
    }
  }
  else
  {
    fullname = NULL;
  }

  return (fullname);
}

void CCTK_FCALL CCTK_FNAME (CCTK_FullName) (int *nchars, int *var, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_PTR (fullnameptr)
  ONE_FORTSTRING_LEN (fullnamelen)
  char *fullname;


  fullname = CCTK_FullName (*var);
  *nchars = CCTK_FortranString (fullname ? fullname : "", fullnameptr,
                                fullnamelen);
  free (fullname);
}


 /*@@
   @routine    CCTK_GroupTypeNumber
   @date       Mon Feb  8 14:44:45 1999
   @author     Tom Goodale
   @desc
               Gets the type number associated with a group.
   @enddesc

   @returntype int
   @returndesc
                the type number of the given group type, or
               -1 if given group type is invalid
   @endreturndesc
@@*/
int CCTK_GroupTypeNumber (const char *type)
{
  int retval;


  if (! strcmp (type, "SCALAR"))
  {
    retval = CCTK_SCALAR;
  }
  else if (! strcmp (type, "GF"))
  {
    retval = CCTK_GF;
  }
  else if (! strcmp (type, "ARRAY"))
  {
    retval = CCTK_ARRAY;
  }
  else
  {
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupTypeNumber)
                           (int *number, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (type)
  *number = CCTK_GroupTypeNumber (type);
  free (type);
}


 /*@@
   @routine    CCTK_VarTypeNumber
   @date       Mon Feb  8 14:44:45 1999
   @author     Tom Goodale
   @desc
               Gets the type number associated with a variable.
   @enddesc

   @returntype int
   @returndesc
                the type number of the given variable type, or
               -1 if given variable type is invalid
   @endreturndesc
@@*/
int CCTK_VarTypeNumber (const char *type)
{
  int retval;


  if (! strcmp (type, "BYTE"))
  {
    retval = CCTK_VARIABLE_BYTE;
  }
  else if (! strcmp (type, "INT"))
  {
    retval = CCTK_VARIABLE_INT;
  }
  else if (! strcmp (type, "INT1"))
  {
    retval = CCTK_VARIABLE_INT1;
  }
  else if (! strcmp (type, "INT2"))
  {
    retval = CCTK_VARIABLE_INT2;
  }
  else if (! strcmp (type, "INT4"))
  {
    retval = CCTK_VARIABLE_INT4;
  }
  else if (! strcmp (type, "INT8"))
  {
    retval = CCTK_VARIABLE_INT8;
  }
  else if (! strcmp (type, "INT16"))
  {
    retval = CCTK_VARIABLE_INT16;
  }
  else if (! strcmp (type, "REAL"))
  {
    retval = CCTK_VARIABLE_REAL;
  }
  else if (! strcmp (type, "REAL4"))
  {
    retval = CCTK_VARIABLE_REAL4;
  }
  else if (! strcmp (type, "REAL8"))
  {
    retval = CCTK_VARIABLE_REAL8;
  }
  else if (! strcmp (type, "REAL16"))
  {
    retval = CCTK_VARIABLE_REAL16;
  }
  else if (! strcmp (type, "COMPLEX"))
  {
    retval = CCTK_VARIABLE_COMPLEX;
  }
  else if (! strcmp (type, "COMPLEX8"))
  {
    retval = CCTK_VARIABLE_COMPLEX8;
  }
  else if (! strcmp (type, "COMPLEX16"))
  {
    retval = CCTK_VARIABLE_COMPLEX16;
  }
  else if (! strcmp (type, "COMPLEX32"))
  {
    retval = CCTK_VARIABLE_COMPLEX32;
  }
  else
  {
    retval = -1;
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_VarTypeNumber)
                           (int *number, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (type)
  *number = CCTK_VarTypeNumber (type);
  free (type);
}


 /*@@
   @routine    CCTK_VarTypeName
   @date       Mon Jan  3 13:50:56 CET 2000
   @author     Gabrielle Allen
   @desc
               Gets the variable type name associated with a variable type.
   @enddesc

   @returntype const char *
   @returndesc
                the type name of the given variable type, or
               -1 if given variable type is invalid
   @endreturndesc
@@*/
const char *CCTK_VarTypeName (int vtype)
{
  const char *retval;


  switch (vtype)
  {
    case CCTK_VARIABLE_BYTE:
      retval = "CCTK_VARIABLE_BYTE";
      break;

    case CCTK_VARIABLE_INT:
      retval = "CCTK_VARIABLE_INT";
      break;

    case CCTK_VARIABLE_INT1:
      retval = "CCTK_VARIABLE_INT1";
      break;

    case CCTK_VARIABLE_INT2:
      retval = "CCTK_VARIABLE_INT2";
      break;

    case CCTK_VARIABLE_INT4:
      retval = "CCTK_VARIABLE_INT4";
      break;

    case CCTK_VARIABLE_INT8:
      retval = "CCTK_VARIABLE_INT8";
      break;

    case CCTK_VARIABLE_REAL:
      retval = "CCTK_VARIABLE_REAL";
      break;

    case CCTK_VARIABLE_REAL4:
      retval = "CCTK_VARIABLE_REAL4";
      break;

    case CCTK_VARIABLE_REAL8:
      retval = "CCTK_VARIABLE_REAL8";
      break;

    case CCTK_VARIABLE_COMPLEX:
      retval = "CCTK_VARIABLE_COMPLEX";
      break;

    case CCTK_VARIABLE_COMPLEX8:
      retval = "CCTK_VARIABLE_COMPLEX8";
      break;

    case CCTK_VARIABLE_COMPLEX16:
      retval = "CCTK_VARIABLE_COMPLEX16";
      break;

    case CCTK_VARIABLE_COMPLEX32:
      retval = "CCTK_VARIABLE_COMPLEX32";
      break;

    case CCTK_VARIABLE_CHAR:
      retval = "CCTK_VARIABLE_CHAR";
      break;

    case CCTK_VARIABLE_STRING:
      retval = "CCTK_VARIABLE_STRING";
      break;

    case CCTK_VARIABLE_POINTER:
      retval = "CCTK_VARIABLE_POINTER";
      break;

    case CCTK_VARIABLE_POINTER_TO_CONST:
      retval = "CCTK_VARIABLE_POINTER_TO_CONST";
      break;

    case CCTK_VARIABLE_FPOINTER:
      retval = "CCTK_VARIABLE_FPOINTER";
      break;

    default:
      retval = NULL;
      break;
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_VarTypeName) (int *nchars, int *var, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_PTR (vartypenameptr)
  ONE_FORTSTRING_LEN (vartypenamelen)
  const char *vartypename;


  vartypename = CCTK_VarTypeName (*var);
  if (!vartypename)
  {
    *nchars = CCTK_FortranString ("", vartypenameptr, vartypenamelen);
  }
  else
  {
    *nchars = CCTK_FortranString (vartypename, vartypenameptr, vartypenamelen);
  }
}


 /*@@
   @routine    CCTK_GroupScopeNumber
   @date       Tuesday June 22 1999
   @author     Gabrielle Allen
   @desc
               Gets the scope number associated with a group.
   @enddesc

   @returntype int
   @returndesc
                the scope number of the given scope type, or
               -1 if given scope type is invalid
   @endreturndesc
@@*/
int CCTK_GroupScopeNumber (const char *type)
{
  int retval;


  if (! strcmp (type, "PRIVATE"))
  {
    retval = CCTK_PRIVATE;
  }
  else if (! strcmp (type, "PROTECTED"))
  {
    retval = CCTK_PROTECTED;
  }
  else if (! strcmp (type, "PUBLIC"))
  {
    retval = CCTK_PUBLIC;
  }
  else
  {
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupScopeNumber)
                           (int *number, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (type)
  *number = CCTK_GroupScopeNumber (type);
  free (type);
}


 /*@@
   @routine    CCTK_GroupDistribNumber
   @date       Tuesday June 22 1999
   @author     Gabrielle Allen
   @desc
               Gets the distribution number associated with a group.
   @enddesc

   @returntype int
   @returndesc
                the distribution number of the given distribution type, or
               -1 if given distribution type is invalid
   @endreturndesc
@@*/
int CCTK_GroupDistribNumber (const char *dtype)
{
  int retval;


  if (! strcmp (dtype, "CONSTANT"))
  {
    retval = CCTK_DISTRIB_CONSTANT;
  }
  else if (! strcmp (dtype, "DEFAULT"))
  {
    retval = CCTK_DISTRIB_DEFAULT;
  }
  else
  {
    retval = -1;
  }

  return (retval);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupDistribNumber)
                           (int *number, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (type)
  *number = CCTK_GroupDistribNumber (type);
  free (type);
}


 /*@@
   @routine    CCTK_GroupData
   @date       Mon Feb  8 15:56:01 1999
   @author     Tom Goodale
   @desc
               For a given group index, gets the group type, variable type,
               and the number of variables
   @enddesc

   @returntype int
   @returndesc
                0 for success
               -1 if given group index is invalid
               -2 if given pointer to store group data is NULL
   @endreturndesc
@@*/
int CCTK_GroupData (int group, cGroup *gp)
{
  int retval;


  retval = (0 <= group && group < n_groups) ? 0 : -1;
  if (! retval)
  {
    if (gp)
    {
      gp->grouptype     = groups[group].gtype;
      gp->vartype       = groups[group].vtype;
      gp->disttype      = groups[group].dtype;
      gp->dim           = groups[group].dim;
      gp->numvars       = groups[group].n_variables;
      gp->vectorlength  = groups[group].vectorlength;
      gp->numtimelevels = groups[group].n_timelevels;
      gp->tagstable     = groups[group].tags_table;

      if(groups[group].vararraysize)
      {
        gp->vectorgroup = 1;
      }
      else
      {
        gp->vectorgroup = 0;
      }
    }
    else
    {
      retval = -2;
    }
  }

  return (retval);
}


 /*@@
   @routine    CCTK_VarName
   @date       Tue Feb  9 15:34:56 1999
   @author     Tom Goodale
   @desc
               Gets the name of a variable.
   @enddesc

   @returntype const char *
   @returndesc
               the name of the given variable, or
               -1 if given variable index is invalid
   @endreturndesc
@@*/
const char *CCTK_VarName (int var)
{
  return ((0 <= var && var < total_variables) ?
          groups[group_of_variable[var]]
            .variables[var-groups[group_of_variable[var]].variables[0].number]
            .name
          : NULL);
}

void CCTK_FCALL CCTK_FNAME (CCTK_VarName) (int *nchars, int *var, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_PTR (varnameptr)
  ONE_FORTSTRING_LEN (varnamelen)
  const char *varname;


  varname = CCTK_VarName (*var);
  if (!varname)
  {
    *nchars = CCTK_FortranString ("", varnameptr, varnamelen);
  }
  else
  {
    *nchars = CCTK_FortranString (varname, varnameptr, varnamelen);
  }
}


 /*@@
   @routine    CCTK_DecomposeName
   @date       Tue Feb  9 15:39:14 1999
   @author     Tom Goodale
   @desc
               Decomposes a full group or variable name of the form imp::name
   @enddesc

   @returntype int
   @returndesc
                0 for success (implementation and name are set to the
                  full name's implementation and name), or
                  negative otherwise (a non-zero error return code of
                  @seeroutine Util_SplitString is translated into one of the
                  following error codes:
               -2 if failed to catch error code from Util_SplitString
               -3 if given full name is in wrong format
               -4 if memory allocation failed
   @endreturndesc
@@*/
int CCTK_DecomposeName (const char *fullname,
                        char **implementation,
                        char **name)
{
  int retval;


  retval = Util_SplitString (implementation, name, fullname, "::");
  if (retval)
  {
    if (retval == 1)
    {
      CCTK_VWarn (8, __LINE__, __FILE__, "Cactus",
                  "CCTK_DecomposeName: Full name '%s' in wrong format",
                  fullname);
      retval = -3;
    }
    else if (retval == 2)
    {
      CCTK_Warn (2, __LINE__, __FILE__, "Cactus",
                 "CCTK_DecomposeName: Memory allocation failed");
      retval = -4;
    }
    else
    {
      CCTK_Warn (1, __LINE__, __FILE__, "Cactus",
                 "CCTK_DecomposeName: Error failed to be caught");
      retval = -2;
    }
  }

  return (retval);
}


 /*@@
   @routine    CCTK_GroupName
   @date       Tue Apr  9 15:39:14 1999
   @author     Gabrielle Allen
   @desc
               Given a group index returns the group name
   @enddesc

   @returntype char *
   @returndesc
               the full name of the given group (which should be freed
               if not needed anymore), or
               -1 if given group index is invalid
   @endreturndesc
@@*/
char *CCTK_GroupName (int group)
{
  char *name;


  name = NULL;
  if (0 <= group && group < n_groups)
  {
    name = malloc (strlen (groups[group].implementation) +
                   strlen (groups[group].name) + 3);
    if (name)
    {
      sprintf (name, "%s::%s",
               groups[group].implementation, groups[group].name);
    }
  }

  return (name);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupName) (int *nchars, int *var, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_PTR (groupnameptr)
  ONE_FORTSTRING_LEN (groupnamelen)
  char *groupname;


  groupname = CCTK_GroupName (*var);
  *nchars = CCTK_FortranString (groupname ? groupname : "", groupnameptr,
                                groupnamelen);
  free (groupname);
}


 /*@@
   @routine    CCTK_GroupTagsTable
   @date       Wed May 22 00:47:58 2002
   @author     Tom Goodale
   @desc
   Returns the table handle of the TAGS table for a group.
   @enddesc
   @var     groupname
   @vdesc   The group name
   @vtype   int
   @vio     in
   @endvar

   @returntype int
   @returndesc
               -1 if group index out of range
   @endreturndesc
 @@*/
int CCTK_GroupTagsTable(const char *groupname)
{
  int retval;
  int group;

  group = CCTK_GroupIndex (groupname);

  retval = CCTK_GroupTagsTableI(group);

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupTagsTable) (int *table, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)


  *table = CCTK_GroupTagsTable (groupname);
  free (groupname);
}


 /*@@
   @routine    CCTK_GroupTagsTableI
   @date       Wed May 22 00:47:58 2002
   @author     Tom Goodale
   @desc
   Returns the table handle of the TAGS table for a group.
   @enddesc
   @var     group
   @vdesc   The group index
   @vtype   int
   @vio     in
   @endvar

   @returntype int
   @returndesc
               -1 if group index out of range
   @endreturndesc
 @@*/
int CCTK_GroupTagsTableI(int group)
{
  int retval;

  if (0 <= group && group < n_groups)
  {
    retval = groups[group].tags_table;
  }
  else
  {
    retval = -1;
  }

  return retval;
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupTagsTableI) (int *table, int *group)
{
  *table = CCTK_GroupTagsTableI (*group);
}

 /*@@
   @routine    CCTK_FirstVarIndexI
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a group index returns the first variable index in the group
   @enddesc

   @returntype int
   @returndesc
               the index of the first variable in the given group, or
               -1 if given group index is invalid
               -2 if given group has no members
   @endreturndesc
@@*/
int CCTK_FirstVarIndexI (int group)
{
  if (! (0 <= group && group < n_groups)) return -1;
  if (groups[group].n_variables == 0) return -2;
  return groups[group].variables[0].number;
}

void CCTK_FCALL CCTK_FNAME (CCTK_FirstVarIndexI)
                           (int *first, const int *group)
{
  *first = CCTK_FirstVarIndexI (*group);
}


 /*@@
   @routine    CCTK_FirstVarIndex
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a group name returns the first variable index in the group
   @enddesc

   @returntype int
   @returndesc
               the index of the first variable in the given group, or
               -1 if given group name is invalid
               -2 if given group has no members
   @endreturndesc
@@*/
int CCTK_FirstVarIndex (const char *groupname)
{
  return CCTK_FirstVarIndexI (CCTK_GroupIndex (groupname));
}

void CCTK_FCALL CCTK_FNAME (CCTK_FirstVarIndex)
                           (int *first, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *first = CCTK_FirstVarIndex (groupname);
  free (groupname);
}


 /*@@
   @routine    CCTK_NumVarsInGroupI
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a group index returns the number of variables in the group
   @enddesc

   @returntype int
   @returndesc
               the number of variables in the given group, or
               -1 if given group index is invalid
   @endreturndesc
@@*/
int CCTK_NumVarsInGroupI (int group)
{
  return ((0 <= group && group < n_groups) ? groups[group].n_variables : -1);
}

void CCTK_FCALL CCTK_FNAME(CCTK_NumVarsInGroupI) (int *numvars, int *group)
{
  *numvars = CCTK_NumVarsInGroupI (*group);
}


 /*@@
   @routine    CCTK_NumVarsInGroup
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a group name returns the number of variables in the group
   @enddesc

   @returntype int
   @returndesc
               return code of @seeroutine CCTK_NumVarsInGroupI
               -1 if given group name is invalid
   @endreturndesc
@@*/
int CCTK_NumVarsInGroup (const char *groupname)
{
  return CCTK_NumVarsInGroupI (CCTK_GroupIndex (groupname));
}

void CCTK_FCALL CCTK_FNAME (CCTK_NumVarsInGroup)
                           (int *num, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (groupname)
  *num = CCTK_NumVarsInGroup (groupname);
  free (groupname);
}


 /*@@
   @routine    CCTK_GroupTypeFromVarI
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a variable index return the type of group
   @enddesc

   @returntype int
   @returndesc
               the type of the given group, or
               -1 if given variable index is invalid
   @endreturndesc
@@*/
int CCTK_GroupTypeFromVarI (int var)
{
  return ((0 <= var && var < total_variables) ?
          groups[group_of_variable[var]].gtype : -1);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupTypeFromVarI)
                           (int *type, const int *var)
{
  *type = CCTK_GroupTypeFromVarI (*var);
}


 /*@@
   @routine    CCTK_GroupTypeI
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a group index return the type of group
   @enddesc

   @returntype int
   @returndesc
               the type of the given group, or
               -1 if given group index is invalid
   @endreturndesc
@@*/
int CCTK_GroupTypeI (int group)
{
  return ((0 <= group && group < n_groups) ? groups[group].gtype : -1);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupTypeI)
                           (int *type, const int *group)
{
  *type = CCTK_GroupTypeI (*group);
}


 /*@@
   @routine    CCTK_VarTypeI
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a variable index return the variable type
   @enddesc

   @returntype int
   @returndesc
               the type of the given variable, or
               -1 if given variable index is invalid
   @endreturndesc
@@*/
int CCTK_VarTypeI (int var)
{
  return ((0 <= var && var < total_variables) ?
          groups[group_of_variable[var]].vtype : -1);
}

void CCTK_FCALL CCTK_FNAME (CCTK_VarTypeI)
                           (int *type, const int *var)
{
  *type = CCTK_VarTypeI (*var);
}


 /*@@
   @routine    CCTK_DeclaredTimeLevelsGI
   @date       July 16 2003
   @author     Gabrielle Allen
   @desc
               Given a group index return the maximum number of timelevels
   @enddesc

   @returntype int
   @returndesc
               the number of timelevels of variables in the group, or
               -1 if given group index is invalid
   @endreturndesc
@@*/
int CCTK_DeclaredTimeLevelsGI (int group)
{
  return ((0 <= group && group < n_groups) ?
          groups[group].n_timelevels : -1);
}

void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevelsGI)
                           (int *num, const int *group)
{
  *num = CCTK_DeclaredTimeLevelsGI (*group);
}


 /*@@
   @routine    CCTK_DeclaredTimeLevelsVI
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a variable index return the maximum number of timelevels
   @enddesc

   @returntype int
   @returndesc
               the number of timelevels of the given variable, or
               -1 if given variable index is invalid
   @endreturndesc
@@*/
int CCTK_DeclaredTimeLevelsVI (int var)
{
  return ((0 <= var && var < total_variables) ?
          groups[group_of_variable[var]].n_timelevels : -1);
}

void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevelsVI)
                           (int *num, const int *var)
{
  *num = CCTK_DeclaredTimeLevelsVI (*var);
}


 /*@@
   @routine    CCTK_DeclaredTimeLevels
   @date       8 June 2003
   @author     Gabrielle Allen
   @desc
               Given a group name return the maximum number of timelevels
   @enddesc

   @returntype int
   @returndesc
               return code of @seeroutine CCTK_DeclaredTimeLevelsI
   @endreturndesc
@@*/
int CCTK_DeclaredTimeLevels (const char *group)
{
  return CCTK_DeclaredTimeLevelsGI (CCTK_GroupIndex (group));
}

void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevels)
                           (int *num, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (group)
  *num = CCTK_DeclaredTimeLevels (group);
  free (group);
}

 /*@@
   @routine    CCTK_DeclaredTimeLevelsGN
   @date       8 June 2003
   @author     Gabrielle Allen
   @desc
               Given a group name return the number of timelevels
   @enddesc

   @returntype int
   @returndesc
               return code of @seeroutine CCTK_DeclaredTimeLevelsGI
   @endreturndesc
@@*/
int CCTK_DeclaredTimeLevelsGN (const char *group)
{
  return CCTK_DeclaredTimeLevelsGI (CCTK_GroupIndex (group));
}

void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevelsGN)
                           (int *num, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (group)
  *num = CCTK_DeclaredTimeLevelsGN (group);
  free (group);
}

 /*@@
   @routine    CCTK_DeclaredTimeLevelsVN
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a variable name return the maximum number of timelevels
   @enddesc

   @returntype int
   @returndesc
               return code of @seeroutine CCTK_DeclaredTimeLevelsVI
   @endreturndesc
@@*/
int CCTK_DeclaredTimeLevelsVN (const char *var)
{
  return CCTK_DeclaredTimeLevelsVI (CCTK_VarIndex (var));
}

void CCTK_FCALL CCTK_FNAME (CCTK_DeclaredTimeLevelsVN)
                           (int *num, ONE_FORTSTRING_ARG)
{
  ONE_FORTSTRING_CREATE (var)
  *num = CCTK_DeclaredTimeLevelsVN (var);
  free (var);
}




 /*@@
   @routine    CCTK_PrintGroup
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a group index print the group name.
               This is for debugging purposes for Fortran routines.
   @enddesc
@@*/
void CCTK_FCALL CCTK_FNAME (CCTK_PrintGroup)
                           (const int *group)
{
  fprintf (stdout, "Group %d is %s\n", *group, CCTK_GroupName (*group));
}


 /*@@
   @routine    CCTK_PrintVar
   @date       3 July 1999
   @author     Gabrielle Allen
   @desc
               Given a group index print the variable name.
               This is for debugging purposes for Fortran.
   @enddesc
@@*/
void CCTK_FCALL CCTK_FNAME (CCTK_PrintVar)
                           (const int *var)
{
  fprintf (stdout, "Variable %d is %s\n", *var, CCTK_VarName (*var));
}


 /*@@
   @routine    CCTK_GroupSizesI
   @date       Sun Nov 28 12:56:44 1999
   @author     Tom Goodale
   @desc
               Returns the size array for a group.
   @enddesc

   @returntype CCTK_INT **
   @returndesc
               the pointer to the size array of the given group, or
               -1 if invalid group index was given
   @endreturndesc
@@*/
CCTK_INT **CCTK_GroupSizesI (int group)
{
  return ((0 <= group && group < n_groups) ? groups[group].size : NULL);
}


 /*@@
   @routine    CCTK_GroupGhostsizesI
   @date       Sun Jan 23 12:56:44 2000
   @author     Gabrielle Allen
   @desc
               Returns the ghostsize array for a group.
   @enddesc

   @returntype CCTK_INT **
   @returndesc
               the pointer to the ghostsize array of the given group, or
               NULL if invalid group index was given or the group exists
               but has no ghostsize information
   @endreturndesc
@@*/
CCTK_INT **CCTK_GroupGhostsizesI (int group)
{
  return ((0 <= group && group < n_groups) ? groups[group].ghostsize : NULL);
}


 /*@@
   @routine    CCTK_VarTypeSize
   @date       Sun Dec  5 10:08:05 1999
   @author     Gabrielle Allen
   @desc
               Returns the size of a given variable type
   @enddesc

   @returntype int
   @returndesc
               positive for the variable type's size (in bytes), or
               -1 if invalid variable type was given
   @endreturndesc
@@*/
int CCTK_VarTypeSize (int vtype)
{
  int var_size;


  switch (vtype)
  {
    case CCTK_VARIABLE_BYTE:
      var_size = sizeof (CCTK_BYTE);
      break;

    case CCTK_VARIABLE_INT:
      var_size = sizeof (CCTK_INT);
      break;

    case CCTK_VARIABLE_REAL:
      var_size = sizeof (CCTK_REAL);
      break;

    case CCTK_VARIABLE_COMPLEX:
      var_size = sizeof (CCTK_COMPLEX);
      break;

#ifdef HAVE_CCTK_INT1
    case CCTK_VARIABLE_INT1:
      var_size = sizeof (CCTK_INT1);
      break;
#endif

#ifdef HAVE_CCTK_INT2
    case CCTK_VARIABLE_INT2:
      var_size = sizeof (CCTK_INT2);
      break;
#endif

#ifdef HAVE_CCTK_INT4
    case CCTK_VARIABLE_INT4:
      var_size = sizeof (CCTK_INT4);
      break;
#endif

#ifdef HAVE_CCTK_INT8
    case CCTK_VARIABLE_INT8:
      var_size = sizeof (CCTK_INT8);
      break;
#endif

#ifdef HAVE_CCTK_INT16
    case CCTK_VARIABLE_INT16:
      var_size = sizeof (CCTK_INT16);
      break;
#endif

#ifdef HAVE_CCTK_REAL4
    case CCTK_VARIABLE_REAL4:
      var_size = sizeof (CCTK_REAL4);
      break;

    case CCTK_VARIABLE_COMPLEX8:
      var_size = sizeof (CCTK_COMPLEX8);
      break;
#endif

#ifdef HAVE_CCTK_REAL8
    case CCTK_VARIABLE_REAL8:
      var_size = sizeof (CCTK_REAL8);
      break;

    case CCTK_VARIABLE_COMPLEX16:
      var_size = sizeof (CCTK_COMPLEX16);
      break;
#endif

#ifdef HAVE_CCTK_REAL16
    case CCTK_VARIABLE_REAL16:
      var_size = sizeof (CCTK_REAL16);
      break;

    case CCTK_VARIABLE_COMPLEX32:
      var_size = sizeof (CCTK_COMPLEX32);
      break;
#endif

    case CCTK_VARIABLE_CHAR:
      var_size = sizeof (CCTK_CHAR);
      break;

    case CCTK_VARIABLE_POINTER:
      var_size = sizeof (CCTK_POINTER);
      break;

    case CCTK_VARIABLE_POINTER_TO_CONST:
      var_size = sizeof (CCTK_POINTER_TO_CONST);
      break;

    case CCTK_VARIABLE_FPOINTER:
      var_size = sizeof (CCTK_FPOINTER);
      break;

    default:
      CCTK_VWarn (4, __LINE__, __FILE__, "Cactus",
                  "CCTK_VarTypeSize: Unknown variable type (%d)", vtype);
      var_size = -1;
  }

  return (var_size);
}

void CCTK_FCALL CCTK_FNAME(CCTK_VarTypeSize) (int *size, int *vtype)
{
  *size = CCTK_VarTypeSize (*vtype);
}


 /*@@
   @routine    CCTK_GroupDimI
   @date       Wed Feb 2 2000
   @author     Gabrielle Allen
   @desc
               Given a group index returns the group dimension
   @enddesc

   @returntype int
   @returndesc
               the dimension of the given group, or
               -1 if given group index is invalid
   @endreturndesc
@@*/
int CCTK_GroupDimI (int group)
{
  return ((0 <= group && group < n_groups) ? groups[group].dim : -1);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupDimI)
                           (int *dim, const int *group)
{
  *dim = CCTK_GroupDimI (*group);
}


 /*@@
   @routine    CCTK_GroupDimFromVarI
   @date       Wed Feb 2 2000
   @author     Gabrielle Allen
   @desc
               Given a variable index returns the group dimension
   @enddesc

   @returntype int
   @returndesc
               the dimension of the variable's group, or
               -1 if given variable index is invalid
   @endreturndesc
@@*/
int CCTK_GroupDimFromVarI (int var)
{
  return ((0 <= var && var < total_variables) ?
          groups[group_of_variable[var]].dim : -1);
}

void CCTK_FCALL CCTK_FNAME (CCTK_GroupDimFromVarI)
                           (int *dim, const int *var)
{
  *dim = CCTK_GroupDimFromVarI (*var);
}


 /*@@
   @routine    CCTK_TraverseString
   @date       Wed 20 Sep 2000
   @author     Thomas Radke
   @desc
               Traverse through all variables and/or groups whose names
               appear in the given string, and call the callback routine
               with those indices and an optional option string appended
               to the variable/group name enclosed in square braces.
               The special keyword "all" in the string can be used to
               indicate that the callback should be called for all
               variables/groups.
   @enddesc
   @var        traverse_string
   @vdesc      list of variable and/or group names
   @vtype      const char *
   @vio        in
   @endvar
   @var        callback
   @vdesc      routine to call for every variable and/or group found
   @vtype      int (*) (int idx, const char *optstring, void *callback_arg)
   @vio        int
   @endvar
   @var        callback_arg
   @vdesc      an arbitrary argument which gets passed to the callback routine
   @vtype      void *
   @vio        in
   @endvar
   @var        selection
   @vdesc      decides whether group and/or variable names are accepted
               in the string
   @vtype      int
   @vio        in
   @endvar

   @returntype int
   @returndesc
               positive for the number of traversed variables, or
               -1 no callback routine was given<BR>
               -2 option string is not associated with a group or variable<BR>
               -3 unterminated option string<BR>
               -4 garbage found at end of option string<BR>
               -5 invalid token in traversed string found
   @endreturndesc
@@*/
int CCTK_TraverseString (const char *traverse_string,
                         void (*callback) (int idx,
                                           const char *optstring,
                                           void *callback_arg),
                         void *callback_arg,
                         int selection)
{
  int retval, nesting, vindex, gindex, first, last, selected_all;
  char delimiter, options_start, options_end;
  char *tmp, *string, *parse_string, *group_var_string, *option_string;
  cGroup gdata;


  if (callback == NULL)
  {
    CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                "CCTK_TraverseString: No callback given");
    return (-1);
  }

  retval = 0;

  /* create a work copy of the string to traverse
     which we can edit in-place */
  parse_string = strdup (traverse_string);

  /* parse it token by token */
  string = parse_string;
  while (string && *string)
  {

    /* skip leading spaces */
    while (*string && isspace (*string))
    {
      string++;
    }
    if (! *string)
    {
      break;
    }

    /* find end of group/varname string (can be either EOS,
       space before next token, or following option string) */
    group_var_string = string;
    options_start = '{'; options_end = '}';
    while (*string)
    {
      if (! *string || isspace (*string) || *string == '{')
      {
        break;
      }

      /* check for an old-style options string (enclosed in square brackets) */
      if (*string == '[')
      {
        gindex = -1;
        /* find matching ']' */
        for (tmp = string+1; *tmp != ']' && *tmp; tmp++);
        if (*tmp == ']')
        {
          delimiter = *++tmp;
          *tmp = 0;
          gindex = CCTK_GroupIndexFromVar (group_var_string);
          *tmp = delimiter;
        }
        /* continue if the current token refers to a valid vector variable name
           otherwise assume the start of an old-style options string */
        if (gindex < 0 || CCTK_GroupData (gindex, &gdata) ||
            ! gdata.vectorgroup)
        {
          options_start = '['; options_end = ']';
          break;
        }
      }
      string++;
    }

    /* mark end of group/varname string */
    delimiter = *string;
    *string = 0;

    /* parse the option string if there is one */
    option_string = delimiter == options_start ? string + 1: NULL;
    if (option_string)
    {
      /* find end of option string (matching bracket) */
      nesting = 1;
      while (*(++string))
      {
        if (*string == options_start)
        {
          nesting++;
        }
        else if (*string == options_end)
        {
          if (--nesting == 0)
          {
            break;
          }
        }
      }
      delimiter = *string;
      *string = 0;
      if (option_string == group_var_string + 1)
      {
        CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                    "CCTK_TraverseString: option string '%s' in traversed "
                    "string '%s' is not associated with a group or variable "
                    "name", option_string, traverse_string);
        retval = -2;
        break;
      }
      else if (! (delimiter == options_end && nesting == 0))
      {
        CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                    "CCTK_TraverseString: unterminated option string '%s' "
                    "in traversed string '%s'", option_string, traverse_string);
        retval = -3;
        break;
      }
      else if (! (string[1] == 0 || isspace (string[1])))
      {
        CCTK_VWarn (2, __LINE__, __FILE__, "Cactus",
                    "CCTK_TraverseString: garbage at end of option string '%s' "
                    "in traversed string '%s'", option_string, traverse_string);
        retval = -4;
        break;
      }
      else if (options_start == '[')
      {
        CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                    "CCTK_TraverseString: enclosing options string '%s' in "
                    "square brackets is deprecated in BETA13. Please use "
                    "curly braces instead !", option_string);
      }
    }

#ifdef DEBUG_GROUPS
    printf ("group/varname is '%s', option string is '%s'\n",
            group_var_string, option_string ? option_string : "(null)");
#endif

    /* Look for the token 'all' */
    selected_all = CCTK_Equals (group_var_string, "all");

    /* See if this name is "<implementation>::<variable>" */
    if (! selected_all &&
        (selection == CCTK_VAR || selection == CCTK_GROUP_OR_VAR))
    {
      first = last = CCTK_VarIndex (group_var_string);
    }
    else
    {
      first = last = -1;
    }
    if (first < 0)
    {

      /* See if this name is "<implementation>::<group>" */
      if (! selected_all &&
          (selection == CCTK_GROUP || selection == CCTK_GROUP_OR_VAR))
      {
        gindex = CCTK_GroupIndex (group_var_string);
      }
      else
      {
        gindex = -1;
      }
      if (gindex >= 0)
      {
        /* We have a group so now need all the variables in the group */
        /* Note: CCTK_FirstVarIndexI is negative if there are zero
           variables in the group */
        first = CCTK_FirstVarIndexI (gindex);
        last = first + CCTK_NumVarsInGroupI (gindex) - 1;
      }
      else if (selected_all)
      {
        first = 0;
        if (selection == CCTK_GROUP)
        {
          last = CCTK_NumGroups () - 1;
        }
        else
        {
          last = CCTK_NumVars () - 1;
        }
      }
      else
      {
        first = last = -1;
      }
    }

    /* invoke the callback */
    if (first >= 0)
    {
      for (vindex = first; vindex <= last; vindex++)
      {
        (*callback) (vindex, option_string, callback_arg);
      }
      if (retval >= 0) retval += last - first + 1;
    }
    /* Only emit an error message if the name is really invalid.  If
       it is a valid group name, but the group has zero variables,
       then it will be first<0 and gindex>=0. */
    else if (gindex < 0)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "CCTK_TraverseString: invalid group/variable name '%s' in "
                  "traversed string '%s'", group_var_string, traverse_string);
      retval = -5;
    }

    /* advance the parse string pointer */
    if (delimiter)
    {
      string++;
    }
  } /* end of while loop over all tokens in parse string */

  /* clean up */
  free (parse_string);

  return (retval);
}


#if 0
/*@@
   @routine    CCTKi_PrintGroupInfo
   @date       Thu Jan 14 15:25:54 1999
   @author     Gerd Lanfermann
   @desc
     Debugging info on the Groups.
   @enddesc
@@*/

void CCTKi_PrintGroupInfo (void)
{
  int group_num;

  for (group_num = 0; group_num < n_groups; group_num++)
  {
    printf ("GROUP INFO: GrpNo./imp_name/name %d   >%s<   >%s<\n",
           group_num,
           groups[group_num].implementation,
           groups[group_num].name);
  }
}
#endif


 /*@@
   @routine    CCTKi_CreateGroup
   @date       Thu Jan 14 15:25:54 1999
   @author     Tom Goodale
   @desc
               Creates a new CCTK group
   @enddesc

   @returntype int
   @returndesc
               0 for success, non-zero otherwise
   @endreturndesc
@@*/
int CCTKi_CreateGroup (const char *gname,
                       const char *thorn,
                       const char *imp,
                       const char *gtype,
                       const char *vtype,
                       const char *gscope,
                       int         dimension,
                       int         ntimelevels,
                       const char *dtype,
                       const char *size,
                       const char *ghostsize,
                       const char *tags,
                       const char *vararraysize,
                       int         n_basevars,
                       ...
                       )
{
  int retval;
  int groupscope;
  int variable;

  int vectorlength;
  int elem;

  va_list ap;

  char *variable_name;

  cGroupDefinition *group;

  vectorlength = 1;

  group = NULL;
  variable_name = NULL;

  retval = 0;

  va_start (ap, n_basevars);

  if (vararraysize)
  {
    vararraysize = Util_Strdup(vararraysize);
    assert (vararraysize);

    vectorlength = CCTKi_ParamExpressionToInt(vararraysize,thorn);

    if(vectorlength < 0)
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "CCTKi_CreateGroup: length of group %s less than 0 !",
                  gname);
    }
  }

  /* Allocate storage for the group */
  groupscope = CCTK_GroupScopeNumber (gscope);
  if (groupscope == CCTK_PUBLIC || groupscope == CCTK_PROTECTED)
  {
    group = CCTKi_SetupGroup (imp, gname,
                              n_basevars * vectorlength, vectorlength);
  }
  else if (groupscope == CCTK_PRIVATE)
  {
    group = CCTKi_SetupGroup (thorn, gname,
                              n_basevars * vectorlength, vectorlength);
  }
  else
  {
    CCTK_Warn (1, __LINE__, __FILE__, "Cactus",
              "CCTKi_CreateGroup: Unrecognised group scope");
  }

  /* Allocate storage for the group and setup some stuff. */
  if (group)
  {
    group->dim          = dimension;
    group->gtype        = CCTK_GroupTypeNumber (gtype);
    group->vtype        = CCTK_VarTypeNumber (vtype);
    group->gscope       = groupscope;
    group->dtype        = CCTK_GroupDistribNumber (dtype);
    group->n_timelevels = ntimelevels;
    group->tags_string  = Util_Strdup(tags);

    group->tags_table   = Util_TableCreateFromString(tags);
    if(group->tags_table < 0)
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "CCTKi_CreateGroup: Failed to create TAGS table for group '%s' from thorn '%s'",
                  gname, thorn);
    }

    /* Extract the variable names from the argument list. */

    group->vararraysize = vararraysize;

    for (variable = 0; variable < n_basevars; variable++)
    {
      variable_name = va_arg (ap, char *);

      if (!vararraysize)
      {
        assert (vectorlength == 1);
        group->variables[variable].name = Util_Strdup(variable_name);
      }
      else
      {
        for (elem = 0; elem < vectorlength; elem++)
        {
          char *name = NULL;
          Util_asprintf(&name, "%s[%d]", variable_name, elem);

          group->variables[variable * vectorlength + elem].name = name;
        }
      }
    }

    va_end (ap);

    if (dimension > maxdim)
    {
      maxdim    = dimension;
    }
    group->size      = CCTKi_ExtractSize (dimension, thorn, size, gname);
    group->ghostsize = CCTKi_ExtractSize (dimension, thorn, ghostsize, gname);

    /* Only typically have GFs in a single dimension */
    if (group->gtype == CCTK_GF)
    {
      if (gfdim > 0)
      {
        if (group->dim != gfdim)
        {
          retval = 1;
        }
      }
      else
      {
        gfdim = group->dim;
      }
    }

  }
  else
  {
    retval = 2;
  }

  if (retval)
  {
    CCTK_Warn (4, __LINE__, __FILE__, "Cactus", "CCTKi_CreateGroup: Error");
  }

  return (retval);
}


 /*@@
   @routine    CCTK_GroupImpI
   @date       20 Oct 2001
   @author     Gabrielle Allen
   @desc
   Return the implementation which created a group
   @enddesc

   @returntype const char *
   @returndesc
   Thorn name
   @endreturndesc
@@*/

const char *CCTK_GroupImplementationI(int group)
{
  const char *imp;

  imp = groups[group].implementation;

  return imp;
}

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

 /*@@
   @routine    CCTKi_SetupGroup
   @date       Thu Jan 14 16:38:40 1999
   @author     Tom Goodale
   @desc
               Stores the data associated with a group.
   @enddesc

   @returntype cGroupDefinition *
   @returndesc
               pointer to the associated group data structure, or
               NULL if out of memory
   @endreturndesc
@@*/
static cGroupDefinition *CCTKi_SetupGroup (const char *implementation,
                                           const char *name,
                                           int n_variables,
                                           int vectorlength)
{
  int *temp_int;
  void *temp;
  cGroupDefinition *returndata;
  int variable;
  int group;


  for (group = 0; group < n_groups; group++)
  {
    if (CCTK_Equals (implementation, groups[group].implementation) &&
        CCTK_Equals (name, groups[group].name))
    {
      break;
    }
  }

  if (group >= n_groups)
  {
    /* Resize the array of groups */
    temp = realloc (groups, (n_groups + 1) * sizeof (cGroupDefinition));
    if (temp)
    {
      groups = (cGroupDefinition *) temp;

      /* Allocate memory to various fields */
      groups[n_groups].implementation = malloc (strlen (implementation)+1);

      groups[n_groups].name = malloc (strlen (name) + 1);

      groups[n_groups].variables = malloc (n_variables *
                                           sizeof (cVariableDefinition));

      /* Resize the array holding correspondence between vars and groups. */
      temp_int = realloc (group_of_variable,
                          (total_variables+n_variables) * sizeof (int));

      if ((n_variables==0 || (groups[n_groups].implementation &&
                              groups[n_groups].name &&
                              groups[n_groups].variables))
          &&
          ((total_variables+n_variables==0) || (temp_int)))
      {
        /* Fill in the data structures. */
        group_of_variable = temp_int;

        strcpy (groups[n_groups].implementation, implementation);
        strcpy (groups[n_groups].name, name);

        groups[n_groups].number       = n_groups;
        groups[n_groups].n_variables  = n_variables;
        groups[n_groups].vectorlength = vectorlength;

        /* Fill in global variable numbers. */
        for (variable = 0; variable < n_variables; variable++)
        {
          groups[n_groups].variables[variable].number = total_variables;

          group_of_variable[total_variables] = n_groups;

          total_variables++;
        }

        n_groups++;
      }
      else
      {
        /* Memory allocation failed, so free any which may have been allocated. */
        free (groups[n_groups].implementation);
        groups[n_groups].implementation = NULL;

        free (groups[n_groups].name);
        groups[n_groups].name = NULL;

        free (groups[n_groups].variables);
        groups[n_groups].variables = NULL;

      }
    }

    /* Return the new group definition structure if successful, otherwise NULL.*/
    if (temp && groups[n_groups-1].name)
    {
      returndata =  &groups[n_groups-1];
    }
    else
    {
      returndata = NULL;
    }
  }
  else
  {
    returndata = &groups[group];
  }

#ifdef DEBUG_GROUPS
  printf ("Setting up group %s::%s\n", implementation, name);
#endif

  return (returndata);
}


 /*@@
   @routine    CCTKi_ExtractSize
   @date       Sun Nov 28 12:38:38 1999
   @author     Tom Goodale
   @desc
               Extracts the size array from a comma-separated list of
               positive integer constants or arithmetical combinations of
               integer parameters.

               FIXME:  This routine originally returned a pointer to the
                       parameter values, which is why it returns **.
                       Now it doesn't, so we probably have a memory leak.
   @enddesc

   @returntype CCTK_INT **
   @returndesc
               pointer to an allocated array of sizes for the given list, or
               NULL if no parameter list is given or out of memory
   @endreturndesc
@@*/
static CCTK_INT **CCTKi_ExtractSize (int dimension,
                                     const char *this_thorn,
                                     const char *sizestring,
                                     const char *gname)
{
  int         dim;
  char       *tmp;
  const char *last_comma, *next_comma;
  CCTK_INT   **size_array;

  if (dimension < 0)
  {
    CCTK_Warn (0, __LINE__, __FILE__, "Cactus","Illegal dimension specified");
  }

  if (strlen (sizestring))
  {
    next_comma = sizestring;

    size_array = malloc (dimension * sizeof (CCTK_INT *));

    if (size_array)
    {
      if (dimension > 0)
      {
        size_array[0] = malloc (dimension * sizeof (CCTK_INT));

        for (dim = 1; dim < dimension; dim++)
        {
          size_array[dim] = size_array[0] + dim;
        }

        for (dim = 0; dim < dimension; dim++)
        {
          if (!next_comma)
          {
             CCTK_VWarn (0, __LINE__, __FILE__, "Cactus","Insufficient dimension size specified in %s::%s -> %s", this_thorn, gname, sizestring);
          }
          /* find the comma as a delimiter for different dimension sizes */
          last_comma = next_comma[0] == ',' ? next_comma+1 : next_comma;
          next_comma = strstr (last_comma, ",");

          /* copy dimension size token into a work string buffer */
          tmp = strdup (last_comma);
          if (next_comma)
          {
            tmp[next_comma-last_comma] = '\0';
          }

          *size_array[dim] = CCTKi_ParamExpressionToInt (tmp, this_thorn);

          free (tmp);
        }
      }
    }
  }
  else
  {
    /* No size specified */
    size_array = NULL;
  }

  return size_array;
}

 /*@@
   @routine    CCTKi_GroupLengthAsPointer
   @date       Sun Oct  7 03:58:44 2001
   @author     Tom Goodale
   @desc
               Get the number of variables in a group,
               or the number of elements in a vector group
   @enddesc

   @var     fullgroupname
   @vdesc   The full name of a GV group
   @vtype   const char *
   @vio     in
   @endvar

   @returntype const int *
   @returndesc
               pointer to an integer containing the number of variables in the group
               NULL if group doesn't exist
   @endreturndesc
@@*/
const int *CCTKi_GroupLengthAsPointer(const char *fullgroupname)
{
  int group;
  const int *retval;
  char *impname, *groupname;


  retval = NULL;
  impname = groupname = NULL;

  if (! CCTK_DecomposeName (fullgroupname, &impname, &groupname))
  {
    for (group = 0; group < n_groups; group++)
    {
      if (CCTK_Equals (impname, groups[group].implementation) &&
          CCTK_Equals (groupname, groups[group].name))
      {
        retval = groups[group].vectorlength ?
                 &groups[group].vectorlength : &groups[group].n_variables;
        break;
      }
    }

    if (retval == NULL)
    {
      CCTK_VWarn (6, __LINE__, __FILE__, "Cactus",
                  "CCTKi_GroupLengthAsPointer: No group named '%s' found",
                  fullgroupname);
    }
  }

  /* Free memory from CCTK_DecomposeName */
  free (impname);
  free (groupname);

  return retval;
}

 /*@@
   @routine    IntParameterEvaluator
   @date       Fri Oct 12 10:01:32 2001
   @author     Tom Goodale
   @desc
   Evaluates integer parameters for the expression parser.
   @enddesc
   @var     nvars
   @vdesc   Number of variables to evaluate
   @vtype   int
   @vio     in
   @endvar
   @var     vars
   @vdesc   an array of parameter names or integers
   @vtype   const char * const *
   @vio     in
   @endvar
   @var     vals
   @vdesc   Output array to hold values
   @vtype   uExpressionValue *
   @vio     out
   @endvar
   @var     data
   @vdesc   Data passed from expression evaluator
   @vtype   const void *
   @vio     in
   @vcomment
     Should be a char * with the thorn name.
   @endvar

   @returntype int
   @returndesc
   0
   @endreturndesc
@@*/
static int IntParameterEvaluator(int nvars,
                                 const char * const *vars,
                                 uExpressionValue *vals,
                                 const void *data)
{
  int i;

  char *tmp;
  char *endptr;
  const CCTK_INT *paramval;
  char *thorn;
  char *param;
  const char *use_thorn;
  const char *use_param;
  int type;

  for(i=0; i < nvars; i++)
  {
    vals[i].type = ival;
    tmp = Util_Strdup(vars[i]);

    vals[i].value.ival = strtol(tmp,&endptr,0);

    if(endptr == tmp)
    {
      if(CCTK_DecomposeName (tmp,&thorn,&param))
      {
        thorn = NULL;
        param = NULL;
        use_thorn = (const char *)data;
        use_param = tmp;
      }
      else
      {
        /* FIXME: need to be a bit more careful here. */
        if (CCTK_IsImplementationActive (thorn))
        {
          use_thorn = CCTK_ActivatingThorn (thorn);
        }
        else
        {
          use_thorn = thorn;
        }
        use_param = param;
      }

      paramval = (const CCTK_INT *) CCTK_ParameterGet (use_param, use_thorn, &type);

      if(paramval && type==PARAMETER_INTEGER)
      {
        vals[i].value.ival = *paramval;
      }
      else if(!paramval)
      {
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "IntParameterEvaluator: '%s::%s' is not a parameter",
                    use_thorn, use_param);
      }
      else
      {
        CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                    "IntParameterEvaluator: '%s::%s' is not an integer parameter",
                    use_thorn, use_param);
      }
      free(thorn);
      free(param);
    }

    free(tmp);
  }

  return 0;
}

 /*@@
   @routine    CCTKi_ParamExpressionToInt
   @date       Fri Oct 12 10:04:18 2001
   @author     Tom Goodale
   @desc
   Parses an arithmetic expression involving integer parameter names and
   returns the final integer
   @enddesc
   @var     expression
   @vdesc   The expression to be parsed
   @vtype   const char *
   @vio     in
   @endvar
   @var     thorn
   @vdesc   The thorn name
   @vtype   const char *
   @vio     in
   @vcomment
     Used for parameters which aren't fully qualified.
   @endvar

   @returntype int
   @returndesc
     The final value of the exppression
   @endreturndesc
@@*/
static int CCTKi_ParamExpressionToInt(const char *expression, const char *thorn)
{
  int retval;
  uExpression parsed_expression;
  uExpressionValue value;
  char *this_thorn;

  this_thorn = Util_Strdup(thorn);

  parsed_expression = Util_ExpressionParse(expression);

  if(parsed_expression)
  {
    /* Evaluate the expression */
    retval = Util_ExpressionEvaluate(parsed_expression,
                                     &value,
                                     IntParameterEvaluator,
                                     this_thorn);

    Util_ExpressionFree(parsed_expression);
  }
  else
  {
    CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                "Unable to parse parameter expression '%s'",
                expression);
    retval = -1;
  }

  free(this_thorn);

  if(retval == 0)
  {
    retval = value.value.ival;
  }
  else
  {
    CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                "Unable to evaluate parameter expression '%s'",
                expression);
  }

  return retval;
}

