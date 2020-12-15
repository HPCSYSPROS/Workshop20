 /*@@
   @header    cctk_Groups.h
   @date      Mon Feb  8 14:47:10 1999
   @author    Tom Goodale
   @desc 
   Prototypes and constants for group functions.
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CCTK_GROUPS_H_
#define _CCTK_GROUPS_H_ 1

typedef struct
{
  int grouptype;
  int vartype;
  int disttype;
  int dim;
  int numvars;
  int numtimelevels;
  int vectorgroup;
  int vectorlength;
  int tagstable;
} cGroup;

/* Prototypes */

#include "cctk_Types.h"

#ifdef __cplusplus 
extern "C" 
{
#endif

int         CCTK_DecomposeName(const char *fullname, 
                               char **implementation, 
                               char **name);

int         CCTK_FirstVarIndex(const char *group);
int         CCTK_FirstVarIndexI(int group);
char      * CCTK_FullName(int var);

int         CCTK_GroupData(int group, cGroup *gp);
int         CCTK_GroupDimI(int group);
int         CCTK_GroupDimFromVarI(int vi);
int         CCTK_GroupDistribNumber(const char *dtype);
CCTK_INT ** CCTK_GroupGhostsizesI(int group);
const char *CCTK_GroupImplementationI(int group);
int         CCTK_GroupIndex(const char *groupname);
int         CCTK_GroupIndexFromVar(const char *var);
int         CCTK_GroupIndexFromVarI(int var);
char      * CCTK_GroupName(int groupnum);
char      * CCTK_GroupNameFromVarI(int varnum);
int         CCTK_GroupScopeNumber(const char *type);
CCTK_INT ** CCTK_GroupSizesI(int group);
int         CCTK_GroupTypeFromVarI(int var);
int         CCTK_GroupTypeNumber(const char *type);
int         CCTK_GroupTypeI(int group);

const char *CCTK_ImpFromVarI(int var);

int         CCTK_MaxDim(void);
int         CCTK_MaxGFDim(void);

int         CCTK_NumGroups(void);

#define HAVE_CCTK_DECLARED_TIMELEVELS
int         CCTK_DeclaredTimeLevels(const char *group);
int         CCTK_DeclaredTimeLevelsVN(const char *var);
int         CCTK_DeclaredTimeLevelsVI(int var);
int         CCTK_DeclaredTimeLevelsGN(const char *group);
int         CCTK_DeclaredTimeLevelsGI(int group);


int         CCTK_NumVars(void);
int         CCTK_NumVarsInGroup(const char *group);
int         CCTK_NumVarsInGroupI(int group);

int         CCTK_VarIndex(const char *variablename);
const char *CCTK_VarName(int varnum);
int         CCTK_VarTypeI(int var);
int         CCTK_VarTypeNumber(const char *type);
const char *CCTK_VarTypeName(int vartype);

int         CCTK_VarTypeSize(int vtype);

const int * CCTKi_GroupLengthAsPointer(const char *fullgroupname);

/* traverse a string of group and/or variable names */
int CCTK_TraverseString (const char *parsestring,
                         void (*callback) (int index,
                                           const char *optstring,
                                           void *callback_arg),
                         void *callback_arg,
                         int selection);

int CCTK_GroupTagsTable(const char *group);
int CCTK_GroupTagsTableI(int group);

#ifdef __cplusplus 
}
#endif


#endif /* _CCTK_GROUPS_H_ */
