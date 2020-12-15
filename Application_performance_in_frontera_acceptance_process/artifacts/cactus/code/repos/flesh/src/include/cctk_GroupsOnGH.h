 /*@@
   @header    cctk_GroupsOnGH.h
   @date      Wed Apr 7 1999
   @author    Gabrielle Allen
   @desc 
   Prototypes and constants for group functions which use GH stucture.
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CCTK_GROUPSONGH_H_
#define _CCTK_GROUPSONGH_H_

typedef struct GROUPDYNAMICDATA
{
  int dim;
  const int *gsh;
  const int *lsh;
  const int *lbnd;
  const int *ubnd;
  const int *ash;
  const int *bbox;
  const int *nghostzones;
  int maxtimelevels;
  int activetimelevels;
} cGroupDynamicData;

/* Prototypes */

#ifdef __cplusplus 
extern "C" {
#endif

void *CCTK_VarDataPtr(const cGH *GH, int timelevel, const char *fullvarname);
void *CCTK_VarDataPtrI(const cGH *GH, int timelevel, int varindex);
void *CCTK_VarDataPtrB(const cGH *GH, int timelevel, int varindex, char *fullvarname);

int CCTK_DisableGroupStorageI(const cGH *GH, int group);
int CCTK_DisableGroupCommI(const cGH *GH, int group);
int CCTK_EnableGroupStorageI(const cGH *GH, int group);
int CCTK_EnableGroupCommI(const cGH *GH, int group);

int CCTK_GrouplbndGN(const cGH *GH, int dim, int *lbnd, const char *groupname);
int CCTK_GrouplbndVN(const cGH *GH, int dim, int *lbnd, const char *varname);
int CCTK_GrouplbndGI(const cGH *GH, int dim, int *lbnd, int groupindex);
int CCTK_GrouplbndVI(const cGH *GH, int dim, int *lbnd, int varindex);

int CCTK_GroupubndGN(const cGH *GH, int dim, int *ubnd, const char *groupname);
int CCTK_GroupubndVN(const cGH *GH, int dim, int *ubnd, const char *varname);
int CCTK_GroupubndGI(const cGH *GH, int dim, int *ubnd, int groupindex);
int CCTK_GroupubndVI(const cGH *GH, int dim, int *ubnd, int varindex);

int CCTK_GrouplshGN(const cGH *GH, int dim, int *lsh, const char *groupname);
int CCTK_GrouplshVN(const cGH *GH, int dim, int *lsh, const char *varname);
int CCTK_GrouplshGI(const cGH *GH, int dim, int *lsh, int groupindex);
int CCTK_GrouplshVI(const cGH *GH, int dim, int *lsh, int varindex);

int CCTK_GroupashGN(const cGH *GH, int dim, int *ash, const char *groupname);
int CCTK_GroupashVN(const cGH *GH, int dim, int *ash, const char *varname);
int CCTK_GroupashGI(const cGH *GH, int dim, int *ash, int groupindex);
int CCTK_GroupashVI(const cGH *GH, int dim, int *ash, int varindex);

int CCTK_GroupgshGN(const cGH *GH, int dim, int *gsh, const char *groupname);
int CCTK_GroupgshVN(const cGH *GH, int dim, int *gsh, const char *varname);
int CCTK_GroupgshGI(const cGH *GH, int dim, int *gsh, int groupindex);
int CCTK_GroupgshVI(const cGH *GH, int dim, int *gsh, int varindex);

int CCTK_GroupbboxGI(const cGH *cctkGH, int size, int *bbox, int groupindex);
int CCTK_GroupbboxGN(const cGH *cctkGH, int size, int *bbox, const char *groupname);
int CCTK_GroupbboxVI(const cGH *cctkGH, int size, int *bbox, int varindex);
int CCTK_GroupbboxVN(const cGH *cctkGH, int size, int *bbox, const char *varname);

int CCTK_GroupnghostzonesGN(const cGH *GH, int dim, int *nghostzones, const char *groupname);
int CCTK_GroupnghostzonesVN(const cGH *GH, int dim, int *nghostzones, const char *varname);
int CCTK_GroupnghostzonesGI(const cGH *GH, int dim, int *nghostzones, int groupindex);
int CCTK_GroupnghostzonesVI(const cGH *GH, int dim, int *nghostzones, int varindex);

int CCTK_ActiveTimeLevels(const cGH *GH, const char *groupname);
int CCTK_ActiveTimeLevelsGN(const cGH *GH, const char *groupname);
int CCTK_ActiveTimeLevelsGI(const cGH *GH, int groupindex);
int CCTK_ActiveTimeLevelsVN(const cGH *GH, const char *varname);
int CCTK_ActiveTimeLevelsVI(const cGH *GH, int varindex);

int CCTK_MaxActiveTimeLevels(const cGH *GH, const char *groupname);
int CCTK_MaxActiveTimeLevelsGN(const cGH *GH, const char *groupname);
int CCTK_MaxActiveTimeLevelsGI(const cGH *GH, int groupindex);
int CCTK_MaxActiveTimeLevelsVN(const cGH *GH, const char *varname);
int CCTK_MaxActiveTimeLevelsVI(const cGH *GH, int varindex);

#ifdef __cplusplus 
}
#endif
 
#endif
