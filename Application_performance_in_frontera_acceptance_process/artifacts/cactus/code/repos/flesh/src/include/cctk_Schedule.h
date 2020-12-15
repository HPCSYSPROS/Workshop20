 /*@@
   @header  cctk_schedule.h
   @date    Thu Sep 16 19:05:27 1999
   @author  Tom Goodale
   @desc
            Routines for creating schedule stuff.
   @enddesc
   @version $Header$
 @@*/

#ifndef _CCTK_SCHEDULE_H_
#define _CCTK_SCHEDULE_H_

#include <stdio.h>

typedef enum {LangNone, LangC, LangFortran} cLanguage;

typedef enum {FunctionNoArgs, FunctionOneArg, FunctionStandard} cFunctionType;

typedef struct
{
  cLanguage language;

  int (*FortranCaller)(cGH *, void *);

  cFunctionType type;

  int n_SyncGroups;
  int *SyncGroups;

  /* Option Flags */

  int meta;
  int meta_early;
  int meta_late;
  int global;
  int global_early;
  int global_late;
  int level;
  int singlemap;
  int local;

  int loop_meta;
  int loop_global;
  int loop_level;
  int loop_singlemap;
  int loop_local;

  int tags;

  /* The last items should be considered volatile and may not stay here */

  int n_TriggerGroups;
  int *TriggerGroups;

#define CACTUS_HAS_READS_CLAUSES
  int n_WritesClauses;
  const char **WritesClauses;
  int n_ReadsClauses;
  const char **ReadsClauses;

  char *where;
  char *routine;
  char *thorn;

} cFunctionData;

#ifdef __cplusplus
extern "C"
{
#endif

int CCTK_CallFunction(void *function,
                      cFunctionData *fdata,
                      void *data);

const cFunctionData *CCTK_ScheduleQueryCurrentFunction(const cGH *GH);

int CCTK_ScheduleTraverse(const char *where,
                          void *GH,
                          int (*CallFunction)(void *, cFunctionData *, void *));

int CCTK_SchedulePrint(const char *where);
int CCTK_SchedulePrintTimes(const char *where);
int CCTK_SchedulePrintTimesToFile(const char *where, FILE *file);

cLanguage CCTK_TranslateLanguage(const char *sval);

#if 0
int CCTK_ScheduleFunction(void *function,
                          const char *name,
                          const char *thorn,
                          const char *implementation,
                          const char *description,
                          const char *where,
                          const char *language,
                          int n_mem_groups,
                          int n_comm_groups,
                          int n_trigger_groups,
                          int n_before,
                          int n_after,
                          int n_writes,
                          int n_reads,
                          int n_while,
                          int n_if,
                          ...);

int CCTK_ScheduleGroup(const char *name,
                       const char *thorn,
                       const char *implementation,
                       const char *description,
                       const char *where,
                       int n_mem_groups,
                       int n_comm_groups,
                       int n_trigger_groups,
                       int n_before,
                       int n_after,
                       int n_while,
                       int n_if,
                       ...);

int CCTK_ScheduleGroupStorage(const char *group);

int CCTK_ScheduleGroupComm(const char *group);

int CCTK_ScheduleTraverse(const char *where,
                          void *GH,
                          int (*calling_function)(void *, void *, void *));

int CCTK_ScheduleGHInit(void *GH);
#endif


#ifdef __cplusplus
}
#endif

#endif /*_CCTK_SCHEDULE_H_*/
