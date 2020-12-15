/*@@
   @file      ScheduleInterface.c
   @date      Thu Sep 16 14:06:21 1999
   @author    Tom Goodale
   @desc
              Routines to interface the main part of Cactus to the schedular.
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "cctk_Flesh.h"
#include "cctk_WarnLevel.h"
#include "cctk_Misc.h"

#include "cctk_Schedule.h"
#include "cctki_ScheduleBindings.h"
#include "cctki_Schedule.h"

#include "cctk_Comm.h"
#include "cctk_Sync.h"

#include "cctk_Constants.h"
#include "cctk_Groups.h"
#include "cctk_GroupsOnGH.h"

#include "cctki_FortranWrappers.h"

#include "cctk_Timers.h"

#include "util_Table.h"
#include <mpi.h>

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_ScheduleInterface_c);


/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

typedef enum {sched_none, sched_group, sched_function} iSchedType;
typedef enum {schedpoint_misc, schedpoint_analysis} iSchedPoint;

/* Timer data for a routine in a given schedule bin */
typedef struct t_timer
{
  struct t_timer *next;
  int timer_handle;
  char *schedule_bin;
  int has_been_output;
} t_timer;

typedef struct
{
  /* Static data */
  char *description;

  /*char *thorn; MOVED TO FunctionData */
  char *implementation;

  iSchedType type;

  cFunctionData FunctionData;

  int n_mem_groups;
  int *mem_groups;
  int *timelevels;

  int n_comm_groups;
  int *comm_groups;

  /* Timer data */
  t_timer *timers;

  /* Dynamic data */
  int *CommOnEntry;
  int *StorageOnEntry;

  int done_entry;
  int synchronised;

} t_attribute;

typedef struct
{
  cGH *GH;
  iSchedPoint schedpoint;
  const char *schedule_bin;
  unsigned int n_functions;

  cTimerData *info;
  cTimerData *total_time;
  int print_headers;
  FILE *file;                   /* output file */

  /* Stuff passed in in user calls */

  int (*CallFunction)(void *, cFunctionData *, void *);

} t_sched_data;


/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int ScheduleTraverse(const char *where,
                            void *GH,
                            int (*CallFunction)(void *, cFunctionData *, void *));

static t_attribute *CreateAttribute(const char *where,
                                    const char *name,
                                    const char *description,
                                    const char *language,
                                    const char *thorn,
                                    const char *implementation,
                                    int n_mem_groups,
                                    int n_comm_groups,
                                    int n_trigger_groups,
                                    int n_sync_groups,
                                    int n_writes,
                                    int n_reads,
                                    int n_options,
                                    const int *timelevels,
                                    va_list *ap);

static int ParseOptionList(int n_items,
                           t_attribute *attribute,
                           va_list *ap);

static int ParseTagsTable(t_attribute *attribute,
                          va_list *ap);

static int InitialiseOptionList(t_attribute *attribute);

static int ParseOption(t_attribute *attribute,
                       const char *option);

static t_sched_modifier *CreateModifiers(int n_before,
                                         int n_after,
                                         int n_while,
                                         int n_if,
                                         va_list *ap);

int ValidateModifiers(t_sched_modifier *modifier);

static int CreateGroupIndexList(int n_items, int *array, va_list *ap);
static int CreateStringList(int n_items, const char **array, va_list *ap);
static t_sched_modifier *CreateTypedModifier(t_sched_modifier *modifier,
                                             const char *type,
                                             int n_items,
                                             va_list *ap);
static cFunctionType TranslateFunctionType(const char *where);

static int SchedulePrint(const char *where);

static int CCTKi_SchedulePrintEntry(t_attribute *attribute, t_sched_data *data);
static int CCTKi_SchedulePrintExit(t_attribute *attribute, t_sched_data *data);
static int CCTKi_SchedulePrintWhile(int n_whiles,
                                    char **whiles,
                                    t_attribute *attribute,
                                    t_sched_data *data,
                                    int first);
static int CCTKi_SchedulePrintIf(int n_if,
                                 char **ifs,
                                 t_attribute *attribute,
                                 t_sched_data *data,
                                 int first);
static int CCTKi_SchedulePrintFunction(void *function, t_attribute *attribute, t_sched_data *data);

static int CCTKi_ScheduleCallEntry(t_attribute *attribute, t_sched_data *data);
static int CCTKi_ScheduleCallExit(t_attribute *attribute, t_sched_data *data);
static int CCTKi_ScheduleCallWhile(int n_whiles,
                                   char **whiles,
                                   t_attribute *attribute,
                                   t_sched_data *data,
                                   int first);
static int CCTKi_ScheduleCallIf(int n_ifs,
                                char **ifs,
                                t_attribute *attribute,
                                t_sched_data *data);
static int CCTKi_ScheduleCallFunction(void *function,
                                      t_attribute *attribute,
                                      t_sched_data *data);

static int SchedulePrintTimes(const char *where, t_sched_data *data);
static int CCTKi_SchedulePrintTimesFunction(void *function,
                                            t_attribute *attribute,
                                            t_sched_data *data);
static void CCTKi_SchedulePrintTimerInfo(cTimerData *info,
                                         cTimerData *total_time,
                                         const char *where,
                                         const char *description,
                                         FILE *file);
static void CCTKi_SchedulePrintTimerHeaders(cTimerData *info,
                                            FILE *file);
static int CCTKi_ScheduleResetTimerOutputFlag(void *function,
                                              t_attribute *attribute,
                                              t_sched_data *data);
static void PrintDelimiterLine (char delimiter,
                                const cTimerData *timer,
                                FILE *file);


/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/* FIXME: these should be put in a header somewhere */

int CCTKi_TriggerSaysGo(cGH *GH, int variable);
int CCTKi_TriggerAction(void *GH, int variable);


/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/


static int indent_level = 0;

static int n_scheduled_comm_groups = 0;
static int *scheduled_comm_groups = NULL;

static int n_scheduled_storage_groups = 0;
static int *scheduled_storage_groups = NULL;
static int *scheduled_storage_groups_timelevels = NULL;

static cTimerData *timerinfo = NULL;
static int total_timer = -1;

static const cFunctionData *current_scheduled_function = NULL;


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTK_CallFunction
   @date       Thu Jan 27 11:29:47 2000
   @author     Tom Goodale
   @desc
   Calls a function depending upon the data passed in the the
   fdata structure.
   @enddesc
   @calls

   @var     function
   @vdesc   pointer to function
   @vtype   void *
   @vio     in
   @endvar
   @var     fdata
   @vdesc   data about the function
   @vtype   cFunctionData *
   @vio     in
   @endvar
   @var     data
   @vdesc   Data to be passed to the function
   @vtype   void *
   @vio     inout
   @endvar

   @returntype int
   @returndesc
   0 - didn't synchronise
   @endreturndesc
@@*/
int CCTK_CallFunction(void *function,
                      cFunctionData *fdata,
                      void *data)
{
  void (*standardfunc)(void *);

  int (*noargsfunc)(void);

  int (*oneargfunc)(void *);

  if(current_scheduled_function != NULL)
  {
    CCTK_VWarn(CCTK_WARN_PICKY, __LINE__, __FILE__, "Cactus",
               "CCTK_CallFunction: recursive call, calling "
               "'%s: %s::%s' while within '%s: %s::%s'",
               fdata->where, fdata->thorn, fdata->routine,
               current_scheduled_function->where,
               current_scheduled_function->thorn,
               current_scheduled_function->routine);
  }
  current_scheduled_function = fdata;
  double start = MPI_Wtime();

  switch(fdata->type)
  {
    case FunctionNoArgs:
      noargsfunc = (int (*)(void))function;
      noargsfunc();
      break;
    case FunctionOneArg:
      oneargfunc = (int (*)(void *))function;
      oneargfunc(data);
      break;
    case FunctionStandard:
      switch(fdata->language)
      {
        case LangC:
          standardfunc = (void (*)(void *))function;
          standardfunc(data);
          break;
        case LangFortran:
          fdata->FortranCaller(data, function);
          break;
        case LangNone:
          /* this should never happen */
          /* fall through */
        default :
          CCTK_Warn(1,__LINE__,__FILE__,"Cactus",
                    "CCTK_CallFunction: Unknown language.");
      }
      break;
    default :
      CCTK_Warn(1,__LINE__,__FILE__,"Cactus",
                "CCTK_CallFunction: Unknown function type.");
  }

  if(CCTK_MyProc(NULL) == 0) {
    //fprintf(stderr, "%s::%s %g %d\n", fdata->thorn, fdata->routine, MPI_Wtime() - start, CCTK_RunTime());
  }

  current_scheduled_function = NULL;

  /* Return 0, meaning didn't synchronise */
  return 0;
}

 /*@@
   @routine    CCTK_ScheduleQueryCurrentFunction
   @date       Fri Apr 20 08:57:49 PDT 2012
   @author     Roland Haas
   @desc
   Returns the cFunctionData of the function currenlty executing via
   CCTK_CallFunction.
   @enddesc
   @calls

   @var     GH
   @vdesc   GH data
   @vtype   const cGH *
   @vio     in
   @endvar

   @returntype cFunctionData *
   @returndesc
   Data about the function. NULL if no function is currently executing via
   CCTK_CallFunction.
   @endreturndesc
@@*/
const cFunctionData *CCTK_ScheduleQueryCurrentFunction(const cGH * CCTK_ATTRIBUTE_UNUSED GH)
{
  return current_scheduled_function;
}

/*@@
   @routine    CCTKi_ScheduleFunction
   @date       Thu Sep 16 18:19:01 1999
   @author     Tom Goodale
   @desc
   Schedules a function.
   @enddesc
   @calls

   @var     function
   @vdesc   function to be scheduled
   @vtype   void *
   @vio     in
   @endvar
   @var     name
   @vdesc   working name of function to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     thorn
   @vdesc   name of thorn providing function to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     implementation
   @vdesc   name of implementation thorn belongs to
   @vtype   const char *
   @vio     in
   @endvar
   @var     description
   @vdesc   desciption of function to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     where
   @vdesc   where to schedule the function
   @vtype   const char *
   @vio     in
   @endvar
   @var     language
   @vdesc   language of function to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     n_mem_groups
   @vdesc   Number of groups needing memory switched on during this function
   @vtype   int
   @vio     in
   @endvar
   @var     n_comm_groups
   @vdesc   Number of groups needing communication switched on during this function
   @vtype   int
   @vio     in
   @endvar
   @var     n_trigger_groups
   @vdesc   Number of groups to trigger this function on
   @vtype   int
   @vio     in
   @endvar
   @var     n_sync_groups
   @vdesc   Number of groups needing synchronisation after this function
   @vtype   int
   @vio     in
   @endvar
   @var     n_writes
   @vdesc   Number of writes clauses
   @vtype   int
   @vio     in
   @endvar
   @var     n_reads
   @vdesc   Number of reads clauses
   @vtype   int
   @vio     in
   @endvar
   @var     n_options
   @vdesc   Number of options for this schedule block
   @vtype   int
   @vio     in
   @endvar
   @var     n_before
   @vdesc   Number of functions/groups to schedule before
   @vtype   int
   @vio     in
   @endvar
   @var     n_after
   @vdesc   Number of functions/groups to schedule after
   @vtype   int
   @vio     in
   @endvar
   @var     n_while
   @vdesc   Number of vars to schedule while
   @vtype   int
   @vio     in
   @endvar
   @var     n_if
   @vdesc   Number of vars to schedule if
   @vtype   int
   @vio     in
   @endvar
   @var     timelevels
   @vdesc   The number of timelevels of the storage groups to enable.
   @vtype   const int *
   @vio     in
   @var     ...
   @vdesc   remaining options
   @vtype   multiple const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   Return val of DoScheduleFunction or
   -1 - memory failure
   @endreturndesc
@@*/
int CCTKi_ScheduleFunction(void *function,
                           const char *name,
                           const char *thorn,
                           const char *implementation,
                           const char *description,
                           const char *where,
                           const char *language,
                           int n_mem_groups,
                           int n_comm_groups,
                           int n_trigger_groups,
                           int n_sync_groups,
                           int n_writes,
                           int n_reads,
                           int n_options,
                           int n_before,
                           int n_after,
                           int n_while,
                           int n_if,
                           const int *timelevels,
                           ...
                           )
{
  int retcode;
  t_attribute *attribute;
  t_sched_modifier *modifier;
  va_list ap;

  va_start(ap, timelevels);

  attribute = CreateAttribute(where,name,description, language, thorn, implementation,
                              n_mem_groups, n_comm_groups, n_trigger_groups,
                              n_sync_groups, n_writes, n_reads,
                              n_options, timelevels, &ap);
  modifier  = CreateModifiers(n_before, n_after, n_while, n_if, &ap);

  va_end(ap);

  ValidateModifiers(modifier);

  if(attribute && (modifier || (n_before == 0 && n_after == 0 &&
                                n_while == 0 && n_if == 0)))
  {
    attribute->FunctionData.type = TranslateFunctionType(where);

    retcode = CCTKi_DoScheduleFunction(where, name, function, modifier, (void *)attribute);

    if(retcode == -2)
    {
      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "Tried to schedule duplicate item '%s' from thorn '%s' in '%s'",
                  name, thorn, where);
    }
#ifdef DEBUG
    fprintf(stderr, "Scheduled %s at %s\n", name, where);
#endif
  }
  else
  {
    fprintf(stderr, "Internal error: Failed to schedule %s at %s!!!\n", name, where);
    exit(2);
    retcode = -1;
  }

  return retcode;
}

/*@@
   @routine    CCTKi_ScheduleGroup
   @date       Thu Sep 16 18:19:18 1999
   @author     Tom Goodale
   @desc
   Schedules a group.
   @enddesc
   @calls

   @var     realname
   @vdesc   real name of group to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     name
   @vdesc   working name of group to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     thorn
   @vdesc   name of thorn providing group to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     implementation
   @vdesc   name of implementation group belongs to
   @vtype   const char *
   @vio     in
   @endvar
   @var     description
   @vdesc   desciption of group to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     where
   @vdesc   where to schedule the group
   @vtype   const char *
   @vio     in
   @endvar
   @var     n_mem_groups
   @vdesc   Number of groups needing memory switched on during this function
   @vtype   int
   @vio     in
   @endvar
   @var     n_comm_groups
   @vdesc   Number of groups needing communication switched on during this function
   @vtype   int
   @vio     in
   @endvar
   @var     n_trigger_groups
   @vdesc   Number of groups to trigger this function on
   @vtype   int
   @vio     in
   @endvar
   @var     n_sync_groups
   @vdesc   Number of groups needing synchronisation after this function
   @vtype   int
   @vio     in
   @endvar
   @var     n_writes
   @vdesc   Number of writes clauses
   @vtype   int
   @vio     in
   @endvar
   @var     n_reads
   @vdesc   Number of reads clauses
   @vtype   int
   @vio     in
   @endvar
   @var     n_options
   @vdesc   Number of options for this schedule block
   @vtype   int
   @vio     in
   @endvar
   @var     n_before
   @vdesc   Number of functions/groups to schedule before
   @vtype   int
   @vio     in
   @endvar
   @var     n_after
   @vdesc   Number of functions/groups to schedule after
   @vtype   int
   @vio     in
   @endvar
   @var     n_while
   @vdesc   Number of vars to schedule while
   @vtype   int
   @vio     in
   @endvar
   @var     n_if
   @vdesc   Number of vars to schedule if
   @vtype   int
   @vio     in
   @endvar
   @var     timelevels
   @vdesc   The number of timelevels of the storage groups to enable.
   @vtype   const int *
   @var     ...
   @vdesc   remaining options
   @vtype   multiple const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   Return val of DoScheduleGroup or
   -1 - memory failure
   @endreturndesc
@@*/
int CCTKi_ScheduleGroup(const char *realname,
                        const char *name,
                        const char *thorn,
                        const char *implementation,
                        const char *description,
                        const char *where,
                        int n_mem_groups,
                        int n_comm_groups,
                        int n_trigger_groups,
                        int n_sync_groups,
                        int n_writes,
                        int n_reads,
                        int n_options,
                        int n_before,
                        int n_after,
                        int n_while,
                        int n_if,
                        const int *timelevels,
                        ...
                        )
{
  int retcode;
  t_attribute *attribute;
  t_sched_modifier *modifier;
  va_list ap;

  va_start(ap, timelevels);

  attribute = CreateAttribute(where,name,description, NULL, thorn, implementation,
                              n_mem_groups, n_comm_groups, n_trigger_groups,
                              n_sync_groups, n_writes, n_reads,
                              n_options, timelevels, &ap);
  modifier  = CreateModifiers(n_before, n_after, n_while, n_if, &ap);

  va_end(ap);

  ValidateModifiers(modifier);

  if(attribute && (modifier || (n_before == 0 && n_after == 0 &&
                                n_writes == 0 && n_reads == 0 &&
                                n_while == 0 && n_if == 0)))
  {
    retcode = CCTKi_DoScheduleGroup(where, name, realname, modifier, (void *)attribute);

    if(retcode == -2)
    {

      CCTK_VWarn (0, __LINE__, __FILE__, "Cactus",
                  "Tried to schedule duplicate item '%s' from thorn '%s' in '%s'",
                  name, thorn, where);
    }

#ifdef DEBUG
    fprintf(stderr, "Scheduled %s at %s\n", name, where);
#endif
  }
  else
  {
#ifdef DEBUG
    fprintf(stderr, "Failed to schedule %s at %s!!!\n", name, where);
#endif
    retcode = -1;
  }

  return retcode;

}


/*@@
   @routine    CCTKi_ScheduleGroupStorage
   @date       Fri Sep 17 18:55:59 1999
   @author     Tom Goodale
   @desc
   Schedules a group for storage when a GH is created.
   @enddesc
   @calls

   @var     group
   @vdesc   group name
   @vtype   const char *
   @vio     in
   @endvar
   @var     timelevels
   @vdesc   The number of timelevels of this group to enable.
   @vtype   int
   @vio     in

   @returntype int
   @returndesc
   Group index or
   -1 - memory failure
   @endreturndesc
@@*/
int CCTKi_ScheduleGroupStorage(const char *group, int timelevels)
{
  int *temp;
  int *temp2;

  temp = realloc(scheduled_storage_groups,
                 (n_scheduled_storage_groups+1) * sizeof(int));
  temp2 = realloc(scheduled_storage_groups_timelevels,
                  (n_scheduled_storage_groups+1) * sizeof(int));

  if(temp && temp2)
  {
    temp[n_scheduled_storage_groups++] = CCTK_GroupIndex(group);
    scheduled_storage_groups = temp;
    scheduled_storage_groups_timelevels = temp2;

    scheduled_storage_groups_timelevels[n_scheduled_storage_groups-1] = timelevels;
  }

  return (temp && temp2  ? temp[n_scheduled_storage_groups-1] : -1);
}


/*@@
   @routine    CCTKi_ScheduleGroupComm
   @date       Fri Sep 17 18:55:59 1999
   @author     Tom Goodale
   @desc
   Schedules a group for communication when a GH is created.
   @enddesc
   @calls

   @var     group
   @vdesc   group name
   @vtype   const char *
   @vio     in
   @endvar
   @returntype int
   @returndesc
   Group index or
   -1 - memory failure
   @endreturndesc
@@*/
int CCTKi_ScheduleGroupComm(const char *group)
{
  int *temp;

  temp = realloc(scheduled_comm_groups,
                 (n_scheduled_comm_groups+1) * sizeof(int));
  if(temp)
  {
    temp[n_scheduled_comm_groups++] = CCTK_GroupIndex(group);
    scheduled_comm_groups = temp;
  }

  return (temp ? temp[n_scheduled_comm_groups-1] : -1);
}


 /*@@
   @routine    CCTK_ScheduleTraverse
   @date       Tue Apr  4 08:05:27 2000
   @author     Tom Goodale
   @desc
   Traverses a schedule point, and its entry and exit points if necessary.
   @enddesc
   @calls

   @var     where
   @vdesc   Schedule point
   @vtype   const char *
   @vio     in
   @endvar
   @var     GH
   @vdesc   GH data
   @vtype   void *
   @vio     inout
   @endvar
   @var     CallFunction
   @vdesc   Function called to call a function
   @vtype   int (*)(void *, cFubctionData, void *)
   @vio     in
   @vcomment
   Set to NULL to use the default
   @endvar

   @returntype int
   @returndesc
   0 - success
   1 - memory failure
   @endreturndesc
@@*/
int CCTK_ScheduleTraverse(const char *where,
                          void *GH,
                          int (*CallFunction)(void *, cFunctionData *, void *))
{
  int retcode;

  int special;
  const char *current;

  static char *current_point = NULL;
  static unsigned int current_length = 0;
  char *temp;

  special=0;

  /* Special entry points have $ in them */
  for(current=where; *current; current++)
  {
    if(*current == '$')
    {
      special = 1;
      break;
    }
  }

  retcode = 0;

  if(special)
  {
    ScheduleTraverse(where, GH, CallFunction);
  }
  else
  {
    if(current_length < strlen(where) + 7)
    {
      current_length = strlen(where)+7;

      temp = realloc(current_point, current_length);

      if(temp)
      {
        current_point = temp;
      }
      else
      {
        retcode = 1;
      }
    }
    if(retcode == 0)
    {
      sprintf(current_point, "%s$ENTRY", where);
      ScheduleTraverse(current_point, GH, CallFunction);

      ScheduleTraverse(where, GH, CallFunction);

      sprintf(current_point, "%s$EXIT", where);
      ScheduleTraverse(current_point, GH, CallFunction);
    }
  }

  return retcode;
}


/*@@
   @routine    CCTKi_ScheduleGHInit
   @date       Fri Sep 17 21:25:13 1999
   @author     Tom Goodale
   @desc
   Does any scheduling stuff setup which requires a GH.
   @enddesc
   @calls

   @var     GH
   @vdesc   GH data
   @vtype   void *
   @vio     inout
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
int CCTKi_ScheduleGHInit(void *GH)
{
  int i;


  /* create and start the CCTK total timer */
  total_timer = CCTK_TimerCreate ("CCTK total time");
  if (total_timer >= 0)
  {
    CCTK_TimerStartI (total_timer);
  }
  else
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "Couldn't create CCTK total timer. "
                "No timing information will be available.");
  }

  if (n_scheduled_storage_groups>0)
  {
    CCTK_GroupStorageIncrease(GH,
                              n_scheduled_storage_groups,
                              scheduled_storage_groups,
                              scheduled_storage_groups_timelevels,
                              NULL);
  }

  for(i = 0; i < n_scheduled_comm_groups; i++)
  {
    CCTK_EnableGroupCommI(GH,scheduled_comm_groups[i]);
  }

  return 0;
}

/*@@
   @routine    CCTK_SchedulePrint
   @date       Fri Sep 17 21:52:44 1999
   @author     Tom Goodale
   @desc
   Prints out the schedule info.
   @enddesc
   @calls

   @var     where
   @vdesc   Schedule point
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
#define PUTS(x) printf("%*s%s\n", indent_level, "", x)
int CCTK_SchedulePrint(const char *where)
{
  if(!where)
  {
    indent_level = 2;
    PUTS("if (recover initial data)");
    PUTS("  Recover parameters");
    PUTS("endif");
    putchar ('\n');
    PUTS("Startup routines");
    PUTS("  [CCTK_STARTUP]");
    SchedulePrint("CCTK_STARTUP");
    putchar ('\n');
    PUTS("Startup routines which need an existing grid hierarchy");
    PUTS("  [CCTK_WRAGH]");
    SchedulePrint("CCTK_WRAGH");
    PUTS("Parameter checking routines");
    PUTS("  [CCTK_PARAMCHECK]");
    SchedulePrint("CCTK_PARAMCHECK");
    putchar ('\n');
    PUTS("Initialisation");
    PUTS("  if (NOT (recover initial data AND recovery_mode is 'strict'))");
    indent_level += 2;
    PUTS("  [CCTK_PREREGRIDINITIAL]");
    SchedulePrint("CCTK_PREREGRIDINITIAL$ENTRY");
    SchedulePrint("CCTK_PREREGRIDINITIAL");
    SchedulePrint("CCTK_PREREGRIDINITIAL$EXIT");
    PUTS("  Set up grid hierarchy");
    PUTS("  [CCTK_POSTREGRIDINITIAL]");
    SchedulePrint("CCTK_POSTREGRIDINITIAL$ENTRY");
    SchedulePrint("CCTK_POSTREGRIDINITIAL");
    SchedulePrint("CCTK_POSTREGRIDINITIAL$EXIT");
    PUTS("  [CCTK_BASEGRID]");
    SchedulePrint("CCTK_BASEGRID$ENTRY");
    SchedulePrint("CCTK_BASEGRID");
    SchedulePrint("CCTK_BASEGRID$EXIT");
    PUTS("  [CCTK_INITIAL]");
    SchedulePrint("CCTK_INITIAL$ENTRY");
    SchedulePrint("CCTK_INITIAL");
    SchedulePrint("CCTK_INITIAL$EXIT");
    PUTS("  [CCTK_POSTINITIAL]");
    SchedulePrint("CCTK_POSTINITIAL$ENTRY");
    SchedulePrint("CCTK_POSTINITIAL");
    SchedulePrint("CCTK_POSTINITIAL$EXIT");
    PUTS("  Initialise finer grids recursively");
    PUTS("  Restrict from finer grids");
    PUTS("  [CCTK_POSTRESTRICTINITIAL]");
    SchedulePrint("CCTK_POSTRESTRICTINITIAL$ENTRY");
    SchedulePrint("CCTK_POSTRESTRICTINITIAL");
    SchedulePrint("CCTK_POSTRESTRICTINITIAL$EXIT");
    PUTS("  [CCTK_POSTPOSTINITIAL]");
    SchedulePrint("CCTK_POSTPOSTINITIAL$ENTRY");
    SchedulePrint("CCTK_POSTPOSTINITIAL");
    SchedulePrint("CCTK_POSTPOSTINITIAL$EXIT");
    PUTS("  [CCTK_POSTSTEP]");
    SchedulePrint("CCTK_POSTSTEP$ENTRY");
    SchedulePrint("CCTK_POSTSTEP");
    SchedulePrint("CCTK_POSTSTEP$EXIT");
    PUTS("endif");
    PUTS("if (recover initial data)");
    PUTS("  [CCTK_BASEGRID]");
    SchedulePrint("CCTK_BASEGRID$ENTRY");
    SchedulePrint("CCTK_BASEGRID");
    SchedulePrint("CCTK_BASEGRID$EXIT");
    PUTS("  [CCTK_RECOVER_VARIABLES]");
    SchedulePrint("CCTK_RECOVER_VARIABLES");
    PUTS("  [CCTK_POST_RECOVER_VARIABLES]");
    SchedulePrint("CCTK_POST_RECOVER_VARIABLES");
    PUTS("endif");
    PUTS("if (checkpoint initial data)");
    PUTS("  [CCTK_CPINITIAL]");
    SchedulePrint("CCTK_CPINITIAL");
    PUTS("endif");
    PUTS("if (analysis)");
    PUTS("  [CCTK_ANALYSIS]");
    SchedulePrint("CCTK_ANALYSIS$ENTRY");
    SchedulePrint("CCTK_ANALYSIS");
    SchedulePrint("CCTK_ANALYSIS$EXIT");
    indent_level -=2;
    PUTS("endif");
    PUTS("Output grid variables");
    putchar ('\n');
    PUTS("do loop over timesteps");
    PUTS("  [CCTK_PREREGRID]");
    SchedulePrint("CCTK_PREREGRID$ENTRY");
    SchedulePrint("CCTK_PREREGRID");
    SchedulePrint("CCTK_PREREGRID$EXIT");
    PUTS("  Change grid hierarchy");
    PUTS("  [CCTK_POSTREGRID]");
    SchedulePrint("CCTK_POSTREGRID$ENTRY");
    SchedulePrint("CCTK_POSTREGRID");
    SchedulePrint("CCTK_POSTREGRID$EXIT");
    PUTS("  Rotate timelevels");
    PUTS("  iteration = iteration+1");
    PUTS("  t = t+dt");
    PUTS("  [CCTK_PRESTEP]");
    SchedulePrint("CCTK_PRESTEP$ENTRY");
    SchedulePrint("CCTK_PRESTEP");
    SchedulePrint("CCTK_PRESTEP$EXIT");
    PUTS("  [CCTK_EVOL]");
    SchedulePrint("CCTK_EVOL$ENTRY");
    SchedulePrint("CCTK_EVOL");
    SchedulePrint("CCTK_EVOL$EXIT");
    PUTS("  Evolve finer grids recursively");
    PUTS("  Restrict from finer grids");
    PUTS("  [CCTK_POSTRESTRICT]");
    SchedulePrint("CCTK_POSTRESTRICT$ENTRY");
    SchedulePrint("CCTK_POSTRESTRICT");
    SchedulePrint("CCTK_POSTRESTRICT$EXIT");
    PUTS("  [CCTK_POSTSTEP]");
    SchedulePrint("CCTK_POSTSTEP$ENTRY");
    SchedulePrint("CCTK_POSTSTEP");
    SchedulePrint("CCTK_POSTSTEP$EXIT");
    PUTS("if (checkpoint)");
    PUTS("  [CCTK_CHECKPOINT]");
    SchedulePrint("CCTK_CHECKPOINT");
    PUTS("endif");
    PUTS("if (analysis)");
    PUTS("  [CCTK_ANALYSIS]");
    SchedulePrint("CCTK_ANALYSIS$ENTRY");
    SchedulePrint("CCTK_ANALYSIS");
    SchedulePrint("CCTK_ANALYSIS$EXIT");
    PUTS("endif");
    PUTS("Output grid variables");
    PUTS("enddo");
    putchar ('\n');
    PUTS("Termination routines");
    PUTS("  [CCTK_TERMINATE]");
    SchedulePrint("CCTK_TERMINATE");
    putchar ('\n');
    PUTS("Shutdown routines");
    PUTS("  [CCTK_SHUTDOWN]");
    SchedulePrint("CCTK_SHUTDOWN");
    putchar ('\n');
    PUTS("Routines run after changing the grid hierarchy:");
    PUTS("  [CCTK_POSTREGRID]");
    SchedulePrint("CCTK_POSTREGRID");
  }
  else
  {
    SchedulePrint(where);
  }

  return (0);
}

/*@@
   @routine    CCTK_SchedulePrintTimes
   @date       2007-01-16
   @author     Erik Schnetter
   @desc
   Prints out the schedule timings.
   @enddesc
   @calls

   @var     where
   @vdesc   Schedule point
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
int CCTK_SchedulePrintTimes(const char *where)
{
  return CCTK_SchedulePrintTimesToFile(where, stdout);
}

/*@@
   @routine    CCTK_SchedulePrintTimesToFile
   @date       Fri Sep 17 21:52:44 1999
   @author     Tom Goodale
   @desc
   Prints out the schedule timings to a file.
   @enddesc
   @calls

   @var     where
   @vdesc   Schedule point
   @vtype   const char *
   @vio     in
   @endvar

   @var     file
   @vdesc   Output file, must be open for writing
   @vtype   FILE *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
int CCTK_SchedulePrintTimesToFile(const char *where, FILE *file)
{
  t_sched_data data;

  if(!timerinfo)
  {
    timerinfo = CCTK_TimerCreateData();
  }

  data.GH = NULL;
  data.schedule_bin = where;
  data.schedpoint = schedpoint_misc;
  data.print_headers = 1;
  data.file = file;
  data.info = timerinfo;
  data.total_time = CCTK_TimerCreateData();

  if(!where)
  {
    SchedulePrintTimes("CCTK_STARTUP", &data);
    SchedulePrintTimes("CCTK_WRAGH", &data);
    SchedulePrintTimes("CCTK_PARAMCHECK", &data);
    SchedulePrintTimes("CCTK_PREREGRIDINITIAL$ENTRY", &data);
    SchedulePrintTimes("CCTK_PREREGRIDINITIAL", &data);
    SchedulePrintTimes("CCTK_PREREGRIDINITIAL$EXIT", &data);
    SchedulePrintTimes("CCTK_POSTREGRIDINITIAL$ENTRY", &data);
    SchedulePrintTimes("CCTK_POSTREGRIDINITIAL", &data);
    SchedulePrintTimes("CCTK_POSTREGRIDINITIAL$EXIT", &data);
    SchedulePrintTimes("CCTK_BASEGRID$ENTRY", &data);
    SchedulePrintTimes("CCTK_BASEGRID", &data);
    SchedulePrintTimes("CCTK_BASEGRID$EXIT", &data);
    SchedulePrintTimes("CCTK_INITIAL$ENTRY", &data);
    SchedulePrintTimes("CCTK_INITIAL", &data);
    SchedulePrintTimes("CCTK_INITIAL$EXIT", &data);
    SchedulePrintTimes("CCTK_POSTINITIAL$ENTRY", &data);
    SchedulePrintTimes("CCTK_POSTINITIAL", &data);
    SchedulePrintTimes("CCTK_POSTINITIAL$EXIT", &data);
    SchedulePrintTimes("CCTK_POSTRESTRICTINITIAL$ENTRY", &data);
    SchedulePrintTimes("CCTK_POSTRESTRICTINITIAL", &data);
    SchedulePrintTimes("CCTK_POSTRESTRICTINITIAL$EXIT", &data);
    SchedulePrintTimes("CCTK_POSTPOSTINITIAL$ENTRY", &data);
    SchedulePrintTimes("CCTK_POSTPOSTINITIAL", &data);
    SchedulePrintTimes("CCTK_POSTPOSTINITIAL$EXIT", &data);
    SchedulePrintTimes("CCTK_RECOVER_VARIABLES", &data);
    SchedulePrintTimes("CCTK_POST_RECOVER_VARIABLES", &data);
    SchedulePrintTimes("CCTK_CPINITIAL", &data);
    SchedulePrintTimes("CCTK_PREREGRID$ENTRY", &data);
    SchedulePrintTimes("CCTK_PREREGRID", &data);
    SchedulePrintTimes("CCTK_PREREGRID$EXIT", &data);
    SchedulePrintTimes("CCTK_POSTREGRID$ENTRY", &data);
    SchedulePrintTimes("CCTK_POSTREGRID", &data);
    SchedulePrintTimes("CCTK_POSTREGRID$EXIT", &data);
    SchedulePrintTimes("CCTK_PRESTEP$ENTRY", &data);
    SchedulePrintTimes("CCTK_PRESTEP", &data);
    SchedulePrintTimes("CCTK_PRESTEP$EXIT", &data);
    SchedulePrintTimes("CCTK_EVOL$ENTRY", &data);
    SchedulePrintTimes("CCTK_EVOL", &data);
    SchedulePrintTimes("CCTK_EVOL$EXIT", &data);
    SchedulePrintTimes("CCTK_POSTRESTRICT$ENTRY", &data);
    SchedulePrintTimes("CCTK_POSTRESTRICT", &data);
    SchedulePrintTimes("CCTK_POSTRESTRICT$EXIT", &data);
    SchedulePrintTimes("CCTK_POSTSTEP$ENTRY", &data);
    SchedulePrintTimes("CCTK_POSTSTEP", &data);
    SchedulePrintTimes("CCTK_POSTSTEP$EXIT", &data);
    SchedulePrintTimes("CCTK_CHECKPOINT", &data);
    SchedulePrintTimes("CCTK_ANALYSIS$ENTRY", &data);
    SchedulePrintTimes("CCTK_ANALYSIS", &data);
    SchedulePrintTimes("CCTK_ANALYSIS$EXIT", &data);
    SchedulePrintTimes("CCTK_TERMINATE", &data);
    SchedulePrintTimes("CCTK_SHUTDOWN", &data);
  }
  else
  {
    SchedulePrintTimes(where, &data);
  }

  CCTK_TimerDestroyData(data.total_time);

  /* also print total time at the bottom */
  if (total_timer >= 0)
  {
    int total_timer_running = CCTK_TimerIsRunningI(total_timer);
    if (total_timer_running)
    {
      CCTK_TimerStopI (total_timer);
    }
    CCTK_TimerI (total_timer, timerinfo);
    CCTKi_SchedulePrintTimerInfo
      (timerinfo, NULL, "", "Total time for simulation", file);

    /* just in case this is not at termination yet ... */
    if (total_timer_running)
    {
      CCTK_TimerStartI (total_timer);
    }
  }

  return 0;
}

/*@@
   @routine    CCTK_TranslateLanguage
   @date       Thu Sep 16 18:18:31 1999
   @author     Tom Goodale
   @desc
   Translates a language string into an internal enum.
   @enddesc
   @calls

   @var     sval
   @vdesc   Language
   @vtype   const char *
   @vio     in
   @endvar

   @returntype cLanguage
   @returndesc
   The language
   @endreturndesc
@@*/
cLanguage CCTK_TranslateLanguage(const char *sval)
{
  cLanguage retcode;

  if(CCTK_Equals(sval, "C"))
  {
    retcode = LangC;
  }
  else if(CCTK_Equals(sval, "Fortran"))
  {
    retcode = LangFortran;
  }
  else
  {
    fprintf(stderr, "Unknown language %s\n", sval);
    retcode = LangNone;
  }

  return retcode;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

/*@@
   @routine    ScheduleTraverse
   @date       Fri Sep 17 21:52:44 1999
   @author     Tom Goodale
   @desc
   Traverses the given schedule point.
   @enddesc
   @calls

   @var     where
   @vdesc   Schedule point
   @vtype   const char *
   @vio     in
   @endvar
   @var     GH
   @vdesc   GH data
   @vtype   void *
   @vio     inout
   @endvar
   @var     CallFunction
   @vdesc   Function called to call a function
   @vtype   int (*)(void *, cFubctionData, void *)
   @vio     in
   @vcomment
   Set to NULL to use the default
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/

static int ScheduleTraverse(const char *where,
                            void *GH,
                            int (*CallFunction)(void *, cFunctionData *, void *))
{
  t_sched_data data;
  int (*calling_function)(void *, t_attribute *, t_sched_data *);

  data.GH = (cGH *)GH;
  data.CallFunction = CallFunction ? CallFunction : CCTK_CallFunction;
  data.schedule_bin = where;
  data.schedpoint = CCTK_Equals(data.schedule_bin, "CCTK_ANALYSIS") ?
                    schedpoint_analysis : schedpoint_misc;
  calling_function = CCTKi_ScheduleCallFunction;

  CCTKi_DoScheduleTraverse(where,
     (int (*)(void *, void *))                    CCTKi_ScheduleCallEntry,
     (int (*)(void *, void *))                    CCTKi_ScheduleCallExit,
     (int (*)(int, char **, void *, void *, int)) CCTKi_ScheduleCallWhile,
     (int (*)(int, char **, void *, void *, int)) CCTKi_ScheduleCallIf,
     (int (*)(void *, void *, void *))            calling_function,
     (void *)&data);

  return 0;
}

 /*@@
   @routine    CreateAttribute
   @date       Thu Sep 16 18:22:48 1999
   @author     Tom Goodale
   @desc
   Creates an attribute structure for a schedule item.
   @enddesc
   @calls

   @var     description
   @vdesc   desciption of function to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     language
   @vdesc   language of function to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     thorn
   @vdesc   name of thorn providing function to be scheduled
   @vtype   const char *
   @vio     in
   @endvar
   @var     implementation
   @vdesc   name of implementation thorn belongs to
   @vtype   const char *
   @vio     in
   @endvar
   @var     n_mem_groups
   @vdesc   Number of groups needing memory switched on during this function
   @vtype   int
   @vio     in
   @endvar
   @var     n_comm_groups
   @vdesc   Number of groups needing communication switched on during this function
   @vtype   int
   @vio     in
   @endvar
   @var     n_trigger_groups
   @vdesc   Number of groups to trigger this function on
   @vtype   int
   @vio     in
   @endvar
   @var     n_sync_groups
   @vdesc   Number of groups needing synchronisation after this function
   @vtype   int
   @vio     in
   @endvar
   @var     n_writes
   @vdesc   Number of writes clauses
   @vtype   int
   @vio     in
   @endvar
   @var     n_reads
   @vdesc   Number of reads clauses
   @vtype   int
   @vio     in
   @endvar
   @var     n_options
   @vdesc   Number of options for this schedule block
   @vtype   int
   @vio     in
   @endvar
   @var     timelevels
   @vdesc   The number of timelevels of the storage groups to enable.
   @vtype   const int *
   @var     ap
   @vdesc   options
   @vtype   va_list of multiple const char *
   @vio     inout
   @vcomment
   This should have as many items as the sum of the above n_* options
   @endvar

   @returntype t_attribute
   @returndesc
   The attribute
   @endreturndesc
@@*/
static t_attribute *CreateAttribute(const char *where,
                                    const char *name,
                                    const char *description,
                                    const char *language,
                                    const char *thorn,
                                    const char *implementation,
                                    int n_mem_groups,
                                    int n_comm_groups,
                                    int n_trigger_groups,
                                    int n_sync_groups,
                                    int n_writes,
                                    int n_reads,
                                    int n_options,
                                    const int *timelevels,
                                    va_list *ap)
{
  t_attribute *this;
  int i;

  this = calloc(1, sizeof(t_attribute));

  if(this)
  {
    this->FunctionData.where = malloc((strlen(where)+1)*sizeof(char));
    this->FunctionData.routine = malloc((strlen(name)+1)*sizeof(char));
    this->description    = malloc((strlen(description)+1)*sizeof(char));
    this->FunctionData.thorn = malloc((strlen(thorn)+1)*sizeof(char));
    this->implementation = malloc((strlen(implementation)+1)*sizeof(char));
    if (n_mem_groups > 0)
    {
      this->mem_groups     = malloc(n_mem_groups*sizeof(int));
      this->timelevels     = malloc(n_mem_groups*sizeof(int));
      this->StorageOnEntry = malloc(n_mem_groups*sizeof(int));
    }
    if (n_trigger_groups > 0)
    {
      this->FunctionData.TriggerGroups = malloc(n_trigger_groups*sizeof(int));
    }
    if (n_sync_groups > 0)
    {
      this->FunctionData.SyncGroups = malloc(n_sync_groups*sizeof(int));
    }
    if (n_writes > 0)
    {
      this->FunctionData.WritesClauses = malloc(n_writes*sizeof(char*));
    }
    if (n_reads > 0)
    {
      this->FunctionData.ReadsClauses = malloc(n_reads*sizeof(char*));
    }
    if (n_comm_groups > 0)
    {
      this->comm_groups    = malloc(n_comm_groups*sizeof(int));
      this->CommOnEntry    = malloc(n_comm_groups*sizeof(int));
    }

    if(this->FunctionData.where &&
       this->FunctionData.routine &&
       this->description     &&
       this->FunctionData.thorn  &&
       this->implementation  &&
       (this->mem_groups || n_mem_groups==0)         &&
       (this->timelevels || n_mem_groups==0)         &&
       (this->StorageOnEntry || n_mem_groups==0)         &&
       (this->comm_groups || n_comm_groups==0)       &&
       (this->FunctionData.TriggerGroups || n_trigger_groups==0) &&
       (this->FunctionData.SyncGroups || n_sync_groups==0) &&
       (this->FunctionData.WritesClauses || n_writes==0) &&
       (this->FunctionData.ReadsClauses || n_reads==0))
    {
      strcpy(this->FunctionData.where,where);
      strcpy(this->FunctionData.routine,name);
      strcpy(this->description,    description);
      strcpy(this->FunctionData.thorn, thorn);
      strcpy(this->implementation, implementation);

      if(language)
      {
        this->type = sched_function;
        this->FunctionData.language = CCTK_TranslateLanguage(language);
        this->FunctionData.FortranCaller = (int (*)(cGH *,void *))CCTKi_FortranWrapper(thorn);
      }
      else
      {
        this->type = sched_group;
      }

      /* Create the lists of indices of groups we're interested in. */
      CreateGroupIndexList(n_mem_groups,     this->mem_groups, ap);
      CreateGroupIndexList(n_comm_groups,    this->comm_groups, ap);
      CreateGroupIndexList(n_trigger_groups, this->FunctionData.TriggerGroups, ap);
      CreateGroupIndexList(n_sync_groups,    this->FunctionData.SyncGroups, ap);
      CreateStringList    (n_writes,       this->FunctionData.WritesClauses, ap);
      CreateStringList    (n_reads,       this->FunctionData.ReadsClauses, ap);

      for(i=0; i< n_mem_groups; i++)
      {
        this->timelevels[i] = timelevels[i];
      }

      /* Check the miscellaneous options */

      InitialiseOptionList(this);
      ParseOptionList(n_options, this, ap);

      /* Check the tags */
      ParseTagsTable(this, ap);

      this->n_mem_groups     = n_mem_groups;
      this->n_comm_groups    = n_comm_groups;
      this->FunctionData.n_TriggerGroups = n_trigger_groups;
      this->FunctionData.n_SyncGroups = n_sync_groups;
      this->FunctionData.n_WritesClauses = n_writes;
      this->FunctionData.n_ReadsClauses = n_reads;

      this->timers = NULL;
    }
    else
    {
      free(this->FunctionData.where);
      free(this->FunctionData.routine);
      free(this->description);
      free(this->comm_groups);
      free(this->FunctionData.TriggerGroups);
      free(this->FunctionData.SyncGroups);
      free(this->FunctionData.WritesClauses);
      free(this->FunctionData.ReadsClauses);
      free(this);
      this = NULL;
    }
  }

  return this;
}

/*@@
   @routine    CreateModifier
   @date       Thu Sep 16 18:23:13 1999
   @author     Tom Goodale
   @desc
   Creates a schedule modifier list.
   @enddesc
   @calls

   @var     n_before
   @vdesc   Number of functions/groups to schedule before
   @vtype   int
   @vio     in
   @endvar
   @var     n_after
   @vdesc   Number of functions/groups to schedule after
   @vtype   int
   @vio     in
   @endvar
   @var     n_while
   @vdesc   Number of vars to schedule while
   @vtype   int
   @vio     in
   @endvar
   @var     n_if
   @vdesc   Number of vars to schedule if
   @vtype   int
   @vio     in
   @endvar
   @var     ap
   @vdesc   options
   @vtype   va_list of multiple const char *
   @vio     inout
   @vcomment
   This should have as many items as the sum of the above n_* options
   @endvar

   @returntype t_sched_modifier *
   @returndesc
   the schedule modifier
   @endreturndesc
@@*/
static t_sched_modifier *CreateModifiers(int n_before,
                                         int n_after,
                                         int n_while,
                                         int n_if,
                                         va_list *ap)
{
  t_sched_modifier *modifier;

  modifier = CreateTypedModifier(NULL, "before", n_before, ap);
  modifier = CreateTypedModifier(modifier, "after", n_after, ap);
  modifier = CreateTypedModifier(modifier, "while", n_while, ap);
  modifier = CreateTypedModifier(modifier, "if", n_if, ap);

  return modifier;
}

/*@@
   @routine    ValidateModifier
   @date       Sat Apr 14 18:28:13 2001
   @author     Gabrielle Allen
   @desc
   Validates a schedule modifier list. At the moment just check that
   the while and if modifiers use a CCTK_INT grid variable.
   @enddesc
   @calls

   @var     modifier
   @vdesc
   @vtype   t_sched_modifier *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   Negative if modifier not valid, zero if modifier is valid.
   @endreturndesc
@@*/
int ValidateModifiers(t_sched_modifier *modifier)
{
  int retval = 0;
  int vindex;
  int type;

  for (;modifier;modifier=modifier->next)
  {
    if (modifier->type != sched_while ||
        modifier->type != sched_if)
    {
      continue;
    }
    vindex = CCTK_VarIndex(modifier->argument);
    if (vindex < 0)
    {
      CCTK_VWarn(0,__LINE__,__FILE__,"Cactus",
                 "While qualifier %s could not be parsed as fully "
                 "specified grid variable name (e.g. MyThorn::MyVar",
                 modifier->argument);
      retval = -1;
    }
    else
    {
      type = CCTK_VarTypeI(vindex);
      if (type != CCTK_VARIABLE_INT)
      {
        if (modifier->type == sched_while)
        {
          CCTK_VWarn(0,__LINE__,__FILE__,"Cactus",
                     "While qualifier %s is not a CCTK_INT grid variable",
                     modifier->argument);
        }
        else if (modifier->type == sched_if)
        {
          CCTK_VWarn(0,__LINE__,__FILE__,"Cactus",
                     "If qualifier %s is not a CCTK_INT grid variable",
                     modifier->argument);
        }
        retval = -1;
      }
    }
  }
  return retval;
}

/*@@
   @routine    CreateGroupIndexList
   @date       Fri Sep 17 21:51:51 1999
   @author     Tom Goodale
   @desc
   Gets the next n_items group names from the variable argument list
   and converts them to indices.
   @enddesc
   @calls

   @var     n_items
   @vdesc   number of items on the list
   @vtype   int
   @vio     in
   @endvar
   @var     array
   @vdesc   array of indices
   @vtype   int *
   @vio     out
   @endvar
   @var     ap
   @vdesc   argument list
   @vtype   va_list of const char *
   @vio     inout
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
static int CreateGroupIndexList(int n_items, int *array, va_list *ap)
{
  int i;
  const char *item;

  for(i=0; i < n_items; i++)
  {
    item = va_arg(*ap, const char *);

    array[i] = CCTK_GroupIndex(item);
  }

  return 0;
}


/*@@
   @routine    CreateStringList
   @date       2007-05-24
   @author     Erik Schnetter
   @desc
   Gets the next n_items group names from the variable argument list
   and converts them to strings.
   @enddesc
   @calls

   @var     n_items
   @vdesc   number of items on the list
   @vtype   int
   @vio     in
   @endvar
   @var     array
   @vdesc   array of strings
   @vtype   char **
   @vio     out
   @endvar
   @var     ap
   @vdesc   argument list
   @vtype   va_list of const char *
   @vio     inout
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
static int CreateStringList(int n_items, const char **array, va_list *ap)
{
  int i;
  const char *item;

  for(i=0; i < n_items; i++)
  {
    item = va_arg(*ap, const char *);

    array[i] = strdup(item);
  }

  return 0;
}


 /*@@
   @routine    ParseOptionList
   @date       Thu Jan 27 20:26:42 2000
   @author     Tom Goodale
   @desc
   Extracts the list of miscellaneous options in a schedule
   group definition.
   @enddesc
   @calls

   @var     n_items
   @vdesc   number of items on the list
   @vtype   int
   @vio     in
   @endvar
   @var     attribute
   @vdesc   attribute list
   @vtype   t_attribute *
   @vio     inout
   @endvar
   @var     ap
   @vdesc   argument list
   @vtype   va_list of const char *
   @vio     inout
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
static int ParseOptionList(int n_items,
                           t_attribute *attribute,
                           va_list *ap)
{
  int i;
  const char *item;

  for(i=0; i < n_items; i++)
  {
    item = va_arg(*ap, const char *);

    ParseOption(attribute, item);
  }

  return 0;
}

 /*@@
   @routine    InitialiseOptionList
   @date       Thu Jan 27 20:36:54 2000
   @author     Tom Goodale
   @desc
   Initialises the miscellaneous option list for a schedule group.
   @enddesc
   @calls

   @var     attribute
   @vdesc   option attribute
   @vtype   t_attribute *
   @vio     out
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
static int InitialiseOptionList(t_attribute *attribute)
{
  attribute->FunctionData.meta = 0;
  attribute->FunctionData.meta_early = 0;
  attribute->FunctionData.meta_late = 0;
  attribute->FunctionData.global = 0;
  attribute->FunctionData.global_early = 0;
  attribute->FunctionData.global_late = 0;
  attribute->FunctionData.level = 0;
  attribute->FunctionData.singlemap = 0;
  attribute->FunctionData.local = 0;

  attribute->FunctionData.loop_meta = 0;
  attribute->FunctionData.loop_global = 0;
  attribute->FunctionData.loop_level = 0;
  attribute->FunctionData.loop_singlemap = 0;
  attribute->FunctionData.loop_local = 0;

  return 0;
}

 /*@@
   @routine    ParseOption
   @date       Thu Jan 27 20:29:36 2000
   @author     Tom Goodale
   @desc
   Parses an individual option to a schedule group.
   @enddesc
   @calls

   @var     attribute
   @vdesc   option attribute
   @vtype   t_attribute *
   @vio     out
   @endvar
   @var     option
   @vdesc   Option
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
static int ParseOption(t_attribute *attribute,
                       const char *option)
{
  if(CCTK_Equals(option, "META"))
  {
    attribute->FunctionData.meta = 1;
  }
  else if(CCTK_Equals(option, "META-EARLY"))
  {
    attribute->FunctionData.meta_early = 1;
  }
  else if(CCTK_Equals(option, "META-LATE"))
  {
    attribute->FunctionData.meta_late = 1;
  }
  else if(CCTK_Equals(option, "GLOBAL"))
  {
    attribute->FunctionData.global = 1;
  }
  else if(CCTK_Equals(option, "GLOBAL-EARLY"))
  {
    attribute->FunctionData.global_early = 1;
  }
  else if(CCTK_Equals(option, "GLOBAL-LATE"))
  {
    attribute->FunctionData.global_late = 1;
  }
  else if(CCTK_Equals(option, "LEVEL"))
  {
    attribute->FunctionData.level = 1;
  }
  else if(CCTK_Equals(option, "SINGLEMAP"))
  {
    attribute->FunctionData.singlemap = 1;
  }
  else if(CCTK_Equals(option, "LOCAL"))
  {
    attribute->FunctionData.local = 1;
  }
  else if(CCTK_Equals(option, "LOOP-META"))
  {
    attribute->FunctionData.loop_meta = 1;
  }
  else if(CCTK_Equals(option, "LOOP-GLOBAL"))
  {
    attribute->FunctionData.loop_global = 1;
  }
  else if(CCTK_Equals(option, "LOOP-LEVEL"))
  {
    attribute->FunctionData.loop_level = 1;
  }
  else if(CCTK_Equals(option, "LOOP-SINGLEMAP"))
  {
    attribute->FunctionData.loop_singlemap = 1;
  }
  else if(CCTK_Equals(option, "LOOP-LOCAL"))
  {
    attribute->FunctionData.loop_local = 1;
  }
  else
  {
    CCTK_VWarn(1,__LINE__,__FILE__,"Cactus",
               "ParseOption: Unknown option \"%s\" for schedule item %s::%s (scheduled at %s)",
               option,
               attribute->FunctionData.thorn, attribute->FunctionData.routine,
               attribute->FunctionData.where);
  }

  return 0;
}


 /*@@
   @routine    ParseTagsTable
   @date       2006-08-03
   @author     Erik Schnetter
   @desc
   Extract the tags table in a schedule group definition.
   @enddesc
   @calls

   @var     attribute
   @vdesc   attribute list
   @vtype   t_attribute *
   @vio     inout
   @endvar
   @var     ap
   @vdesc   argument list
   @vtype   va_list of const char *
   @vio     inout
   @endvar

   @returntype int
   @returndesc
   0 - success
   @endreturndesc
@@*/
static int ParseTagsTable(t_attribute *attribute,
                          va_list *ap)
{
  const char *item;
  int table;

  item = va_arg(*ap, const char *);
  table = Util_TableCreateFromString (item);
  if (table < 0)
  {
    return table;
  }
  attribute->FunctionData.tags = table;

  return 0;
}


/*@@
   @routine    CreateTypedModifier
   @date       Fri Sep 17 21:50:59 1999
   @author     Tom Goodale
   @desc
   Adds the next n_items items from the variable argument list
   onto the modifer.
   @enddesc
   @calls

   @var     modifier
   @vdesc   base schedule modifier
   @vtype   t_sched_modifier
   @vio     inout
   @vcomment
   This is a list which gets expanded by this function
   @endvar
   @var     type
   @vdesc   modifier type
   @vtype   const char *
   @vio     in
   @vcomment
   before, after, while, if
   @endvar
   @var     n_items
   @vdesc   Number of items on list
   @vtype   int
   @vio     in
   @endvar
   @var     ap
   @vdesc   argument list
   @vtype   va_list of const char *
   @vio     inout
   @endvar

   @returntype t_sched_modifier *
   @returndesc
   modifier list
   @endreturndesc

@@*/
static t_sched_modifier *CreateTypedModifier(t_sched_modifier *modifier,
                                             const char *type,
                                             int n_items,
                                             va_list *ap)
{
  int i;
  const char *item;

  for(i=0; i < n_items; i++)
  {
    item = va_arg(*ap, const char *);

    modifier = CCTKi_ScheduleAddModifier(modifier, type, item);
  }

  return modifier;
}

/*@@
   @routine    TranslateFunctionType
   @date       Mon Jan 24 16:52:06 2000
   @author     Tom Goodale
   @desc
   Translates a string saying what schedule point
   a function is registered at into the appropriate
   function type.
   @enddesc
   @calls

   @var     where
   @vdesc   schedule point
   @vtype   const char *
   @vio     in
   @endvar

   @returntype cFunctionType
   @returndesc
   The function type
   @endreturndesc
@@*/
static cFunctionType TranslateFunctionType(const char *where)
{
  cFunctionType retcode;

  int special;
  const char *current;

  special = 0;

  /* Special entry points have $ in them */
  for(current=where; *current; current++)
  {
    if(*current == '$')
    {
      special = 1;
      break;
    }
  }

  if(special)
  {
    retcode = FunctionOneArg;
  }
  else if(CCTK_Equals(where, "CCTK_STARTUP"))
  {
    retcode = FunctionNoArgs;
  }
  else if(CCTK_Equals(where, "CCTK_SHUTDOWN"))
  {
    retcode = FunctionNoArgs;
  }
  else
  {
    retcode = FunctionStandard;
  }

  return retcode;
}

/*@@
   @routine    SchedulePrint
   @date       Sun Sep 19 13:31:23 1999
   @author     Tom Goodale
   @desc
   Traverses the schedule data for a particular entry point and
   prints out the data.
   @enddesc
   @calls

   @var     where
   @vdesc   Schedule point
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   return code of DoScheduleTraverse or
   0 - where is NULL
   @endreturndesc
@@*/
static int SchedulePrint(const char *where)
{
  int retcode;
  t_sched_data data;

  data.GH = NULL;
  data.schedule_bin = where;
  data.schedpoint = schedpoint_misc;

  if(where)
  {
    retcode = CCTKi_DoScheduleTraverse(where,
       (int (*)(void *, void *))                    CCTKi_SchedulePrintEntry,
       (int (*)(void *, void *))                    CCTKi_SchedulePrintExit,
       (int (*)(int, char **, void *, void *, int)) CCTKi_SchedulePrintWhile,
       (int (*)(int, char **, void *, void *, int)) CCTKi_SchedulePrintIf,
       (int (*)(void *, void *, void *))            CCTKi_SchedulePrintFunction,
       (void *)&data);
  }
  else
  {
    retcode = 0;
  }

  return retcode;
}

/*@@
   @routine    SchedulePrintTimes
   @date       Fri Oct 22 12:35:06 1999
   @author     Tom Goodale
   @desc
   Prints the times for a particular schedule entry point.
   @enddesc
   @calls

   @var     where
   @vdesc   Schedule point
   @vtype   const char *
   @vio     in
   @endvar
   @var     data
   @vdesc   schedule data
   @vtype   t_sched_data
   @vio     in
   @endvar
   @var     file
   @vdesc   Output file, must be open for writing
   @vtype   FILE *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   return code of DoScheduleTraverse or
   0 - where is NULL
   @endreturndesc
@@*/
static int SchedulePrintTimes(const char *where, t_sched_data *data)
{
  int i;
  int retcode;
  char *description;

  if(where)
  {
    data->schedule_bin = where;
    memset (data->total_time->vals, 0,
            data->total_time->n_vals * sizeof (data->total_time->vals[0]));

    retcode = CCTKi_DoScheduleTraverse(where, NULL, NULL, NULL, NULL,
       (int (*)(void *, void *, void *)) CCTKi_ScheduleResetTimerOutputFlag,
       NULL);

    data->n_functions = 0;
    retcode = CCTKi_DoScheduleTraverse(where, NULL, NULL, NULL, NULL,
       (int (*)(void *, void *, void *)) CCTKi_SchedulePrintTimesFunction,
       (void *)data);

    /* print total time for this schedule bin
       if any routines have been called in it */
    if (retcode >= 0 && data->n_functions > 0)
    {
      for (i = 0; i < data->total_time->n_vals; i++)
      {
        data->total_time->vals[i].type = data->info->vals[i].type;
        data->total_time->vals[i].units = data->info->vals[i].units;
        data->total_time->vals[i].heading = data->info->vals[i].heading;
      }
      description = malloc (strlen (where) + 16);
      sprintf (description, "Total time for %s", where);
      CCTKi_SchedulePrintTimerInfo(data->total_time, NULL, "", description,
                                   data->file);
      free (description);
    }
  }
  else
  {
    retcode = 0;
  }

  return retcode;
}

/********************************************************************
 *********************     Printing Routines   **********************
 ********************************************************************/


 /*@@
   @routine    CCTKi_SchedulePrintEntry
   @date       Sun Sep 19 13:31:23 1999
   @author     Tom Goodale
   @desc
   Routine called on entry to a group when traversing for printing.
   @enddesc
   @calls

   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - schedule item is inactive
   1 - schedule item is active
   @endreturndesc
@@*/
static int CCTKi_SchedulePrintEntry(t_attribute *attribute,
                                    t_sched_data *data)
{
  /* prevent compiler warnings about unused parameters */
  data = data;

  indent_level += 2;

  if (attribute && attribute->type == sched_group)
  {
    printf("%*s %s: %s\n", indent_level + 5, "GROUP",
           attribute->FunctionData.routine,attribute->description);
  }

  return 1;
}

/*@@
   @routine    CCTKi_SchedulePrintExit
   @date       Sun Sep 19 13:31:23 1999
   @author     Tom Goodale
   @desc
   Routine called on exit to a group when traversing for printing.
   @enddesc
   @calls

   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - this has no meaning
   @endreturndesc
@@*/
static int CCTKi_SchedulePrintExit(t_attribute *attribute,
                                   t_sched_data *data)
{
  /* prevent compiler warnings about unused parameters */
  data       = data;
  attribute  = attribute;

  indent_level -= 2;

  return 1;
}

/*@@
   @routine    CCTKi_SchedulePrintWhile
   @date       Sun Sep 19 13:31:23 1999
   @author     Tom Goodale
   @desc
   Routine called for while of a group when traversing for printing.
   @enddesc
   @calls

   @var     n_whiles
   @vdesc   number of while statements
   @vtype   int
   @vio     in
   @endvar
   @var     whiles
   @vdesc   while statements
   @vtype   char **
   @vio     in
   @endvar
   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar
   @var     first
   @vdesc   flag - is this the first time we are checking while on this schedule item
   @vtype   int
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - schedule item is inactive
   1 - schedule item is active
   @endreturndesc
@@*/
static int CCTKi_SchedulePrintWhile(int n_whiles,
                                    char **whiles,
                                    t_attribute *attribute,
                                    t_sched_data *data,
                                    int first)
{
  int i;

  /* prevent compiler warnings about unused parameters */
  attribute = attribute;
  data = data;

  if(first)
  {
    printf("%*s", indent_level + 2 + 7, "while (");

    for(i = 0; i < n_whiles; i++)
    {
      if(i > 0)
      {
        printf(" && ");
      }

      printf("%s", whiles[i]);
    }
    printf(")\n");
    indent_level += 2;
  }
  else
  {
    indent_level -= 2;
    printf("%*s\n", indent_level + 9, "end while");
  }

  return first;
}

/*@@
   @routine    CCTKi_SchedulePrintIf
   @date       Dec 27, 2005
   @author     Erik Schnetter
   @desc
   Routine called for if of a group when traversing for printing.
   @enddesc
   @calls

   @var     n_ifs
   @vdesc   number of if statements
   @vtype   int
   @vio     in
   @endvar
   @var     ifs
   @vdesc   if statements
   @vtype   char **
   @vio     in
   @endvar
   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar
   @var     first
   @vdesc   flag - 1 before the item is entered, 0 afterwards
   @vtype   int
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - schedule item is inactive
   1 - schedule item is active
   @endreturndesc
@@*/
static int CCTKi_SchedulePrintIf(int n_ifs,
                                 char **ifs,
                                 t_attribute *attribute,
                                 t_sched_data *data,
                                 int first)
{
  int i;

  /* prevent compiler warnings about unused parameters */
  attribute = attribute;
  data = data;

  if(first)
  {
    printf("%*s", indent_level + 2 + 4, "if (");
    for(i = 0; i < n_ifs; i++)
    {
      if(i > 0)
      {
        printf(" && ");
      }

      printf("%s", ifs[i]);
    }
    printf(")\n");
    indent_level += 2;
  }
  else
  {
    indent_level -= 2;
  }

  return first;
}

/*@@
   @routine    CCTKi_SchedulePrintFunction
   @date       Sun Sep 19 13:36:25 1999
   @author     Tom Goodale
   @desc
   Function which actually prints out data about a group or a function.
   @enddesc
   @calls

   @var     function
   @vdesc   the function to be called
   @vtype   void *
   @vio     in
   @endvar
   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - this has no meaning
   @endreturndesc
@@*/
static int CCTKi_SchedulePrintFunction(void *function,
                                       t_attribute *attribute,
                                       t_sched_data *data)
{
  const char* mode;
  const char* loop_mode;
  /* prevent compiler warnings about unused parameters */
  function = function;
  data = data;

  if (attribute->FunctionData.meta ||
      attribute->FunctionData.meta_early ||
      attribute->FunctionData.meta_late)
  {
    mode = "[meta] ";
  }
  else if (attribute->FunctionData.global ||
           attribute->FunctionData.global_early ||
           attribute->FunctionData.global_late)
  {
    mode = "[global] ";
  }
  else if (attribute->FunctionData.level)
  {
    mode = "[level] ";
  }
  else if (attribute->FunctionData.singlemap)
  {
    mode = "[singlemap] ";
  }
  else if (attribute->FunctionData.local)
  {
    mode = "[local] ";
  }
  else
  {
    mode = "";
  }

  if (attribute->FunctionData.loop_meta)
  {
    loop_mode = "[loop-meta] ";
  }
  else if (attribute->FunctionData.loop_global)
  {
    loop_mode = "[loop-global] ";
  }
  else if (attribute->FunctionData.loop_level)
  {
    loop_mode = "[loop-level] ";
  }
  else if (attribute->FunctionData.loop_singlemap)
  {
    loop_mode = "[loop-singlemap] ";
  }
  else if (attribute->FunctionData.loop_local)
  {
    loop_mode = "[loop-local] ";
  }
  else
  {
    loop_mode = "";
  }

  if (indent_level > 0)
  {
    printf ("%*s", indent_level, " ");
  }
  printf("%s::%s: %s%s%s\n",
         attribute->FunctionData.thorn, attribute->FunctionData.routine,
         mode, loop_mode,
         attribute->description);

  return 1;
}


/********************************************************************
 *********************     Calling Routines   ***********************
 ********************************************************************/


 /*@@
   @routine    CCTKi_ScheduleCallEntry
   @date       Sun Sep 19 13:24:06 1999
   @author     Tom Goodale
   @desc
   Routine called when a schedule group is entered.
   @enddesc
   @calls

   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - schedule item is inactive
   1 - schedule item is active
   @endreturndesc
@@*/
static int CCTKi_ScheduleCallEntry(t_attribute *attribute,
                                   t_sched_data *data)
{
  int i;
  int indx;
  int last;
  int go;

  if(attribute)
  {
    if(data->schedpoint == schedpoint_analysis)
    {
      /* In analysis, so check triggers */
      if(attribute->FunctionData.n_TriggerGroups == 0)
      {
        go = 1;
      }
      else
      {
        /* Check if it is now being triggered */
        go = 0;
        for (i = 0; i < attribute->FunctionData.n_TriggerGroups ; i++)
        {
          indx = CCTK_FirstVarIndexI(attribute->FunctionData.TriggerGroups[i]);
          last  = indx + CCTK_NumVarsInGroupI(attribute->FunctionData.TriggerGroups[i]) -1;
          for(; indx <= last ; indx++)
          {
            go = go || CCTKi_TriggerSaysGo(data->GH, indx);
          }
        }
      }
    }
    else
    {
      go = 1;
    }

    if(go)
    {
      /* Switch on storage for groups */
/*       for(i = 0; i < attribute->n_mem_groups; i++) */
/*       { */
/*         attribute->StorageOnEntry[i] = CCTK_EnableGroupStorageI(data->GH,attribute->mem_groups[i]); */
/*       } */

      if(attribute->n_mem_groups > 0)
      {
        CCTK_GroupStorageIncrease(data->GH,
                                  attribute->n_mem_groups,
                                  attribute->mem_groups,
                                  attribute->timelevels,
                                  attribute->StorageOnEntry);
      }

      /* Switch on communication for groups. */
      for(i = 0; i < attribute->n_comm_groups; i++)
      {
        attribute->CommOnEntry[i] = CCTK_EnableGroupCommI(data->GH,attribute->comm_groups[i]);
      }
    }

    /* Initialise the synchronised flag. */
    attribute->synchronised = 0;

    /* Remember if we have switched on storage and comm or not. */
    attribute->done_entry = go;
  }
  else
  {
    go = 1;
  }

  return go;
}

/*@@
   @routine    CCTKi_ScheduleCallExit
   @date       Sun Sep 19 13:25:24 1999
   @author     Tom Goodale
   @desc
   Routine called on exit from a schedule group.
   @enddesc
   @calls

   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - this has no meaning
   @endreturndesc
@@*/
static int CCTKi_ScheduleCallExit(t_attribute *attribute,
                                  t_sched_data *data)
{
  int i;
  int vindex;
  int last;

  /* Only do this if the entry routine did stuff. */
  if(attribute && attribute->done_entry)
  {

    /* Synchronise variable groups associated with this schedule group. */
    if(attribute->FunctionData.n_SyncGroups > 0 && ! attribute->synchronised)
    {
      CCTK_SyncGroupsI(data->GH,
                       attribute->FunctionData.n_SyncGroups,
                       attribute->FunctionData.SyncGroups);
      attribute->synchronised = 0;
    }

    if(data->schedpoint == schedpoint_analysis)
    {
      /* In analysis, so do any trigger actions. */
      for (i = 0; i < attribute->FunctionData.n_TriggerGroups ; i++)
      {
        vindex = CCTK_FirstVarIndexI(attribute->FunctionData.TriggerGroups[i]);
        last  = vindex + CCTK_NumVarsInGroupI(attribute->FunctionData.TriggerGroups[i]) - 1;
        for(; vindex <= last ; vindex++)
        {
          CCTKi_TriggerAction(data->GH, vindex);
        }
      }
    }

    /* Switch off communication if it was done in entry. */
    for(i = 0; i < attribute->n_comm_groups; i++)
    {
      if(!attribute->CommOnEntry[i])
      {
        CCTK_DisableGroupCommI(data->GH,attribute->comm_groups[i]);
      }
    }

    /* Switch off storage if it was switched on in entry. */
/*     for(i = 0; i < attribute->n_mem_groups; i++) */
/*     { */
/*       if(!attribute->StorageOnEntry[i]) */
/*       { */
/*         CCTK_DisableGroupStorageI(data->GH,attribute->mem_groups[i]); */
/*       } */
/*     } */

    if(attribute->n_mem_groups > 0)
    {
      CCTK_GroupStorageDecrease(data->GH,
                                attribute->n_mem_groups,
                                attribute->mem_groups,
                                attribute->StorageOnEntry,
                                NULL);
    }
  }

  return 1;
}

/*@@
   @routine    CCTKi_ScheduleCallWhile
   @date       Sun Sep 19 13:27:53 1999
   @author     Tom Goodale
   @desc
   Routine called to check variables to see if a group or function should be executed.
   @enddesc
   @calls

   @var     n_whiles
   @vdesc   number of while statements
   @vtype   int
   @vio     in
   @endvar
   @var     whiles
   @vdesc   while statements
   @vtype   char **
   @vio     in
   @endvar
   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar
   @var     first
   @vdesc   flag - is this the first time we are checking while on this schedule item
   @vtype   int
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - schedule item is inactive
   1 - schedule item is active
   @endreturndesc
@@*/
static int CCTKi_ScheduleCallWhile(int n_whiles,
                                   char **whiles,
                                   t_attribute *attribute,
                                   t_sched_data *data,
                                   int first)
{
  int i;
  int retcode;

  /* prevent compiler warnings about unused parameters */
  attribute = attribute;
  first = first;

  retcode = 1;

  /* FIXME - should do a lot of validation either here or on registration */
  for(i = 0; i < n_whiles; i++)
  {
    retcode = retcode && *((CCTK_INT *)CCTK_VarDataPtr(data->GH, 0, whiles[i]));
  }

  return retcode;
}

/*@@
   @routine    CCTKi_ScheduleCallIf
   @date       Dec 2007, 2005
   @author     Erik Schnetter
   @desc
   Routine called to check variables to see if a group or function should be executed.
   @enddesc
   @calls

   @var     n_ifs
   @vdesc   number of if statements
   @vtype   int
   @vio     in
   @endvar
   @var     ifs
   @vdesc   if statements
   @vtype   char **
   @vio     in
   @endvar
   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0 - schedule item is inactive
   1 - schedule item is active
   @endreturndesc
@@*/
static int CCTKi_ScheduleCallIf(int n_ifs,
                                char **ifs,
                                t_attribute *attribute,
                                t_sched_data *data)
{
  int i;
  int retcode;

  /* prevent compiler warnings about unused parameters */
  attribute = attribute;

  retcode = 1;

  /* FIXME - should do a lot of validation either here or on registration */
  for(i = 0; i < n_ifs; i++)
  {
    retcode = retcode && *((CCTK_INT *)CCTK_VarDataPtr(data->GH, 0, ifs[i]));
  }

  return retcode;
}

/*@@
   @routine    CCTKi_ScheduleCallFunction
   @date       Sun Sep 19 13:29:14 1999
   @author     Tom Goodale
   @desc
   The routine which actually calls a function.
   @enddesc
   @calls

   @var     function
   @vdesc   the function to be called
   @vtype   void *
   @vio     in
   @endvar
   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - this has no meaning
   @endreturndesc
@@*/
static int CCTKi_ScheduleCallFunction(void *function,
                                      t_attribute *attribute,
                                      t_sched_data *data)
{
  /* find the timer for this function and this schedule bin */
  t_timer *timer = attribute->timers;
  while (timer && strcmp(timer->schedule_bin, data->schedule_bin))
  {
    timer = timer->next;
  }

  /* create the timer if it doesn't exist yet */
  if (! timer)
  {
    /* build a unique name for this timer */
    const char *thorn = attribute->FunctionData.thorn;
    const char *name = attribute->FunctionData.routine;
    const char *where = data->schedule_bin;
    char *timername = malloc (strlen (thorn) + strlen (name) +
                              strlen (where) + 20);
    if (! timername)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                  "Could not allocate memory for timer name");
    }
    else
    {
      static int timernum = 0;
      sprintf (timername, "[%04d] %s: %s in %s", timernum++, thorn, name,where);
      const int timer_handle = CCTK_TimerCreate(timername);
      if (timer_handle < 0)
      {
        CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                    "Could not create timer with name '%s'", timername);
      }
      else
      {
        timer = malloc (sizeof (*timer));
        if (! timer)
        {
          CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                      "Could not allocate memory for timer with name '%s'",
                      timername);
          CCTK_TimerDestroy (timername);
        }
        else
        {
          timer->timer_handle = timer_handle;
          timer->schedule_bin = strdup (where);
          timer->next = attribute->timers;
          attribute->timers = timer;
        }
      }
      free (timername);
    }
  }

  if (timer)
  {
    CCTK_TimerStartI(timer->timer_handle);
  }

  /* Use whatever has been chosen as the calling function for this
   * function.
   */
  attribute->synchronised = data->CallFunction(function, &(attribute->FunctionData), data->GH);

  if (timer)
  {
    CCTK_TimerStopI(timer->timer_handle);
  }

  return 1;
}

/********************************************************************
 ****************     Timer Printing Routines   *********************
 ********************************************************************/

/*@@
   @routine    CCTKi_SchedulePrintTimesFunction
   @date       Fri Oct 22 12:26:26 1999
   @author     Tom Goodale
   @desc
   Function which actually prints out data about a group or a function.
   @enddesc
   @calls

   @var     function
   @vdesc   the function to be called
   @vtype   void *
   @vio     in
   @endvar
   @var     attribute
   @vdesc   schedule item attributes
   @vtype   t_attribute *
   @vio     in
   @endvar
   @var     data
   @vdesc   data associated with schedule item
   @vtype   t_sched_data
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - this has no meaning
   @endreturndesc
@@*/
static int CCTKi_SchedulePrintTimesFunction(void *function,
                                            t_attribute *attribute,
                                            t_sched_data *data)
{
  /* prevent compiler warnings about unused parameters */
  function = function;

  /* find the timer for this function and this schedule bin */
  t_timer *timer = attribute->timers;
  while (timer && strcmp(timer->schedule_bin, data->schedule_bin))
  {
    timer = timer->next;
  }
  if (timer && ! timer->has_been_output)
  {
    timer->has_been_output = 1;
    CCTK_TimerI(timer->timer_handle, data->info);

    if(data->print_headers)
    {
      CCTKi_SchedulePrintTimerHeaders(data->info, data->file);

      data->print_headers = 0;
    }

    CCTKi_SchedulePrintTimerInfo(data->info, data->total_time,
                                 attribute->FunctionData.thorn,
                                 attribute->description,
                                 data->file);
    data->n_functions++;
  }

  return 1;
}

static void CCTKi_SchedulePrintTimerInfo(cTimerData *timer,
                                         cTimerData *total_time,
                                         const char *where,
                                         const char *description,
                                         FILE *file)
{
  int i, j;


  /* print delimiter line */
  if (*where == 0)
  {
    PrintDelimiterLine ('-', timer, file);
  }

  /* print the timer description */
  fprintf (file, "%-16.16s| %-40.40s", where, description);

  /* print the actual timer values */
  for (i = 0; i < timer->n_vals; i++)
  {
    j = strlen (timer->vals[i].heading) + strlen (timer->vals[i].units) + 3;

    fprintf (file, "| %*.8f ", j, timer->vals[i].seconds);
    if (total_time)
    {
      total_time->vals[i].seconds += timer->vals[i].seconds;
    }
  }

  fprintf (file, "\n");

  /* print delimiter line */
  if (*where == 0)
  {
    PrintDelimiterLine ('=', timer, file);
  }
}


static void CCTKi_SchedulePrintTimerHeaders (cTimerData *timer, FILE *file)
{
  int i;


  PrintDelimiterLine ('=', timer, file);

  fprintf (file,
           "%-16.16s| %-40.40s", "Thorn", "Scheduled routine in time bin");
  for (i = 0; i < timer->n_vals; i++)
  {
    fprintf (file, "| %s [%s] ", timer->vals[i].heading, timer->vals[i].units);
  }
  fprintf (file, "\n");

  PrintDelimiterLine ('=', timer, file);
}


static int CCTKi_ScheduleResetTimerOutputFlag(void *function,
                                              t_attribute *attribute,
                                              t_sched_data *data)
{
  /* prevent compiler warnings about unused parameters */
  function = function;
  data = data;

  /* find the timer for this function and this schedule bin */
  t_timer *timer = attribute->timers;
  while (timer)
  {
    timer->has_been_output = 0;
    timer = timer->next;
  }

  return 1;
}


static void PrintDelimiterLine (char delimiter,
                                const cTimerData *timer,
                                FILE *file)
{
  int i, len;


  len = 58;
  for (i = 0; i < timer->n_vals; i++)
  {
    len += strlen (timer->vals[i].heading) + strlen (timer->vals[i].units) + 6;
  }
  for (i = 0; i < len; i++)
  {
    fprintf (file, "%c", delimiter);
  }
  fprintf (file, "\n");
}
