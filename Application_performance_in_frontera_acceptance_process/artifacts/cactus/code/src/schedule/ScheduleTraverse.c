 /*@@
   @file      ScheduleTraverse.c
   @date      Thu Sep 16 08:58:37 1999
   @author    Tom Goodale
   @desc 
   Routins to traverse schedule groups.
   @enddesc 
   @version $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>

#include "cctk_Flesh.h"
#include "cctk_WarnLevel.h"

#include "cctki_Schedule.h"
#include "StoreHandledData.h"
#include "Schedule.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(schedule_ScheduleTraverse_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int ScheduleTraverseGroup(cHandledData *schedule_groups, 
                                 t_sched_group *group,
                                 void *attributes,
                                 int n_whiles,
                                 char **whiles,
                                 int n_ifs,
                                 char **ifs,
                                 int (*item_entry)(void *, void *),
                                 int (*item_exit)(void *, void *),
                                 int (*while_check)(int, char **, void *, void *, int),
                                 int (*if_check)(int, char **, void *, void *, int),
                                 int (*function_process)(void *, void *, void *),
                                 void *data);

static int ScheduleTraverseFunction(void *function,
                                    void *attributes,
                                    int n_whiles,
                                    char **whiles,
                                    int n_ifs,
                                    char **ifs,
                                    int (*item_entry)(void *, void *),
                                    int (*item_exit)(void *, void *),
                                    int (*while_check)(int, char **, void *, void *, int),
                                    int (*if_check)(int, char **, void *, void *, int),
                                    int (*function_process)(void *, void *, void *),
                                    void *data);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

 /*@@
   @routine    CCTKi_DoScheduleTraverse
   @date       Thu Sep 16 09:05:23 1999
   @author     Tom Goodale
   @desc 
   Traverses the group with the given name.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     group_name
   @vdesc   group to traverse
   @vtype   const char *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     item_entry
   @vdesc   function to be called on entry to an item
   @vtype   int (*)(void *, void *)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     item_exit
   @vdesc   function to be called on exit from an item
   @vtype   int (*)(void *, void *)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     while_check
   @vdesc   function to be called to check a while statement
   @vtype   int (*)(int, char **, void *, void *, int)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     if_check
   @vdesc   function to be called to check an if statement
   @vtype   int (*)(int, char **, void *, void *, int)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     function_process
   @vdesc   function to be called on any function
   @vtype   int (*)(void *, void *, void *)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     data
   @vdesc   data to be passed to the functions
   @vtype   void *
   @vio     inout
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc 
   Return code of @seeroutine ScheduleTraverseGroup or return value of 
   @seeroutine Util_GetHandle if that fails.
   @endreturndesc
@@*/
int CCTKi_DoScheduleTraverse(const char *group_name,
                             int (*item_entry)(void *, void *),
                             int (*item_exit)(void *, void *),
                             int (*while_check)(int, char **, void *, void *, int),
                             int (*if_check)(int, char **, void *, void *, int),
                             int (*function_process)(void *, void *, void *),
                             void *data)
{
  cHandledData *schedule_groups;
  t_sched_group *group;
  int handle;
  int retcode;

  schedule_groups = CCTKi_DoScheduleGetGroups();

  handle = Util_GetHandle(schedule_groups, group_name, (void *)&group);

  if(handle >= 0)
  {
    retcode = ScheduleTraverseGroup(schedule_groups, 
                                    group, 
                                    NULL,
                                    0,
                                    NULL,
                                    0,
                                    NULL,
                                    item_entry, 
                                    item_exit, 
                                    while_check, 
                                    if_check, 
                                    function_process,
                                    data);
  }
  else
  {
    retcode = handle;
  }

  return retcode;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    ScheduleTraverseGroup
   @date       Thu Sep 16 09:07:44 1999
   @author     Tom Goodale
   @desc 
   Traverses the given schedule group.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     schedule_groups
   @vdesc   the schedule groups
   @vtype   cHandledData
   @vio     in
   @vcomment 
 
   @endvar 
   @var     group
   @vdesc   the group to traverse
   @vtype   t_sched_group
   @vio     in
   @vcomment 
 
   @endvar 
   @var     attributes
   @vdesc   group attributes
   @vtype   void *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     n_whiles
   @vdesc   number of whiles
   @vtype   int
   @vio     in
   @vcomment 
 
   @endvar 
   @var     whiles
   @vdesc   array of while strings
   @vtype   char **
   @vio     in
   @vcomment 
 
   @endvar 
   @var     n_ifs
   @vdesc   number of ifs
   @vtype   int
   @vio     in
   @vcomment 
 
   @endvar 
   @var     ifs
   @vdesc   array of if strings
   @vtype   char **
   @vio     in
   @vcomment 
 
   @endvar 
   @var     item_entry
   @vdesc   function to be called on entry to an item
   @vtype   int (*)(void *, void *)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     item_exit
   @vdesc   function to be called on exit from an item
   @vtype   int (*)(void *, void *)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     while_check
   @vdesc   function to be called to check a while statement
   @vtype   int (*)(int, char **, void *, void *, int)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     if_check
   @vdesc   function to be called to check a if statement
   @vtype   int (*)(int, char **, void *, void *, int)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     function_process
   @vdesc   function to be called on any function
   @vtype   int (*)(void *, void *, void *)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     data
   @vdesc   data to be passed to the functions
   @vtype   void *
   @vio     inout
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc 
   0 - success
   @endreturndesc
@@*/
static int ScheduleTraverseGroup(cHandledData *schedule_groups, 
                                 t_sched_group *group,
                                 void *attributes,
                                 int n_whiles,
                                 char **whiles,
                                 int n_ifs,
                                 char **ifs,
                                 int (*item_entry)(void *, void *),
                                 int (*item_exit)(void *, void *),
                                 int (*while_check)(int, char **, void *, void *, int),
                                 int (*if_check)(int, char **, void *, void *, int),
                                 int (*function_process)(void *, void *, void *),
                                 void *data)
{
  int item;
  int doit;
  int called_item_entry;
  t_sched_group *newgroup;

  doit = 1;

  /* If there is an if-list associated with this item, check if the
   * group should be exectuted at all.
   */

  if(n_ifs > 0 && if_check)
  {
    doit = doit && if_check(n_ifs, ifs, attributes, data, 1);
  }

  /* If there is a while-list associated with this item, check if the
   * group should be exectuted at all.
   *
   * If there are both an if-list and a while-list, then both are
   * checked when the group is entered, but only the while-list is
   * checked to determine whether the loop should continue.  This
   * means that the if-statement is outside the while-statement, as in
   *
   * IF (...)
   *   WHILE (...)
   *     schedule stuff
   *   END WHILE
   * END IF
   */

  if(n_whiles > 0 && while_check)
  {
    doit = doit && while_check(n_whiles, whiles, attributes, data, 1);
  }

  /* Call a item entry function if it is defined. */
  if(doit)
  {
    called_item_entry = 1;

    if(item_entry)
    {
      doit = item_entry(attributes, data);
    }
  }
  else
  {
    called_item_entry = 0;
  }

  /* Now traverse the group. */
  while(doit)
  {
      
    /* Traverse in the sorted order - assumes group has been sorted ! */
    for(item = 0 ; item < group->n_scheditems; item++)
    {
      switch(group->scheditems[group->order[item]].type)
      {
        case sched_function :
          ScheduleTraverseFunction(group->scheditems[group->order[item]].function, 
                                   group->scheditems[group->order[item]].attributes,
                                   group->scheditems[group->order[item]].n_whiles,
                                   group->scheditems[group->order[item]].whiles,
                                   group->scheditems[group->order[item]].n_ifs,
                                   group->scheditems[group->order[item]].ifs,
                                   item_entry,
                                   item_exit,
                                   while_check,
                                   if_check,
                                   function_process,
                                   data);
          break;
        case sched_group :
          newgroup = (t_sched_group *)Util_GetHandledData(schedule_groups, 
                                                          group->scheditems[group->order[item]].group);
          ScheduleTraverseGroup(schedule_groups,
                                newgroup, 
                                group->scheditems[group->order[item]].attributes,
                                group->scheditems[group->order[item]].n_whiles,
                                group->scheditems[group->order[item]].whiles,
                                group->scheditems[group->order[item]].n_ifs,
                                group->scheditems[group->order[item]].ifs,
                                item_entry, 
                                item_exit, 
                                while_check, 
                                if_check, 
                                function_process,
                                data);
          break;
        case sched_item_none :
          /* FALLTHROUGH */
        default :
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                     "Unknown schedule item type %d\n", group->scheditems[group->order[item]].type);
          break;
      }
    }

    /* Check the while_list again. */
    if(n_whiles > 0 && while_check)
    {
      doit = while_check(n_whiles, whiles, attributes, data, 0);
    }
    else
    {
      doit = 0;
    }
  }

  /* Call the group_exit function if it's defined. */
  if(called_item_entry)
  {
    if(item_exit)
    {
      item_exit(attributes, data);
    }
  }

  /* this is used for SchedulePrint to reduce the indentation level */
  if(n_ifs > 0 && if_check)
  {
    if_check(n_ifs, ifs, attributes, data, 0);
  }

  return 0;
}

 /*@@
   @routine    ScheduleTraverseFunction
   @date       Thu Sep 16 11:51:58 1999
   @author     Tom Goodale
   @desc 
   Deals with a function in the schedule list
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     function
   @vdesc   the function to be called
   @vtype   void *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     attributes
   @vdesc   function attributes
   @vtype   void *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     n_whiles
   @vdesc   number of whiles
   @vtype   int
   @vio     in
   @vcomment 
 
   @endvar 
   @var     whiles
   @vdesc   array of while strings
   @vtype   char **
   @vio     in
   @vcomment 
 
   @endvar 
   @var     n_ifs
   @vdesc   number of ifs
   @vtype   int
   @vio     in
   @vcomment 
 
   @endvar 
   @var     ifs
   @vdesc   array of if strings
   @vtype   char **
   @vio     in
   @vcomment 
 
   @endvar 
   @var     item_entry
   @vdesc   function to be called on entry to an item
   @vtype   int (*)(void *, void *)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     item_exit
   @vdesc   function to be called on exit from an item
   @vtype   int (*)(void *, void *)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     while_check
   @vdesc   function to be called to check a while statement
   @vtype   int (*)(int, char **, void *, void *, int)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     if_check
   @vdesc   function to be called to check a if statement
   @vtype   int (*)(int, char **, void *, void *, int)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     function_process
   @vdesc   function to be called on any function
   @vtype   int (*)(void *, void *, void *)
   @vio     in
   @vcomment 
 
   @endvar 
   @var     data
   @vdesc   data to be passed to the functions
   @vtype   void *
   @vio     inout
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc 
   0 - success
   @endreturndesc
@@*/
static int ScheduleTraverseFunction(void *function,
                                    void *attributes,
                                    int n_whiles,
                                    char **whiles,
                                    int n_ifs,
                                    char **ifs,
                                    int (*item_entry)(void *, void *),
                                    int (*item_exit)(void *, void *),
                                    int (*while_check)(int, char **, void *, void *, int),
                                    int (*if_check)(int, char **, void *, void *, int),
                                    int (*function_process)(void *, void *, void *),
                                    void *data)
{
  int doit;
  int called_item_entry;

  doit = 1;

  /* If there is an if-list associated with this function, check if
   * the function should be executed at all.
   */

  if(n_ifs > 0 && if_check)
  {
    doit = doit && if_check(n_ifs, ifs, attributes, data, 1);
  }

  /* If there is a while-list associated with this function, check if
   * the function should be executed at all.
   */

  if(n_whiles > 0 && while_check)
  {
    doit = doit && while_check(n_whiles, whiles, attributes, data, 1);
  }

  /* Call a item entry function if it is defined. */
  if(doit)
  {
    called_item_entry = 1;

    if(item_entry)
    {
      doit = item_entry(attributes, data);
    }
  }
  else
  {
    called_item_entry = 0;
  }

  /* Now traverse the . */
  while(doit)
  {
    
    /* Now actually do something with the function. */
    function_process(function, attributes, data);

    /* Check the while_list again. */
    if(n_whiles > 0 && while_check)
    {
      doit = while_check(n_whiles, whiles, attributes, data, 0);
    }
    else
    {
      doit = 0;
    }
  }

  /* Call the item_exit function if it's defined. */
  if(called_item_entry)
  {
    if(item_exit)
    {
      item_exit(attributes, data);
    }
  }

  /* this is used for SchedulePrint to reduce the indentation level */
  if(n_ifs > 0 && if_check)
  {
    if_check(n_ifs, ifs, attributes, data, 0);
  }

  return 0;
}

/********************************************************************
 ********************************************************************
 ********************************************************************/

#ifdef TEST_SCHEDULETRAVERSE

#define func_x(x) \
int func_ ## x (void) { return printf("I'm func " #x "\n"); }

func_x(a)
func_x(b)
func_x(c)

int fprocess(void *function, void *attributes, void *data)
{
  int (*func)(void);

  func = (int (*)(void)) function;

  func();

  return 1;
}

int main(int argc, char *argv[])
{
  t_sched_modifier *modifier;

  modifier = CCTKi_DoScheduleAddModifer(NULL, "before", "c");
  modifier = CCTKi_DoScheduleAddModifer(modifier, "after",  "a");

  CCTKi_DoScheduleFunction("group_a", "c", func_c, NULL, NULL);
  CCTKi_DoScheduleFunction("group_a", "b", func_b, modifier, NULL);
  CCTKi_DoScheduleFunction("group_a", "a", func_a, NULL, NULL);
  CCTKi_DoScheduleFunction("group_b", "a", func_a, NULL, NULL);
  CCTKi_DoScheduleFunction("group_b", "b", func_b, NULL, NULL);
  CCTKi_DoScheduleGroup("group_a", "group_b", modifier, NULL);

  CCTKi_DoScheduleSortAllGroups();

  CCTKi_DoScheduleTraverse("group_a", NULL, NULL, NULL, fprocess, NULL);

  return 0;
}
#endif
