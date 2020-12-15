 /*@@
   @file      ScheduleCreater.c
   @date      Tue Aug 31 12:46:08 1999
   @author    Tom Goodale
   @desc 
   
   @enddesc 
   @version $Header$
 @@*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cctk.h"
#include "cctk_WarnLevel.h"
#include "cctk_Parameters.h"
#include "cctk_Flesh.h"

#include "util_String.h"

#include "cctki_Schedule.h"
#include "StoreHandledData.h"
#include "Schedule.h"

static const char *rcsid="$Header$";

CCTK_FILEVERSION(schedule_ScheduleCreater_c);

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

static int ScheduleCreateGroup(const char *name);

static t_sched_item *ScheduleCreateItem(const char *name, 
                                        t_sched_modifier *modifiers, 
                                        void *attributes);

static int ScheduleAddItem(int ghandle, t_sched_item *item);

static int ScheduleSortGroup(t_sched_group *group);

static t_sched_modifier_type ScheduleTranslateModifierType(const char *modifier);

static const char * ScheduleModifierTypeName(t_sched_modifier_type type);

static int ScheduleItemNumber(t_sched_group *group, 
                              const char *name);

static t_sort_order ScheduleTranslateSortOrder(const char *order);

static int ScheduleCompareAscending(const void *a, const void *b);
static int ScheduleCompareDescending(const void *a, const void *b);

static int ScheduleSetupWhiles(t_sched_item *item);
static int ScheduleSetupIfs(t_sched_item *item);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

static int n_schedule_groups = 0;
static cHandledData *schedule_groups = NULL;


/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

 /*@@
   @routine    CCTKi_ScheduleAddModifer
   @date       Thu Sep  9 21:45:25 1999
   @author     Tom Goodale
   @desc 
   Adds a schedule modifier to a modifier list.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     orig
   @vdesc   original schedule modifier list
   @vtype   t_sched_modifier *
   @vio     inout
   @vcomment 
 
   @endvar 
   @var     modifier
   @vdesc   new modifier
   @vtype   const char *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     argument
   @vdesc   modifier argument
   @vtype   const char *
   @vio     in
   @vcomment 
 
   @returntype t_sched_modifier *
   @returndesc
   New schedule modifier list or NULL
   @endreturndesc
@@*/
t_sched_modifier *CCTKi_ScheduleAddModifier(t_sched_modifier *orig, 
                                            const char *modifier, 
                                            const char *argument)
{
  t_sched_modifier *this;

  this = (t_sched_modifier *)malloc(sizeof(t_sched_modifier));

  if(this)
  {
    this->argument = (char *)malloc((strlen(argument)+1)*sizeof(char));
    if(this->argument)
    {
      strcpy(this->argument, argument);

      this->type = ScheduleTranslateModifierType(modifier);

      this->next = orig;
    }
    else
    {
      free(this);
      this = NULL;
    }
  }

  return this;
}

 /*@@
   @routine    CCTKi_DoScheduleFunction
   @date       Thu Sep  9 21:42:58 1999
   @author     Tom Goodale
   @desc 
   Adds a function to a schedule group.  Creates the group if necessary.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     gname
   @vdesc   name of group to schedule function in
   @vtype   const char *
   @vio     in
   @vcomment 
   @var     fname
   @vdesc   name of function to be scheduled (may be aliased)
   @vtype   const char *
   @vio     in
   @vcomment 
   @var     function
   @vdesc   function to be scheduled
   @vtype   void *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     modifiers
   @vdesc   moodifier list
   @vtype   t_sched_modifier *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     attributes
   @vdesc   function attributes
   @vtype   void *
   @vio     in
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc 
    0 - success
   -1 - failure
   @endreturndesc
@@*/
int CCTKi_DoScheduleFunction(const char *gname, 
                             const char *fname, 
                             void *func, 
                             t_sched_modifier *modifiers, 
                             void *attributes)
{
  int retcode;
  int handle;
  t_sched_group *this_group;
  t_sched_item *newitem;

  handle = Util_GetHandle(schedule_groups, gname, (void **)&this_group);

  if(handle < 0)
  {
    handle = ScheduleCreateGroup(gname);
  }

  if(handle < 0)
  {
    retcode = -1;
  }
  else
  {
    newitem = ScheduleCreateItem(fname, modifiers, attributes);

    if(newitem)
    {
      newitem->type = sched_function;
      newitem->function = func;
      retcode = ScheduleAddItem(handle, newitem);
    }
    else
    {
      retcode = -1;
    }
  }

  return retcode;
}
    
 /*@@
   @routine    CCTKi_DoScheduleGroup
   @date       Thu Sep  9 21:43:44 1999
   @author     Tom Goodale
   @desc 
   Adds a group to a schedule group.  Creates the group if necessary.   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     gname
   @vdesc   name of group to schedule group in
   @vtype   const char *
   @vio     in
   @var     thisname
   @vdesc   name of group to be scheduled (possibly an alias)
   @vtype   const char *
   @vio     in
   @endvar 
   @var     realname
   @vdesc   actual name of group to be scheduled
   @vtype   const char *
   @vio     in
   @endvar 
   @var     modifiers
   @vdesc   moodifier list
   @vtype   t_sched_modifier *
   @vio     in
   @endvar 
   @var     attributes
   @vdesc   function attributes
   @vtype   void *
   @vio     in
   @endvar 

   @returntype int
   @returndesc 
    0 - success
   -1 - failure
   @endreturndesc
@@*/
int CCTKi_DoScheduleGroup(const char *gname, 
                          const char *thisname,
                          const char *realname,
                          t_sched_modifier *modifiers, 
                          void *attributes)
{
  int retcode;
  int handle;
  int thishandle;
  t_sched_group *this_group;
  t_sched_item *newitem;

  /* Find the group within which to schedule this group */
  handle = Util_GetHandle(schedule_groups, gname, (void **)&this_group);

  if(handle < 0)
  {
    handle = ScheduleCreateGroup(gname);
  }

  /* Find this group */
  thishandle = Util_GetHandle(schedule_groups, realname, (void **)&this_group);

  if(thishandle < 0)
  {
    thishandle = ScheduleCreateGroup(realname);
  }

  if(handle < 0 || thishandle < 0)
  {
    retcode = -1;
  }
  else
  {
    newitem = ScheduleCreateItem(thisname, modifiers, attributes);

    if(newitem)
    {
      newitem->type = sched_group;
      newitem->group = thishandle;
      retcode = ScheduleAddItem(handle, newitem);
    }
    else
    {
      retcode = -1;
    }
  }

  return retcode;
}

 /*@@
   @routine    CCTKi_DoScheduleSortAllGroups
   @date       Wed Sep 15 22:37:49 1999
   @author     Tom Goodale
   @desc 
   Sorts all the schedule groups.
   @enddesc 
   @calls   ScheduleSortGroups  
   @calledby   
   @history 
    @hdate    Mon Mar 21 16:33:07 CET 2005
    @hauthor  Jonathan Thornburg <jthorn@aei.mpg.de>
    @nddesc   change error-handling strategy from
              "print a message to stderr, then return error code"
              to CCTK_VWarn(CCTK_WARN_ABORT, ...), i.e.
              "print a message to stderr, then abort the Cactus run";
              this is a partial fix for bug Cactus/1908
   @endhistory 

   @returntype int
   @returndesc 
    0  - success
   -ve - -1* number of errors
   @endreturndesc
@@*/
int CCTKi_DoScheduleSortAllGroups(void)
{
  int group;
  t_sched_group *gdata;
  int errcode;
  int n_errors;
  
  n_errors  = 0;

  for(group = 0; group < n_schedule_groups; group++)
  {
    if((gdata = (t_sched_group *)Util_GetHandledData(schedule_groups, group)))
    {
      errcode = ScheduleSortGroup(gdata);

      if(errcode)
      {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, "Cactus",
                "Error while sorting group '%s' - %d remaining unsorted routines.\n", 
                gdata->name,
                -errcode);                                      /*NOTREACHED*/

        n_errors += -errcode;
      }
    }
  }
   
  return -n_errors;
}

 /*@@
   @routine    CCTKi_DoScheduleGetGroups
   @date       Wed Sep 15 22:37:49 1999
   @author     Tom Goodale
   @desc 
   Gets the schedule groups
   @enddesc 
   @calls   ScheduleSortGroups  
   @calledby   
   @history 
 
   @endhistory 

   @returntype cHandledData
   @returndesc 
   The scheduled groups.
   @endreturndesc
@@*/
cHandledData *CCTKi_DoScheduleGetGroups(void)
{
  return schedule_groups;
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

 /*@@
   @routine    ScheduleCreateGroup
   @date       Wed Sep  8 11:15:32 1999
   @author     Tom Goodale
   @desc 
   Creates a schedule group.
   @enddesc 
   @calls     Util_GetHandle
   @calledby   
   @history 
 
   @endhistory 
   @var     name
   @vdesc   name of the group
   @vtype   const char *
   @vio     in
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc 
    0 - success
   -1 - group already exists
   -2 - memory failure
   @endreturndesc
@@*/
static int ScheduleCreateGroup(const char *name)
{
  int retcode;
  int handle;

  t_sched_group *this_group;

  handle = Util_GetHandle(schedule_groups, name, (void **)&this_group);

  if(handle > -1)
  {
    /* Group already exists */
    retcode = -1;
  }
  else
  {
    this_group = (t_sched_group *)malloc(sizeof(t_sched_group));

    if(this_group)
    {
      this_group->name = (char *)malloc((strlen(name)+1)*sizeof(char));

      if(this_group->name)
      {
        strcpy(this_group->name, name);

        this_group->order = NULL;
        this_group->n_scheditems = 0;
        this_group->scheditems = NULL;
        retcode = Util_NewHandle(&schedule_groups, name, (void *)this_group);
        n_schedule_groups++;
      }
      else
      {
        free(this_group);

        retcode = -2;
      }
    }
    else
    {
      retcode = -2;
    }
  }

  return retcode;
}

 /*@@
   @routine    ScheduleCreateItem
   @date       Thu Sep  9 21:44:17 1999
   @author     Tom Goodale
   @desc 
   Creates a schedule item to be scheduled.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     name
   @vdesc   name of item
   @vtype   const char *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     modifiers
   @vdesc   modifier list
   @vtype   t_sched_modifier *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     attributes
   @vdesc   item attributes
   @vtype   void *
   @vio     in
   @vcomment 
 
   @endvar 

   @returntype t_sched_item *
   @returndesc 
   The new schedule item or NULL.
   @endreturndesc
@@*/
static t_sched_item *ScheduleCreateItem(const char *name, t_sched_modifier *modifiers, void *attributes)
{
  t_sched_item *this;

  this = (t_sched_item *)malloc(sizeof(t_sched_item));

  if(this)
  {
    this->name = (char *)malloc((strlen(name)+1)*sizeof(char));

    if(this->name)
    {
      strcpy(this->name, name);

      this->type     = sched_item_none;
      this->function = NULL;
      this->group    = -1;
      this->modifiers = modifiers;

      this->n_whiles = 0;
      this->whiles = NULL;
      this->n_ifs = 0;
      this->ifs = NULL;

      ScheduleSetupWhiles(this);
      ScheduleSetupIfs(this);

      this->attributes = attributes;

#ifdef DEBUG_SCHEDULAR
      printf("Created Schedule item %s\n", this->name);
#endif

    }
    else
    {
      free(this);
      this = NULL;
    }
  }

  return this;
}

 /*@@
   @routine    ScheduleAddItem
   @date       Thu Sep  9 21:45:03 1999
   @author     Tom Goodale
   @desc 
   Adds a schedule item to a group.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     ghandle
   @vdesc   The handle of the group
   @vtype   int
   @vio     in
   @vcomment 
 
   @endvar 
   @var     item
   @vdesc   The schedule item
   @vtype   t_sched_item *
   @vio     in
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc 
    0 - success
   -1 - memory failure
   -2 - duplicate item
   @endreturndesc
@@*/
static int ScheduleAddItem(int ghandle, t_sched_item *item)
{
  int retcode;
  t_sched_group *this_group;
  t_sched_item *temp;

  this_group = (t_sched_group *)Util_GetHandledData(schedule_groups, ghandle);

  if(ScheduleItemNumber(this_group, item->name) == -1)
  {
    this_group->n_scheditems++;

    temp = (t_sched_item *)realloc(this_group->scheditems, this_group->n_scheditems*sizeof(t_sched_item));

    if(temp)
    {
      this_group->scheditems = temp;
      this_group->scheditems[this_group->n_scheditems-1] = *item;

#ifdef DEBUG_SCHEDULAR
      printf("Added item '%s' to group '%s'\n", item->name, this_group->name);
#endif

      free(item);

      retcode = 0;
    }
    else
    {
      this_group->n_scheditems--;
      retcode = -1;
    }
  }
  else
  {
    /* Item already existed. */
    retcode = -2;
  }

  return retcode;
}
  


 /*@@
   @routine    ScheduleTranslateModifierType
   @date       Thu Sep  9 21:45:56 1999
   @author     Tom Goodale
   @desc 
   Translates a modifier type from a string to an enum.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     modifier
   @vdesc   The modifier string
   @vtype   const char *
   @vio     in
   @vcomment 
 
   @endvar 

   @returntype t_sched_modifier_type
   @returndesc 
   The enumerated schedule modifier type.
   @endreturndesc
@@*/
static t_sched_modifier_type ScheduleTranslateModifierType(const char *modifier)
{
  /* FIXME */
  
  t_sched_modifier_type retval;

  retval = sched_mod_none;

  if(!strcmp(modifier, "before"))
  {
    retval = sched_before;
  }
  else if(!strcmp(modifier, "after"))
  {
    retval = sched_after;
  }
  else if(!strcmp(modifier, "while"))
  {
    retval = sched_while;
  }
  else if(!strcmp(modifier, "if"))
  {
    retval = sched_if;
  }

#ifdef DEBUG_SCHEDULAR
  printf("Translated modifier type %s to %d\n", modifier, retval);
#endif

  return retval;
}

 /*@@
   @routine    ScheduleTranslateSortOrder
   @date       Sun Jul 15 20:29:29 PDT 2012
   @author     Roland Haas
   @desc 
   Translates a sorting order type from a string to an enum.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     order
   @vdesc   The order string
   @vtype   const char *
   @vio     in
   @vcomment 
 
   @endvar 

   @returntype t_sort_order
   @returndesc 
   The enumerated sort order.
   @endreturndesc
@@*/
static t_sort_order ScheduleTranslateSortOrder(const char *order)
{
  t_sort_order retval;

  if(CCTK_Equals(order, "none"))
  {
    retval = sort_order_none;
  }
  else if(CCTK_Equals(order, "ascending"))
  {
    retval = sort_order_ascending;
  }
  else if(CCTK_Equals(order, "descending"))
  {
    retval = sort_order_descending;
  }
  else
  {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, "Cactus",
        "Unknown sort ordering type '%s'", order);
    retval = -1; /* NOTREACHED */
  }

#ifdef DEBUG_SCHEDULAR
  printf("Translated sort order %s to %d\n", order, retval);
#endif

  return retval;
}


 /*@@
   @routine    ScheduleModifierTypeName
   @date       Sun Jul 15 18:42:48 PDT 2012
   @author     Roland Haas
   @desc 
   Translates a modifier type from an enum to a string.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     type
   @vdesc   The modifier type
   @vtype   t_sched_modifier_type
   @vio     in
   @vcomment 
 
   @endvar 

   @returntype const char *
   @returndesc 
   The named schedule modifier type.
   @endreturndesc
@@*/
static const char *ScheduleModifierTypeName(t_sched_modifier_type type)
{
  switch(type)
  {
    case sched_before:
      return "before";
      break;
    case sched_after:
      return "after";
      break;
    case sched_while:
      return "while";
      break;
    case sched_if:
      return "if";
      break;
    case sched_mod_none:
      /* FALLTHROUGH */
    default:
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, "Cactus",
          "Internal error: Unknown schedule modifier type %d",
          (int)type);
      break;
  }
  return NULL; /* NOTREACHED */
}

 /*@@
   @routine    ScheduleSortGroup
   @date       Mon Sep 13 11:30:19 1999
   @author     Tom Goodale
   @desc 
   Sorts the routines in a group.
   @enddesc 
   @calls     CCTKi_ScheduleCreateArray CCTKi_ScheduleCreateIVec
              ScheduleItemNumber CCTKi_ScheduleAddRow
              CCTKi_ScheduleSort
              CCTKi_ScheduleDestroyArray CCTKi_ScheduleDestroyIVec
   @calledby   
   @history 
    @hdate    Mon Mar 21 16:33:07 CET 2005
    @hauthor  Jonathan Thornburg <jthorn@aei.mpg.de>
    @nddesc   change error-handling strategy from
              "print a message to stderr, then return error code"
              to CCTK_VWarn(CCTK_WARN_ABORT, ...), i.e.
              "print a message to stderr, then abort the Cactus run"
              this is a partial fix for bug Cactus/1908
   @endhistory 
   @var     group
   @vdesc   The schedule group
   @vtype   t_sched_group *
   @vio     inout
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc 
    0 - success
    Or number of scheduling errors
   @endreturndesc
@@*/
static int ScheduleSortGroup(t_sched_group *group)
{
  DECLARE_CCTK_PARAMETERS;

  int item;
  int *order;
  int *thisorders;
  t_sched_modifier *modifier;
  signed char **array;
  int number;
  int mod;
  int i;
  int errcode;
  int sort_order;

#ifdef DEBUG_SCHEDULAR
  int j;
#endif

  /* Create the data staructures */
  array      = CCTKi_ScheduleCreateArray(group->n_scheditems);
  order      = CCTKi_ScheduleCreateIVec(group->n_scheditems);
  thisorders = CCTKi_ScheduleCreateIVec(group->n_scheditems);

  sort_order = ScheduleTranslateSortOrder(schedule_sort_mode);
  if(sort_order == sort_order_ascending)
  {
    qsort(group->scheditems, group->n_scheditems,
        sizeof(*group->scheditems), ScheduleCompareAscending);
  }
  else if(sort_order == sort_order_descending)
  {
    qsort(group->scheditems, group->n_scheditems,
        sizeof(*group->scheditems), ScheduleCompareDescending);
  }
  
  for(item=0; item < group->n_scheditems; item++)
  {
#ifdef DEBUG_SCHEDULAR
    printf("Scheduling item %d '%s'\n", item, group->scheditems[item].name);
#endif
    for(modifier = group->scheditems[item].modifiers; modifier; modifier = modifier->next)
    {
      if(modifier->type == sched_while || modifier->type == sched_if)
      {
        continue;
      }
      number = ScheduleItemNumber(group, modifier->argument);

      if(number < 0 && schedule_sort_warnings)
      {
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
            "Modifier '%s' of schedule item '%s' in group '%s' refers to non-existing item '%s'.",
            ScheduleModifierTypeName(modifier->type), group->scheditems[item].name,
            group->name, modifier->argument);
      }

      if(number >= 0 && number < group->n_scheditems)
      {
#ifdef DEBUG_SCHEDULAR
        printf("Scheduling against item %d '%s' - mod-type %d\n", number, group->scheditems[number].name, modifier->type);
#endif
        switch(modifier->type)
        {
          case sched_before : mod = -1;  break;
          case sched_after  : mod = 1; break;
          case sched_while: /* FALLTHROUGH */
          case sched_if: /* FALLTHROUGH */
          case sched_mod_none: /* FALLTHROUGH */
          default :
            mod = 0;
        }
#ifdef DEBUG_SCHEDULAR
        printf("Modifier is %d\n", mod);
#endif

        thisorders[number] = mod;
      }
    }

#ifdef DEBUG_SCHEDULAR
    printf("Orderlist for item %d is...\n", item);
    for(i=0; i < group->n_scheditems; i++)
    {
      printf("  %d", thisorders[i]);
    }
    printf("\n");
#endif

    errcode = CCTKi_ScheduleAddRow(group->n_scheditems, array, order, item, thisorders);
    if(errcode)
    {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, "Cactus",
      "Adding item %s to group %s failed due to a circular dependency on %s.\n",
      group->scheditems[item].name, group->name,
      group->scheditems[-errcode-1].name );  /*NOTREACHED*/
    }

    /* Clear the array for the next item. */
    for(i=0; i < group->n_scheditems; i++)
    {
      thisorders[i] = 0;
    }
  }

#ifdef DEBUG_SCHEDULAR
  printf("Initial array is...\n");
  for(i=0; i < group->n_scheditems; i++)
  {
    for(j=0; j < group->n_scheditems; j++)
    {
      printf("  %d", (int)array[i][j]);
    }

    printf("\n");
  }

  printf("Initial order is...\n");
  for(i=0; i < group->n_scheditems; i++)
  {
    printf("  %d", order[i]);
  }
  printf("\n");
  
  printf("Sorting array...\n");
#endif

  errcode = CCTKi_ScheduleSort(group->n_scheditems, array, order);

  if(errcode)
  {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, "Cactus",
    "Schedule sort of group %s failed with error code %d\n",
    group->name, errcode);                                  /*NOTREACHED*/
  }

#ifdef DEBUG_SCHEDULAR
  printf("Final array is...\n");
  for(i=0; i < group->n_scheditems; i++)
  {
    for(j=0; j < group->n_scheditems; j++)
    {
      printf("  %d", (int)array[i][j]);
    }

    printf("\n");
  }

  printf("Final order is...\n");
  for(i=0; i < group->n_scheditems; i++)
  {
    printf("  %d", order[i]);
  }
  printf("\n");
#endif

  /* Free memory */
  CCTKi_ScheduleDestroyIVec(group->n_scheditems,thisorders);
  CCTKi_ScheduleDestroyArray(group->n_scheditems, array);

  group->order = order;

  return errcode;
}

 /*@@
   @routine    ScheduleItemNumber
   @date       Mon Sep 13 11:30:49 1999
   @author     Tom Goodale
   @desc 
   Returns the number of a specific item in the array of schedule items.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     group
   @vdesc   schedule group
   @vtype   t_sched_group *
   @vio     in
   @vcomment 
 
   @endvar 
   @var     name
   @vdesc   name of schedule item
   @vtype   const char *
   @vio     in
   @vcomment 
 
   @endvar 

   @returntype int
   @returndesc 
   The number of the schedule item in the group
   @endreturndesc
@@*/
static int ScheduleItemNumber(t_sched_group *group, const char *name)
{
  int retval;
  int i;
  
  retval = -1;
  /* This is the quick-to-write way to do this.  Should really put into a
   * hash table or a tree, but this will do for the moment.
   */
  for(i = 0 ; i < group->n_scheditems; i++)
  {
    if(!strcmp(group->scheditems[i].name, name))
    {
      retval = i;
      break;
    }
  }

  return retval;
}


 /*@@
   @routine    ScheduleSetupWhiles
   @date       Wed Sep 15 20:10:28 1999
   @author     Tom Goodale
   @desc 
   Make an array of all the whiles in the modifier list for a schedule item.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     item
   @vdesc   The schedule item to work on
   @vtype   t_sched_item *
   @vio     inout
   @vcomment 
 
   @endvar 
   @returntype int
   @returndesc 
   Number of whiles
   @endreturndesc
@@*/
static int ScheduleSetupWhiles(t_sched_item *item)
{
  int retval;
  t_sched_modifier *modifier;
  char **temp;

  retval = 0;

  for(modifier = item->modifiers; modifier; modifier = modifier->next)
  {
    if(modifier->type == sched_while)
    {
      item->n_whiles++;
      temp = (char **)realloc(item->whiles, item->n_whiles*sizeof(char *));

      if(temp)
      {
        item->whiles = temp;

        temp[item->n_whiles-1] = (char *)malloc((strlen(modifier->argument)+1)*sizeof(char));
        if(temp[item->n_whiles-1])
        {
          strcpy(temp[item->n_whiles-1], modifier->argument);
        }
        else
        {
          item->n_whiles--;
          retval--;
        }
      }
      else
      {
        retval--;
      }
    }
  }

  return retval;
}

 /*@@
   @routine    ScheduleSetupIfs
   @date       Dec 27, 2005
   @author     Erik Schnetter
   @desc 
   Make an array of all the ifs in the modifier list for a schedule item.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 
   @var     item
   @vdesc   The schedule item to work on
   @vtype   t_sched_item *
   @vio     inout
   @vcomment 
 
   @endvar 
   @returntype int
   @returndesc 
   Number of whiles
   @endreturndesc
@@*/
static int ScheduleSetupIfs(t_sched_item *item)
{
  int retval;
  t_sched_modifier *modifier;
  char **temp;

  retval = 0;

  for(modifier = item->modifiers; modifier; modifier = modifier->next)
  {
    if(modifier->type == sched_if)
    {
      item->n_ifs++;
      temp = (char **)realloc(item->ifs, item->n_ifs*sizeof(char *));

      if(temp)
      {
        item->ifs = temp;

        temp[item->n_ifs-1] = (char *)malloc((strlen(modifier->argument)+1)*sizeof(char));
        if(temp[item->n_ifs-1])
        {
          strcpy(temp[item->n_ifs-1], modifier->argument);
        }
        else
        {
          item->n_ifs--;
          retval--;
        }
      }
      else
      {
        retval--;
      }
    }
  }

  return retval;
}

 /*@@
   @routine    ScheduleCompareAscending
   @date       Mon Jul 16 07:27:39 PDT 2012
   @author     Roland Haas
   @desc
   Compare two schedule items based on their thorn and routine names
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
   @var     a
   @vdesc   The first schedule items to work on
   @vtype   t_sched_item *
   @vio     ino
   @vcomment

   @var     b
   @vdesc   The first schedule items to work on
   @vtype   t_sched_item *
   @vio     ino
   @vcomment

   @returntype int
   @returndesc
   Returns +strcmpi(a->name, b->name)
   @endreturndesc
@@*/
static int ScheduleCompareAscending(const void *a_, const void *b_)
{
  const t_sched_item *a = a_;
  const t_sched_item *b = b_;
  return +Util_StrCmpi(a->name, b->name);
}

 /*@@
   @routine    ScheduleCompareDescending
   @date       Mon Jul 16 07:27:39 PDT 2012
   @author     Roland Haas
   @desc
   Compare (in reverse alphabetical order) two schedule items based on their
   thorn and routine names
   @enddesc
   @calls
   @calledby
   @history

   @endhistory
   @var     a
   @vdesc   The first schedule items to work on
   @vtype   t_sched_item *
   @vio     ino
   @vcomment

   @var     b
   @vdesc   The first schedule items to work on
   @vtype   t_sched_item *
   @vio     ino
   @vcomment

   @returntype int
   @returndesc
   Returns -strcmpi(a->name, b->name)
   @endreturndesc
@@*/
static int ScheduleCompareDescending(const void *a_, const void *b_)
{
  const t_sched_item *a = a_;
  const t_sched_item *b = b_;
  return -Util_StrCmpi(a->name, b->name);
}

/********************************************************************
 ********************************************************************
 ********************************************************************/

#ifdef TEST_SCHEDULECREATOR

#include <stdarg.h>

#define func_x(x) \
int func_ ## x (void) { return printf("I'm func " #x "\n"); }

int CCTK_VWarn(int level, int line, const char *file, 
    const char *thorn, const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);

  fprintf(stdout, "WARNING[L%d,P0] (%s) in line %d of %s: ", level, thorn, line, file);
  vfprintf(stdout, fmt, ap);
  fputc('\n', stdout);

  va_end(ap);

  if(level == 0)
    exit(1);

  return 0;
}

func_x(a)
func_x(b)
func_x(c)

#define func_x_proto(x) int func_ ## x (void);

int main(int argc, char *argv[])
{
  t_sched_modifier *modifier;

  modifier = CCTKi_ScheduleAddModifier(NULL, "before", "c");
  modifier = CCTKi_ScheduleAddModifier(modifier, "after",  "a");

  CCTKi_DoScheduleFunction("group_a", "c", func_c, NULL, NULL);
  CCTKi_DoScheduleFunction("group_a", "b", func_b, modifier, NULL);
  CCTKi_DoScheduleFunction("group_a", "a", func_a, NULL, NULL);
  CCTKi_DoScheduleFunction("group_b", "a", func_a, NULL, NULL);
  CCTKi_DoScheduleFunction("group_b", "b", func_b, NULL, NULL);
  CCTKi_DoScheduleGroup("group_a", "group_b", modifier, NULL);

  CCTKi_DoScheduleSortAllGroups();

  return 0;
}
#endif
