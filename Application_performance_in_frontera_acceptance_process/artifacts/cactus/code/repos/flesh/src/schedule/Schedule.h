 /*@@
   @header    Schedule.h
   @date      Mon Sep 13 12:24:48 1999
   @author    Tom Goodale
   @desc 
   Header file for Schedule routines, etc.
   @enddesc 
   @version $Header$
 @@*/

#include "cctki_Schedule.h"
#include "StoreHandledData.h"

#ifndef _SCHEDULE_H_
#define _SCHEDULE_H_

#ifdef __cplusplus
extern "C" {
#endif


int *CCTKi_ScheduleCreateIVec(int size);
void CCTKi_ScheduleDestroyIVec(int size, int *vector);
signed char **CCTKi_ScheduleCreateArray(int size);
void CCTKi_ScheduleDestroyArray(int size, signed char **array);

int CCTKi_ScheduleAddRow(int size, 
                         signed char **array, 
                         int *order, 
                         int item, 
                         int *thisorders);

int CCTKi_ScheduleSort(int size, signed char **array, int *order);

cHandledData *CCTKi_ScheduleGetGroups(void);

cHandledData *CCTKi_DoScheduleGetGroups(void);

#ifdef __cplusplus
}
#endif

/* Internal type data */

typedef enum {sched_item_none, sched_group, sched_function} t_sched_item_type;

typedef enum {sort_order_none, sort_order_ascending, sort_order_descending} t_sort_order;

typedef struct
{
  char *name;
  
  t_sched_item_type type;

  void *function;
  int group;

  int n_whiles;
  char **whiles;
  int n_ifs;
  char **ifs;

  void *attributes;

  t_sched_modifier *modifiers;
} t_sched_item;

typedef struct 
{
  char *name;
  int *order;
  
  int n_scheditems;
  
  t_sched_item *scheditems;

} t_sched_group;


#endif
