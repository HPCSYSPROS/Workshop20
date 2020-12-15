/*@@
   @file      DefaultTimers.c
   @date      Wed Oct 20 16:17:42 1999
   @author    Tom Goodale
   @desc
              Default Cactus timers
   @enddesc
   @version   $Id$
 @@*/

#include <stdio.h>
#include <stdlib.h>

#include "cctk_Config.h"
#include "cctk_Flesh.h"

#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if HAVE_SYS_TIME_H
#  include <sys/time.h>
# elif HAVE_TIME_H
#  include <time.h>
# endif
#endif

#ifdef HAVE_TIME_GETTIMEOFDAY
#include <unistd.h>
#endif

#ifdef HAVE_TIME_GETRUSAGE
#include <sys/resource.h>
#endif

#include "cctk_Timers.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_DefaultTimers_c);


/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/
#ifdef HAVE_TIME_GETTIMEOFDAY
/* A structure to hold the relevent data */
typedef struct
{
  struct timeval total;
  struct timeval last;
} t_GetTimeOfDayTimer;
#endif

#ifdef HAVE_TIME_GETRUSAGE
/* A structure to hold the relevent data */
typedef struct
{
  struct timeval total;
  struct timeval last;
} t_GetrUsageTimer;
#endif


/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/
/* Prototypes for registration functions */
#ifdef HAVE_TIME_GETTIMEOFDAY
static void CCTKi_RegisterTimersGetTimeOfDay(void);
#endif

#ifdef HAVE_TIME_GETRUSAGE
static void CCTKi_RegisterTimersGetrUsage(void);
#endif


/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/
#ifdef HAVE_TIME_GETTIMEOFDAY
static char *GetTimeOfDayHeading = "gettimeofday";
static char *GetTimeOfDayUnits   = "secs";
#endif

#ifdef HAVE_TIME_GETRUSAGE
static char *GetrUsageHeading = "getrusage";
static char *GetrUsageUnits   = "secs";
#endif


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/
int CCTKi_RegisterDefaultTimerFunctions(void);


 /*@@
   @routine    CCTKi_RegisterDefaultTimerFunctions
   @date       Wed Oct 20 18:27:20 1999
   @author     Tom Goodale
   @desc
               Master flesh timer registration function.
   @enddesc
   @calls      CCTKi_RegisterTimersGetTimeOfDay
               CCTKi_RegisterTimersGetrUsage

   @returntype int
   @returndesc
               0
   @endreturndesc
@@*/
int CCTKi_RegisterDefaultTimerFunctions(void)
{
#ifdef HAVE_TIME_GETTIMEOFDAY
  CCTKi_RegisterTimersGetTimeOfDay();
#endif

#ifdef HAVE_TIME_GETRUSAGE
  CCTKi_RegisterTimersGetrUsage();
#endif

  return (0);
}


/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/

/*********************************************************************
 ****************   gettimeofday based timer     *********************
 *********************************************************************/

#ifdef HAVE_TIME_GETTIMEOFDAY

 /*@@
   @routine CCTKi_TimerGetTimeOfDayCreate
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Create the timer structure for use with the gettimeofday function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar

   @returntype void *
   @returndesc
               timer structure, or NULL if allocation failed
   @endreturndesc
@@*/
static void *CCTKi_TimerGetTimeOfDayCreate(int timernum)
{
  t_GetTimeOfDayTimer *retval;


  (void) (timernum + 0);

  retval = calloc (1, sizeof (t_GetTimeOfDayTimer));

  return (retval);
}


 /*@@
   @routine CCTKi_TimerGetTimeOfDayDestroy
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Destroy the timer structure for use with the gettimeofday function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     data
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
@@*/
static void CCTKi_TimerGetTimeOfDayDestroy(int timernum, void *data)
{
  (void) (timernum + 0);
  if(data)
  {
    free(data);
  }
}


 /*@@
   @routine CCTKi_TimerGetTimeOfDayStart
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Start the timer with the gettimeofday function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
@@*/
static void CCTKi_TimerGetTimeOfDayStart(int timernum, void *idata)
{
  struct timeval tp;
  t_GetTimeOfDayTimer *data;


  (void) (timernum + 0);

  gettimeofday(&tp, NULL);

  data = (t_GetTimeOfDayTimer *) idata;
  data->last = tp;
}


 /*@@
   @routine CCTKi_TimerGetTimeOfDayStop
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Stop the timer with the gettimeofday function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
@@*/
static void CCTKi_TimerGetTimeOfDayStop(int timernum, void *idata)
{
  struct timeval tp;
  t_GetTimeOfDayTimer *data;


  (void) (timernum + 0);

  gettimeofday(&tp, NULL);

  data = (t_GetTimeOfDayTimer *) idata;
  data->total.tv_sec  += tp.tv_sec  - data->last.tv_sec;
  data->total.tv_usec += tp.tv_usec - data->last.tv_usec;
}


 /*@@
   @routine CCTKi_TimerGetTimeOfDayReset
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Reset the timer with the gettimeofday function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
@@*/
static void CCTKi_TimerGetTimeOfDayReset(int timernum, void *idata)
{
  t_GetTimeOfDayTimer *data;


  (void) (timernum + 0);

  data = (t_GetTimeOfDayTimer *) idata;
  data->last.tv_sec  = data->last.tv_usec  = 0;
  data->total.tv_sec = data->total.tv_usec = 0;
}


 /*@@
   @routine CCTKi_TimerGetTimeOfDayGet
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Get the time recorded with the gettimeofday function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
   @var     vals
   @vdesc   Timer value structure
   @vtype   cTimerVal
   @vio     out
   @endvar
@@*/
static void CCTKi_TimerGetTimeOfDayGet(int timernum, void *idata, cTimerVal *vals)
{
  t_GetTimeOfDayTimer *data;


  (void) (timernum + 0);

  data = (t_GetTimeOfDayTimer *) idata;
  vals[0].type    = val_double;
  vals[0].heading = GetTimeOfDayHeading;
  vals[0].units   = GetTimeOfDayUnits;
  vals[0].val.d   = data->total.tv_sec + (double)data->total.tv_usec/1000000.0;
  vals[0].seconds = vals[0].val.d;
  vals[0].resolution = 1.e-6;
}


 /*@@
   @routine CCTKi_TimerGetTimeOfDaySet
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Set the time for a gettimeofday function based timer.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
   @var     vals
   @vdesc   Timer value structure
   @vtype   cTimerVal
   @vio     in
   @endvar
@@*/
static void CCTKi_TimerGetTimeOfDaySet(int timernum, void *idata, cTimerVal *vals)
{
  t_GetTimeOfDayTimer *data;


  (void) (timernum + 0);

  data = (t_GetTimeOfDayTimer *) idata;
  /* Note: the struct timeval fields tv_sec and tv_usec are defined in SYSV
   * as time_t and suseconds_t.  But these types don't exist on all systems.
   * They are signed integral types for which long should be enough. */
  data->total.tv_sec  = (long)vals[0].val.d;
  data->total.tv_usec = (long)(1000000*vals[0].val.d)
                                            - (long)data->total.tv_sec;
}


 /*@@
   @routine CCTKi_RegisterTimersGetTimeOfDay
   @date    Wed Oct 20 18:32:17 1999
   @author  Tom Goodale
   @desc
            Register all the timer functions associated with the
            gettimeofday function.
   @enddesc
   @calls   CCTK_ClockRegister
@@*/
static void CCTKi_RegisterTimersGetTimeOfDay(void)
{
  cClockFuncs functions;


  functions.n_vals  = 1;
  functions.create  = CCTKi_TimerGetTimeOfDayCreate;
  functions.destroy = CCTKi_TimerGetTimeOfDayDestroy;
  functions.start   = CCTKi_TimerGetTimeOfDayStart;
  functions.stop    = CCTKi_TimerGetTimeOfDayStop;
  functions.reset   = CCTKi_TimerGetTimeOfDayReset;
  functions.get     = CCTKi_TimerGetTimeOfDayGet;
  functions.set     = CCTKi_TimerGetTimeOfDaySet;

  CCTK_ClockRegister("GetTimeOfDay", &functions);
}

#endif /* HAVE_TIME_GETTIMEOFDAY */


/*********************************************************************
 ****************       getrusage based timer      *******************
 *********************************************************************/

#ifdef HAVE_TIME_GETRUSAGE

 /*@@
   @routine CCTKi_TimerGetrUsageCreate
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Create the timer structure for use with the getrusage function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar

   @returntype void *
   @returndesc
               timer structure, or NULL if allocation failed
   @endreturndesc
@@*/
static void *CCTKi_TimerGetrUsageCreate(int timernum)
{
  t_GetrUsageTimer *retval;


  (void) (timernum + 0);

  retval = calloc (1, sizeof (t_GetrUsageTimer));

  return (retval);
}


 /*@@
   @routine CCTKi_TimerGetrUsageDestroy
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Destroy the timer structure for use with the getrusage function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     data
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
@@*/
static void CCTKi_TimerGetrUsageDestroy(int timernum, void *data)
{
  (void) (timernum + 0);
  if(data)
  {
    free(data);
  }
}


 /*@@
   @routine CCTKi_TimerGetrUsageStart
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Start the timer with the getrusage function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
@@*/
static void CCTKi_TimerGetrUsageStart(int timernum, void *idata)
{
  struct rusage ru;
  t_GetrUsageTimer *data;


  (void) (timernum + 0);

  getrusage(RUSAGE_SELF, &ru);

  data = (t_GetrUsageTimer *) idata;
  data->last = ru.ru_utime;
}


 /*@@
   @routine CCTKi_TimerGetrUsageStop
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Stop the timer with the getrusage function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
@@*/
static void CCTKi_TimerGetrUsageStop(int timernum, void *idata)
{
  struct rusage ru;
  t_GetrUsageTimer *data;


  (void) (timernum + 0);

  getrusage(RUSAGE_SELF, &ru);

  data = (t_GetrUsageTimer *) idata;
  data->total.tv_sec  += ru.ru_utime.tv_sec  - data->last.tv_sec;
  data->total.tv_usec += ru.ru_utime.tv_usec - data->last.tv_usec;
}


 /*@@
   @routine CCTKi_TimerGetrUsageReset
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Reset the timer with the getrusage function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
@@*/
static void CCTKi_TimerGetrUsageReset(int timernum, void *idata)
{
  t_GetrUsageTimer *data;


  (void) (timernum + 0);

  data = (t_GetrUsageTimer *) idata;
  data->total.tv_sec = data->total.tv_usec = 0;
}


 /*@@
   @routine CCTKi_TimerGetrUsageGet
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Get the time recorded with the getrusage function.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
   @var     vals
   @vdesc   Timer value structure
   @vtype   cTimerVal
   @vio     out
   @endvar
@@*/
static void CCTKi_TimerGetrUsageGet(int timernum, void *idata, cTimerVal *vals)
{
  t_GetrUsageTimer *data;


  (void) (timernum + 0);

  data = (t_GetrUsageTimer *) idata;
  vals[0].type    = val_double;
  vals[0].heading = GetrUsageHeading;
  vals[0].units   = GetrUsageUnits;
  vals[0].val.d   = data->total.tv_sec + (double)data->total.tv_usec/1000000.0;
  vals[0].seconds = vals[0].val.d;
  vals[0].resolution = 1.e-6;
}


 /*@@
   @routine CCTKi_TimerGetrUsageSet
   @date    Wed Oct 20 18:28:19 1999
   @author  Tom Goodale
   @desc
            Set the time for a getrusage function based timer.
   @enddesc

   @var     timernum
   @vdesc   timer number
   @vtype   int
   @vio     unused
   @endvar
   @var     idata
   @vdesc   timer data
   @vtype   void *
   @vio     inout
   @endvar
   @var     vals
   @vdesc   Timer value structure
   @vtype   cTimerVal *
   @vio     in
   @endvar
@@*/
static void CCTKi_TimerGetrUsageSet(int timernum, void *idata, cTimerVal *vals)
{
  t_GetrUsageTimer *data;


  (void) (timernum + 0);

  data = (t_GetrUsageTimer *) idata;
  data->total.tv_sec  = (long)vals[0].val.d;
  data->total.tv_usec = (long)(1000000*vals[0].val.d)
                                            - (long)data->total.tv_sec;
}


/*@@
   @routine CCTKi_RegisterTimersGetrUsage
   @date    Wed Oct 20 18:32:17 1999
   @author  Tom Goodale
   @desc
            Register all the timer functions associated with the
            getrusage function.
   @enddesc
   @calls   CCTK_ClockRegister
@@*/
static void CCTKi_RegisterTimersGetrUsage(void)
{
  cClockFuncs functions;

  functions.n_vals  = 1;
  functions.create  = CCTKi_TimerGetrUsageCreate;
  functions.destroy = CCTKi_TimerGetrUsageDestroy;
  functions.start   = CCTKi_TimerGetrUsageStart;
  functions.stop    = CCTKi_TimerGetrUsageStop;
  functions.reset   = CCTKi_TimerGetrUsageReset;
  functions.get     = CCTKi_TimerGetrUsageGet;
  functions.set     = CCTKi_TimerGetrUsageSet;

  CCTK_ClockRegister("GetrUsage", &functions);
}

#endif /* HAVE_TIME_GETRUSAGE */
