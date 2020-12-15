 /*@@
   @file      Time.c
   @date      Wed Sep 17 2000
   @author    Gabrielle Allen
   @desc 
   Date and time routines
   @enddesc 
   @version $Header$
 @@*/

/* #define DEBUG_TIME */

#include <ctype.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "cctk_Flesh.h"
#include "cctk_Misc.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(util_Time_c);


/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    Util_CurrentTime
   @date       Tue Sep 19
   @author     Gabrielle Allen
   @desc 
   Fills string with current local time, returning the number
   of characters filled.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int Util_CurrentTime(int len, char *now)
{
  int retval;
  time_t timep;
  const char *fmt = "%X%z";

  timep = time(NULL);
  strftime(now, len, fmt, localtime(&timep));

  retval = strlen(now);
  retval=retval > len ? 0 : retval;

#ifdef DEBUG_TIME
  printf("CurrentTime = %s\n",now);
#endif

  return retval;
}


 /*@@
   @routine    Util_CurrentDate
   @date       Tue Sep 19
   @author     Gabrielle Allen
   @desc 
   Fills string with current local date, returning the number of
   characters filled.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

int Util_CurrentDate(int len, char *now)
{
  int retval;
  time_t timep;
  const char *fmt = "%b %d %Y";

  timep = time(NULL);
  strftime(now, 50, fmt, localtime(&timep));

  retval = strlen(now);
  retval=retval > len ? 0 : retval;

#ifdef DEBUG_TIME
  printf("CurrentDate = %s\n",thedate);
#endif

  return retval;
}



 /*@@
   @routine    Util_CurrentDateTime
   @date       Wed 19 July 2006
   @author     Thomas Radke
   @desc
               Returns the current datetime in a machine-processable format
               as defined in ISO 8601 chapter 5.4.
   @enddesc

   @returntype char *
   @returndesc
               pointer to an allocated string buffer containing the datetime,
               must be freed by the user
   @endreturndesc
@@*/
char *Util_CurrentDateTime(void)
{
  char *buffer;
  time_t timep;
  const int len = sizeof ("YYYY-MM-DDThh:mm:ss+hh:mm");
  const char *fmt = "%Y-%m-%dT%H:%M:%S%z";

  buffer = calloc (1, len);
  if (buffer)
  {
    timep = time (NULL);
    strftime (buffer, len, fmt, localtime (&timep));

    /* if the timezone part is returned as "(+|-)hhmm"
       then turn it into "(+|-)hh:mm" */
    if ((buffer[len-7] == '+' || buffer[len-7] == '-') &&
        isdigit (buffer[len-6]) && isdigit (buffer[len-5]) &&
        isdigit (buffer[len-4]) && isdigit (buffer[len-3]) &&
        buffer[len-2] == 0)
    {
      buffer[len-2] = buffer[len-3];
      buffer[len-3] = buffer[len-4];
      buffer[len-4] = ':';
    }
  }

  return (buffer);
}
