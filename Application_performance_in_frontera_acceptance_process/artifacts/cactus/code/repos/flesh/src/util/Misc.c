 /*@@
   @file      Misc.c
   @date      Wed Jan 20 10:06:35 1999
   @author    Tom Goodale
   @desc
              Miscellaneous routines.
   @enddesc
   @version   $Id$
 @@*/

/*#define DEBUG_MISC*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include "cctk_Flesh.h"
#include "cctk_Misc.h"
#include "cctk_FortranString.h"
#include "cctk_WarnLevel.h"
#include "util_String.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(util_Misc_c);


/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

int CCTK_SetStringInRegexList(char **data, const char *value,
                              int n_elements, ...);
void CCTK_PrintString(const char *data);

int CCTK_FCALL CCTK_FNAME(CCTK_Equals)
                         (const char **arg1, ONE_FORTSTRING_ARG);

CCTK_POINTER CCTK_FCALL CCTK_FNAME(CCTK_PointerTo)
                                  (CCTK_POINTER var);
CCTK_POINTER CCTK_FCALL CCTK_FNAME(CCTK_NullPointer)
                                  (void);

void CCTK_FCALL CCTK_FNAME(CCTK_PrintString)
                          (const char **arg1);

void CCTK_FCALL CCTK_FNAME (CCTK_FortranString)
                           (CCTK_INT *nchars,
                            const char *const *cstring,
                            ONE_FORTSTRING_ARG);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

 /*@@
   @routine    CCTK_Equals
   @date       Wed Jan 20 10:25:30 1999
   @author     Tom Goodale
   @desc
   Does a case independent comparison of strings.
   Returns true if they are equal.
   @enddesc

   @var     string1
   @vdesc   First string in comparison
   @vtype   const char *
   @vio     in
   @endvar
   @var     string2
   @vdesc   Second string in comparison
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - equal
   0 - not equal
   @endreturndesc
@@*/
int CCTK_Equals(const char *string1, const char *string2)
{
  int retval;

  retval = 1;

  /* Check that string1 isn't null */
  if (!string1 || !string2)
  {
    retval = 0;
    if (!string1 && string2)
    {
      CCTK_VError(__LINE__,__FILE__,"Cactus",
                  "CCTK_Equals: First string is null (2nd is \"%s\")",string2);

    }
    else if (string1 && !string2)
    {
      CCTK_VError(__LINE__,__FILE__,"Cactus",
                  "CCTK_Equals: Second string is null (1st is \"%s\")",string1);
    }
    else
    {
      CCTK_Error(__LINE__,__FILE__,"Cactus",
                 "CCTK_Equals: Both strings are null");
    }
  }
  else
  {
    if (Util_StrCmpi(string1,string2))
    {
      retval = 0;
    }
  }
  return retval;
}

int CCTK_FCALL CCTK_FNAME(CCTK_Equals)
     (const char **arg1,ONE_FORTSTRING_ARG)
{
  int retval;
  const char *string1 = *arg1;
  ONE_FORTSTRING_PTR(string2);
  ONE_FORTSTRING_LEN(len2);

  retval = 1;

  /* Check that string1 isn't null */
  if (!string1)
  {
    retval = 0;
    CCTK_VWarn(0,__LINE__,__FILE__,"Cactus",
               "CCTK_Equals: First string null (2nd is %s)",
               Util_NullTerminateString(string2,len2));
  }
  else
  {
    if (Util_StrMemCmpi(string1,string2,len2))
    {
      retval = 0;
    }
  }
  return retval;
}


 /*@@
   @routine    CCTK_PointerTo
   @date       Tue 10 Dec 2002
   @author     Thomas Radke
   @desc
               Returns the pointer to a variable which is passed in
               (by reference) from a fortran routine.
   @enddesc

   @var        var
   @vdesc      variable passed in from fortran by reference
   @vtype      CCTK_POINTER
   @vio        in

   @returntype CCTK_POINTER
   @returndesc
               pointer to the variable
   @endreturndesc
@@*/
CCTK_POINTER CCTK_FCALL CCTK_FNAME (CCTK_PointerTo) (CCTK_POINTER var)
{
  return (var);
}

 /*@@
   @routine    CCTK_NullPointer
   @date       Wed 12 Feb 2003
   @author     Thomas Radke
   @desc
               Returns a NULL pointer to be used for the initialization of
               CCTK_POINTER variables in the calling fortran routine.
   @enddesc

   @returntype CCTK_POINTER
   @returndesc
               a C NULL pointer
   @endreturndesc
@@*/
CCTK_POINTER CCTK_FCALL CCTK_FNAME (CCTK_NullPointer) (void)
{
  return (NULL);
}

 /*@@
   @routine Util_NullTerminateString
   @author Paul Walker
   @desc
   Null terminates a fortran string. Need to free the
   string it returns
   @enddesc

   @var     instring
   @vdesc   String to null terminate
   @vtype   const char *
   @vio     in
   @endvar
   @var     len
   @vdesc   Length of string to be null terminated
   @vtype   unsigned int
   @vio     in
   @endvar

   @returntype char *
   @returndesc
   The null terminated string.
   @endreturndesc
@@*/

char *Util_NullTerminateString(const char *instring, unsigned int len)
{
  char *outstring;
  unsigned int i;
  unsigned int position;

  if (len > 100000)
  {
    CCTK_VWarn(CCTK_WARN_ALERT,__LINE__,__FILE__,"Cactus",
               "Null-terminating a string with length %d; "
               "this is probably an error in calling a C routine from Fortran",
               len);
  }

#ifdef DEBUG_MISC
  printf("Util_NullTerminateString: -%s-, (%u)\n",instring,len);
#endif

  position = len;
  while (position > 0 && instring[position-1] == ' ')
  {
    position--;
  }

  outstring = (char *)malloc((position+1)*sizeof(char));

  if (outstring)
  {
    for (i=0;i<position;i++)
    {
      outstring[i] = instring[i];
    }
    outstring[position] = '\0';
  }
  return(outstring);
}


 /*@@
   @routine    Util_InList
   @date       Wed Jan 20 10:31:25 1999
   @author     Tom Goodale
   @desc
   Determines if a string is in a list of other strings.
   @enddesc

   @var     string1
   @vdesc   The string to search for
   @vtype   const char *
   @vio     in
   @endvar
   @var     n_elements
   @vdesc   The number of elements in the list
   @vtype   int
   @vio     in
   @endvar
   @var     ...
   @vdesc   List of strings to search.
   @vtype   multiple const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - In list
   0 - Not in list
   @endreturndesc
@@*/
int Util_InList(const char *string1, int n_elements, ...)
{
  int retval;
  int arg;
  va_list ap;

  char *element;

  retval = 0;

  /* Walk through the element list. */
  va_start(ap, n_elements);

  for(arg = 0; arg < n_elements; arg++)
  {
    element = va_arg(ap, char *);

    if(CCTK_Equals(string1, element))
    {
      retval = 1;
      break;
    }
  }

  va_end(ap);

  return retval;

}


 /*@@
   @routine    Util_IntInRange
   @date       Wed Jan 20 10:32:36 1999
   @author     Tom Goodale
   @desc
   This routine will determine if an integer is in the range specified
   in the range string.
   @enddesc

   @var     inval
   @vdesc   The value to check
   @vtype   int
   @vio     in
   @endvar
   @var     range
   @vdesc   The range to look in
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - in range
   0 - not in range
   @endreturndesc
@@*/
int Util_IntInRange(int inval, const char *range)
{
  int retval = 0;
  regmatch_t pmatch[8];
  int start_closed, end_closed;
  int start, end, step;
  int matched;

  /* Valid ranges are of the form start:end:step
   * possibly preceeded by a [ or ( and ended by a ) or ] to indicate
   * closure.  The end and step are optional.  A * can also be used
   * to indicate any value.
   *
   * The following regular expression may match five substrings:
   *
   * 1 - [ or (
   * 2 - start
   * 3 - whether end is present
   * 4 - end
   * 5 - whether step is present
   * 6 - step
   * 7 - ) or ]
   */

#define R_BEGIN "(\\[|\\()?"
#define R_VALUE "([^]):]*)"
#define R_SEP   ":"
#define R_END   "(\\]|\\))?"

#define R_MAYBE(x) "(" x ")?"

  const char pattern[] =
    R_BEGIN
    R_VALUE
    R_MAYBE(R_SEP
            R_VALUE
            R_MAYBE(R_SEP
                    R_VALUE))
    R_END;

#undef R_BEGIN
#undef R_VALUE
#undef R_SEP
#undef R_END

#undef R_MAYBE

  matched = CCTK_RegexMatch(range, pattern, 8, pmatch);
  if(matched > 0)
  {
    /* First work out if the range is closed at the lower end. */
    if(pmatch[1].rm_so != -1)
    {
      switch(range[pmatch[1].rm_so])
      {
        case '(' : start_closed = 0; break;
        case '[' :
        default  : start_closed = 1;
      }
    }
    else
    {
      start_closed = 1;
    }

    /* Next find the start of the range */
    if(pmatch[2].rm_so != -1 &&
       (pmatch[2].rm_eo-pmatch[2].rm_so > 0) &&
       range[pmatch[2].rm_so] != '*')
    {
      char *endptr;
      start = strtol(range+pmatch[2].rm_so, &endptr, 10);
      if (endptr == range+pmatch[2].rm_so)
      {
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                   "Invalid start of integer parameter range; range descriptor is \"%s\" (value is %d)",
                   range, inval);
        start = INT_MIN;
      }
    }
    else
    {
      /* No start range given, so use the smallest integer available. */
      start = INT_MIN;
    }

    /* Next find whether end of the range is present */
    if(pmatch[3].rm_so != -1)
    {
      /* Next find the end of the range */
      if(pmatch[4].rm_so != -1 &&
         (pmatch[4].rm_eo-pmatch[4].rm_so > 0) &&
         range[pmatch[4].rm_so] != '*')
      {
        char *endptr;
        end = strtol(range+pmatch[4].rm_so, &endptr, 10);
        if (endptr == range+pmatch[4].rm_so)
        {
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                     "Invalid end of integer parameter range; range descriptor is \"%s\" (value is %d)",
                     range, inval);
          end = INT_MAX;
        }
      }
      else
      {
        /* No end range given, so use the largest integer available. */
        end = INT_MAX;
      }
    }
    else
    {
      /* End range not given, so interval has either zero or infinite length. */
      end = start == INT_MIN ? INT_MAX : start;
    }

    /* Next find whether step of the range is present */
    if(pmatch[5].rm_so != -1)
    {
      /* Next find the step of the range */
      if(pmatch[6].rm_so != -1 && (pmatch[6].rm_eo-pmatch[6].rm_so > 0))
      {
        char *endptr;
        step = strtol(range+pmatch[6].rm_so, &endptr, 10);
        if (endptr == range+pmatch[6].rm_so)
        {
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                     "Invalid step of integer parameter range; range descriptor is \"%s\" (value is %d)",
                    range, inval);
          step = 1;
        }
      }
      else
      {
        /* No step given, so default to 1. */
        step = 1;
      }
    }
    else
    {
      /* Step not given, so default to 1. */
      step = 1;
    }
    if (step <= 0)
    {
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                 "Non-positive step of integer parameter range; range descriptor is \"%s\" (value is %d)",
                 range, inval);
      step = 1;
    }

    /* Finally work out if the range is closed at the upper end. */
    if(pmatch[7].rm_so != -1)
    {
      switch(range[pmatch[7].rm_so])
      {
        case ')' : end_closed = 0; break;
        case ']' :
        default  : end_closed = 1;
      }
    }
    else
    {
      end_closed = 1;
    }

    /* Cast to unsigned long, because the subtraction (inval - start)
       can overflow, and wrap-around is legal in C only for unsigned
       types.  */
    if(inval >= start + !start_closed &&
       inval <= end   - !end_closed   &&
       ! (((unsigned long)inval-(unsigned long)start) % (unsigned long)step))
    {
      retval = 1;
    }

  }
  else if(!matched)
  {
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
               "Invalid range of integer parameter range; range descriptor is \"%s\" (value is %d)",
               range, inval);
  }
  else
  {
    CCTK_VError(__LINE__, __FILE__, "Cactus",
                "Invalid patten '%s' used to parse range", pattern);
  }

  return retval;
}

 /*@@
   @routine    Util_DoubleInRange
   @date       Wed Jan 20 10:32:36 1999
   @author     Tom Goodale
   @desc
   This routine will determine if a double is in the range specified
   in the range string.

   Valid ranges are of the form start:end:step
   possibly preceeded by a [ or ( and ended by a ) or ] to indicate
   closure.  The end and step are optional.  A * can also be used
   to indicate any value.
   @enddesc

   @var     inval
   @vdesc   The value to check
   @vtype   double
   @vio     in
   @endvar
   @var     range
   @vdesc   The range to look in
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - in range
   0 - not in range
   @endreturndesc
@@*/
int Util_DoubleInRange(double inval, const char *range)
{
  int retval = 0;
  regmatch_t pmatch[8];
  int start_closed, end_closed;
  double start, end, step;

  /*
   * The following regular expression may match five substrings:
   *
   * 1 - [ or (
   * 2 - start
   * 3 - whether end is present
   * 4 - end
   * 5 - whether step is present
   * 6 - step
   * 7 - ) or ]
   */

#define R_BEGIN "(\\[|\\()?"
#define R_VALUE "([^]):]*)"
#define R_SEP   ":"
#define R_END   "(\\]|\\))?"

#define R_MAYBE(x) "(" x ")?"

  const char pattern[] =
    R_BEGIN
    R_VALUE
    R_MAYBE(R_SEP
            R_VALUE
            R_MAYBE(R_SEP
                    R_VALUE))
    R_END;

#undef R_BEGIN
#undef R_VALUE
#undef R_SEP
#undef R_END

#undef R_MAYBE

  if(CCTK_RegexMatch(range, pattern, 8, pmatch) > 0)
  {
    /* First work out if the range is closed at the lower end. */
    if(pmatch[1].rm_so != -1)
    {
      switch(range[pmatch[1].rm_so])
      {
        case '(' : start_closed = 0; break;
        case '[' :
        default  : start_closed = 1;
      }
    }
    else
    {
      start_closed = 1;
    }

    /* Next find the start of the range */
    if(pmatch[2].rm_so != -1 &&
       (pmatch[2].rm_eo-pmatch[2].rm_so > 0) &&
       range[pmatch[2].rm_so] != '*')
    {
      char *endptr;
      start = strtod(range+pmatch[2].rm_so, &endptr);
      if (endptr == range+pmatch[2].rm_so)
      {
        CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                   "Invalid start of real-valued parameter range; range descriptor is \"%s\" (value is %.17g)",
                   range, inval);
        start = -HUGE_VAL;
      }
    }
    else
    {
      /* No start range given, so use the smallest double available. */
      start = -HUGE_VAL;
    }

    /* Next find whether end of the range is present */
    if(pmatch[3].rm_so != -1)
    {
      /* Next find the end of the range */
      if(pmatch[4].rm_so != -1 &&
         (pmatch[4].rm_eo-pmatch[4].rm_so > 0) &&
         range[pmatch[4].rm_so] != '*')
      {
        char *endptr;
        end = strtod(range+pmatch[4].rm_so, &endptr);
        if (endptr == range+pmatch[4].rm_so)
        {
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                     "Invalid end of real-valued parameter range; range descriptor is \"%s\" (value is %.17g)",
                     range, inval);
          end = HUGE_VAL;
        }
      }
      else
      {
        /* No end range given, so use the largest double available. */
        end = HUGE_VAL;
      }
    }
    else
    {
      /* End range not given, so interval has either zero or infinite length. */
      end = start == -HUGE_VAL ? HUGE_VAL : start;
    }

    /* Next find whether step of the range is present */
    if(pmatch[5].rm_so != -1)
    {
      /* Next find the step of the range */
      if(pmatch[6].rm_so != -1 && (pmatch[6].rm_eo-pmatch[6].rm_so > 0))
      {
        char *endptr;
        step = strtod(range+pmatch[6].rm_so, &endptr);
        if (endptr == range+pmatch[6].rm_so)
        {
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                     "Invalid step of real-valued parameter range; range descriptor is \"%s\" (value is %.17g)",
                     range, inval);
          step = 0;
        }
      }
      else
      {
        /* No step given, so default to 0. */
        step = 0;
      }
    }
    else
    {
      /* Step not given, so default to 0. */
      step = 0;
    }
    if (step < 0)
    {
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
                 "Negative step of real-valued parameter range; range descriptor is \"%s\" (value is %.17g)",
                 range, inval);
      step = 0;
    }

    /* Finally work out if the range is closed at the upper end. */
    if(pmatch[7].rm_so != -1)
    {
      switch(range[pmatch[7].rm_so])
      {
        case ')' : end_closed = 0; break;
        case ']' :
        default  : end_closed = 1;
      }
    }
    else
    {
      end_closed = 1;
    }

    if((start_closed ? inval >= start : inval > start) &&
       (end_closed   ? inval <= end   : inval < end  ) &&
       (step         ? ! fmod (inval - start, step) : 1))
    {
      retval = 1;
    }

  }
  else
  {
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, "Cactus",
               "Invalid range of real-valued parameter range; range descriptor is \"%s\" (value is %.17g)",
               range, inval);
  }

  return retval;
}


 /*@@
   @routine    Util_IntInRangeList
   @date       Wed Jan 20 10:36:31 1999
   @author     Tom Goodale
   @desc
   Determines if an integer is in a given list of ranges.
   @enddesc

   @var     inval
   @vdesc   The value to check
   @vtype   int
   @vio     in
   @endvar
   @var     n_elements
   @vdesc   The number of elements in the list
   @vtype   int
   @vio     in
   @endvar

   @var     ...
   @vdesc   The list of ranges to look in
   @vtype   multiple const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - in range in list
   0 - not in range in list
   @endreturndesc
@@*/
int Util_IntInRangeList(int inval, int n_elements, ...)
{
  int retval;
  int arg;
  va_list ap;

  char *element;

  retval = 0;

  /* Walk through the element list. */
  va_start(ap, n_elements);

  for(arg = 0; arg < n_elements; arg++)
  {
    element = va_arg(ap, char *);

    if(Util_IntInRange(inval, element))
    {
      retval = 1;
      break;
    }
  }

  va_end(ap);

  return retval;

}


 /*@@
   @routine    Util_DoubleInRangeList
   @date       Wed Jan 20 10:36:31 1999
   @author     Tom Goodale
   @desc
   Determines if a double is in a given list of ranges.
   @enddesc

   @var     inval
   @vdesc   The value to check
   @vtype   double
   @vio     in
   @endvar
   @var     n_elements
   @vdesc   The number of elements in the list
   @vtype   int
   @vio     in
   @endvar

   @var     ...
   @vdesc   The list of ranges to look in
   @vtype   multiple const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - in range in list
   0 - not in range in list
   @endreturndesc
@@*/
int Util_DoubleInRangeList(double inval, int n_elements, ...)
{
  int retval;
  int arg;
  va_list ap;

  char *element;

  retval = 0;

  /* Walk through the element list. */
  va_start(ap, n_elements);

  for(arg = 0; arg < n_elements; arg++)
  {
    element = va_arg(ap, char *);

    if(Util_DoubleInRange(inval, element))
    {
      retval = 1;
      break;
    }
  }

  va_end(ap);

  return retval;

}


 /*@@
   @routine    CCTK_SetDoubleInRangeList
   @date       Thu Jan 21 09:41:21 1999
   @author     Tom Goodale
   @desc
   Sets the value of a double if the desired value is in one of
   the specified ranges.
   @enddesc

   @var     data
   @vdesc   Pointer to the value to set
   @vtype   CCTK_REAL *
   @vio     out
   @endvar
   @var     value
   @vdesc   The value to check
   @vtype   const char *
   @vio     in
   @endvar
   @var     n_elements
   @vdesc   The number of elements in the list
   @vtype   int
   @vio     in
   @endvar

   @var     ...
   @vdesc   The list of ranges to look in
   @vtype   multiple const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - in range in list
   0 - not in range in list
   @endreturndesc
@@*/
int CCTK_SetDoubleInRangeList(CCTK_REAL *data, const char *value,
                              int n_elements, ...)
{
  int retval;
  char temp[1001];
  unsigned int p;
  int arg;
  va_list ap;

  char *element;
  char *endptr;

  CCTK_REAL inval;

  retval = 1;

  /* Convert the value string to a double.
   * Allow various formats.
   */
  strncpy(temp, value, sizeof(temp) - 1);
  temp[sizeof(temp) - 1] = 0;

  for (p=0;temp[p];p++)
  {
    if (temp[p] == 'E' ||
        temp[p] == 'd' ||
        temp[p] == 'D')
    {
      temp[p] = 'e';
      break;
    }
  }

  inval = strtod(temp, &endptr);
  if (endptr == temp)
  {
    retval = 0;
    return retval;
  }

  /* Walk through the element list. */
  va_start(ap, n_elements);

  for(arg = 0; arg < n_elements; arg++)
  {
    element = va_arg(ap, char *);

    if(Util_DoubleInRange(inval, element))
    {
      retval = 0;
      *data = inval;
      break;
    }
  }

  va_end(ap);

  return retval;
}

 /*@@
   @routine    CCTK_SetIntInRangeList
   @date       Thu Jan 21 10:27:26 1999
   @author     Tom Goodale
   @desc
   Sets the value of an integer if the desired value is in one of
   the specified ranges.
   @enddesc

   @var     data
   @vdesc   Pointer to the value to set
   @vtype   CCTK_INT *
   @vio     out
   @endvar
   @var     value
   @vdesc   The value to check
   @vtype   const char *
   @vio     in
   @endvar
   @var     n_elements
   @vdesc   The number of elements in the list
   @vtype   int
   @vio     in
   @endvar

   @var     ...
   @vdesc   The list of ranges to look in
   @vtype   multiple const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - in range in list
   0 - not in range in list
   @endreturndesc
@@*/
int CCTK_SetIntInRangeList(CCTK_INT *data, const char *value,
                           int n_elements, ...)
{
  int retval;
  int arg;
  va_list ap;

  char *element;
  char *endptr;

  CCTK_INT inval;

  retval = 1;

  /* Convert the value string to an int.*/

  inval = strtol(value, &endptr, 10);
  if (endptr == value)
  {
    retval = 0;
    return retval;
  }

  /* Walk through the element list. */
  va_start(ap, n_elements);

  for(arg = 0; arg < n_elements; arg++)
  {
    element = va_arg(ap, char *);

    if(Util_IntInRange(inval, element))
    {
      retval = 0;
      *data = inval;
      break;
    }
  }

  va_end(ap);

  return retval;
}

 /*@@
   @routine    CCTK_SetKeywordInRangeList
   @date       Thu Jan 21 10:28:00 1999
   @author     Tom Goodale
   @desc
   Sets the value of a keyword if the desired value is in one of
   the specified ranges.
   @enddesc

   @var     data
   @vdesc   Pointer to the value to set
   @vtype   char **
   @vio     out
   @endvar
   @var     value
   @vdesc   The value to check
   @vtype   const char *
   @vio     in
   @endvar
   @var     n_elements
   @vdesc   The number of elements in the list
   @vtype   int
   @vio     in
   @endvar

   @var     ...
   @vdesc   The list of ranges to look in
   @vtype   multiple const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - in range in list
   0 - not in range in list
   @endreturndesc
@@*/
int CCTK_SetKeywordInRangeList(char **data, const char *value,
                               int n_elements, ...)
{
  int retval;
  int arg;
  va_list ap;

  char *element;

  retval = 1;

  /* Walk through the element list. */
  va_start(ap, n_elements);

  for(arg = 0; arg < n_elements; arg++)
  {
    element = va_arg(ap, char *);

    if(CCTK_Equals(value, element))
    {
      if(*data) free(*data);
      *data = (char *)malloc((strlen(value)+1)*sizeof(char));
      if(*data)
      {
        strcpy(*data, value);
        retval = 0;
      }
      else
      {
        retval =-1;
      }
      break;
    }
  }

  va_end(ap);

  return retval;
}


 /*@@
   @routine    CCTK_SetStringInRegexList
   @date       Fri Apr 16 08:37:02 1999
   @author     Tom Goodale
   @desc
   Sets the value of a string if it matches any of the given regular
   expressions.
   @enddesc

   @var     data
   @vdesc   Pointer to the value to set
   @vtype   char **
   @vio     out
   @endvar
   @var     value
   @vdesc   The value to check
   @vtype   const char *
   @vio     in
   @endvar
   @var     n_elements
   @vdesc   The number of elements in the list
   @vtype   int
   @vio     in
   @endvar

   @var     ...
   @vdesc   The list of ranges to look in
   @vtype   multiple const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   1 - in range in list
   0 - not in range in list
   @endreturndesc
@@*/
int CCTK_SetStringInRegexList(char **data, const char *value,
                               int n_elements, ...)
{
  int retval;
  int arg;
  va_list ap;

  char *element;

  retval = 1;

  /* Walk through the element list. */
  va_start(ap, n_elements);

  for(arg = 0; arg < n_elements; arg++)
  {
    element = va_arg(ap, char *);

    if(CCTK_RegexMatch(value, element, 0, NULL) > 0)
    {
      retval = CCTK_SetString(data, value);
      break;
    }
  }

  va_end(ap);

  return retval;
}

 /*@@
   @routine    CCTK_SetString
   @date       Thu Jan 21 10:28:27 1999
   @author     Tom Goodale
   @desc
   Sets the value of a string
   @enddesc

   @var     data
   @vdesc   Pointer to the value to set
   @vtype   char **
   @vio     out
   @endvar
   @var     value
   @vdesc   The value to check
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0  - success
   -1 - out of memory
   @endreturndesc
@@*/
int CCTK_SetString(char **data, const char *value)
{
  int retval;

  retval = 1;

  if(*data) free(*data);
  *data = (char *)malloc((strlen(value)+1)*sizeof(char));
  if(*data)
  {
    strcpy(*data, value);
    retval = 0;
  }
  else
  {
    retval = -1;
  }

  return retval;
}


 /*@@
   @routine    CCTK_SetBoolean
   @date       Thu Jan 21 10:35:11 1999
   @author     Tom Goodale
   @desc
   Sets the value of a boolean to true or false according to
   the value of the value string.
   @enddesc

   @var     data
   @vdesc   Pointer to the value to set
   @vtype   CCTK_INT *
   @vio     out
   @endvar
   @var     value
   @vdesc   The value to check
   @vtype   const char *
   @vio     in
   @endvar

   @returntype int
   @returndesc
   0  - success
   -1 - out of memory
   @endreturndesc
@@*/
int CCTK_SetBoolean(CCTK_INT *data, const char *value)
{
  int retval = 0;

  if(Util_InList(value, 5, "true", "t", "yes", "y", "1"))
  {
    *data = 1;
  }
  else if(Util_InList(value, 5, "false", "f", "no", "n", "0"))
  {
    *data = 0;
  }
  else
  {
    retval = -1;
  }

  return retval;
}

 /*@@
   @routine    CCTK_RegexMatch
   @date       Fri Apr 16 08:40:14 1999
   @author     Tom Goodale
   @desc
   Perform a regular expression match of string against pattern.
   Also returns the specified number of matched substrings as
   give by regexec.
   This is a modified form of the example routine given in the SGI
   man page for regcomp.
   @enddesc

   @var     string
   @vdesc   String to match against
   @vtype   const char *
   @vio     in
   @endvar
   @var     pattern
   @vdesc   Regex pattern
   @vtype   const char *
   @vio     in
   @endvar
   @var     nmatch
   @vdesc   The size of the pmatch array
   @vtype   const int
   @vio     in
   @endvar
   @var     pmatch
   @vdesc   Array in which to place the matches
   @vtype   regmatch_t
   @vio     out
   @endvar

   @returntype int
   @returndesc
   1 - pattern matches
   0 - pattern doesn't match
   <0 - pattern is invalid
   @endreturndesc
@@*/

int CCTK_RegexMatch(const char *string,
                    const char *pattern,
                    const int nmatch,
                    regmatch_t *pmatch)
{
  int retval;
  int status;
  regex_t re;

  /* BSD says: an empty string is not a legal regular expression.
     Handle this case specially.  */
  if (strcmp(pattern, "") == 0)
  {
    retval = 1;                 /* report success */
  }
  else
  {
    if (regcomp(&re, pattern, REG_EXTENDED) == 0)
    {
      status = regexec(&re, string, (size_t)nmatch, pmatch, 0);
      regfree(&re);
      if (status != 0)
      {
        retval = 0;      /* report error */
      }
      else
      {
        retval = 1;
      }
    }
    else
    {
      retval = -1;       /* repost error */
    }
  }

  return retval;
}

 /*@@
   @routine    CCTK_PrintString
   @date       Fri Apr 1 1999
   @author     Gabrielle Allen
   @desc
   Prints the value of a string (this is for fortran)
   @enddesc

   @var     data
   @vdesc   string to print
   @vtype   const char *
   @vio     in
   @endvar

@@*/
void CCTK_PrintString(const char *data)
{
  printf("%s\n",data);
}

void CCTK_FCALL CCTK_FNAME(CCTK_PrintString)
     (const char **arg1)
{
  CCTK_PrintString(*arg1);
}


/*@@
   @routine    CCTK_FortranString
   @date       Thu Jan 22 14:44:39 1998
   @author     Paul Walker
   @desc
               Copies a C string into a Fortran string.
   @enddesc

   @var        nchars
   @vdesc      Number of characters in the C string
   @vtype      CCTK_INT *
   @vio        out
   @vcomment
               It will copy only as many characters as fit into the fortran
               string. You should check for truncation by comparing 'nchars'
               against the length of your fortran string.
   @endvar
   @var        c_string
   @vdesc      C string to be copied
   @vtype      const char *const *
   @vio        in
   @endvar
   @var        ONE_FORTSTRING_ARG
   @vdesc      Fortran string
   @vtype      FORTRAN string macro
   @vio        out
   @endvar
@@*/
int CCTK_FortranString (const char *c_string,
                        char *fortran_string,
                        int fortran_length)
{
  int nchars, c_strlen;


  nchars = c_strlen = strlen (c_string);
  if (c_strlen > fortran_length)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_FortranString: Fortran string buffer is too short to "
                "hold C string '%s, string will be truncated", c_string);
    c_strlen = fortran_length;
  }

  /* copy up to the size of the fortran string
     and pad remaining chars in the fortran string with spaces */
  memcpy (fortran_string, c_string, c_strlen);
  memset (fortran_string + c_strlen, ' ', fortran_length - c_strlen);

  return nchars;
}

void CCTK_FCALL CCTK_FNAME (CCTK_FortranString)
                           (CCTK_INT *nchars,
                            const char *const *c_string,
                            ONE_FORTSTRING_ARG)
{
  size_t c_strlen;
  ONE_FORTSTRING_PTR (fortran_string)
  ONE_FORTSTRING_LEN (fortran_length)


  *nchars = c_strlen = strlen (*c_string);
  if (c_strlen > fortran_length)
  {
    CCTK_VWarn (1, __LINE__, __FILE__, "Cactus",
                "CCTK_FortranString: Fortran string buffer is too short to "
                "hold C string '%s', string will be truncated", *c_string);
    c_strlen = fortran_length;
  }

  /* copy up to the size of the fortran string
     and pad remaining chars in the fortran string with spaces */
  memcpy (fortran_string, *c_string, c_strlen);
  memset (fortran_string + c_strlen, ' ', fortran_length - c_strlen);
}
