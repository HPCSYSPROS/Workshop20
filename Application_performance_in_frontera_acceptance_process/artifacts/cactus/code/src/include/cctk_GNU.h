 /*@@
   @header    cctk_GNU.h
   @date      Tue Feb 15 16:37:50 2000
   @author    Tom Goodale
   @desc 
   Includes various gnu header files
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CCTK_GNU_H_
#define _CCTK_GNU_H_ 1

#include "cctk_Config.h"

#ifdef HAVE_GETOPT_H
#ifdef HAVE_GETOPT_LONG_ONLY
#include <getopt.h>
#else
#include "../gnu/getopt.h"
#endif
#else
#include "../gnu/getopt.h"
#endif

#ifdef HAVE_REGEX_H
#include <regex.h>
#else
#include "../gnu/regex.h"
#endif

#endif /* _CCTK_GNU_H_ */
