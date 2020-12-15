 /*@@
   @header    CactusCommFunctions.h
   @date      Thu Jan 14 18:37:58 1999
   @author    Tom Goodale
   @desc 
   Overloadable communication functions
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CACTUSCOMMFUNCTIONS_H_
#define _CACTUSCOMMFUNCTIONS_H_

#include <stdarg.h>

#include "OverloadMacros.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The functions. */

#define OVERLOADABLE(name) OVERLOADABLE_PROTOTYPE(name)

#include "CommOverloadables.h"

#undef OVERLOADABLE

#ifdef __cplusplus
}
#endif

#endif
