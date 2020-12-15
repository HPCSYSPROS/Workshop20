 /*@@
   @header    CactusMainFunctions.h
   @date      Thu Jan 14 18:34:52 1999
   @author    Tom Goodale
   @desc 
   Main registerable functions.
   @enddesc 
   @version $Header$
 @@*/

#ifndef _CACTUSMAINFUNCTIONS_H_
#define _CACTUSMAINFUNCTIONS_H_

#include "OverloadMacros.h"

/* Function prototypes */

#ifdef __cplusplus
extern "C" {
#endif

#define OVERLOADABLE(name) OVERLOADABLE_PROTOTYPE(name)

#include "MainOverloadables.h"

#undef OVERLOADABLE

#ifdef __cplusplus
}
#endif

#endif
