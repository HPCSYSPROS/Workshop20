#ifndef FORTRAN_H
#define FORTRAN_H

#include "cctk.h"

#if FORTRAN_CPP_ANSI
/* For an ANSI-like cpp */

#define CONCAT(a,b) a##b
#define STRINGIFY(a) STRINGIFY_(a)
#define STRINGIFY_(a) #a

#else
/* For a traditional cpp */

#define CONCAT(a,b) a/**/b
#define STRINGIFY(a) "a"

#endif

#endif /* #ifdef FORTRAN_H */
