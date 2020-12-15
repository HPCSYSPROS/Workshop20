 /*@@
   @header    cctk_types.h
   @date      Mon Jun 21 21:03:27 1999
   @author    Tom Goodale
   @desc 
   Defines the appropriate types based upon the precision.
   Should be included by cctk_Config.h .
   @enddesc 
 @@*/

#ifndef _CCTK_TYPES_H_
#define _CCTK_TYPES_H_

/* Make sure that cctk_config.h is available in case someone includes
 * this by hand.
 */
#include "cctk_Config.h"

/* Define stuff for C. */
#ifdef CCODE

typedef void *CCTK_POINTER;
typedef const void *CCTK_POINTER_TO_CONST;
typedef void (*CCTK_FPOINTER)(void);
#define HAVE_CCTK_POINTER 1
#define HAVE_CCTK_POINTER_TO_CONST 1
#define HAVE_CCTK_FPOINTER 1

/* Character types */
typedef char CCTK_CHAR;
typedef const char * CCTK_STRING;
#define HAVE_CCTK_CHAR 1
#define HAVE_CCTK_STRING 1

/* Floating point types */

#ifdef HAVE_CCTK_REAL16
typedef CCTK_REAL16_TYPE CCTK_REAL16;
#endif
#ifdef HAVE_CCTK_REAL8
typedef CCTK_REAL8_TYPE CCTK_REAL8;
#endif
#ifdef HAVE_CCTK_REAL4
typedef CCTK_REAL4_TYPE CCTK_REAL4;
#endif

/* Declarations for complex types */

#ifdef __cplusplus
#  include <complex>
#endif

#ifdef HAVE_CCTK_REAL16
#  define HAVE_CCTK_COMPLEX32 1
#  ifdef __cplusplus
typedef std::complex<CCTK_REAL16> CCTK_COMPLEX32;
#  else
typedef long double _Complex CCTK_COMPLEX32;
#  endif
#endif

#ifdef HAVE_CCTK_REAL8
#  define HAVE_CCTK_COMPLEX16 1
#  ifdef __cplusplus
typedef std::complex<CCTK_REAL8> CCTK_COMPLEX16;
#  else
typedef double _Complex CCTK_COMPLEX16;
#  endif
#endif

#ifdef HAVE_CCTK_REAL4
#  define HAVE_CCTK_COMPLEX8 1
#  ifdef __cplusplus
typedef std::complex<CCTK_REAL4> CCTK_COMPLEX8;
#  else
typedef float _Complex CCTK_COMPLEX8;
#  endif
#endif

/* Small positive integer type */
typedef unsigned char CCTK_BYTE;
#define HAVE_CCTK_BYTE 1

/* Integer types */

#ifndef __CUDACC__
#ifdef HAVE_CCTK_INT16
typedef CCTK_INT16_TYPE CCTK_INT16;
#endif
#endif
#ifdef HAVE_CCTK_INT8
typedef CCTK_INT8_TYPE CCTK_INT8;
#endif
#ifdef HAVE_CCTK_INT4
typedef CCTK_INT4_TYPE CCTK_INT4;
#endif
#ifdef HAVE_CCTK_INT2
typedef CCTK_INT2_TYPE CCTK_INT2;
#endif
#ifdef HAVE_CCTK_INT1
typedef CCTK_INT1_TYPE CCTK_INT1;
#endif

#endif /* CCODE */

/* Define stuff for Fortran. */
#ifdef FCODE

#define CCTK_POINTER          integer*SIZEOF_CHAR_P
#define CCTK_POINTER_TO_CONST integer*SIZEOF_CHAR_P
/* TODO: add autoconf for determining the size of function pointers */
#define CCTK_FPOINTER         integer*SIZEOF_CHAR_P
#define HAVE_CCTK_POINTER 1
#define HAVE_CCTK_POINTER_TO_CONST 1
#define HAVE_CCTK_FPOINTER 1

/* Character types */
/* A single character does not exist in Fortran; in Fortran, all
   character types are strings.  Hence we do not define CCTK_CHAR.  */
/* #define CCTK_CHAR   CHARACTER */
/* #define HAVE_CCTK_CHAR 1 */
/* This is a C-string, i.e., only a pointer */
#define CCTK_STRING CCTK_POINTER_TO_CONST
#define HAVE_CCTK_STRING 1

#ifdef HAVE_CCTK_INT16
#  define CCTK_INT16 INTEGER*16
#endif
#ifdef HAVE_CCTK_INT8
#  define CCTK_INT8 INTEGER*8
#endif
#ifdef HAVE_CCTK_INT4
#  define CCTK_INT4 INTEGER*4
#endif
#ifdef HAVE_CCTK_INT2
#  define CCTK_INT2 INTEGER*2
#endif
#ifdef HAVE_CCTK_INT1
#  define CCTK_INT1 INTEGER*1
#endif

#ifdef HAVE_CCTK_REAL16
#  define CCTK_REAL16 REAL*CCTK_REAL16_KIND
#  define HAVE_CCTK_COMPLEX32 1
#  define CCTK_COMPLEX32  COMPLEX*CCTK_COMPLEX32_KIND
#endif

#ifdef HAVE_CCTK_REAL8
#  define CCTK_REAL8  REAL*8
#  define HAVE_CCTK_COMPLEX16 1
#  define CCTK_COMPLEX16  COMPLEX*16
#endif

#ifdef HAVE_CCTK_REAL4
#  define CCTK_REAL4  REAL*4
#  define HAVE_CCTK_COMPLEX8 1
#  define CCTK_COMPLEX8   COMPLEX*8
#endif

/* Should be unsigned, but Fortran doesn't have that */
#define CCTK_BYTE INTEGER*1
#define HAVE_CCTK_BYTE 1

#endif /*FCODE */

/* Now pick the types based upon the precision variable. */

/* Floating point precision */
#ifdef CCTK_REAL_PRECISION_16
#  define CCTK_REAL_PRECISION 16
#  define CCTK_REAL CCTK_REAL16
#endif

#ifdef CCTK_REAL_PRECISION_8
#  define CCTK_REAL_PRECISION 8
#  define CCTK_REAL CCTK_REAL8
#endif

#ifdef CCTK_REAL_PRECISION_4
#  define CCTK_REAL_PRECISION 4
#  define CCTK_REAL CCTK_REAL4
#endif

/* Integer precision */

#ifdef CCTK_INTEGER_PRECISION_16
#  define CCTK_INTEGER_PRECISION 16
#  define CCTK_INT CCTK_INT16
#endif

#ifdef CCTK_INTEGER_PRECISION_8
#  define CCTK_INTEGER_PRECISION 8
#  define CCTK_INT CCTK_INT8
#endif

#ifdef CCTK_INTEGER_PRECISION_4
#  define CCTK_INTEGER_PRECISION 4
#  define CCTK_INT CCTK_INT4
#endif

#ifdef CCTK_INTEGER_PRECISION_2
#  define CCTK_INTEGER_PRECISION 2
#  define CCTK_INT CCTK_INT2
#endif

#ifdef CCTK_INTEGER_PRECISION_1
#  define CCTK_INTEGER_PRECISION 1
#  define CCTK_INT CCTK_INT1
#endif

/* Complex precision */
#ifdef CCTK_REAL_PRECISION_16
#  define CCTK_COMPLEX_PRECISION 32
#  define CCTK_COMPLEX CCTK_COMPLEX32
#endif

#ifdef CCTK_REAL_PRECISION_8
#  define CCTK_COMPLEX_PRECISION 16
#  define CCTK_COMPLEX CCTK_COMPLEX16
#endif

#ifdef CCTK_REAL_PRECISION_4
#  define CCTK_COMPLEX_PRECISION 8
#  define CCTK_COMPLEX CCTK_COMPLEX8
#endif

/* Determine whether we have a traditional or an ANSI cpp. */
#ifdef FCODE
/* The empty
   comment in the definition of CCTK_ANSI_FPP will either turn into
   nothing or into white space.  There must not be any add spaces
   around this empty comment.

   A traditional cpp will turn it into nothing, an ANSI cpp will turn
   it into white space.  Depending on this, CCTK_ANSI_FPP will either
   turn into a single separate token (which lead to the value 0), or
   into two separate tokens (which lead to the value 1).

   This is magic.  */
#define CCTKi_FPP_A
#define CCTKi_FPP_B 1
#define CCTKi_FPP_ACCTKi_FPP_B 0
#define CCTK_ANSI_FPP CCTKi_FPP_A/**/CCTKi_FPP_B
#endif

/* Handle 'unused' function arguments */
#ifdef FCODE
/* Declare a variable and tell the compiler that it may be unused.
   This is used for CCTK_ARGUMENTS.

   The macro CCTK_DECLARE (typ, nam, dim) is used with
   typ: a type, used to declare the variable (e.g. "CCTK_REAL")
   nam: the variable name (e.g. "x")
   dim: optional array dimensions, (e.g. "(10,10)")
*/
#ifdef F90CODE
/* Declare it, and use it for a dummy operation  */

#if CCTK_ANSI_FPP
#define CCTK_DECLARE(typ,nam,dim)                       \
  typ nam dim &&                                        \
  integer, parameter :: cctki_use_##nam = kind(nam)
#else
#define CCTK_DECLARE(typ,nam,dim)                       \
  typ nam dim &&                                        \
  integer, parameter :: cctki_use_/**/nam = kind(nam)
#endif

#else /* #ifndef F90CODE */

/* Just declare it; FORTRAN 77 has no good way of marking it as used
   within a block of declarations  */
#define CCTK_DECLARE(typ,nam,dim)               \
  typ nam dim

#endif /* F90CODE */
#endif /* FCODE */

#ifdef CCODE
/* Declare and initialise a variable and tell the compiler that it may
   be unused.  This is used for CCTK_PARAMETERS and CCTK_ARGUMENTS.

   The macro CCTK_DECLARE_INIT (typ, nam, val) is used with
   typ: a type, used to declare the variable (e.g. "const int")
   nam: the variable name (e.g. "x")
   val: the value used to initialise it (e.g. "42")
*/
#if  (! defined(__cplusplus) && defined(HAVE_CCTK_C_ATTRIBUTE_UNUSED  )) \
  || (  defined(__cplusplus) && defined(HAVE_CCTK_CXX_ATTRIBUTE_UNUSED))

/* We have __attribute__((unused)), so use it */
#define CCTK_DECLARE_INIT(typ,nam,val)          \
  typ nam CCTK_ATTRIBUTE_UNUSED = (val);

#else

/* Some fallback, bound to fool most compilers */
#define CCTK_DECLARE_INIT(typ,nam,val)                                  \
  typ nam = (val);                                                      \
  enum cctki_use_##nam { cctki_use0_##nam = sizeof nam };

#endif /* HAVE_..._ATTRIBUTE_UNUSED */

#endif /* CCODE */

#endif /*_CCTK_TYPES_H_ */
