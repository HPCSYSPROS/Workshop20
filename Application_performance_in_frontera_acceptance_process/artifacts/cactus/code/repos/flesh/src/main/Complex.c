 /*@@
   @file      Complex.c
   @date      Tue Dec 14 12:09:43 1999
   @author    Tom Goodale
   @desc
              Complex variable stuff
   @enddesc
   @version   $Id$
 @@*/

#include <complex.h>
#include <math.h>

#include "cctk_Flesh.h"

#ifndef DEFINE_CCTK_COMPLEX_INLINE_FUNCTIONS
#  define DEFINE_CCTK_COMPLEX_EXTERN_FUNCTIONS
#  include "cctk_Complex.h"
#  undef DEFINE_CCTK_COMPLEX_EXTERN_FUNCTIONS
#endif


#ifndef DEFINE_CCTK_COMPLEX_INLINE_FUNCTIONS

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(main_Complex_c);

#endif


/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

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
   @routine    CCTK_Cmplx
   @date       Tue Dec 14 12:16:01 1999
   @author     Tom Goodale
   @desc
               Turns two reals into a complex number
   @enddesc
 
   @var        Re
   @vdesc      Real part
   @vtype      CCTK_REAL
   @vio        in
   @endvar
   @var        Im
   @vdesc      Imaginary part
   @vtype      CCTK_REAL
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX(CCTK_Cmplx, cctk_real, cctk_complex, type)    \
cctk_complex CCTK_Cmplx (cctk_real Re, cctk_real Im)                    \
{                                                                       \
  return Re + _Complex_I * Im;                                          \
}


 /*@@
   @routine    CCTK_CmplxReal
   @date       Tue Dec 14 12:22:41 1999
   @author     Tom Goodale
   @desc
               Returns the real part of a complex number.
   @enddesc

   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The real part
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_REAL(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_real CCTK_Cmplx##Real (cctk_complex x)                             \
{                                                                       \
  return creal##type(x);                                                \
}


 /*@@
   @routine    CCTK_CmplxImag
   @date       Tue Dec 14 12:22:41 1999
   @author     Tom Goodale
   @desc
               Returns the imaginary part of a complex number.
   @enddesc

   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The imaginary part
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_IMAG(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_real CCTK_Cmplx##Imag (cctk_complex x)                             \
{                                                                       \
  return cimag##type(x);                                                \
}


 /*@@
   @routine    CCTK_CmplxNeg
   @date       2006-06-07
   @author     Erik Schnetter
   @desc
               Returns the negative of a complex number.
   @enddesc
 
   @var       x
   @vdesc     The complex number
   @vtype     CCTK_COMPLEX
   @vio       in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The negative
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_NEG(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Neg (cctk_complex x)                           \
{                                                                       \
  return -x;                                                            \
}


 /*@@
   @routine    CCTK_CmplxConjg
   @date       Tue Dec 14 12:22:41 1999
   @author     Tom Goodale
   @desc
               Returns the complex conjugate of a complex number.
   @enddesc
 
   @var       x
   @vdesc     The complex number
   @vtype     CCTK_COMPLEX
   @vio       in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The complex conjugate
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_CONJG(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Conjg (cctk_complex x)                         \
{                                                                       \
  return conj##type(x);                                                 \
}


 /*@@
   @routine    CCTK_CmplxAbs
   @date       Tue Dec 14 12:26:33 1999
   @author     Tom Goodale
   @desc
               Return the absolute value of a complex number.
   @enddesc
 
   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The absolute value of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_ABS(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_real CCTK_Cmplx##Abs (cctk_complex x)                              \
{                                                                       \
  return cabs##type(x);                                                 \
}


 /*@@
   @routine    CCTK_CmplxNorm
   @date       2006-07-06
   @author     Erik Schnetter
   @desc
               Return the absolute value squared of a complex number.
   @enddesc
 
   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The absolute value squared of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_NORM(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_real CCTK_Cmplx##Norm (cctk_complex x)                             \
{                                                                       \
  return creal##type(x)*creal##type(x) + cimag##type(x)*cimag##type(x); \
}


 /*@@
   @routine    CCTK_CmplxArg
   @date       2005-11-17
   @author     Erik Schnetter
   @desc
               Return the argument of a complex number.
   @enddesc
 
   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_REAL
   @returndesc
               The argument of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_ARG(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_real CCTK_Cmplx##Arg (cctk_complex x)                              \
{                                                                       \
  return carg##type(x);                                                 \
}


 /*@@
   @routine    CCTK_CmplxAdd
   @date       Sat Dec 4 12:11:04 1999
   @author     Gabrielle Allen
   @desc
               Adds two complex numbers
   @enddesc
 
   @var        a
   @vdesc      First summand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      Second summand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The sum of a and b
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_ADD(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Add (cctk_complex a, cctk_complex b)           \
{                                                                       \
  return a+b;                                                           \
}


 /*@@
   @routine    CCTK_CmplxSub
   @date       Sat Dec 4 12:11:04 1999
   @author     Gabrielle Allen
   @desc
               Subtracts two complex numbers
   @enddesc
 
   @var        a
   @vdesc      First operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      Second operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The difference
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_SUB(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Sub (cctk_complex a, cctk_complex b)           \
{                                                                       \
  return a-b;                                                           \
}


 /*@@
   @routine    CCTK_CmplxMul
   @date       Sat Dec 4 12:11:04 1999
   @author     Gabrielle Allen
   @desc
               Multiplies two complex numbers
   @enddesc
 
   @var        a
   @vdesc      First operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      Second operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The product
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_MUL(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Mul (cctk_complex a, cctk_complex b)           \
{                                                                       \
  return a*b;                                                           \
}


 /*@@
   @routine    CCTK_CmplxDiv
   @date       Sat Dec 4 12:11:04 1999
   @author     Gabrielle Allen
   @desc
               Divide two complex numbers.
   @enddesc
 
   @var        a
   @vdesc      First operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      Second operand
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The quotient
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_DIV(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Div (cctk_complex a, cctk_complex b)           \
{                                                                       \
  return a/b;                                                           \
}


 /*@@
   @routine    CCTK_CmplxCPow
   @date       2005-11-17
   @author     Erik Schnetter
   @desc
               Raises a complex number to a given power
   @enddesc
 
   @var        a
   @vdesc      The base
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        b
   @vdesc      The exponent
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               a to the power of b
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_CPOW(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##CPow (cctk_complex a, cctk_complex b)          \
{                                                                       \
  return cpow##type(a, b);                                              \
}


 /*@@
   @routine    CCTK_CmplxSin
   @date       Wed 12 Dec 2001
   @author     Thomas Radke
   @desc
               Returns the sine of a complex number.
   @enddesc
 
   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The sine
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_SIN(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Sin (cctk_complex x)                           \
{                                                                       \
  return csin##type(x);                                                 \
}


 /*@@
   @routine    CCTK_CmplxCos
   @date       Wed 12 Dec 2001
   @author     Thomas Radke
   @desc
               Returns the cosine of a complex number.
   @enddesc
 
   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The cosine
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_COS(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Cos (cctk_complex x)                           \
{                                                                       \
  return ccos##type(x);                                                 \
}


 /*@@
   @routine    CCTK_CmplxExp
   @date       Wed 12 Dec 2001
   @author     Thomas Radke
   @desc
               Returns the exponential of a complex number.
   @enddesc
 
   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The exponential
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_EXP(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Exp (cctk_complex x)                           \
{                                                                       \
  return cexp##type(x);                                                 \
}


 /*@@
   @routine    CCTK_CmplxLog
   @date       2005-11-17
   @author     Erik Schnetter
   @desc
               Returns the logarithm of a complex number.
   @enddesc
 
   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The logarithm
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_LOG(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Log (cctk_complex x)                           \
{                                                                       \
  return clog##type(x);                                                 \
}


 /*@@
   @routine    CCTK_CmplxSqrt
   @date       Wed 12 Dec 2001
   @author     Thomas Radke
   @desc
               Returns the square root of a complex number.
   @enddesc

   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The square root
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_SQRT(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Sqrt (cctk_complex x)                          \
{                                                                       \
  return csqrt##type(x);                                                \
}


 /*@@
   @routine    CCTK_CmplxPow
   @date       
   @author     Yaakoub Y El Khamra
   @desc
               Raises a complex number to a given power
   @enddesc
 
   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        w
   @vdesc      The exponent
   @vtype      CCTK_REAL
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The power of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_POW(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##Pow (cctk_complex x, cctk_real w)              \
{                                                                       \
  return cpow##type(x, w);                                              \
}


 /*@@
   @routine    CCTK_CmplxIPow
   @date       
   @author     Erik Schnetter
   @desc
               Raises a complex number to a given power
   @enddesc
 
   @var        x
   @vdesc      The complex number
   @vtype      CCTK_COMPLEX
   @vio        in
   @endvar
   @var        w
   @vdesc      The exponent
   @vtype      int
   @vio        in
   @endvar

   @returntype CCTK_COMPLEX
   @returndesc
               The power of the complex number
   @endreturndesc
@@*/
#define DEFINE_CCTK_CMPLX_IPOW(CCTK_Cmplx, cctk_real, cctk_complex, type) \
cctk_complex CCTK_Cmplx##IPow (cctk_complex x, int w)                   \
{                                                                       \
  int w0 = w;                                                           \
  cctk_complex result = 1;                                              \
  while (w)                                                             \
  {                                                                     \
    if (w % 2)                                                          \
    {                                                                   \
      result *= x;                                                      \
    }                                                                   \
    w /= 2;                                                             \
    x *= x;                                                             \
  }                                                                     \
  if (w0 < 0)                                                           \
  {                                                                     \
    result = 1 / result;                                                \
  }                                                                     \
  return result;                                                        \
}



/* macro to define a set of complex functions for a given precision */
#define DEFINE_CMPLX_FUNCTIONS(CCTK_Cmplx, cctk_real, cctk_complex, type) \
DEFINE_CCTK_CMPLX       (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_REAL  (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_IMAG  (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_NEG   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_CONJG (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_ABS   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_NORM  (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_ARG   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_ADD   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_SUB   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_MUL   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_DIV   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_CPOW  (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_SIN   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_COS   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_EXP   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_LOG   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_SQRT  (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_POW   (CCTK_Cmplx, cctk_real, cctk_complex, type)     \
DEFINE_CCTK_CMPLX_IPOW  (CCTK_Cmplx, cctk_real, cctk_complex, type)


/* define complex functions for all available precisions */
#ifdef HAVE_CCTK_REAL4
  #define KIND_SUFFIX f
  DEFINE_CMPLX_FUNCTIONS (CCTK_Cmplx8, CCTK_REAL4, CCTK_COMPLEX8, KIND_SUFFIX)
  #undef KIND_SUFFIX
#endif

#ifdef HAVE_CCTK_REAL8
  #define KIND_SUFFIX
  DEFINE_CMPLX_FUNCTIONS (CCTK_Cmplx16, CCTK_REAL8, CCTK_COMPLEX16, KIND_SUFFIX)
  #undef KIND_SUFFIX
#endif

#ifdef HAVE_CCTK_REAL16
  #define KIND_SUFFIX l
  DEFINE_CMPLX_FUNCTIONS (CCTK_Cmplx32, CCTK_REAL16, CCTK_COMPLEX32, KIND_SUFFIX)
  #undef KIND_SUFFIX
#endif
