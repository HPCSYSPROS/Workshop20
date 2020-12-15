 /*@@
   @header    cctk_Complex.h
   @date      Tue Dec 14 12:28:05 1999
   @author    Tom Goodale
   @desc
              Prototypes for complex numbers.
   @enddesc
   @version   $Header$
 @@*/

#ifndef _CCTK_COMPLEX_H_
#define _CCTK_COMPLEX_H_

#include "cctk_Config.h"

#ifdef __cplusplus

#include <complex>
/* Macro to declare a set of complex functions for a given precision */
/* Note: We do not declare these as extern "C" since that they are
   different from their C counterparts */
#define DECLARE_CMPLX_FUNCTIONS(CCTK_Cmplx, cctk_real, cctk_complex)          \
static inline cctk_complex CCTK_Cmplx(cctk_real Re, cctk_real Im) {           \
  return cctk_complex(Re,Im);                                                 \
}                                                                             \
static inline cctk_real    CCTK_Cmplx##Real(cctk_complex a) {                 \
  return std::real(a);                                                        \
}                                                                             \
static inline cctk_real    CCTK_Cmplx##Imag(cctk_complex a) {                 \
  return std::imag(a);                                                        \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Neg(cctk_complex a) {                  \
  return -a;                                                                  \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Conjg(cctk_complex a) {                \
  return std::conj(a);                                                        \
}                                                                             \
static inline cctk_real    CCTK_Cmplx##Abs(cctk_complex a) {                  \
  return std::abs(a);                                                         \
}                                                                             \
static inline cctk_real    CCTK_Cmplx##Arg(cctk_complex a) {                  \
  return std::arg(a);                                                         \
}                                                                             \
static inline cctk_real    CCTK_Cmplx##Norm(cctk_complex a) {                 \
  return std::norm(a);                                                        \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Add(cctk_complex a, cctk_complex b) {  \
  return a+b;                                                                 \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Sub(cctk_complex a, cctk_complex b) {  \
  return a-b;                                                                 \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Mul(cctk_complex a, cctk_complex b) {  \
  return a*b;                                                                 \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Div(cctk_complex a, cctk_complex b) {  \
  return a/b;                                                                 \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##CPow(cctk_complex a, cctk_complex b) { \
  return pow(a,b);                                                            \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Sin(cctk_complex a) {                  \
  return std::sin(a);                                                         \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Cos(cctk_complex a) {                  \
  return std::cos(a);                                                         \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Exp(cctk_complex a) {                  \
  return std::exp(a);                                                         \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Log(cctk_complex a) {                  \
  return std::log(a);                                                         \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Sqrt(cctk_complex a) {                 \
  return std::sqrt(a);                                                        \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##Pow(cctk_complex a, cctk_real b) {     \
  return std::pow(a,b);                                                       \
}                                                                             \
static inline cctk_complex CCTK_Cmplx##IPow(cctk_complex a, int b) {          \
  return std::pow(a,(cctk_real)b);                                            \
}

#else

/* Macro to declare a set of complex functions for a given precision */
/* Note: We do not provide inline implementations here, since this
   would require including <complex.h>, which pushes the identifier
   "I" into the global namespace. */
#define DECLARE_CMPLX_FUNCTIONS(CCTK_Cmplx, cctk_real, cctk_complex)    \
cctk_complex CCTK_Cmplx        (cctk_real Re, cctk_real Im);            \
cctk_real    CCTK_Cmplx##Real  (cctk_complex a);                        \
cctk_real    CCTK_Cmplx##Imag  (cctk_complex a);                        \
cctk_complex CCTK_Cmplx##Neg   (cctk_complex a);                        \
cctk_complex CCTK_Cmplx##Conjg (cctk_complex a);                        \
cctk_real    CCTK_Cmplx##Abs   (cctk_complex a);                        \
cctk_real    CCTK_Cmplx##Arg   (cctk_complex a);                        \
cctk_real    CCTK_Cmplx##Norm  (cctk_complex a);                        \
cctk_complex CCTK_Cmplx##Add   (cctk_complex a, cctk_complex b);        \
cctk_complex CCTK_Cmplx##Sub   (cctk_complex a, cctk_complex b);        \
cctk_complex CCTK_Cmplx##Mul   (cctk_complex a, cctk_complex b);        \
cctk_complex CCTK_Cmplx##Div   (cctk_complex a, cctk_complex b);        \
cctk_complex CCTK_Cmplx##CPow  (cctk_complex a, cctk_complex b);        \
cctk_complex CCTK_Cmplx##Sin   (cctk_complex a);                        \
cctk_complex CCTK_Cmplx##Cos   (cctk_complex a);                        \
cctk_complex CCTK_Cmplx##Exp   (cctk_complex a);                        \
cctk_complex CCTK_Cmplx##Log   (cctk_complex a);                        \
cctk_complex CCTK_Cmplx##Sqrt  (cctk_complex a);                        \
cctk_complex CCTK_Cmplx##Pow   (cctk_complex a, cctk_real b);           \
cctk_complex CCTK_Cmplx##IPow  (cctk_complex a, int b);
#endif


/* declare complex functions for all available precisions */
#ifdef HAVE_CCTK_REAL4
DECLARE_CMPLX_FUNCTIONS (CCTK_Cmplx8, CCTK_REAL4, CCTK_COMPLEX8)
#endif

#ifdef HAVE_CCTK_REAL8
DECLARE_CMPLX_FUNCTIONS (CCTK_Cmplx16, CCTK_REAL8, CCTK_COMPLEX16)
#endif

#ifdef HAVE_CCTK_REAL16
DECLARE_CMPLX_FUNCTIONS (CCTK_Cmplx32, CCTK_REAL16, CCTK_COMPLEX32)
#endif


/* declare the default precision complex functions as #define'd macros */
#ifdef CCTK_REAL_PRECISION_4
#define CCTK_Cmplx      CCTK_Cmplx8
#define CCTK_CmplxReal  CCTK_Cmplx8Real
#define CCTK_CmplxImag  CCTK_Cmplx8Imag
#define CCTK_CmplxNeg   CCTK_Cmplx8Neg
#define CCTK_CmplxConjg CCTK_Cmplx8Conjg
#define CCTK_CmplxAbs   CCTK_Cmplx8Abs
#define CCTK_CmplxArg   CCTK_Cmplx8Arg
#define CCTK_CmplxNorm  CCTK_Cmplx8Norm
#define CCTK_CmplxAdd   CCTK_Cmplx8Add
#define CCTK_CmplxSub   CCTK_Cmplx8Sub
#define CCTK_CmplxMul   CCTK_Cmplx8Mul
#define CCTK_CmplxDiv   CCTK_Cmplx8Div
#define CCTK_CmplxCPow  CCTK_Cmplx8CPow
#define CCTK_CmplxSin   CCTK_Cmplx8Sin
#define CCTK_CmplxCos   CCTK_Cmplx8Cos
#define CCTK_CmplxExp   CCTK_Cmplx8Exp
#define CCTK_CmplxLog   CCTK_Cmplx8Log
#define CCTK_CmplxSqrt  CCTK_Cmplx8Sqrt
#define CCTK_CmplxPow   CCTK_Cmplx8Pow
#define CCTK_CmplxIPow  CCTK_Cmplx8IPow
#elif CCTK_REAL_PRECISION_8
#define CCTK_Cmplx      CCTK_Cmplx16
#define CCTK_CmplxReal  CCTK_Cmplx16Real
#define CCTK_CmplxImag  CCTK_Cmplx16Imag
#define CCTK_CmplxNeg   CCTK_Cmplx16Neg
#define CCTK_CmplxConjg CCTK_Cmplx16Conjg
#define CCTK_CmplxAbs   CCTK_Cmplx16Abs
#define CCTK_CmplxArg   CCTK_Cmplx16Arg
#define CCTK_CmplxNorm  CCTK_Cmplx16Norm
#define CCTK_CmplxAdd   CCTK_Cmplx16Add
#define CCTK_CmplxSub   CCTK_Cmplx16Sub
#define CCTK_CmplxMul   CCTK_Cmplx16Mul
#define CCTK_CmplxDiv   CCTK_Cmplx16Div
#define CCTK_CmplxCPow  CCTK_Cmplx16CPow
#define CCTK_CmplxSin   CCTK_Cmplx16Sin
#define CCTK_CmplxCos   CCTK_Cmplx16Cos
#define CCTK_CmplxExp   CCTK_Cmplx16Exp
#define CCTK_CmplxLog   CCTK_Cmplx16Log
#define CCTK_CmplxSqrt  CCTK_Cmplx16Sqrt
#define CCTK_CmplxPow   CCTK_Cmplx16Pow
#define CCTK_CmplxIPow  CCTK_Cmplx16IPow
#elif CCTK_REAL_PRECISION_16
#define CCTK_Cmplx      CCTK_Cmplx32
#define CCTK_CmplxReal  CCTK_Cmplx32Real
#define CCTK_CmplxImag  CCTK_Cmplx32Imag
#define CCTK_CmplxNeg   CCTK_Cmplx32Neg
#define CCTK_CmplxConjg CCTK_Cmplx32Conjg
#define CCTK_CmplxAbs   CCTK_Cmplx32Abs
#define CCTK_CmplxArg   CCTK_Cmplx32Arg
#define CCTK_CmplxNorm  CCTK_Cmplx32Norm
#define CCTK_CmplxAdd   CCTK_Cmplx32Add
#define CCTK_CmplxSub   CCTK_Cmplx32Sub
#define CCTK_CmplxMul   CCTK_Cmplx32Mul
#define CCTK_CmplxDiv   CCTK_Cmplx32Div
#define CCTK_CmplxCPow  CCTK_Cmplx32CPow
#define CCTK_CmplxSin   CCTK_Cmplx32Sin
#define CCTK_CmplxCos   CCTK_Cmplx32Cos
#define CCTK_CmplxExp   CCTK_Cmplx32Exp
#define CCTK_CmplxLog   CCTK_Cmplx32Log
#define CCTK_CmplxSqrt  CCTK_Cmplx32Sqrt
#define CCTK_CmplxPow   CCTK_Cmplx32Pow
#define CCTK_CmplxIPow  CCTK_Cmplx32IPow
#endif

#endif /* _CCTK_COMPLEX_H_ */
