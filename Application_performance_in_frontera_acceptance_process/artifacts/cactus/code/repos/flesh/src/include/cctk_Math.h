 /*@@
   @header    cctk_Math.h
   @date      2012-10-17
   @author    Erik Schnetter
   @desc
              Miscellaneous math routines, providing fallback C
              implementations for broken C++ compilers, and providing
              dummy implementations for broken C compilers.
   @enddesc
 @@*/

#ifndef _CCTK_MATH_H_
#define _CCTK_MATH_H_

#ifdef __cplusplus
/* Note: Some C++ implementations differ in their behaviour (how isnan
   etc. can be accessed) depending on the order (!) in which <cmath>
   and <math.h> is included, if both are included.
   
   We therefore require: A file including this file may include
   <cmath> at will, and is allowed to include <math.h> if <cmath> is
   included before. We do this because <cmath> seems "more standard"
   than <math.h> for C++ these days. */
#  include <cmath>
#  include <math.h>
#  include <cctk_Config.h>
#  include <cctk_Types.h>
#endif



#ifdef __cplusplus
extern "C" {
#endif

double CCTK_copysign(double x, double y);
int CCTK_fpclassify(double x);
int CCTK_isfinite(double x);
int CCTK_isinf(double x);
int CCTK_isnan(double x);
int CCTK_isnormal(double x);
int CCTK_signbit(double x);

  /* int CCTK_IEEE_fpclassify(double x); */
int CCTK_IEEE_isfinite(double x);
int CCTK_IEEE_isinf(double x);
int CCTK_IEEE_isnan(double x);
int CCTK_IEEE_isnormal(double x);
int CCTK_IEEE_signbit(double x);

#ifdef __cplusplus
}   
#endif



#if defined __cplusplus && ! defined __CUDACC__

/* If necessary, provide fallback C implementations for C++, using the
   routines declared above. */

#ifndef HAVE_CCTK_CXX_COPYSIGN
#  define HAVE_CCTK_CXX_COPYSIGN 1
#  define CCTK_CXX_COPYSIGN CCTK_copysign
#  define HAVE_CCTK_COPYSIGN HAVE_CCTK_CXX_COPYSIGN
#  define CCTK_COPYSIGN CCTK_CXX_COPYSIGN
#endif

#ifndef HAVE_CCTK_CXX_FPCLASSIFY
#  define HAVE_CCTK_CXX_FPCLASSIFY 1
#  define CCTK_CXX_FPCLASSIFY CCTK_fpclassify
#  define HAVE_CCTK_FPCLASSIFY HAVE_CCTK_CXX_FPCLASSIFY
#  define CCTK_FPCLASSIFY CCTK_CXX_FPCLASSIFY
#endif

#ifndef HAVE_CCTK_CXX_ISFINITE
#  define HAVE_CCTK_CXX_ISFINITE 1
#  define CCTK_CXX_ISFINITE CCTK_isfinite
#  define HAVE_CCTK_ISFINITE HAVE_CCTK_CXX_ISFINITE
#  define CCTK_ISFINITE CCTK_CXX_ISFINITE
#endif

#ifndef HAVE_CCTK_CXX_ISINF
#  define HAVE_CCTK_CXX_ISINF 1
#  define CCTK_CXX_ISINF CCTK_isinf
#  define HAVE_CCTK_ISINF HAVE_CCTK_CXX_ISINF
#  define CCTK_ISINF CCTK_CXX_ISINF
#endif

#ifndef HAVE_CCTK_CXX_ISNAN
#  define HAVE_CCTK_CXX_ISNAN 1
#  define CCTK_CXX_ISNAN CCTK_isnan
#  define HAVE_CCTK_ISNAN HAVE_CCTK_CXX_ISNAN
#  define CCTK_ISNAN CCTK_CXX_ISNAN
#endif

#ifndef HAVE_CCTK_CXX_ISNORMAL
#  define HAVE_CCTK_CXX_ISNORMAL 1
#  define CCTK_CXX_ISNORMAL CCTK_isnormal
#  define HAVE_CCTK_ISNORMAL HAVE_CCTK_CXX_ISNORMAL
#  define CCTK_ISNORMAL CCTK_CXX_ISNORMAL
#endif

#ifndef HAVE_CCTK_CXX_SIGNBIT
#  define HAVE_CCTK_CXX_SIGNBIT 1
#  define CCTK_CXX_SIGNBIT CCTK_signbit
#  define HAVE_CCTK_SIGNBIT HAVE_CCTK_CXX_SIGNBIT
#  define CCTK_SIGNBIT CCTK_CXX_SIGNBIT
#endif



/* Provide macros, so that these routines can be accessed in a simple
   manner. */

/* First, we define known-to-be-good implementations in a separate
   namespace "std::Cactus". These routines can be called from anywhere
   and will work fine, independent of whether they are implemented as
   macros, or in which (other) namespace they actually live (could be
   std or the global namespace).
   
   Then we define our own macros, making sure to undefine any possibly
   previously existing macros. To avoid infinite recursion when
   expanding the macros, we use a prefix "good_" for the functions in
   the Cactus namespace. */
  
#define IMPLEMENT_FUNCTIONS(T)                  \
                                                \
  inline T good_copysign(T x, T y)              \
  {                                             \
    return CCTK_COPYSIGN(x, y);                 \
  }                                             \
                                                \
  inline int good_fpclassify(T x)               \
  {                                             \
    return CCTK_FPCLASSIFY(x);                  \
  }                                             \
                                                \
  inline int good_isfinite(T x)                 \
  {                                             \
    return CCTK_ISFINITE((double)x);            \
  }                                             \
                                                \
  inline int good_isinf(T x)                    \
  {                                             \
    return CCTK_ISINF(x);                       \
  }                                             \
                                                \
  inline int good_isnan(T x)                    \
  {                                             \
    return CCTK_ISNAN(x);                       \
  }                                             \
                                                \
  inline int good_isnormal(T x)                 \
  {                                             \
    return CCTK_ISNORMAL(x);                    \
  }                                             \
                                                \
  inline int good_signbit(T x)                  \
  {                                             \
    return CCTK_SIGNBIT(x);                     \
  }
  
/* Define the functions in std::Cactus, so that the macro isnan (see
   below) can be called as std::isnan. */
namespace std {
  namespace Cactus {
    
#ifdef HAVE_CCTK_REAL4
  IMPLEMENT_FUNCTIONS(CCTK_REAL4)
#endif
#ifdef HAVE_CCTK_REAL8
  IMPLEMENT_FUNCTIONS(CCTK_REAL8)
#endif
#ifdef HAVE_CCTK_REAL16
  IMPLEMENT_FUNCTIONS(CCTK_REAL16)
#endif
    
  }
}

#undef IMPLEMENT_FUNCTIONS
  
/* Define the macros after the functions above, so that the functions
   above can refer to the standard (non-Cactus) definitions. */

#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#define copysign   Cactus::good_copysign
#define fpclassify Cactus::good_fpclassify
#define isfinite   Cactus::good_isfinite
#define isinf      Cactus::good_isinf
#define isnan      Cactus::good_isnan
#define isnormal   Cactus::good_isnormal
#define signbit    Cactus::good_signbit

#endif  /* #if defined __cplusplus && ! defined __CUDACC__ */



#endif  /* #ifndef _CCTK_MATH_H_ */
