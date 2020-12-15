// -*-C++-*-

#ifndef VECTORS_H
#define VECTORS_H

#include <cctk.h>



#define vec_static_assert(x) namespace { typedef int vsa[(x) ? 1 : -1]; }



#if VECTORISE

#  if defined __AVX__ && !defined DISABLE_AVX // Intel AVX
#    include "vectors-4-AVX.h"
#  elif defined __SSE__         // Intel SSE
#    include "vectors-4-SSE.h"
#  elif defined __ALTIVEC__     // Power Altivec
#    include "vectors-4-Altivec.h"
#  endif

#  if defined __MIC__           // Intel MIC
#    include "vectors-8-MIC.h"
#  elif defined __AVX__ && !defined DISABLE_AVX // Intel AVX
#    include "vectors-8-AVX.h"
#  elif defined __SSE2__        // Intel SSE2
#    include "vectors-8-SSE2.h"
#  elif defined __bgq__ && defined __VECTOR4DOUBLE__ // Blue Gene/Q QPX
#    include "vectors-8-QPX.h"
#  elif defined __ALTIVEC__ && defined _ARCH_PWR7 // Power VSX
#    include "vectors-8-VSX.h"
#  elif defined _ARCH_450D      // Blue Gene/P Double Hummer
#    include "vectors-8-DoubleHummer.h"
#  endif

#endif

// Default implementation, do not vectorise
#ifndef CCTK_REAL4_VEC_SIZE
#  include "vectors-4-default.h"
#endif
#ifndef CCTK_REAL8_VEC_SIZE
#  include "vectors-8-default.h"
#endif



// Operation counters
#ifndef VEC_COUNT
#  define VEC_COUNT(x)
#endif
// This expects variables declared as
//    ptrdiff_t vec_op_counter, vec_mem_counter;
#define vec_op_inc  ((void)(VEC_COUNT(vec_op_counter+=CCTK_REAL_VEC_SIZE)+0))
#define vec_mem_inc ((void)(VEC_COUNT(vec_mem_counter+=CCTK_REAL_VEC_SIZE)+0))



// Define macros for CCTK_REAL

#if defined CCTK_REAL_PRECISION_4 

#  define vec_architecture vec4_architecture

#  define CCTK_REAL_VEC      CCTK_REAL4_VEC
#  define CCTK_REAL_VEC_SIZE CCTK_REAL4_VEC_SIZE
#  define CCTK_INTEGER       CCTK_INTEGER4
#  define CCTK_BOOLEAN       CCTK_BOOLEAN4
#  define CCTK_INTEGER_VEC   CCTK_INTEGER4_VEC
#  define CCTK_BOOLEAN_VEC   CCTK_BOOLEAN4_VEC

#  define vec_set1 vec4_set1
#  define vec_set  vec4_set

#  define vec_elt  vec4_elt
#  define vec_elti vec4_elti
#  define vec_eltb vec4_eltb

#  define vec_load(p)            (vec_mem_inc, vec4_load(p))
#  define vec_loadu(p)           (vec_mem_inc, vec4_loadu(p))
#  define vec_loadu_maybe(off,p) (vec_mem_inc, vec4_loadu_maybe(off,p))
#  define vec_loadu_maybe3(off1,off2,off3,p)            \
  (vec_mem_inc, vec4_loadu_maybe3(off1,off2,off3,p))
#  define vec_store(p,x)             (vec_mem_inc, vec4_store(p,x))
#  define vec_store_nta(p,x)         (vec_mem_inc, vec4_store_nta(p,x))
#  define vec_store_partial_prepare  vec4_store_partial_prepare
#  define vec_store_nta_partial(p,x) (vec_mem_inc, vec4_store_nta_partial(p,x))
#  define vec_storeu_partial(p,x)    (vec_mem_inc, vec4_storeu_partial(p,x))
#  define vec_store_nta_partial_lo   vec4_store_nta_partial_lo
#  define vec_store_nta_partial_hi   vec4_store_nta_partial_hi
#  define vec_store_nta_partial_mid  vec4_store_nta_partial_mid

#  define kneg(x) (vec_op_inc, k4neg(x))

#  define kadd(x,y) (vec_op_inc, k4add(x,y))
#  define ksub(x,y) (vec_op_inc, k4sub(x,y))
#  define kmul(x,y) (vec_op_inc, k4mul(x,y))
#  define kdiv(x,y) (vec_op_inc, k4div(x,y))

#  define kmadd(x,y,z)  (vec_op_inc, vec_op_inc, k4madd(x,y,z))
#  define kmsub(x,y,z)  (vec_op_inc, vec_op_inc, k4msub(x,y,z))
#  define knmadd(x,y,z) (vec_op_inc, vec_op_inc, k4nmadd(x,y,z))
#  define knmsub(x,y,z) (vec_op_inc, vec_op_inc, k4nmsub(x,y,z))

#  define kacos     k4acos
#  define kacosh    k4acosh
#  define kasin     k4asin
#  define kasinh    k4asinh
#  define katan     k4atan
#  define katan2    k4atan2
#  define katanh    k4atanh
#  define kcopysign(x,y) (vec_op_inc, k4copysign(x,y))
#  define kcos      k4cos
#  define kcosh     k4cosh
#  define kexp      k4exp
#  define kfabs(x)   (vec_op_inc, k4fabs(x))
#  define kfmax(x,y) (vec_op_inc, k4fmax(x,y))
#  define kfmin(x,y) (vec_op_inc, k4fmin(x,y))
#  define kfmod(x,y) (vec_op_inc, k4fmod(x,y))
#  define kfnabs(x)  (vec_op_inc, k4fnabs(x))
#  define klog      k4log
#  define kpow      k4pow
#  define ksignbit(x) (vec_op_inc, k4signbit(x))
#  define ksin      k4sin
#  define ksinh     k4sinh
#  define ksgn      k4sgn
#  define ksqrt     k4sqrt
#  define ktan      k4tan
#  define ktanh     k4tanh

#  define klfalse k4lfalse
#  define kltrue  k4ltrue
#  define klnot   k4lnot
#  define kland   k4land
#  define klor    k4lor
#  define klxor   k4lxor
#  define kifthen k4ifthen

#  define kcmpeq k4cmpeq
#  define kcmpne k4cmpne
#  define kcmpgt k4cmpgt
#  define kcmpge k4cmpge
#  define kcmplt k4cmplt
#  define kcmple k4cmple

#elif defined CCTK_REAL_PRECISION_8 

#  define vec_architecture vec8_architecture

#  define CCTK_REAL_VEC      CCTK_REAL8_VEC
#  define CCTK_REAL_VEC_SIZE CCTK_REAL8_VEC_SIZE
#  define CCTK_INTEGER       CCTK_INTEGER8
#  define CCTK_BOOLEAN       CCTK_BOOLEAN8
#  define CCTK_INTEGER_VEC   CCTK_INTEGER8_VEC
#  define CCTK_BOOLEAN_VEC   CCTK_BOOLEAN8_VEC

#  define vec_set1 vec8_set1
#  define vec_set  vec8_set

#  define vec_elt  vec8_elt
#  define vec_elti vec8_elti
#  define vec_eltb vec8_eltb

#  define vec_load(p)                (vec_mem_inc, vec8_load(p))
#  define vec_loadu(p)               (vec_mem_inc, vec8_loadu(p))
#  define vec_loadu_maybe(off,p)     (vec_mem_inc, vec8_loadu_maybe(off,p))
#  define vec_loadu_maybe3(off1,off2,off3,p)            \
  (vec_mem_inc, vec8_loadu_maybe3(off1,off2,off3,p))
#  define vec_store(p,x)             (vec_mem_inc, vec8_store(p,x))
#  define vec_store_nta(p,x)         (vec_mem_inc, vec8_store_nta(p,x))
#  define vec_store_partial_prepare  vec8_store_partial_prepare
#  define vec_store_partial_prepare_fixed       \
  vec8_store_partial_prepare_fixed
#  define vec_store_nta_partial(p,x) (vec_mem_inc, vec8_store_nta_partial(p,x))
#  define vec_storeu_partial(p,x)    (vec_mem_inc, vec8_storeu_partial(p,x))
#  define vec_store_nta_partial_lo   vec8_store_nta_partial_lo
#  define vec_store_nta_partial_hi   vec8_store_nta_partial_hi
#  define vec_store_nta_partial_mid  vec8_store_nta_partial_mid

#  define kneg(x) (vec_op_inc, k8neg(x))

#  define kadd(x,y) (vec_op_inc, k8add(x,y))
#  define ksub(x,y) (vec_op_inc, k8sub(x,y))
#  define kmul(x,y) (vec_op_inc, k8mul(x,y))
#  define kdiv(x,y) (vec_op_inc, k8div(x,y))

#  define kmadd(x,y,z)  (vec_op_inc, vec_op_inc, k8madd(x,y,z))
#  define kmsub(x,y,z)  (vec_op_inc, vec_op_inc, k8msub(x,y,z))
#  define knmadd(x,y,z) (vec_op_inc, vec_op_inc, k8nmadd(x,y,z))
#  define knmsub(x,y,z) (vec_op_inc, vec_op_inc, k8nmsub(x,y,z))

#  define kacos     k8acos
#  define kacosh    k8acosh
#  define kasin     k8asin
#  define kasinh    k8asinh
#  define katan     k8atan
#  define katan2    k8atan2
#  define katanh    k8atanh
#  define kcopysign(x,y) (vec_op_inc, k8copysign(x,y))
#  define kcos      k8cos
#  define kcosh     k8cosh
#  define kexp      k8exp
#  define kfabs(x)   (vec_op_inc, k8fabs(x))
#  define kfmax(x,y) (vec_op_inc, k8fmax(x,y))
#  define kfmin(x,y) (vec_op_inc, k8fmin(x,y))
#  define kfmod(x,y) (vec_op_inc, k8fmod(x,y))
#  define kfnabs(x)  (vec_op_inc, k8fnabs(x))
#  define klog      k8log
#  define kpow      k8pow
#  define ksignbit(x) (vec_op_inc, k8signbit(x))
#  define ksin      k8sin
#  define ksinh     k8sinh
#  define ksgn      k8sgn
#  define ksqrt     k8sqrt
#  define ktan      k8tan
#  define ktanh     k8tanh

#  define klfalse k8lfalse
#  define kltrue  k8ltrue
#  define klnot   k8lnot
#  define kland   k8land
#  define klor    k8lor
#  define klxor   k8lxor
#  define kifthen k8ifthen

#  define kcmpeq k8cmpeq
#  define kcmpne k8cmpne
#  define kcmpgt k8cmpgt
#  define kcmpge k8cmpge
#  define kcmplt k8cmplt
#  define kcmple k8cmple

#else

#  error "Unknown CCTK_REAL_PRECISION"

#endif



// Deprecated
#define kifmsb(a,b,c) kifthen(a,b,c)
#define kifneg(a,b,c) kifmsb(a,b,c)
#define kifpos(a,b,c) kifmsb(a,c,b)

#define kisgn(a) (-42424242)



#if CCTK_REAL_VEC_SIZE == 1
#  define vec_index vec_set(0)
#elif CCTK_REAL_VEC_SIZE == 2
#  define vec_index vec_set(0,1)
#elif CCTK_REAL_VEC_SIZE == 4
#  define vec_index vec_set(0,1,2,3)
#elif CCTK_REAL_VEC_SIZE == 8
#  define vec_index vec_set(0,1,2,3,4,5,6,7)
#elif CCTK_REAL_VEC_SIZE == 16
#  define vec_index vec_set(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)
#else
#  error "Unsupported vector size"
#endif

  

// Define a class template for easier access from C++

#ifdef __cplusplus

#include <cmath>
#include <cstdlib>

template<typename T>
struct vecprops {
  typedef T scalar_t;
  typedef T vector_t;
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  std::size_t size()
  {
    return 1;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t load (scalar_t const& a)
  {
    return a;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t loadu (scalar_t const& a)
  {
    return a;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t set1 (scalar_t const& a)
  {
    return a;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  scalar_t elt (vector_t const& x, std::ptrdiff_t const d)
  {
    return x;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t neg (vector_t const& x)
  {
    return -x;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t add (vector_t const& x, vector_t const& y)
  {
    return x+y;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t sub (vector_t const& x, vector_t const& y)
  {
    return x-y;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t mul (vector_t const& x, vector_t const& y)
  {
    return x*y;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t div (vector_t const& x, vector_t const& y)
  {
    return x/y;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t madd (vector_t const& x, vector_t const& y, vector_t const& z)
  {
    return x*y+z;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t sqrt (vector_t const& x)
  {
    return std::sqrt(x);
  }
};

template<>
struct vecprops<CCTK_REAL4> {
  typedef CCTK_REAL4     scalar_t;
  typedef CCTK_REAL4_VEC vector_t;
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  int size()
  {
    return CCTK_REAL4_VEC_SIZE;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t load (scalar_t const& a)
  {
    return vec4_load(a);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
   vector_t loadu (scalar_t const& a)
  {
    return vec4_loadu(a);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t set1 (scalar_t const& a)
  {
    return vec4_set1(a);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  scalar_t elt (vector_t const& x, int const d)
  {
    return vec4_elt(x,d);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t neg (vector_t const& x)
  {
    return k4neg(x);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t add (vector_t const& x, vector_t const& y)
  {
    return k4add(x,y);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t sub (vector_t const& x, vector_t const& y)
  {
    return k4sub(x,y);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t mul (vector_t const& x, vector_t const& y)
  {
    return k4mul(x,y);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t div (vector_t const& x, vector_t const& y)
  {
    return k4div(x,y);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t madd (vector_t const& x, vector_t const& y, vector_t const& z)
  {
    return k4madd(x,y,z);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t sqrt (vector_t const& x)
  {
    return k4sqrt(x);
  }
};

template<>
struct vecprops<CCTK_REAL8> {
  typedef CCTK_REAL8     scalar_t;
  typedef CCTK_REAL8_VEC vector_t;
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  int size()
  {
    return CCTK_REAL8_VEC_SIZE;
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t load (scalar_t const& a)
  {
    return vec8_load(a);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t loadu (scalar_t const& a)
  {
    return vec8_loadu(a);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t set1 (scalar_t const& a)
  {
    return vec8_set1(a);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  scalar_t elt (vector_t const& x, int const d)
  {
    return vec8_elt(x,d);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t neg (vector_t const& x)
  {
    return k8neg(x);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t add (vector_t const& x, vector_t const& y)
  {
    return k8add(x,y);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t sub (vector_t const& x, vector_t const& y)
  {
    return k8sub(x,y);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t mul (vector_t const& x, vector_t const& y)
  {
    return k8mul(x,y);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t div (vector_t const& x, vector_t const& y)
  {
    return k8div(x,y);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t madd (vector_t const& x, vector_t const& y, vector_t const& z)
  {
    return k8madd(x,y,z);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vector_t sqrt (vector_t const& x)
  {
    return k8sqrt(x);
  }
};

template<typename T>
class vectype {
  typedef vecprops<T> props;
public:
  typedef typename props::vector_t vector_t;
  typedef typename props::scalar_t scalar_t;
  vector_t v;
  vectype() { }
  vectype(vectype const& x): v(x.v) { }
  vectype(vector_t const& x): v(x) { }
  explicit vectype(scalar_t const& a): v(props::set1(a)) { }
  operator vector_t() const { return v; }
  vectype& operator=(vectype const& x) { v=x.v; return *this; }
  
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  std::size_t size() const {
    return props::size();
  }
  
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype load(scalar_t const& a)
  {
    return props::load(a);
  }
  static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype loadu(scalar_t const& a)
  {
    return props::loadu(a);
  }
  
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  scalar_t elt(std::ptrdiff_t const d) const
  {
    return props::elt(*this, d);
  }
  
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype operator+() const
  {
    return *this;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype operator-() const
  {
    return props::neg(*this);
  }
  
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype operator+(vectype const& x) const
  {
    return props::add(*this, x);
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype operator-(vectype const& x) const
  {
    return props::sub(*this, x);
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype operator*(vectype const& x) const
  {
    return props::mul(*this, x);
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype operator/(vectype const& x) const
  {
    return props::div(*this, x);
  }
  
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype& operator+=(vectype const& x)
  {
    return *this = *this+x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype& operator-=(vectype const& x)
  {
    return *this = *this-x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype& operator*=(vectype const& x)
  {
    return *this = *this*x;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE
  vectype& operator/=(vectype const& x)
  {
    return *this = *this/x;
  }
};

template<typename T>
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
vectype<T> sqrt(vectype<T> const& x)
{
  return vecprops<T>::sqrt(x);
}

#endif



// Cache information

// Size of a a cache line in bytes
#ifndef CCTK_CACHELINE_SIZE
// TODO: Determine this properly
#  define CCTK_CACHELINE_SIZE 64
#endif

// Number of CCTK_REALs in a cache line
#define CCTK_REAL_CACHELINE_SIZE (CCTK_CACHELINE_SIZE / CCTK_REAL_PRECISION)
// If this fails, something is most likely wrong -- this would be a
// very weird (and inefficient?) architecture indeed
#if CCTK_REAL_CACHELINE_SIZE % CCTK_REAL_VEC_SIZE != 0
#  error "The cache line size is not a multiple of sizeof(CCTK_REAL_VEC)"
#endif



// For Kranc

#ifdef KRANC_C

#  undef KRANC_DIFF_FUNCTIONS
#  if ! VECTORISE_INLINE
#    define KRANC_DIFF_FUNCTIONS
#  endif

#  undef ToReal
#  define ToReal(x) (vec_set1(CCTK_REAL(x)))

#  undef IfThen
#  if (defined __PGI ||                                                 \
       defined _ARCH_450D ||                                            \
       (defined __ALTIVEC__ && defined _ARCH_PWR7))
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL_VEC vec_IfThen(CCTK_BOOLEAN x, CCTK_REAL_VEC y, CCTK_REAL_VEC z)
{
  if (x) return y; else return z;
}
#    define IfThen(x,y,z) vec_IfThen(x,y,z)
#  else
#    define IfThen(x,y,z) ((x) ? CCTK_REAL_VEC(y) : CCTK_REAL_VEC(z))
#  endif

#  undef KRANC_GFOFFSET3D
#  define KRANC_GFOFFSET3D(var,i,j,k)                                   \
  vec_loadu_maybe3((i),(j),(k),                                         \
                   *(CCTK_REAL const*)&                                 \
                   ((char const*)(var))[cdi*(i)+cdj*(j)+cdk*(k)])

#endif  // KRANC_C

#endif  // #ifndef VECTORS_H
