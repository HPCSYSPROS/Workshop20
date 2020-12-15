// -*-C++-*-
// Vectorise using IBM's Altivec VSX (Power)

// Use the type vector double directly, without introducing a wrapper class
// Use macros instead of inline functions

// See <http://pic.dhe.ibm.com/infocenter/comphelp/v111v131/index.jsp>



#include <altivec.h>
#include <math.h>



#define vec8_architecture "VSX"

// Vector type corresponding to CCTK_REAL
typedef vector double           CCTK_REAL8_VEC;
typedef vector signed long long CCTK_INTEGER8_VEC;
typedef vector bool long long   CCTK_BOOLEAN8_VEC;

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL8_VEC_SIZE 2

vec_static_assert(sizeof(CCTK_REAL8_VEC) ==
                  sizeof(CCTK_REAL8) * CCTK_REAL8_VEC_SIZE);

// Integer and boolean types corresponding to this real type
typedef long long          CCTK_INTEGER8;
typedef unsigned long long CCTK_BOOLEAN8;



// Create vectors, extract vector elements

#define vec8_set1(a)  (vec_splats(a))
#define vec8_set1i(a) (vec_splats(a))
#define vec8_set(a,b)                           \
  ({                                            \
    CCTK_REAL8_VEC x;                           \
    x[0]=(a);                                   \
    x[1]=(b);                                   \
    x;                                          \
  })

#define vec8_elt0(x)   ((x)[0])
#define vec8_elt1(x)   ((x)[1])
// #define vec8_elt(x,d)  ((x)[d])
// #define vec8_elti(x,d) ((x)[d])
// #define vec8_eltb(x,d) ((x)[d])
static inline CCTK_REAL8 vec8_elt(CCTK_REAL8_VEC x, int d) { return x[d]; }
static inline CCTK_INTEGER8 vec8_elti(CCTK_INTEGER8_VEC x, int d) { return x[d]; }
static inline CCTK_BOOLEAN8 vec8_eltb(CCTK_BOOLEAN8_VEC x,int d) { return x[d]; }



// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
#define vec8_load(p)  (*(CCTK_REAL8_VEC const*)&(p))
#define vec8_loadu(p) (*(CCTK_REAL8_VEC const*)&(p))

// Load a vector from memory that may or may not be aligned, as
// decided by the offset and the vector size
#define vec8_loadu_maybe(off,p)             (vec8_loadu(p))
#define vec8_loadu_maybe3(off1,off2,off3,p) (vec8_loadu(p))

// Store a vector to memory (aligned and non-temporal); this stores to
// a reference to a scalar
#define vec8_store(p,x)  (*(CCTK_REAL8_VEC*)&(p)=(x))
#define vec8_storeu(p,x) (*(CCTK_REAL8_VEC*)&(p)=(x))
// stvxl instruction doesn't exist for double precision
#define vec8_store_nta(p,x) vec8_store(p,x)

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare(i,imin,imax)                 \
  bool const v8stp_lo = (i)>=(imin);                            \
  bool const v8stp_hi = (i)+CCTK_REAL8_VEC_SIZE-1<(imax)
#define vec8_store_nta_partial(p_,x_)                           \
  ({                                                            \
    CCTK_REAL8& p__=(p_);                                       \
    CCTK_REAL8& p=p__;                                          \
    CCTK_REAL8_VEC const x__=(x_);                              \
    CCTK_REAL8_VEC const x=x__;                                 \
    if (CCTK_BUILTIN_EXPECT(v8stp_lo and v8stp_hi, true)) {     \
      vec8_store(p,x);                                          \
    } else if (v8stp_lo) {                                      \
      (&p)[0]=vec8_elt0(x);                                     \
    } else if (v8stp_hi) {                                      \
      (&p)[1]=vec8_elt1(x);                                     \
    }                                                           \
  })

// Store a lower or higher partial vector (aligned and non-temporal);
// the non-temporal hint is probably ignored
#define vec8_store_nta_partial_lo(p,x,n) ((&(p))[0]=(x)[0])
#define vec8_store_nta_partial_hi(p,x,n) ((&(p))[1]=(x)[1])
#define vec8_store_nta_partial_mid(p,x,nlo,nhi) (assert(0))



// Functions and operators

// Operators
#define k8neg(x) (-(x))

#define k8add(x,y) ((x)+(y))
#define k8sub(x,y) ((x)-(y))
#define k8mul(x,y) ((x)*(y))
#define k8div(x,y) ((x)/(y))

// Fused multiply-add, defined as [+-]x*y[+-]z
#define k8madd(x,y,z)  (vec_madd(x,y,z))
#define k8msub(x,y,z)  (vec_msub(x,y,z))
#define k8nmadd(x,y,z) (vec_nmadd(x,y,z))
#define k8nmsub(x,y,z) (vec_nmsub(x,y,z))

// Cheap functions
#define k8copysign(x,y) (vec_cpsgn(y,x))
#define k8fabs(x)       (vec_abs(x))
#define k8fmax(x,y)     (vec_max(x,y))
#define k8fmin(x,y)     (vec_min(x,y))
#define k8fnabs(x)      (vec_nabs(x))
#define k8sgn(x_)                                                       \
  ({                                                                    \
    CCTK_REAL8_VEC x__=(x_);                                            \
    CCTK_REAL8_VEC x=x__;                                               \
    CCTK_BOOLEAN8_VEC iszero = k8cmpeq(x,vec8_set1((CCTK_REAL8)0.0));   \
    CCTK_REAL8_VEC signedone = k8copysign(vec8_set1((CCTK_REAL8)1.0),x); \
    k8ifthen(iszero, vec8_set1((CCTK_REAL8)0.0), signedone);            \
  })
#define k8sqrt(x)       (vec_sqrt(x))

// Expensive functions
#define K8REPL(f,x_)                            \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    vec8_set(f(vec8_elt0(x)),                   \
             f(vec8_elt1(x)));                  \
  })
#define K8REPL2S(f,x_,a_)                       \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8     const a__=(a_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    CCTK_REAL8     const a=a__;                 \
    vec8_set(f(vec8_elt0(x),a),                 \
             f(vec8_elt1(x),a));                \
  })
#define K8REPL2(f,x_,y_)                        \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const y__=(y_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    CCTK_REAL8_VEC const y=y__;                 \
    vec8_set(f(vec8_elt0(x),vec8_elt0(y)),      \
             f(vec8_elt1(x),vec8_elt1(y)));     \
  })

#define k8acos(x)    K8REPL(acos,x)
#define k8acosh(x)   K8REPL(acosh,x)
#define k8asin(x)    K8REPL(asin,x)
#define k8asinh(x)   K8REPL(asinh,x)
#define k8atan(x)    K8REPL(atan,x)
#define k8atan2(x,y) K8REPL2(atan2,x,y)
#define k8atanh(x)   K8REPL(atanh,x)
#define k8cos(x)     K8REPL(cos,x)
#define k8cosh(x)    K8REPL(cosh,x)
#define k8exp(x)     K8REPL(exp,x)
#define k8fmod(x,y)  K8REPL2(fmod,x,y)
#define k8log(x)     K8REPL(log,x)
#define k8pow(x,a)   K8REPL2S(pow,x,a)
#define k8sin(x)     K8REPL(sin,x)
#define k8sinh(x)    K8REPL(sinh,x)
#define k8tan(x)     K8REPL(tan,x)
#define k8tanh(x)    K8REPL(tanh,x)

// canonical true is -1LL, canonical false is 0LL
// truth values are interpreted bit-wise
// #define k8lfalse        ({ CCTK_BOOLEAN8_VEC dummy; vec_xor(dummy,dummy); })
// #define k8ltrue         (k8lnot(k8lfalse))
static inline CCTK_BOOLEAN8_VEC k8lfalse1()
{
  CCTK_BOOLEAN8_VEC dummy;
  return vec_xor(dummy,dummy);
}
#define k8lfalse (k8lfalse1())
#define k8ltrue  (k8lnot(k8lfalse))

// #define k8lnot(x_)                              \
//   ({                                            \
//     CCTK_BOOLEAN8_VEC x__=(x_);                 \
//     CCTK_BOOLEAN8_VEC x=x__;                    \
//     vec_nor(x,x);                               \
//   })
static inline CCTK_BOOLEAN8_VEC k8lnot(CCTK_BOOLEAN8_VEC x)
{
  return vec_nor(x,x);
}

#define k8land(x,y)     (vec_and(x,y))
#define k8lor(x,y)      (vec_or(x,y))
#define k8lxor(x,y)     (vec_xor(x,y))
#define k8ifthen(x,y,z) (vec_sel(z,y,x))

#define k8cmpeq(x,y) (vec_cmpeq(x,y))
#define k8cmpne(x,y) (k8lnot(vec_cmpeq(x,y)))
#define k8cmpgt(x,y) (vec_cmpgt(x,y))
#define k8cmpge(x,y) (vec_cmpge(x,y))
#define k8cmplt(x,y) (vec_cmplt(x,y))
#define k8cmple(x,y) (vec_cmple(x,y))
