// -*-C++-*-
// Vectorise using IBM's Blue Gene/P Double Hummer (Power)

// Use the type double _Complex directly, without introducing a wrapper class
// Use macros instead of inline functions

// See <http://publib.boulder.ibm.com/infocenter/compbgpl/v9v111/index.jsp>



#include <assert.h>

#ifdef __cplusplus
#  include <builtins.h>
#endif



#define vec8_architecture "Double Hummer"

// Vector type corresponding to CCTK_REAL
#define CCTK_REAL8_VEC    double _Complex
#define CCTK_INTEGER8_VEC CCTK_REAL8_VEC
#define CCTK_BOOLEAN8_VEC CCTK_REAL8_VEC

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL8_VEC_SIZE 2

// Integer and boolean types corresponding to this real type
#define CCTK_INTEGER8     CCTK_REAL8
#define CCTK_BOOLEAN8     CCTK_REAL8



union k8const_t {
  double             f[2];
  unsigned long long i[2];
  CCTK_REAL8_VEC     vf;
};



// Create vectors, extract vector elements

#define vec8_set1(a)  (__cmplx(a,a))
#define vec8_set(a,b) (__cmplx(a,b))

#define vec8_elt0(x) (__creal(x))
#define vec8_elt1(x) (__cimag(x))
#define vec8_elt(x_,d)                          \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    CCTK_REAL8 a;                               \
    switch (d) {                                \
    case 0: a=vec8_elt0(x); break;              \
    case 1: a=vec8_elt1(x); break;              \
    }                                           \
    a;                                          \
  })
#define vec8_elti(x,d) vec8_elt(x,d)
#define vec8_eltb(x,d) vec8_elt(x,d)



// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
#define vec8_load(p)  (__lfpd((CCTK_REAL8 *)&(p)))
#if ! VECTORISE_ALWAYS_USE_ALIGNED_LOADS
#  define vec8_load_off1(p_)                    \
  ({                                            \
    CCTK_REAL8 const& p__=(p_);                 \
    CCTK_REAL8 const& p=p__;                    \
    vec8_set((&p)[0],(&p)[1]);                  \
  })
#else
#if 0
#  define vec8_load_off1(p_)                                    \
  ({                                                            \
    CCTK_REAL8 const& p__=(p_);                                 \
    CCTK_REAL8 const& p=p__;                                    \
    CCTK_REAL8_VEC const lo = __lfxd((CCTK_REAL8 *)(&p-1));     \
    CCTK_REAL8_VEC const hi = __lfxd((CCTK_REAL8 *)(&p+1));     \
    __fpsel(vec8_set(-1.0,+1.0),lo,hi);                         \
  })
#endif
#  define vec8_load_off1(p_)                            \
  ({                                                    \
    CCTK_REAL8 const& p__=(p_);                         \
    CCTK_REAL8 const& p=p__;                            \
    CCTK_REAL8_VEC const lo = vec8_load((&p)[-1]);      \
    CCTK_REAL8_VEC const hi = vec8_load((&p)[+1]);      \
    __fxmr(__fpsel(vec8_set(+1.0,-1.0),lo,hi));         \
  })
#endif
#define vec8_loadu(p_)                          \
  ({                                            \
    CCTK_REAL8 const& p__=(p_);                 \
    CCTK_REAL8 const& p=p__;                    \
    int const off = (ptrdiff_t)&p & 0xf;        \
    off==0 ? vec8_load(p) : vec8_load_off1(p);  \
  })

// Load a vector from memory that may or may not be aligned, as
// decided by the offset and the vector size
#if VECTORISE_ALWAYS_USE_UNALIGNED_LOADS
// Implementation: Always use unaligned load
#  define vec8_loadu_maybe(off,p)             vec8_loadu(p)
#  define vec8_loadu_maybe3(off1,off2,off3,p) vec8_loadu(p)
#else
#  define vec8_loadu_maybe(off,p_)              \
  ({                                            \
    CCTK_REAL8 const& p__=(p_);                 \
    CCTK_REAL8 const& p=p__;                    \
    (off) % CCTK_REAL8_VEC_SIZE == 0 ?          \
      vec8_load(p) :                            \
      vec8_load_off1(p);                        \
  })
#  if VECTORISE_ALIGNED_ARRAYS
// Assume all array x sizes are multiples of the vector size
#    define vec8_loadu_maybe3(off1,off2,off3,p) vec8_loadu_maybe(off1,p)
#  else
#    define vec8_loadu_maybe3(off1,off2,off3,p_)        \
  ({                                                    \
    CCTK_REAL8 const& p__=(p_);                         \
    CCTK_REAL8 const& p=p__;                            \
    ((off2) % CCTK_REAL8_VEC_SIZE != 0 or               \
     (off3) % CCTK_REAL8_VEC_SIZE != 0) ?               \
      vec8_loadu(p) :                                   \
      vec8_loadu_maybe(off1,p);                         \
  })
#  endif
#endif

// Store a vector to memory (aligned and non-temporal); this stores to
// a reference to a scalar
#define vec8_store(p,x)     (__stfpd(&(p),x))
#define vec8_storeu(p,x)    (__stfpd(&(p),x)) // this may not work
#define vec8_store_nta(p,x) (__stfpd(&(p),x)) // this doesn't avoid the cache

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare(i,imin,imax)                 \
  bool const v8stp_lo = (i)>=(imin);                            \
  bool const v8stp_hi = (i)+CCTK_REAL_VEC_SIZE-1<(imax)
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
#define vec8_store_nta_partial_lo(p,x,n) ((&(p))[0]=vec8_elt0(x))
#define vec8_store_nta_partial_hi(p,x,n) ((&(p))[1]=vec8_elt1(x))
#define vec8_store_nta_partial_mid(p,x,nlo,nhi) assert(0)



// Functions and operators

// Operators
#define k8neg(x) (__fpneg(x))

#define k8add(x,y) (__fpadd(x,y))
#define k8sub(x,y) (__fpsub(x,y))
#define k8mul(x,y) (__fpmul(x,y))
// Estimate for reciprocal
#define k8inv_init(x) (__fpre(x))
// One Newton iteration for reciprocal
#define k8inv_iter(x_,r_)                               \
  ({                                                    \
    CCTK_REAL8_VEC const x__=(x_);                      \
    CCTK_REAL8_VEC const r__=(r_);                      \
    CCTK_REAL8_VEC const x=x__;                         \
    CCTK_REAL8_VEC const r=r__;                         \
    /* r + r * (1 - x*r) */                             \
    k8madd(r, k8nmsub(x, r, vec8_set1(1.0)), r);        \
  })
// Reciprocal: First estimate, then apply two Newton iterations
#define k8inv(x_)                               \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    CCTK_REAL8_VEC const r0 = k8inv_init(x);    \
    CCTK_REAL8_VEC const r1 = k8inv_iter(x,r0); \
    CCTK_REAL8_VEC const r2 = k8inv_iter(x,r1); \
    r2;                                         \
  })
#define k8div(x,y) (__fpmul(x,k8inv(y)))

// Fused multiply-add, defined as [+-]x*y[+-]z
#define k8madd(x,y,z)  (__fpmadd(z,x,y))
#define k8msub(x,y,z)  (__fpmsub(z,x,y))
#define k8nmadd(x,y,z) (__fpnmadd(z,x,y))
#define k8nmsub(x,y,z) (__fpnmsub(z,x,y))

// Cheap functions
// TODO: handle -0 correctly
#define k8copysign(x,y) (k8ifthen(y,k8fabs(x),k8fnabs(x)))
#define k8fabs(x)       (__fpabs(x))
#define k8fmax(x_,y_)                           \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const y__=(y_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    CCTK_REAL8_VEC const y=y__;                 \
    __fpsel(k8sub(y,x),x,y);                    \
  })
#define k8fmin(x_,y_)                           \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const y__=(y_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    CCTK_REAL8_VEC const y=y__;                 \
    __fpsel(k8sub(x,y),x,y);                    \
  })
#define k8fnabs(x) (__fpnabs(x))
static const k8const_t k8zero = {{  0.0,  0.0, }};
static const k8const_t k8one  = {{ +1.0, +1.0, }};
static const k8const_t k8mone = {{ -1.0, -1.0, }};
#define k8sgn(x_)                                                       \
  ({                                                                    \
    CCTK_REAL_VEC x__=(x_);                                             \
    CCTK_REAL_VEC x=x__;                                                \
    CCTK_REAL_VEC iszero = k8fnabs(x);                                  \
    CCTK_REAL_VEC signedone = k8ifthen(x, k8one.vf, k8mone.vf);         \
    k8ifthen(iszero, k8zero.vf, signedone);                             \
  })
// Estimate for reciprocal square root
#define k8rsqrt_init(x) (__fprsqrte(x))
// One Newton iteration for reciprocal square root
#define k8rsqrt_iter(x_,rs_)                                    \
  ({                                                            \
    CCTK_REAL8_VEC const x__=(x_);                              \
    CCTK_REAL8_VEC const rs__=(rs_);                            \
    CCTK_REAL8_VEC const x=x__;                                 \
    CCTK_REAL8_VEC const rs=rs__;                               \
    /* rs (3/2 - x/2 rs^2) */                                   \
    k8mul(rs, k8msub(vec8_set1(1.5), x2, k8mul(rs, rs)));       \
  })
// Reciprocal square root: First estimate, then apply two Newton iterations
#define k8rsqrt(x_)                                     \
  ({                                                    \
    CCTK_REAL8_VEC const x__=(x_);                      \
    CCTK_REAL8_VEC const x=x__;                         \
    CCTK_REAL8_VEC const rs0 = k8rsqrt_init(x);         \
    CCTK_REAL8_VEC const x2  = k8mul(vec8_set1(0.5),x); \
    CCTK_REAL8_VEC const rs1 = k8rsqrt_iter(x2,rs0);    \
    CCTK_REAL8_VEC const rs2 = k8rsqrt_iter(x2,rs1);    \
    rs2;                                                \
  })
#define k8sqrt(x_)                              \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    k8mul(x, k8rsqrt(x));                       \
  })

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

// canonical true is +1.0, canonical false is -1.0
// >=0 is true, -0 is true (?), nan is false (?)
static const k8const_t k8lfalse_ = {{ -1.0, -1.0, }};
static const k8const_t k8ltrue_  = {{ +1.0, +1.0, }};
#define k8lfalse        (k8lfalse_.vf)
#define k8ltrue         (k8ltrue_.vf)
#define k8lnot(x)       (k8ifthen(x,k8lfalse,k8ltrue))
#define k8land(x,y)     (k8ifthen(x,y,k8lfalse))
#define k8lor(x,y)      (k8ifthen(x,k8ltrue,y))
#define k8lxor(x,y)     (k8ifthen(x,k8lnot(y),y))
#define k8ifthen(x,y,z) (__fpsel(x,z,y))

#define k8cmpeq(x,y) (__fpnabs((x)-(y)))
#define k8cmpne(x,y) (k8lnot(__fpnabs((x)-(y))))
#define k8cmpgt(x,y) (k8lnot(k8cmple(x,y)))
#define k8cmpge(x,y) ((x)-(y))
#define k8cmplt(x,y) (k8lnot(k8cmpge(x,y)))
#define k8cmple(x,y) ((y)-(x))
