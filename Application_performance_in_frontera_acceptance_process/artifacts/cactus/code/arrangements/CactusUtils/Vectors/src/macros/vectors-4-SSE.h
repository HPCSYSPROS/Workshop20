// Vectorise using Intel's or AMD's SSE

// Use the type __m128 directly, without introducing a wrapper class
// Use macros instead of inline functions



#include <assert.h>
#include <math.h>

#include <xmmintrin.h>
#ifdef __SSE4_1__
// Intel's SSE 4.1
#  include <smmintrin.h>
#endif
#ifdef __SSE4A__
// AMD's SSE 4a
#  include <ammintrin.h>
#endif
#ifdef __FMA4__
#  include <fma4intrin.h>
#endif



#ifdef __SSE4_1__
#  define vec4_architecture_SSE4_1 "+SSE4.1"
#else
#  define vec4_architecture_SSE4_1 ""
#endif
#ifdef __SSE4A__
#  define vec4_architecture_SSE4a "+SSE4A"
#else
#  define vec4_architecture_SSE4a ""
#endif
#ifdef __FMA4__
#  define vec4_architecture_FMA4 "+FMA4"
#else
#  define vec4_architecture_FMA4 ""
#endif
#define vec4_architecture "SSE" vec4_architecture_SSE4_1 vec4_architecture_SSE4a vec4_architecture_FMA4 " (32-bit precision)"



// Vector type corresponding to CCTK_REAL
#define CCTK_REAL4_VEC __m128

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL4_VEC_SIZE 4

// Integer and boolean types corresponding to this real type
#define CCTK_INTEGER4     CCTK_REAL4
#define CCTK_BOOLEAN4     CCTK_REAL4
#define CCTK_INTEGER4_VEC CCTK_REAL4_VEC
#define CCTK_BOOLEAN4_VEC CCTK_REAL4_VEC



union k4const_t {
  unsigned i[4];
  float    f[4];
  __m128i  vi;
  __m128   vf;
};

#define K4_ZERO 0x00000000UL
#define K4_IMIN 0x80000000UL
#define K4_IMAX 0x7fffffffUL



// Create vectors, extract vector elements

#define vec4_set1(a)      (_mm_set1_ps(a))
#define vec4_set(a,b,c,d) (_mm_set_ps(d,c,b,a)) // note reversed arguments

// original order is 0123
#define vec4_swap1032(x_)                       \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    _mm_shuffle_ps(x,x, _MM_SHUFFLE(2,3,0,1));  \
  })
#define vec4_swap2301(x_)                       \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    _mm_shuffle_ps(x,x, _MM_SHUFFLE(1,0,3,2));  \
  })
#define vec4_swap3210(x_)                       \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    _mm_shuffle_ps(x,x, _MM_SHUFFLE(0,1,2,3));  \
  })

#if defined(__PGI)
// _mm_cvtss_f32 does not exist on PGI compilers
#  define vec4_elt0(x)                          \
  ({                                            \
    CCTK_REAL4 a;                               \
    asm ("" : "=x" (a) : "0" (x));              \
    a;                                          \
  })
#else
#  define vec4_elt0(x) (_mm_cvtss_f32(x)) // this is a no-op
#endif
#define vec4_elt1(x) vec4_elt0(vec4_swap1032(x))
#define vec4_elt2(x) vec4_elt0(vec4_swap2301(x))
#define vec4_elt3(x) vec4_elt0(vec4_swap3210(x))
#if defined(__PGI)
#  define vec4_elt(x_,d)                        \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    CCTK_REAL4 a;                               \
    if (d==0)      a=vec4_elt0(x);              \
    else if (d==1) a=vec4_elt1(x);              \
    else if (d==2) a=vec4_elt2(x);              \
    else if (d==3) a=vec4_elt3(x);              \
    a;                                          \
  })
#else
#  define vec4_elt(x_,d)                        \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    CCTK_REAL4 a;                               \
    switch (d) {                                \
    case 0: a=vec4_elt0(x); break;              \
    case 1: a=vec4_elt1(x); break;              \
    case 2: a=vec4_elt2(x); break;              \
    case 3: a=vec4_elt3(x); break;              \
    }                                           \
    a;                                          \
  })
#endif



// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
#define vec4_load(p)  (_mm_load_ps(&(p)))
#define vec4_loadu(p) (_mm_loadu_ps(&(p)))
#if ! VECTORISE_ALWAYS_USE_ALIGNED_LOADS
#  define vec4_load_off1(p) vec_loadu(p)
#  define vec4_load_off2(p) vec_loadu(p)
#  define vec4_load_off3(p) vec_loadu(p)
#else
#  define vec4_load_off1(p_)                                            \
  ({                                                                    \
    CCTK_REAL4 const& p__=(p_);                                         \
    CCTK_REAL4 const& p=p__;                                            \
    CCTK_REAL4_VEC const lo=vec4_load((&p)[-1]);                        \
    CCTK_REAL4_VEC const hi=vec4_load((&p)[+3]);                        \
    assert(0);                                                          \
    CCTK_REAL4_VEC const hi2=_mm_shuffle_ps(lo,hi, _MM_SHUFFLE(0,1,2,3)); \
    _mm_shuffle_ps(lo,hi2, _MM_SHUFFLE(2,1,3,0));                       \
  })
#  define vec4_load_off2(p_)                            \
  ({                                                    \
    CCTK_REAL4 const& p__=(p_);                         \
    CCTK_REAL4 const& p=p__;                            \
    CCTK_REAL4_VEC const lo=vec4_load((&p)[-2]);        \
    CCTK_REAL4_VEC const hi=vec4_load((&p)[+2]);        \
    _mm_shuffle_ps(lo,hi, _MM_SHUFFLE(1,0,3,2));        \
  })
#  define vec4_load_off3(p_)                                            \
  ({                                                                    \
    CCTK_REAL4 const& p__=(p_);                                         \
    CCTK_REAL4 const& p=p__;                                            \
    CCTK_REAL4_VEC const lo=vec4_load((&p)[-1]);                        \
    CCTK_REAL4_VEC const hi=vec4_load((&p)[+3]);                        \
    assert(0);                                                          \
    CCTK_REAL4_VEC const lo2=_mm_shuffle_ps(lo,hi, _MM_SHUFFLE(0,1,2,3)); \
    _mm_shuffle_ps(lo2,hi, _MM_SHUFFLE(3,0,2,1));                       \
  })
#endif

// Load a vector from memory that may or may not be aligned, as
// decided by the offset off and the vector size
#if VECTORISE_ALWAYS_USE_UNALIGNED_LOADS
// Implementation: Always use unaligned load
#  define vec4_loadu_maybe(off,p)             vec4_loadu(p)
#  define vec4_loadu_maybe3(off1,off2,off3,p) vec4_loadu(p)
#else
#  define vec4_loadu_maybe(off,p_)              \
  ({                                            \
    CCTK_REAL4 const& p__=(p_);                 \
    CCTK_REAL4 const& p=p__;                    \
    (off) % CCTK_REAL4_VEC_SIZE == 0 ?          \
      vec4_load(p) :                            \
      vec4_loadu(p);                            \
  })
#  if VECTORISE_ALIGNED_ARRAYS
// Assume all array x sizes are multiples of the vector size
#    define vec4_loadu_maybe3(off1,off2,off3,p) \
  vec4_loadu_maybe(off1,p)
#  else
#    define vec4_loadu_maybe3(off1,off2,off3,p) \
  vec4_loadu_maybe((off1)|(off2)|(off3),p)
#  endif
#endif

// Store a vector to memory (aligned and non-temporal); this stores to
// a reference to a scalar
#define vec4_store(p,x)  (_mm_store_ps(&(p),x))
#define vec4_storeu(p,x) (_mm_storeu_ps(&(p),x))
#if ! VECTORISE_STREAMING_STORES
#  define vec4_store_nta(p,x) vec4_store(p,x)
#else
#  define vec4_store_nta(p,x) (_mm_stream_ps(&(p),x))
#endif

// Store a partial vector (aligned and non-temporal)
#define vec4_store_partial_prepare(i,imin,imax)                         \
  int v4stp_lo_skip = (imin)-(i);                                       \
  int v4stp_hi_skip = (i)+CCTK_REAL_VEC_SIZE-(imax);                    \
  if (CCTK_BUILTIN_EXPECT(v4stp_lo_skip < 0, true)) v4stp_lo_skip = 0;  \
  if (CCTK_BUILTIN_EXPECT(v4stp_hi_skip < 0, true)) v4stp_hi_skip = 0;
// Ignoring VECTORISE_STREAMING_STORES for partial stores
#define vec4_store_nta_partial(p_,x_)                                   \
  ({                                                                    \
    CCTK_REAL4& p__=(p_);                                               \
    CCTK_REAL4& p=p__;                                                  \
    CCTK_REAL4_VEC const x__=(x_);                                      \
    CCTK_REAL4_VEC const x=x__;                                         \
    if (CCTK_BUILTIN_EXPECT(v4stp_lo_skip==0 and v4stp_hi_skip==0, true)) { \
      vec4_store_nta(p,x);                                              \
    } else {                                                            \
      /* these cases fall through */                                    \
      switch (v4stp_lo_skip) {                                          \
      case 0:                                                           \
        (&p)[0] = vec4_elt0(x);                                         \
      case 1:                                                           \
        if (v4stp_hi_skip>=3) break;                                    \
        (&p)[1] = vec4_elt1(x);                                         \
      case 2:                                                           \
        if (v4stp_hi_skip>=2) break;                                    \
        (&p)[2] = vec4_elt2(x);                                         \
      case 3:                                                           \
        if (v4stp_hi_skip>=1) break;                                    \
        (&p)[3] = vec4_elt3(x);                                         \
      }                                                                 \
    }                                                                   \
  })

// Ignoring VECTORISE_STREAMING_STORES for partial stores
#define vec4_store_nta_partial_lo(p_,x_,n)      \
  ({                                            \
    CCTK_REAL4          & p__=(p_);             \
    CCTK_REAL4_VEC const  x__=(x_);             \
    CCTK_REAL4          & p=p__;                \
    CCTK_REAL4_VEC const  x=x__;                \
    /* these cases fall through */              \
    switch (n) {                                \
    case 3: (&p)[2] = vec4_elt2(x);             \
    case 2: (&p)[1] = vec4_elt1(x);             \
    case 1: (&p)[0] = vec4_elt0(x);             \
    }                                           \
  })
#define vec4_store_nta_partial_hi(p_,x_,n)      \
  ({                                            \
    CCTK_REAL4          & p__=(p_);             \
    CCTK_REAL4_VEC const  x__=(x_);             \
    CCTK_REAL4          & p=p__;                \
    CCTK_REAL4_VEC const  x=x__;                \
    /* these cases fall through */              \
    switch (n) {                                \
    case 3: (&p)[1]=vec4_elt1(x);               \
    case 2: (&p)[2]=vec4_elt2(x);               \
    case 1: (&p)[3]=vec4_elt3(x);               \
    }                                           \
  })
#define vec4_store_nta_partial_mid(p_,x_,nlo,nhi)               \
  ({                                                            \
    CCTK_REAL4          & p__=(p_);                             \
    CCTK_REAL4_VEC const  x__=(x_);                             \
    CCTK_REAL4          & p=p__;                                \
    CCTK_REAL4_VEC const  x=x__;                                \
    /* these cases fall through */                              \
    switch (nhi) {                                              \
    case 3: if (nlo<2) break; (&p)[1] = vec4_elt1(x);           \
    case 2: if (nlo<3) break; (&p)[2] = vec4_elt2(x);           \
    }                                                           \
  })



// Functions and operators

static const k4const_t k4sign_mask = {{ K4_IMIN, K4_IMIN, K4_IMIN, K4_IMIN, }};

// Operators
#define k4neg(x) (_mm_xor_ps(k4sign_mask.vf,x))
// #define k4inv(x)
// TODO: provide k4inv via rcp and Newton-Raphson
// This is described in AMD's publication 47414.
// This should apply for AVX as well.

#define k4add(x,y) (_mm_add_ps(x,y))
#define k4sub(x,y) (_mm_sub_ps(x,y))
#define k4mul(x,y) (_mm_mul_ps(x,y))
// TODO: use k4inv and k4mul instead
#define k4div(x,y) (_mm_div_ps(x,y))

// Fused multiply-add, defined as [+-]x*y[+-]z
#ifdef __FMA4__
#  define k4madd(x,y,z)  (_mm_macc_ps(x,y,z))
#  define k4msub(x,y,z)  (_mm_msub_ps(x,y,z))
#  define k4nmadd(x,y,z) (_mm_nmsub_ps(x,y,z))
#  define k4nmsub(x,y,z) (_mm_nmacc_ps(x,y,z))
#else
#  define k4madd(x,y,z)  (k4add(k4mul(x,y),z))
#  define k4msub(x,y,z)  (k4sub(k4mul(x,y),z))
#  define k4nmadd(x,y,z) (k4sub(k4neg(z),k4mul(x,y)))
#  define k4nmsub(x,y,z) (k4sub(z,k4mul(x,y)))
#endif

// Cheap functions
#define k4copysign(x,y)                                         \
  (_mm_or_ps(_mm_andnot_ps(k4sign_mask.vf,x),                   \
             _mm_and_ps(k4sign_mask.vf,y)))
#define k4fabs(x)   (_mm_andnot_ps(k4sign_mask.vf,x))
#define k4fmax(x,y) (_mm_max_ps(x,y))
#define k4fmin(x,y) (_mm_min_ps(x,y))
#define k4fnabs(x)  (_mm_or_ps(k4sign_mask.vf,x))
#define k4sgn(x_)                                                       \
  ({                                                                    \
    CCTK_REAL_VEC const x__=(x_);                                       \
    CCTK_REAL_VEC const x=x__;                                          \
    CCTK_REAL_VEC const iszero = _mm_cmpeq_ps(vec4_set1(0.0f), x);      \
    CCTK_REAL_VEC const sign = _mm_and_ps(k4sign_mask.vf, x);           \
    CCTK_REAL_VEC const signedone = _mm_or_ps(vec4_set1(1.0f), sign);   \
    k4ifthen(iszero, vec4_set1(0.0f), signedone);                       \
  })
// TODO: maybe use rsqrt and Newton-Raphson
#define k4sqrt(x)   (_mm_sqrt_ps(x))

// Expensive functions
#define K4REPL(f,x_)                            \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    vec4_set(f(vec4_elt0(x)),                   \
             f(vec4_elt1(x)),                   \
             f(vec4_elt2(x)),                   \
             f(vec4_elt3(x)));                  \
  })
#define K4REPL2S(f,x_,a_)                       \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4     const a__=(a_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    CCTK_REAL4     const a=a__;                 \
    vec4_set(f(vec4_elt0(x),a),                 \
             f(vec4_elt1(x),a),                 \
             f(vec4_elt2(x),a),                 \
             f(vec4_elt3(x),a));                \
  })
#define K4REPL2(f,x_,y_)                        \
  ({                                            \
    CCTK_REAL4_VEC const x__=(x_);              \
    CCTK_REAL4_VEC const y__=(y_);              \
    CCTK_REAL4_VEC const x=x__;                 \
    CCTK_REAL4_VEC const y=y__;                 \
    vec4_set(f(vec4_elt0(x),vec4_elt0(y)),      \
             f(vec4_elt1(x),vec4_elt1(y)),      \
             f(vec4_elt2(x),vec4_elt2(y)),      \
             f(vec4_elt3(x),vec4_elt3(y)));     \
  })

#define k4acos(x)    K4REPL(acosf,x)
#define k4acosh(x)   K4REPL(acoshf,x)
#define k4asin(x)    K4REPL(asinf,x)
#define k4asinh(x)   K4REPL(asinhf,x)
#define k4atan(x)    K4REPL(atanf,x)
#define k4atan2(x,y) K4REPL2(atan2f,x,y)
#define k4atanh(x)   K4REPL(atanhf,x)
#define k4cos(x)     K4REPL(cosf,x)
#define k4cosh(x)    K4REPL(coshf,x)
#define k4exp(x)     K4REPL(expf,x)
#define k4log(x)     K4REPL(logf,x)
#define k4pow(x,a)   K4REPL2S(powf,x,a)
#define k4sin(x)     K4REPL(sinf,x)
#define k4sinh(x)    K4REPL(sinhf,x)
#define k4tan(x)     K4REPL(tanf,x)
#define k4tanh(x)    K4REPL(tanhf,x)

static const k4const_t k4lfalse_ = {{  0U,  0U,  0U,  0U, }};
static const k4const_t k4ltrue_  = {{ ~0U, ~0U, ~0U, ~0U, }};
#define k4lfalse (k4lfalse_.vf)
#define k4ltrue  (k4ltrue_.vf)
#define k4lnot(x)   (_mm_xor_ps(k4ltrue,x))
#define k4land(x,y) (_mm_and_ps(x,y))
#define k4lor(x,y)  (_mm_or_ps(x,y))
#define k4lxor(x,y) (_mm_xor_ps(x,y))

#ifdef __SSE4_1__
#  define k4ifthen(x,y,z) (_mm_blendv_ps(z,y,x))
#elif 0
#  ifdef __cplusplus
#    define k4signbit(x) ({ using namespace std; signbit(x); })
#  else
#    define k4signbit(x) (signbitf(x))
#  endif
#  define k4ifthen(x,y,z)                                               \
  ({                                                                    \
    CCTK_REAL4_VEC const x__=(x_);                                      \
    CCTK_REAL4_VEC const y__=(y_);                                      \
    CCTK_REAL4_VEC const z__=(z_);                                      \
    CCTK_REAL4_VEC const x=x__;                                         \
    CCTK_REAL4_VEC const y=y__;                                         \
    CCTK_REAL4_VEC const z=z__;                                         \
    vec4_set(k4signbit(vec4_elt0(x)) ? vec4_elt0(y) : vec4_elt0(z),     \
             k4signbit(vec4_elt1(x)) ? vec4_elt1(y) : vec4_elt1(z),     \
             k4signbit(vec4_elt2(x)) ? vec4_elt2(y) : vec4_elt2(z),     \
             k4signbit(vec4_elt3(x)) ? vec4_elt3(y) : vec4_elt3(z));    \
  })
#elif 0
// We don't need to shift -- the condition (mask) will be either all
// zeros or all ones
#  define k4ifthen(x_,y_,z_)                                    \
  ({                                                            \
    CCTK_REAL4_VEC const x__=(x_);                              \
    CCTK_REAL4_VEC const y__=(y_);                              \
    CCTK_REAL4_VEC const z__=(z_);                              \
    CCTK_REAL4_VEC const x=x__;                                 \
    CCTK_REAL4_VEC const y=y__;                                 \
    CCTK_REAL4_VEC const z=z__;                                 \
    CCTK_REAL4_VEC const mask =                                 \
      (__m128)_mm_srai_epi32((__m128i)x, 31);                   \
    /* (z & ~mask) | (y & mask) */                              \
    _mm_or_ps(_mm_andnot_ps(mask, z), _mm_and_ps(mask, y));     \
  })
#else
#  define k4ifthen(x_,y_,z_)                                    \
  ({                                                            \
    CCTK_REAL4_VEC const x__=(x_);                              \
    CCTK_REAL4_VEC const y__=(y_);                              \
    CCTK_REAL4_VEC const z__=(z_);                              \
    CCTK_REAL4_VEC const x=x__;                                 \
    CCTK_REAL4_VEC const y=y__;                                 \
    CCTK_REAL4_VEC const z=z__;                                 \
    /* (z & ~mask) | (y & mask)   where imask = ~mask */        \
    _mm_or_ps(_mm_and_ps(x, y), _mm_andnot_ps(x, z));           \
  })
#endif

#define k4cmpeq(x,y) (_mm_cmpeq_ps(x,y))
#define k4cmpne(x,y) (_mm_cmpneq_ps(x,y))
#define k4cmpgt(x,y) (_mm_cmpgt_ps(x,y))
#define k4cmpge(x,y) (_mm_cmpge_ps(x,y))
#define k4cmplt(x,y) (_mm_cmplt_ps(x,y))
#define k4cmple(x,y) (_mm_cmple_ps(x,y))
