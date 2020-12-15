// Vectorise using Intel's or AMD's SSE2

// Use the type __m128d directly, without introducing a wrapper class
// Use macros instead of inline functions



#include <assert.h>
#include <math.h>

#include <emmintrin.h>
#ifdef __SSE4_1__
// Intel's SSE 4.1
#  include <smmintrin.h>
#endif
#ifdef __SSE4A__
// AMD's SSE 4a
#  include <ammintrin.h>

// Intel compilers don't support SSE 4a. Here is how we can implement
// these instructions in assembler instead:

// inline void __attribute__((__always_inline__))
//   _mm_stream_sd (double *p, __m128d x)
// {
//   asm ("movntsd %[x],%[p]" : "=m" (*p) : [p] "m" (*p), [x] "x" (x));
// }

#endif
#ifdef __FMA4__
#  include <fma4intrin.h>
#endif



#ifdef __SSE4_1__
#  define vec8_architecture_SSE4_1 "+SSE4.1"
#else
#  define vec8_architecture_SSE4_1 ""
#endif
#ifdef __SSE4A__
#  define vec8_architecture_SSE4a "+SSE4A"
#else
#  define vec8_architecture_SSE4a ""
#endif
#ifdef __FMA4__
#  define vec8_architecture_FMA4 "+FMA4"
#else
#  define vec8_architecture_FMA4 ""
#endif
#define vec8_architecture "SSE2" vec8_architecture_SSE4_1 vec8_architecture_SSE4a vec8_architecture_FMA4 " (64-bit precision)"



// Vector type corresponding to CCTK_REAL
#define CCTK_REAL8_VEC __m128d

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL8_VEC_SIZE 2

// Integer and boolean types corresponding to this real type
#define CCTK_INTEGER8     CCTK_REAL8
#define CCTK_BOOLEAN8     CCTK_REAL8
#define CCTK_INTEGER8_VEC CCTK_REAL8_VEC
#define CCTK_BOOLEAN8_VEC CCTK_REAL8_VEC



union k8const_t {
  long long i[2];
  double    f[2];
  __m128i   vi;
  __m128d   vf;
};

#define K8_IMIN ((long long)0x8000000000000000ULL)



// Create vectors, extract vector elements

#define vec8_set1(a)  (_mm_set1_pd(a))
#define vec8_set(a,b) (_mm_set_pd(b,a)) // note reversed arguments

// original order is 01
#define vec8_swap10(x_)                         \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    _mm_shuffle_pd(x,x, _MM_SHUFFLE2(0,1));     \
  })

#define vec8_elt0(x) (((CCTK_REAL8 const*)&(x))[0])
#define vec8_elt1(x) (((CCTK_REAL8 const*)&(x))[1])
#define vec8_elt(x,d) (((CCTK_REAL8 const*)&(x))[d])



// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
#define vec8_load(p)  (_mm_load_pd(&(p)))
#define vec8_loadu(p) (_mm_loadu_pd(&(p)))
#if ! VECTORISE_ALWAYS_USE_ALIGNED_LOADS
#  define vec8_load_off1(p) vec_loadu(p)
#else
#  define vec8_load_off1(p_)                                    \
  ({                                                            \
    CCTK_REAL8 const& p__=(p_);                                 \
    CCTK_REAL8 const& p=p__;                                    \
    _mm_shuffle_pd(vec8_load((&p)[-1]),                         \
                   vec8_load((&p)[+1]), _MM_SHUFFLE2(0,1));     \
  })
#endif

// Load a vector from memory that may or may not be aligned, as
// decided by the offset off and the vector size
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
#    define vec8_loadu_maybe3(off1,off2,off3,p) \
  vec8_loadu_maybe(off1,p)
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
#define vec8_store(p,x)  (_mm_store_pd(&(p),x))
#define vec8_storeu(p,x) (_mm_storeu_pd(&(p),x))
#if ! VECTORISE_STREAMING_STORES
#  define vec8_store_nta(p,x) vec8_store(p,x)
#else
#  define vec8_store_nta(p,x) (_mm_stream_pd(&(p),x))
#endif

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare(i,imin,imax)                 \
  bool const v8stp_lo = (i)>=(imin);                            \
  bool const v8stp_hi = (i)+CCTK_REAL_VEC_SIZE-1<(imax)
#if VECTORISE_STREAMING_STORES && defined(__SSE4A__)
#  define vec8_store_nta_partial(p_,x_)                         \
  ({                                                            \
    CCTK_REAL8& p__=(p_);                                       \
    CCTK_REAL8& p=p__;                                          \
    CCTK_REAL8_VEC const x__=(x_);                              \
    CCTK_REAL8_VEC const x=x__;                                 \
    if (CCTK_BUILTIN_EXPECT(v8stp_lo and v8stp_hi, true)) {     \
      vec8_store_nta(p,x);                                      \
    } else if (v8stp_lo) {                                      \
      _mm_stream_sd(&p,x);                                      \
    } else if (v8stp_hi) {                                      \
      _mm_stream_sd(&p+1, vec8_swap10(x));                      \
    }                                                           \
  })
#else
#  define vec8_store_nta_partial(p_,x_)                         \
  ({                                                            \
    CCTK_REAL8& p__=(p_);                                       \
    CCTK_REAL8& p=p__;                                          \
    CCTK_REAL8_VEC const x__=(x_);                              \
    CCTK_REAL8_VEC const x=x__;                                 \
    if (CCTK_BUILTIN_EXPECT(v8stp_lo and v8stp_hi, true)) {     \
      vec8_store_nta(p,x);                                      \
    } else if (v8stp_lo) {                                      \
      _mm_storel_pd(&p,x);                                      \
    } else if (v8stp_hi) {                                      \
      _mm_storeh_pd(&p+1,x);                                    \
    }                                                           \
  })
#endif

// Store a lower or higher partial vector (aligned and non-temporal)
#if ! VECTORISE_STREAMING_STORES
#  define vec8_store_nta_partial_lo(p,x,n) (_mm_storel_pd(&(p),x))
#  define vec8_store_nta_partial_hi(p,x,n) (_mm_storeh_pd(&(p)+1,x))
#else
#  if defined(__SSE4A__)
#    define vec8_store_nta_partial_lo(p,x,n) (_mm_stream_sd(&(p),x))
#    define vec8_store_nta_partial_hi(p,x,n)    \
  (_mm_stream_sd(&(p)+1, vec8_swap10(x)))
#  else
// TODO: use clflush once a whole cache line has been written (cache
// lines are usually larger than the CPU vector size)
#    define vec8_store_nta_partial_lo(p_,x,n)   \
  ({                                            \
    CCTK_REAL8& p__=(p_);                       \
    CCTK_REAL8& p=p__;                          \
    _mm_storel_pd(&p,x);                        \
    /* _mm_clflush(&p); */                      \
  })
#    define vec8_store_nta_partial_hi(p_,x,n)   \
  ({                                            \
    CCTK_REAL8& p__=(p_);                       \
    CCTK_REAL8& p=p__;                          \
    _mm_storeh_pd(&p+1,x);                      \
    /* _mm_clflush(&p+1); */                    \
  })
#  endif
#endif
#if 0
// This is slower; we would need a non-temporal read
#define vec8_store_nta_partial_lo(p,x,n)        \
  vec8_store_nta(p, _mm_loadh_pd(x,&(p)+1))
#define vec8_store_nta_partial_hi(p,x,n)        \
  vec8_store_nta(p, _mm_loadl_pd(x,&(p)))
#endif
#define vec8_store_nta_partial_mid(p,x,nlo,nhi) assert(0)



// Functions and operators

static const k8const_t k8sign_mask = {{ K8_IMIN, K8_IMIN, }};

// Operators

// #define k8inot(x) (_mm_xor_si128(k8all_mask,x))
// 
// #define k8iand(x,y) (_mm_and_si128(x,y))
// #define k8ior(x,y)  (_mm_or_si128(x,y))
// #define k8ixor(x,y) (_mm_xor_si128(x,y))
// 
// #define k8ineg(x) (_mm_xor_pd(k8sign_mask,x))
// 
// #define k8iadd(x,y) (_mm_add_epi64(x,y))
// #define k8isub(x,y) (_mm_sub_epi64(x,y))
// 
// #define k8not(x) (_mm_xor_pd(k8all_mask,x))
// 
// #define k8and(x,y) (_mm_and_pd(x,y))
// #define k8or(x,y)  (_mm_or_pd(x,y))
// #define k8xor(x,y) (_mm_xor_pd(x,y))

#define k8neg(x) (_mm_xor_pd(k8sign_mask.vf,x))

#define k8add(x,y) (_mm_add_pd(x,y))
#define k8sub(x,y) (_mm_sub_pd(x,y))
#define k8mul(x,y) (_mm_mul_pd(x,y))
#define k8div(x,y) (_mm_div_pd(x,y))

// Fused multiply-add, defined as [+-]x*y[+-]z
#ifdef __FMA4__
#  define k8madd(x,y,z)  (_mm_macc_pd(x,y,z))
#  define k8msub(x,y,z)  (_mm_msub_pd(x,y,z))
#  define k8nmadd(x,y,z) (_mm_nmsub_pd(x,y,z))
#  define k8nmsub(x,y,z) (_mm_nmacc_pd(x,y,z))
#else
#  define k8madd(x,y,z)  (k8add(k8mul(x,y),z))
#  define k8msub(x,y,z)  (k8sub(k8mul(x,y),z))
#  define k8nmadd(x,y,z) (k8sub(k8neg(z),k8mul(x,y)))
#  define k8nmsub(x,y,z) (k8sub(z,k8mul(x,y)))
#endif

// Cheap functions
#define k8copysign(x,y)                                         \
  (_mm_or_pd(_mm_andnot_pd(k8sign_mask.vf,x),                   \
             _mm_and_pd(k8sign_mask.vf,y)))
#define k8fabs(x)   (_mm_andnot_pd(k8sign_mask.vf,x))
#define k8fmax(x,y) (_mm_max_pd(x,y))
#define k8fmin(x,y) (_mm_min_pd(x,y))
#define k8fnabs(x)  (_mm_or_pd(k8sign_mask.vf,x))
#define k8sgn(x_)                                                       \
  ({                                                                    \
    CCTK_REAL_VEC const x__=(x_);                                       \
    CCTK_REAL_VEC const x=x__;                                          \
    CCTK_REAL_VEC const iszero = _mm_cmpeq_pd(vec8_set1(0.0), x);       \
    CCTK_REAL_VEC const sign = _mm_and_pd(k8sign_mask.vf, x);           \
    CCTK_REAL_VEC const signedone = _mm_or_pd(vec8_set1(1.0), sign);    \
    k8ifthen(iszero, vec8_set1(0.0), signedone);                        \
  })
#define k8sqrt(x)   (_mm_sqrt_pd(x))

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
#define k8log(x)     K8REPL(log,x)
#define k8pow(x,a)   K8REPL2S(pow,x,a)
#define k8sin(x)     K8REPL(sin,x)
#define k8sinh(x)    K8REPL(sinh,x)
#define k8tan(x)     K8REPL(tan,x)
#define k8tanh(x)    K8REPL(tanh,x)

static const k8const_t k8lfalse_ = {{ +0LL, +0LL, }};
static const k8const_t k8ltrue_  = {{ -1LL, -1LL, }};
#define k8lfalse (k8lfalse_.vf)
#define k8ltrue  (k8ltrue_.vf)
#define k8lnot(x)   (_mm_xor_pd(k8ltrue,x))
#define k8land(x,y) (_mm_and_pd(x,y))
#define k8lor(x,y)  (_mm_or_pd(x,y))
#define k8lxor(x,y) (_mm_xor_pd(x,y))

#ifdef __SSE4_1__
#  define k8ifthen(x,y,z) (_mm_blendv_pd(z,y,x))
#elif 0
// This is slow (but this is what Intel/PGI produce by themselves)
#  define k8ifthen(x_,y_,z_)                    \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const y__=(y_);              \
    CCTK_REAL8_VEC const z__=(z_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    CCTK_REAL8_VEC const y=y__;                 \
    CCTK_REAL8_VEC const z=z__;                 \
    int const m = _mm_movemask_pd(x);           \
    CCTK_REAL8_VEC r;                           \
    switch (m) {                                \
    case 0: r = y; break;                       \
    case 1: r = _mm_move_sd(y,z); break;        \
    case 2: r = _mm_move_sd(z,y); break;        \
    case 3: r = z; break;                       \
    }                                           \
    r;                                          \
  })
#elif 0
#  ifdef __cplusplus
#    define k8signbit(x) ({ using namespace std; signbit(x); })
#  else
#    define k8signbit(x) (signbit(x))
#  endif
#  define k8ifthen(x_,y_,z_)                                            \
  ({                                                                    \
    CCTK_REAL8_VEC const x__=(x_);                                      \
    CCTK_REAL8_VEC const y__=(y_);                                      \
    CCTK_REAL8_VEC const z__=(z_);                                      \
    CCTK_REAL8_VEC const x=x__;                                         \
    CCTK_REAL8_VEC const y=y__;                                         \
    CCTK_REAL8_VEC const z=z__;                                         \
    vec8_set(k8signbit(vec8_elt0(x)) ? vec8_elt0(y) : vec8_elt0(z),     \
             k8signbit(vec8_elt1(x)) ? vec8_elt1(y) : vec8_elt1(z));    \
  })
#elif 0
// We don't need to shift -- the condition (mask) will be either all
// zeros or all ones
static const k8const_t k8ione  = {{ 0x1ULL, 0x1ULL, }};
#  define k8ifthen(x_,y_,z_)                                            \
  ({                                                                    \
    CCTK_REAL8_VEC const x__=(x_);                                      \
    CCTK_REAL8_VEC const y__=(y_);                                      \
    CCTK_REAL8_VEC const z__=(z_);                                      \
    CCTK_REAL8_VEC const x=x__;                                         \
    CCTK_REAL8_VEC const y=y__;                                         \
    CCTK_REAL8_VEC const z=z__;                                         \
    /* there is no _mm_srai_epi64(x, 63); we therefore calculate srli(x)-1 */ \
    __m128i const x_int = *(__m128i const*)&x;                          \
    __m128i const imask_int =                                           \
      _mm_sub_epi64(_mm_srli_epi64(x_int, 63), k8ione.vi);              \
    CCTK_REAL8_VEC const imask = *(CCTK_REAL8_VEC const*)&imask_int;    \
    /* (z & ~mask) | (y & mask)   where imask = ~mask */                \
    _mm_or_pd(_mm_and_pd(imask, z), _mm_andnot_pd(imask, y));           \
  })
#else
#  define k8ifthen(x_,y_,z_)                                    \
  ({                                                            \
    CCTK_REAL8_VEC const x__=(x_);                              \
    CCTK_REAL8_VEC const y__=(y_);                              \
    CCTK_REAL8_VEC const z__=(z_);                              \
    CCTK_REAL8_VEC const x=x__;                                 \
    CCTK_REAL8_VEC const y=y__;                                 \
    CCTK_REAL8_VEC const z=z__;                                 \
    /* (z & ~mask) | (y & mask)   where imask = ~mask */        \
    _mm_or_pd(_mm_and_pd(x, y), _mm_andnot_pd(x, z));           \
  })
#endif

#define k8cmpeq(x,y) (_mm_cmpeq_pd(x,y))
#define k8cmpne(x,y) (_mm_cmpneq_pd(x,y))
#define k8cmpgt(x,y) (_mm_cmpgt_pd(x,y))
#define k8cmpge(x,y) (_mm_cmpge_pd(x,y))
#define k8cmplt(x,y) (_mm_cmplt_pd(x,y))
#define k8cmple(x,y) (_mm_cmple_pd(x,y))
