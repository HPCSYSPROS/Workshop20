// Vectorise using Intel's or AMD's AVX

// Use the type __m256d directly, without introducing a wrapper class
// Use macros instead of inline functions



#if VECTORISE_EMULATE_AVX
#  include "avxintrin_emu.h"
#else
#  include <immintrin.h>
#endif
#ifdef __FMA4__
#  include <fma4intrin.h>
#endif



#ifdef __FMA4__
#  define vec8_architecture_FMA4 "+FMA4"
#else
#  define vec8_architecture_FMA4 ""
#endif
#define vec8_architecture "AVX" vec8_architecture_FMA4 " (64-bit precision)"



// Vector type corresponding to CCTK_REAL
#define CCTK_REAL8_VEC __m256d

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL8_VEC_SIZE 4

// Integer and boolean types corresponding to this real type
#define CCTK_INTEGER8     CCTK_REAL8
#define CCTK_BOOLEAN8     CCTK_REAL8
#define CCTK_INTEGER8_VEC CCTK_REAL8_VEC
#define CCTK_BOOLEAN8_VEC CCTK_REAL8_VEC



union k8const_t {
  unsigned long long i[4];
  double             f[4];
  __m256i            vi;
  __m256d            vf;
};

#define K8_ZERO    0x0000000000000000ULL
#define K8_NOTZERO 0xffffffffffffffffULL
#define K8_IMIN    0x8000000000000000ULL
#define K8_IMAX    0x7fffffffffffffffULL



// Create vectors, extract vector elements

#define vec8_set1(a)      (_mm256_set1_pd(a))
#define vec8_set(a,b,c,d) (_mm256_set_pd(d,c,b,a)) // note reversed arguments

#define vec8_elt0(x) (((CCTK_REAL8 const*)&(x))[0])
#define vec8_elt1(x) (((CCTK_REAL8 const*)&(x))[1])
#define vec8_elt2(x) (((CCTK_REAL8 const*)&(x))[2])
#define vec8_elt3(x) (((CCTK_REAL8 const*)&(x))[3])
#define vec8_elt(x,d) (((CCTK_REAL8 const*)&(x))[d])



// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
#define vec8_load(p)  (_mm256_load_pd(&(p)))
#define vec8_loadu(p) (_mm256_loadu_pd(&(p)))
#if ! VECTORISE_ALWAYS_USE_ALIGNED_LOADS
#  define vec8_load_off1(p) vec_loadu(p)
#else
#  error "VECTORISE_ALWAYS_USE_ALIGNED_LOADS not yet supported"
#endif

// Load a vector from memory that may or may not be aligned, as
// decided by the offset off and the vector size
#if VECTORISE_ALWAYS_USE_UNALIGNED_LOADS
// Implementation: Always use unaligned load
#  define vec8_loadu_maybe(off,p)             (vec8_loadu(p))
#  define vec8_loadu_maybe3(off1,off2,off3,p) (vec8_loadu(p))
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
#define vec8_store(p,x)  (_mm256_store_pd(&(p),x))
#define vec8_storeu(p,x) (_mm256_storeu_pd(&(p),x))
#if ! VECTORISE_STREAMING_STORES
#  define vec8_store_nta(p,x) (vec8_store(p,x))
#else
#  define vec8_store_nta(p,x) (_mm256_stream_pd(&(p),x))
#endif

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare(i,imin_,imax_)                       \
  bool v8stp_all;                                                       \
  __m256i v8stp_mask;                                                   \
  ({                                                                    \
    ptrdiff_t const imin__=(imin_);                                     \
    ptrdiff_t const imin=imin__;                                        \
    ptrdiff_t const imax__=(imax_);                                     \
    ptrdiff_t const imax=imax__;                                        \
                                                                        \
    v8stp_all = i>=imin and i+CCTK_REAL_VEC_SIZE-1<imax;                \
                                                                        \
    if (not CCTK_BUILTIN_EXPECT(v8stp_all, true)) {                     \
      /*                                                                \
        __m256i const v8stp_mask =                                      \
          _mm256_andnot_pd(_mm256_add_epi64(_mm256_set1_epi64x(i-imin), \
                                            vec_index),                 \
                           _mm256_add_epi64(_mm256_set1_epi64x(i-imax), \
                                            vec_index));                \
      */                                                                \
      __m128i const termlo0 =                                           \
        _mm_add_epi64(_mm_set1_epi64x(i-imin), _mm_set_epi64x(1, 0));   \
      __m128i const termup0 =                                           \
        _mm_add_epi64(_mm_set1_epi64x(i-imax), _mm_set_epi64x(1, 0));   \
      __m128i const term0 = _mm_andnot_si128(termlo0, termup0);         \
      __m128i const termlo1 =                                           \
        _mm_add_epi64(_mm_set1_epi64x(i-imin), _mm_set_epi64x(3, 2));   \
      __m128i const termup1 =                                           \
        _mm_add_epi64(_mm_set1_epi64x(i-imax), _mm_set_epi64x(3, 2));   \
      __m128i const term1 = _mm_andnot_si128(termlo1, termup1);         \
      v8stp_mask =                                                      \
        _mm256_insertf128_si256(_mm256_castsi128_si256(term0), term1, 1); \
    }                                                                   \
  })

#define vec8_store_nta_partial(p,x)             \
  ({                                            \
    if (CCTK_BUILTIN_EXPECT(v8stp_all, true)) { \
      vec8_store_nta(p,x);                      \
    } else {                                    \
      _mm256_maskstore_pd(&p,v8stp_mask,x);     \
    }                                           \
  })

// Store a lower or higher partial vector (aligned and non-temporal);
// the non-temporal hint is probably ignored
// Masks indicating which vector element should be stored:
static const k8const_t k8store_lo[5] =
  {
    {{ K8_ZERO   , K8_ZERO   , K8_ZERO   , K8_ZERO   , }},
    {{ K8_NOTZERO, K8_ZERO   , K8_ZERO   , K8_ZERO   , }},
    {{ K8_NOTZERO, K8_NOTZERO, K8_ZERO   , K8_ZERO   , }},
    {{ K8_NOTZERO, K8_NOTZERO, K8_NOTZERO, K8_ZERO   , }},
    {{ K8_NOTZERO, K8_NOTZERO, K8_NOTZERO, K8_NOTZERO, }},
  };
static const k8const_t k8store_hi[5] =
  {
    {{ K8_ZERO   , K8_ZERO   , K8_ZERO   , K8_ZERO   , }},
    {{ K8_ZERO   , K8_ZERO   , K8_ZERO   , K8_NOTZERO, }},
    {{ K8_ZERO   , K8_ZERO   , K8_NOTZERO, K8_NOTZERO, }},
    {{ K8_ZERO   , K8_NOTZERO, K8_NOTZERO, K8_NOTZERO, }},
    {{ K8_NOTZERO, K8_NOTZERO, K8_NOTZERO, K8_NOTZERO, }},
  };
#if !defined(__INTEL_COMPILER) && defined(__GNUC__) && __GNUC__==4 && __GNUC_MINOR__<=4
// gcc 4.4 uses a wrong prototype for _mm256_maskstore_pd
#  define vec8_store_nta_partial_lo(p,x,n)                              \
  (_mm256_maskstore_pd(&(p),_mm256_castsi256_pd(k8store_lo[n].vi),x))
#  define vec8_store_nta_partial_hi(p,x,n)                              \
  (_mm256_maskstore_pd(&(p),_mm256_castsi256_pd(k8store_hi[n].vi),x))
#  define vec8_store_nta_partial_mid(p,x,nlo,nhi)                       \
  (_mm256_maskstore_pd                                                  \
   (&(p),                                                               \
    _mm256_castsi256_pd(k8store_lo[nlo].vi & k8store_hi[nhi].vi),       \
    x))
#else
#  define vec8_store_nta_partial_lo(p,x,n)              \
  (_mm256_maskstore_pd(&(p),k8store_lo[n].vi,x))
#  define vec8_store_nta_partial_hi(p,x,n)              \
  (_mm256_maskstore_pd(&(p),k8store_hi[n].vi,x))
#  define vec8_store_nta_partial_mid(p,x,nlo,nhi)               \
  (_mm256_maskstore_pd                                          \
   (&(p),                                                       \
    _mm256_castpd_si256(_mm256_and_pd(k8store_lo[nlo].vf,       \
                                      k8store_hi[nhi].vf)),     \
    x))
#endif



// Functions and operators

static const k8const_t k8sign_mask = {{ K8_IMIN, K8_IMIN, K8_IMIN, K8_IMIN, }};

// Operators
#define k8neg(x) (_mm256_xor_pd(x,k8sign_mask.vf))

#define k8add(x,y) (_mm256_add_pd(x,y))
#define k8sub(x,y) (_mm256_sub_pd(x,y))
#define k8mul(x,y) (_mm256_mul_pd(x,y))
#define k8div(x,y) (_mm256_div_pd(x,y))

// Fused multiply-add, defined as [+-]x*y[+-]z
#ifdef __FMA4__
#  define k8madd(x,y,z)  (_mm256_macc_pd(x,y,z))
#  define k8msub(x,y,z)  (_mm256_msub_pd(x,y,z))
#  define k8nmadd(x,y,z) (_mm256_nmsub_pd(x,y,z))
#  define k8nmsub(x,y,z) (_mm256_nmacc_pd(x,y,z))
#else
#  define k8madd(x,y,z)  (k8add(k8mul(x,y),z))
#  define k8msub(x,y,z)  (k8sub(k8mul(x,y),z))
#  define k8nmadd(x,y,z) (k8sub(k8neg(z),k8mul(x,y)))
#  define k8nmsub(x,y,z) (k8sub(z,k8mul(x,y)))
#endif

// Cheap functions
#define k8copysign(x,y)                                 \
  (_mm256_or_pd(_mm256_andnot_pd(k8sign_mask.vf,x),     \
                _mm256_and_pd(k8sign_mask.vf,y)))
#define k8fabs(x)   (_mm256_andnot_pd(k8sign_mask.vf,x))
#define k8fmax(x,y) (_mm256_max_pd(x,y))
#define k8fmin(x,y) (_mm256_min_pd(x,y))
#define k8fnabs(x)  (_mm256_or_pd(x,k8sign_mask.vf))
static const k8const_t k8zero = { f: { 0.0, 0.0, 0.0, 0.0, }};
static const k8const_t k8one  = { f: { 1.0, 1.0, 1.0, 1.0, }};
#define k8sgn(x_)                                                       \
  ({                                                                    \
    CCTK_REAL_VEC x__=(x_);                                             \
    CCTK_REAL_VEC x=x__;                                                \
    CCTK_REAL_VEC iszero = _mm256_cmp_pd(x, k8zero.vf, _CMP_EQ_OQ);     \
    CCTK_REAL_VEC sign = _mm256_and_pd(k8sign_mask.vf, x);              \
    CCTK_REAL_VEC signedone = _mm256_or_pd(sign, k8one.vf);             \
    k8ifthen(iszero, k8zero.vf, signedone);                             \
  })
#define k8sqrt(x)   (_mm256_sqrt_pd(x))

// Expensive functions
#define K8REPL(f,x_)                            \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    vec8_set(f(vec8_elt0(x)),                   \
             f(vec8_elt1(x)),                   \
             f(vec8_elt2(x)),                   \
             f(vec8_elt3(x)));                  \
  })
#define K8REPL2S(f,x_,a_)                       \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8     const a__=(a_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    CCTK_REAL8     const a=a__;                 \
    vec8_set(f(vec8_elt0(x),a),                 \
             f(vec8_elt1(x),a),                 \
             f(vec8_elt2(x),a),                 \
             f(vec8_elt3(x),a));                \
  })
#define K8REPL2(f,x_,y_)                        \
  ({                                            \
    CCTK_REAL8_VEC const x__=(x_);              \
    CCTK_REAL8_VEC const y__=(y_);              \
    CCTK_REAL8_VEC const x=x__;                 \
    CCTK_REAL8_VEC const y=y__;                 \
    vec8_set(f(vec8_elt0(x),vec8_elt0(y)),      \
             f(vec8_elt1(x),vec8_elt1(y)),      \
             f(vec8_elt2(x),vec8_elt2(y)),      \
             f(vec8_elt3(x),vec8_elt3(y)));     \
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

static const k8const_t k8lfalse_ =
  {{ K8_ZERO, K8_ZERO, K8_ZERO, K8_ZERO, }};
static const k8const_t k8ltrue_  =
  {{ K8_NOTZERO, K8_NOTZERO, K8_NOTZERO, K8_NOTZERO, }};
#define k8lfalse        (k8lfalse_.vf)
#define k8ltrue         (k8ltrue_.vf)
#define k8lnot(x)       (_mm256_xor_pd(k8ltrue,x))
#define k8land(x,y)     (_mm256_and_pd(x,y))
#define k8lor(x,y)      (_mm256_or_pd(x,y))
#define k8lxor(x,y)     (_mm256_xor_pd(x,y))
#define k8ifthen(x,y,z) (_mm256_blendv_pd(z,y,x))

#define k8cmpeq(x,y) (_mm256_cmp_pd(x,y,_CMP_EQ_OQ))
#define k8cmpne(x,y) (_mm256_cmp_pd(x,y,_CMP_NEQ_OQ))
#define k8cmpgt(x,y) (_mm256_cmp_pd(x,y,_CMP_GT_OQ))
#define k8cmpge(x,y) (_mm256_cmp_pd(x,y,_CMP_GE_OQ))
#define k8cmplt(x,y) (_mm256_cmp_pd(x,y,_CMP_LT_OQ))
#define k8cmple(x,y) (_mm256_cmp_pd(x,y,_CMP_LE_OQ))
