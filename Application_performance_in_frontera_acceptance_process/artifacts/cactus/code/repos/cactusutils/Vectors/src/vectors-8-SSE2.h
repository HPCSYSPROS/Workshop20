// -*-C++-*-
// Vectorise using Intel's or AMD's SSE2

// Use the type __m128d directly, without introducing a wrapper class

#include <cassert>
#include <cmath>



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
#ifdef __AVX__
#  include <immintrin.h>
#endif
#ifdef __FMA4__
#  include <x86intrin.h>
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
// Note: some boolean masks (e.g. ~0) correspond to nan when
// interpreted as floating point number. gcc 4.8 is clever enough to
// optimize away such constants with fast-math. We therefore need to
// handle this constant as integer number.
typedef __m128d CCTK_REAL8_VEC;
typedef __m128i CCTK_INTEGER8_VEC;
typedef __m128i CCTK_BOOLEAN8_VEC;

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL8_VEC_SIZE 2

vec_static_assert(sizeof(CCTK_REAL8_VEC) ==
                  sizeof(CCTK_REAL8) * CCTK_REAL8_VEC_SIZE);

// Integer and boolean types corresponding to this real type
typedef CCTK_INT8 CCTK_INTEGER8;
typedef CCTK_INT8 CCTK_BOOLEAN8;



// These macros are undefined at the end of this file -- use them only
// within functions, not within macros that are exported
#define I2R(x) _mm_castsi128_pd(x)
#define R2I(x) _mm_castpd_si128(x)



union k8const_t {
  CCTK_INTEGER8     i[CCTK_REAL8_VEC_SIZE];
  CCTK_INTEGER8_VEC vi;
};

#define k8sign    (vec8_set1i(  (CCTK_INTEGER8)(1ULL << 63ULL)))
#define k8notsign (vec8_set1i(~ (CCTK_INTEGER8)(1ULL << 63ULL)))



// Create vectors, extract vector elements

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_set1(CCTK_REAL8 const a)
{
  return _mm_set1_pd(a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_INTEGER8_VEC vec8_set1i(CCTK_INT8 const a)
{
#if defined(__INTEL_COMPILER)
  // Intel 11.1 does not support _mm_set1_epi64x
  return R2I(_mm_set1_pd(*(CCTK_REAL8 const*)&a));
#else
  return _mm_set1_epi64x(a);
#endif
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_set(CCTK_REAL8 const a, CCTK_REAL8 const b)
{
  return _mm_set_pd(b,a);       // note reversed arguments
}

// original order is 01
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_swap10(CCTK_REAL8_VEC const x)
{
  return _mm_shuffle_pd(x,x, _MM_SHUFFLE2(0,1));
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8 vec8_elt(CCTK_REAL8_VEC const x, std::ptrdiff_t const d)
{
  CCTK_REAL8 e;
  std::memcpy(&e, &((char const*)&x)[d*sizeof e], sizeof e);
  return e;
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_INTEGER8 vec8_elti(CCTK_INTEGER8_VEC const x, std::ptrdiff_t const d)
{
  CCTK_INTEGER8 e;
  std::memcpy(&e, &((char const*)&x)[d*sizeof e], sizeof e);
  return e;
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8 vec8_eltb(CCTK_BOOLEAN8_VEC const x, std::ptrdiff_t const d)
{
  CCTK_BOOLEAN8 e;
  std::memcpy(&e, &((char const*)&x)[d*sizeof e], sizeof e);
  return e;
}



// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_load(CCTK_REAL8 const& p)
{
  return _mm_load_pd(&p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_loadu(CCTK_REAL8 const& p)
{
  return _mm_loadu_pd(&p);
}
#if ! VECTORISE_ALWAYS_USE_ALIGNED_LOADS
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_load_off1(CCTK_REAL8 const& p)
{
  return vec8_loadu(p);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_load_off1(CCTK_REAL8 const& p)
{
  return _mm_shuffle_pd(vec8_load((&p)[-1]),
                        vec8_load((&p)[+1]), _MM_SHUFFLE2(0,1));
}
#endif

// Load a vector from memory that may or may not be aligned, as
// decided by the offset off and the vector size
#if VECTORISE_ALWAYS_USE_UNALIGNED_LOADS
// Implementation: Always use unaligned load
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_loadu_maybe(std::ptrdiff_t const off, CCTK_REAL8 const& p)
{
  return vec8_loadu(p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_loadu_maybe3(std::ptrdiff_t const off1,
                                 std::ptrdiff_t const off2,
                                 std::ptrdiff_t const off3,
                                 CCTK_REAL8 const& p)
{
  return vec8_loadu(p);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_loadu_maybe(std::ptrdiff_t const off, CCTK_REAL8 const& p)
{
  // Note: This function is mis-translated by some versions of the
  // Intel compiler. The symptom is a segfault during the Vectors
  // selftest in the vec_loadu_maybe test.
  // return off % CCTK_REAL8_VEC_SIZE == 0 ? vec8_load(p) : vec8_load_off1(p);
  if (off % CCTK_REAL8_VEC_SIZE == 0) return vec8_load(p);
  return vec8_load_off1(p);
}
#  if VECTORISE_ALIGNED_ARRAYS
// Assume all array x sizes are multiples of the vector size
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_loadu_maybe3(std::ptrdiff_t const off1,
                                 std::ptrdiff_t const off2,
                                 std::ptrdiff_t const off3,
                                 CCTK_REAL8 const& p)
{
  return vec8_loadu_maybe(off1, p);
}
#  else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_loadu_maybe3(std::ptrdiff_t const off1,
                                 std::ptrdiff_t const off2,
                                 std::ptrdiff_t const off3,
                                 CCTK_REAL8 const& p)
{
  return
    off2 % CCTK_REAL8_VEC_SIZE != 0 or
    off3 % CCTK_REAL8_VEC_SIZE != 0 ?
    vec8_loadu(p) :
    vec8_loadu_maybe(off1, p);
}
#  endif
#endif

// Store a vector to memory (aligned and non-temporal); this stores to
// a reference to a scalar
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store(CCTK_REAL8& p, CCTK_REAL8_VEC const x)
{
  _mm_store_pd(&p, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_storeu(CCTK_REAL8& p, CCTK_REAL8_VEC const x)
{
  _mm_storeu_pd(&p, x);
}
#if ! VECTORISE_STREAMING_STORES
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta(CCTK_REAL8& p, CCTK_REAL8_VEC const x)
{
  vec8_store(p, x);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta(CCTK_REAL8& p, CCTK_REAL8_VEC const x)
{
  // TODO: requires _mm_sfence() afterwards? requires _mm_lfence() in
  // readers afterwards? maybe better just an _mm_mfence() afterwards?
  _mm_stream_pd(&p, x);
}
#endif

// TODO: Use _mm_maskstore_pd if AVX is available?

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare(i, imin,imax)                        \
  bool v8stp_lo, v8stp_hi;                                              \
  vec8_store_partial_prepare_(v8stp_lo, v8stp_hi, i, imin, imax);
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_partial_prepare_(bool& lo, bool& hi,
                                 std::ptrdiff_t const i,
                                 std::ptrdiff_t const imin,
                                 std::ptrdiff_t const imax)
{
  lo = i >= imin;
  hi = i+CCTK_REAL8_VEC_SIZE-1 < imax;
}

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare_fixed(i, imin,imax)  \
  vec8_store_partial_prepare(i, imin,imax)

#define vec8_store_nta_partial(p, x)                    \
  vec8_store_nta_partial_(v8stp_lo, v8stp_hi, p, x)
#if VECTORISE_STREAMING_STORES && defined(__SSE4A__)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_(bool const lo, bool const hi,
                             CCTK_REAL8& p,
                             CCTK_REAL8_VEC const x)
{
  if (CCTK_BUILTIN_EXPECT(lo and hi, true)) {
    vec8_store_nta(p, x);
  } else if (lo) {
    _mm_stream_sd(&p, x);
  } else if (hi) {
    _mm_stream_sd(&p+1, vec8_swap10(x));
  }
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_(bool const lo, bool const hi,
                             CCTK_REAL8& p,
                             CCTK_REAL8_VEC const x)
{
  if (CCTK_BUILTIN_EXPECT(lo and hi, true)) {
    vec8_store_nta(p, x);
  } else if (lo) {
    _mm_storel_pd(&p, x);
  } else if (hi) {
    _mm_storeh_pd(&p+1, x);
  }
}
#endif

#define vec8_storeu_partial(p, x)                       \
  vec8_storeu_partial_(v8stp_lo, v8stp_hi, p, x)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_storeu_partial_(bool const lo, bool const hi,
                          CCTK_REAL8& p,
                          CCTK_REAL8_VEC const x)
{
  if (CCTK_BUILTIN_EXPECT(lo and hi, true)) {
    vec8_storeu(p, x);
  } else if (lo) {
    _mm_storel_pd(&p, x);
  } else if (hi) {
    _mm_storeh_pd(&p+1, x);
  }
}

// Store a lower or higher partial vector (aligned and non-temporal)
#if ! VECTORISE_STREAMING_STORES
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_lo(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm_storel_pd(&p, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_hi(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm_storeh_pd(&p+1, x);
}
#else
#  if defined(__SSE4A__)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_lo(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm_stream_sd(&p, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_hi(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm_stream_sd(&p+1, vec8_swap10(x));
}
#  else
// TODO: use clflush once a whole cache line has been written (cache
// lines are usually larger than the CPU vector size)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_lo(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm_storel_pd(&p, x);
  // _mm_clflush(&p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_hi(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm_storeh_pd(&p+1, x);
  // _mm_clflush(&p+1);
}
#  endif
#endif
#if 0
// This is slower; we would need a non-temporal read
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_lo(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  vec8_store_nta(p, _mm_loadh_pd(x, &p+1));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_hi(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  vec8_store_nta(p, _mm_loadl_pd(x, &p));
}
#endif
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_mid(CCTK_REAL8& p,
                                CCTK_REAL8_VEC const x,
                                ptrdiff_t const nlo,
                                ptrdiff_t const nhi)
{
  assert(0);
}



// Functions and operators

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

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8neg(CCTK_REAL8_VEC const x)
{
  return _mm_xor_pd(I2R(k8sign), x);
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8add(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm_add_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sub(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm_sub_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8mul(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm_mul_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8div(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm_div_pd(x, y);
}

// Fused multiply-add, defined as [+-]x*y[+-]z
#ifdef __FMA4__
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8madd(CCTK_REAL8_VEC const x,
                      CCTK_REAL8_VEC const y,
                      CCTK_REAL8_VEC const z)
{
  return _mm_macc_pd(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8msub(CCTK_REAL8_VEC const x,
                      CCTK_REAL8_VEC const y,
                      CCTK_REAL8_VEC const z)
{
  return _mm_msub_pd(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8nmadd(CCTK_REAL8_VEC const x,
                       CCTK_REAL8_VEC const y,
                       CCTK_REAL8_VEC const z)
{
  return _mm_nmsub_pd(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8nmsub(CCTK_REAL8_VEC const x,
                       CCTK_REAL8_VEC const y,
                       CCTK_REAL8_VEC const z)
{
  return _mm_nmacc_pd(x, y, z);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8madd(CCTK_REAL8_VEC const x,
                      CCTK_REAL8_VEC const y,
                      CCTK_REAL8_VEC const z)
{
  return k8add(k8mul(x, y), z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8msub(CCTK_REAL8_VEC const x,
                      CCTK_REAL8_VEC const y,
                      CCTK_REAL8_VEC const z)
{
  return k8sub(k8mul(x, y), z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8nmadd(CCTK_REAL8_VEC const x,
                       CCTK_REAL8_VEC const y,
                       CCTK_REAL8_VEC const z)
{
  return k8sub(k8neg(z), k8mul(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8nmsub(CCTK_REAL8_VEC const x,
                       CCTK_REAL8_VEC const y,
                       CCTK_REAL8_VEC const z)
{
  return k8sub(z, k8mul(x, y));
}
#endif

// Cheap functions
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8copysign(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm_or_pd(_mm_and_pd(I2R(k8notsign), x),
                   _mm_and_pd(I2R(k8sign   ), y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fabs(CCTK_REAL8_VEC const x)
{
  return _mm_and_pd(I2R(k8notsign), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmax(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm_max_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmin(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm_min_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fnabs(CCTK_REAL8_VEC const x)
{
  return _mm_or_pd(I2R(k8sign), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8signbit(CCTK_REAL8_VEC const x)
{
  return R2I(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sqrt(CCTK_REAL8_VEC const x)
{
  return _mm_sqrt_pd(x);
}



// Expensive functions

#define K8REPL(f,x)                             \
  vec8_set(f(vec8_elt(x,0)),                    \
           f(vec8_elt(x,1)));
#define K8REPL2S(f,x,a)                         \
  vec8_set(f(vec8_elt(x,0),a),                  \
           f(vec8_elt(x,1),a));
#define K8REPL2(f,x,y)                          \
  vec8_set(f(vec8_elt(x,0),vec8_elt(y,0)),      \
           f(vec8_elt(x,1),vec8_elt(y,1)));

#if defined __ICC
// The Intel compiler provides intrinsics for these

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8acos(CCTK_REAL8_VEC const x)
{
  return _mm_acos_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8acosh(CCTK_REAL8_VEC const x)
{
  return _mm_acosh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8asin(CCTK_REAL8_VEC const x)
{
  return _mm_asin_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8asinh(CCTK_REAL8_VEC const x)
{
  return _mm_asinh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atan(CCTK_REAL8_VEC const x)
{
  return _mm_atan_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atan2(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm_atan2_pd(x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atanh(CCTK_REAL8_VEC const x)
{
  return _mm_atanh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8cos(CCTK_REAL8_VEC const x)
{
  return _mm_cos_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8cosh(CCTK_REAL8_VEC const x)
{
  return _mm_cosh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8exp(CCTK_REAL8_VEC const x)
{
  return _mm_exp_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmod(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
#if __ICC > 1400
  return _mm_fmod_pd(x,y);
#else
  return K8REPL2(fmod,x,y);
#endif
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8log(CCTK_REAL8_VEC const x)
{
  return _mm_log_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8pow(CCTK_REAL8_VEC const x, CCTK_REAL8 const a)
{
  return _mm_pow_pd(x, _mm_set1_pd(a));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sin(CCTK_REAL8_VEC const x)
{
  return _mm_sin_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sinh(CCTK_REAL8_VEC const x)
{
  return _mm_sinh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8tan(CCTK_REAL8_VEC const x)
{
  return _mm_tan_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8tanh(CCTK_REAL8_VEC const x)
{
  return _mm_tanh_pd(x);
}

#else

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8acos(CCTK_REAL8_VEC const x)
{
  return K8REPL(acos,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8acosh(CCTK_REAL8_VEC const x)
{
  return K8REPL(acosh,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8asin(CCTK_REAL8_VEC const x)
{
  return K8REPL(asin,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8asinh(CCTK_REAL8_VEC const x)
{
  return K8REPL(asinh,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atan(CCTK_REAL8_VEC const x)
{
  return K8REPL(atan,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atan2(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return K8REPL2(atan2,x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atanh(CCTK_REAL8_VEC const x)
{
  return K8REPL(atanh,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8cos(CCTK_REAL8_VEC const x)
{
  return K8REPL(cos,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8cosh(CCTK_REAL8_VEC const x)
{
  return K8REPL(cosh,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8exp(CCTK_REAL8_VEC const x)
{
  return K8REPL(exp,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmod(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return K8REPL2(fmod,x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8log(CCTK_REAL8_VEC const x)
{
  return K8REPL(log,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8pow(CCTK_REAL8_VEC const x, CCTK_REAL8 const a)
{
  return K8REPL2S(pow,x,a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sin(CCTK_REAL8_VEC const x)
{
  return K8REPL(sin,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sinh(CCTK_REAL8_VEC const x)
{
  return K8REPL(sinh,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8tan(CCTK_REAL8_VEC const x)
{
  return K8REPL(tan,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8tanh(CCTK_REAL8_VEC const x)
{
  return K8REPL(tanh,x);
}

#endif



#define k8lfalse (vec8_set1i( 0))
#define k8ltrue  (vec8_set1i(~0))
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8lnot(CCTK_BOOLEAN8_VEC const x)
{
  return R2I(_mm_xor_pd(I2R(k8ltrue), I2R(x)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8land(CCTK_BOOLEAN8_VEC const x, CCTK_BOOLEAN8_VEC const y)
{
  return R2I(_mm_and_pd(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8lor(CCTK_BOOLEAN8_VEC const x, CCTK_BOOLEAN8_VEC const y)
{
  return R2I(_mm_or_pd(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8lxor(CCTK_BOOLEAN8_VEC const x, CCTK_BOOLEAN8_VEC const y)
{
  return R2I(_mm_xor_pd(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8ifthen(CCTK_BOOLEAN8_VEC const x,
                        CCTK_REAL8_VEC const y,
                        CCTK_REAL8_VEC const z)
{
#ifdef __SSE4_1__
  return _mm_blendv_pd(z,y,I2R(x));
#elif 0
  // This is slow (but this is what Intel/PGI produce by themselves)
  int const m = _mm_movemask_pd(x);
  switch (m) {
  case 0: return y;
  case 1: return _mm_move_sd(y,z);
  case 2: return _mm_move_sd(z,y);
  }
  return z;
#elif 0
  return vec8_set(std::signbit(vec8_elt(x,0)) ? vec8_elt(y,0) : vec8_elt(z,0),
                  std::signbit(vec8_elt(x,1)) ? vec8_elt(y,1) : vec8_elt(z,1));
#elif 0
  // We don't need to shift -- the condition (mask) will be either all
  // zeros or all ones
  k8const_t const k8ione  = {{ 1, 1, }};
  // there is no _mm_srai_epi64(x, 63); we therefore calculate srli(x)-1
  __m128i const x_int = *(__m128i const*)&x;
  __m128i const imask_int = _mm_sub_epi64(_mm_srli_epi64(x_int, 63), k8ione.vi);
  CCTK_REAL8_VEC const imask = *(CCTK_REAL8_VEC const*)&imask_int;
  // (z & ~mask) | (y & mask)   where imask = ~mask
  return _mm_or_pd(_mm_and_pd(imask, z), _mm_andnot_pd(imask, y));
#else
  // This assumes that all logical operations always return either
  // lfalse or ltrue, and nothing "in between"
  // (z & ~mask) | (y & mask)   where imask = ~mask
  return _mm_or_pd(_mm_and_pd(I2R(x), y), _mm_andnot_pd(I2R(x), z));
#endif
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpeq(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm_cmpeq_pd(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpne(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm_cmpneq_pd(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpgt(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm_cmpgt_pd(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpge(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm_cmpge_pd(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmplt(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm_cmplt_pd(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmple(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm_cmple_pd(x, y));
}



static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sgn(CCTK_REAL8_VEC const x)
{
  CCTK_BOOLEAN8_VEC const iszero    = k8cmpeq(x, vec8_set1(0.0));
  CCTK_REAL8_VEC    const sign      = _mm_and_pd(I2R(k8sign), x);
  CCTK_REAL8_VEC    const signedone = _mm_or_pd(sign, vec8_set1(1.0));
  return k8ifthen(iszero, vec8_set1(0.0), signedone);
}



#undef I2R
#undef R2I
