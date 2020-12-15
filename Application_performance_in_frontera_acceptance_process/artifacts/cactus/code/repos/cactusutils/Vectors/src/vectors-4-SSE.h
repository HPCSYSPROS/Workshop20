// -*-C++-*-
// Vectorise using Intel's or AMD's SSE

// Use the type __m128 directly, without introducing a wrapper class

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>



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
#  include <x86intrin.h>
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
// Note: some boolean masks (e.g. ~0) correspond to nan when
// interpreted as floating point number. gcc 4.8 is clever enough to
// optimize away such constants with fast-math. We therefore need to
// handle this constant as integer number.
typedef __m128  CCTK_REAL4_VEC;
typedef __m128i CCTK_INTEGER4_VEC;
typedef __m128i CCTK_BOOLEAN4_VEC;

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL4_VEC_SIZE 4

vec_static_assert(sizeof(CCTK_REAL4_VEC) ==
                  sizeof(CCTK_REAL4) * CCTK_REAL4_VEC_SIZE);

// Integer and boolean types corresponding to this real type
typedef CCTK_INT4  CCTK_INTEGER4;
typedef CCTK_REAL4 CCTK_BOOLEAN4;



// These macros are undefined at the end of this file -- use them only
// within functions, not within macros that are exported
#define I2R(x) _mm_castsi128_ps(x)
#define R2I(x) _mm_castps_si128(x)



#define k4sign    (vec4_set1i(  (CCTK_INTEGER4)(1UL << 31UL)))
#define k4notsign (vec4_set1i(~ (CCTK_INTEGER4)(1UL << 31UL)))



// Create vectors, extract vector elements

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_set1(CCTK_REAL4 const a)
{
  return _mm_set1_ps(a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_INTEGER4_VEC vec4_set1i(CCTK_INT4 const a)
{
  return _mm_set1_epi32(a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_set(CCTK_REAL4 const a,
                        CCTK_REAL4 const b,
                        CCTK_REAL4 const c,
                        CCTK_REAL4 const d)
{
  return _mm_set_ps(d,c,b,a);   // note reversed arguments
}

// original order is 0123
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_swap1032(CCTK_REAL4_VEC const x)
{
  return _mm_shuffle_ps(x, x, _MM_SHUFFLE(2,3,0,1));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_swap2301(CCTK_REAL4_VEC const x)
{
  return _mm_shuffle_ps(x, x, _MM_SHUFFLE(1,0,3,2));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_swap3210(CCTK_REAL4_VEC const x)
{
  return _mm_shuffle_ps(x, x, _MM_SHUFFLE(0,1,2,3));
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4 vec4_elt(CCTK_REAL4_VEC const x, std::ptrdiff_t const d)
{
  CCTK_REAL4 e;
  std::memcpy(&e, &((char const*)&x)[d*sizeof e], sizeof e);
  return e;
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_INTEGER4 vec4_elti(CCTK_INTEGER4_VEC const x, std::ptrdiff_t const d)
{
  CCTK_INTEGER4 e;
  std::memcpy(&e, &((char const*)&x)[d*sizeof e], sizeof e);
  return e;
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4 vec4_eltb(CCTK_BOOLEAN4_VEC const x, std::ptrdiff_t const d)
{
  CCTK_BOOLEAN4 e;
  std::memcpy(&e, &((char const*)&x)[d*sizeof e], sizeof e);
  return e;
}



// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_load(CCTK_REAL4 const& p)
{
  return _mm_load_ps(&p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_loadu(CCTK_REAL4 const& p)
{
  return _mm_loadu_ps(&p);
}
#if ! VECTORISE_ALWAYS_USE_ALIGNED_LOADS
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_load_off1(CCTK_REAL4 const& p)
{
  return vec4_loadu(p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_load_off2(CCTK_REAL4 const& p)
{
  return vec4_loadu(p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_load_off3(CCTK_REAL4 const& p)
{
  return vec4_loadu(p);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_load_off1(CCTK_REAL4 const& p)
{
  CCTK_REAL4_VEC const lo = vec4_load((&p)[-1]);
  CCTK_REAL4_VEC const hi = vec4_load((&p)[+3]);
  CCTK_REAL4_VEC const hi2 = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(0,1,2,3));
  return _mm_shuffle_ps(lo, hi2, _MM_SHUFFLE(2,1,3,0));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_load_off2(CCTK_REAL4 const& p)
{
  CCTK_REAL4_VEC const lo = vec4_load((&p)[-2]);
  CCTK_REAL4_VEC const hi = vec4_load((&p)[+2]);
  return _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(1,0,3,2));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_load_off3(CCTK_REAL4 const& p)
{
  CCTK_REAL4_VEC const lo = vec4_load((&p)[-1]);
  CCTK_REAL4_VEC const hi = vec4_load((&p)[+3]);
  CCTK_REAL4_VEC const lo2 = _mm_shuffle_ps(lo, hi, _MM_SHUFFLE(0,1,2,3));
  return _mm_shuffle_ps(lo2, hi, _MM_SHUFFLE(3,0,2,1));
}
#endif

// Load a vector from memory that may or may not be aligned, as
// decided by the offset off and the vector size
#if VECTORISE_ALWAYS_USE_UNALIGNED_LOADS
// Implementation: Always use unaligned load
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_loadu_maybe(std::ptrdiff_t const off, CCTK_REAL4 const& p)
{
  return vec4_loadu(p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_loadu_maybe3(std::ptrdiff_t const off1,
                                 std::ptrdiff_t const off2,
                                 std::ptrdiff_t const off3,
                                 CCTK_REAL4 const& p)
{
  return vec4_loadu(p);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_loadu_maybe(std::ptrdiff_t const off, CCTK_REAL4 const& p)
{
  // The :? operator probably breaks with the Intel compiler
  //return off % CCTK_REAL4_VEC_SIZE == 0 ? vec4_load(p) : vec4_loadu(p);
  if (off % CCTK_REAL4_VEC_SIZE == 0) return vec4_load(p);
  return vec4_loadu(p);
}
#  if VECTORISE_ALIGNED_ARRAYS
// Assume all array x sizes are multiples of the vector size
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_loadu_maybe3(std::ptrdiff_t const off1,
                                 std::ptrdiff_t const off2,
                                 std::ptrdiff_t const off3,
                                 CCTK_REAL4 const& p)
{
  return vec4_loadu_maybe(off1, p);
}
#  else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_loadu_maybe3(std::ptrdiff_t const off1,
                                 std::ptrdiff_t const off2,
                                 std::ptrdiff_t const off3,
                                 CCTK_REAL4 const& p)
{
  return
    off2 % CCTK_REAL4_VEC_SIZE != 0 or
    off3 % CCTK_REAL4_VEC_SIZE != 0 ?
    vec4_loadu(p) :
    vec4_loadu_maybe(off1, p);
}
#  endif
#endif

// Store a vector to memory (aligned and non-temporal); this stores to
// a reference to a scalar
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store(CCTK_REAL4& p, CCTK_REAL4_VEC const x)
{
  return _mm_store_ps(&p, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_storeu(CCTK_REAL4& p, CCTK_REAL4_VEC const x)
{
  return _mm_storeu_ps(&p, x);
}
#if ! VECTORISE_STREAMING_STORES
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta(CCTK_REAL4& p, CCTK_REAL4_VEC const x)
{
  return vec4_store(p, x);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta(CCTK_REAL4& p, CCTK_REAL4_VEC const x)
{
  return _mm_stream_ps(&p, x);
}
#endif

// Store a partial vector (aligned and non-temporal)
// We ignoring VECTORISE_STREAMING_STORES for partial stores
#define vec4_store_partial_prepare(i, imin, imax)                       \
  std::ptrdiff_t v4stp_lo_skip, v4stp_hi_skip;                          \
  vec4_store_partial_prepare_(v4stp_lo_skip, v4stp_hi_skip, i, imin, imax)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_partial_prepare_(std::ptrdiff_t& lo_skip,
                                 std::ptrdiff_t& hi_skip,
                                 std::ptrdiff_t const i,
                                 std::ptrdiff_t const imin,
                                 std::ptrdiff_t const imax)
{
  lo_skip = std::max(std::ptrdiff_t(0), imin - i);
  hi_skip = std::max(std::ptrdiff_t(0), i+CCTK_REAL4_VEC_SIZE - imax);
}
#define vec4_store_nta_partial(p, x)                            \
  vec4_store_nta_partial_(v8stp_lo_skip, v8stp_hi_skip, p, x)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_(std::ptrdiff_t const lo_skip,
                             std::ptrdiff_t const hi_skip,
                             CCTK_REAL4& p,
                             CCTK_REAL4_VEC const x)
{
  if (CCTK_BUILTIN_EXPECT(lo_skip==0 and hi_skip==0, true)) {
    vec4_store_nta(p, x);
  } else {
    // these cases fall through
    switch (lo_skip) {
    case 0:
      (&p)[0] = vec4_elt(x, 0);
    case 1:
      if (hi_skip>=3) break;
      (&p)[1] = vec4_elt(x, 1);
    case 2:
      if (hi_skip>=2) break;
      (&p)[2] = vec4_elt(x, 2);
    case 3:
      if (hi_skip>=1) break;
      (&p)[3] = vec4_elt(x, 3);
    }
  }
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_lo(CCTK_REAL4& p,
                               CCTK_REAL4_VEC const x,
                               std::ptrdiff_t const n)
{
  // these cases fall through
  switch (n) {
  case 3: (&p)[2] = vec4_elt(x, 2);
  case 2: (&p)[1] = vec4_elt(x, 1);
  case 1: (&p)[0] = vec4_elt(x, 0);
  }
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_hi(CCTK_REAL4& p,
                               CCTK_REAL4_VEC const x,
                               std::ptrdiff_t const n)
{
  // these cases fall through
  switch (n) {
  case 3: (&p)[1]=vec4_elt(x, 1);
  case 2: (&p)[2]=vec4_elt(x, 2);
  case 1: (&p)[3]=vec4_elt(x, 3);
  }
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_hi(CCTK_REAL4& p,
                               CCTK_REAL4_VEC const x,
                               std::ptrdiff_t const nlo,
                               std::ptrdiff_t const nhi)
{
  // these cases fall through
  switch (nhi) {
  case 3:
    if (nlo<2) break;
    (&p)[1] = vec4_elt(x, 1);
  case 2:
    if (nlo<3) break;
    (&p)[2] = vec4_elt(x, 2);
  }
}



// Functions and operators

// Operators
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4neg(CCTK_REAL4_VEC const x)
{
  return _mm_xor_ps(I2R(k4sign), x);
}
// #define k4inv(x)
// TODO: provide k4inv via rcp and Newton-Raphson
// This is described in AMD's publication 47414.
// This should apply to AVX as well.

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4add(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm_add_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sub(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm_sub_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4mul(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm_mul_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4div(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  // TODO: maybe use k4inv and k4mul instead
  return _mm_div_ps(x, y);
}

// Fused multiply-add, defined as [+-]x*y[+-]z
#ifdef __FMA4__
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4madd(CCTK_REAL4_VEC const x,
                      CCTK_REAL4_VEC const y,
                      CCTK_REAL4_VEC const z)
{
  return _mm_macc_ps(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4msub(CCTK_REAL4_VEC const x,
                      CCTK_REAL4_VEC const y,
                      CCTK_REAL4_VEC const z)
{
  return _mm_msub_ps(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4nmadd(CCTK_REAL4_VEC const x,
                       CCTK_REAL4_VEC const y,
                       CCTK_REAL4_VEC const z)
{
  return _mm_nmsub_ps(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4nmsub(CCTK_REAL4_VEC const x,
                       CCTK_REAL4_VEC const y,
                       CCTK_REAL4_VEC const z)
{
  return _mm_nmacc_ps(x, y, z);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4madd(CCTK_REAL4_VEC const x,
                      CCTK_REAL4_VEC const y,
                      CCTK_REAL4_VEC const z)
{
  return k4add(k4mul(x, y), z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4msub(CCTK_REAL4_VEC const x,
                      CCTK_REAL4_VEC const y,
                      CCTK_REAL4_VEC const z)
{
  return k4sub(k4mul(x, y), z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4nmadd(CCTK_REAL4_VEC const x,
                       CCTK_REAL4_VEC const y,
                       CCTK_REAL4_VEC const z)
{
  return k4sub(k4neg(z), k4mul(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4nmsub(CCTK_REAL4_VEC const x,
                       CCTK_REAL4_VEC const y,
                       CCTK_REAL4_VEC const z)
{
  return k4sub(z, k4mul(x, y));
}
#endif

// Cheap functions
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4copysign(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm_or_ps(_mm_and_ps(I2R(k4notsign), x),
                   _mm_and_ps(I2R(k4sign   ), y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fabs(CCTK_REAL4_VEC const x)
{
  return _mm_and_ps(I2R(k4notsign), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fmax(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm_max_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fmin(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm_min_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fnabs(CCTK_REAL4_VEC const x)
{
  return _mm_or_ps(I2R(k4sign), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4signbit(CCTK_REAL4_VEC const x)
{
  return R2I(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sqrt(CCTK_REAL4_VEC const x)
{
  // TODO: maybe use rsqrt and Newton-Raphson
  return _mm_sqrt_ps(x);
}



// Expensive functions

#define K4REPL(f,x)                             \
  vec4_set(f(vec4_elt(x,0)),                    \
           f(vec4_elt(x,1)),                    \
           f(vec4_elt(x,2)),                    \
           f(vec4_elt(x,3)));
#define K4REPL2S(f,x,a)                         \
  vec4_set(f(vec4_elt(x,0),a),                  \
           f(vec4_elt(x,1),a),                  \
           f(vec4_elt(x,2),a),                  \
           f(vec4_elt(x,3),a));
#define K4REPL2(f,x,y)                          \
  vec4_set(f(vec4_elt(x,0),vec4_elt(y,0)),      \
           f(vec4_elt(x,1),vec4_elt(y,1)),      \
           f(vec4_elt(x,2),vec4_elt(y,2)),      \
           f(vec4_elt(x,3),vec4_elt(y,3)));

#if defined __ICC
// The Intel compiler provides intrinsics for these

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4acos(CCTK_REAL4_VEC const x)
{
  return _mm_acos_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4acosh(CCTK_REAL4_VEC const x)
{
  return _mm_acosh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4asin(CCTK_REAL4_VEC const x)
{
  return _mm_asin_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4asinh(CCTK_REAL4_VEC const x)
{
  return _mm_asinh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4atan(CCTK_REAL4_VEC const x)
{
  return _mm_atan_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4atan2(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm_atan2_ps(x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4atanh(CCTK_REAL4_VEC const x)
{
  return _mm_atanh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4cos(CCTK_REAL4_VEC const x)
{
  return _mm_cos_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4cosh(CCTK_REAL4_VEC const x)
{
  return _mm_cosh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4exp(CCTK_REAL4_VEC const x)
{
  return _mm_exp_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fmod(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
#if __ICC > 1400
  return _mm_fmod_ps(x,y);
#else
  return K4REPL2(fmodf,x,y);
#endif
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4log(CCTK_REAL4_VEC const x)
{
  return _mm_log_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4pow(CCTK_REAL4_VEC const x, CCTK_REAL4 const a)
{
  return _mm_pow_ps(x, _mm_set1_ps(a));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sin(CCTK_REAL4_VEC const x)
{
  return _mm_sin_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sinh(CCTK_REAL4_VEC const x)
{
  return _mm_sinh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4tan(CCTK_REAL4_VEC const x)
{
  return _mm_tan_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4tanh(CCTK_REAL4_VEC const x)
{
  return _mm_tanh_ps(x);
}

#else

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4acos(CCTK_REAL4_VEC const x)
{
  return K4REPL(acosf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4acosh(CCTK_REAL4_VEC const x)
{
  return K4REPL(acoshf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4asin(CCTK_REAL4_VEC const x)
{
  return K4REPL(asinf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4asinh(CCTK_REAL4_VEC const x)
{
  return K4REPL(asinhf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4atan(CCTK_REAL4_VEC const x)
{
  return K4REPL(atanf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4atan2(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return K4REPL2(atan2f,x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4atanh(CCTK_REAL4_VEC const x)
{
  return K4REPL(atanhf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4cos(CCTK_REAL4_VEC const x)
{
  return K4REPL(cosf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4cosh(CCTK_REAL4_VEC const x)
{
  return K4REPL(coshf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4exp(CCTK_REAL4_VEC const x)
{
  return K4REPL(expf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fmod(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return K4REPL2(fmodf,x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4log(CCTK_REAL4_VEC const x)
{
  return K4REPL(logf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4pow(CCTK_REAL4_VEC const x, CCTK_REAL4 const a)
{
  return K4REPL2S(powf,x,a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sin(CCTK_REAL4_VEC const x)
{
  return K4REPL(sinf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sinh(CCTK_REAL4_VEC const x)
{
  return K4REPL(sinhf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4tan(CCTK_REAL4_VEC const x)
{
  return K4REPL(tanf,x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4tanh(CCTK_REAL4_VEC const x)
{
  return K4REPL(tanhf,x);
}

#endif



#define k4lfalse (vec4_set1i( 0))
#define k4ltrue  (vec4_set1i(~0))
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4lnot(CCTK_BOOLEAN4_VEC const x)
{
  return R2I(_mm_xor_ps(I2R(k4ltrue), I2R(x)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4land(CCTK_BOOLEAN4_VEC const x, CCTK_BOOLEAN4_VEC const y)
{
  return R2I(_mm_and_ps(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4lor(CCTK_BOOLEAN4_VEC const x, CCTK_BOOLEAN4_VEC const y)
{
  return R2I(_mm_or_ps(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4lxor(CCTK_BOOLEAN4_VEC const x, CCTK_BOOLEAN4_VEC const y)
{
  return R2I(_mm_xor_ps(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4ifthen(CCTK_BOOLEAN4_VEC const x,
                        CCTK_REAL4_VEC const y,
                        CCTK_REAL4_VEC const z)
{
#ifdef __SSE4_1__
  return _mm_blendv_ps(z,y,I2R(x));
#elif 0
  return vec4_set(std::signbit(vec4_elt(x,0)) ? vec4_elt(y,0) : vec4_elt(z,0),
                  std::signbit(vec4_elt(x,1)) ? vec4_elt(y,1) : vec4_elt(z,1),
                  std::signbit(vec4_elt(x,2)) ? vec4_elt(y,2) : vec4_elt(z,2),
                  std::signbit(vec4_elt(x,3)) ? vec4_elt(y,3) : vec4_elt(z,3));
#elif 0
  // We don't need to shift -- the condition (mask) will be either all
  // zeros or all ones
  CCTK_REAL4_VEC const mask = (__m128)_mm_srai_epi32((__m128i)x, 31);
  // (z & ~mask) | (y & mask)
  return _mm_or_ps(_mm_andnot_ps(mask, z), _mm_and_ps(mask, y));
#else
  // This assumes that all logical operations always return either
  // lfalse or ltrue, and nothing "in between"
  // (z & ~mask) | (y & mask)   where imask = ~mask
  return _mm_or_ps(_mm_and_ps(I2R(x), y), _mm_andnot_ps(I2R(x), z));
#endif
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmpeq(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm_cmpeq_ps(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmpne(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm_cmpneq_ps(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmpgt(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm_cmpgt_ps(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmpge(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm_cmpge_ps(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmplt(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm_cmplt_ps(x, y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmple(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm_cmple_ps(x, y));
}



static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sgn(CCTK_REAL4_VEC const x)
{
  CCTK_BOOLEAN4_VEC const iszero    = k4cmpeq(x, vec4_set1(0.0));
  CCTK_REAL4_VEC    const sign      = _mm_and_ps(I2R(k4sign), x);
  CCTK_REAL4_VEC    const signedone = _mm_or_ps(sign, vec4_set1(1.0));
  return k4ifthen(iszero, vec4_set1(0.0), signedone);
}



#undef I2R
#undef R2I
