// -*-C++-*-
// Vectorise using Intel's or AMD's AVX

// Use the type __m256 directly, without introducing a wrapper class



#include <cstdlib>
#include <cstring>



#include <immintrin.h>
#ifdef __FMA4__
#  include <x86intrin.h>
#endif



#ifdef __FMA4__
#  define vec4_architecture_FMA4 "+FMA4"
#else
#  define vec4_architecture_FMA4 ""
#endif
#define vec4_architecture "AVX" vec4_architecture_FMA4 " (32-bit precision)"



// Vector type corresponding to CCTK_REAL
// Note: some boolean masks (e.g. ~0) correspond to nan when
// interpreted as floating point number. gcc 4.8 is clever enough to
// optimize away such constants with fast-math. We therefore need to
// handle this constant as integer number.
typedef __m256  CCTK_REAL4_VEC;
typedef __m256i CCTK_INTEGER4_VEC;
typedef __m256i CCTK_BOOLEAN4_VEC;

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL4_VEC_SIZE 8

vec_static_assert(sizeof(CCTK_REAL4_VEC) ==
                  sizeof(CCTK_REAL4) * CCTK_REAL4_VEC_SIZE);

// Integer and boolean types corresponding to this real type
typedef CCTK_INT4 CCTK_INTEGER4;
typedef CCTK_INT4 CCTK_BOOLEAN4;



// These macros are undefined at the end of this file -- use them only
// within functions, not within macros that are exported
#define I2R(x) _mm256_castsi256_ps(x)
#define R2I(x) _mm256_castps_si256(x)



union k4const_t {
  CCTK_INTEGER4     i[CCTK_REAL4_VEC_SIZE];
  CCTK_INTEGER4_VEC vi;
};

#define k4sign    (vec4_set1i(  (CCTK_INTEGER4)(1UL << 31UL)))
#define k4notsign (vec4_set1i(~ (CCTK_INTEGER4)(1UL << 31UL)))



// Create vectors, extract vector elements

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_set1(CCTK_REAL4 const a)
{
  return _mm256_set1_ps(a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_INTEGER4_VEC vec4_set1i(CCTK_INT4 const a)
{
  return _mm256_set1_epi32(a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_set(CCTK_REAL4 const a,
                        CCTK_REAL4 const b,
                        CCTK_REAL4 const c,
                        CCTK_REAL4 const d,
                        CCTK_REAL4 const e,
                        CCTK_REAL4 const f,
                        CCTK_REAL4 const g,
                        CCTK_REAL4 const h)
{
  return _mm256_set_ps(h,g,f,e,d,c,b,a); // note reversed arguments
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
  return _mm256_load_ps(&p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC vec4_loadu(CCTK_REAL4 const& p)
{
  return _mm256_loadu_ps(&p);
}
#if VECTORISE_ALWAYS_USE_ALIGNED_LOADS
#  error "VECTORISE_ALWAYS_USE_ALIGNED_LOADS not yet supported"
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
  return off % CCTK_REAL4_VEC_SIZE == 0 ? vec4_load(p) : vec4_loadu(p);
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
  return _mm256_store_ps(&p, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_storeu(CCTK_REAL4& p, CCTK_REAL4_VEC const x)
{
  return _mm256_storeu_ps(&p, x);
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
  return _mm256_stream_ps(&p, x);
}
#endif

// Store a partial vector (aligned and non-temporal)
#define vec4_store_partial_prepare(i,imin,imax)                         \
  bool v4stp_all;                                                       \
  __m256i v4stp_mask;                                                   \
  vec4_store_partial_prepare_(v4stp_all, v4stp_mask, i, imin, imax);
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_partial_prepare_(bool& all, __m256i& mask,
                                 std::ptrdiff_t const i,
                                 std::ptrdiff_t const imin,
                                 std::ptrdiff_t const imax)
{
  all = i>=imin and i+CCTK_REAL4_VEC_SIZE-1<imax;
  
  if (not CCTK_BUILTIN_EXPECT(all, true)) {
    /* __m256i const v4stp_mask = */
    /*   _mm256_andnot_ps(_mm256_add_epi64(_mm256_set1_epi64x(i-imin), */
    /*                                     vec_index), */
    /*                    _mm256_add_epi64(_mm256_set1_epi64x(i-imax), */
    /*                                     vec_index)); */
    __m128i const termlo0123 =
      _mm_add_epi32(_mm_set1_epi32(i-imin), _mm_set_epi32(3, 2, 1, 0));
    __m128i const termup0123 =
      _mm_add_epi32(_mm_set1_epi32(i-imax), _mm_set_epi32(3, 2, 1, 0));
    __m128i const term0123 = _mm_andnot_si128(termlo0123, termup0123);
    __m128i const termlo4567 =
      _mm_add_epi32(_mm_set1_epi32(i-imin), _mm_set_epi32(7, 6, 5, 4));
    __m128i const termup4567 =
      _mm_add_epi32(_mm_set1_epi32(i-imax), _mm_set_epi32(7, 6, 5, 4));
    __m128i const term4567 = _mm_andnot_si128(termlo4567, termup4567);
    mask =
      _mm256_insertf128_si256(_mm256_castsi128_si256(term0123), term4567, 1);
  }
}

#define vec4_store_nta_partial(p, x)                    \
  vec4_store_nta_partial_(v4stp_all, v4stp_mask, p, x)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_(bool const all, __m256i const mask,
                             CCTK_REAL4& p, CCTK_REAL4_VEC const x)
{
  if (CCTK_BUILTIN_EXPECT(all, true)) {
    vec4_store_nta(p, x);
  } else {
    _mm256_maskstore_ps(&p, mask, x);
  }
}

// Store a lower or higher partial vector (aligned and non-temporal);
// the non-temporal hint is probably ignored
// Masks indicating which vector element should be stored:
static k4const_t const k4store_lo[9] =
  {
    {{  0,  0,  0,  0,  0,  0,  0,  0, }},
    {{ ~0,  0,  0,  0,  0,  0,  0,  0, }},
    {{ ~0, ~0,  0,  0,  0,  0,  0,  0, }},
    {{ ~0, ~0, ~0,  0,  0,  0,  0,  0, }},
    {{ ~0, ~0, ~0, ~0,  0,  0,  0,  0, }},
    {{ ~0, ~0, ~0, ~0, ~0,  0,  0,  0, }},
    {{ ~0, ~0, ~0, ~0, ~0, ~0,  0,  0, }},
    {{ ~0, ~0, ~0, ~0, ~0, ~0, ~0,  0, }},
    {{ ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, }},
  };
static k4const_t const k4store_hi[9] =
  {
    {{  0,  0,  0,  0,  0,  0,  0,  0, }},
    {{  0,  0,  0,  0,  0,  0,  0, ~0, }},
    {{  0,  0,  0,  0,  0,  0, ~0, ~0, }},
    {{  0,  0,  0,  0,  0, ~0, ~0, ~0, }},
    {{  0,  0,  0,  0, ~0, ~0, ~0, ~0, }},
    {{  0,  0,  0, ~0, ~0, ~0, ~0, ~0, }},
    {{  0,  0, ~0, ~0, ~0, ~0, ~0, ~0, }},
    {{  0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, }},
    {{ ~0, ~0, ~0, ~0, ~0, ~0, ~0, ~0, }},
  };
#if !defined(__INTEL_COMPILER) && defined(__GNUC__) && __GNUC__==4 && __GNUC_MINOR__<=4
// gcc 4.4 uses a wrong prototype for _mm256_maskstore_ps
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_lo(CCTK_REAL4& p,
                               CCTK_REAL4_VEC const x,
                               ptrdiff_t const n)
{
  _mm256_maskstore_ps(&p, I2R(k4store_lo[n].vi), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_hi(CCTK_REAL4& p,
                               CCTK_REAL4_VEC const x,
                               ptrdiff_t const n)
{
  _mm256_maskstore_ps(&p, I2R(k4store_hi[n].vi), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_mid(CCTK_REAL4& p,
                                CCTK_REAL4_VEC const x,
                                ptrdiff_t const nlo,
                                ptrdiff_t const nhi)
{
  _mm256_maskstore_ps(&p, I2R(k4store_lo[nlo].vi & k4store_hi[nhi].vi), x);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_lo(CCTK_REAL4& p,
                               CCTK_REAL4_VEC const x,
                               ptrdiff_t const n)
{
  _mm256_maskstore_ps(&p, k4store_lo[n].vi, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_hi(CCTK_REAL4& p,
                               CCTK_REAL4_VEC const x,
                               ptrdiff_t const n)
{
  _mm256_maskstore_ps(&p, k4store_hi[n].vi, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec4_store_nta_partial_mid(CCTK_REAL4& p,
                                CCTK_REAL4_VEC const x,
                                ptrdiff_t const nlo,
                                ptrdiff_t const nhi)
{
  _mm256_maskstore_ps
    (&p,
     R2I(_mm256_and_ps(I2R(k4store_lo[nlo].vi), I2R(k4store_hi[nhi].vi))),
     x);
}
#endif



// Functions and operators

// Operators
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4neg(CCTK_REAL4_VEC const x)
{
  return _mm256_xor_ps(x, I2R(k4sign));
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4add(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm256_add_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sub(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm256_sub_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4mul(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm256_mul_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4div(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm256_div_ps(x, y);
}

// Fused multiply-add, defined as [+-]x*y[+-]z
#ifdef __FMA4__
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4madd(CCTK_REAL4_VEC const x,
                      CCTK_REAL4_VEC const y,
                      CCTK_REAL4_VEC const z)
{
  return _mm256_macc_ps(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4msub(CCTK_REAL4_VEC const x,
                      CCTK_REAL4_VEC const y,
                      CCTK_REAL4_VEC const z)
{
  return _mm256_msub_ps(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4nmadd(CCTK_REAL4_VEC const x,
                       CCTK_REAL4_VEC const y,
                       CCTK_REAL4_VEC const z)
{
  return _mm256_nmsub_ps(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4nmsub(CCTK_REAL4_VEC const x,
                       CCTK_REAL4_VEC const y,
                       CCTK_REAL4_VEC const z)
{
  return _mm256_nmacc_ps(x, y, z);
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
  return _mm256_or_ps(_mm256_and_ps(I2R(k4notsign), x),
                      _mm256_and_ps(I2R(k4sign   ), y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fabs(CCTK_REAL4_VEC const x)
{
  return _mm256_and_ps(I2R(k4notsign), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fmax(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm256_max_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fmin(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm256_min_ps(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fnabs(CCTK_REAL4_VEC const x)
{
  return _mm256_or_ps(x, I2R(k4sign));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4signbit(CCTK_REAL4_VEC const x)
{
  return R2I(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sqrt(CCTK_REAL4_VEC const x)
{
  return _mm256_sqrt_ps(x);
}



// Expensive functions

#define K4REPL(f,x)                             \
  vec4_set(f(vec4_elt(x,0)),                    \
           f(vec4_elt(x,1)),                    \
           f(vec4_elt(x,2)),                    \
           f(vec4_elt(x,3)),                    \
           f(vec4_elt(x,4)),                    \
           f(vec4_elt(x,5)),                    \
           f(vec4_elt(x,6)),                    \
           f(vec4_elt(x,7)));
#define K4REPL2S(f,x,a)                         \
  vec4_set(f(vec4_elt(x,0),a),                  \
           f(vec4_elt(x,1),a),                  \
           f(vec4_elt(x,2),a),                  \
           f(vec4_elt(x,3),a),                  \
           f(vec4_elt(x,4),a),                  \
           f(vec4_elt(x,5),a),                  \
           f(vec4_elt(x,6),a),                  \
           f(vec4_elt(x,7),a));
#define K4REPL2(f,x,y)                          \
  vec4_set(f(vec4_elt(x,0),vec4_elt(y,0)),      \
           f(vec4_elt(x,1),vec4_elt(y,1)),      \
           f(vec4_elt(x,2),vec4_elt(y,2)),      \
           f(vec4_elt(x,3),vec4_elt(y,3)),      \
           f(vec4_elt(x,4),vec4_elt(y,4)),      \
           f(vec4_elt(x,5),vec4_elt(y,5)),      \
           f(vec4_elt(x,6),vec4_elt(y,6)),      \
           f(vec4_elt(x,7),vec4_elt(y,7)));

#if defined __ICC
// The Intel compiler provides intrinsics for these

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4acos(CCTK_REAL4_VEC const x)
{
  return _mm256_acos_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4acosh(CCTK_REAL4_VEC const x)
{
  return _mm256_acosh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4asin(CCTK_REAL4_VEC const x)
{
  return _mm256_asin_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4asinh(CCTK_REAL4_VEC const x)
{
  return _mm256_asinh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4atan(CCTK_REAL4_VEC const x)
{
  return _mm256_atan_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4atan2(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return _mm256_atan2_ps(x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4atanh(CCTK_REAL4_VEC const x)
{
  return _mm256_atanh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4cos(CCTK_REAL4_VEC const x)
{
  return _mm256_cos_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4cosh(CCTK_REAL4_VEC const x)
{
  return _mm256_cosh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4exp(CCTK_REAL4_VEC const x)
{
  return _mm256_exp_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4fmod(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
#if __ICC > 1400
  return _mm256_fmod_ps(x,y);
#else
  return K4REPL2(fmodf,x,y);
#endif
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4log(CCTK_REAL4_VEC const x)
{
  return _mm256_log_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4pow(CCTK_REAL4_VEC const x, CCTK_REAL4 const a)
{
  return _mm256_pow_ps(x, _mm256_set1_ps(a));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sin(CCTK_REAL4_VEC const x)
{
  return _mm256_sin_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sinh(CCTK_REAL4_VEC const x)
{
  return _mm256_sinh_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4tan(CCTK_REAL4_VEC const x)
{
  return _mm256_tan_ps(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4tanh(CCTK_REAL4_VEC const x)
{
  return _mm256_tanh_ps(x);
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
  return R2I(_mm256_xor_ps(I2R(k4ltrue), I2R(x)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4land(CCTK_BOOLEAN4_VEC const x, CCTK_BOOLEAN4_VEC const y)
{
  return R2I(_mm256_and_ps(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4lor(CCTK_BOOLEAN4_VEC const x, CCTK_BOOLEAN4_VEC const y)
{
  return R2I(_mm256_or_ps(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4lxor(CCTK_BOOLEAN4_VEC const x, CCTK_BOOLEAN4_VEC const y)
{
  return R2I(_mm256_xor_ps(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4ifthen(CCTK_BOOLEAN4_VEC const x,
                        CCTK_REAL4_VEC const y,
                        CCTK_REAL4_VEC const z)
{
  return _mm256_blendv_ps(z, y, I2R(x));
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmpeq(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm256_cmp_ps(x, y, _CMP_EQ_OQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmpne(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm256_cmp_ps(x, y, _CMP_NEQ_UQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmpgt(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm256_cmp_ps(x, y, _CMP_GT_OQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmpge(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm256_cmp_ps(x, y, _CMP_GE_OQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmplt(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm256_cmp_ps(x, y, _CMP_LT_OQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN4_VEC k4cmple(CCTK_REAL4_VEC const x, CCTK_REAL4_VEC const y)
{
  return R2I(_mm256_cmp_ps(x, y, _CMP_LE_OQ));
}



static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL4_VEC k4sgn(CCTK_REAL4_VEC const x)
{
  CCTK_BOOLEAN4_VEC const iszero    = k4cmpeq(x, vec4_set1(0.0));
  CCTK_REAL4_VEC    const sign      = _mm256_and_ps(I2R(k4sign), x);
  CCTK_REAL4_VEC    const signedone = _mm256_or_ps(sign, vec4_set1(1.0));
  return k4ifthen(iszero, vec4_set1(0.0), signedone);
}



#undef I2R
#undef R2I
