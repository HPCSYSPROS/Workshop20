// -*-C++-*-
// Vectorise using Intel's or AMD's AVX

// Use the type __m256d directly, without introducing a wrapper class



#include <cstdlib>
#include <cstring>



#include <immintrin.h>
#ifdef __FMA4__
#  include <x86intrin.h>
#endif



#ifdef __FMA4__
#  define vec8_architecture_FMA4 "+FMA4"
#else
#  define vec8_architecture_FMA4 ""
#endif
#define vec8_architecture "AVX" vec8_architecture_FMA4 " (64-bit precision)"



// Vector type corresponding to CCTK_REAL
// Note: some boolean masks (e.g. ~0) correspond to nan when
// interpreted as floating point number. gcc 4.8 is clever enough to
// optimize away such constants with fast-math. We therefore need to
// handle this constant as integer number.
typedef __m256d CCTK_REAL8_VEC;
typedef __m256i CCTK_INTEGER8_VEC;
typedef __m256i CCTK_BOOLEAN8_VEC;

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL8_VEC_SIZE 4

vec_static_assert(sizeof(CCTK_REAL8_VEC) ==
                  sizeof(CCTK_REAL8) * CCTK_REAL8_VEC_SIZE);

// Integer and boolean types corresponding to this real type
typedef CCTK_INT8 CCTK_INTEGER8;
typedef CCTK_INT8 CCTK_BOOLEAN8;



// These macros are undefined at the end of this file -- use them only
// within functions, not within macros that are exported
#define I2R(x) _mm256_castsi256_pd(x)
#define R2I(x) _mm256_castpd_si256(x)



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
  return _mm256_set1_pd(a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_INTEGER8_VEC vec8_set1i(CCTK_INT8 const a)
{
  return _mm256_set1_epi64x(a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_set(CCTK_REAL8 const a,
                        CCTK_REAL8 const b,
                        CCTK_REAL8 const c,
                        CCTK_REAL8 const d)
{
  return _mm256_set_pd(d,c,b,a); // note reversed arguments
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_INTEGER8_VEC vec8_seti(CCTK_INTEGER8 const a,
                            CCTK_INTEGER8 const b,
                            CCTK_INTEGER8 const c,
                            CCTK_INTEGER8 const d)
{
  return _mm256_set_epi64x(d,c,b,a); // note reversed arguments
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
  return _mm256_load_pd(&p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_loadu(CCTK_REAL8 const& p)
{
  return _mm256_loadu_pd(&p);
}
#if VECTORISE_ALWAYS_USE_ALIGNED_LOADS
#  error "VECTORISE_ALWAYS_USE_ALIGNED_LOADS not yet supported"
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
  return off % CCTK_REAL8_VEC_SIZE == 0 ? vec8_load(p) : vec8_loadu(p);
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
    off2 % CCTK_REAL8_VEC_SIZE != 0 or off3 % CCTK_REAL8_VEC_SIZE != 0 ?
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
  _mm256_store_pd(&p, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_storeu(CCTK_REAL8& p, CCTK_REAL8_VEC const x)
{
  _mm256_storeu_pd(&p, x);
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
  _mm256_stream_pd(&p, x);
}
#endif

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare(i, imin,imax)                        \
  bool v8stp_all;                                                       \
  __m256i v8stp_mask;                                                   \
  vec8_store_partial_prepare_(v8stp_all, v8stp_mask, i, imin, imax)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_partial_prepare_(bool& all, __m256i& mask,
                                 std::ptrdiff_t const i,
                                 std::ptrdiff_t const imin,
                                 std::ptrdiff_t const imax)
{
  all = i>=imin and i+CCTK_REAL8_VEC_SIZE-1<imax;
  
  if (not CCTK_BUILTIN_EXPECT(all, true)) {
    /* __m256i const v8stp_mask = */
    /*   _mm256_andnot_pd(_mm256_add_epi64(_mm256_set1_epi64x(i-imin), */
    /*                                     vec_index), */
    /*                    _mm256_add_epi64(_mm256_set1_epi64x(i-imax), */
    /*                                     vec_index)); */
    __m128i const termlo01 =
      _mm_add_epi64(_mm_set1_epi64x(i-imin), _mm_set_epi64x(1, 0));
    __m128i const termup01 =
      _mm_add_epi64(_mm_set1_epi64x(i-imax), _mm_set_epi64x(1, 0));
    __m128i const term01 = _mm_andnot_si128(termlo01, termup01);
    __m128i const termlo23 =
      _mm_add_epi64(_mm_set1_epi64x(i-imin), _mm_set_epi64x(3, 2));
    __m128i const termup23 =
      _mm_add_epi64(_mm_set1_epi64x(i-imax), _mm_set_epi64x(3, 2));
    __m128i const term23 = _mm_andnot_si128(termlo23, termup23);
    mask = _mm256_insertf128_si256(_mm256_castsi128_si256(term01), term23, 1);
  }
}

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare_fixed(i, imin,imax)                  \
  bool v8stp_all;                                                       \
  __m256i v8stp_mask;                                                   \
  vec8_store_partial_prepare_fixed_(v8stp_all, v8stp_mask, i, imin, imax)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_partial_prepare_fixed_(bool& all, __m256i& mask,
                                       std::ptrdiff_t const i,
                                       std::ptrdiff_t const imin,
                                       std::ptrdiff_t const imax)
{
  all = i>=imin and i+CCTK_REAL8_VEC_SIZE-1<imax;
  
  if (not CCTK_BUILTIN_EXPECT(all, true)) {
    mask = vec8_seti(i+0>=imin and i+0<imax ? ~0 : 0,
                     i+1>=imin and i+1<imax ? ~0 : 0,
                     i+2>=imin and i+2<imax ? ~0 : 0,
                     i+3>=imin and i+3<imax ? ~0 : 0);
  }
}

#define vec8_store_nta_partial(p, x)                    \
  vec8_store_nta_partial_(v8stp_all, v8stp_mask, p, x)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_(bool const all, __m256i const mask,
                             CCTK_REAL8& p,
                             CCTK_REAL8_VEC const x)
{
  if (CCTK_BUILTIN_EXPECT(all, true)) {
    vec8_store_nta(p, x);
  } else {
    _mm256_maskstore_pd(&p, mask, x);
  }
}

#define vec8_storeu_partial(p, x)                       \
  vec8_storeu_partial_(v8stp_all, v8stp_mask, p, x)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_storeu_partial_(bool const all, __m256i const mask,
                          CCTK_REAL8& p,
                          CCTK_REAL8_VEC const x)
{
  if (CCTK_BUILTIN_EXPECT(all, true)) {
    vec8_storeu(p, x);
  } else {
    _mm256_maskstore_pd(&p, mask, x);
  }
}

// Store a lower or higher partial vector (aligned and non-temporal);
// the non-temporal hint is probably ignored
// Masks indicating which vector element should be stored:
/*static*/ k8const_t const k8store_lo[5] =
  {
    {{  0,  0,  0,  0, }},
    {{ ~0,  0,  0,  0, }},
    {{ ~0, ~0,  0,  0, }},
    {{ ~0, ~0, ~0,  0, }},
    {{ ~0, ~0, ~0, ~0, }},
  };
/*static*/ k8const_t const k8store_hi[5] =
  {
    {{  0,  0,  0,  0, }},
    {{  0,  0,  0, ~0, }},
    {{  0,  0, ~0, ~0, }},
    {{  0, ~0, ~0, ~0, }},
    {{ ~0, ~0, ~0, ~0, }},
  };
#if !defined(__INTEL_COMPILER) && defined(__GNUC__) && __GNUC__==4 && __GNUC_MINOR__<=4
// gcc 4.4 uses a wrong prototype for _mm256_maskstore_pd
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_lo(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm256_maskstore_pd(&p, I2R(k8store_lo[n].vi), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_hi(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm256_maskstore_pd(&p, I2R(k8store_hi[n].vi), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_mid(CCTK_REAL8& p,
                                CCTK_REAL8_VEC const x,
                                ptrdiff_t const nlo,
                                ptrdiff_t const nhi)
{
  _mm256_maskstore_pd(&p, I2R(k8store_lo[nlo].vi & k8store_hi[nhi].vi), x);
}
#else
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_lo(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm256_maskstore_pd(&p, k8store_lo[n].vi, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_hi(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm256_maskstore_pd(&p, k8store_hi[n].vi, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_mid(CCTK_REAL8& p,
                                CCTK_REAL8_VEC const x,
                                ptrdiff_t const nlo,
                                ptrdiff_t const nhi)
{
  _mm256_maskstore_pd
    (&p,
     R2I(_mm256_and_pd(I2R(k8store_lo[nlo].vi), I2R(k8store_hi[nhi].vi))),
     x);
}
#endif



// Functions and operators

// Operators
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8neg(CCTK_REAL8_VEC const x)
{
  return _mm256_xor_pd(I2R(k8sign), x);
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8add(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm256_add_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sub(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm256_sub_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8mul(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm256_mul_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8div(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm256_div_pd(x, y);
}

// Fused multiply-add, defined as [+-]x*y[+-]z
#ifdef __FMA4__
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8madd(CCTK_REAL8_VEC const x,
                      CCTK_REAL8_VEC const y,
                      CCTK_REAL8_VEC const z)
{
  return _mm256_macc_pd(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8msub(CCTK_REAL8_VEC const x,
                      CCTK_REAL8_VEC const y,
                      CCTK_REAL8_VEC const z)
{
  return _mm256_msub_pd(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8nmadd(CCTK_REAL8_VEC const x,
                       CCTK_REAL8_VEC const y,
                       CCTK_REAL8_VEC const z)
{
  return _mm256_nmsub_pd(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8nmsub(CCTK_REAL8_VEC const x,
                       CCTK_REAL8_VEC const y,
                       CCTK_REAL8_VEC const z)
{
  return _mm256_nmacc_pd(x, y, z);
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
  return _mm256_or_pd(_mm256_and_pd(I2R(k8notsign), x),
                      _mm256_and_pd(I2R(k8sign   ), y));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fabs(CCTK_REAL8_VEC const x)
{
  return _mm256_and_pd(I2R(k8notsign), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmax(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm256_max_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmin(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm256_min_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fnabs(CCTK_REAL8_VEC const x)
{
  return _mm256_or_pd(I2R(k8sign), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8signbit(CCTK_REAL8_VEC const x)
{
  return R2I(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sqrt(CCTK_REAL8_VEC const x)
{
  return _mm256_sqrt_pd(x);
}



// Expensive functions

#define K8REPL(f,x)                             \
  vec8_set(f(vec8_elt(x,0)),                    \
           f(vec8_elt(x,1)),                    \
           f(vec8_elt(x,2)),                    \
           f(vec8_elt(x,3)));
#define K8REPL2S(f,x,a)                         \
  vec8_set(f(vec8_elt(x,0),a),                  \
           f(vec8_elt(x,1),a),                  \
           f(vec8_elt(x,2),a),                  \
           f(vec8_elt(x,3),a));
#define K8REPL2(f,x,y)                          \
  vec8_set(f(vec8_elt(x,0),vec8_elt(y,0)),      \
           f(vec8_elt(x,1),vec8_elt(y,1)),      \
           f(vec8_elt(x,2),vec8_elt(y,2)),      \
           f(vec8_elt(x,3),vec8_elt(y,3)));

#if defined __ICC
// The Intel compiler provides intrinsics for these

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8acos(CCTK_REAL8_VEC const x)
{
  return _mm256_acos_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8acosh(CCTK_REAL8_VEC const x)
{
  return _mm256_acosh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8asin(CCTK_REAL8_VEC const x)
{
  return _mm256_asin_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8asinh(CCTK_REAL8_VEC const x)
{
  return _mm256_asinh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atan(CCTK_REAL8_VEC const x)
{
  return _mm256_atan_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atan2(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm256_atan2_pd(x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atanh(CCTK_REAL8_VEC const x)
{
  return _mm256_atanh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8cos(CCTK_REAL8_VEC const x)
{
  return _mm256_cos_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8cosh(CCTK_REAL8_VEC const x)
{
  return _mm256_cosh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8exp(CCTK_REAL8_VEC const x)
{
  return _mm256_exp_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmod(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
#if __ICC > 1400
  return _mm256_fmod_pd(x, y);
#else
  return K8REPL2(fmod,x,y);
#endif
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8log(CCTK_REAL8_VEC const x)
{
  return _mm256_log_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8pow(CCTK_REAL8_VEC const x, CCTK_REAL8 const a)
{
  return _mm256_pow_pd(x, _mm256_set1_pd(a));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sin(CCTK_REAL8_VEC const x)
{
  return _mm256_sin_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sinh(CCTK_REAL8_VEC const x)
{
  return _mm256_sinh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8tan(CCTK_REAL8_VEC const x)
{
  return _mm256_tan_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8tanh(CCTK_REAL8_VEC const x)
{
  return _mm256_tanh_pd(x);
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
  return R2I(_mm256_xor_pd(I2R(k8ltrue), I2R(x)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8land(CCTK_BOOLEAN8_VEC const x, CCTK_BOOLEAN8_VEC const y)
{
  return R2I(_mm256_and_pd(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8lor(CCTK_BOOLEAN8_VEC const x, CCTK_BOOLEAN8_VEC const y)
{
  return R2I(_mm256_or_pd(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8lxor(CCTK_BOOLEAN8_VEC const x, CCTK_BOOLEAN8_VEC const y)
{
  return R2I(_mm256_xor_pd(I2R(x), I2R(y)));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8ifthen(CCTK_BOOLEAN8_VEC const x,
                        CCTK_REAL8_VEC const y,
                        CCTK_REAL8_VEC const z)
{
  return _mm256_blendv_pd(z, y, I2R(x));
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpeq(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm256_cmp_pd(x, y, _CMP_EQ_OQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpne(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm256_cmp_pd(x, y, _CMP_NEQ_UQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpgt(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm256_cmp_pd(x, y, _CMP_GT_OQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpge(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm256_cmp_pd(x, y, _CMP_GE_OQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmplt(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm256_cmp_pd(x, y, _CMP_LT_OQ));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmple(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return R2I(_mm256_cmp_pd(x, y, _CMP_LE_OQ));
}



static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sgn(CCTK_REAL8_VEC const x)
{
  CCTK_BOOLEAN8_VEC const iszero    = k8cmpeq(x, vec8_set1(0.0));
  CCTK_REAL8_VEC    const sign      = _mm256_and_pd(I2R(k8sign), x);
  CCTK_REAL8_VEC    const signedone = _mm256_or_pd(sign, vec8_set1(1.0));
  return k8ifthen(iszero, vec8_set1(0.0), signedone);
}



#undef I2R
#undef R2I
