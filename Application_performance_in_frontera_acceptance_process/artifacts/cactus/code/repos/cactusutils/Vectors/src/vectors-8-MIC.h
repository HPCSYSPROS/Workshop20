// -*-C++-*-
// Vectorise using Intel's MIC

// Use the type __m512d directly, without introducing a wrapper class

// See
// <http://software.intel.com/sites/products/documentation/doclib/stdxe/2013/composerxe/compiler/cpp-lin/index.htm#GUID-B8DF6000-6872-47B4-AA64-D47A38AF21BD.htm>
// and
// <http://software.intel.com/sites/default/files/forum/278102/327364001en.pdf>.



#include <cstdlib>
#include <cstring>

#include <immintrin.h>



#define vec8_architecture "MIC (64-bit precision)"



// Vector type corresponding to CCTK_REAL
typedef __m512d  CCTK_REAL8_VEC;
typedef __m512i  CCTK_INTEGER8_VEC;
typedef __mmask8 CCTK_BOOLEAN8_VEC;

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL8_VEC_SIZE 8

vec_static_assert(sizeof(CCTK_REAL8_VEC) ==
                  sizeof(CCTK_REAL8) * CCTK_REAL8_VEC_SIZE);

// Integer and boolean types corresponding to this real type
typedef CCTK_INT8     CCTK_INTEGER8;
typedef bool          CCTK_BOOLEAN8;



#define k8sign    (vec8i_set1i(  (CCTK_INTEGER8)(1ULL << 63ULL)))
#define k8notsign (vec8i_set1i(~ (CCTK_INTEGER8)(1ULL << 63ULL)))



// Create vectors, extract vector elements

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_set1(CCTK_REAL8 const a)
{
  return _mm512_set1_pd(a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_INTEGER8_VEC vec8i_set1i(CCTK_INT8 const a)
{
  return _mm512_set1_epi64(a);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_set(CCTK_REAL8 const a0,
                        CCTK_REAL8 const a1,
                        CCTK_REAL8 const a2,
                        CCTK_REAL8 const a3,
                        CCTK_REAL8 const a4,
                        CCTK_REAL8 const a5,
                        CCTK_REAL8 const a6,
                        CCTK_REAL8 const a7)
{
  return _mm512_set_pd(a7,a6,a5,a4,a3,a2,a1,a0); // note reversed arguments
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8 vec8_elt(CCTK_REAL8_VEC const x, std::ptrdiff_t const d)
{
  CCTK_REAL8 e;
  std::memcpy(&e, &((char const*)&x)[d*sizeof e], sizeof e);
  return e;
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8 vec8_eltb(CCTK_BOOLEAN8_VEC const x, std::ptrdiff_t const d)
{
  return _mm512_mask2int(x) & (1 << d);
}



// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_load(CCTK_REAL8 const& p)
{
  return _mm512_load_pd(&p);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC vec8_loadu(CCTK_REAL8 const& p)
{
  CCTK_REAL8_VEC x = _mm512_undefined_pd();
  x = _mm512_loadunpacklo_pd(x, &p);
  x = _mm512_loadunpackhi_pd(x, &p+8);
  return x;
}
#if VECTORISE_ALWAYS_USE_ALIGNED_LOADS
#  error "VECTORISE_ALWAYS_USE_ALIGNED_LOADS is not yet supported"
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
  _mm512_store_pd(&p, x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_storeu(CCTK_REAL8& p, CCTK_REAL8_VEC const x)
{
  // TODO: Intel erratum suggests that hi should come before lo
  _mm512_packstorehi_pd(&p+8, x);
  _mm512_packstorelo_pd(&p  , x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta(CCTK_REAL8& p, CCTK_REAL8_VEC const x)
{
#if VECTORISE_STREAMING_STORES
  // non-temporal hint:
  // _mm512_extstore_pd(&p, x, _MM_DOWNCONV_PD_NONE, _MM_HINT_NT);
  // no-read hint:
  _mm512_storenr_pd(&p, x);
  _mm_clevict(&p, _MM_HINT_T1);
  // no-read hint, not globally ordered (requires fence?):
  // _mm512_storenrngo_pd(&p, x);
  // _mm_clevict(&p, _MM_HINT_T1);

#else
  _mm512_store_pd(&p, x);
#endif
}

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare(i, imin,imax)                        \
  __mmask8 v8stp_mask;                                                  \
  vec8_store_partial_prepare_(v8stp_mask, i, imin, imax)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_partial_prepare_(__mmask8& mask,
                                 std::ptrdiff_t const i,
                                 std::ptrdiff_t const imin,
                                 std::ptrdiff_t const imax)
{
  unsigned char m = 255;
  if (i < imin) {
    /* clear lower imin-i bits */
    m &= 255 << (imin-i);
  }
  if (i+CCTK_REAL8_VEC_SIZE > imax) {
    /* clear upper i+CCTK_REAL8_VEC_SIZE-imax bits */
    m &= 255 >> (i+CCTK_REAL8_VEC_SIZE-imax);
  }
  mask = _mm512_int2mask(m);
}

#define vec8_store_nta_partial(p, x)            \
  vec8_store_nta_partial_(v8stp_mask, p, x)
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_(__mmask8 const mask,
                             CCTK_REAL8& p,
                             CCTK_REAL8_VEC const x)
{
  // TODO: use vec8_store_nta(p, x) if all=true?
  _mm512_mask_store_pd(&p, mask, x);
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_lo(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm512_mask_store_pd(&p, _mm512_int2mask(255 >> (8-n)), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_hi(CCTK_REAL8& p,
                               CCTK_REAL8_VEC const x,
                               ptrdiff_t const n)
{
  _mm512_mask_store_pd(&p, _mm512_int2mask(255 << (8-n)), x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
void vec8_store_nta_partial_mid(CCTK_REAL8& p,
                                CCTK_REAL8_VEC const x,
                                ptrdiff_t const nlo,
                                ptrdiff_t const nhi)
{
  _mm512_mask_store_pd
    (&p, _mm512_int2mask((255 >> (8-nlo)) & (255 << (8-nhi))), x);
}



// Functions and operators

// Operators
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8neg(CCTK_REAL8_VEC const x)
{
  // Could also multiply by -1
  // Could also invert sign bit
  return _mm512_sub_pd(_mm512_set1_pd(0.0), x);
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8add(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_add_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sub(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_sub_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8mul(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_mul_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8div(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_div_pd(x, y);
}

// Fused multiply-add, defined as [+-]x*y[+-]z
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8madd(CCTK_REAL8_VEC const x,
                      CCTK_REAL8_VEC const y,
                      CCTK_REAL8_VEC const z)
{
  return _mm512_fmadd_pd(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8msub(CCTK_REAL8_VEC const x,
                      CCTK_REAL8_VEC const y,
                      CCTK_REAL8_VEC const z)
{
  return _mm512_fmsub_pd(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8nmadd(CCTK_REAL8_VEC const x,
                       CCTK_REAL8_VEC const y,
                       CCTK_REAL8_VEC const z)
{
  return _mm512_fnmsub_pd(x, y, z);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8nmsub(CCTK_REAL8_VEC const x,
                       CCTK_REAL8_VEC const y,
                       CCTK_REAL8_VEC const z)
{
  return _mm512_fnmadd_pd(x, y, z);
}

// Cheap functions
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8copysign(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  CCTK_INTEGER8_VEC ix = _mm512_castpd_si512(x);
  CCTK_INTEGER8_VEC iy = _mm512_castpd_si512(y);
  CCTK_INTEGER8_VEC ir = _mm512_or_epi64(_mm512_and_epi64(k8notsign, ix),
                                         _mm512_and_epi64(k8sign   , iy));
  return _mm512_castsi512_pd(ir);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fabs(CCTK_REAL8_VEC const x)
{
  // Could also do k8fmax(x, k8neg(x))
  CCTK_INTEGER8_VEC ix = _mm512_castpd_si512(x);
  CCTK_INTEGER8_VEC ir = _mm512_and_epi64(k8notsign, ix);
  return _mm512_castsi512_pd(ir);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmax(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_gmax_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmin(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_gmin_pd(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fnabs(CCTK_REAL8_VEC const x)
{
  // Could also do k8fmin(x, k8neg(x))
  CCTK_INTEGER8_VEC ix = _mm512_castpd_si512(x);
  CCTK_INTEGER8_VEC ir = _mm512_or_epi64(k8sign, ix);
  return _mm512_castsi512_pd(ir);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sqrt(CCTK_REAL8_VEC const x)
{
  return _mm512_sqrt_pd(x);
}



// Expensive functions

// These implementations are very expensive
#define K8REPL(f,x)                             \
  vec8_set(f(vec8_elt(x,0)),                    \
           f(vec8_elt(x,1)),                    \
           f(vec8_elt(x,2)),                    \
           f(vec8_elt(x,3)),                    \
           f(vec8_elt(x,4)),                    \
           f(vec8_elt(x,5)),                    \
           f(vec8_elt(x,6)),                    \
           f(vec8_elt(x,7)));
#define K8REPL2S(f,x,a)                         \
  vec8_set(f(vec8_elt(x,0),a),                  \
           f(vec8_elt(x,1),a),                  \
           f(vec8_elt(x,2),a),                  \
           f(vec8_elt(x,3),a),                  \
           f(vec8_elt(x,4),a),                  \
           f(vec8_elt(x,5),a),                  \
           f(vec8_elt(x,6),a),                  \
           f(vec8_elt(x,7),a));
#define K8REPL2(f,x,y)                          \
  vec8_set(f(vec8_elt(x,0),vec8_elt(y,0)),      \
           f(vec8_elt(x,1),vec8_elt(y,1)),      \
           f(vec8_elt(x,2),vec8_elt(y,2)),      \
           f(vec8_elt(x,3),vec8_elt(y,3)),      \
           f(vec8_elt(x,4),vec8_elt(y,4)),      \
           f(vec8_elt(x,5),vec8_elt(y,5)),      \
           f(vec8_elt(x,6),vec8_elt(y,6)),      \
           f(vec8_elt(x,7),vec8_elt(y,7)));

#if defined __ICC
// The Intel compiler provides intrinsics for these

// These implementations lead to an ICE with icpc 13.0.1
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8acos(CCTK_REAL8_VEC const x)
{
  return _mm512_acos_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8acosh(CCTK_REAL8_VEC const x)
{
  return _mm512_acosh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8asin(CCTK_REAL8_VEC const x)
{
  return _mm512_asin_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8asinh(CCTK_REAL8_VEC const x)
{
  return _mm512_asinh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atan(CCTK_REAL8_VEC const x)
{
  return _mm512_atan_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atan2(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_atan2_pd(x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8atanh(CCTK_REAL8_VEC const x)
{
  return _mm512_atanh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8cos(CCTK_REAL8_VEC const x)
{
  return _mm512_cos_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8cosh(CCTK_REAL8_VEC const x)
{
  return _mm512_cosh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8exp(CCTK_REAL8_VEC const x)
{
  return _mm512_exp_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8fmod(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return K8REPL2(fmod,x,y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8log(CCTK_REAL8_VEC const x)
{
  return _mm512_log_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8pow(CCTK_REAL8_VEC const x, CCTK_REAL8 const a)
{
  return _mm512_pow_pd(x, _mm512_set1_pd(a));
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sin(CCTK_REAL8_VEC const x)
{
  return _mm512_sin_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sinh(CCTK_REAL8_VEC const x)
{
  return _mm512_sinh_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8tan(CCTK_REAL8_VEC const x)
{
  return _mm512_tan_pd(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8tanh(CCTK_REAL8_VEC const x)
{
  return _mm512_tanh_pd(x);
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



// TODO: try k8lxor(x,x) and k8lxnor(x,x)
#define k8lfalse (_mm512_int2mask( 0))
#define k8ltrue  (_mm512_int2mask(~0))
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8lnot(CCTK_BOOLEAN8_VEC const x)
{
  return _mm512_knot(x);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8land(CCTK_BOOLEAN8_VEC const x, CCTK_BOOLEAN8_VEC const y)
{
  return _mm512_kand(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8lor(CCTK_BOOLEAN8_VEC const x, CCTK_BOOLEAN8_VEC const y)
{
  return _mm512_kor(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8lxor(CCTK_BOOLEAN8_VEC const x, CCTK_BOOLEAN8_VEC const y)
{
  return _mm512_kxor(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8ifthen(CCTK_BOOLEAN8_VEC const x,
                        CCTK_REAL8_VEC const y,
                        CCTK_REAL8_VEC const z)
{
  // This leads to an ICE
  // return _mm512_mask_blend_pd(x, z, y);
#if 0
  // This works:
  return _mm512_mask_mov_pd(z, x, y);
#endif
  // Intel suggests this:
  return x==0 ? z : _mm512_mask_blend_pd(x, z, y);
}

static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpeq(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_cmpeq_pd_mask(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpne(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_cmpneq_pd_mask(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpgt(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_cmpnle_pd_mask(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmpge(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_cmpnlt_pd_mask(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmplt(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_cmplt_pd_mask(x, y);
}
static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_BOOLEAN8_VEC k8cmple(CCTK_REAL8_VEC const x, CCTK_REAL8_VEC const y)
{
  return _mm512_cmple_pd_mask(x, y);
}



static inline CCTK_ATTRIBUTE_ALWAYS_INLINE
CCTK_REAL8_VEC k8sgn(CCTK_REAL8_VEC const x)
{
  CCTK_BOOLEAN8_VEC const iszero = k8cmpeq(x, vec8_set1(0.0));
  CCTK_BOOLEAN8_VEC const isneg  = k8cmplt(x, vec8_set1(0.0));
  return k8ifthen(iszero, vec8_set1(0.0),
                  k8ifthen(isneg, vec8_set1(-1.0), vec8_set1(+1.0)));
}
