// -*-C++-*-
// Vectorise using IBM's Blue Gene/Q QPX (Power)

#ifndef __clang__

// Use the type vector4double directly, without introducing a wrapper class
// Use macros instead of inline functions

// Note: bgxlC_r does not like const declarations, so we need to cast
// them away and/or omit them everywhere

// See <http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp>

#include <assert.h>

#define vec8_assert(x) ((void)0)
// #define vec8_assert(x) assert(x)

// Argonne
#ifdef __cplusplus
#include <builtins.h>
#endif
#include <mass_simd.h>

#define vec8_architecture "QPX"

// Vector type corresponding to CCTK_REAL
// We use a struct to avoid the "const" issue
// #define CCTK_REAL8_VEC vector4double
struct CCTK_REAL8_VEC {
  vector4double v;
  CCTK_REAL8_VEC() {}
  CCTK_REAL8_VEC(CCTK_REAL8_VEC const &x) : v(x.v) {}
  CCTK_REAL8_VEC(vector4double v_) : v(v_) {}
  operator vector4double() const { return v; }
};

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL8_VEC_SIZE 4

// Integer and boolean types corresponding to this real type
#define CCTK_INTEGER8 CCTK_INT8
#define CCTK_BOOLEAN8 CCTK_REAL8
#define CCTK_INTEGER8_VEC CCTK_REAL8_VEC
#define CCTK_BOOLEAN8_VEC CCTK_REAL8_VEC

// Create vectors, extract vector elements

#define vec8_set1(a) (vec_splats(a))
#if 0
#define vec8_set(a, b, c, d)                                                   \
  (vec_insert(                                                                 \
      d,                                                                       \
      vec_insert(c, vec_insert(b, vec_insert(a, CCTK_REAL8_VEC(), 0), 1), 2),  \
      3))
#endif
#define vec8_set(a_, b_, c_, d_)                                               \
  ({                                                                           \
    CCTK_REAL8 const a__ = (a_);                                               \
    CCTK_REAL8 const b__ = (b_);                                               \
    CCTK_REAL8 const c__ = (c_);                                               \
    CCTK_REAL8 const d__ = (d_);                                               \
    CCTK_REAL8 const a = a__;                                                  \
    CCTK_REAL8 const b = b__;                                                  \
    CCTK_REAL8 const c = c__;                                                  \
    CCTK_REAL8 const d = d__;                                                  \
    CCTK_REAL8_VEC x;                                                          \
    ((CCTK_REAL8 *)&x)[0] = a;                                                 \
    ((CCTK_REAL8 *)&x)[1] = b;                                                 \
    ((CCTK_REAL8 *)&x)[2] = c;                                                 \
    ((CCTK_REAL8 *)&x)[3] = d;                                                 \
    x;                                                                         \
  })

#define vec8_b2r(b) ((b) ? +1.0 : -1.0)
#define vec8_setb(a, b, c, d)                                                  \
  (vec8_set(vec8_b2r(a), vec8_b2r(b), vec8_b2r(c), vec8_b2r(d)))

#define vec8_elt0(x) (vec_extract(x, 0))
#define vec8_elt1(x) (vec_extract(x, 1))
#define vec8_elt2(x) (vec_extract(x, 2))
#define vec8_elt3(x) (vec_extract(x, 3))
#define vec8_elt(x, d) (vec_extract(x, d))
#define vec8_elts(x, a, b, c, d)                                               \
  ({                                                                           \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8_VEC x = x__;                                                    \
    a = ((CCTK_REAL8 *)&x)[0];                                                 \
    b = ((CCTK_REAL8 *)&x)[1];                                                 \
    c = ((CCTK_REAL8 *)&x)[2];                                                 \
    d = ((CCTK_REAL8 *)&x)[3];                                                 \
  })

#define vec8_r2b(x) ((x) >= 0.0)
#define vec8_eltb(x, d) (vec8_r2b(vec8_elt(x, d)))

// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
#define vec8_load(p) (vec_lda(0, (CCTK_REAL8 *)&(p)))
#define vec8_loadu(p_)                                                         \
  ({                                                                           \
    CCTK_REAL8 const &p__ = (p_);                                              \
    CCTK_REAL8 &p = *(CCTK_REAL8 *)&p__;                                       \
    vector4double v1, v2, vp;                                                  \
    /* code taken from IBM's compiler documentation */                         \
    v1 = vec_ld(0, &p);   /* load the left part of the vector */               \
    v2 = vec_ld(31, &p);  /* load the right part of the vector */              \
    vp = vec_lvsl(0, &p); /* generate control value */                         \
    vec_perm(v1, v2, vp); /* generate the aligned vector */                    \
  })
#define vec8_loadu_off(off_, p_)                                               \
  ({                                                                           \
    int const off__ = (off_);                                                  \
    CCTK_REAL8 const &p__ = (p_);                                              \
    int off = off__;                                                           \
    CCTK_REAL8 &p = *(CCTK_REAL8 *)&p__;                                       \
    vector4double v1, v2;                                                      \
    off &= CCTK_REAL8_VEC_SIZE - 1;                                            \
    v1 = vec_lda(0, &p - off);                                                 \
    v2 = vec_lda(0, &p - off + CCTK_REAL8_VEC_SIZE);                           \
    off == 1 ? vec_sldw(v1, v2, 1) : off == 2                                  \
                                         ? vec_sldw(v1, v2, 2)                 \
                                         : off == 3 ? vec_sldw(v1, v2, 3)      \
                                                    : (vec8_assert(0), v1);    \
  })

// Load a vector from memory that may or may not be aligned, as
// decided by the offset and the vector size
#if VECTORISE_ALWAYS_USE_UNALIGNED_LOADS
// Implementation: Always use unaligned load
#define vec8_loadu_maybe(off, p) vec8_loadu(p)
#define vec8_loadu_maybe3(off1, off2, off3, p) vec8_loadu(p)
#else
#define vec8_loadu_maybe(off_, p_)                                             \
  ({                                                                           \
    CCTK_REAL8 const &p__ = (p_);                                              \
    int const off__ = (off_);                                                  \
    CCTK_REAL8 const &p = p__;                                                 \
    int const off = off__;                                                     \
    off % CCTK_REAL8_VEC_SIZE == 0 ? vec8_load(p) : vec8_loadu_off(off, p);    \
  })
#if VECTORISE_ALIGNED_ARRAYS
// Assume all array x sizes are multiples of the vector size
#define vec8_loadu_maybe3(off1, off2, off3, p) vec8_loadu_maybe(off1, p)
#else
#define vec8_loadu_maybe3(off1, off2, off3, p_)                                \
  ({                                                                           \
    CCTK_REAL8 const &p__ = (p_);                                              \
    CCTK_REAL8 const &p = p__;                                                 \
    ((off2) % CCTK_REAL8_VEC_SIZE != 0 or (off3) % CCTK_REAL8_VEC_SIZE != 0)   \
        ? vec8_loadu(p)                                                        \
        : vec8_loadu_maybe(off1, p);                                           \
  })
#endif
#endif

// Store a vector to memory (aligned and non-temporal); this stores to
// a reference to a scalar
#define vec8_store(p, x) (vec_sta(x, 0, &(p)))
#define vec8_storeu(p_, x_)                                                    \
  ({                                                                           \
    CCTK_REAL8 &p__ = (p_);                                                    \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8 &p = p__;                                                       \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC v1, v2, v3, vp, m1, m2, m3;                                 \
    /* code taken from IBM's compiler documentation */                         \
    /* generate insert masks */                                                \
    vp = vec_lvsr(0, &p);                                                      \
    m1 = k8lfalse;                                                             \
    m2 = k8ltrue;                                                              \
    m3 = vec_perm(m1, m2, vp);                                                 \
    v3 = vec_perm(x, x, vp);                                                   \
    _Pragma("tm_atomic") {                                                     \
      /* get existing data */                                                  \
      v1 = vec_ld(0, &p);                                                      \
      v2 = vec_ld(31, &p);                                                     \
      /* permute and insert */                                                 \
      v1 = vec_sel(v1, v3, m3);                                                \
      v2 = vec_sel(v3, v2, m3);                                                \
      /* store data back */                                                    \
      vec_st(v1, 0, &p);                                                       \
      vec_st(v2, 31, &p);                                                      \
    }                                                                          \
  })
#define vec8_store_nta(p, x)                                                   \
  (vec_sta(x, 0, &(p))) // this doesn't avoid the cache

#ifdef _OPENMP
#if VECTORISE_ALIGNED_ARRAYS
// Arrays are aligned; wrap-around is not an issue
#define vec8_store_omp
#else
// Need to protect partial stores, as they may wrap around to the
// beginning of the next line in the array
#define vec8_store_omp _Pragma("tm_atomic")
#endif
#else
#define vec8_store_omp
#endif

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare(i, imin_, imax_)                            \
  bool v8stp_all;                                                              \
  CCTK_BOOLEAN8_VEC v8stp_mask;                                                \
  bool v8stp_mask0, v8stp_mask1, v8stp_mask2, v8stp_mask3;                     \
  ({                                                                           \
    ptrdiff_t const imin__ = (imin_);                                          \
    ptrdiff_t const imax__ = (imax_);                                          \
    ptrdiff_t const imin = imin__;                                             \
    ptrdiff_t const imax = imax__;                                             \
                                                                               \
    v8stp_all = i >= imin and i + CCTK_REAL8_VEC_SIZE - 1 < imax;              \
                                                                               \
    if (not CCTK_BUILTIN_EXPECT(v8stp_all, true)) {                            \
      CCTK_INTEGER8_VEC vp_lo, vp_hi;                                          \
      CCTK_BOOLEAN8_VEC mask_lo, mask_hi;                                      \
      /* this is correct but slow */                                           \
      /*                                                                       \
      mask_lo = vec8_setb(i+0>=imin, i+1>=imin, i+2>=imin, i+3>=imin);         \
      mask_hi = vec8_setb(i+0<imax, i+1<imax, i+2<imax, i+3<imax);             \
      */                                                                       \
      /* Note: vec_lvsl(i,p) =  &p[i] / 8 % 4                                  \
         Note: vec_lvsr(i,p) = -&p[i] / 8 % 4                                  \
         /8: 8 bytes per double                                                \
         %4: 4 doubles per vector                                              \
      */                                                                       \
      /* We assume p[i] is aligned */                                          \
      /* Ensure at least one vector element is inside the active region */     \
      vec8_assert(i - imin >= -(CCTK_REAL8_VEC_SIZE - 1));                     \
      vp_lo = vec_lvsl(8 * (i - imin), (CCTK_REAL *)0);                        \
      mask_lo =                                                                \
          (i - imin >= 0 ? k8ltrue : vec_perm(k8lfalse, k8ltrue, vp_lo));      \
      /* Ensure at least one vector element is inside the active region */     \
      vec8_assert(i < imax);                                                   \
      vp_hi = vec_lvsl(8 * (i - imax), (CCTK_REAL *)0);                        \
      mask_hi = (i - imax < -(CCTK_REAL8_VEC_SIZE - 1)                         \
                     ? k8ltrue                                                 \
                     : vec_perm(k8ltrue, k8lfalse, vp_hi));                    \
      v8stp_mask = k8land(mask_lo, mask_hi);                                   \
      v8stp_mask0 = vec8_eltb(v8stp_mask, 0);                                  \
      v8stp_mask1 = vec8_eltb(v8stp_mask, 1);                                  \
      v8stp_mask2 = vec8_eltb(v8stp_mask, 2);                                  \
      v8stp_mask3 = vec8_eltb(v8stp_mask, 3);                                  \
    }                                                                          \
  })
#define vec8_store_nta_partial(p_, x_)                                         \
  ({                                                                           \
    CCTK_REAL8 &p__ = (p_);                                                    \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8 &p = p__;                                                       \
    CCTK_REAL8_VEC x = x__;                                                    \
    if (CCTK_BUILTIN_EXPECT(v8stp_all, true)) {                                \
      vec8_store(p, x);                                                        \
    } else {                                                                   \
      /*                                                                       \
      vec8_store_omp                                                           \
        vec8_store(p, k8ifthen(v8stp_mask, x, vec8_load(p)));                  \
      */                                                                       \
      if (VECTORISE_ALIGNED_ARRAYS) {                                          \
        vec8_store(p, k8ifthen(v8stp_mask, x, vec8_load(p)));                  \
      } else {                                                                 \
        if (v8stp_mask0)                                                       \
          (&p)[0] = vec8_elt0(x);                                              \
        if (v8stp_mask1)                                                       \
          (&p)[1] = vec8_elt1(x);                                              \
        if (v8stp_mask2)                                                       \
          (&p)[2] = vec8_elt2(x);                                              \
        if (v8stp_mask3)                                                       \
          (&p)[3] = vec8_elt3(x);                                              \
      }                                                                        \
    }                                                                          \
  })

// Store a lower or higher partial vector (aligned and non-temporal);
// the non-temporal hint is probably ignored
#define vec8_store_nta_partial_lo(p_, x_, n)                                   \
  ({                                                                           \
    CCTK_REAL8 &p__ = (p_);                                                    \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8 &p = p__;                                                       \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC vp, mask;                                                   \
    /* Ensure at least one and but all vector elements are active */           \
    vec8_assert(n > 0 and n < CCTK_REAL8_VEC_SIZE - 1);                        \
    vp = vec_lvsl(-8 * n, (CCTK_REAL *)0);                                     \
    mask = vec_perm(k8ltrue, k8lfalse, vp);                                    \
    vec8_store_omp vec8_store(p, k8ifthen(mask, x, vec8_load(p)));             \
  })
#define vec8_store_nta_partial_hi(p_, x_, n)                                   \
  ({                                                                           \
    CCTK_REAL8 &p__ = (p_);                                                    \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8 &p = p__;                                                       \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC vp, mask;                                                   \
    /* Ensure at least one but not all vector elements are active */           \
    vec8_assert(n > 0 and n < CCTK_REAL8_VEC_SIZE - 1);                        \
    vp = vec_lvsl(8 * n, (CCTK_REAL *)0);                                      \
    mask = vec_perm(k8lfalse, k8ltrue, vp);                                    \
    vec8_store_omp vec8_store(p, k8ifthen(mask, x, vec8_load(p)));             \
  })
#define vec8_store_nta_partial_mid(p_, x_, nlo, nhi)                           \
  ({                                                                           \
    CCTK_REAL8 &p__ = (p_);                                                    \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8 p = p__;                                                        \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC vp_lo, mask_lo;                                             \
    /* Ensure at least one but not all vector elements are active */           \
    vec8_assert(nlo > 0 and nlo < CCTK_REAL8_VEC_SIZE - 1);                    \
    vp_lo = vec_lvsl(-8 * nlo, (CCTK_REAL *)0);                                \
    mask_lo = vec_perm(k8lfalse, k8ltrue, vp_lo);                              \
    CCTK_REAL8_VEC vp_hi, mask_hi;                                             \
    /* Ensure at least one but not all vector elements are active */           \
    vec8_assert(nhi > 0 and nhi < CCTK_REAL8_VEC_SIZE - 1);                    \
    vp_hi = vec_lvsl(8 * nhi, (CCTK_REAL *)0);                                 \
    mask_hi = vec_perm(k8lfalse, k8ltrue, vp_hi);                              \
    CCTK_REAL8_VEC mask;                                                       \
    mask = vec_and(mask_lo, mask_hi);                                          \
    vec8_store_omp vec8_store(p, k8ifthen(mask, x, vec8_load(p)));             \
  })

// Functions and operators

// Operators
#define k8neg(x) (vec_neg(x))

#define k8add(x, y) (vec_add(x, y))
#define k8sub(x, y) (vec_sub(x, y))
#define k8mul(x, y) (vec_mul(x, y))
#define k8div(x, y) (vec_swdiv_nochk(x, y))

// Fused multiply-add, defined as [+-]x*y[+-]z
#define k8madd(x, y, z) (vec_madd(x, y, z))
#define k8msub(x, y, z) (vec_msub(x, y, z))
#define k8nmadd(x, y, z) (vec_nmadd(x, y, z))
#define k8nmsub(x, y, z) (vec_nmsub(x, y, z))

// Cheap functions
#define k8copysign(x, y) (vec_cpsgn(y, x))
#define k8fabs(x) (vec_abs(x))
#define k8fmax(x_, y_)                                                         \
  ({                                                                           \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8_VEC y__ = (y_);                                                 \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC y = y__;                                                    \
    k8ifthen(k8cmplt(x, y), y, x);                                             \
  })
#define k8fmin(x_, y_)                                                         \
  ({                                                                           \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8_VEC y__ = (y_);                                                 \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC y = y__;                                                    \
    k8ifthen(k8cmpgt(x, y), y, x);                                             \
  })
#define k8fnabs(x) (vec_nabs(x))
#define k8sgn(x_)                                                              \
  ({                                                                           \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC one, zero, iszero;                                          \
    one = k8ltrue;                                                             \
    zero = k8sub(one, one);                                                    \
    iszero = k8cmpeq(x, zero);                                                 \
    k8ifthen(iszero, zero, k8copysign(one, x));                                \
  })
#define k8sqrt(x) (vec_swsqrt_nochk(x))

// Expensive functions

#define K8REPL(f, x_)                                                          \
  ({                                                                           \
    CCTK_REAL8_VEC const x__ = (x_);                                           \
    CCTK_REAL8_VEC const x = x__;                                              \
    vec8_set(f(vec8_elt0(x)), f(vec8_elt1(x)), f(vec8_elt2(x)),                \
             f(vec8_elt3(x)));                                                 \
  })
#define K8REPL2S(f, x_, a_)                                                    \
  ({                                                                           \
    CCTK_REAL8_VEC const x__ = (x_);                                           \
    CCTK_REAL8 const a__ = (a_);                                               \
    CCTK_REAL8_VEC const x = x__;                                              \
    CCTK_REAL8 const a = a__;                                                  \
    vec8_set(f(vec8_elt0(x), a), f(vec8_elt1(x), a), f(vec8_elt2(x), a),       \
             f(vec8_elt3(x), a));                                              \
  })
#define K8REPL2(f, x_, y_)                                                     \
  ({                                                                           \
    CCTK_REAL8_VEC const x__ = (x_);                                           \
    CCTK_REAL8_VEC const y__ = (y_);                                           \
    CCTK_REAL8_VEC const x = x__;                                              \
    CCTK_REAL8_VEC const y = y__;                                              \
    vec8_set(f(vec8_elt0(x), vec8_elt0(y)), f(vec8_elt1(x), vec8_elt1(y)),     \
             f(vec8_elt2(x), vec8_elt2(y)), f(vec8_elt3(x), vec8_elt3(y)));    \
  })

#define k8acos(x) acosd4(x)
#define k8acosh(x) acoshd4(x)
#define k8asin(x) asind4(x)
#define k8asinh(x) asinhd4(x)
#define k8atan(x) atand4(x)
#define k8atan2(x, y) atan2d4(x, y)
#define k8atanh(x) atanhd4(x)
#define k8cos(x) cosd4(x)
#define k8cosh(x) coshd4(x)
#define k8exp(x) expd4(x)
// #define k8fmod(x,y)  fmodd4(x,y)
#define k8fmod(x, y) K8REPL2(fmod, x, y)
#define k8log(x) logd4(x)
#define k8pow(x, a) powd4(x, vec_set1(a))
#define k8sin(x) sind4(x)
#define k8sinh(x) sinhd4(x)
#define k8tan(x) tand4(x)
#define k8tanh(x) tanhd4(x)

// canonical true is +1.0, canonical false is -1.0
// >=0 is true, -0 is true, nan is false
#define k8lfalse                                                               \
  ({                                                                           \
    CCTK_REAL8_VEC dummy;                                                      \
    vec_logical(dummy, dummy, 0x0);                                            \
  })
#define k8ltrue                                                                \
  ({                                                                           \
    CCTK_REAL8_VEC dummy;                                                      \
    vec_logical(dummy, dummy, 0xf);                                            \
  })
#define k8lnot(x) (vec_not(x))
#define k8land(x, y) (vec_and(x, y))
#define k8lor(x, y) (vec_or(x, y))
#define k8lxor(x, y) (vec_xor(x, y))
#define k8ifthen(x, y, z) (vec_sel(z, y, x))

#define k8cmpeq(x, y) (vec_cmpeq(x, y))
#define k8cmpne(x, y) (k8lnot(k8cmpeq(x, y)))
#define k8cmpgt(x, y) (vec_cmpgt(x, y))
#define k8cmpge(x, y) (k8lnot(k8cmplt(x, y)))
#define k8cmplt(x, y) (vec_cmplt(x, y))
#define k8cmple(x, y) (k8lnot(k8cmpgt(x, y)))

#else // #ifdef __clang__

// Use the type vector4double directly, without introducing a wrapper class
// Use macros instead of inline functions

// See <http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp>

#include <assert.h>

#define vec8_assert(x) ((void)0)
// #define vec8_assert(x) assert(x)

#include <qpxmath.h>

#define vec8_architecture "QPX"

// Vector type corresponding to CCTK_REAL
typedef vector4double CCTK_REAL8_VEC;

// Number of vector elements in a CCTK_REAL_VEC
#define CCTK_REAL8_VEC_SIZE 4

// Integer and boolean types corresponding to this real type
typedef CCTK_INT8 CCTK_INTEGER8;
typedef CCTK_REAL8 CCTK_BOOLEAN8;
typedef CCTK_REAL8_VEC CCTK_INTEGER8_VEC;
typedef CCTK_REAL8_VEC CCTK_BOOLEAN8_VEC;

// canonical true is +1.0, canonical false is -1.0
// >=0 is true, -0 is true, nan is false
inline CCTK_BOOLEAN8_VEC k8lfalse_() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
  CCTK_BOOLEAN8_VEC dummy;
  return vec_logical(dummy, dummy, 0x0);
#pragma clang diagnostic pop
}
inline CCTK_BOOLEAN8_VEC k8ltrue_() {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
  CCTK_BOOLEAN8_VEC dummy;
  return vec_logical(dummy, dummy, 0xf);
#pragma clang diagnostic pop
}
#define k8lfalse k8lfalse_()
#define k8ltrue k8ltrue_()

// Create vectors, extract vector elements

inline CCTK_REAL8_VEC vec8_set1(CCTK_REAL8 a) { return vec_splats(a); }
inline CCTK_REAL8_VEC vec8_set(CCTK_REAL8 a, CCTK_REAL8 b, CCTK_REAL8 c,
                               CCTK_REAL8 d) {
  return (CCTK_REAL8_VEC){a, b, c, d};
}

inline CCTK_REAL8 vec8_b2r(bool b) { return b ? +1.0 : -1.0; }
inline CCTK_BOOLEAN8_VEC vec8_setb(bool a, bool b, bool c, bool d) {
  return vec8_set(vec8_b2r(a), vec8_b2r(b), vec8_b2r(c), vec8_b2r(d));
}

inline CCTK_REAL8 vec8_elt0(CCTK_REAL8_VEC x) { return vec_extract(x, 0); }
inline CCTK_REAL8 vec8_elt1(CCTK_REAL8_VEC x) { return vec_extract(x, 1); }
inline CCTK_REAL8 vec8_elt2(CCTK_REAL8_VEC x) { return vec_extract(x, 2); }
inline CCTK_REAL8 vec8_elt3(CCTK_REAL8_VEC x) { return vec_extract(x, 3); }
inline CCTK_REAL8 vec8_elt(CCTK_REAL8_VEC x, int d) {
  return vec_extract(x, d);
}
/*
inline void vec8_elts(CCTK_REAL8_VEC x, CCTK_REAL8 &a, CCTK_REAL8 &b,
                      CCTK_REAL8 &c, CCTK_REAL8 &d) {
  a = x[0];
  b = x[1];
  c = x[2];
  d = x[3];
}
*/

inline bool vec8_r2b(CCTK_REAL8 x) { return x >= 0.0; }
inline bool vec8_eltb(CCTK_BOOLEAN8_VEC x, int d) {
  return vec8_r2b(vec8_elt(x, d));
}

// Load and store vectors

// Load a vector from memory (aligned and unaligned); this loads from
// a reference to a scalar
inline CCTK_REAL8_VEC vec8_load(const CCTK_REAL8 &p) {
  return vec_lda(0, (CCTK_REAL8 *)&p);
}
inline CCTK_REAL8_VEC vec8_loadu(const CCTK_REAL8 &p) {
  vector4double v1, v2, vp; // code taken from IBM's compiler documentation
  v1 = vec_ld(0, (CCTK_REAL8 *)&p);   // load the left part of the vector
  v2 = vec_ld(31, (CCTK_REAL8 *)&p);  // load the right part of the vector
  vp = vec_lvsl(0, (CCTK_REAL8 *)&p); // generate control value
  return vec_perm(v1, v2, vp);        // generate the aligned vector
}
inline CCTK_REAL8_VEC vec8_loadu_off(ptrdiff_t off, const CCTK_REAL8 &p) {
  vector4double v1, v2;
  off &= CCTK_REAL8_VEC_SIZE - 1;
  v1 = vec_lda(0, (CCTK_REAL8 *)&p - off);
  v2 = vec_lda(0, (CCTK_REAL8 *)&p - off + CCTK_REAL8_VEC_SIZE);
  return off == 1 ? vec_sldw(v1, v2, 1) : off == 2
                                              ? vec_sldw(v1, v2, 2)
                                              : off == 3 ? vec_sldw(v1, v2, 3)
                                                         : (vec8_assert(0), v1);
}

// Load a vector from memory that may or may not be aligned, as
// decided by the offset and the vector size
#if VECTORISE_ALWAYS_USE_UNALIGNED_LOADS
// Implementation: Always use unaligned load
inline CCTK_REAL8_VEC vec8_loadu_maybe(ptrdiff_t off, const CCTK_REAL8 &p) {
  return vec8_loadu(p);
}
inline CCTK_REAL8_VEC vec8_loadu_maybe3(ptrdiff_t off1, ptrdiff_t off2,
                                        ptrdiff_t off3, const CCTK_REAL8 &p) {
  return vec8_loadu(p);
}
#else
inline CCTK_REAL8_VEC vec8_loadu_maybe(ptrdiff_t off, const CCTK_REAL8 &p) {
  return off % CCTK_REAL8_VEC_SIZE == 0 ? vec8_load(p) : vec8_loadu_off(off, p);
}
#if VECTORISE_ALIGNED_ARRAYS
// Assume all array x sizes are multiples of the vector size
inline CCTK_REAL8_VEC vec8_loadu_maybe3(ptrdiff_t off1, ptrdiff_t off2,
                                        ptrdiff_t off3, const CCTK_REAL8 &p) {
  return vec8_loadu_maybe(off1, p);
}
#else
inline CCTK_REAL8_VEC vec8_loadu_maybe3(ptrdiff_t off1, ptrdiff_t off2,
                                        ptrdiff_t off3, const CCTK_REAL8 &p) {
  return off2 % CCTK_REAL8_VEC_SIZE != 0 or off3 % CCTK_REAL8_VEC_SIZE != 0
             ? vec8_loadu(p)
             : vec8_loadu_maybe(off1, p);
}
#endif
#endif

// Store a vector to memory (aligned and non-temporal); this stores to
// a reference to a scalar
inline void vec8_store(CCTK_REAL8 &p, CCTK_REAL8_VEC x) { vec_sta(x, 0, &p); }
inline void vec8_storeu(CCTK_REAL8 &p, CCTK_REAL8_VEC x) {
  CCTK_REAL8_VEC v1, v2, v3, vp, m1, m2, m3;
  // code taken from IBM's compiler documentation
  // generate insert masks
  vp = vec_lvsr(0, &p);
  m1 = k8lfalse;
  m2 = k8ltrue;
  m3 = vec_perm(m1, m2, vp);
  v3 = vec_perm(x, x, vp);
#pragma tm_atomic
  {
    // get existing data
    v1 = vec_ld(0, &p);
    v2 = vec_ld(31, &p);
    // permute and insert
    v1 = vec_sel(v1, v3, m3);
    v2 = vec_sel(v3, v2, m3);
    // store data back
    vec_st(v1, 0, &p);
    vec_st(v2, 31, &p);
  }
}
inline void vec8_store_nta(CCTK_REAL8 &p, CCTK_REAL8_VEC x) {
  // this doesn't avoid the cache
  return vec_sta(x, 0, &p);
}

#ifdef _OPENMP
#if VECTORISE_ALIGNED_ARRAYS
// Arrays are aligned; wrap-around is not an issue
#define vec8_store_omp
#else
// Need to protect partial stores, as they may wrap around to the
// beginning of the next line in the array
#define vec8_store_omp _Pragma("tm_atomic")
#endif
#else
#define vec8_store_omp
#endif

// Store a partial vector (aligned and non-temporal)
#define vec8_store_partial_prepare(i, imin_, imax_)                            \
  bool v8stp_all;                                                              \
  CCTK_BOOLEAN8_VEC v8stp_mask;                                                \
  bool v8stp_mask0, v8stp_mask1, v8stp_mask2, v8stp_mask3;                     \
  ({                                                                           \
    ptrdiff_t const imin__ = (imin_);                                          \
    ptrdiff_t const imax__ = (imax_);                                          \
    ptrdiff_t const imin = imin__;                                             \
    ptrdiff_t const imax = imax__;                                             \
                                                                               \
    v8stp_all = i >= imin and i + CCTK_REAL8_VEC_SIZE - 1 < imax;              \
                                                                               \
    if (not CCTK_BUILTIN_EXPECT(v8stp_all, true)) {                            \
      CCTK_INTEGER8_VEC vp_lo, vp_hi;                                          \
      CCTK_BOOLEAN8_VEC mask_lo, mask_hi;                                      \
      /* this is correct but slow */                                           \
      /*                                                                       \
      mask_lo = vec8_setb(i+0>=imin, i+1>=imin, i+2>=imin, i+3>=imin);         \
      mask_hi = vec8_setb(i+0<imax, i+1<imax, i+2<imax, i+3<imax);             \
      */                                                                       \
      /* Note: vec_lvsl(i,p) =  &p[i] / 8 % 4                                  \
         Note: vec_lvsr(i,p) = -&p[i] / 8 % 4                                  \
         /8: 8 bytes per double                                                \
         %4: 4 doubles per vector                                              \
      */                                                                       \
      /* We assume p[i] is aligned */                                          \
      /* Ensure at least one vector element is inside the active region */     \
      vec8_assert(i - imin >= -(CCTK_REAL8_VEC_SIZE - 1));                     \
      vp_lo = vec_lvsl(8 * (i - imin), (CCTK_REAL *)0);                        \
      mask_lo =                                                                \
          (i - imin >= 0 ? k8ltrue : vec_perm(k8lfalse, k8ltrue, vp_lo));      \
      /* Ensure at least one vector element is inside the active region */     \
      vec8_assert(i < imax);                                                   \
      vp_hi = vec_lvsl(8 * (i - imax), (CCTK_REAL *)0);                        \
      mask_hi = (i - imax < -(CCTK_REAL8_VEC_SIZE - 1)                         \
                     ? k8ltrue                                                 \
                     : vec_perm(k8ltrue, k8lfalse, vp_hi));                    \
      v8stp_mask = k8land(mask_lo, mask_hi);                                   \
      v8stp_mask0 = vec8_eltb(v8stp_mask, 0);                                  \
      v8stp_mask1 = vec8_eltb(v8stp_mask, 1);                                  \
      v8stp_mask2 = vec8_eltb(v8stp_mask, 2);                                  \
      v8stp_mask3 = vec8_eltb(v8stp_mask, 3);                                  \
    }                                                                          \
  })
#define vec8_store_nta_partial(p_, x_)                                         \
  ({                                                                           \
    CCTK_REAL8 &p__ = (p_);                                                    \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8 &p = p__;                                                       \
    CCTK_REAL8_VEC x = x__;                                                    \
    if (CCTK_BUILTIN_EXPECT(v8stp_all, true)) {                                \
      vec8_store(p, x);                                                        \
    } else {                                                                   \
      /*                                                                       \
      vec8_store_omp                                                           \
        vec8_store(p, k8ifthen(v8stp_mask, x, vec8_load(p)));                  \
      */                                                                       \
      if (VECTORISE_ALIGNED_ARRAYS) {                                          \
        vec8_store(p, k8ifthen(v8stp_mask, x, vec8_load(p)));                  \
      } else {                                                                 \
        if (v8stp_mask0)                                                       \
          (&p)[0] = vec8_elt0(x);                                              \
        if (v8stp_mask1)                                                       \
          (&p)[1] = vec8_elt1(x);                                              \
        if (v8stp_mask2)                                                       \
          (&p)[2] = vec8_elt2(x);                                              \
        if (v8stp_mask3)                                                       \
          (&p)[3] = vec8_elt3(x);                                              \
      }                                                                        \
    }                                                                          \
  })

// Store a lower or higher partial vector (aligned and non-temporal);
// the non-temporal hint is probably ignored
#define vec8_store_nta_partial_lo(p_, x_, n)                                   \
  ({                                                                           \
    CCTK_REAL8 &p__ = (p_);                                                    \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8 &p = p__;                                                       \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC vp, mask;                                                   \
    /* Ensure at least one and but all vector elements are active */           \
    vec8_assert(n > 0 and n < CCTK_REAL8_VEC_SIZE - 1);                        \
    vp = vec_lvsl(-8 * n, (CCTK_REAL *)0);                                     \
    mask = vec_perm(k8ltrue, k8lfalse, vp);                                    \
    vec8_store_omp vec8_store(p, k8ifthen(mask, x, vec8_load(p)));             \
  })
#define vec8_store_nta_partial_hi(p_, x_, n)                                   \
  ({                                                                           \
    CCTK_REAL8 &p__ = (p_);                                                    \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8 &p = p__;                                                       \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC vp, mask;                                                   \
    /* Ensure at least one but not all vector elements are active */           \
    vec8_assert(n > 0 and n < CCTK_REAL8_VEC_SIZE - 1);                        \
    vp = vec_lvsl(8 * n, (CCTK_REAL *)0);                                      \
    mask = vec_perm(k8lfalse, k8ltrue, vp);                                    \
    vec8_store_omp vec8_store(p, k8ifthen(mask, x, vec8_load(p)));             \
  })
#define vec8_store_nta_partial_mid(p_, x_, nlo, nhi)                           \
  ({                                                                           \
    CCTK_REAL8 &p__ = (p_);                                                    \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8 p = p__;                                                        \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC vp_lo, mask_lo;                                             \
    /* Ensure at least one but not all vector elements are active */           \
    vec8_assert(nlo > 0 and nlo < CCTK_REAL8_VEC_SIZE - 1);                    \
    vp_lo = vec_lvsl(-8 * nlo, (CCTK_REAL *)0);                                \
    mask_lo = vec_perm(k8lfalse, k8ltrue, vp_lo);                              \
    CCTK_REAL8_VEC vp_hi, mask_hi;                                             \
    /* Ensure at least one but not all vector elements are active */           \
    vec8_assert(nhi > 0 and nhi < CCTK_REAL8_VEC_SIZE - 1);                    \
    vp_hi = vec_lvsl(8 * nhi, (CCTK_REAL *)0);                                 \
    mask_hi = vec_perm(k8lfalse, k8ltrue, vp_hi);                              \
    CCTK_REAL8_VEC mask;                                                       \
    mask = vec_and(mask_lo, mask_hi);                                          \
    vec8_store_omp vec8_store(p, k8ifthen(mask, x, vec8_load(p)));             \
  })

// Functions and operators

// Operators
#define k8neg(x) (vec_neg(x))

#define k8add(x, y) (vec_add(x, y))
#define k8sub(x, y) (vec_sub(x, y))
#define k8mul(x, y) (vec_mul(x, y))
#define k8div(x, y) (vec_swdiv_nochk(x, y))

// Fused multiply-add, defined as [+-]x*y[+-]z
#define k8madd(x, y, z) (vec_madd(x, y, z))
#define k8msub(x, y, z) (vec_msub(x, y, z))
#define k8nmadd(x, y, z) (vec_nmadd(x, y, z))
#define k8nmsub(x, y, z) (vec_nmsub(x, y, z))

// Cheap functions
#define k8copysign(x, y) (vec_cpsgn(y, x))
#define k8fabs(x) (vec_abs(x))
#define k8fmax(x_, y_)                                                         \
  ({                                                                           \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8_VEC y__ = (y_);                                                 \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC y = y__;                                                    \
    k8ifthen(k8cmplt(x, y), y, x);                                             \
  })
#define k8fmin(x_, y_)                                                         \
  ({                                                                           \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8_VEC y__ = (y_);                                                 \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC y = y__;                                                    \
    k8ifthen(k8cmpgt(x, y), y, x);                                             \
  })
#define k8fnabs(x) (vec_nabs(x))
#define k8sgn(x_)                                                              \
  ({                                                                           \
    CCTK_REAL8_VEC x__ = (x_);                                                 \
    CCTK_REAL8_VEC x = x__;                                                    \
    CCTK_REAL8_VEC one, zero, iszero;                                          \
    one = k8ltrue;                                                             \
    zero = k8sub(one, one);                                                    \
    iszero = k8cmpeq(x, zero);                                                 \
    k8ifthen(iszero, zero, k8copysign(one, x));                                \
  })
#define k8sqrt(x) (vec_swsqrt_nochk(x))

// Expensive functions

#define K8REPL(f, x_)                                                          \
  ({                                                                           \
    CCTK_REAL8_VEC const x__ = (x_);                                           \
    CCTK_REAL8_VEC const x = x__;                                              \
    vec8_set(f(vec8_elt0(x)), f(vec8_elt1(x)), f(vec8_elt2(x)),                \
             f(vec8_elt3(x)));                                                 \
  })
#define K8REPL2S(f, x_, a_)                                                    \
  ({                                                                           \
    CCTK_REAL8_VEC const x__ = (x_);                                           \
    CCTK_REAL8 const a__ = (a_);                                               \
    CCTK_REAL8_VEC const x = x__;                                              \
    CCTK_REAL8 const a = a__;                                                  \
    vec8_set(f(vec8_elt0(x), a), f(vec8_elt1(x), a), f(vec8_elt2(x), a),       \
             f(vec8_elt3(x), a));                                              \
  })
#define K8REPL2(f, x_, y_)                                                     \
  ({                                                                           \
    CCTK_REAL8_VEC const x__ = (x_);                                           \
    CCTK_REAL8_VEC const y__ = (y_);                                           \
    CCTK_REAL8_VEC const x = x__;                                              \
    CCTK_REAL8_VEC const y = y__;                                              \
    vec8_set(f(vec8_elt0(x), vec8_elt0(y)), f(vec8_elt1(x), vec8_elt1(y)),     \
             f(vec8_elt2(x), vec8_elt2(y)), f(vec8_elt3(x), vec8_elt3(y)));    \
  })

#define k8acos(x) xacos(x)
#define k8acosh(x) xacosh(x)
#define k8asin(x) xasin(x)
#define k8asinh(x) xasinh(x)
#define k8atan(x) xatan(x)
#define k8atan2(x, y) xatan2(x, y)
#define k8atanh(x) xatanh(x)
#define k8cos(x) xcos(x)
#define k8cosh(x) xcosh(x)
#define k8exp(x) xexp(x)
#define k8fmod(x, y) K8REPL2(fmod, x, y)
#define k8log(x) xlog(x)
#define k8pow(x, a) xpow(x, vec_set1(a))
#define k8sin(x) xsin(x)
#define k8sinh(x) xsinh(x)
#define k8tan(x) xtan(x)
#define k8tanh(x) xtanh(x)

#define k8lnot(x) (vec_not(x))
#define k8land(x, y) (vec_and(x, y))
#define k8lor(x, y) (vec_or(x, y))
#define k8lxor(x, y) (vec_xor(x, y))
#define k8ifthen(x, y, z) (vec_sel(z, y, x))

#define k8cmpeq(x, y) (vec_cmpeq(x, y))
#define k8cmpne(x, y) (k8lnot(k8cmpeq(x, y)))
#define k8cmpgt(x, y) (vec_cmpgt(x, y))
#define k8cmpge(x, y) (k8lnot(k8cmplt(x, y)))
#define k8cmplt(x, y) (vec_cmplt(x, y))
#define k8cmple(x, y) (k8lnot(k8cmpgt(x, y)))

#endif
