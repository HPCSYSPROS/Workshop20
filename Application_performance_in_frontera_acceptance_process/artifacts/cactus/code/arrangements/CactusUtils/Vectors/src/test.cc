#include "vectors.h"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>

using namespace std;

inline CCTK_REAL my_sgn(CCTK_REAL const x) {
  return x == (CCTK_REAL)0.0 ? (CCTK_REAL)0.0 : copysign((CCTK_REAL)1.0, x);
}

#define SCALARTEST(testname, vecexpr, scalarexpr)                              \
  do {                                                                         \
    if (verbose) {                                                             \
      CCTK_VInfo(CCTK_THORNSTRING, "Test %s...", testname);                    \
      fflush(stdout);                                                          \
    }                                                                          \
    CCTK_REAL const res = (scalarexpr);                                        \
    CCTK_REAL const vecres = (vecexpr);                                        \
    CCTK_REAL const eps = numeric_limits<CCTK_REAL>::epsilon();                \
    if ((fabs(vecres - res) <= 10 * eps) or                                    \
        (CCTK_isnan(vecres) and CCTK_isnan(res))) {                            \
      passed++;                                                                \
    } else {                                                                   \
      CCTK_VParamWarn(CCTK_THORNSTRING,                                        \
                      "Failed test %s: expected %.17g, received %.17g",        \
                      testname, (double)res, (double)vecres);                  \
    }                                                                          \
    numtests++;                                                                \
  } while (0)

#define VECTEST(testname, vecexpr, scalarexpr)                                 \
  do {                                                                         \
    if (verbose) {                                                             \
      CCTK_VInfo(CCTK_THORNSTRING, "Test %s...", testname);                    \
      fflush(stdout);                                                          \
    }                                                                          \
    CCTK_REAL_VEC rv = (vecexpr);                                              \
    for (int i = 0; i < CCTK_REAL_VEC_SIZE; i++) {                             \
      CCTK_REAL res = (scalarexpr);                                            \
      CCTK_REAL vecres = vec_elt(rv, i);                                       \
      CCTK_REAL eps = numeric_limits<CCTK_REAL>::epsilon();                    \
      if ((fabs(vecres - res) <= 10 * eps) or                                  \
          (CCTK_isnan(vecres) and CCTK_isnan(res))) {                          \
        passed++;                                                              \
      } else {                                                                 \
        CCTK_VParamWarn(CCTK_THORNSTRING,                                      \
                        "Failed test %s: "                                     \
                        "for element %d, expected %.17g, received %.17g",      \
                        testname, i, (double)res, (double)vecres);             \
      }                                                                        \
      numtests++;                                                              \
    }                                                                          \
  } while (0)

#define VECBITTEST(testname, vecexpr, scalarexpr)                              \
  do {                                                                         \
    if (verbose) {                                                             \
      CCTK_VInfo(CCTK_THORNSTRING, "Test %s...", testname);                    \
      fflush(stdout);                                                          \
    }                                                                          \
    CCTK_BOOLEAN_VEC rv = (vecexpr);                                           \
    for (int i = 0; i < CCTK_REAL_VEC_SIZE; i++) {                             \
      CCTK_BOOLEAN res = (scalarexpr);                                         \
      CCTK_BOOLEAN vecres = vec_eltb(rv, i);                                   \
      if (memcmp(&vecres, &res, sizeof vecres) == 0) {                         \
        passed++;                                                              \
      } else {                                                                 \
        CCTK_INTEGER ires, ivecres;                                            \
        memcpy(&ires, &res, sizeof ires);                                      \
        memcpy(&ivecres, &vecres, sizeof ivecres);                             \
        CCTK_VParamWarn(CCTK_THORNSTRING,                                      \
                        "Failed test %s: "                                     \
                        "for element %d, expected %lld, received %lld",        \
                        testname, i, (long long)ires, (long long)ivecres);     \
      }                                                                        \
      numtests++;                                                              \
    }                                                                          \
  } while (0)

extern "C" void Vectors_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INFO("Testing vectorisation... [errors may result in segfaults]");
  fflush(stdout);

  char testname[100];

  int passed = 0, numtests = 0;

  CCTK_REAL a[CCTK_REAL_VEC_SIZE];
  CCTK_REAL b[CCTK_REAL_VEC_SIZE];
  CCTK_REAL c[CCTK_REAL_VEC_SIZE];
  CCTK_REAL x[CCTK_REAL_VEC_SIZE]; // small
  CCTK_REAL y[CCTK_REAL_VEC_SIZE]; // small and positive
  CCTK_REAL z[CCTK_REAL_VEC_SIZE]; // >=1

  for (int i = 0; i < CCTK_REAL_VEC_SIZE; i++) {
    a[i] = (i + 1) * 1.23456789;
    b[i] = -(i + 1) * 9.87654321;
    c[i] = (i + 1) * 1.01010101;
    x[i] = (i + 1) * 0.00123 * (i % 2 ? +1 : -1);
    y[i] = (i + 1) * 0.00234;
    z[i] = (i)*0.00234 + 1.0;
  }

  CCTK_REAL_VEC av = vec_loadu(a[0]);
  CCTK_REAL_VEC bv = vec_loadu(b[0]);
  CCTK_REAL_VEC cv = vec_loadu(c[0]);
  CCTK_REAL_VEC xv = vec_loadu(x[0]);
  CCTK_REAL_VEC yv = vec_loadu(y[0]);
  CCTK_REAL_VEC zv = vec_loadu(z[0]);

  /* l and lv are similar to a and av, except that it is larger, and
     guaranteed to be aligned */
  CCTK_REAL_VEC lv[4];
  lv[0] = vec_loadu(a[0]);
  lv[1] = vec_loadu(b[0]);
  lv[2] = vec_loadu(c[0]);
  lv[3] = vec_loadu(a[0]);
  CCTK_REAL *const l = (CCTK_REAL *)&lv[0];

  /* s and sv are similar to a and av, but are aligned and not
     initialised */
  CCTK_REAL_VEC sv;
  CCTK_REAL *const s = (CCTK_REAL *)&sv;

  VECTEST("vec_set1", vec_set1(a[0]), a[0]);
#if CCTK_REAL_VEC_SIZE == 1
  VECTEST("vec_set", vec_set(a[0]), a[i]);
#elif CCTK_REAL_VEC_SIZE == 2
  VECTEST("vec_set", vec_set(a[0], a[1]), a[i]);
#elif CCTK_REAL_VEC_SIZE == 4
  VECTEST("vec_set", vec_set(a[0], a[1], a[2], a[3]), a[i]);
#elif CCTK_REAL_VEC_SIZE == 8
  VECTEST("vec_set", vec_set(a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]),
          a[i]);
#else
#error "Unsupported vector size"
#endif
  for (int d = 0; d < CCTK_REAL_VEC_SIZE; ++d) {
    snprintf(testname, sizeof testname, "vec_elt[%d]", d);
    SCALARTEST(testname, vec_elt(av, d), a[d]);
  }

  /* These tests will probably fail with a segfault, if they fail */
  VECTEST("vec_load", vec_load(*l), l[i]);
  for (int d = 0; d < CCTK_REAL_VEC_SIZE; ++d) {
    snprintf(testname, sizeof testname, "vec_loadu[%d]", d);
    VECTEST(testname, vec_loadu(l[d]), l[i + d]);
  }
  for (int d = 0; d < CCTK_REAL_VEC_SIZE; ++d) {
    snprintf(testname, sizeof testname, "vec_loadu_maybe[%d]", d);
    VECTEST(testname, vec_loadu_maybe(d, l[d]), l[i + d]);
  }
  for (int d1 = 0; d1 < CCTK_REAL_VEC_SIZE; ++d1) {
    for (int d2 = 0; d2 < CCTK_REAL_VEC_SIZE; ++d2) {
      for (int d3 = 0; d3 < CCTK_REAL_VEC_SIZE; ++d3) {
        if (not(VECTORISE and VECTORISE_ALIGNED_ARRAYS) or
            (d2 == 0 and d3 == 0)) {
          snprintf(testname, sizeof testname, "vec_loadu_maybe3[%d,%d,%d]", d1,
                   d2, d3);
          VECTEST(testname, vec_loadu_maybe3(d1, d2, d3, l[d1 + d2 + d3]),
                  l[i + d1 + d2 + d3]);
        }
      }
    }
  }

  /* These tests may fail with a segfault, if they fail */
  sv = av;
  vec_store(*s, bv);
  VECTEST("vec_store", sv, b[i]);
  sv = av;
  vec_store_nta(*s, bv);
  VECTEST("vec_store_nta", sv, b[i]);
  int const arenaitems = 3;
  int const arenasize = arenaitems * CCTK_REAL_VEC_SIZE;
  for (int ilo = 0; ilo < arenasize; ++ilo) {
    for (int ihi = 0; ihi < arenasize; ++ihi) {
      // Initialise arena to a everywhere
      for (int n = 0; n < arenaitems; ++n) {
        lv[n] = av;
      }
      // Set arena to b in the "interior"
      for (int i = (ilo & -CCTK_REAL_VEC_SIZE); i < ihi;
           i += CCTK_REAL_VEC_SIZE) {
        vec_store_partial_prepare(i, ilo, ihi);
        vec_store_nta_partial(l[i], bv);
      }
      for (int i = 0; i < arenasize; ++i) {
        snprintf(testname, sizeof testname, "vec_store_nta_partial[%d..%d][%d]",
                 ilo, ihi, i);
        SCALARTEST(
            testname, l[i],
            (i >= ilo and i < ihi ? b : a)[i & (CCTK_REAL_VEC_SIZE - 1)]);
      }
    }
  }
  /* The partial stores are not implemented for d==0 and
     d==CCTK_REAL_VEC_SIZE (because these are trivial) */
  for (int d = 1; d < CCTK_REAL_VEC_SIZE - 1; ++d) {
    sv = av;
    vec_store_nta_partial_lo(*s, bv, d);
    snprintf(testname, sizeof testname, "vec_store_nta_partial_lo[%d]", d);
    VECTEST(testname, sv, i < d ? b[i] : a[i]);
  }
  for (int d = 1; d < CCTK_REAL_VEC_SIZE - 1; ++d) {
    sv = av;
    vec_store_nta_partial_hi(*s, bv, d);
    snprintf(testname, sizeof testname, "vec_store_nta_partial_hi[%d]", d);
    VECTEST(testname, sv, i >= CCTK_REAL_VEC_SIZE - d ? b[i] : a[i]);
  }
  for (int dlo = 1; dlo < CCTK_REAL_VEC_SIZE - 1; ++dlo) {
    for (int dhi = 1; dhi < CCTK_REAL_VEC_SIZE - 1; ++dhi) {
      sv = av;
      vec_store_nta_partial_mid(*s, bv, dlo, dhi);
      snprintf(testname, sizeof testname, "vec_store_nta_partial_mid[%d,%d]",
               dlo, dhi);
      VECTEST(testname, sv,
              i < dlo and i >= CCTK_REAL_VEC_SIZE - dhi ? b[i] : a[i]);
    }
  }

  VECTEST("kneg", kneg(av), -a[i]);

  VECTEST("kadd", kadd(av, bv), a[i] + b[i]);
  VECTEST("ksub", ksub(av, bv), a[i] - b[i]);
  VECTEST("kmul", kmul(av, bv), a[i] * b[i]);
  VECTEST("kdiv", kdiv(av, bv), a[i] / b[i]);

  VECTEST("kmadd", kmadd(av, bv, cv), +a[i] * b[i] + c[i]);
  VECTEST("kmsub", kmsub(av, bv, cv), +a[i] * b[i] - c[i]);
  VECTEST("knmadd", knmadd(av, bv, cv), -a[i] * b[i] - c[i]);
  VECTEST("knmsub", knmsub(av, bv, cv), -a[i] * b[i] + c[i]);

  VECTEST("kacos", kacos(xv), acos(x[i]));
  VECTEST("kacosh", kacosh(zv), acosh(z[i]));
  VECTEST("kasin", kasin(xv), asin(x[i]));
  VECTEST("kasinh", kasinh(xv), asinh(x[i]));
  VECTEST("katan", katan(xv), atan(x[i]));
  VECTEST("katan2", katan2(xv, yv), atan2(x[i], y[i]));
  VECTEST("katanh", katanh(xv), atanh(x[i]));
  VECTEST("kcopysign", kcopysign(xv, yv), copysign(x[i], y[i]));
  VECTEST("kcos", kcos(xv), cos(x[i]));
  VECTEST("kcosh", kcosh(xv), cosh(x[i]));
  VECTEST("kexp", kexp(xv), exp(x[i]));
  VECTEST("kfabs", kfabs(xv), fabs(x[i]));
  VECTEST("kfmax", kfmax(xv, yv), fmax(x[i], y[i]));
  VECTEST("kfmin", kfmin(xv, yv), fmin(x[i], y[i]));
  VECTEST("kfmod", kfmod(xv, yv), fmod(x[i], y[i]));
  VECTEST("kfnabs", kfnabs(xv), -fabs(x[i]));
  VECTEST("klog", klog(yv), log(y[i]));
  VECTEST("kpow", kpow(yv, x[0]), pow(y[i], x[0]));
  VECTEST("ksin", ksin(xv), sin(x[i]));
  VECTEST("ksinh", ksinh(xv), sinh(x[i]));
  VECTEST("ksgn", ksgn(xv), my_sgn(x[i]));
  VECTEST("ksqrt", ksqrt(yv), sqrt(y[i]));
  VECTEST("ktan", ktan(xv), tan(x[i]));
  VECTEST("ktanh", ktanh(xv), tanh(x[i]));

#if 0
  VECTEST("kifpos positive",
          kifpos(av, bv, cv), CCTK_SIGNBIT(a[i]) ? c[i] : b[i]);
  VECTEST("kifpos negative",
          kifpos(bv, bv, cv), CCTK_SIGNBIT(b[i]) ? c[i] : b[i]);
  VECTEST("kifpos 0",  kifpos(vec_set1(0.),bv,cv),  b[i]);
  VECTEST("kifpos -0", kifpos(vec_set1(-0.),bv,cv), c[i]);

  VECTEST("kifneg positive",
          kifneg(av, bv, cv), CCTK_SIGNBIT(a[i]) ? b[i] : c[i]);
  VECTEST("kifneg negative",
          kifneg(bv, bv, cv), CCTK_SIGNBIT(b[i]) ? b[i] : c[i]);
  VECTEST("kifneg 0",  kifneg(vec_set1(0.),bv,cv),  c[i]);
  VECTEST("kifneg -0", kifneg(vec_set1(-0.),bv,cv), b[i]);
#endif

  CCTK_BOOLEAN klfalse1 = vec_eltb(klfalse, 0);
  CCTK_BOOLEAN kltrue1 = vec_eltb(kltrue, 0);
  VECBITTEST("constant F", klfalse, klfalse1);
  VECBITTEST("constant T", kltrue, kltrue1);
  VECBITTEST("klnot F", klnot(klfalse), kltrue1);
  VECBITTEST("klnot T", klnot(kltrue), klfalse1);
  VECBITTEST("kland FF", kland(klfalse, klfalse), klfalse1);
  VECBITTEST("kland FT", kland(klfalse, kltrue), klfalse1);
  VECBITTEST("kland TF", kland(kltrue, klfalse), klfalse1);
  VECBITTEST("kland TT", kland(kltrue, kltrue), kltrue1);
  VECBITTEST("klor FF", klor(klfalse, klfalse), klfalse1);
  VECBITTEST("klor FT", klor(klfalse, kltrue), kltrue1);
  VECBITTEST("klor TF", klor(kltrue, klfalse), kltrue1);
  VECBITTEST("klor TT", klor(kltrue, kltrue), kltrue1);
  VECBITTEST("klxor FF", klxor(klfalse, klfalse), klfalse1);
  VECBITTEST("klxor FT", klxor(klfalse, kltrue), kltrue1);
  VECBITTEST("klxor TF", klxor(kltrue, klfalse), kltrue1);
  VECBITTEST("klxor TT", klxor(kltrue, kltrue), klfalse1);
  VECTEST("kifthen F", kifthen(klfalse, av, bv), b[i]);
  VECTEST("kifthen T", kifthen(kltrue, av, bv), a[i]);

  VECBITTEST("kcmpeq", kcmpeq(av, bv), a[i] == b[i] ? kltrue1 : klfalse1);
  VECBITTEST("kcmpne", kcmpne(av, bv), a[i] != b[i] ? kltrue1 : klfalse1);
  VECBITTEST("kcmpgt", kcmpgt(av, bv), a[i] > b[i] ? kltrue1 : klfalse1);
  VECBITTEST("kcmpge", kcmpge(av, bv), a[i] >= b[i] ? kltrue1 : klfalse1);
  VECBITTEST("kcmplt", kcmplt(av, bv), a[i] < b[i] ? kltrue1 : klfalse1);
  VECBITTEST("kcmple", kcmple(av, bv), a[i] <= b[i] ? kltrue1 : klfalse1);

  CCTK_VInfo(CCTK_THORNSTRING, "%d/%d tests passed ", passed, numtests);
  fflush(stdout);
  if (passed != numtests) {
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Failed %d correctness tests", numtests - passed);
  }
}
