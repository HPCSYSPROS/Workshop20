#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

#ifdef _OPENMP
#include <omp.h>
#else
static int omp_get_num_threads(void) { return 1; }
#endif

static double get_time(void) {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec / 1.0e+6;
}

/* Generate a matrix with random entries between 0.0 and 1.0 */
static void dmatgen(const ptrdiff_t m, const ptrdiff_t n,
                    CCTK_REAL *restrict const A, const ptrdiff_t ldA) {
  const CCTK_REAL x = rand() / (RAND_MAX + 1.0);

  for (ptrdiff_t j = 0; j < n; ++j) {
    for (ptrdiff_t i = 0; i < m; ++i) {
      A[i + ldA * j] = x;
    }
  }
}

/* DGEMM (for transposed A and un-transposed B) */
static void dgemm(const ptrdiff_t m, const ptrdiff_t n, const ptrdiff_t l,
                  const CCTK_REAL *restrict const A, const ptrdiff_t ldA,
                  const CCTK_REAL *restrict const B, const ptrdiff_t ldB,
                  CCTK_REAL *restrict const C, const ptrdiff_t ldC) {
  for (ptrdiff_t j = 0; j < n; ++j) {
    for (ptrdiff_t i = 0; i < m; ++i) {
      CCTK_REAL s = C[i + ldC * j];
      for (ptrdiff_t k = 0; k < l; ++k) {
        s += A[k + ldA * i] * B[k + ldB * j];
      }
      C[i + ldC * j] = s;
    }
  }
}

/* "Use" a matrix */
static CCTK_REAL dmatuse(const ptrdiff_t m, const ptrdiff_t n,
                         const CCTK_REAL *restrict const A,
                         const ptrdiff_t ldA) {
  CCTK_REAL s = 0.0;

  for (ptrdiff_t j = 0; j < n; ++j) {
    for (ptrdiff_t i = 0; i < m; ++i) {
      s += A[i + ldA * j];
    }
  }

  return s;
}

static void run_dgemm(const ptrdiff_t N) {
  DECLARE_CCTK_PARAMETERS;

#pragma omp barrier
  const double t0 = get_time();

  CCTK_REAL *restrict const A = malloc(N * N * sizeof *A);
  CCTK_REAL *restrict const B = malloc(N * N * sizeof *B);
  CCTK_REAL *restrict const C = malloc(N * N * sizeof *C);

  dmatgen(N, N, A, N);
  dmatgen(N, N, B, N);
  dmatgen(N, N, C, N);

  const double t1 = get_time();
#pragma omp barrier
  dgemm(N, N, N, A, N, B, N, C, N);
#pragma omp barrier
  const double t2 = get_time();
  const int num_threads = omp_get_num_threads();
  const double Gflop = num_threads * 2.0 * N * N * N / 1.0e+9;

  const volatile CCTK_REAL s CCTK_ATTRIBUTE_UNUSED = dmatuse(N, N, C, N);

  free(A);
  free(B);
  free(C);

#pragma omp barrier
  const double t3 = get_time();

  const double Gflop_sec_pure = Gflop / (t2 - t1);
  const double Gflop_sec_overall = Gflop / (t3 - t0);

#pragma omp master
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Floating point operations: %g Gflop over %d threads", Gflop,
               num_threads);
    CCTK_VInfo(CCTK_THORNSTRING, "Pure performance:          %g Gflop/sec",
               Gflop_sec_pure);
    CCTK_VInfo(CCTK_THORNSTRING, "Overall performance:       %g Gflop/sec",
               Gflop_sec_overall);
  }
}

void PAPI_dgemm(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const ptrdiff_t N = dgemm_N;

  CCTK_VInfo(CCTK_THORNSTRING,
             "Running local single-threaded DGEMM benchmark with N=%td...", N);
  run_dgemm(N);

  CCTK_VInfo(CCTK_THORNSTRING,
             "Running local multi-threaded DGEMM benchmark with N=%td...", N);
#pragma omp parallel
  { run_dgemm(N); }

  CCTK_INFO("NOTE: These performance numbers do not have to be large, but they "
            "need to agree with what PAPI reports");
}
