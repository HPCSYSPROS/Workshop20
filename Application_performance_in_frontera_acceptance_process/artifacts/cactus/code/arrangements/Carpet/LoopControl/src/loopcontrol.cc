#include "loopcontrol.h"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#else

// Simple non-OpenMP implementations
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_num_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }

#endif

#if defined HAVE_CAPABILITY_CYCLECLOCK
// We have a fast, accurate clock

#include <cycleclock.h>

#elif defined _OPENMP
// We use the OpenMP clock

typedef double ticks;
static inline ticks getticks() { return omp_get_wtime(); }
static inline double elapsed(ticks t1, ticks t0) { return t1 - t0; }
static inline double seconds_per_tick() { return 1.0; }

#else
// We use gettimeofday as fallback

#include <sys/time.h>
typedef timeval ticks;
static inline ticks getticks() {
  timeval tp;
  gettimeofday(&tp, NULL);
  return tp;
}
static inline double elapsed(ticks t1, ticks t0) {
  return 1.0e+6 * (t1.tv_sec - t0.tv_sec) + (t1.tv_usec - t0.tv_usec);
}
static inline double seconds_per_tick() { return 1.0e-6; }

#endif

using namespace std;

static bool lc_do_explore_eagerly = false;
static bool lc_do_settle = false;

struct lc_thread_info_t {
  volatile int idx;            // linear index of next coarse thread block
} CCTK_ATTRIBUTE_ALIGNED(128); // align to prevent sharing cache lines

struct lc_fine_thread_comm_t {
  volatile int state;          // waiting threads
  volatile int value;          // broadcast value
} CCTK_ATTRIBUTE_ALIGNED(128); // align to prevent sharing cache lines

// One object per coarse thread: shared between fine threads
// Note: Since we use a vector, the individual elements may not
// actually be aligned, but they will still be spaced apart and thus
// be placed into different cache lines.
static vector<lc_fine_thread_comm_t> lc_fine_thread_comm;

// Statistics
struct lc_stats_t {
  double points, threads;
  double count, sum, sum2, min, max;
  int init_count;
  lc_stats_t()
      : points(0.0), threads(0.0), count(0.0), sum(0.0), sum2(0.0),
        min(numeric_limits<double>::max()), max(0.0), init_count(2) {}
  void add(const int npoints, const int nthreads, const double elapsed_time) {
    if (init_count > 0) {
      --init_count;
      // Reset statistics after the first iteration
      if (init_count == 0) {
        points = 0.0;
        threads = 0.0;
        count = 0.0;
        sum = 0.0;
        sum2 = 0.0;
        min = numeric_limits<double>::max();
        max = 0.0;
      }
    }
    points += double(npoints);
    threads += double(nthreads);
    count += 1.0;
    sum += elapsed_time;
    sum2 += pow(elapsed_time, 2.0);
    min = fmin(min, elapsed_time);
    max = fmax(max, elapsed_time);
  }
  double avg_thread() const { return sum / count; }
  double avg_point() const { return sum * threads / (count * points); }
};

// Parameters that determine how a loop is traversed. This corresponds
// to the choices one can make to optimize. This loosely corresponds
// to parameters of thorn LoopControl.
struct lc_params_key_t {
  lc_vec_t tilesize;
  lc_vec_t loopsize;

  bool operator==(const lc_params_key_t &x) const {
    return memcmp(this, &x, sizeof *this) == 0;
  }
  bool operator<(const lc_params_key_t &x) const {
    return memcmp(this, &x, sizeof *this) < 0;
  }
};

// Unique identifier for an iteration setup, where differing setups
// need to be optimized separately. This corresponds to the
// information passed to lc_control_init.
struct lc_setup_key_t {
  lc_vec_t min, max, ash;
  int num_coarse_threads, num_fine_threads;

  bool operator==(const lc_setup_key_t &x) const {
    return memcmp(this, &x, sizeof *this) == 0;
  }
  bool operator<(const lc_setup_key_t &x) const {
    return memcmp(this, &x, sizeof *this) < 0;
  }
};

struct lc_setup_t;
struct lc_params_t {
  lc_setup_t &setup;   // setup
  lc_params_key_t key; // copy of params

  lc_params_t(lc_setup_t &setup_, lc_params_key_t &key_)
      : setup(setup_), key(key_) {}

  lc_stats_t stats; // statistics for these params
};

struct lc_descr_t;
struct lc_setup_t {
  lc_descr_t &descr;  // descriptor
  lc_setup_key_t key; // copy of setup

  lc_setup_t(lc_descr_t &descr_, lc_setup_key_t &key_)
      : descr(descr_), key(key_), default_params(0), best_params(0),
        current_params(0) {}

  typedef map<lc_params_key_t, lc_params_t *> params_map_t;
  params_map_t params;
  lc_params_t *default_params;
  lc_params_t *best_params;

  lc_params_t *current_params;

  lc_stats_t stats; // statistics for all params for this setup
};

struct lc_descr_t {
  string name;
  string file;
  int line;

  lc_descr_t(const char *name_, const char *file_, int line_)
      : name(name_), file(file_), line(line_), current_setup(0),
        current_params(0) {}

  typedef map<lc_setup_key_t, lc_setup_t *> setup_map_t;
  setup_map_t setups;

  lc_setup_t *current_setup;   // current setup
  lc_params_t *current_params; // current params

  lc_stats_t stats; // global statistics for all setups
  ticks start_time; // current start time
};

// The Intel compiler keeps increasing the amount of "free" memory
// allocated by libc. Work around this by allocating memory in large
// batches.
namespace {

template <typename T> class mempool {
  T *next;
  int nleft;

public:
  mempool() : nleft(0) {}
  void *allocate() {
    if (nleft < 1) {
      nleft = 1000000 / sizeof(T);
      next = (T *)new char[nleft * sizeof(T)];
    }
    assert(nleft >= 1);
    return nleft--, next++;
  }
};

mempool<lc_params_t> params_mempool;
mempool<lc_setup_t> setup_mempool;
}

extern "C" CCTK_FCALL void
    CCTK_FNAME(lc_get_fortran_type_sizes)(ptrdiff_t *type_sizes);

namespace {

typedef vector<lc_descr_t *> all_descrs_t;
all_descrs_t all_descrs;

struct descr_comp_name {
  bool operator()(const lc_descr_t *a, const lc_descr_t *b) const {
    return a->name < b->name;
  }
};

struct params_comp_time {
  bool operator()(const lc_params_t *a, const lc_params_t *b) const {
    return a->stats.avg_point() < b->stats.avg_point();
  }
};

void check_fortran_type_sizes() {
  ptrdiff_t type_sizes[3];
  CCTK_FNAME(lc_get_fortran_type_sizes)(type_sizes);
  assert(type_sizes[0] == sizeof(lc_vec_t));
  assert(type_sizes[1] == sizeof(lc_space_t));
  assert(type_sizes[2] == sizeof(lc_control_t));
}

template <typename T> T divup(const T i, const T j) {
  assert(i >= 0 and j > 0);
  return (i + j - 1) / j;
}

template <typename T> T divdown(const T i, const T j) {
  assert(i >= 0 and j > 0);
  return i / j;
}

template <typename T> T divexact(const T i, const T j) {
  assert(i % j == 0);
  return i / j;
}

template <typename T> T moddown(const T i, const T j) {
  assert(i >= 0 and j > 0);
  return i % j;
}

template <typename T> T alignup(const T i, const T j) {
  return divup(i, j) * j;
}

template <typename T> T aligndown(const T i, const T j) {
  return divdown(i, j) * j;
}

// random uniform integer
template <typename T> T randomui(const T imin, const T imax, const T istr = 1) {
  assert(imin < imax);
  // const T res =
  //   imin + istr * floor(rand() / (RAND_MAX + 1.0) * (imax - imin) / istr);
  const T res =
      imin +
      istr * llrint(floor(random() / (RAND_MAX + 1.0) * (imax - imin) / istr));
  assert(res >= imin and res < imax and (res - imin) % istr == 0);
  return res;
}

ostream &operator<<(ostream &os, const lc_vec_t &x) {
  os << "[";
  for (int d = 0; d < LC_DIM; ++d) {
    if (d > 0)
      os << ",";
    os << x.v[d];
  }
  os << "]";
  return os;
}

ostream &operator<<(ostream &os, const lc_space_t &s) {
  os << "{"
     << "min:" << s.min << ","
     << "max:" << s.max << ","
     << "step:" << s.step << ","
     << "pos:" << s.pos << ","
     << "count:" << s.count << ","
     << "idx:" << s.idx << "}";
  return os;
}

ostream &operator<<(ostream &os, const lc_control_t &c) {
  os << "lc_control{\n"
     << "   ash:" << c.ash << ",\n"
     << "   overall:" << c.overall << ",\n"
     << "   coarse_thread:" << c.coarse_thread << ",\n"
     << "   coarse_loop:" << c.coarse_loop << ",\n"
     << "   fine_loop:" << c.fine_loop << "\n"
     << "   fine_thread:" << c.fine_thread << "\n"
     << "}\n";
  return os;
}

ptrdiff_t prod(const lc_vec_t &x) {
  ptrdiff_t r = 1;
  for (int d = 0; d < LC_DIM; ++d) {
    assert(x.v[d] >= 0);
    r *= x.v[d];
  }
  return r;
}

ptrdiff_t ind(const lc_vec_t &shape, const lc_vec_t &pos) {
  ptrdiff_t r = 0;
  ptrdiff_t f = 1;
  for (int d = 0; d < LC_DIM; ++d) {
    assert(pos.v[d] >= 0 and pos.v[d] < shape.v[d]);
    r += f * pos.v[d];
    assert(shape.v[d] >= 0);
    f *= shape.v[d];
  }
  return r;
}

ptrdiff_t ind(const lc_vec_t &shape, const ptrdiff_t i, const ptrdiff_t j,
              const ptrdiff_t k) {
  const lc_vec_t pos = {{i, j, k}};
  return ind(shape, pos);
}

void space_set_count(lc_space_t &space) {
  for (int d = 0; d < LC_DIM; ++d) {
    space.count.v[d] = divup(space.max.v[d] - space.min.v[d], space.step.v[d]);
  }
}

void space_idx2pos(lc_space_t &space) {
  for (int d = 0; d < LC_DIM; ++d) {
    space.pos.v[d] = space.min.v[d] + space.idx.v[d] * space.step.v[d];
  }
}

bool space_global2local(lc_space_t &space, ptrdiff_t gidx) {
  assert(gidx >= 0);
  for (int d = 0; d < LC_DIM; ++d) {
    if (space.count.v[d] > 0) {
      space.idx.v[d] = moddown(gidx, space.count.v[d]);
      gidx = divdown(gidx, space.count.v[d]);
    } else {
      space.idx.v[d] = 0;
    }
  }
  return gidx != 0;
}

int space_local2global(const lc_space_t &space) {
  int gidx = 0;
  int fact = 1;
  for (int d = 0; d < LC_DIM; ++d) {
    assert(space.idx.v[d] >= 0 and space.idx.v[d] < space.count.v[d]);
    gidx += fact * space.idx.v[d];
    fact *= space.count.v[d];
  }
  return gidx;
}

static int num_smt_threads = 0;

int get_num_fine_threads() {
  DECLARE_CCTK_PARAMETERS;
  if (not use_smt_threads)
    return 1;
  if (omp_get_num_threads() == 1)
    return 1;
  return num_smt_threads;
}

int get_fine_thread_num() {
  DECLARE_CCTK_PARAMETERS;
  if (not use_smt_threads)
    return 0;
  if (omp_get_num_threads() == 1)
    return 0;
  const int thread_num = omp_get_thread_num();
  const int num_fine_threads = get_num_fine_threads();
  return moddown(thread_num, num_fine_threads);
}

int get_num_coarse_threads() {
  const int num_threads = omp_get_num_threads();
  const int num_fine_threads = get_num_fine_threads();
  return divexact(num_threads, num_fine_threads);
}

int get_coarse_thread_num() {
  const int thread_num = omp_get_thread_num();
  const int num_fine_threads = get_num_fine_threads();
  return divdown(thread_num, num_fine_threads);
}

// Wait until *ptr is different from old_value
void thread_wait(volatile int *const ptr, const int old_value) {
  while (*ptr == old_value) {
#pragma omp flush
    (void)0; // PGI compiler needs this
  }
}

int fine_thread_broadcast(lc_fine_thread_comm_t *const comm, int value) {
  const int num_fine_threads = get_num_fine_threads();
  if (num_fine_threads == 1)
    return value;
  assert(num_fine_threads < 8 * int(sizeof comm->state));
  const int fine_thread_num = get_fine_thread_num();
  const int master_mask = 1;

  // Assume comm->count == 0 initially
  if (fine_thread_num == 0) { // if on master

    const int all_threads_mask = (1 << num_fine_threads) - 1;
    if (comm->state != 0) {
// wait until everybody has acknowledged the previous value
#pragma omp flush
      for (;;) {
        const int state = comm->state;
        if (state == all_threads_mask)
          break;
        thread_wait(&comm->state, state);
      }
      // mark the value as invalid
      comm->state = 0;
#pragma omp flush
    }
    // publish value
    comm->value = value;
#pragma omp flush
    // mark the value as valid
    comm->state = master_mask;
#pragma omp flush

  } else { // if not on master

    // wait until the value is valid, and it is a new value
    const int thread_mask = 1 << fine_thread_num;
#pragma omp flush
    for (;;) {
      const int state = comm->state;
      if ((state & (master_mask | thread_mask)) == master_mask)
        break;
      thread_wait(&comm->state, state);
    }
    // read value
    value = comm->value;
#pragma omp flush
    (void)0; // PGI compiler needs this
             // acknowledge the value
#pragma omp atomic
    comm->state |= thread_mask;
#pragma omp flush

  } // if not on master

  return value;
}

} // namespace

void lc_descr_init(lc_descr_t **const descr_ptr, const char *const name,
                   const char *const file, const int line) {
  if (CCTK_BUILTIN_EXPECT(*descr_ptr != 0, true))
    return;

#pragma omp barrier
#pragma omp master
  {
    lc_descr_t *const descr = new lc_descr_t(name, file, line);
    all_descrs.push_back(descr);
    *descr_ptr = descr;

    // Determine number of SMT threads
    if (CCTK_BUILTIN_EXPECT(num_smt_threads == 0, false)) {
      if (CCTK_IsFunctionAliased("GetNumSMTThreads")) {
        num_smt_threads = GetNumSMTThreads();
      } else {
        num_smt_threads = 1;
      }
    }

    // Allocate fine thread communicators
    if (CCTK_BUILTIN_EXPECT(lc_fine_thread_comm.empty(), false)) {
      lc_fine_thread_comm.resize(omp_get_max_threads());
    }
  }
#pragma omp barrier
}

void lc_control_init(lc_control_t *restrict const control,
                     lc_descr_t *const descr, ptrdiff_t imin, ptrdiff_t jmin,
                     ptrdiff_t kmin, ptrdiff_t imax, ptrdiff_t jmax,
                     ptrdiff_t kmax, ptrdiff_t iash, ptrdiff_t jash,
                     ptrdiff_t kash, ptrdiff_t istr) {
  DECLARE_CCTK_PARAMETERS;

  // Get cache line size
  static ptrdiff_t max_cache_linesize = -1;
  if (CCTK_BUILTIN_EXPECT(max_cache_linesize < 0, false)) {
#pragma omp barrier
#pragma omp master
    {
      max_cache_linesize = 1;
      if (CCTK_IsFunctionAliased("GetCacheInfo1")) {
        const int num_levels =
            GetCacheInfo1(NULL, NULL, NULL, NULL, NULL, NULL, 0);
        vector<CCTK_INT> types(num_levels);
        vector<CCTK_INT> linesizes(num_levels);
        vector<CCTK_INT> strides(num_levels);
        GetCacheInfo1(NULL, &types[0], NULL, &linesizes[0], &strides[0], NULL,
                      num_levels);
        for (int level = 0; level < num_levels; ++level) {
          if (types[level] == 0) { // if this is a cache
            max_cache_linesize =
                max(max_cache_linesize, ptrdiff_t(linesizes[level]));
          }
        }
      }
    }
#pragma omp barrier
  }

  ptrdiff_t tilesize_alignment = 1;
  if (align_with_cachelines) {
    tilesize_alignment =
        divup(max_cache_linesize, ptrdiff_t(sizeof(CCTK_REAL)));
    tilesize_alignment = alignup(tilesize_alignment, istr);
  }

#pragma omp barrier
#pragma omp master
  {
    // Start timing
    descr->start_time = getticks();

    // Capture loop setup key
    lc_setup_key_t setup_key;
    setup_key.min.v[0] = imin;
    setup_key.min.v[1] = jmin;
    setup_key.min.v[2] = kmin;
    setup_key.max.v[0] = imax;
    setup_key.max.v[1] = jmax;
    setup_key.max.v[2] = kmax;
    setup_key.ash.v[0] = iash;
    setup_key.ash.v[1] = jash;
    setup_key.ash.v[2] = kash;
    setup_key.num_coarse_threads = get_num_coarse_threads();
    setup_key.num_fine_threads = get_num_fine_threads();

    // Determine loop setup
    {
      const pair<lc_descr_t::setup_map_t::iterator, bool> res =
          descr->setups.insert(
              make_pair(setup_key, static_cast<lc_setup_t *>(0)));
      const lc_descr_t::setup_map_t::iterator setup_i = res.first;
      lc_setup_t *&setup_p = setup_i->second;
      const bool isnew = res.second;
      assert(isnew == not setup_p);
      if (isnew) {
        void *ptr = setup_mempool.allocate();
        setup_p = new (ptr) lc_setup_t(*descr, setup_key);
      }
      assert(not descr->current_setup);
      descr->current_setup = setup_p;
    }

    // Choose loop params

    lc_setup_t &setup = *descr->current_setup;

    enum choices_t {
      choice_set_default,
      choice_keep_current,
      choice_use_best,
      choice_random_jump
    };
    choices_t choice = choice_set_default;

    if (setup.current_params) {
      choice = choice_keep_current;

      if (setup.current_params->stats.avg_point() >
          very_expensive_factor * setup.best_params->stats.avg_point()) {
        // Bail out if this params setting is too expensive
        choice = choice_use_best;
      }
      if (setup.current_params->stats.count >= double(tryout_iterations)) {
        // Switch if we tried this setting for some time
        choice = choice_use_best;
      }
    }
    if (choice == choice_use_best) {
      // Make a random jump every so often
      if (not lc_do_settle and
          (lc_do_explore_eagerly or
           random() / (RAND_MAX + 1.0) < random_jump_probability)) {
        choice = choice_random_jump;
      }
    }

    lc_params_key_t params_key;
    switch (choice) {
    case choice_set_default:
      // Set default
      params_key.tilesize.v[0] =
          alignup(ptrdiff_t(tilesize_i), tilesize_alignment);
      params_key.tilesize.v[1] = tilesize_j;
      params_key.tilesize.v[2] = tilesize_k;
      params_key.loopsize.v[0] =
          alignup(ptrdiff_t(loopsize_i), params_key.tilesize.v[0]);
      params_key.loopsize.v[1] =
          alignup(ptrdiff_t(loopsize_j), params_key.tilesize.v[1]);
      params_key.loopsize.v[2] =
          alignup(ptrdiff_t(loopsize_k), params_key.tilesize.v[2]);
      break;
    case choice_keep_current:
      params_key = setup.current_params->key;
      break;
    case choice_use_best:
      params_key = setup.best_params->key;
      break;
    case choice_random_jump: {
      const ptrdiff_t tilesizes[LC_DIM] = {tilesize_i, tilesize_j, tilesize_k};
      const ptrdiff_t loopsizes[LC_DIM] = {loopsize_i, loopsize_j, loopsize_k};
      for (int d = 0; d < LC_DIM; ++d) {
        const ptrdiff_t align = d == 0 ? int(tilesize_alignment) : 1;
        params_key.tilesize.v[d] = randomui(
            align, alignup(max_size_factor * tilesizes[d], align), align);
        const ptrdiff_t tilesize = params_key.tilesize.v[d];
        params_key.loopsize.v[d] = randomui(
            tilesize, alignup(max_size_factor * loopsizes[d], tilesize),
            tilesize);
        assert(moddown(params_key.tilesize.v[d], align) == 0);
        assert(moddown(params_key.loopsize.v[d], params_key.tilesize.v[d]) ==
               0);
      }
      break;
    }
    default:
      assert(0);
    }

    // Determine loop params
    {
      const pair<lc_setup_t::params_map_t::iterator, bool> res =
          setup.params.insert(
              make_pair(params_key, static_cast<lc_params_t *>(0)));
      const lc_setup_t::params_map_t::iterator params_i = res.first;
      lc_params_t *&params_p = params_i->second;
      const bool isnew = res.second;
      assert(isnew == not params_p);
      if (isnew) {
        void *ptr = params_mempool.allocate();
        params_p = new (ptr) lc_params_t(setup, params_key);
      }
      assert(not descr->current_params);
      descr->current_params = params_p;
      setup.current_params = descr->current_params;
      if (not setup.default_params) {
        setup.default_params = setup.current_params;
      }
      if (not setup.best_params) {
        setup.best_params = setup.current_params;
      }
    }
  }
#pragma omp barrier

  // Ensure thread counts are consistent
  assert(get_num_coarse_threads() * get_num_fine_threads() ==
         omp_get_num_threads());

#if 0
  {
    lc_descr_t* global_descr;
    int global_num_threads;
    int global_num_coarse_threads;
    int global_num_fine_threads;
    static int is_inconsistent;
#pragma omp single copyprivate(global_descr, global_num_threads,               \
                               global_num_coarse_threads,                      \
                               global_num_fine_threads)
    {
      global_descr = descr;
      global_num_threads = omp_get_num_threads();
      global_num_coarse_threads = get_num_coarse_threads();
      global_num_fine_threads = get_num_fine_threads();
      is_inconsistent = 0;
    }
#pragma omp atomic
    is_inconsistent |=
      global_descr != descr or
      global_num_threads != omp_get_num_threads() or
      global_num_coarse_threads != get_num_coarse_threads() or
      global_num_fine_threads != get_num_fine_threads();
#pragma omp barrier
    if (is_inconsistent) {
#pragma omp critical
      cout << "thread: " << omp_get_thread_num() << "\n"
           << "   loop name: " << descr->name << "\n"
           << "   file: " << descr->file << ":" << descr->line << "\n" << flush;
#pragma omp barrier
#pragma omp critical
      CCTK_ERROR("Thread inconsistency");
    }
  }
#endif

  // Initialize everything with a large, bogus value
  memset(control, 123, sizeof *control);

  // Parameters (all in units of grid points)
  const ptrdiff_t tilesize[LC_DIM] = {
      descr->current_params->key.tilesize.v[0],
      descr->current_params->key.tilesize.v[1],
      descr->current_params->key.tilesize.v[2],
  };
  const ptrdiff_t loopsize[LC_DIM] = {
      descr->current_params->key.loopsize.v[0],
      descr->current_params->key.loopsize.v[1],
      descr->current_params->key.loopsize.v[2],
  };
  ptrdiff_t smt_size[LC_DIM] = {1, 1, 1};
  {
    const int num_fine_threads = get_num_fine_threads();
    // If possible, stagger fine threads in the i direction, so that
    // they share cache lines
    if (istr * num_fine_threads <= loopsize[0]) {
      smt_size[0] = num_fine_threads;
    } else if (num_fine_threads <= loopsize[1]) {
      smt_size[1] = num_fine_threads;
    } else {
      smt_size[2] = num_fine_threads;
    }
  }

  // Arguments
  const ptrdiff_t loop_min[LC_DIM] = {imin, jmin, kmin};
  const ptrdiff_t loop_max[LC_DIM] = {imax, jmax, kmax};
  const ptrdiff_t ash[LC_DIM] = {iash, jash, kash};
  const ptrdiff_t vect_size[LC_DIM] = {istr, 1, 1};

  // Copy ash arguments
  for (int d = 0; d < LC_DIM; ++d) {
    control->ash.v[d] = ash[d];
  }

  // Set up multithreading state
  {
    lc_thread_info_t *thread_info_ptr;
#pragma omp single copyprivate(thread_info_ptr)
    { thread_info_ptr = new lc_thread_info_t; }
    control->coarse_thread_info_ptr = thread_info_ptr;
  }

  // Set loop sizes
  for (int d = 0; d < LC_DIM; ++d) {
    // Overall loop: as specified
    control->overall.min.v[d] = loop_min[d];
    control->overall.max.v[d] = loop_max[d];
// Thread loop
#if VECTORISE && VECTORISE_ALIGNED_ARRAYS
    // Move start to be aligned with vector size
    control->coarse_thread.min.v[d] =
        aligndown(control->overall.min.v[d], vect_size[d]);
#else
    control->coarse_thread.min.v[d] = control->overall.min.v[d];
#endif
    control->coarse_thread.max.v[d] = loop_max[d];
    // Fine threads
    control->fine_thread.count.v[d] = smt_size[d];
  }
  {
    const int fine_thread_num = get_fine_thread_num();
    const bool outside =
        space_global2local(control->fine_thread, fine_thread_num);
    assert(not outside);
  }

  // Set loop step sizes
  if (CCTK_EQUALS(initial_setup, "legacy")) {
    // Like a non-LoopControl loop: no loop tiling (i.e. do not use
    // coarse loops), parallelise only in k direction (i.e. assign
    // equal k ranges to threads)
    for (int d = 0; d < LC_DIM; ++d) {
      assert(smt_size[d] == 1); // TODO: implement this
      control->fine_thread.step.v[d] = vect_size[d];
      control->fine_loop.step.v[d] = vect_size[d];
      const ptrdiff_t npoints =
          control->overall.max.v[d] - control->overall.min.v[d];
      const ptrdiff_t nthreads = d != LC_DIM - 1 ? 1 : get_num_coarse_threads();
      control->coarse_loop.step.v[d] =
          alignup(divup(npoints, nthreads), control->fine_loop.step.v[d]);
      control->coarse_thread.step.v[d] =
          alignup(npoints, control->coarse_loop.step.v[d]);
    }
  } else if (CCTK_EQUALS(initial_setup, "tiled")) {
    // Basic LoopControl setup
    for (int d = 0; d < LC_DIM; ++d) {
      control->fine_thread.step.v[d] = vect_size[d];
      control->fine_loop.step.v[d] =
          smt_size[d] * control->fine_thread.step.v[d];
      control->coarse_loop.step.v[d] =
          alignup(tilesize[d], control->fine_loop.step.v[d]);
      control->coarse_thread.step.v[d] =
          alignup(loopsize[d], control->coarse_loop.step.v[d]);
    }
  } else {
    assert(0);
  }

  if (veryverbose) {
#pragma omp master
    CCTK_VInfo(
        CCTK_THORNSTRING,
        "Loop %s (%s:%d): imin=[%td,%td,%td] imax=[%td,%td,%td]\n"
        "   threads=%d coarse_threads=%d fine_threads=%d\n"
        "   fine_thread.step=[%td,%td,%td] fine_loop.step=[%td,%td,%td] "
        "coarse_loop.step=[%td,%td,%td] coarse_thread.step=[%td,%td,%td]",
        descr->name.c_str(), descr->file.c_str(), descr->line,
        control->overall.min.v[0], control->overall.min.v[1],
        control->overall.min.v[2], control->overall.max.v[0],
        control->overall.max.v[1], control->overall.max.v[2],
        omp_get_num_threads(), get_num_coarse_threads(), get_num_fine_threads(),
        control->fine_thread.step.v[0], control->fine_thread.step.v[1],
        control->fine_thread.step.v[2], control->fine_loop.step.v[0],
        control->fine_loop.step.v[1], control->fine_loop.step.v[2],
        control->coarse_loop.step.v[0], control->coarse_loop.step.v[1],
        control->coarse_loop.step.v[2], control->coarse_thread.step.v[0],
        control->coarse_thread.step.v[1], control->coarse_thread.step.v[2]);
  }

  // Initialise selftest
  if (selftest) {
    unsigned char *selftest_array;
#pragma omp single copyprivate(selftest_array)
    {
      const ptrdiff_t npoints = prod(control->ash);
      selftest_array = new unsigned char[npoints];
      memset(selftest_array, 0, npoints * sizeof *selftest_array);
    }
    control->selftest_array = selftest_array;
  } else {
    control->selftest_array = NULL;
  }
}

void lc_control_finish(lc_control_t *restrict const control,
                       lc_descr_t *const descr) {
  DECLARE_CCTK_PARAMETERS;

  // Finish selftest
  if (selftest) {
    assert(control->selftest_array);
#pragma omp barrier
#pragma omp master
    {
      ptrdiff_t nfailed = 0;
      for (ptrdiff_t k = 0; k < control->ash.v[2]; ++k) {
        for (ptrdiff_t j = 0; j < control->ash.v[1]; ++j) {
          for (ptrdiff_t i = 0; i < control->ash.v[0]; ++i) {
            const bool inside = i >= control->overall.min.v[0] and
                                j >= control->overall.min.v[1] and
                                k >= control->overall.min.v[2] and
                                i < control->overall.max.v[0] and
                                j < control->overall.max.v[1] and
                                k < control->overall.max.v[2];
            const ptrdiff_t ipos = ind(control->ash, i, j, k);
            nfailed += control->selftest_array[ipos] != inside;
          }
        }
      }
      if (nfailed > 0) {
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "LoopControl self-test failed");
      }
      delete[] control->selftest_array;
    }
    control->selftest_array = NULL;
  }

#pragma omp barrier
#pragma omp master
  {
    // Finish timing
    const ticks end_time = getticks();
    const double elapsed_time =
        seconds_per_tick() * elapsed(end_time, descr->start_time);
    ptrdiff_t npoints = 1;
    for (int d = 0; d < LC_DIM; ++d) {
      npoints *= control->overall.max.v[d] - control->overall.min.v[d];
    }

    // Collect statistics
    const double old_avg = descr->current_params->stats.avg_point();
    descr->current_params->stats.add(npoints, omp_get_num_threads(),
                                     elapsed_time);
    const double new_avg = descr->current_params->stats.avg_point();
    descr->current_setup->stats.add(npoints, omp_get_num_threads(),
                                    elapsed_time);
    descr->stats.add(npoints, omp_get_num_threads(), elapsed_time);
    if (veryverbose) {
      if (descr->stats.count == 0.0) {
        const double time_point =
            elapsed_time * omp_get_num_threads() / npoints;
        CCTK_VInfo(CCTK_THORNSTRING, "Loop %s: time=%g, time/point=%g s",
                   descr->name.c_str(), elapsed_time, time_point);
      } else {
        CCTK_VInfo(CCTK_THORNSTRING,
                   "Loop %s: count=%g, avg/thread=%g s, avg/point=%g s",
                   descr->name.c_str(), descr->stats.count,
                   descr->stats.avg_thread(), descr->stats.avg_point());
      }
    }

    lc_setup_t *const setup = descr->current_setup;
    if (setup->current_params == setup->best_params and new_avg > old_avg) {
      // The current best params just became worse, so forget it
      setup->best_params = NULL;
    } else if (setup->current_params != setup->best_params and
               new_avg < setup->best_params->stats.avg_point()) {
      // We found a new best params
      setup->best_params = setup->current_params;
    }
    if (not setup->best_params) {
      // We don't know which params is best, so find it
      // TODO: This is expensive -- maintain a tree instead?
      double best_avg = -1.0;
      for (lc_setup_t::params_map_t::iterator params_i = setup->params.begin(),
                                              params_end = setup->params.end();
           params_i != params_end; ++params_i) {
        lc_params_t *const params = params_i->second;
        const double avg = params->stats.avg_point();
        if (best_avg < 0.0 or avg < best_avg) {
          setup->best_params = params;
          best_avg = avg;
        }
      }
    }
    assert(setup->best_params);

    descr->current_setup = NULL;
    descr->current_params = NULL;

    // Tear down multithreading state
    delete control->coarse_thread_info_ptr;
    control->coarse_thread_info_ptr = NULL;
  }
#pragma omp barrier
}

void lc_thread_init(lc_control_t *restrict const control) {
  space_set_count(control->coarse_thread);
#pragma omp single
  { control->coarse_thread_info_ptr->idx = get_num_coarse_threads(); }
  control->coarse_thread_done =
      space_global2local(control->coarse_thread, get_coarse_thread_num());
  space_idx2pos(control->coarse_thread);
}

int lc_thread_done(const lc_control_t *restrict const control) {
  return control->coarse_thread_done;
}

void lc_thread_step(lc_control_t *restrict const control) {
  // Get next thread block
  int new_global_idx = -1;
  if (get_fine_thread_num() == 0) {
#pragma omp critical(LoopControl_lc_thread_step)
    { new_global_idx = control->coarse_thread_info_ptr->idx++; }
  }
  new_global_idx = fine_thread_broadcast(
      &lc_fine_thread_comm[get_coarse_thread_num()], new_global_idx);
  control->coarse_thread_done =
      space_global2local(control->coarse_thread, new_global_idx);
  space_idx2pos(control->coarse_thread);
}

void lc_selftest_set(const lc_control_t *restrict control, const ptrdiff_t imin,
                     const ptrdiff_t imax, const ptrdiff_t istr,
                     const ptrdiff_t i0, const ptrdiff_t j, const ptrdiff_t k) {
  DECLARE_CCTK_PARAMETERS;
  assert(selftest);
  assert(imin >= 0 and imin < imax and imax <= control->ash.v[0]);
  assert(istr > 0);
  assert(j >= 0 and j < control->ash.v[1]);
  assert(k >= 0 and k < control->ash.v[2]);
  assert(i0 + istr - 1 >= control->overall.min.v[0] and
         i0 < control->overall.max.v[0]);
  if (imin > control->overall.min.v[0]) {
    const ptrdiff_t ipos_imin = ind(control->ash, imin, j, k);
    assert(ipos_imin % istr == 0);
  }
  if (imax < control->overall.max.v[0]) {
    const ptrdiff_t ipos_imax = ind(control->ash, imax, j, k);
    assert(ipos_imax % istr == 0);
  }
  assert(j >= control->overall.min.v[1] and j < control->overall.max.v[1]);
  assert(k >= control->overall.min.v[2] and k < control->overall.max.v[2]);
  for (ptrdiff_t i = i0; i < i0 + istr; ++i) {
    if (i >= imin and i < imax) {
      assert(i >= 0 and i < control->ash.v[0]);
      assert(i >= control->overall.min.v[0] and i < control->overall.max.v[0]);
      const ptrdiff_t ipos = ind(control->ash, i, j, k);
      unsigned char &elt = control->selftest_array[ipos];
#ifdef _CRAYC
// Cray C++ compiler 8.1.2 segfaults on atomic
#pragma omp critical(lc_selftest_set)
      ++elt;
#else
#pragma omp atomic
      ++elt;
#endif
      if (elt != 1) {
#pragma omp critical
        {
          fflush(stdout);
          fprintf(stderr, "thread=%d/%d fine_thread=%d/%d ijk=[%td,%td,%td]\n",
                  get_coarse_thread_num(), get_num_coarse_threads(),
                  get_fine_thread_num(), get_num_fine_threads(), i, j, k);
          assert(0);
        }
      }
    }
  }
}

int lc_setup(void) {
  check_fortran_type_sizes();
  return 0;
}

void lc_steer(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  lc_do_settle = settle_after_iteration >= 0 and
                 cctkGH->cctk_iteration >= settle_after_iteration;

  lc_do_explore_eagerly =
      not lc_do_settle and
      cctkGH->cctk_iteration < explore_eagerly_before_iteration;
}

void lc_statistics(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  {
    CCTK_INFO("LoopControl statistics:");
    const size_t nloops = all_descrs.size();
    size_t nsetups = 0, nparams = 0;
    double time_default = 0.0, time_best = 0.0, time_actual = 0.0;
    for (all_descrs_t::const_iterator idescr = all_descrs.begin();
         idescr != all_descrs.end(); ++idescr) {
      const lc_descr_t &descr = **idescr;
      nsetups += descr.setups.size();
      for (lc_descr_t::setup_map_t::const_iterator
               setup_i = descr.setups.begin(),
               setup_end = descr.setups.end();
           setup_i != setup_end; ++setup_i) {
        const lc_setup_t &setup = *setup_i->second;
        nparams += setup.params.size();
        const double setup_count = setup.stats.count * setup.stats.points;
        time_default += setup_count * setup.default_params->stats.avg_point();
        time_best += setup_count * setup.best_params->stats.avg_point();
        time_actual += setup_count * setup.stats.avg_point();
      }
    }
    const double ratio_unopt = time_default == 0.0 && time_actual == 0.0
                                   ? 0.0
                                   : time_default / time_actual - 1.0;
    const double ratio_ideal = time_default == 0.0 && time_actual == 0.0
                                   ? 0.0
                                   : time_best / time_actual - 1.0;
    const size_t nbytes =
        (nloops * (sizeof(lc_descr_t *) + sizeof(lc_descr_t))) +
        (nsetups *
         (sizeof(lc_setup_key_t) + sizeof(lc_setup_t *) + sizeof(lc_setup_t))) +
        (nparams * (sizeof(lc_params_key_t) + sizeof(lc_params_t *) +
                    sizeof(lc_params_t)));
    CCTK_VInfo(CCTK_THORNSTRING, "  Loops traversed:    %td", nloops);
    CCTK_VInfo(CCTK_THORNSTRING, "  Setups encountered: %td", nsetups);
    CCTK_VInfo(CCTK_THORNSTRING, "  Params explored:    %td", nparams);
    CCTK_VInfo(CCTK_THORNSTRING, "    Actual time spent:                %g s",
               time_actual);
    CCTK_VInfo(CCTK_THORNSTRING,
               "    Unoptimized time would have been: %g s   (%+.1f%%)",
               time_default, 100.0 * ratio_unopt);
    CCTK_VInfo(CCTK_THORNSTRING,
               "    Ideal time could have been:       %g s   (%+.1f%%)",
               time_best, 100.0 * ratio_ideal);
    CCTK_VInfo(CCTK_THORNSTRING, "  Memory allocated: %g MB", nbytes / 1.0e+6);
  }

  if (strcmp(statistics_filename, "") == 0)
    return;

  static bool did_truncate = false;
  const bool do_truncate = IO_TruncateOutputFiles(cctkGH);
  const char *const mode = do_truncate and not did_truncate ? "w" : "a";
  did_truncate = true;

  char filename[10000];
  snprintf(filename, sizeof filename, "%s/%s.%06d.txt", out_dir,
           statistics_filename, CCTK_MyProc(cctkGH));
  FILE *const descrfile = fopen(filename, mode);

  fprintf(descrfile, "LoopControl statistics:\n");
  vector<lc_descr_t *> all_descrs_sorted;
  all_descrs_sorted = all_descrs;
  sort(all_descrs_sorted.begin(), all_descrs_sorted.end(), descr_comp_name());
  for (vector<lc_descr_t *>::const_iterator descr_i = all_descrs_sorted.begin(),
                                            descr_e = all_descrs_sorted.end();
       descr_i != descr_e; ++descr_i) {
    const lc_descr_t &descr = **descr_i;
    fprintf(descrfile, "   Loop %s (%s:%d):\n", descr.name.c_str(),
            descr.file.c_str(), descr.line);
    // TODO: sort setups?
    for (lc_descr_t::setup_map_t::const_iterator setup_i = descr.setups.begin(),
                                                 setup_end = descr.setups.end();
         setup_i != setup_end; ++setup_i) {
      const lc_setup_t &setup = *setup_i->second;
      fprintf(descrfile,
              "      setup=[%d,%d,%d]:[%d,%d,%d]/[%d,%d,%d] nt=%d/%d\n",
              int(setup.key.min.v[0]), int(setup.key.min.v[1]),
              int(setup.key.min.v[2]), int(setup.key.max.v[0]),
              int(setup.key.max.v[1]), int(setup.key.max.v[2]),
              int(setup.key.ash.v[0]), int(setup.key.ash.v[1]),
              int(setup.key.ash.v[2]), int(setup.key.num_coarse_threads),
              int(setup.key.num_fine_threads));
      double best_avg = numeric_limits<double>::max(), worst_avg = 0.0;
      vector<lc_params_t *> params_sorted;
      params_sorted.reserve(setup.params.size());
      for (lc_setup_t::params_map_t::const_iterator
               params_i = setup.params.begin(),
               params_e = setup.params.end();
           params_i != params_e; ++params_i) {
        params_sorted.push_back(params_i->second);
      }
      sort(params_sorted.begin(), params_sorted.end(), params_comp_time());
      for (vector<lc_params_t *>::const_iterator
               params_i = params_sorted.begin(),
               params_e = params_sorted.end();
           params_i != params_e; ++params_i) {
        const lc_params_t &params = **params_i;
        fprintf(descrfile, "         tilesize=[%d,%d,%d] loopsize=[%d,%d,%d]\n",
                int(params.key.tilesize.v[0]), int(params.key.tilesize.v[1]),
                int(params.key.tilesize.v[2]), int(params.key.loopsize.v[0]),
                int(params.key.loopsize.v[1]), int(params.key.loopsize.v[2]));
        const lc_stats_t &stats = params.stats;
        fprintf(descrfile,
                "            count=%g, avg/thread=%g s, avg/point=%g s%s%s\n",
                stats.count, stats.avg_thread(), stats.avg_point(),
                &params == setup.default_params ? " (DEFAULT)" : "",
                &params == setup.best_params ? " (BEST)" : "");
        best_avg = min(best_avg, stats.avg_point());
        worst_avg = max(worst_avg, stats.avg_point());
      }
      const double default_avg = setup.default_params->stats.avg_point();
      fprintf(descrfile, "         best(avg/point)=%g s, worst(avg/point)=%g "
                         "s, default(avg/point)=%g s\n",
              best_avg, worst_avg, default_avg);
      const lc_stats_t &stats = setup.stats;
      fprintf(descrfile, "         count=%g, avg/thread=%g s, avg/point=%g s\n",
              stats.count, stats.avg_thread(), stats.avg_point());
    }
    const lc_stats_t &stats = descr.stats;
    fprintf(descrfile, "      count=%g, avg/thread=%g s, avg/point=%g s\n",
            stats.count, stats.avg_thread(), stats.avg_point());
  }
  fprintf(descrfile, "\n");
  fclose(descrfile);
}

void lc_statistics_analysis(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static double last_output = 0.0;
  const double run_time = CCTK_RunTime();

  if (veryverbose || (statistics_every_seconds >= 0.0 &&
                      run_time >= last_output + statistics_every_seconds)) {
    lc_statistics(CCTK_PASS_CTOC);
    last_output = run_time;
  }
}

void lc_statistics_terminate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (verbose or veryverbose or statistics_every_seconds >= 0.0) {
    lc_statistics(CCTK_PASS_CTOC);
  }
}

extern "C" CCTK_FCALL void CCTK_FNAME(lc_descr_init)(CCTK_POINTER &descr,
                                                     int &line,
                                                     TWO_FORTSTRINGS_ARGS) {
  TWO_FORTSTRINGS_CREATE(file, name);
  lc_descr_init((lc_descr_t **)&descr, name, file, line);
  free(name);
  free(file);
}

extern "C" CCTK_FCALL void CCTK_FNAME(lc_control_init)(
    lc_control_t &restrict control, CCTK_POINTER &descr, const int &imin,
    const int &jmin, const int &kmin, const int &imax, const int &jmax,
    const int &kmax, const int &iash, const int &jash, const int &kash,
    const int &istr) {
  lc_control_init(&control, (lc_descr_t *)descr, imin, jmin, kmin, imax, jmax,
                  kmax, iash, jash, kash, istr);
}

extern "C" CCTK_FCALL void
    CCTK_FNAME(lc_control_finish)(lc_control_t &restrict control,
                                  CCTK_POINTER &descr) {
  lc_control_finish(&control, (lc_descr_t *)descr);
}

extern "C" CCTK_FCALL void CCTK_FNAME(lc_thread_init)(lc_control_t &control) {
  lc_thread_init(&control);
}

extern "C" CCTK_FCALL int
    CCTK_FNAME(lc_thread_done)(const lc_control_t &control) {
  return lc_thread_done(&control);
}

extern "C" CCTK_FCALL void CCTK_FNAME(lc_thread_step)(lc_control_t &control) {
  lc_thread_step(&control);
}
