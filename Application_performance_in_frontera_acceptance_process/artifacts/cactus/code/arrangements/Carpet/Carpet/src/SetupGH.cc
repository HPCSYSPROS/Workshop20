#define _GNU_SOURCE                                                            \
  1 // needed for sched_getaffinity, best at the top to avoid inconsistent
    // includes

#include <algorithm>
#include <cassert>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stack>
#include <string>
#include <vector>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_ErrorCodes.h>
#include <util_Network.h>
#include <util_Table.h>

#ifdef HAVE_SCHED_H
#include <sched.h>
#endif
#include <unistd.h>

#include <Requirements.hh>

#include <Timer.hh>

#include <bbox.hh>
#include <defs.hh>
#include <dist.hh>
#include <ggf.hh>
#include <gh.hh>
#include <mpi_string.hh>
#include <region.hh>
#include <vect.hh>

#include <carpet.hh>

namespace Carpet {

using namespace std;

static void setup_model_information(cGH *cctkGH);
static void setup_multigrid_information(cGH *cctkGH);
static void setup_refinement_information();
static void setup_map_information();
static void setup_time_information();
static void setup_domain_extents(cGH const *cctkGH);
static void allocate_grid_hierarchy(cGH const *cctkGH, int m);
static void allocate_data_hierarchy(cGH const *cctkGH, int m);
static void allocate_time_hierarchy(cGH const *cctkGH);
static void setup_grid_hierarchy(cGH const *cctkGH);
static void set_base_extent(int m, vector<vector<region_t> > &regss);

static void allocate_group_data(cGH const *cctkGH);
static void allocate_group_hierarchies(int group, ivect const &sizes,
                                       vector<i2vect> const &ghosts);
static void setup_group_grid_hierarchy(cGH const *cctkGH, int group,
                                       cGroup const &gdata,
                                       ivect const &convpowers,
                                       ivect const &convoffsets);
static void initialise_group_info(cGH const *cctkGH, int group,
                                  cGroup const &gdata);
static void set_state(cGH *cctkGH);
static void enable_storage_for_all_groups(cGH const *cctkGH);

static vector<i2vect> get_ghostzones();
static vector<int> get_prolongation_orders_space();
static ivect get_npoints();
static void get_boundary_specification(cGH const *cctkGH, int m,
                                       vector<i2vect> const &ghosts,
                                       i2vect &nboundaryzones,
                                       b2vect &is_internal,
                                       b2vect &is_staggered, i2vect &shiftout);
static void get_domain_specification(cGH const *cctkGH, int m,
                                     ivect const &npoints, rvect &physical_min,
                                     rvect &physical_max, rvect &base_spacing);
static void adapt_domain_specification(int m, rvect const &physical_min,
                                       rvect const &physical_max,
                                       rvect const &base_spacing,
                                       rvect &exterior_min, rvect &exterior_max,
                                       rvect &spacing);
static void calculate_grid_points(int m, vector<i2vect> const &ghosts,
                                  rvect const &exterior_min,
                                  rvect const &exterior_max,
                                  rvect const &spacing, ivect &npoints);
static void
calculate_base_extents(ibbox const &baseextent, centering refcentering,
                       centering mgcentering, i2vect const &nboundaryzones,
                       b2vect const &is_internal, b2vect const &is_staggered,
                       i2vect const &shiftout,
                       vector<vector<ibbox> > &baseextents);
static void find_processor_decomposition(
    cGH const *cctkGH, vector<vector<vector<region_t> > > &superregsss,
    vector<vector<vector<vector<region_t> > > > &regssss);

static void get_group_size(int group, cGroup const &gdata, ivect &sizes,
                           vector<i2vect> &ghosts);
static void adapt_group_size_mglevel(int group, cGroup const &gdata,
                                     ivect &sizes, ivect &convpowers,
                                     ivect &convoffsets);
static void get_convergence_options(int group, cGroup const &gdata,
                                    ivect &convpowers, ivect &convoffsets);
static void adapt_group_size_disttype(cGH const *cctkGH, int group,
                                      cGroup const &gdata, ivect &sizes,
                                      vector<i2vect> const &ghosts);
static void output_group_statistics(cGH const *cctkGH);
static operator_type get_transport_operator(cGH const *cctkGH, int group,
                                            cGroup const &gdata);
static bool can_transfer_variable_type(cGH const *cctkGH, int group,
                                       cGroup const &gdata);

static void ensure_CartGrid3D_type();
#if 0
  static void
  ensure_CartGrid3D_domain ();  // UNUSED
#endif
static void ensure_CartGrid3D_avoid_origin();
static void ensure_ReflectionSymmetry_avoid_origin(centering refcentering);
static void ensure_ghostzones(int m, vector<i2vect> const &ghosts);
static void ensure_group_options(int group, cGroup const &gdata);

void *SetupGH(tFleshConfig *const fc, int const convLevel, cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // Some statistics
  {
    int const nprocs = dist::size();
    int const myproc = dist::rank();
    dist::set_num_threads(num_threads);
    int const mynthreads = dist::num_threads();
    int const nthreads_total = dist::total_num_threads();
    // Can't store this in a variable because the value is different
    // for each thread
    // int const mythreadnum = dist::thread_num();
    char const *const CACTUS_NUM_PROCS = getenv("CACTUS_NUM_PROCS");
    int const cactus_num_procs = CACTUS_NUM_PROCS ? atoi(CACTUS_NUM_PROCS) : 0;
#ifdef CCTK_MPI
    CCTK_VInfo(CCTK_THORNSTRING, "MPI is enabled");
    CCTK_VInfo(CCTK_THORNSTRING, "Carpet is running on %d processes", nprocs);
    if (not CACTUS_NUM_PROCS) {
      CCTK_VWarn(CCTK_WARN_COMPLAIN, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Although MPI is enabled, the environment variable "
                 "CACTUS_NUM_PROCS is not set.");
    } else {
      if (cactus_num_procs != nprocs) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "The environment variable CACTUS_NUM_PROCS is set to %d, "
                   "but there are %d MPI processes. This may indicate a severe "
                   "problem with the MPI startup mechanism.",
                   cactus_num_procs, nprocs);
      }
    }
    CCTK_VInfo(CCTK_THORNSTRING, "This is process %d", myproc);
#else
    CCTK_VInfo(CCTK_THORNSTRING, "MPI is disabled");
    if (CACTUS_NUM_PROCS and cactus_num_procs != nprocs) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Although MPI is disabled, the environment variable "
                 "CACTUS_NUM_PROCS is set to %d. This may indicate a severe "
                 "problem with the Cactus startup mechanism.",
                 cactus_num_procs);
    }
#endif
    char const *const OMP_NUM_THREADS = getenv("OMP_NUM_THREADS");
    char const *const CACTUS_NUM_THREADS = getenv("CACTUS_NUM_THREADS");
    int const cactus_num_threads =
        CACTUS_NUM_THREADS ? atoi(CACTUS_NUM_THREADS) : 0;
#ifdef _OPENMP
    CCTK_VInfo(CCTK_THORNSTRING, "OpenMP is enabled");
    if (not OMP_NUM_THREADS and num_threads == -1) {
      CCTK_VWarn(CCTK_WARN_COMPLAIN, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Although OpenMP is enabled, neither the environment variable "
                 "OMP_NUM_THREADS nor the parameter Carpet::num_threads are "
                 "set. A system-specific default value is used instead.");
    }
    CCTK_VInfo(CCTK_THORNSTRING,
               "This process contains %d threads, this is thread %d",
               mynthreads, dist::thread_num());
    if (not CACTUS_NUM_THREADS) {
      CCTK_VWarn(CCTK_WARN_COMPLAIN, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Although OpenMP is enabled, the environment variable "
                 "CACTUS_NUM_THREADS is not set.");
    } else {
      if (cactus_num_threads != mynthreads) {
        CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "The environment variable CACTUS_NUM_THREADS is set to %d, "
                   "but there are %d threads on this process. This may "
                   "indicate a severe problem with the OpenMP startup "
                   "mechanism.",
                   cactus_num_threads, mynthreads);
      }
    }
    CCTK_VInfo(CCTK_THORNSTRING, "There are %d threads in total",
               nthreads_total);
#ifdef CCTK_MPI
    CCTK_VInfo(CCTK_THORNSTRING, "There are %g threads per process",
               1.0 * nthreads_total / nprocs);
#endif
#else
    CCTK_VInfo(CCTK_THORNSTRING, "OpenMP is disabled");
    int const omp_num_threads = OMP_NUM_THREADS ? atoi(OMP_NUM_THREADS) : 0;
    if (omp_num_threads > 0) {
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Although OpenMP is disabled, the environment variable "
                 "OMP_NUM_THREADS is set to %d. It will be ignored.",
                 omp_num_threads);
    }
    if (num_threads > 0) {
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Although OpenMP is disabled, the parameter "
                 "Carpet::num_threads is set to %d. It will be ignored.",
                 num_threads);
    }
    if (CACTUS_NUM_THREADS and cactus_num_threads != mynthreads) {
      CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Although OpenMP is disabled, the environment variable "
                 "CACTUS_NUM_THREADS is set to %d. This may indicate a severe "
                 "problem with the Cactus startup mechanism.",
                 cactus_num_threads);
    }
#endif

    char hostnamebuf[MPI_MAX_PROCESSOR_NAME];
    int hostnamelen;
    MPI_Get_processor_name(hostnamebuf, &hostnamelen);
    string const hostname(hostnamebuf);
    DetermineHosts(hostname, verbose);
#if HAVE_GETPID
    int const mypid = static_cast<int>(getpid());
#else
    int const mypid = -1;
#endif
    // Output
    CCTK_VInfo(CCTK_THORNSTRING, "This process runs on host %s, pid=%d",
               hostname.c_str(), mypid);
    if (verbose or veryverbose) {
      // Collect process ids
      vector<int> pids(nprocs);
      MPI_Allgather(const_cast<int *>(&mypid), 1, MPI_INT, &pids.front(), 1,
                    MPI_INT, dist::comm());
      // Collect number of threads
      vector<int> nthreads(nprocs);
      MPI_Allgather(const_cast<int *>(&mynthreads), 1, MPI_INT,
                    &nthreads.front(), 1, MPI_INT, dist::comm());
      // Output
      CCTK_VInfo(CCTK_THORNSTRING, "Running on the following hosts:");
      for (int n = 0; n < nprocs; ++n) {
        int const host_id = HostId(n);
        CCTK_VInfo(CCTK_THORNSTRING,
                   "   %6d: hid=%d (%s), pid=%d, num_threads=%d", n, host_id,
                   HostName(host_id).c_str(), pids.AT(n), nthreads.AT(n));
      }
    }

    if (set_cpu_affinity) {
#ifdef HAVE_SCHED_GETAFFINITY
      // Cores available to this process
      assert(mynthreads <= CPU_SETSIZE);
      vector<bool> mask(CPU_SETSIZE, false);
#pragma omp parallel
      {
        cpu_set_t cpumask;
        int const ierr = sched_getaffinity(0, sizeof cpumask, &cpumask);
        assert(not ierr);
        for (int n = 0; n < CPU_SETSIZE; ++n) {
          if (CPU_ISSET(n, &cpumask)) {
#pragma omp critical
            mask.at(n) = true;
          }
        }
      }
      int const host_id = HostId(myproc);
      vector<int> const host_procs = HostProcs(host_id);
      int const num_host_procs = host_procs.size();
      assert(num_host_procs > 0);
#if 0
        // Collect information from all processes on this host
        vector<char> sendbuf(num_host_procs * CPU_SETSIZE, 0);
        vector<char> recvbuf(num_host_procs * CPU_SETSIZE);
        for (int i=0; i<num_host_procs; ++i) {
          for (int n=0; n<CPU_SETSIZE; ++n) {
            sendbuf.at(i * CPU_SETSIZE + n) = mask.at(n);
          }
        }
        vector<int> sendcount(nprocs, 0);
        vector<int> recvcount(nprocs, 0);
        for (int i=0; i<num_host_procs; ++i) {
          int const p = host_procs.at(i);
          sendcount.at(p) = CPU_SETSIZE;
          recvcount.at(p) = CPU_SETSIZE;
        }
        vector<int> senddispls(nprocs);
        vector<int> recvdispls(nprocs);
        senddispls.at(0) = 0;
        recvdispls.at(0) = 0;
        for (int p=1; p<nprocs; ++p) {
          senddispls.at(p) = senddispls.at(p-1) + sendcount.at(p-1);
          recvdispls.at(p) = recvdispls.at(p-1) + recvcount.at(p-1);
        }
        MPI_Alltoallv(&sendbuf[0], &sendcounts[0], &senddispls[0], MPI_CHAR,
                      &recvbuf[0], &recvcounts[0], &recvdispls[0], MPI_CHAR,
                      MPI_COMM_WORLD);
        // Cores available on this host
        vector<bool> othermasks(CPU_SETSIZE, false);
        for (int i=0; i<num_host_procs; ++i) {
          int const p = host_procs.at(i);
          if (p != myproc) {
            for (int n=0; n<CPU_SETSIZE; ++n) {
              othermasks.at(n) =
                othermasks.at(n) or recvbuf.at(i * CPU_SETSIZE + n);
            }
          }
        }
#endif
      int num_cores = 0;
      for (int n = 0; n < CPU_SETSIZE; ++n) {
        if (mask.at(n))
          ++num_cores;
      }
      if (num_cores == 0) {
        CCTK_WARN(CCTK_WARN_ALERT, "Cannot select core set for this process "
                                   "(no cores seem available to this process)");
      } else {
        bool const other_procs_need_cores =
            num_cores >= mynthreads * num_host_procs;
        int my_host_index = 0;
        if (other_procs_need_cores) {
          // It seems that all cores are available to all processes
          // on this host. Split the cores evenly across the
          // processes.
          for (int i = 0; i < num_host_procs; ++i) {
            int const p = host_procs.at(i);
            if (p == myproc) {
              my_host_index = i;
              break;
            }
          }
          if (num_cores % num_host_procs != 0) {
            CCTK_WARN(CCTK_WARN_ALERT, "The number of processes on this host "
                                       "does not evenly divide the number of "
                                       "cores on this host -- leaving some "
                                       "cores unassigned");
          }
          num_cores /= num_host_procs;
        }
        if (num_cores % mynthreads != 0) {
          CCTK_WARN(CCTK_WARN_ALERT, "The number of threads in this process "
                                     "does not evenly divide the number of "
                                     "cores for this process -- leaving some "
                                     "cores unassigned");
        }
        int const num_cores_per_thread = num_cores / mynthreads;
        // Choose cores for this process
        int const first_core_index =
            my_host_index * mynthreads * num_cores_per_thread;
        vector<int> available_cores;
        for (int n = 0; n < CPU_SETSIZE; ++n) {
          if (mask.at(n))
            available_cores.push_back(n);
        }
        int const n0 = available_cores.at(first_core_index);
        CCTK_VInfo(CCTK_THORNSTRING,
                   "Selecting cores %d (and following) for this process", n0);
#pragma omp parallel
        {
          cpu_set_t cpumask;
          CPU_ZERO(&cpumask);
          int const i0 =
              first_core_index + dist::thread_num() * num_cores_per_thread;
          for (int i = 0; i < num_cores_per_thread; ++i) {
            CPU_SET(available_cores[i0 + i], &cpumask);
          }
          int const ierr = sched_setaffinity(0, sizeof(cpumask), &cpumask);
          assert(not ierr);
        }
      }
#else
      CCTK_WARN(CCTK_WARN_ALERT, "Cannot select core set for this process; "
                                 "sched_getaffinity is not available on this "
                                 "system");
#endif
    }

#ifdef HAVE_SCHED_GETAFFINITY
    {
      vector<bool> mask(CPU_SETSIZE, false);
#pragma omp parallel
      {
        cpu_set_t cpumask;
        int const ierr = sched_getaffinity(0, sizeof cpumask, &cpumask);
        assert(not ierr);
        for (int n = 0; n < CPU_SETSIZE; ++n) {
          if (CPU_ISSET(n, &cpumask)) {
#pragma omp critical
            mask.at(n) = true;
          }
        }
      }

      ostringstream buf;
      int num_cores = 0;
      bool isfirst = true;
      int first_active = -1;
      for (int n = 0; n < CPU_SETSIZE; ++n) {
        if (mask.at(n))
          ++num_cores;
        if (first_active == -1 and mask.at(n)) {
          if (not isfirst)
            buf << ", ";
          isfirst = false;
          buf << n;
          first_active = n;
        } else if (first_active >= 0 and not mask.at(n)) {
          if (n - 1 > first_active)
            buf << "-" << n - 1;
          first_active = -1;
        }
      }
      CCTK_VInfo(CCTK_THORNSTRING, "This process runs on %d core%s: %s",
                 num_cores, num_cores == 1 ? "" : "s", buf.str().c_str());
      if (mynthreads > num_cores) {
        CCTK_WARN(CCTK_WARN_ALERT, "The number of threads for this process is "
                                   "larger its number of cores. This may "
                                   "indicate a performance problem.");
      }
    }

    for (int thread = 0; thread < mynthreads; ++thread) {
      vector<bool> mask(CPU_SETSIZE, false);
      cpu_set_t cpumask;
#pragma omp parallel
      if (thread == dist::thread_num()) {
        int const ierr = sched_getaffinity(0, sizeof cpumask, &cpumask);
        assert(not ierr);
      }
      for (int n = 0; n < CPU_SETSIZE; ++n) {
        if (CPU_ISSET(n, &cpumask)) {
          mask.at(n) = true;
        }
      }

      ostringstream buf;
      int num_cores = 0;
      bool isfirst = true;
      int first_active = -1;
      for (int n = 0; n < CPU_SETSIZE; ++n) {
        if (mask.at(n))
          ++num_cores;
        if (first_active == -1 and mask.at(n)) {
          if (not isfirst)
            buf << ", ";
          isfirst = false;
          buf << n;
          first_active = n;
        } else if (first_active >= 0 and not mask.at(n)) {
          if (n - 1 > first_active)
            buf << "-" << n - 1;
          first_active = -1;
        }
      }
      CCTK_VInfo(CCTK_THORNSTRING, "Thread %d runs on %d core%s: %s", thread,
                 num_cores, num_cores == 1 ? "" : "s", buf.str().c_str());
    }
#else
    CCTK_INFO("Cannot determine core affinity");
#endif
  }

  // Initialise current position (must be the very first thing,
  // before the first output)
  mglevel = -1;
  reflevel = -1;
  mc_grouptype = -1;
  map = -1;
  component = -1;
  local_component = -1;
  timelevel = -1;
  timelevel_offset = -1;
#ifdef HAVE_CGH_CCTK_MODE
  cctkGH->cctk_mode = CCTK_MODE_META;
#endif

  // Say hello

  Timers::Timer timer("CarpetStartup");
  timer.start();

  Waypoint("Setting up the grid hierarchy");

  // Check arguments:
  CCTK_VInfo(CCTK_THORNSTRING, "This simulation is running in %d dimensions",
             cctkGH->cctk_dim);
  // Only a specific number of dimensions is supported
  assert(cctkGH->cctk_dim == dim);
  // Not sure what to do with that
  assert(convLevel == 0);

  // Set up descriptions from user parameters
  setup_model_information(cctkGH);
  setup_multigrid_information(cctkGH);
  setup_refinement_information();
  setup_map_information();
  setup_time_information();

  // Calculate domain extents for each map
  setup_domain_extents(cctkGH);

  // Set up grid hierarchy
  setup_grid_hierarchy(cctkGH);

  // Allocate space for group descriptions
  allocate_group_data(cctkGH);

  // Set times, origin, and spacings, and go to meta mode
  set_state(cctkGH);

  // Enable prolongating
  do_allow_past_timelevels = true;
  do_prolongate = true;
  do_taper = false;
  do_warn_about_storage = false; // This is enabled later
  in_analysis_bin = false;

  if (enable_all_storage) {
    if (not enable_no_storage) {
      enable_storage_for_all_groups(cctkGH);
    }
  }

  Waypoint("Done with setting up the grid hierarchy");
  timer.stop();

  return &carpetGH;
}

//////////////////////////////////////////////////////////////////////////////
// Routines which perform actions ////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void setup_model_information(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  cctkGH->identity = strdup(model);
}

void setup_multigrid_information(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  basemglevel = convergence_level;
  mglevels = num_convergence_levels;
  mgfact = convergence_factor;
  maxmglevelfact = ipow(mgfact, mglevels - 1);
  cctkGH->cctk_convfac = mgfact;
}

void setup_refinement_information() {
  DECLARE_CCTK_PARAMETERS;

  // Set maximum number of refinement levels
  maxreflevels = max_refinement_levels;

  // Set the per-level spatial refinement factors
  if (CCTK_EQUALS(space_refinement_factors, "")) {
    // Calculate them from the default refinement factor
    spacereffacts.resize(maxreflevels);
    assert(ipow(2, 0) == 1);
    assert(ipow(3, 1) == 3);
    assert(ipow(4, 2) == 16);
    for (int rl = 0; rl < maxreflevels; ++rl) {
      spacereffacts.AT(rl) = ivect(ipow(int(refinement_factor), rl));
    }
  } else {
    // Read them from the parameter
    try {
      istringstream srf(space_refinement_factors);
      srf >> spacereffacts;
    } catch (input_error) {
      CCTK_ERROR("Could not parse parameter \"space_refinement_factors\"");
    }
  }
  // TODO: turn these into real error messages
  assert(int(spacereffacts.size()) >= maxreflevels);
  assert(all(spacereffacts.AT(0) == 1));
  for (int rl = 1; rl < maxreflevels; ++rl) {
    assert(all(spacereffacts.AT(rl) >= spacereffacts.AT(rl - 1)));
    assert(all(spacereffacts.AT(rl) % spacereffacts.AT(rl - 1) == 0));
  }
  spacereffacts.resize(maxreflevels);

  // Set the per-level temporal refinement factors
  if (CCTK_EQUALS(time_refinement_factors, "")) {
    // Calculate them from the default refinement factor
    timereffacts.resize(maxreflevels);
    for (int rl = 0; rl < maxreflevels; ++rl) {
      timereffacts.AT(rl) = ipow(int(refinement_factor), rl);
    }
  } else {
    // Read them from the parameter
    try {
      istringstream trf(time_refinement_factors);
      trf >> timereffacts;
    } catch (input_error) {
      CCTK_ERROR("Could not parse parameter \"time_refinement_factors\"");
    }
  }
  // TODO: turn these into real error messages
  assert(int(timereffacts.size()) >= maxreflevels);
  assert(timereffacts.AT(0) == 1);
  for (int rl = 1; rl < maxreflevels; ++rl) {
    assert(timereffacts.AT(rl) >= timereffacts.AT(rl - 1));
    assert(timereffacts.AT(rl) % timereffacts.AT(rl - 1) == 0);
  }
  timereffacts.resize(maxreflevels);

  // Calculate the maximum refinement factors
  maxtimereflevelfact = timereffacts.AT(maxreflevels - 1);
  maxspacereflevelfact = spacereffacts.AT(maxreflevels - 1);
}

void setup_map_information() {
  DECLARE_CCTK_PARAMETERS;

  if (domain_from_multipatch) {
    assert(num_maps == 1); // must be the default to avoid confusion
    assert(CCTK_IsFunctionAliased("MultiPatch_GetSystemSpecification"));
    CCTK_INT maps1;
    check(not MultiPatch_GetSystemSpecification(&maps1));
    maps = maps1;
  } else {
    maps = num_maps;
  }
  carpetGH.maps = maps;
#ifdef REQUIREMENTS_HH
  Requirements::Setup(maps);
#endif
}

void setup_time_information() {
  DECLARE_CCTK_PARAMETERS;

  // Set maximum number of time levels
  if (max_timelevels < 0) {
    // Set automatically (backward compatibility)
    maxtimelevels = prolongation_order_time + 1;
  } else {
    maxtimelevels = max_timelevels;
  }
  if (maxtimelevels < prolongation_order_time + 1) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "There are not enough time levels for this time prolongation "
                "order: max_timelevels=%d, prolongation_order_time=%d",
                int(max_timelevels), int(prolongation_order_time));
  }
}

void setup_domain_extents(cGH const *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  for (int m = 0; m < maps; ++m) {

    // Allocate hierarchies
    allocate_grid_hierarchy(cctkGH, m);
    allocate_data_hierarchy(cctkGH, m);

  } // for m

  allocate_time_hierarchy(cctkGH);
}

void allocate_grid_hierarchy(cGH const *const cctkGH, int const m) {
  DECLARE_CCTK_PARAMETERS;

  // Centering
  centering refcentering;
  int reffactdenom;
  if (CCTK_EQUALS(refinement_centering, "vertex")) {
    refcentering = vertex_centered;
    reffactdenom = 1;
  } else if (CCTK_EQUALS(refinement_centering, "cell")) {
    refcentering = cell_centered;
    reffactdenom = 2;
  } else {
    assert(0);
  }

  // Number of grid points
  ivect npoints = get_npoints();

  // Number of ghost zones
  vector<i2vect> const ghosts = get_ghostzones();

  // Boundary description
  i2vect nboundaryzones;
  b2vect is_internal;
  b2vect is_staggered;
  i2vect shiftout;
  get_boundary_specification(cctkGH, m, ghosts, nboundaryzones, is_internal,
                             is_staggered, shiftout);

  // Grid size
  rvect physical_min, physical_max;
  rvect base_spacing;
  get_domain_specification(cctkGH, m, npoints, physical_min, physical_max,
                           base_spacing);

  if (maxreflevels > 1) {
    // Ensure that ReflectionSymmetry::avoid_origin is set correctly
    ensure_ReflectionSymmetry_avoid_origin(refcentering);
  }

  // Adapt domain specification for convergence level
  rvect exterior_min, exterior_max;
  rvect spacing;
  adapt_domain_specification(m, physical_min, physical_max, base_spacing,
                             exterior_min, exterior_max, spacing);

  // Calculate global number of grid points
  calculate_grid_points(m, ghosts, exterior_min, exterior_max, spacing,
                        npoints);

  centering const mgcentering = refcentering;

  // Base grid extent
  ivect const str(maxspacereflevelfact * reffactdenom);
  ivect const lb(0);
  ivect const ub((npoints - 1) * str);
  ibbox const baseext(lb, ub, str);

  vector<vector<ibbox> > baseexts;
  calculate_base_extents(baseext, refcentering, mgcentering, nboundaryzones,
                         is_internal, is_staggered, shiftout, baseexts);

  // Allocate grid hierarchy

  Timers::Timer timer("AllocateGridHierarchy");
  timer.start();
  vhh.resize(maps);
  vhh.AT(m) = new gh(spacereffacts, refcentering, convergence_factor,
                     mgcentering, baseexts, nboundaryzones);
  timer.stop();
}

void allocate_data_hierarchy(cGH const *const cctkGH, int const m) {
  DECLARE_CCTK_PARAMETERS;

  // Number of ghost zones
  vector<i2vect> const ghosts = get_ghostzones();

  int buffer_factor;
  if (use_buffer_zones) {
    if (num_integrator_substeps == -1) {
      assert(CCTK_IsFunctionAliased("MoLNumIntegratorSubsteps"));
      buffer_factor = MoLNumIntegratorSubsteps();
    } else {
      buffer_factor = num_integrator_substeps;
    }
  } else {
    buffer_factor = 1;
  }
  assert(buffer_factor > 0);

  int const taper_factor = use_tapered_grids ? refinement_factor : 1;
  for (int rl = 0; rl < maxreflevels; ++rl) {
    assert(all(all(
        buffer_factor * ghosts.AT(rl) + int(additional_buffer_zones) >= 0)));
  }

  const streamsize oldprecision = cout.precision();
  const ios_base::fmtflags oldflags = cout.flags();
  cout.setf(ios::fixed);
  CCTK_INFO("Buffer zone counts (excluding ghosts):");
  vector<i2vect> buffers(maxreflevels);
  for (int rl = 0; rl < maxreflevels; ++rl) {
    buffers.AT(rl) =
        rl == 0 ? i2vect(0) : taper_factor * (buffer_factor * ghosts.AT(rl) +
                                              int(additional_buffer_zones)) -
                                  ghosts.AT(rl);
    cout << "   [" << rl << "]: " << buffers.AT(rl) << "\n";
    assert(all(all(buffers.AT(rl) >= 0)));
  }
  CCTK_INFO("Overlap zone counts:");
  vector<i2vect> overlaps(maxreflevels);
  for (int rl = 0; rl < maxreflevels; ++rl) {
    gh const &hh = *vhh.AT(m);
    overlaps.AT(rl) =
        rl == 0 ? i2vect(0)
                : (use_overlap_zones
                       ? hh.reffacts.AT(rl) / hh.reffacts.AT(rl - 1) *
                             (ghosts.AT(rl) + int(additional_overlap_zones))
                       : i2vect(0));
    cout << "   [" << rl << "]: " << overlaps.AT(rl) << "\n";
    assert(all(all(overlaps.AT(rl) >= 0)));
  }
  cout.precision(oldprecision);
  cout.setf(oldflags);

  vector<int> const my_prolongation_orders_space =
      get_prolongation_orders_space();

  vdd.resize(maps);
  vdd.AT(m) = new dh(*vhh.AT(m), ghosts, buffers, overlaps,
                     my_prolongation_orders_space);

  if (maxreflevels > 1) {
    ensure_ghostzones(m, ghosts);
  }
}

void allocate_time_hierarchy(cGH const *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // We are using the gh of the first map.  This works because all
  // maps have the same number of refinement levels.
  tt = new th(*vhh.AT(0), maxtimelevels, timereffacts,
              time_interpolation_during_regridding);
}

void setup_grid_hierarchy(cGH const *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  vector<vector<vector<region_t> > > superregsss(maps);
  for (int m = 0; m < maps; ++m) {
    set_base_extent(m, superregsss.AT(m));
  }

  vector<vector<vector<vector<region_t> > > > regssss;
  find_processor_decomposition(cctkGH, superregsss, regssss);

  for (int m = 0; m < maps; ++m) {

    // Check the regions
    CheckRegions(regssss.AT(m));

    // Recompose grid hierarchy
    vhh.AT(m)->regrid(superregsss.AT(m), regssss.AT(m), false);
    int const rl = 0;
    vhh.AT(m)->recompose(rl, false);
    vhh.AT(m)->regrid_free(false);

    // Output grid structure to screen but not to files since
    // IO_TruncateFiles might not be present yet
    OutputGrids(cctkGH, m, *vhh.AT(m), *vdd.AT(m));

  } // for m

  regridding_epoch = 0;
  level_regridding_epochs.resize(1);
  level_regridding_epochs.AT(0) = 0;

  if (verbose or veryverbose) {
    CCTK_INFO("Grid structure (grid points):");
    for (int ml = 0; ml < mglevels; ++ml) {
      int const rl = 0;
      for (int m = 0; m < maps; ++m) {
        for (int c = 0; c < vhh.AT(m)->components(rl); ++c) {
          ibbox const ext = vhh.AT(m)->extent(ml, rl, c);
          ivect const lower = ext.lower();
          ivect const upper = ext.upper();
          int const convfact = ipow(mgfact, ml);
          ibbox const &base = vhh.AT(m)->baseextent(ml, 0);
          ivect const &bstride = base.stride();
          assert(all(lower % bstride == 0));
          assert(all(upper % bstride == 0));
          assert(all(((upper - lower) / bstride) % convfact == 0));
          cout << "   [" << ml << "][" << rl << "][" << m << "][" << c << "]"
               << "   exterior: "
               << "proc " << vhh.AT(m)->processor(rl, c) << "   "
               << lower / bstride << " : " << upper / bstride << "   ("
               << (upper - lower) / bstride / convfact + 1 << ") "
               << prod((upper - lower) / bstride / convfact + 1) << endl;
        }
      }
    }
  }

  // Assert that all maps have one refinement level
  reflevels = 1;
  for (int m = 0; m < maps; ++m) {
    assert(vhh.AT(m)->reflevels() == reflevels);
  }

  timelevels = maxtimelevels;
}

void set_base_extent(int const m, vector<vector<region_t> > &regss) {
  DECLARE_CCTK_PARAMETERS;

  // Create one refinement level
  int const rl = 0;
  regss.resize(1);
  vector<region_t> &regs = regss.AT(rl);

  if (CCTK_EQUALS(base_extents, "")) {

    // Default: one grid component covering everything
    region_t reg;
    reg.extent = vhh.AT(m)->baseextents.AT(0).AT(0);
    reg.outer_boundaries = b2vect(bvect(true));
    reg.map = m;
    regs.push_back(reg);

  } else {

    // Read explicit grid components
    // TODO: invent something for the other convergence levels
    vector<ibbox> exts;
    istringstream ext_str(base_extents);
    try {
      ext_str >> exts;
    } catch (input_error) {
      CCTK_ERROR("Could not parse parameter \"base_extents\"");
    }
    CCTK_VInfo(CCTK_THORNSTRING, "Using %d grid components", int(exts.size()));
    if (exts.size() == 0) {
      CCTK_ERROR("Cannot evolve with zero grid components");
    }

    vector<bbvect> obs;
    istringstream ob_str(base_outerbounds);
    try {
      ob_str >> obs;
    } catch (input_error) {
      CCTK_ERROR("Could not parse parameter \"base_outerbounds\"");
    }
    assert(obs.size() == exts.size());

    for (size_t n = 0; n < exts.size(); ++n) {
      region_t reg;
      reg.extent = exts.AT(n);
      reg.outer_boundaries = xpose(obs.AT(n));
      reg.map = m;
      regs.push_back(reg);
    }
  }
}

void allocate_group_data(cGH const *const cctkGH) {
  groupdata.resize(CCTK_NumGroups());
  arrdata.resize(CCTK_NumGroups());

  for (int group = 0; group < CCTK_NumGroups(); ++group) {

    cGroup gdata;
    check(not CCTK_GroupData(group, &gdata));

    // Check for compact, contiguous, and staggered groups
    ensure_group_options(group, gdata);

    switch (gdata.grouptype) {

    // Grid functions
    case CCTK_GF: {

      // All grid function groups must have the standard rank
      assert(gdata.dim == dim);

      // Set up one refinement level
      groupdata.AT(group).activetimelevels.resize(mglevels);
      for (int ml = 0; ml < mglevels; ++ml) {
        groupdata.AT(group).activetimelevels.AT(ml).resize(1);
      }

      // Initial maximum number of timelevels from interface.ccl
      groupdata.AT(group).info.maxtimelevels = gdata.numtimelevels;

      // Grid function groups use the global grid descriptors
      arrdata.AT(group).resize(maps);
      for (int m = 0; m < maps; ++m) {
        arrdata.AT(group).AT(m).hh = vhh.AT(m);
        arrdata.AT(group).AT(m).dd = vdd.AT(m);
        arrdata.AT(group).AT(m).tt = tt;
      }

      break;
    }

    // Grid scalars and grid arrays are treated in the same way
    case CCTK_SCALAR:
    case CCTK_ARRAY: {

      // All grid variables must have at most the standard rank
      assert(gdata.dim >= 0 and gdata.dim <= dim);

      // Use only one refinement level for grid arrays
      groupdata.AT(group).activetimelevels.resize(mglevels);
      for (int ml = 0; ml < mglevels; ++ml) {
        groupdata.AT(group).activetimelevels.AT(ml).resize(1);
      }

      // Initial maximum number of timelevels from interface.ccl
      groupdata.AT(group).info.maxtimelevels = gdata.numtimelevels;

      // Use only one map for grid arrays
      arrdata.AT(group).resize(1);

      ivect sizes;
      vector<i2vect> ghosts;
      get_group_size(group, gdata, sizes, ghosts);

      // Adapt group sizes for convergence level
      ivect convpowers, convoffsets;
      adapt_group_size_mglevel(group, gdata, sizes, convpowers, convoffsets);
      // Adapt group sizes for disttype
      adapt_group_size_disttype(cctkGH, group, gdata, sizes, ghosts);

      allocate_group_hierarchies(group, sizes, ghosts);

      setup_group_grid_hierarchy(cctkGH, group, gdata, convpowers, convoffsets);

      break;
    } // case scalar or array

    default:
      assert(0);
    } // switch grouptype

    initialise_group_info(cctkGH, group, gdata);

  } // for groups

  output_group_statistics(cctkGH);
}

void allocate_group_hierarchies(int const group, ivect const &sizes,
                                vector<i2vect> const &ghosts) {
  DECLARE_CCTK_PARAMETERS;

  // Calculate base extent
  ivect const str(1);
  ivect const lb(0);
  ivect const ub((sizes - 1) * str);
  ibbox const baseext(lb, ub, str);
  vector<vector<ibbox> > baseexts(1);
  baseexts.AT(0).resize(1);
  baseexts.AT(0).AT(0) = baseext;
  i2vect const nboundaryzones(0);

  // One refinement level
  vector<int> grouptimereffacts(1);
  grouptimereffacts.AT(0) = 1;
  vector<ivect> groupspacereffacts(1);
  groupspacereffacts.AT(0) = ivect(1);

  // There is only one map
  int const m = 0;

  arrdata.AT(group).AT(m).hh =
      new gh(groupspacereffacts, vertex_centered, convergence_factor,
             vertex_centered, baseexts, nboundaryzones);

  vector<i2vect> const buffers(1, i2vect(0));
  vector<i2vect> const overlaps(1, i2vect(0));
  vector<int> const my_prolongation_orders_space(1, 0);
  arrdata.AT(group).AT(m).dd =
      new dh(*arrdata.AT(group).AT(m).hh, ghosts, buffers, overlaps,
             my_prolongation_orders_space);

  arrdata.AT(group).AT(m).tt =
      new th(*arrdata.AT(group).AT(m).hh, maxtimelevels, grouptimereffacts,
             time_interpolation_during_regridding);
}

void setup_group_grid_hierarchy(cGH const *const cctkGH, int const group,
                                cGroup const &gdata, ivect const &convpowers,
                                ivect const &convoffsets) {
#if 0
    // Do not set up anything for groups that have zero variables
    // TODO: To do this, need to change modes.cc to provide illegal
    // entries for GroupDynamicData
    int const nvars = CCTK_NumVarsInGroupI (group);
    assert (nvars >= 0);
    if (nvars == 0) return;
#endif

  // Set refinement structure for scalars and arrays
  vector<region_t> superregs(1);
  int const m = 0;
  int const rl = 0;
  {
    int const c = 0;
    superregs.AT(c).extent =
        arrdata.AT(group).AT(m).hh->baseextents.AT(rl).AT(c);
    superregs.AT(c).outer_boundaries = b2vect(true);
    superregs.AT(c).map = m;
  }
  vector<region_t> regs;

  // Split it into components, one for each processor
  switch (gdata.disttype) {
  case CCTK_DISTRIB_DEFAULT: {
    CCTK_INT no_split_directions[dim];
    bvect no_split_dims(false);
    int const nvals = Util_TableGetIntArray(
        gdata.tagstable, dim, no_split_directions, "no_split_directions");
    assert((0 <= nvals && nvals <= dim) ||
           nvals == UTIL_ERROR_TABLE_NO_SUCH_KEY);
    for (int i = 0; i < nvals; ++i) {
      assert(no_split_directions[i] < dim);
      no_split_dims[no_split_directions[i]] = true;
    }
    SplitRegions_Automatic(cctkGH, superregs, regs, no_split_dims);
    break;
  }
  case CCTK_DISTRIB_CONSTANT: {
    int const d = min(dim - 1, gdata.dim == 0 ? 1 : gdata.dim);
    SplitRegions_AlongDir(cctkGH, superregs, regs, d);
    break;
  }
  default:
    assert(0);
  }

  // Add empty regions if there are fewer regions than processors
  {
    int const nprocs = CCTK_nProcs(cctkGH);
    int const oldsize = regs.size();
    if (oldsize < nprocs) {
      // Ensure that there are at least nprocs components

      // Construct an empty bbox that is just to the right of the
      // last bbox, and which has the correct stride
      assert(not regs.empty());
      ibbox const &lbox = (*regs.rbegin()).extent;
      ibbox const ebox(lbox.upper() + lbox.stride(), lbox.upper(),
                       lbox.stride());

      regs.resize(nprocs);
      for (int c = oldsize; c < nprocs; ++c) {
        region_t empty;
        empty.extent = ebox;
        empty.processor = c;
        regs.AT(c) = empty;
      }
    }
  }

  // Check processor distribution
  {
    int const nprocs = CCTK_nProcs(cctkGH);
    vector<bool> used(nprocs, false);
    for (int c = 0; c < nprocs; ++c) {
      int const p = regs.AT(c).processor;
      assert(p >= 0 and p < nprocs);
      assert(not used.AT(p));
      used.AT(p) = true;
    }
    for (int c = 0; c < nprocs; ++c) {
      assert(used.AT(c));
    }
  }

  // Only one refinement level
  vector<vector<region_t> > superregss(1);
  superregss.AT(rl) = superregs;
  vector<vector<region_t> > regss(1);
  regss.AT(rl) = regs;

  // Create all multigrid levels
  vector<vector<vector<region_t> > > regsss(mglevels);
  ivect mgfact1;
  i2vect offset;
  for (int d = 0; d < dim; ++d) {
    mgfact1[d] = ipow(mgfact, convpowers[d]);
    offset[0][d] = 0;
    offset[1][d] = convoffsets[d];
  }
  regsss.AT(0) = regss;

  for (int ml = 1; ml < mglevels; ++ml) {
    for (int c = 0; c < int(regss.AT(rl).size()); ++c) {
      // this base
      ivect const baselo = ivect(0);
      ivect const baselo1 = baselo;
      // next finer grid
      ivect const flo = regsss.AT(ml - 1).AT(rl).AT(c).extent.lower();
      ivect const fhi = regsss.AT(ml - 1).AT(rl).AT(c).extent.upper();
      ivect const fstr = regsss.AT(ml - 1).AT(rl).AT(c).extent.stride();
      // this grid
      ivect const str = fstr * mgfact1;
      ivect const lo =
          flo + either(regsss.AT(ml - 1).AT(rl).AT(c).outer_boundaries[0],
                       +(offset[0] - mgfact1 * offset[0]) * fstr, ivect(0));
      ivect const hi =
          fhi + either(regsss.AT(ml - 1).AT(rl).AT(c).outer_boundaries[1],
                       -(offset[1] - mgfact1 * offset[1]) * fstr, ivect(0));
      ivect const lo1 = baselo1 + (lo - baselo1 + str - 1) / str * str;
      ivect const hi1 = lo1 + (hi - lo1) / str * str;
      regsss.AT(ml).AT(rl).AT(c) = regsss.AT(ml - 1).AT(rl).AT(c);
      regsss.AT(ml).AT(rl).AT(c).extent = ibbox(lo1, hi1, str);
    }
  }

  // Recompose for this map
  {
    char *const groupname = CCTK_GroupName(group);
    assert(groupname);
    Checkpoint("Recomposing grid array group \"%s\"...", groupname);
    arrdata.AT(group).AT(0).hh->regrid(superregss, regsss, false);
    arrdata.AT(group).AT(0).hh->recompose(0, false);
    arrdata.AT(group).AT(0).hh->regrid_free(false);
    Checkpoint("Done recomposing grid array group \"%s\".", groupname);
    free(groupname);
  }
}

void initialise_group_info(cGH const *const cctkGH, int const group,
                           cGroup const &gdata) {
  DECLARE_CCTK_PARAMETERS;

  // Initialise group information
  groupdata.AT(group).info.dim = gdata.dim;
  groupdata.AT(group).info.gsh = new int[dim];
  groupdata.AT(group).info.lsh = new int[dim];
  groupdata.AT(group).info.ash = new int[dim];
  groupdata.AT(group).info.lbnd = new int[dim];
  groupdata.AT(group).info.ubnd = new int[dim];
  groupdata.AT(group).info.bbox = new int[2 * dim];
  groupdata.AT(group).info.nghostzones = new int[dim];

  groupdata.AT(group).transport_operator =
      get_transport_operator(cctkGH, group, gdata);

  groupdata.AT(group).info.activetimelevels = 0;

  // Initialise group variables
  for (size_t m = 0; m < arrdata.AT(group).size(); ++m) {

    arrdata.AT(group).AT(m).data.resize(CCTK_NumVarsInGroupI(group));
    for (size_t var = 0; var < arrdata.AT(group).AT(m).data.size(); ++var) {
      arrdata.AT(group).AT(m).data.AT(var) = 0;
    }
  }
}

void set_state(cGH *const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // // Allocate level times
  // leveltimes.resize (mglevels);
  // for (int ml=0; ml<mglevels; ++ml) {
  //   leveltimes.AT(ml).resize (1);
  // }

  // Allocate orgin and spacings
  origin_space.resize(maps);
  delta_space.resize(maps);
  for (int m = 0; m < maps; ++m) {
    origin_space.AT(m).resize(mglevels);
  }

  // Current state
  mglevelfact = 1;
  cctkGH->cctk_time = 0.0;
  cctkGH->cctk_delta_time = 1.0;
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_origin_space[d] = 0.0;
    cctkGH->cctk_delta_space[d] = 1.0;
  }

  // Set up things as if in local mode
  mglevel = 0;
  reflevel = 0;
  mc_grouptype = CCTK_GF;
  map = 0;
  component = 0;
  local_component = -1;
  timelevel = 0;
  timelevel_offset = 0;

  // Leave everything, so that everything is set up correctly
  Timers::Timer("meta mode", 1).start();
  Timers::Timer("global mode", 1).start();
  Timers::Timer("level(0)", 1).start();
  if (include_maps_in_mode_timer_tree) {
    Timers::Timer("map(0)", 1).start();
  }
  if (include_local_mode_in_mode_timer_tree) {
    Timers::Timer("local", 1).start();
  }
  leave_local_mode(cctkGH);
  leave_singlemap_mode(cctkGH);
  leave_level_mode(cctkGH);
  leave_global_mode(cctkGH);
}

void enable_storage_for_all_groups(cGH const *const cctkGH) {
  BEGIN_MGLEVEL_LOOP(cctkGH) {
    BEGIN_REFLEVEL_LOOP(cctkGH) {

      for (int group = 0; group < CCTK_NumGroups(); ++group) {
        char *const groupname = CCTK_GroupName(group);
        EnableGroupStorage(cctkGH, groupname);
        free(groupname);
      }
    }
    END_REFLEVEL_LOOP;
  }
  END_MGLEVEL_LOOP;
}

//////////////////////////////////////////////////////////////////////////////
// Routines which do not change state ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

vector<i2vect> get_ghostzones() {
  DECLARE_CCTK_PARAMETERS;

  vector<i2vect> ghostzones;
  // Decide which parameters to use
  if (CCTK_EQUALS(ghost_sizes, "")) {
    i2vect ghostzones1;
    if (ghost_size == -1) {
      ghostzones1 = i2vect(ivect(ghost_size_x, ghost_size_y, ghost_size_z),
                           ivect(ghost_size_x, ghost_size_y, ghost_size_z));
    } else {
      ghostzones1 = i2vect(ivect(ghost_size, ghost_size, ghost_size),
                           ivect(ghost_size, ghost_size, ghost_size));
    }
    ghostzones.resize(maxreflevels, ghostzones1);
  } else {
    // Read them from the parameter
    vector<int> ghostzones1;
    try {
      istringstream gs(ghost_sizes);
      gs >> ghostzones1;
    } catch (input_error) {
      CCTK_WARN(CCTK_WARN_ABORT, "Could not parse parameter \"ghost_sizes\"");
    }
    assert(int(ghostzones1.size()) >= maxreflevels);
    for (int rl = 0; rl < maxreflevels; ++rl) {
      ghostzones.AT(rl) = i2vect(ghostzones1.AT(rl));
    }
  }
  for (int rl = 0; rl < maxreflevels; ++rl) {
    assert(all(all(ghostzones.AT(rl) >= 0)));
  }
  return ghostzones;
}

vector<int> get_prolongation_orders_space() {
  DECLARE_CCTK_PARAMETERS;

  vector<int> orders;
  // Decide which parameters to use
  if (CCTK_EQUALS(prolongation_orders_space, "")) {
    orders.resize(maxreflevels, prolongation_order_space);
  } else {
    // Read them from the parameter
    try {
      istringstream pos(prolongation_orders_space);
      pos >> orders;
    } catch (input_error) {
      CCTK_WARN(CCTK_WARN_ABORT,
                "Could not parse parameter \"prolongation_orders_space\"");
    }
    assert(int(orders.size()) >= maxreflevels);
    for (int rl = 0; rl < maxreflevels; ++rl) {
      assert(orders.AT(rl) >= 0);
    }
  }
  return orders;
}

ivect get_npoints() {
  DECLARE_CCTK_PARAMETERS;

  // Decide which parameters to use
  ivect npoints;
  if (global_nsize == -1) {
    npoints = ivect(global_nx, global_ny, global_nz);
  } else {
    npoints = ivect(global_nsize, global_nsize, global_nsize);
  }

  // Modify npoints for benchmarks
  if (constant_load_per_processor) {
    if (CCTK_EQUALS(processor_topology, "manual")) {
      // Enlarge the domain so that each processor has the specified
      // number of grid points, using the specified processor
      // topology
      assert(processor_topology_3d_x >= 1);
      assert(processor_topology_3d_y >= 1);
      assert(processor_topology_3d_z >= 1);
      ivect const topo = ivect(processor_topology_3d_x, processor_topology_3d_y,
                               processor_topology_3d_z);
      npoints *= topo;
    } else if (CCTK_EQUALS(processor_topology, "automatic")) {
      // Enlarge the domain in a smart way so that each processor
      // has the specified number of grid points
      int const nprocs = dist::total_num_threads();
      // Factorise the number of processors, placing the smallest
      // factors in the tail
      stack<int> factors;
      for (int procsleft = nprocs; procsleft > 1;) {
        for (int divisor = 2; divisor <= procsleft; ++divisor) {
          while (procsleft % divisor == 0) {
            factors.push(divisor);
            procsleft /= divisor;
          }
        }
      }
      // Distribute the factors greedily onto the directions,
      // starting with the largest factor, and preferring to enlarge
      // the x direction
      while (not factors.empty()) {
        int const mindir = minloc(npoints);
        assert(mindir >= 0 and mindir < dim);
        int const factor = factors.top();
        factors.pop();
        npoints[mindir] *= factor;
      }
    } else {
      // TODO: handle processor_topology values "along-z" and "along-dir"
      CCTK_ERROR("Unsupported value of parameter processor_topology");
    }
  } // if constant_load_per_processor
  return npoints;
}

void get_boundary_specification(cGH const *const cctkGH, int const m,
                                vector<i2vect> const &ghosts,
                                i2vect &nboundaryzones, b2vect &is_internal,
                                b2vect &is_staggered, i2vect &shiftout) {
  DECLARE_CCTK_PARAMETERS;

  if (domain_from_multipatch or domain_from_coordbase) {

    jjvect nboundaryzones_, is_internal_, is_staggered_, shiftout_;
    if (CCTK_IsFunctionAliased("MultiPatch_GetBoundarySpecification")) {
      check(not MultiPatch_GetBoundarySpecification(
          m, 2 * dim, &nboundaryzones_[0][0], &is_internal_[0][0],
          &is_staggered_[0][0], &shiftout_[0][0]));
    } else {
      check(not GetBoundarySpecification(
          2 * dim, &nboundaryzones_[0][0], &is_internal_[0][0],
          &is_staggered_[0][0], &shiftout_[0][0]));
    }
    nboundaryzones = xpose(nboundaryzones_);
    is_internal = xpose(is_internal_);
    is_staggered = xpose(is_staggered_);
    shiftout = xpose(shiftout_);

  } else {
    // Legacy code

    // Assume that there are 0 boundary points at outer boundaries
    // and nghostzones boundary points at symmetry boundaries

    // Ensure that all levels have the same number of ghost zones
    for (size_t rl = 1; rl < ghosts.size(); ++rl) {
      assert(all(all(ghosts.AT(rl) == ghosts.AT(0))));
    }

#if 0
      // We cannot call GetSymmetryBoundaries boundaries here, since
      // the symmetry boundaries have not yet been registered.
      jjvect symbnd_;
      check (not GetSymmetryBoundaries (cctkGH, 2*dim, &symbnd_[0][0]));
      b2vect const symbnd = xpose (symbnd_);
#else
    b2vect const symbnd = b2vect(true);
#endif

    for (int f = 0; f < 2; ++f) {
      for (int d = 0; d < dim; ++d) {
        if (symbnd[f][d]) {
          nboundaryzones[f][d] = ghosts.AT(0)[f][d];
          is_internal[f][d] = false;
          // TODO: Look at what CartGrid3D's avoid_origin
          // TODO: Take cell centring into account
          is_staggered[f][d] = false;
          shiftout[f][d] = 1;
        } else {
          // TODO: This is wrong
          nboundaryzones[f][d] = 0;
          is_internal[f][d] = false;
          is_staggered[f][d] = false;
          shiftout[f][d] = 0;
        }
      }
    }
  }

  ostringstream buf;
  buf << "Boundary specification for map " << m << ":" << endl
      << "   nboundaryzones: " << nboundaryzones << endl
      << "   is_internal   : " << is_internal << endl
      << "   is_staggered  : " << is_staggered << endl
      << "   shiftout      : " << shiftout;
  Output(buf.str().c_str());

  if (max_refinement_levels > 1) {
    // Ensure that the boundary is not staggered
    for (int d = 0; d < dim; ++d) {
      for (int f = 0; f < 2; ++f) {
        if (CCTK_EQUALS(refinement_centering, "vertex")) {
          if (is_staggered[f][d]) {
            CCTK_ERROR("The parameters CoordBase::boundary_staggered specify a "
                       "staggered boundary.  Carpet does not support staggered "
                       "boundaries when Carpet::max_refinement_levels > 1 with "
                       "Carpet::centering = \"vertex\"");
          }
        } else if (CCTK_EQUALS(refinement_centering, "cell")) {
          if (not is_staggered[f][d]) {
            CCTK_ERROR("The parameters CoordBase::boundary_staggered specify a "
                       "non-staggered boundary.  Carpet does not support "
                       "non-staggered boundaries when "
                       "Carpet::max_refinement_levels > 1 with "
                       "Carpet::centering = \"cell\"");
          }
        } else {
          assert(0);
        }
      }
    }
  }
}

void get_domain_specification(cGH const *cctkGH, int const m,
                              ivect const &npoints, rvect &physical_min,
                              rvect &physical_max, rvect &base_spacing) {
  DECLARE_CCTK_PARAMETERS;

  rvect interior_min, interior_max;
  rvect exterior_min, exterior_max;

  if (domain_from_multipatch) {
    assert(not domain_from_coordbase);

    // TODO: handle CartGrid3D: either add parameter
    // type=multipatch, and make it handle map numbers, or ignore it
    // altogether, maybe creating a new thorn

    assert(CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification"));
    check(not MultiPatch_GetDomainSpecification(
        m, dim, &physical_min[0], &physical_max[0], &interior_min[0],
        &interior_max[0], &exterior_min[0], &exterior_max[0],
        &base_spacing[0]));

  } else if (domain_from_coordbase) {

    assert(not CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification"));

    // Ensure that CartGrid3D::type = "coordbase"
    ensure_CartGrid3D_type();

    check(not GetDomainSpecification(dim, &physical_min[0], &physical_max[0],
                                     &interior_min[0], &interior_max[0],
                                     &exterior_min[0], &exterior_max[0],
                                     &base_spacing[0]));

  } else {
    // Legacy code

    assert(not CCTK_IsFunctionAliased("MultiPatch_GetDomainSpecification"));

    if (max_refinement_levels > 1) {
      // Ensure that CartGrid3D::avoid_origin = no
      ensure_CartGrid3D_avoid_origin();
    } // if max_refinement_levels > 1

    ostringstream buf;
    buf << "Standard grid specification for map " << m << ":" << endl
        << "   number of grid points: " << npoints;
    Output(buf.str().c_str());

    // Reduce to physical domain
    // TODO: This is not the true domain specification.  However, it
    // is later written to the domainspec, and it is used by Carpet
    // for screen output.
    exterior_min = 0.0;
    exterior_max = rvect(npoints - 1);
    base_spacing = 1.0;
    check(not ConvertFromExteriorBoundary(dim, &physical_min[0],
                                          &physical_max[0], &interior_min[0],
                                          &interior_max[0], &exterior_min[0],
                                          &exterior_max[0], &base_spacing[0]));

  } // if legacy domain specification

  ostringstream buf;
  buf << "CoordBase domain specification for map " << m << ":" << endl
      << "   physical extent: " << physical_min << " : " << physical_max
      << "   (" << physical_max - physical_min << ")" << endl
      << "   interior extent: " << interior_min << " : " << interior_max
      << "   (" << interior_max - interior_min << ")" << endl
      << "   exterior extent: " << exterior_min << " : " << exterior_max
      << "   (" << exterior_max - exterior_min << ")" << endl
      << "   base_spacing   : " << base_spacing;
  Output(buf.str().c_str());
}

void adapt_domain_specification(int const m, rvect const &physical_min,
                                rvect const &physical_max,
                                rvect const &base_spacing, rvect &exterior_min,
                                rvect &exterior_max, rvect &spacing) {
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL const baseconvfact =
      ipow(static_cast<CCTK_REAL>(convergence_factor), basemglevel);
  spacing = base_spacing * baseconvfact;

  rvect interior_min, interior_max;

  if (domain_from_multipatch and
      CCTK_IsFunctionAliased("MultiPatch_ConvertFromPhysicalBoundary")) {
    assert(not domain_from_coordbase);
    check(not MultiPatch_ConvertFromPhysicalBoundary(
        m, dim, &physical_min[0], &physical_max[0], &interior_min[0],
        &interior_max[0], &exterior_min[0], &exterior_max[0], &spacing[0]));
  } else {
    check(not ConvertFromPhysicalBoundary(
        dim, &physical_min[0], &physical_max[0], &interior_min[0],
        &interior_max[0], &exterior_min[0], &exterior_max[0], &spacing[0]));
  }

  ostringstream buf;
  buf << "Adapted domain specification for map " << m << ":" << endl
      << "   convergence factor: " << convergence_factor << endl
      << "   convergence level : " << basemglevel << endl
      << "   physical extent   : " << physical_min << " : " << physical_max
      << "   (" << physical_max - physical_min << ")" << endl
      << "   interior extent   : " << interior_min << " : " << interior_max
      << "   (" << interior_max - interior_min << ")" << endl
      << "   exterior extent   : " << exterior_min << " : " << exterior_max
      << "   (" << exterior_max - exterior_min << ")" << endl
      << "   spacing           : " << spacing;
  Output(buf.str().c_str());
}

void calculate_grid_points(int const m, vector<i2vect> const &ghosts,
                           rvect const &exterior_min, rvect const &exterior_max,
                           rvect const &spacing, ivect &npoints) {
  DECLARE_CCTK_PARAMETERS;

  rvect const real_npoints = either(
      spacing, (exterior_max - exterior_min) / spacing + rvect(1), rvect(1));

  ostringstream buf;
  buf << "Base grid specification for map " << m << ":" << endl
      << "   number of grid points             : " << real_npoints << endl
      << "   number of coarse grid ghost points: " << ghosts.AT(0);
  Output(buf.str().c_str());

  npoints = floor(real_npoints + static_cast<CCTK_REAL>(0.5));

  // Check domain size
  if (any(fabs(rvect(npoints) - real_npoints) >
          static_cast<CCTK_REAL>(0.001))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "The domain size for map %d scaled for convergence level %d "
                "with convergence factor %d is not integer",
                m, int(basemglevel), int(convergence_factor));
  }

  // Sanity check
  assert(all(npoints <= INT_MAX));
#if 0
    int max = INT_MAX;
    for (int d=0; d<dim; ++d) {
      assert (npoints[d] <= max);
      max /= npoints[d];
    }
#endif
  {
    CCTK_REAL const total_npoints = prod(rvect(npoints));
    CCTK_REAL const size_max = numeric_limits<size_type>::max();
    if (total_npoints > size_max) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "The domain for map %d contains %g grid points.  This number "
                  "is larger than the maximum number supported by Carpet (%g).",
                  m, double(total_npoints), double(size_max));
    }
    CCTK_REAL const int_max = numeric_limits<int>::max();
    if (total_npoints > int_max) {
      if (dist::rank() == 0) {
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "The domain for map %d contains %g grid points.  This "
                   "number is larger than the maximum number that can be "
                   "represented as an integer (%g).  This can lead to strange "
                   "problems in thorns which try to calculate the total number "
                   "of grid points.",
                   m, double(total_npoints), double(int_max));
      }
    }
  }

  // Save domain specification
  domainspecs.resize(maps);
  domainspecs.AT(m).exterior_min = exterior_min;
  domainspecs.AT(m).exterior_max = exterior_max;
  domainspecs.AT(m).npoints = npoints;
}

void calculate_base_extents(ibbox const &baseextent,
                            centering const refcentering,
                            centering const mgcentering,
                            i2vect const &nboundaryzones,
                            b2vect const &is_internal,
                            b2vect const &is_staggered, i2vect const &shiftout,
                            vector<vector<ibbox> > &baseextents) {
  DECLARE_CCTK_PARAMETERS;

  assert(baseextents.empty());
  baseextents.resize(num_convergence_levels);
  for (int ml = 0; ml < num_convergence_levels; ++ml) {
    baseextents.AT(ml).resize(maxreflevels);
    for (int rl = 0; rl < maxreflevels; ++rl) {
      if (ml == 0) {
        if (rl == 0) {
          // Use base extent
          baseextents.AT(ml).AT(rl) = baseextent;
        } else {
          // Refine next coarser refinement level
          switch (refcentering) {
          case vertex_centered: {
            ibbox const &cbox = baseextents.AT(ml).AT(rl - 1);
            assert(not any(any(is_staggered)));
            i2vect const bnd_shift =
                (nboundaryzones - 1) * i2vect(not is_internal) + shiftout;
            ibbox const cbox_phys = cbox.expand(-bnd_shift);
            assert(all(baseextent.stride() % spacereffacts.AT(rl - 1) == 0));
            ivect const fstride = baseextent.stride() / spacereffacts.AT(rl);
            ibbox const fbox_phys =
                ibbox(cbox_phys.lower(), cbox_phys.upper(), fstride);
            ibbox const fbox = fbox_phys.expand(bnd_shift);
            baseextents.AT(ml).AT(rl) = fbox;
            break;
          }
          case cell_centered: {
            ibbox const &cbox = baseextents.AT(ml).AT(rl - 1);
            assert(all(all(is_staggered)));
            ivect const cstride = cbox.stride();
            assert(all(cstride % 2 == 0));
            i2vect const bnd_shift_cstride =
                +(cstride / 2 - i2vect(is_internal) * cstride) +
                (nboundaryzones - 1) * i2vect(not is_internal) * cstride +
                shiftout * cstride;
            ibbox const cbox_phys(cbox.lower() - (-bnd_shift_cstride[0]),
                                  cbox.upper() + (-bnd_shift_cstride[1]),
                                  cbox.stride());
            assert(all(baseextent.stride() % spacereffacts.AT(rl - 1) == 0));
            ivect const fstride = baseextent.stride() / spacereffacts.AT(rl);
            ibbox const fbox_phys =
                ibbox(cbox_phys.lower(), cbox_phys.upper(), fstride);
            assert(all(cstride % fstride == 0));
            assert(all(all(bnd_shift_cstride % (cstride / fstride) == 0)));
            i2vect const bnd_shift_fstride =
                bnd_shift_cstride / (cstride / fstride);
            ibbox const fbox(fbox_phys.lower() - bnd_shift_fstride[0],
                             fbox_phys.upper() + bnd_shift_fstride[1],
                             fbox_phys.stride());
            baseextents.AT(ml).AT(rl) = fbox;
            break;
          }
          default:
            assert(0);
          }
        }
      } else {
        // Coarsen next finer convergence level
        assert(mgcentering == vertex_centered);
        assert(not any(any(is_staggered)));
        i2vect const bnd_shift =
            (nboundaryzones - 1) * i2vect(not is_internal) + shiftout;
        ibbox const &fbox = baseextents.AT(ml - 1).AT(rl);
        ibbox const fbox_phys = fbox.expand(-bnd_shift);
        ivect const cstride = fbox.stride() * int(convergence_factor);
        ibbox const cbox_phys =
            ibbox(fbox_phys.lower(), fbox_phys.upper(), cstride);
        ibbox const cbox = cbox_phys.expand(bnd_shift);
        baseextents.AT(ml).AT(rl) = cbox;
      }
    }
  }
}

void find_processor_decomposition(
    cGH const *const cctkGH, vector<vector<vector<region_t> > > &superregsss,
    vector<vector<vector<vector<region_t> > > > &regssss) {
  DECLARE_CCTK_PARAMETERS;

  assert(regssss.empty());
  regssss.resize(maps);

  if (not regrid_in_level_mode) {
    // Distribute each map independently

    for (int m = 0; m < maps; ++m) {
      int const rl = 0;

      vector<vector<region_t> > regss(1);

      // Distribute onto the processors
      SplitRegions(cctkGH, superregsss.AT(m).AT(rl), regss.AT(rl));

      // Create all multigrid levels
      MakeMultigridBoxes(cctkGH, m, regss, regssss.AT(m));
    } // for m

  } else {
    // Distribute all maps at the same time

    int const rl = 0;

    vector<vector<region_t> > superregss(maps);
    for (int m = 0; m < maps; ++m) {
      superregss.AT(m) = superregsss.AT(m).AT(rl);
    }

    vector<vector<region_t> > regss(maps);
    SplitRegionsMaps(cctkGH, superregss, regss);

    vector<vector<vector<region_t> > > regsss(maps);
    for (int m = 0; m < maps; ++m) {
      superregsss.AT(m).AT(rl) = superregss.AT(m);
      regsss.AT(m).resize(1);
      regsss.AT(m).AT(rl) = regss.AT(m);
    }

    // Create all multigrid levels
    MakeMultigridBoxesMaps(cctkGH, regsss, regssss);

  } // if
}

void get_group_size(int const group, cGroup const &gdata, ivect &sizes,
                    vector<i2vect> &ghosts) {
  // Default values
  sizes = 1;
  ghosts.resize(1, i2vect(ivect(0)));

  switch (gdata.grouptype) {

  case CCTK_SCALAR:
    // treat scalars as DIM=0, DISTRIB=const arrays
    assert(gdata.dim == 0);
    assert(gdata.disttype == CCTK_DISTRIB_CONSTANT);
    break;

  case CCTK_ARRAY: {
    assert(gdata.dim >= 1 or gdata.dim <= dim);
    CCTK_INT const *const *const sz = CCTK_GroupSizesI(group);
    CCTK_INT const *const *const gsz = CCTK_GroupGhostsizesI(group);
    // Decode group sizes
    for (int d = 0; d < gdata.dim; ++d) {
      if (sz) {
        sizes[d] = *sz[d];
      }
      if (gsz) {
        ghosts.AT(0)[0][d] = *gsz[d];
        ghosts.AT(0)[1][d] = *gsz[d];
      }
    }
    break;
  }

  default:
    assert(0);
  } // switch grouptype
}

void adapt_group_size_mglevel(int const group, cGroup const &gdata,
                              ivect &sizes, ivect &convpowers,
                              ivect &convoffsets) {
  DECLARE_CCTK_PARAMETERS;

  // Adapt array sizes for convergence level
  get_convergence_options(group, gdata, convpowers, convoffsets);

  ivect baseconvpowers = convpowers * int(basemglevel);
  rvect real_sizes = (((rvect(sizes) - rvect(convoffsets)) /
                       ipow(rvect(convergence_factor), baseconvpowers)) +
                      rvect(convoffsets));
  // Do not modify extra dimensions
  for (int d = gdata.dim; d < dim; ++d) {
    real_sizes[d] = sizes[d];
  }

  // Round group sizes
  sizes = floor(real_sizes + static_cast<CCTK_REAL>(0.5));

  if (any(sizes < 0)) {
    char *const groupname = CCTK_GroupName(group);
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "The shape of group \"%s\" scaled for convergence level %d "
                "with convergence factor %d is negative",
                groupname, int(basemglevel), int(convergence_factor));
    free(groupname);
  }

  if (any(fabs(rvect(sizes) - real_sizes) > static_cast<CCTK_REAL>(1.0e-8))) {
    char *const groupname = CCTK_GroupName(group);
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "The shape of group \"%s\" scaled for convergence level %d "
                "with convergence factor %d is not integer",
                groupname, int(basemglevel), int(convergence_factor));
    free(groupname);
  }
}

void get_convergence_options(int const group, cGroup const &gdata,
                             ivect &convpowers, ivect &convoffsets) {
  if (gdata.tagstable >= 0) {
    {
      jvect convpowers1;
      int const status = Util_TableGetIntArray(
          gdata.tagstable, gdata.dim, &convpowers1[0], "convergence_power");
      if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        // use default: independent of convergence level
        convpowers = 0;
      } else if (status == 1) {
        // a scalar was given
        convpowers = convpowers1[0];
      } else if (status == gdata.dim) {
        convpowers = convpowers1;
      } else {
        char *const groupname = CCTK_GroupName(group);
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "The key \"convergence_power\" in the tags table of group "
                    "\"%s\" is wrong",
                    groupname);
        free(groupname);
      }
      assert(all(convpowers >= 0));
    }

    {
      jvect convoffsets1;
      int const status = Util_TableGetIntArray(
          gdata.tagstable, gdata.dim, &convoffsets1[0], "convergence_offset");
      if (status == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
        // use default: offset is 0
        convoffsets = 0;
      } else if (status == 1) {
        // a scalar was given

      } else if (status == gdata.dim) {
        convoffsets = convoffsets1;
      } else {
        char *const groupname = CCTK_GroupName(group);
        CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                    "The key \"convergence_offset\" in the tags table of group "
                    "\"%s\" is wrong",
                    groupname);
        free(groupname);
      }
    }
  }
}

void adapt_group_size_disttype(cGH const *const cctkGH, int const group,
                               cGroup const &gdata, ivect &sizes,
                               vector<i2vect> const &ghosts) {
  switch (gdata.disttype) {

  case CCTK_DISTRIB_DEFAULT: {
    // do nothing
    break;
  }

  case CCTK_DISTRIB_CONSTANT: {
    if (not all(all(ghosts.AT(0) == 0))) {
      char *const groupname = CCTK_GroupName(group);
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "The group \"%s\" has DISTRIB=constant, but its "
                  "ghostsize is not 0",
                  groupname);
      free(groupname);
    }
    assert(all(all(ghosts.AT(0) == 0)));

    int const nprocs = CCTK_nProcs(cctkGH);

    // Find dimension which should be extended
    int const d = min(dim - 1, gdata.dim == 0 ? 1 : gdata.dim);
    // Extend group sizes
    assert(sizes[d] < numeric_limits<int>::max() / nprocs);
    sizes[d] *= nprocs;
    assert(sizes[d] >= 0);
    break;
  }

  default:
    assert(0);
  } // switch disttype
}

void output_group_statistics(cGH const *const cctkGH) {
  int num_gf_groups = 0;
  int num_gf_vars = 0;
  vector<int> num_array_groups(dim + 1), num_array_vars(dim + 1);
  for (int d = 0; d <= dim; ++d) {
    num_array_groups.AT(d) = 0;
    num_array_vars.AT(d) = 0;
  }

  for (int group = 0; group < CCTK_NumGroups(); ++group) {

    cGroup gdata;
    check(not CCTK_GroupData(group, &gdata));

    switch (gdata.grouptype) {
    case CCTK_GF:
      num_gf_groups += 1;
      num_gf_vars += gdata.numvars * gdata.numtimelevels;
      break;
    case CCTK_SCALAR:
    case CCTK_ARRAY:
      assert(gdata.dim <= dim);
      num_array_groups.AT(gdata.dim) += 1;
      num_array_vars.AT(gdata.dim) += gdata.numvars * gdata.numtimelevels;
      break;
    default:
      assert(0);
    }
  } // for group

  CCTK_INFO("Group and variable statistics:");
  CCTK_VInfo(CCTK_THORNSTRING, "   There are %d grid functions in %d groups",
             num_gf_vars, num_gf_groups);
  CCTK_VInfo(CCTK_THORNSTRING, "   There are %d grid scalars in %d groups",
             num_array_vars.AT(0), num_array_groups.AT(0));
  for (int d = 1; d <= dim; ++d) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "   There are %d %d-dimensional grid arrays in %d groups",
               num_array_vars.AT(d), d, num_array_groups.AT(d));
  }
  CCTK_VInfo(CCTK_THORNSTRING,
             "   (The number of variables counts all time levels)");
}

operator_type get_transport_operator(cGH const *const cctkGH, int const group,
                                     cGroup const &gdata) {
  assert(group >= 0 and group < CCTK_NumGroups());

  if (gdata.grouptype != CCTK_GF) {
    // Ignore everything but true grid functions
    return op_sync; // was: op_copy -- why?
  }

  bool const can_transfer = can_transfer_variable_type(cctkGH, group, gdata);

  // Catch a common error (using the tag "prolongate" instead of
  // "prolongation")
  {
    int const prolong_length =
        Util_TableGetString(gdata.tagstable, 0, NULL, "Prolongate");
    if (prolong_length >= 0) {
      char *const groupname = CCTK_GroupName(group);
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" contains the illegal tag \"Prolongate\".  (Did "
                  "you mean \"Prolongation\" instead?)",
                  groupname);
      free(groupname);
    } else if (prolong_length == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      // good -- do nothing
    } else {
      assert(0);
    }
  }

  // Get prolongation method
  char prolong_string[1000];
  bool have_prolong_string = false;
  {
    int const prolong_length = Util_TableGetString(
        gdata.tagstable, sizeof prolong_string, prolong_string, "Prolongation");
    if (prolong_length >= 0) {
      have_prolong_string = true;
    } else if (prolong_length == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      // do nothing
    } else {
      assert(0);
    }
  }

  // Get prolongation parameter name
  char prolong_param_string[1000];
  bool have_prolong_param_string = false;
  {
    int const prolong_param_length =
        Util_TableGetString(gdata.tagstable, sizeof prolong_param_string,
                            prolong_param_string, "ProlongationParameter");
    if (prolong_param_length >= 0) {
      have_prolong_param_string = true;
    } else if (prolong_param_length == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
      // do nothing
    } else {
      assert(0);
    }
  }

  // Complain if both are given
  if (have_prolong_string and have_prolong_param_string) {
    char *const groupname = CCTK_GroupName(group);
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Group \"%s\" has both the tags \"Prolongation\" and "
                "\"ProlongationParameter\".  This is not possible.",
                groupname);
    free(groupname);
  }

  // Map the parameter name
  if (have_prolong_param_string) {
    char *thorn;
    char *name;
    int const ierr = CCTK_DecomposeName(prolong_param_string, &thorn, &name);
    if (ierr < 0) {
      char *const groupname = CCTK_GroupName(group);
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" has the \"ProlongationParameter\" tag \"%s\".  "
                  "This is not a valid parameter name.",
                  groupname, prolong_param_string);
      free(groupname);
    }
    int type;
    char const *const *const value = (static_cast<char const *const *>(
        CCTK_ParameterGet(name, thorn, &type)));
    if (not value or not*value) {
      char *const groupname = CCTK_GroupName(group);
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" has the \"ProlongationParameter\" tag \"%s\".  "
                  "This parameter does not exist.",
                  groupname, prolong_param_string);
      free(groupname);
    }
    if (type != PARAMETER_KEYWORD and type != PARAMETER_STRING) {
      char *const groupname = CCTK_GroupName(group);
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "Group \"%s\" has the \"ProlongationParameter\" tag \"%s\".  "
                  "This parameter has the wrong type; it must be either "
                  "KEYWORD or STRING.",
                  groupname, prolong_param_string);
      free(groupname);
    }
    free(thorn);
    free(name);
    assert(strlen(*value) < sizeof prolong_string);
    strcpy(prolong_string, *value);
    have_prolong_string = true;
  }

  // Select a default, if necessary
  if (not have_prolong_string) {
    if (can_transfer) {
// Use the default
#if 0
        if (gdata.numtimelevels == 1) {
          // Only one time level:
          char * const groupname = CCTK_GroupName (group);
          CCTK_VWarn (2, __LINE__, __FILE__, CCTK_THORNSTRING,
                      "Group \"%s\" has only one time level; therefore it will not be prolongated or restricted.",
                      groupname);
          free (groupname);
          return op_sync;
        } else {
          // Several time levels: use the default
          return op_Lagrange;
        }
#endif
      return op_Lagrange;
    } else {
      if (gdata.grouptype == CCTK_GF) {
        char *const groupname = CCTK_GroupName(group);
        CCTK_VWarn(2, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Group \"%s\" has the variable type \"%s\" which cannot be "
                   "prolongated or restricted.",
                   groupname, CCTK_VarTypeName(gdata.vartype));
        free(groupname);
        return op_sync;
      } else {
        return op_error;
      }
    }
  }

  // Select the prolongation method
  assert(have_prolong_string);
  if (CCTK_Equals(prolong_string, "none")) {
    // This would surprise too many people
    // return op_none;
    return op_sync;
  } else if (CCTK_Equals(prolong_string, "sync")) {
    return op_sync;
  } else if (CCTK_Equals(prolong_string, "restrict")) {
    return op_restrict;
  } else if (CCTK_Equals(prolong_string, "copy")) {
    return op_copy;
  } else if (CCTK_Equals(prolong_string, "Lagrange")) {
    return op_Lagrange;
  } else if (CCTK_Equals(prolong_string, "ENO")) {
    return op_ENO;
  } else if (CCTK_Equals(prolong_string, "WENO")) {
    return op_WENO;
  } else if (CCTK_Equals(prolong_string, "TVD")) {
    return op_TVD;
  } else if (CCTK_Equals(prolong_string, "Lagrange_monotone")) {
    return op_Lagrange_monotone;
  } else if (CCTK_Equals(prolong_string, "STAGGER011")) {
    return op_STAGGER011;
  } else if (CCTK_Equals(prolong_string, "STAGGER101")) {
    return op_STAGGER101;
  } else if (CCTK_Equals(prolong_string, "STAGGER110")) {
    return op_STAGGER110;
  } else if (CCTK_Equals(prolong_string, "STAGGER111")) {
    return op_STAGGER111;
  } else {
    char *const groupname = CCTK_GroupName(group);
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "Group \"%s\" has the unknown prolongation method \"%s\".",
                groupname, prolong_string);
    free(groupname);
    return op_error;
  }
  return op_error;
}

bool can_transfer_variable_type(cGH const *const cctkGH, int const group,
                                cGroup const &gdata) {
// Find out which types correspond to the default types
#if CCTK_INTEGER_PRECISION_1
#define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT1
#elif CCTK_INTEGER_PRECISION_2
#define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT2
#elif CCTK_INTEGER_PRECISION_4
#define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT4
#elif CCTK_INTEGER_PRECISION_8
#define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT8
#elif CCTK_INTEGER_PRECISION_16
#define CCTK_DEFAULT_INTEGER_TYPE CCTK_VARIABLE_INT16
#else
#error "Unsupported default integer type"
#endif

#if CCTK_REAL_PRECISION_4
#define CCTK_DEFAULT_REAL_TYPE CCTK_VARIABLE_REAL4
#define CCTK_DEFAULT_COMPLEX_TYPE CCTK_VARIABLE_COMPLEX8
#elif CCTK_REAL_PRECISION_8
#define CCTK_DEFAULT_REAL_TYPE CCTK_VARIABLE_REAL8
#define CCTK_DEFAULT_COMPLEX_TYPE CCTK_VARIABLE_COMPLEX16
#elif CCTK_REAL_PRECISION_16
#define CCTK_DEFAULT_REAL_TYPE CCTK_VARIABLE_REAL16
#define CCTK_DEFAULT_COMPLEX_TYPE CCTK_VARIABLE_COMPLEX32
#else
#error "Unsupported default real type"
#endif

  int const type0 = gdata.vartype;
  int type1;

  switch (type0) {
  case CCTK_VARIABLE_INT:
    type1 = CCTK_DEFAULT_INTEGER_TYPE;
    break;
  case CCTK_VARIABLE_REAL:
    type1 = CCTK_DEFAULT_REAL_TYPE;
    break;
  case CCTK_VARIABLE_COMPLEX:
    type1 = CCTK_DEFAULT_COMPLEX_TYPE;
    break;
  default:
    type1 = type0;
  }
  switch (type1) {

#ifdef HAVE_CCTK_REAL8
  case CCTK_VARIABLE_REAL8:
    // This type is supported.
    return true;
#endif

#ifdef HAVE_CCTK_REAL4
  case CCTK_VARIABLE_REAL4:
#endif
#ifdef HAVE_CCTK_REAL16
  case CCTK_VARIABLE_REAL16:
#endif
#ifdef HAVE_CCTK_COMPLEX8
  case CCTK_VARIABLE_COMPLEX8:
#endif
#ifdef HAVE_CCTK_COMPLEX16
  case CCTK_VARIABLE_COMPLEX16:
#endif
#ifdef HAVE_CCTK_COMPLEX32
  case CCTK_VARIABLE_COMPLEX32:
#endif
    // This type is not supported, but could be.
    return false;

  case CCTK_VARIABLE_BYTE:
#ifdef HAVE_CCTK_INT1
  case CCTK_VARIABLE_INT1:
#endif
#ifdef HAVE_CCTK_INT2
  case CCTK_VARIABLE_INT2:
#endif
#ifdef HAVE_CCTK_INT4
  case CCTK_VARIABLE_INT4:
#endif
#ifdef HAVE_CCTK_INT8
  case CCTK_VARIABLE_INT8:
#endif
#ifdef HAVE_CCTK_INT16
  case CCTK_VARIABLE_INT16:
#endif
    // This type is not supported, and cannot be.
    return false;

  default: {
    CCTK_VError(
        __LINE__, __FILE__, CCTK_THORNSTRING,
        "Internal error: encountered variable type %d (%s) for group %d (%s)",
        type1, CCTK_VarTypeName(type1), group, CCTK_GroupName(group));
  }
  }

  // not reached
  abort();
  return false;
}

//////////////////////////////////////////////////////////////////////////////
// Parameter and consistency checking ////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

// Ensure that CartGrid3D::type = "coordbase"
void ensure_CartGrid3D_type() {
  if (CCTK_IsThornActive("CartGrid3D")) {
    int type;
    void const *const ptr = CCTK_ParameterGet("type", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_KEYWORD);
    char const *const coordtype = *static_cast<char const *const *>(ptr);
    if (not CCTK_EQUALS(coordtype, "coordbase")) {
      CCTK_ERROR("When Carpet::domain_from_coordbase = yes, and when thorn "
                 "CartGrid3D is active, then you also have to set "
                 "CartGrid3D::type = \"coordbase\"");
    }
  }
}

#if 0
  // UNUSED:
  // Ensure that CartGrid3D doesn't apply symmetries
  void
  ensure_CartGrid3D_domain ()
  {
    if (CCTK_IsThornActive ("CartGrid3D")) {
      int type;
      void const * ptr;
      
      ptr = CCTK_ParameterGet ("domain", "CartGrid3D", & type);
      assert (ptr != 0);
      assert (type == PARAMETER_KEYWORD);
      char const * const domain
        = * static_cast<char const * const *> (ptr);
      if (not CCTK_EQUALS (domain, "full")) {
        CCTK_ERROR ("When Carpet::domain_from_coordbase = no, and when Carpet::max_refinement_levels > 1, then thorn CartGrid3D cannot provide symmetry boundaries");
      }
    }
  }
#endif

// Ensure that CartGrid3D::avoid_origin = no
void ensure_CartGrid3D_avoid_origin() {
  if (CCTK_IsThornActive("CartGrid3D")) {
    int type;
    void const *ptr;

    ptr = CCTK_ParameterGet("no_origin", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const no_origin = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("no_originx", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const no_originx = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("no_originy", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const no_originy = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("no_originz", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const no_originz = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("avoid_origin", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const avoid_origin = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("avoid_originx", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const avoid_originx = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("avoid_originy", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const avoid_originy = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("avoid_originz", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const avoid_originz = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("domain", "CartGrid3D", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_KEYWORD);
    char const *const domain = *static_cast<char const *const *>(ptr);

    bvect stag;
    stag[0] = no_origin and no_originx and avoid_origin and avoid_originx;
    stag[1] = no_origin and no_originy and avoid_origin and avoid_originy;
    stag[2] = no_origin and no_originz and avoid_origin and avoid_originz;

    // TODO: Check only if there is actually a symmetry boundary
    // TODO: Take cell centring into account
    if (not CCTK_EQUALS(domain, "full") and any(stag)) {
      CCTK_ERROR("When Carpet::domain_from_coordbase = no, when "
                 "Carpet::max_refinement_levels > 1, and when thorn CartGrid3D "
                 "provides symmetry boundaries, then you have to set "
                 "CartGrid3D::avoid_origin = no");
    }
  }
}

// Ensure that ReflectionSymmetry::avoid_origin = no (if vertex
// centred)
void ensure_ReflectionSymmetry_avoid_origin(centering const refcentering) {
  if (CCTK_IsThornActive("ReflectionSymmetry")) {
    int type;
    void const *ptr;

    ptr = CCTK_ParameterGet("avoid_origin_x", "ReflectionSymmetry", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const avoid_origin_x = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("avoid_origin_y", "ReflectionSymmetry", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const avoid_origin_y = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("avoid_origin_z", "ReflectionSymmetry", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const avoid_origin_z = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("reflection_x", "ReflectionSymmetry", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const reflection_x = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("reflection_y", "ReflectionSymmetry", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const reflection_y = *static_cast<CCTK_INT const *>(ptr);

    ptr = CCTK_ParameterGet("reflection_z", "ReflectionSymmetry", &type);
    assert(ptr != 0);
    assert(type == PARAMETER_BOOLEAN);
    CCTK_INT const reflection_z = *static_cast<CCTK_INT const *>(ptr);

    if (refcentering == vertex_centered) {
      if ((reflection_x and avoid_origin_x) or
          (reflection_y and avoid_origin_y) or
          (reflection_z and avoid_origin_z)) {
        CCTK_ERROR("When Carpet::max_refinement_levels > 1, and when "
                   "ReflectionSymmetry::symmetry_[xyz] = yes, then you also "
                   "have to set ReflectionSymmetry::avoid_origin_[xyz] = no");
      }
    } else if (refcentering == cell_centered) {
      if ((reflection_x and not avoid_origin_x) or
          (reflection_y and not avoid_origin_y) or
          (reflection_z and not avoid_origin_z)) {
        CCTK_ERROR("When Carpet::max_refinement_levels > 1, and when "
                   "ReflectionSymmetry::symmetry_[xyz] = yes, then you have to "
                   "set ReflectionSymmetry::avoid_origin_[xyz] = yes");
      }
    } else {
      assert(0);
    }
  }
}

void ensure_ghostzones(int const m, vector<i2vect> const &ghosts) {
  DECLARE_CCTK_PARAMETERS;

  int const rls = vhh.AT(m)->reffacts.size();
  for (int rl = 0; rl < rls; ++rl) {
    int const my_prolongation_order_space =
        vdd.AT(m)->prolongation_orders_space.AT(rl);
    int const prolongation_stencil_size =
        vdd.AT(m)->prolongation_stencil_size(rl);
    int const min_nghosts =
        ((prolongation_stencil_size + refinement_factor - 1) /
         (refinement_factor - 1));
    int const min_nghosts_restrict = restriction_order_space / 2;
    if (any(any(ghosts.AT(rl) < i2vect(min_nghosts)))) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "There are not enough ghost zones for the desired spatial "
                  "prolongation order on map %d, refinement level %d.  With a "
                  "spatial prolongation order of %d, you need at least %d "
                  "ghost zones.",
                  m, rl, my_prolongation_order_space, min_nghosts);
    }
    if (use_higher_order_restriction and
        any(any(ghosts.AT(rl) < i2vect(min_nghosts_restrict)))) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                  "There are not enough ghost zones for the desired "
                  "restriction order on map %d, refinement level %d.  With a "
                  "restriction order of %d, you need at least %d ghost zones.",
                  m, rl, int(restriction_order_space), min_nghosts_restrict);
    }
  }
}

void ensure_group_options(int const group, cGroup const &gdata) {
#ifdef CCTK_HAVE_COMPACT_GROUPS
  if (gdata.compact) {
    char *const groupname = CCTK_GroupName(group);
    CCTK_VError(
        __LINE__, __FILE__, CCTK_THORNSTRING,
        "The group \"%s\" has COMPACT=1.  Compact groups are not yet supported",
        groupname);
    free(groupname);
  }
#endif

#ifdef CCTK_HAVE_CONTIGUOUS_GROUPS
  if (gdata.contiguous) {
    char *const groupname = CCTK_GroupName(group);
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
                "The group \"%s\" has CONTIGUOUS=1.  Contiguous groups are not "
                "yet supported",
                groupname);
    free(groupname);
  }
#endif
}

} // namespace Carpet
