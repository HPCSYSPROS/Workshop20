#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstring>
#include <list>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#ifdef HAVE_CAPABILITY_MPI
#include <mpi.h>

#if defined __bgq__
// The processor names on a Blue Gene/Q include the MPI rank, and thus
// do not uniquely identify the host name. We therefore roll our own.
// See <https://wiki.alcf.anl.gov/parts/index.php/Blue_Gene/Q>.
#include <mpix.h>
namespace {
void MPI_Get_processor_name1(char *name, int *resultlen) {
  MPIX_Hardware_t hw;
  MPIX_Hardware(&hw);
  *resultlen =
      snprintf(name, MPI_MAX_PROCESSOR_NAME, "(%u,%u,%u,%u,%u)", hw.Coords[0],
               hw.Coords[1], hw.Coords[2], hw.Coords[3], hw.Coords[4]);
  // ignoring hw.Coords[5], which is the core number inside a node
}
}
#define MPI_Get_processor_name MPI_Get_processor_name1
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#else
#include <sys/time.h>
int omp_get_max_threads() { return 1; }
int omp_get_num_threads() { return 1; }
int omp_get_thread_num() { return 0; }
#endif

#include <hwloc.h>
// On a Blue Gene/Q, hwloc reports per-process hardware information
// instead of per-node hardware information. We need to correct for
// this.
#ifdef __bgq__
#define HWLOC_PER_PROCESS
#endif

using namespace std;

namespace {
int divup(int a, int b) {
  assert(a >= 0);
  assert(b > 0);
  return (a + b - 1) / b;
}

bool is_pow2(int a) {
  if (a <= 0)
    return false;
  while (a != 1) {
    if (a % 2)
      return false;
    a /= 2;
  }
  return true;
}
}

namespace {

// Check that OpenMP counts and numbers threads as expected
void check_openmp() {
  bool found_inconsistency = false;

  // Count OpenMP threads
  int num_threads_direct = 0;
#pragma omp parallel reduction(+ : num_threads_direct)
  { ++num_threads_direct; }
  int num_threads_omp = -1;
#pragma omp parallel
  {
#pragma omp master
    { num_threads_omp = omp_get_num_threads(); }
  }
  if (num_threads_direct != num_threads_omp) {
    found_inconsistency = true;
    CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Number of OpenMP threads is inconsistent: counting %d threads, "
               "but OpenMP run-time reports %d threads",
               num_threads_direct, num_threads_omp);
  }

  // Check OpenMP thread numbers
  vector<int> thread_nums;
#pragma omp parallel
  {
#pragma omp critical
    { thread_nums.push_back(omp_get_thread_num()); }
  }
  // Prevent insanity
  assert(int(thread_nums.size()) == num_threads_direct);
  int max_thread_num = -1;
  for (int i = 0; i < int(thread_nums.size()); ++i) {
    max_thread_num = max(max_thread_num, thread_nums.at(i));
  }
  vector<int> thread_counts(max_thread_num + 1, 0);
  for (int i = 0; i < int(thread_nums.size()); ++i) {
    ++thread_counts.at(thread_nums.at(i));
  }
  int num_threads_direct_again = 0;
  for (size_t i = 0; i < thread_counts.size(); ++i) {
    num_threads_direct_again += thread_counts.at(i);
  }
  // Prevent insanity
  assert(num_threads_direct_again == num_threads_direct);
  bool thread_counts_bad = int(thread_counts.size()) < num_threads_direct;
  for (int i = 0; i < int(thread_counts.size()); ++i) {
    thread_counts_bad =
        thread_counts_bad or thread_counts.at(i) != (i < num_threads_direct);
  }
  if (thread_counts_bad) {
    found_inconsistency = true;
    printf("OpenMP thread numbers:");
    for (int i = 0; i < int(thread_counts.size()); ++i) {
      for (int j = 0; j < thread_counts.at(i); ++j) {
        printf(" %d", i);
      }
    }
    printf("\n");
    CCTK_WARN(CCTK_WARN_ALERT, "OpenMP threads are numbered inconsistently");
  }

  if (found_inconsistency) {
    CCTK_ERROR("Severe OpenMP inconsistency detected -- aborting");
  }
}
}

namespace {

struct mpi_host_mapping_t {
  int mpi_num_procs, mpi_proc_num;
  int mpi_num_hosts, mpi_host_num;
  int mpi_num_procs_on_host, mpi_proc_num_on_host;

  void load();
};

mpi_host_mapping_t *mpi_host_mapping = NULL;

#ifdef HAVE_CAPABILITY_MPI

void mpi_host_mapping_t::load() {
  CCTK_INFO("MPI process-to-host mapping:");

  MPI_Comm comm = MPI_COMM_WORLD;
  int const root = 0;

  MPI_Comm_size(comm, &mpi_num_procs);
  MPI_Comm_rank(comm, &mpi_proc_num);
  printf("This is MPI process %d of %d\n", mpi_proc_num, mpi_num_procs);

  char procname[MPI_MAX_PROCESSOR_NAME];
  int procnamelen;
  MPI_Get_processor_name(procname, &procnamelen);
  vector<char> procnames;
  if (mpi_proc_num == root) {
    procnames.resize(MPI_MAX_PROCESSOR_NAME * mpi_num_procs);
  }
  MPI_Gather(procname, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, &procnames[0],
             MPI_MAX_PROCESSOR_NAME, MPI_CHAR, root, comm);
  vector<int> host_byproc;
  map<int, list<int> > host2procs1;
  if (mpi_proc_num == root) {
    map<string, int> hostname2host;
    vector<string> hostnames;
    hostnames.reserve(mpi_num_procs);
    host_byproc.resize(mpi_num_procs);
    mpi_num_hosts = 0;
    for (int proc = 0; proc < mpi_num_procs; ++proc) {
      string hostname(&procnames[MPI_MAX_PROCESSOR_NAME * proc]);
      if (hostname2host.count(hostname) == 0) {
        hostname2host[hostname] = mpi_num_hosts;
        hostnames.push_back(hostname);
        ++mpi_num_hosts;
      }
      int host = hostname2host[hostname];
      host_byproc[proc] = host;
      host2procs1[host].push_back(proc);
    }
    printf("MPI hosts:\n");
    for (int host = 0; host < mpi_num_hosts; ++host) {
      printf("  %d: %s\n", host, hostnames[host].c_str());
    }
  }
  MPI_Bcast(&mpi_num_hosts, 1, MPI_INT, root, comm);
  MPI_Scatter(&host_byproc[0], 1, MPI_INT, &mpi_host_num, 1, MPI_INT, root,
              comm);
  printf("This MPI process runs on host %d of %d\n", mpi_host_num,
         mpi_num_hosts);

  vector<int> num_procs_on_host_byproc;
  vector<int> proc_num_on_host_byproc;
  if (mpi_proc_num == root) {
    vector<vector<int> > host2procs(mpi_num_hosts);
    vector<int> num_procs_on_host_byhost(mpi_num_hosts);
    num_procs_on_host_byproc.resize(mpi_num_procs);
    proc_num_on_host_byproc.resize(mpi_num_procs);
    for (int host = 0; host < mpi_num_hosts; ++host) {
      list<int> const &host_procs1 = host2procs1[host];
      vector<int> &host_procs = host2procs[host];
      host_procs.reserve(host_procs1.size());
      for (list<int>::const_iterator iproc = host_procs1.begin();
           iproc != host_procs1.end(); ++iproc) {
        host_procs.push_back(*iproc);
      }
      sort(host_procs.begin(), host_procs.end());
      int num_procs_on_host = host_procs.size();
      num_procs_on_host_byhost[host] = num_procs_on_host;
      for (int proc_num_on_host = 0; proc_num_on_host < num_procs_on_host;
           ++proc_num_on_host) {
        int proc = host_procs[proc_num_on_host];
        num_procs_on_host_byproc[proc] = num_procs_on_host;
        proc_num_on_host_byproc[proc] = proc_num_on_host;
      }
    }
  }
  MPI_Scatter(&num_procs_on_host_byproc[0], 1, MPI_INT, &mpi_num_procs_on_host,
              1, MPI_INT, root, comm);
  MPI_Scatter(&proc_num_on_host_byproc[0], 1, MPI_INT, &mpi_proc_num_on_host, 1,
              MPI_INT, root, comm);
  printf("On this host, this is MPI process %d of %d\n", mpi_proc_num_on_host,
         mpi_num_procs_on_host);
}

#else

void mpi_host_mapping_t::load() {
  mpi_num_procs = 1;
  mpi_proc_num = 0;
  mpi_num_hosts = 1;
  mpi_host_num = 0;
  mpi_num_procs_on_host = 1;
  mpi_proc_num_on_host = 0;
}

#endif
}

// Inspired by code in hwloc's documentation

namespace {

#define OUTPUT_SUPPORT(FIELD)                                                  \
  printf("  %-41s: %s\n", #FIELD, topology_support->FIELD ? "yes" : "no")

void output_support(hwloc_topology_t topology) {
  CCTK_INFO("Topology support:");
  hwloc_topology_support const *topology_support =
      hwloc_topology_get_support(topology);
  printf("Discovery support:\n");
  OUTPUT_SUPPORT(discovery->pu);
  printf("CPU binding support:\n");
  OUTPUT_SUPPORT(cpubind->set_thisproc_cpubind);
  OUTPUT_SUPPORT(cpubind->get_thisproc_cpubind);
  OUTPUT_SUPPORT(cpubind->set_proc_cpubind);
  OUTPUT_SUPPORT(cpubind->get_proc_cpubind);
  OUTPUT_SUPPORT(cpubind->set_thisthread_cpubind);
  OUTPUT_SUPPORT(cpubind->get_thisthread_cpubind);
  OUTPUT_SUPPORT(cpubind->set_thread_cpubind);
  OUTPUT_SUPPORT(cpubind->get_thread_cpubind);
  OUTPUT_SUPPORT(cpubind->get_thisproc_last_cpu_location);
  OUTPUT_SUPPORT(cpubind->get_proc_last_cpu_location);
  OUTPUT_SUPPORT(cpubind->get_thisthread_last_cpu_location);
  printf("Memory binding support:\n");
  OUTPUT_SUPPORT(membind->set_thisproc_membind);
  OUTPUT_SUPPORT(membind->get_thisproc_membind);
  OUTPUT_SUPPORT(membind->set_proc_membind);
  OUTPUT_SUPPORT(membind->get_proc_membind);
  OUTPUT_SUPPORT(membind->set_thisthread_membind);
  OUTPUT_SUPPORT(membind->get_thisthread_membind);
  OUTPUT_SUPPORT(membind->set_area_membind);
  OUTPUT_SUPPORT(membind->get_area_membind);
  OUTPUT_SUPPORT(membind->alloc_membind);
  OUTPUT_SUPPORT(membind->firsttouch_membind);
  OUTPUT_SUPPORT(membind->bind_membind);
  OUTPUT_SUPPORT(membind->interleave_membind);
  OUTPUT_SUPPORT(membind->replicate_membind);
  OUTPUT_SUPPORT(membind->nexttouch_membind);
  OUTPUT_SUPPORT(membind->migrate_membind);
}

void output_object(hwloc_topology_t topology, hwloc_obj_t obj, int depth) {
  char type_buf[1000], attr_buf[1000];
  hwloc_obj_type_snprintf(type_buf, sizeof type_buf, obj, 1);
  hwloc_obj_attr_snprintf(attr_buf, sizeof attr_buf, obj, ", ", 1);
  printf("%*s%s L#%d: (P#%d%s%s)\n", 2 * depth, "", type_buf,
         obj->logical_index, obj->os_index, strlen(attr_buf) == 0 ? "" : ", ",
         attr_buf);

  for (unsigned i = 0; i < obj->arity; ++i) {
    // TODO: output index (physical, logical?)
    output_object(topology, obj->children[i], depth + 1);
  }
}

void output_objects(hwloc_topology_t topology) {
  CCTK_INFO("Hardware objects in this node:");
  output_object(topology, hwloc_get_root_obj(topology), 0);
}

void output_bindings(hwloc_topology_t topology,
                     mpi_host_mapping_t const &host_mapping) {
  hwloc_topology_support const *topology_support =
      hwloc_topology_get_support(topology);
  if (not topology_support->cpubind->get_thisthread_cpubind) {
    CCTK_INFO("Cannot determine thread CPU bindings");
    return;
  }

  CCTK_INFO("Thread CPU bindings:");
  // Output all information about host 0
  int const root = 0;
  int const host = 0;
  if (host_mapping.mpi_proc_num == root) {
    if (host_mapping.mpi_host_num != host or
        host_mapping.mpi_proc_num_on_host != 0) {
      CCTK_ERROR("Unexpected host numbering -- root process is not process 0 "
                 "on host 0");
    }
  }
  if (host_mapping.mpi_host_num == host) {
    ostringstream buf;
    buf << "  "
        << "MPI process " << host_mapping.mpi_proc_num << " "
        << "on host " << host_mapping.mpi_host_num << " "
        << "(process " << host_mapping.mpi_proc_num_on_host << " "
        << "of " << host_mapping.mpi_num_procs_on_host << " "
        << "on this host)\n";

    int const pu_depth = hwloc_get_type_or_below_depth(topology, HWLOC_OBJ_PU);
    assert(pu_depth >= 0);
    int const num_pus = hwloc_get_nbobjs_by_depth(topology, pu_depth);
    assert(num_pus > 0);
#pragma omp parallel
    {
      int const num_threads = omp_get_num_threads();
      for (int thread = 0; thread < num_threads; ++thread) {
        if (thread == omp_get_thread_num()) {
          hwloc_cpuset_t cpuset = hwloc_bitmap_alloc();
          if (cpuset) {
            int const ierr =
                hwloc_get_cpubind(topology, cpuset, HWLOC_CPUBIND_THREAD);
            if (not ierr) {
              hwloc_cpuset_t lcpuset = hwloc_bitmap_alloc();
              if (lcpuset) {
                for (int pu_num = 0; pu_num < num_pus; ++pu_num) {
                  hwloc_obj_t const pu_obj =
                      hwloc_get_obj_by_depth(topology, pu_depth, pu_num);
                  if (hwloc_bitmap_isset(cpuset, pu_obj->os_index)) {
                    hwloc_bitmap_set(lcpuset, pu_num);
                  }
                }
                char lcpuset_buf[1000];
                hwloc_bitmap_list_snprintf(lcpuset_buf, sizeof lcpuset_buf,
                                           lcpuset);
                hwloc_bitmap_free(lcpuset);
                char cpuset_buf[1000];
                hwloc_bitmap_list_snprintf(cpuset_buf, sizeof cpuset_buf,
                                           cpuset);
                // printf("OpenMP thread %d: PU set L#{%s} P#{%s}\n",
                //        thread, lcpuset_buf, cpuset_buf);
                buf << "    "
                    << "OpenMP thread " << thread << ": "
                    << "PU set L#{" << lcpuset_buf << "} "
                    << "P#{" << cpuset_buf << "}\n";
              } else {
                CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__,
                           CCTK_THORNSTRING,
                           "Could not allocate bitmap for CPU bindings");
              }
            } else {
              CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                         "Could not obtain CPU binding for thread %d", thread);
            }
            hwloc_bitmap_free(cpuset);
          } else {
            CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                       "Could not allocate bitmap for CPU bindings");
          }
        }
#pragma omp barrier
      }
    }

    // Collect all output
    string const bufstr = buf.str();
    // Output for root process
    printf("%s", bufstr.c_str());
#ifdef HAVE_CAPABILITY_MPI
    // Output for other processes
    if (host_mapping.mpi_proc_num == root) {
      // Receive
      for (int proc = 1; proc < host_mapping.mpi_num_procs_on_host; ++proc) {
        int rbuflen;
        MPI_Recv(&rbuflen, 1, MPI_INT, MPI_ANY_SOURCE, proc, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        vector<char> rbufstr(rbuflen + 1);
        MPI_Recv(&rbufstr[0], rbuflen, MPI_CHAR, MPI_ANY_SOURCE, proc,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        rbufstr[rbuflen] = '\0';
        printf("%s", &rbufstr[0]);
      }
    } else {
      // Send
      int const buflen = bufstr.size();
      MPI_Send(const_cast<int *>(&buflen), 1, MPI_INT, root,
               host_mapping.mpi_proc_num_on_host, MPI_COMM_WORLD);
      MPI_Send(const_cast<char *>(bufstr.c_str()), buflen, MPI_CHAR, root,
               host_mapping.mpi_proc_num_on_host, MPI_COMM_WORLD);
    }
#endif

  } // if on root host
}

void set_bindings(hwloc_topology_t topology,
                  mpi_host_mapping_t const &host_mapping) {
  DECLARE_CCTK_PARAMETERS;

  hwloc_topology_support const *topology_support =
      hwloc_topology_get_support(topology);
  if (not topology_support->cpubind->set_thisthread_cpubind) {
    CCTK_INFO("Cannot set thread CPU bindings");
    return;
  }

  bool dense_layout;
  if (CCTK_EQUALS(thread_layout, "dense")) {
    dense_layout = true;
  } else if (CCTK_EQUALS(thread_layout, "loose")) {
    dense_layout = false;
  } else {
    assert(0);
  }

  // TODO: set memory binding policy as well

  // TODO: use hwloc_distribute instead
  CCTK_INFO("Setting thread CPU bindings:");
#pragma omp parallel
  {
    // All quantities are per host
    int const core_depth =
        hwloc_get_type_or_below_depth(topology, HWLOC_OBJ_CORE);
    assert(core_depth >= 0);
    int const num_cores = hwloc_get_nbobjs_by_depth(topology, core_depth);
    assert(num_cores > 0);
    int const pu_depth = hwloc_get_type_or_below_depth(topology, HWLOC_OBJ_PU);
    assert(pu_depth >= 0);
    int const num_pus = hwloc_get_nbobjs_by_depth(topology, pu_depth);
    assert(num_pus > 0);
    assert(num_pus % num_cores == 0);
    int const max_smt_threads = num_pus / num_cores;
    int const num_threads_in_proc = omp_get_num_threads();
    int const num_procs = host_mapping.mpi_num_procs_on_host;
    int const num_threads = num_threads_in_proc * num_procs;
    int const num_smt_threads = divup(num_threads, num_cores);
    int const core_spacing =
        dense_layout ? 1 : num_cores / (num_threads / num_smt_threads);
    int const proc_num = host_mapping.mpi_proc_num_on_host;
    int const thread_offset = num_threads_in_proc * proc_num;
    // Bind thread to exactly one PU
    for (int thread_num_in_proc = 0; thread_num_in_proc < num_threads_in_proc;
         ++thread_num_in_proc) {
      if (thread_num_in_proc == omp_get_thread_num()) {
        int const thread_num = thread_offset + thread_num_in_proc;
        // Map requested threads to existing cores, oversubscribing
        // or spacing if necessary
        int const core_num = thread_num / num_smt_threads * core_spacing;
        assert(core_num < num_cores);
        // Map requested SMT threads to existing PUs,
        // oversubscribing if necessary
        int const pu_offset =
            thread_num % num_smt_threads * max_smt_threads / num_smt_threads;
        assert(pu_offset < max_smt_threads);
        int const pu_num = core_num * max_smt_threads + pu_offset;
        assert(pu_num < num_pus);
        hwloc_obj_t pu_obj = hwloc_get_obj_by_depth(topology, pu_depth, pu_num);
        assert(pu_obj);
        // hwloc_cpuset_t cpuset = hwloc_bitmap_dup(pu_obj->cpuset);
        hwloc_cpuset_t cpuset = pu_obj->cpuset;
        int ierr;
        ierr = hwloc_set_cpubind(topology, cpuset,
                                 HWLOC_CPUBIND_THREAD | HWLOC_CPUBIND_STRICT);
        if (ierr) {
          ierr = hwloc_set_cpubind(topology, cpuset, HWLOC_CPUBIND_THREAD);
        }
        if (ierr) {
          CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                     "Could not set CPU binding for thread %d",
                     thread_num_in_proc);
        }
        // hwloc_bitmap_free(cpuset);
      }
#pragma omp barrier
    }
  }
}
}

namespace {

enum memory_t { mem_cache = 0, mem_local = 1, mem_global = 2 };
struct node_topology_info_t {
  int num_smt_threads; // threads per core
  int max_smt_threads; // pus per core
  struct cache_info_t {
    char const *name;
    memory_t type;
    ptrdiff_t size; // data cache size in bytes (0 if unknown)
    int linesize;   // data cache line size in bytes (0 if unknown)
    int stride;     // data cache stride in bytes (0 if unknown)
    int num_pus;    // number of PUs which share this cache
  };
  vector<cache_info_t> cache_info;

  void load(hwloc_topology_t const &topology,
            mpi_host_mapping_t const &mpi_host_mapping);
};

node_topology_info_t *node_topology_info = NULL;

void node_topology_info_t::load(hwloc_topology_t const &topology,
                                mpi_host_mapping_t const &host_mapping) {
  CCTK_INFO("Extracting CPU/cache/memory properties:");
  // All quantities are per host
  int const core_depth =
      hwloc_get_type_or_below_depth(topology, HWLOC_OBJ_CORE);
  assert(core_depth >= 0);
  int num_cores = hwloc_get_nbobjs_by_depth(topology, core_depth);
#ifdef HWLOC_PER_PROCESS
  num_cores *= host_mapping.mpi_num_procs_on_host;
#endif
  assert(num_cores > 0);
  int const pu_depth = hwloc_get_type_or_below_depth(topology, HWLOC_OBJ_PU);
  assert(pu_depth >= 0);
  int num_pus = hwloc_get_nbobjs_by_depth(topology, pu_depth);
#ifdef HWLOC_PER_PROCESS
  num_pus *= host_mapping.mpi_num_procs_on_host;
#endif
  assert(num_pus > 0);
  assert(num_pus % num_cores == 0);
  max_smt_threads = num_pus / num_cores;
  printf("  There are %d PUs per core (aka hardware SMT threads)\n",
         max_smt_threads);
  int const num_threads_in_proc = omp_get_max_threads();
  int const num_procs = host_mapping.mpi_num_procs_on_host;
  int const num_threads = num_threads_in_proc * num_procs;
  // TODO: calculate this instead by looking at logical core numbers
  // for each thread
  num_smt_threads = divup(num_threads, num_cores);
  printf("  There are %d threads per core (aka SMT threads used)\n",
         num_smt_threads);
  if (num_smt_threads > max_smt_threads) {
    printf("WARNING: This is larger than the number of hardware SMT threads\n");
  }
  if (num_threads_in_proc % num_smt_threads != 0) {
    printf("WARNING: This does not evenly divide the number of threads per "
           "process\n");
  }
  assert(num_smt_threads > 0);

  assert(cache_info.empty());
  for (int cache_level = 1; true; ++cache_level) {
    int const cache_depth =
        hwloc_get_cache_type_depth(topology, cache_level, HWLOC_OBJ_CACHE_DATA);
    if (cache_depth < 0)
      break;
    int const num_caches = hwloc_get_nbobjs_by_depth(topology, cache_depth);
    assert(num_caches > 0);
    int const cache_num = 0; // just look at first cache
    hwloc_obj_t const cache_obj =
        hwloc_get_obj_by_depth(topology, cache_depth, cache_num);
    assert(cache_obj->type == HWLOC_OBJ_CACHE);
    hwloc_obj_attr_u::hwloc_cache_attr_s const &cache_attr =
        cache_obj->attr->cache;
    char const *const cache_type_str =
        cache_attr.type == HWLOC_OBJ_CACHE_UNIFIED
            ? "unified"
            : cache_attr.type == HWLOC_OBJ_CACHE_DATA
                  ? "data"
                  : cache_attr.type == HWLOC_OBJ_CACHE_INSTRUCTION
                        ? "instruction"
                        : "UNKNOWN";
    ostringstream namebuf;
    switch (cache_attr.type) {
    case HWLOC_OBJ_CACHE_UNIFIED:
      namebuf << "L";
      break;
    case HWLOC_OBJ_CACHE_DATA:
      namebuf << "D";
      break;
    case HWLOC_OBJ_CACHE_INSTRUCTION:
      namebuf << "I";
      break;
    default:
      namebuf << "?";
      break;
    }
    namebuf << cache_attr.depth << " cache";
    string const name = namebuf.str();
    int const cache_stride = (cache_attr.associativity == 0
                                  ? 0
                                  : cache_attr.size / cache_attr.associativity);
    printf("  Cache %s has type \"%s\" depth %u\n"
           "    size %td linesize %u associativity %d stride %d, "
           "for %d PUs\n",
           cache_obj->name ? cache_obj->name : "(unknown name)", cache_type_str,
           cache_attr.depth, (ptrdiff_t)cache_attr.size, cache_attr.linesize,
           cache_attr.associativity, cache_stride, num_pus / num_caches);
    if (cache_attr.type != HWLOC_OBJ_CACHE_INSTRUCTION) {
      // assert(cache_attr.linesize >= 0); // linesize is unsigned
      assert(cache_attr.linesize == 0 or is_pow2(cache_attr.linesize));
      assert(cache_stride >= 0);
      // Cache strides may not be powers of two
      // assert(cache_stride==0 or is_pow2(cache_stride));
      cache_info_t new_cache_info;
      new_cache_info.name = strdup(name.c_str());
      new_cache_info.type = mem_cache;
      new_cache_info.size = cache_attr.size;
      new_cache_info.linesize = cache_attr.linesize;
      new_cache_info.stride = cache_stride;
      new_cache_info.num_pus = num_pus / num_caches; // TODO
      cache_info.push_back(new_cache_info);
    }
  }

  // Describe (socket-local) memory as well, creating fake cache
  // entries
  for (int type_idx = 0; type_idx < 2; ++type_idx) {
    hwloc_obj_type_t const obj_type =
        type_idx == 0 ? HWLOC_OBJ_NODE : HWLOC_OBJ_MACHINE;
    int const node_depth = hwloc_get_type_depth(topology, obj_type);
    if (node_depth >= 0) {
      int const num_nodes = hwloc_get_nbobjs_by_depth(topology, node_depth);
      assert(num_nodes > 0);
      int const num_memory_levels = num_nodes == 1 ? 1 : 2;
      int const node_num = 0; // just look at first node
      hwloc_obj_t const node_obj =
          hwloc_get_obj_by_depth(topology, node_depth, node_num);
      assert(node_obj->type == obj_type);
      hwloc_obj_memory_s const &memory_attr = node_obj->memory;
      for (int memory_level = 0; memory_level < num_memory_levels;
           ++memory_level) {
        int const num_memories = memory_level == 0 ? 1 : num_nodes;
        char const *const name =
            memory_level == 0 ? "local memory" : "global memory";
        ptrdiff_t const memory_size = memory_attr.local_memory;
        if (memory_size > 0) {
          ptrdiff_t page_size;
          if (memory_attr.page_types_len > 0) {
            // Use smallest page size
            page_size = memory_attr.page_types[0].size;
          } else {
            page_size = 0;
          }
          printf("  Memory has type \"%s\" depth %d\n"
                 "    size %td pagesize %td, for %d PUs\n",
                 memory_level == 0 ? "local" : "global", node_depth,
                 memory_size * num_memories, page_size,
                 num_pus * num_memories / num_nodes);
          cache_info_t new_cache_info;
          new_cache_info.name = name;
          new_cache_info.type = memory_level == 0 ? mem_local : mem_global;
          new_cache_info.size = memory_size * num_memories;
          new_cache_info.linesize = page_size;
          new_cache_info.stride = 0;
          new_cache_info.num_pus = num_pus * num_memories / num_nodes; // TODO
          cache_info.push_back(new_cache_info);
        }
      }
    }
  }
}
}

extern "C" CCTK_INT ST_GetNumSMTThreads() {
  return node_topology_info->num_smt_threads;
}

extern "C" CCTK_INT ST_GetMaxSMTThreads() {
  return node_topology_info->max_smt_threads;
}

extern "C" CCTK_INT ST_GetCacheInfo(CCTK_POINTER_TO_CONST *restrict const names,
                                    CCTK_INT *restrict const types,
                                    CCTK_POINTER_TO_CONST *restrict const sizes,
                                    CCTK_INT *restrict const linesizes,
                                    CCTK_INT *restrict const strides,
                                    CCTK_INT *restrict const num_puss,
                                    CCTK_INT const max_num_cache_levels) {
  vector<node_topology_info_t::cache_info_t> const &cache_info =
      node_topology_info->cache_info;

  int const num_levels = min(int(max_num_cache_levels), int(cache_info.size()));
  for (int level = 0; level < num_levels; ++level) {
    if (names)
      names[level] = cache_info[level].name;
    if (types)
      types[level] = cache_info[level].type;
    if (sizes)
      sizes[level] = (void *)(cache_info[level].size);
    if (linesizes)
      linesizes[level] = cache_info[level].linesize;
    if (strides)
      strides[level] = cache_info[level].stride;
    if (num_puss)
      num_puss[level] = cache_info[level].num_pus;
  }
  return cache_info.size();
}

extern "C" CCTK_INT
ST_GetMPIProcessInfo(CCTK_INT *restrict const mpi_num_procs,
                     CCTK_INT *restrict const mpi_proc_num,
                     CCTK_INT *restrict const mpi_num_hosts,
                     CCTK_INT *restrict const mpi_host_num,
                     CCTK_INT *restrict const mpi_num_procs_on_host,
                     CCTK_INT *restrict const mpi_proc_num_on_host) {
  *mpi_num_procs = mpi_host_mapping->mpi_num_procs;
  *mpi_proc_num = mpi_host_mapping->mpi_proc_num;
  *mpi_num_hosts = mpi_host_mapping->mpi_num_hosts;
  *mpi_host_num = mpi_host_mapping->mpi_host_num;
  *mpi_num_procs_on_host = mpi_host_mapping->mpi_num_procs_on_host;
  *mpi_proc_num_on_host = mpi_host_mapping->mpi_proc_num_on_host;
  return 0;
}

extern "C" int ST_system_topology() {
  DECLARE_CCTK_PARAMETERS;

  // Check OpenMP consistency
  check_openmp();

  // Determine MPI (host/process) mapping
  mpi_host_mapping = new mpi_host_mapping_t;
  mpi_host_mapping->load();

  // Determine node topology
  hwloc_topology_t topology;
  hwloc_topology_init(&topology);
  hwloc_topology_load(topology);

  output_support(topology);
  output_objects(topology);
  output_bindings(topology, *mpi_host_mapping);
  // TODO: output distance matrix

  bool do_set_thread_bindings;
  if (CCTK_EQUALS(set_thread_bindings, "yes")) {
    do_set_thread_bindings = true;
  } else if (CCTK_EQUALS(set_thread_bindings, "no")) {
    do_set_thread_bindings = false;
  } else if (CCTK_EQUALS(set_thread_bindings, "auto")) {
// TODO: may want to handle some systems specially here
#ifdef HWLOC_PER_PROCESS
    // The calculation is probably wrong, so don't act on it
    do_set_thread_bindings = false;
#else
    do_set_thread_bindings = true;
#endif
  } else {
    CCTK_BUILTIN_UNREACHABLE();
  }
  if (do_set_thread_bindings) {
    set_bindings(topology, *mpi_host_mapping);
    output_bindings(topology, *mpi_host_mapping);
  }

  // Capture some information for later use
  node_topology_info = new node_topology_info_t;
  node_topology_info->load(topology, *mpi_host_mapping);

  hwloc_topology_destroy(topology);

  return 0;
}
