#include <cctk.h>

#include <cassert>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#ifdef CCTK_MPI
#include <mpi.h>
#else
#include "nompi.h"
#endif

#include <dist.hh>
#include <functions.hh>
#include <mpi_string.hh>

namespace Carpet {

using namespace std;
using namespace CarpetLib;

vector<string> host_names;       // Host id to host name
std::map<string, int> host_map;  // Host name to host id
vector<int> host_ids;            // Process to host id
vector<vector<int> > host_procs; // Host id to processes

vector<string> const &HostNames() { return host_names; }
std::map<string, int> const &HostMap() { return host_map; }
vector<int> const &HostIds() { return host_ids; }
vector<vector<int> > const &HostProcs() { return host_procs; }

string HostName(int const id) { return host_names.AT(id); }

int HostId(string const name) {
  if (host_map.find(name) != host_map.end()) {
    return host_map[name];
  } else {
    return -1;
  }
}

int HostId(int const proc) { return host_ids.AT(proc); }

vector<int> const &HostProcs(int const id) { return host_procs.AT(id); }

void DetermineHosts(string const host, bool const verbose) {
  // Get the total number of processes
  int num_procs;
  MPI_Comm_size(dist::comm(), &num_procs);
  int my_proc;
  MPI_Comm_rank(dist::comm(), &my_proc);

  // Gather all host names
  vector<string> const hosts(allgather_string(dist::comm(), host));

  // Map host strings to small integers
  int num_hosts = 0;
  host_ids.resize(num_procs);
  host_map.clear();
  for (int n = 0; n < num_procs; ++n) {
    if (host_map.find(hosts.AT(n)) != host_map.end()) {
      host_ids.AT(n) = host_map[hosts.AT(n)];
    } else {
      host_map[hosts.AT(n)] = num_hosts;
      host_ids.AT(n) = num_hosts;
      ++num_hosts;
    }
  }

  // Determine processes per host
  vector<int> num_host_procs(num_hosts, 0);
  for (int n = 0; n < num_procs; ++n) {
    ++num_host_procs.AT(host_ids.AT(n));
  }

  host_names.resize(num_hosts);
  host_procs.resize(num_hosts);
  for (int m = 0; m < num_hosts; ++m) {
    host_procs.AT(m).reserve(num_host_procs.AT(m));
  }
  for (int n = 0; n < num_procs; ++n) {
    host_names.AT(host_ids.AT(n)) = hosts.AT(n);
    host_procs.AT(host_ids.AT(n)).push_back(n);
  }
  for (int m = 0; m < num_hosts; ++m) {
    assert(static_cast<int>(host_procs.AT(m).size()) == num_host_procs.AT(m));
  }

  if (verbose) {
    CCTK_INFO("Host listing:");
    for (int m = 0; m < num_hosts; ++m) {
      cout << "   host " << m << ": \"" << host_names.AT(m) << "\"" << endl;
    }
    CCTK_INFO("Host/process mapping:");
    for (int n = 0; n < num_procs; ++n) {
      int const m = host_ids.AT(n);
      bool const same_host_as_prev = n - 1 >= 0 and host_ids.AT(n - 1) == m;
      bool const same_host_as_next =
          n + 1 < num_procs and host_ids.AT(n + 1) == m;
      if (same_host_as_next) {
        if (same_host_as_prev) {
          // Output nothing
        } else {
          // This process has the same host as the next one:
          // output only a partial line
          cout << "   processes " << n << "-";
        }
      } else {
        if (same_host_as_prev) {
          // This process has the same host as the previous one:
          // finish a partial line
          cout << n << ": "
               << "host " << m << " \"" << host_names.AT(m) << "\"" << endl;
        } else {
          cout << "   process " << n << ": "
               << "host " << m << " \"" << host_names.AT(m) << "\"" << endl;
        }
      }
    }
    int const my_host = host_ids.AT(my_proc);
    CCTK_VInfo(CCTK_THORNSTRING,
               "Host mapping: This is process %d, host %d \"%s\"", my_proc,
               my_host, host.c_str());
  }
}

} // namespace Carpet
