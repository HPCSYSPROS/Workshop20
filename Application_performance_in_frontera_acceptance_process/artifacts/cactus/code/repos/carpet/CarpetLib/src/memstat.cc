#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <sys/resource.h>
#include <sys/time.h>

#include "mem.hh"

#include "dh.hh"
#include "gdata.hh"
#include "ggf.hh"
#include "gh.hh"
#include "th.hh"

struct mstat {
  // Carpet object statistics
  double total_bytes;
  double total_objects;
  double max_bytes;
  double max_objects;
  // Carpet administrative data structure statistics
  double total_admin_bytes;
  // malloc statistics
  double malloc_used_bytes;
  double malloc_free_bytes;
};
int const mstat_entries = sizeof(mstat) / sizeof(double);

extern "C" void CarpetLib_printmemstats(CCTK_ARGUMENTS);

void CarpetLib_printmemstats(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const ioproc = 0;

  static int next_output = 0;
  if ((print_memstats_every == 0 and cctk_iteration == 0) or
      (print_memstats_every > 0 and cctk_iteration >= next_output)) {
    next_output = cctk_iteration + print_memstats_every;

    mstat mybuf;
    mybuf.total_bytes = gmem::total_allocated_bytes;
    mybuf.total_objects = gmem::total_allocated_objects;
    mybuf.max_bytes = gmem::max_allocated_bytes;
    mybuf.max_objects = gmem::max_allocated_objects;
    mybuf.total_admin_bytes = gh::allmemory() + dh::allmemory() +
                              th::allmemory() + ggf::allmemory() +
                              gdata::allmemory();
#ifdef HAVE_MALLINFO
    // NOTE: struct mallinfo returns byte-counts as int, which can
    // overflow.  In this case, the information is incorrect.
    struct mallinfo const minfo = mallinfo();
    mybuf.malloc_used_bytes = minfo.uordblks;
    mybuf.malloc_free_bytes = minfo.fordblks;
#else
    mybuf.malloc_used_bytes = 0;
    mybuf.malloc_free_bytes = 0;
#endif

    cout << "Memory statistics from CarpetLib:" << eol
         << "   Current number of objects: "
         << size_t(gmem::total_allocated_objects) << eol
         << "   Current allocated memory:  " << setprecision(3)
         << gmem::total_allocated_bytes / gmem::MEGA << " MB" << eol
         << "   Maximum number of objects: "
         << size_t(gmem::max_allocated_objects) << eol
         << "   Maximum allocated memory:  " << setprecision(3)
         << gmem::max_allocated_bytes / gmem::MEGA << " MB" << eol
         << "   Current administrative memory: " << setprecision(3)
         << mybuf.total_admin_bytes / gmem::MEGA << " MB" << eol
         << "   Total allocated used system memory: " << setprecision(3)
         << mybuf.malloc_used_bytes / gmem::MEGA << " MB" << eol
         << "   Total allocated free system memory: " << setprecision(3)
         << mybuf.malloc_free_bytes / gmem::MEGA << " MB" << endl;

    // TODO: improve this message
    cout << "   gh::allmemory:    " << gh::allmemory() << eol
         << "   dh::allmemory:    " << dh::allmemory() << eol
         << "   th::allmemory:    " << th::allmemory() << eol
         << "   ggf::allmemory:   " << ggf::allmemory() << eol
         << "   gdata::allmemory: " << gdata::allmemory() << endl;

    if (strcmp(memstat_file, "") != 0) {
      vector<mstat> allbuf(dist::size());
      MPI_Gather(&mybuf, mstat_entries, MPI_DOUBLE, &allbuf.front(),
                 mstat_entries, MPI_DOUBLE, ioproc, dist::comm());

      if (dist::rank() == ioproc) {

        double max_total_bytes = 0;
        double avg_total_bytes = 0;
        double cnt_total_bytes = 0;
        double max_max_bytes = 0;
        double avg_max_bytes = 0;
        double cnt_max_bytes = 0;
        double max_admin_bytes = 0;
        double avg_admin_bytes = 0;
        double cnt_admin_bytes = 0;
        double max_used_bytes = 0;
        double avg_used_bytes = 0;
        double cnt_used_bytes = 0;
        double max_free_bytes = 0;
        double avg_free_bytes = 0;
        double cnt_free_bytes = 0;
        for (size_t n = 0; n < allbuf.size(); ++n) {
          max_total_bytes = max(max_total_bytes, allbuf[n].total_bytes);
          avg_total_bytes += allbuf[n].total_bytes;
          ++cnt_total_bytes;
          max_max_bytes = max(max_max_bytes, allbuf[n].max_bytes);
          avg_max_bytes += allbuf[n].max_bytes;
          ++cnt_max_bytes;
          max_admin_bytes = max(max_admin_bytes, allbuf[n].total_admin_bytes);
          avg_admin_bytes += allbuf[n].total_admin_bytes;
          ++cnt_admin_bytes;
          max_used_bytes = max(max_used_bytes, allbuf[n].malloc_used_bytes);
          avg_used_bytes += allbuf[n].malloc_used_bytes;
          ++cnt_used_bytes;
          max_free_bytes = max(max_free_bytes, allbuf[n].malloc_free_bytes);
          avg_free_bytes += allbuf[n].malloc_free_bytes;
          ++cnt_free_bytes;
        }
        avg_total_bytes /= cnt_total_bytes;
        avg_max_bytes /= cnt_max_bytes;
        avg_admin_bytes /= cnt_admin_bytes;
        avg_used_bytes /= cnt_used_bytes;
        avg_free_bytes /= cnt_free_bytes;

        ostringstream filenamebuf;
        filenamebuf << out_dir << "/" << memstat_file;
        string const filename = filenamebuf.str();
        ofstream file;
        static bool did_truncate = false;
        if (not did_truncate) {
          did_truncate = true;
          file.open(filename.c_str(), ios::out | ios::trunc);
          if (CCTK_IsFunctionAliased("UniqueBuildID")) {
            char const *const build_id =
                static_cast<char const *>(UniqueBuildID(cctkGH));
            file << "# Build ID: " << build_id << eol;
          }
          if (CCTK_IsFunctionAliased("UniqueSimulationID")) {
            char const *const job_id =
                static_cast<char const *>(UniqueSimulationID(cctkGH));
            file << "# Simulation ID: " << job_id << eol;
          }
          file << "# Running on " << dist::size() << " processes" << eol;
          file << "#" << eol;
          file << "# iteration   maxtotalbytes avgtotalbytes   maxmaxbytes "
                  "avgm avgfreebytes"
               << eol;
        } else {
          file.open(filename.c_str(), ios::out | ios::app);
        }

        file << cctk_iteration << "\t " << max_total_bytes << " "
             << avg_total_bytes << "\t " << max_max_bytes << " "
             << avg_max_bytes << "\t " << max_admin_bytes << " "
             << avg_admin_bytes << "\t " << max_used_bytes << " "
             << avg_used_bytes << "\t " << max_free_bytes << " "
             << avg_free_bytes << eol;

        file.close();

      } // if on root process
    }   // if output to file
  }
}
