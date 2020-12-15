/*@@
  @file      Output.c
  @date      July 6 2003
  @author    Gabrielle Allen
  @desc
             Functions to report the timers
  @enddesc
@@*/

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Schedule.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

#ifdef CCTK_MPI
#include <mpi.h>
#endif

namespace TimerReport {

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

struct timer_stats {
  int ntimers;
  vector<string> names;       // timer names
  vector<CCTK_REAL> secs_avg; // global average
  vector<CCTK_REAL> secs_min; // global minimum
  vector<CCTK_REAL> secs_max; // global maximum
  timer_stats() : ntimers(0) {}
};

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

void Output(CCTK_ARGUMENTS);
void PrintTimes(CCTK_ARGUMENTS);
void OutputAllTimers(CCTK_ARGUMENTS);
void OutputAllTimersTogether(CCTK_ARGUMENTS);
void OutputAllTimersReadable(CCTK_ARGUMENTS);
void PrintTopTimers(CCTK_ARGUMENTS);

int CollectTimerInfo(const cGH *restrict cctkGH, timer_stats &timers);

string QuoteForCSV(const string &);
string QuoteForTSV(const string &);

/********************************************************************
 *********************  Scheduled Routine Prototypes  ***************
 ********************************************************************/

extern "C" void TimerReport_OutputEvery(CCTK_ARGUMENTS);
extern "C" void TimerReport_OutputTerminate(CCTK_ARGUMENTS);
extern "C" void TimerReport_Checkpoint(CCTK_ARGUMENTS);

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

const string sep(70, '=');

/********************************************************************
 ********************    External Routines   ************************
 ********************************************************************/

/*@@
  @routine    TimerReport_OutputEvery
  @date       2008-11-12
  @author     Erik Schnetter
  @desc
  Output the timer table periodically
  @enddesc
@@*/
void TimerReport_OutputEvery(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (next || out_at == cctk_iteration ||
      (out_every && cctk_iteration % out_every == 0)) {
    CCTK_VInfo(CCTK_THORNSTRING, "Timer Report at iteration %d time %g",
               cctk_iteration, double(cctk_time));
    Output(CCTK_PASS_CTOC);

    if (next) {
      CCTK_ParameterSet("next", CCTK_THORNSTRING, "no");
    }
  }
}

/*@@
  @routine    TimerReport_OutputTerminate
  @date       2008-11-12
  @author     Erik Schnetter
  @desc
  Output the timer table periodically
  @enddesc
@@*/
void TimerReport_OutputTerminate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VInfo(CCTK_THORNSTRING,
             "Timer Report before terminating at iteration %d time %g",
             cctk_iteration, double(cctk_time));
  Output(CCTK_PASS_CTOC);
}

/*@@
  @routine    TimerReport_Checkpoint
  @date       April 10 2004
  @author     Erik Schnetter
  @desc
  Output the timer table if before_checkpoint is set
  @enddesc
@@*/
void TimerReport_Checkpoint(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (before_checkpoint &&
      (checkpoint_every && cctk_iteration % checkpoint_every == 0)) {
    CCTK_VInfo(CCTK_THORNSTRING,
               "Timer Report before checkpointing at iteration %d, time %g",
               cctk_iteration, double(cctk_time));
    Output(CCTK_PASS_CTOC);
  }
}

/********************************************************************
 ********************    Internal Routines   ************************
 ********************************************************************/

void Output(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (output_schedule_timers) {
    PrintTimes(CCTK_PASS_CTOC);
  }

  if (output_all_timers) {
    OutputAllTimers(CCTK_PASS_CTOC);
  }

  if (output_all_timers_together) {
    OutputAllTimersTogether(CCTK_PASS_CTOC);
  }

  if (output_all_timers_readable) {
    OutputAllTimersReadable(CCTK_PASS_CTOC);
  }

  if (n_top_timers > 0) {
    PrintTopTimers(CCTK_PASS_CTOC);
  }
}

void PrintTimes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(out_filename, "")) {
    // Print to stdout.
    CCTK_SchedulePrintTimes(NULL);
  } else {
    // Print to a file.
    ostringstream filename_buf;
    int myproc = CCTK_MyProc(cctkGH);
    filename_buf << out_dir << "/" << out_filename << "." << setw(6)
                 << setfill('0') << myproc << ".txt";
    string filename = filename_buf.str();
    // truncate or append
    static bool first_time = true;
    string flags = first_time && IO_TruncateOutputFiles(cctkGH) ? "w" : "a";
    first_time = false;

    FILE *file = fopen(filename.c_str(), flags.c_str());
    if (file) {
      // Print the schedule to the file
      fprintf(file, "Timer Report at iteration %d time %g:\n\n", cctk_iteration,
              double(cctk_time));
      CCTK_SchedulePrintTimesToFile(NULL, file);
      fprintf(file, "\n********************************************************"
                    "************************\n");
      fclose(file);
    } else {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Could not open timer report output file \"%s\"",
                 filename.c_str());
    }
  }
}

void OutputAllTimers(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  static bool first_time = true;
  static int last_ntimers = -1;

  int ntimers = CCTK_NumTimers();
  assert(ntimers >= 0);
  cTimerData *td = CCTK_TimerCreateData();

  ostringstream filename_buf;
  filename_buf << out_dir << "/"
               << "AllTimers"
               << "." << setw(6) << setfill('0') << CCTK_MyProc(cctkGH)
               << ".txt";
  string filename = filename_buf.str();

  // truncate or append
  string flags = first_time && IO_TruncateOutputFiles(cctkGH) ? "w" : "a";

  FILE *file = fopen(filename.c_str(), flags.c_str());
  if (file) {
    if (first_time) {
      fprintf(file, "# Clock %s\n", all_timers_clock);
      fprintf(file, "# Unit seconds\n");
    }
    // If the number of timers has changed, output the header containing the
    // timer names again. It would be better to track creation and deletion
    // of timers in the flesh, but this always works as long as timers are
    // never deleted.
    if (first_time || last_ntimers != ntimers) {
      fprintf(file, "# Column 1 iteration\n");
      fprintf(file, "# Column 2 simulation time\n");
      for (int i = 0; i < ntimers; i++) {
        const char *name = CCTK_TimerName(i);

        if (name == nullptr)
          name = "";

        fprintf(file, "# Column %d %s\n", i + 3, name);
      }
    }

    fprintf(file, "%d\t%.15g", cctk_iteration, double(cctk_time));

    for (int i = 0; i < ntimers; i++) {
      CCTK_TimerI(i, td);

      const cTimerVal *tv = CCTK_GetClockValue(all_timers_clock, td);

      double timer_secs;
      if (tv != nullptr) {
        timer_secs = CCTK_TimerClockSeconds(tv);
      } else {
        const char *name = CCTK_TimerName(i);
        if (name == nullptr)
          name = "(null)";
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Clock \"%s\" not found for timer #%d \"%s\"",
                   all_timers_clock, i, name);
        timer_secs = -1;
      }

      fprintf(file, "\t%.15g", timer_secs);
    }
    fprintf(file, "\n");
    fclose(file);
  } else {
    CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not open timer report output file \"%s\"",
               filename.c_str());
  }
  CCTK_TimerDestroyData(td);
  first_time = false;
  last_ntimers = ntimers;
}

void OutputAllTimersTogether(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  timer_stats timers;
  if (!CollectTimerInfo(cctkGH, timers))
    return;

  if (CCTK_MyProc(cctkGH) == 0) {
    static bool first_time = true;
    static int last_ntimers = -1;

    ostringstream filename_buf, filename_csv_buf, filename_tsv_buf;
    filename_buf << out_dir << "/"
                 << "AllTimers"
                 << ".txt";
    filename_csv_buf << out_dir << "/"
                     << "AllTimers"
                     << ".csv";
    filename_tsv_buf << out_dir << "/"
                     << "AllTimers"
                     << ".tsv";
    string filename = filename_buf.str();
    string filename_csv = filename_csv_buf.str();
    string filename_tsv = filename_tsv_buf.str();

    // truncate or append
    string flags = first_time && IO_TruncateOutputFiles(cctkGH) ? "w" : "a";

    FILE *file = fopen(filename.c_str(), flags.c_str());
    FILE *file_csv = fopen(filename_csv.c_str(), flags.c_str());
    FILE *file_tsv = fopen(filename_tsv.c_str(), flags.c_str());
    if (file) {
      if (first_time) {
        fprintf(file, "# Clock %s\n", all_timers_clock);
        fprintf(file, "# Unit seconds\n");
        string all_timers_clock_csv = QuoteForCSV(all_timers_clock);
        fprintf(file_csv, "\"Clock %s\",", all_timers_clock_csv.c_str());
        fprintf(file_csv, "\"Unit seconds\"\n");
        string all_timers_clock_tsv = QuoteForTSV(all_timers_clock);
        fprintf(file_tsv, "Clock %s\t", all_timers_clock_tsv.c_str());
        fprintf(file_tsv, "Unit seconds\n");
      }
      // If the number of timers has changed, output the header containing the
      // timer names again.  It would be better to track creation and deletion
      // of timers in the flesh, but this method here works as long as timers
      // are never deleted.
      if (last_ntimers != timers.ntimers) {
        fprintf(file, "# Column 1 iteration\n");
        fprintf(file, "# Column 2 simulation time\n");
        fprintf(file, "# For all following columns:\n");
        fprintf(file, "#    Column 3n   average of all processes\n");
        fprintf(file, "#    Column 3n+1 minimum of all processes\n");
        fprintf(file, "#    Column 3n+2 maximum of all processes\n");
        fprintf(file_csv, "\"iteration\",");
        fprintf(file_csv, "\"simulation time\"");
        fprintf(file_tsv, "iteration\t");
        fprintf(file_tsv, "simulation time");
        for (int i = 0; i < timers.ntimers; i++) {
          const string &name = timers.names[i];
          fprintf(file, "# Column %d %s\n", 3 * i + 3, name.c_str());
          const string &name_csv = QuoteForCSV(name);
          fprintf(file_csv,
                  ",\"%s (average)\",\"%s (minimum)\",\"%s (maximum)\"",
                  name_csv.c_str(), name_csv.c_str(), name_csv.c_str());
          const string &name_tsv = QuoteForTSV(name);
          fprintf(file_tsv, "\t%s (average)\t%s (minimum)\t%s (maximum)",
                  name_tsv.c_str(), name_tsv.c_str(), name_tsv.c_str());
        }
        fprintf(file_csv, "\n");
        fprintf(file_tsv, "\n");
      }

      fprintf(file, "%d\t%.15g", cctk_iteration, double(cctk_time));
      fprintf(file_csv, "%d,%.15g", cctk_iteration, double(cctk_time));
      fprintf(file_tsv, "%d\t%.15g", cctk_iteration, double(cctk_time));

      for (int i = 0; i < timers.ntimers; i++) {
        fprintf(file, "\t%.15g %.15g %.15g", double(timers.secs_avg[i]),
                double(timers.secs_min[i]), double(timers.secs_max[i]));
        fprintf(file_csv, ",%.15g,%.15g,%.15g", double(timers.secs_avg[i]),
                double(timers.secs_min[i]), double(timers.secs_max[i]));
        fprintf(file_tsv, "\t%.15g\t%.15g\t%.15g", double(timers.secs_avg[i]),
                double(timers.secs_min[i]), double(timers.secs_max[i]));
      }
      fprintf(file, "\n");
      fprintf(file_csv, "\n");
      fprintf(file_tsv, "\n");
      fclose(file);
      fclose(file_csv);
      fclose(file_tsv);
    } else {
      CCTK_VWarn(
          1, __LINE__, __FILE__, CCTK_THORNSTRING,
          "Could not open timer report output files \"%s\", \"%s\", and \"%s\"",
          filename.c_str(), filename_csv.c_str(), filename_tsv.c_str());
    }

    first_time = false;
    last_ntimers = timers.ntimers;

  } // if root processor
}

void OutputAllTimersReadable(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  timer_stats timers;
  if (!CollectTimerInfo(cctkGH, timers))
    return;

  if (CCTK_MyProc(cctkGH) == 0) {
    static bool first_time = true;
    static int last_ntimers = -1;

    ostringstream filename_buf;
    filename_buf << out_dir << "/"
                 << "AllTimersReadable"
                 << ".txt";
    string filename = filename_buf.str();

    // truncate or append
    string flags = first_time && IO_TruncateOutputFiles(cctkGH) ? "w" : "a";

    FILE *file = fopen(filename.c_str(), flags.c_str());
    if (file) {
      if (first_time) {
        fprintf(file, "# Clock %s\n", all_timers_clock);
        fprintf(file, "# Unit seconds\n");
      }
      // If the number of timers has changed, output the header
      // containing the timer names again. It would be better to track
      // creation and deletion of timers in the flesh, but this method
      // here works as long as timers are never deleted.
      if (last_ntimers != timers.ntimers) {
        fprintf(file, "# Column 1 iteration\n");
        fprintf(file, "# Column 2 simulation time\n");
        fprintf(file, "# Column 3 timer number\n");
        fprintf(
            file,
            "# Column 4,5,6 average, minimum, maximum over all processors\n");
        fprintf(file, "# Column 7+ timer name\n");
      }
      for (int i = 0; i < timers.ntimers; i++) {
        const string &name = timers.names[i];
        fprintf(file, "%d %.15g\t%d\t%.15g %.15g %.15g\t%s\n", cctk_iteration,
                double(cctk_time), i, double(timers.secs_avg[i]),
                double(timers.secs_min[i]), double(timers.secs_max[i]),
                name.c_str());
      }
      fprintf(file, "\n");
      fclose(file);
    } else {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Could not open timer report output file \"%s\"",
                 filename.c_str());
    }

    first_time = false;
    last_ntimers = timers.ntimers;

  } // if root processor
}

void PrintTopTimers(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Collect timing information from all processes
  timer_stats timers;
  if (!CollectTimerInfo(cctkGH, timers))
    return;

  // Output the times only on the root process (because they are not reduced
  // onto the other processes anyway)
  if (CCTK_MyProc(cctkGH) != 0)
    return;

  // Sort timers (by average), largest times first
  vector<int> idx(timers.ntimers);
  for (int i = 0; i < timers.ntimers; ++i)
    idx[i] = i;
  assert(int(timers.secs_avg.size()) == timers.ntimers);
  sort(idx.begin(), idx.end(), [&](int ia, int ib) {
    // compare average times
    return timers.secs_avg[ia] > timers.secs_avg[ib];
  });

  CCTK_VInfo(CCTK_THORNSTRING, "Top timers at iteration %d time %g",
             cctk_iteration, double(cctk_time));

  printf("%s\n", sep.c_str());
  printf("%5s   %7s %7s %7s   %s (%s)\n", "%", "Time/s", "Min/s", "Max/s",
         "Timer", all_timers_clock);
  printf("%s\n", sep.c_str());

  // This should be the "CCTK total time" timer
  int total_idx = 0;
  assert(timers.ntimers > total_idx);

  CCTK_REAL max_time = timers.secs_max[idx[total_idx]];
  int digits = lrint(floor(log10(max_time) + 1.0e-4)) + 1;

  // Output timing results
  for (int i = 0; i < min(timers.ntimers, n_top_timers); ++i) {
    double percent =
        100.0 * timers.secs_avg[idx[i]] / timers.secs_avg[idx[total_idx]];

    // field widths: 5+3 + 7+1 + 7+1 + 7+3 + n = 80
    printf("%5.1f   %7.*f %7.*f %7.*f   %s\n", percent, 6 - digits,
           double(timers.secs_avg[idx[i]]), 6 - digits,
           double(timers.secs_min[idx[i]]), 6 - digits,
           double(timers.secs_max[idx[i]]), timers.names[idx[i]].c_str());
  }
  printf("%s\n", sep.c_str());
}

// Note: Timer names are truncated to 100 characters for simplicity
const int TIMERNAME_LENGTH = 101; // this includes the terminating NUL character
typedef array<char, TIMERNAME_LENGTH> timername_t;

// Collect timer information onto the root processor
int CollectTimerInfo(const cGH *restrict cctkGH, timer_stats &timers) {
  DECLARE_CCTK_PARAMETERS;

  // Gather number of timers from each process
  int myproc = CCTK_MyProc(cctkGH);
  int nprocs = CCTK_nProcs(cctkGH);

  int my_ntimers = CCTK_NumTimers();
  vector<int> all_ntimers;
  if (myproc == 0)
    all_ntimers.resize(nprocs);
#ifdef CCTK_MPI
  MPI_Gather(&my_ntimers, 1, MPI_INT, &all_ntimers[0], 1, MPI_INT, 0,
             MPI_COMM_WORLD);
#else
  assert(myproc == 0);
  assert(nprocs == 1);
  all_ntimers[0] = my_ntimers;
#endif
  int total_ntimers = 0;
  if (myproc == 0) {
    for (int p = 0; p < nprocs; ++p)
      total_ntimers += all_ntimers[p];
  }

  // Determine local timer names and their values
  vector<timername_t> my_timernames(my_ntimers);
  for (int n = 0; n < my_ntimers; ++n) {
    const char *name = CCTK_TimerName(n);
    if (name == nullptr)
      snprintf(&my_timernames[n][0], TIMERNAME_LENGTH, "DESTROYED TIMER %5d",
               n);
    else
      snprintf(&my_timernames[n][0], TIMERNAME_LENGTH, "%s", name);
  }
  vector<double> my_timervalues(my_ntimers);
  {
    cTimerData *td = CCTK_TimerCreateData();
    for (int n = 0; n < my_ntimers; ++n) {
      CCTK_TimerI(n, td);
      const cTimerVal *tv = CCTK_GetClockValue(all_timers_clock, td);
      if (!tv) {
        const char *name = CCTK_TimerName(n);
        if (name == nullptr)
          name = "(null)";
        CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Clock \"%s\" not found for timer #%d \"%s\"",
                   all_timers_clock, n, name);
        my_timervalues[n] = -1;
      } else {
        my_timervalues[n] = CCTK_TimerClockSeconds(tv);
      }
    }
    CCTK_TimerDestroyData(td);
  }

  // Gather timer names and values from each process
  vector<timername_t> all_timernames;
  vector<double> all_timervalues;
  if (myproc == 0) {
    all_timernames.resize(total_ntimers);
    all_timervalues.resize(total_ntimers);
  }
  vector<int> name_displacements, value_displacements, name_counts;
  if (myproc == 0) {
    name_displacements.resize(nprocs);
    value_displacements.resize(nprocs);
    name_counts.resize(nprocs);
    name_displacements[0] = 0;
    value_displacements[0] = 0;
    name_counts[0] = all_ntimers[0] * TIMERNAME_LENGTH;
    for (int p = 1; p < nprocs; ++p) {
      name_displacements[p] =
          name_displacements[p - 1] + all_ntimers[p - 1] * TIMERNAME_LENGTH;
      value_displacements[p] = value_displacements[p - 1] + all_ntimers[p - 1];
      name_counts[p] = all_ntimers[p] * TIMERNAME_LENGTH;
    }
  }
#ifdef CCTK_MPI
  MPI_Gatherv(&my_timernames[0], my_ntimers * TIMERNAME_LENGTH, MPI_CHAR,
              &all_timernames[0], &name_counts[0], &name_displacements[0],
              MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Gatherv(&my_timervalues[0], my_ntimers, MPI_DOUBLE, &all_timervalues[0],
              &all_ntimers[0], &value_displacements[0], MPI_DOUBLE, 0,
              MPI_COMM_WORLD);
#else
  assert(all_timernames.size() == my_timernames.size());
  assert(all_timervalues.size() == my_timervalues.size());
  copy(my_timernames.begin(), my_timernames.end(), all_timernames.begin());
  copy(my_timervalues.begin(), my_timervalues.end(), all_timervalues.begin());
#endif

  // Continue only on the root process
  if (myproc != 0) {
    return 1;
  }

  // Construct global list of timers: sort, then unique
  // TODO: sort the processes' timers separately (and in parallel), then merge
  // them
  vector<int> sort_index(total_ntimers);
  for (int i = 0; i < total_ntimers; ++i)
    sort_index[i] = i;
  sort(sort_index.begin(), sort_index.end(),
       [&](int ia, int ib) { return all_timernames[ia] < all_timernames[ib]; });
  sort_index.erase(
      unique(sort_index.begin(), sort_index.end(), [&](int ia, int ib) {
        return all_timernames[ia] == all_timernames[ib];
      }), sort_index.end());
  int unique_timers = sort_index.size();

  // Allocate timer data structure
  timers.ntimers = unique_timers;
  timers.names.resize(timers.ntimers);
  timers.secs_avg.resize(timers.ntimers);
  timers.secs_min.resize(timers.ntimers);
  timers.secs_max.resize(timers.ntimers);
  for (int n = 0; n < timers.ntimers; ++n) {
    timers.names[n] = &all_timernames[sort_index[n]][0];

    // Reduce timer values
    CCTK_REAL count = 0.0;
    CCTK_REAL sum = 0.0;
    CCTK_REAL minval = numeric_limits<CCTK_REAL>::infinity();
    CCTK_REAL maxval = 0.0;
    // Reduce over all processes
    for (int p = 0; p < nprocs; ++p) {
      int name_offset = name_displacements[p] / TIMERNAME_LENGTH;
      // Look for this timer
      // TODO: use a map
      for (int i = 0; i < all_ntimers[p]; ++i) {
        if (timers.names[n] == &all_timernames[name_offset + i][0]) {
          // Found the timer
          CCTK_REAL value = all_timervalues[value_displacements[p] + i];
          count += 1;
          sum += value;
          minval = min(minval, value);
          maxval = max(maxval, value);
          break;
        }
      }
      // Ignore timers that do not exist on this process
    }
    assert(count > 0);
    timers.secs_avg[n] = sum / count;
    timers.secs_min[n] = minval;
    timers.secs_max[n] = maxval;
  }

  return 1;
}

// Quote a string so that it can be output as CSV entry.
string QuoteForCSV(const string &str) {
  ostringstream buf;
  // Begin with a quote
  buf << '"';
  // Copy string into result, quoting as necessary
  for (string::const_iterator it = str.begin(); it != str.end(); ++it) {
    char ch = *it;
    // Quote all quote characters by doubling them (so <"hi" there> becomes
    // <""hi"" there>
    if (ch == '"')
      buf << '"';
    buf << ch;
  }
  // End with a quote
  buf << '"';
  return buf.str();
}

// Quote a string so that it can be output as TSV entry.
string QuoteForTSV(const string &str) {
  ostringstream buf;
  // Copy string into result, replacing tabs as necessary
  for (string::const_iterator it = str.begin(); it != str.end(); ++it) {
    char ch = *it;
    if (ch == '\t') {
      buf << ' ';
    } else {
      buf << ch;
    }
  }
  return buf.str();
}
}
