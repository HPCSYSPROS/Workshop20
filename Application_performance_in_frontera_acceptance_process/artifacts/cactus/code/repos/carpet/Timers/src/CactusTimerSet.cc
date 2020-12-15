#include <cassert>
#include <cstdio>
#include <cstring>
#include <list>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <util_String.h>

#if HAVE_UNISTD_H
#include <fcntl.h>
#include <unistd.h>
#endif

#include <defs.hh>

#include <CactusTimer.hh>
#include <CactusTimerSet.hh>
#include <Timer.hh>
#include <TimerTree.hh>

namespace Timers {

using namespace std;

// A global timer set
CactusTimerSet timerSet;

// Add a timer
void CactusTimerSet::add(CactusTimer *const timer) {
  DECLARE_CCTK_PARAMETERS;

  if (disable_cactus_timer_set)
    return;
  timers.insert(timer);
}

// Remove a timer
void CactusTimerSet::remove(CactusTimer *const timer) {
  DECLARE_CCTK_PARAMETERS;

  if (disable_cactus_timer_set)
    return;
  timers.erase(timer);
}

// Print all timer names
void CactusTimerSet::printNames() const {
  printf("Timer names:\n");
  int n = 0;
  for (timers_t::const_iterator itimer = timers.begin(); itimer != timers.end();
       ++itimer) {
    printf("   [%4d] %s\n", n, (*itimer)->name().c_str());
    ++n;
  }
}

// Print all timer data
void CactusTimerSet::printData() {
  for (timers_t::const_iterator itimer = timers.begin(); itimer != timers.end();
       ++itimer) {
    (*itimer)->printData();
  }
  printf("\n");
}

// Print all timer data
void CactusTimerSet::writeData(const cGH *const cctkGH,
                               const char *const filename) {
  const int oldfd = redirect(cctkGH, filename);
#if 0
    printf("********************************************************************************\n");
#endif
  printf("# Carpet timing information at iteration %d time %g:\n",
         cctkGH->cctk_iteration, double(cctkGH->cctk_time));
  timerSet.printData();
  unredirect(oldfd);
}

// If filename is not empty, then redirect stdout to a file
int CactusTimerSet::redirect(const cGH *const cctkGH,
                             const char *const filename) {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(filename, "")) {
    return -1;
  }

#ifndef HAVE_UNISTD_H
  CCTK_WARN(CCTK_WARN_ALERT, "Cannot redirect timer output to a file; the "
                             "operating system does not support this");
  return -1;
#else

  const int myproc = CCTK_MyProc(cctkGH);
  char fullname[10000];
  Util_snprintf(fullname, sizeof fullname, "%s/%s.%04d.txt", out_dir, filename,
                myproc);

  int flags = O_WRONLY | O_CREAT | O_APPEND; // append
  static bool first_time = true;
  if (first_time) {
    first_time = false;
    if (IO_TruncateOutputFiles(cctkGH)) {
      flags = O_WRONLY | O_CREAT | O_TRUNC; // truncate
    }
  }

  // Temporarily redirect stdout
  fflush(stdout);
  const int oldfd = dup(1); // fd 1 is stdout
  const int mode = 0644;    // rw-r--r--, or a+r u+w
  const int fdfile = open(fullname, flags, mode);
  if (fdfile < 0) {
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Could not open timer output file \"%s\"", fullname);
    close(oldfd);
    return -1;
  }
  // close(1);
  // const int fd = dup(fdfile); // dup to 1, i.e., stdout again
  const int fd = dup2(fdfile, 1); // dup to 1, i.e., stdout again
  assert(fd == 1);
  close(fdfile);
  return oldfd;
#endif
}

// Redirect stdout back
void CactusTimerSet::unredirect(const int oldfd) {
  if (oldfd < 0)
    return;

#ifdef HAVE_UNISTD_H
  fflush(stdout);
  // close(1);
  // const int fd = dup(oldfd);
  const int fd = dup2(oldfd, 1);
  if (not(fd == 1)) {
    fprintf(stderr, "oldfd=%d fd=%d\n", oldfd, fd);
  }
  assert(fd == 1);
  close(oldfd);
#endif
}

#if 0
  /// Reduce each timer in the set across all processes and update
  /// each timer with the reduction information.
  void CactusTimerSet::reduce()
  {
    // Collect timer names that each process has
    
    // Construct union of all timer names, sort canonically and assign
    // integer identifiers
    
    // For each timer, identify which processes have that timer
    
    // Reduce the timer across all those processes (return to root proc only)
    
    serialise(cout);
  }
  
  ostream& CactusTimerSet::serialise(ostream &os)
  {
    for (timers_t::const_iterator
           itimer = timers.begin(); itimer != timers.end(); ++itimer)
    {
      (*itimer)->serialise(os);
      os << endl;
    }
    return os;
  }
#endif

} // namespace Carpet
