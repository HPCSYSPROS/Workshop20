#ifndef CACTUSTIMERSET_HH
#define CACTUSTIMERSET_HH

#include <iostream>
#include <set>

#include <cctk.h>
#include "CactusTimer.hh"

namespace Timers {

class CactusTimerSet;
extern CactusTimerSet timerSet;

using namespace std;

// A set of timers
class CactusTimerSet {

  typedef set<CactusTimer *> timers_t;
  timers_t timers;

public:
  // Add a timer
  void add(CactusTimer *timer);

  // Remove a timer
  void remove(CactusTimer *timer);

  // Print all timer names
  void printNames() const;

  // Print all timer data
  void printData();

  // Write all timer data of the global timer set to a file
  static void writeData(cGH const *cctkGH, char const *filename);

#if 0
    // Reduce each timer in the set across all processes and update
    // each timer with the reduction information.
    void reduce();
    
    ostream& serialise(ostream &os);
#endif

private:
  // If filename is not empty, then redirect stdout to a file,
  // returning the old stdout file descriptor
  static int redirect(cGH const *cctkGH, char const *filename);

  // Redirect stdout back
  static void unredirect(int oldfd);
}; // class CactusTimerSet

} // namespace Timers

#endif // CACTUSTIMERSET_HH
