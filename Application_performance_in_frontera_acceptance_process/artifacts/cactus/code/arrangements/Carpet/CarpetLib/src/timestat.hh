#ifndef TIMESTAT_HH
#define TIMESTAT_HH

#include <iostream>
#include <list>
#include <string>

#include <cctk.h>

#include <cycleclock.h>
#ifndef HAVE_TICK_COUNTER
typedef int ticks;
inline ticks getticks() { return 0; }
inline double elapsed(ticks, ticks) { return 0.0; }
inline double seconds_per_tick(void) { return 0.0; }
inline void measure_tick(void) {}
#endif

namespace CarpetLib {

using namespace std;

class Timer;

// A set of timers
class TimerSet {

  list<Timer *> timers;

public:
  // Add a timer
  void add(Timer *timer);

  // Remove a timer
  void remove(Timer *timer);

  // Output all timer names
  void outputNames(ostream &os) const;

  // Output all timer data
  void outputData(ostream &os) const;

}; // class TimerSet

inline ostream &operator<<(ostream &os, TimerSet const &timerSet) {
  timerSet.outputData(os);
  return os;
}

// A timer, which counts time (in seconds) spent in and amount (in
// bytes) used in various operations
class Timer {

  string timername;

public:
  // Create a new timer with the given name
  Timer(char const *timername_);

  // Destroy a timer
  ~Timer();

private:
  // Reset the statistics
  void resetstats();

  // Add statistics of a timing operation
  void addstat(double t, double b);

private:
  double wtime;
  double wtime2;
  double wmin;
  double wmax;

  double bytes;
  double bytes2;
  double bmin;
  double bmax;

  double count;

  bool running;
  ticks starttime;

public:
  // Start the timer
  void start();

  // Stop the timer
  void stop(double b);

  // Reset the timer
  void reset();

  // Timer name
  string name() const CCTK_MEMBER_ATTRIBUTE_PURE;

  // Print timer data
  void outputData(ostream &os) const;
};

inline ostream &operator<<(ostream &os, Timer const &timer) {
  timer.outputData(os);
  return os;
}

} // namespace CarpetLib

#endif // TIMESTAT_HH
