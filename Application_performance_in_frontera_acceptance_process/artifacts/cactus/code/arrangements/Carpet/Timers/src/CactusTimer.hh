#ifndef CACTUSTIMER_HH
#define CACTUSTIMER_HH

#include <cctk.h>

#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>

namespace Timers {

/** The CactusTimer class wraps the Cactus timer mechanism. All
    times are returned as doubles for now. */

class CactusTimer {
  int handle;
  bool running;

public:
  /// Create a new Cactus timer with the given name
  CactusTimer(std::string timername);

  /// Destruct the timer
  ~CactusTimer();

  /// Start the timer
  void start();

  /// Stop the timer
  void stop();

  /// Reset the timer
  void reset();

  /// Timer name
  std::string name() const;

  /// Return the current time of the timer in seconds as a double
  double getTime();

  /// Return the average and maximum current time over all MPI processes
  void getGlobalTime(double &avg, double &max);

  /// Return all clock names and their units
  std::vector<std::pair<std::string, std::string> > getAllTimerNames() const;

  /// Return all clock values of the timer as double
  std::vector<double> getAllTimerValues();

  /// Print timer data
  void printData();

  std::ostream &serialise(std::ostream &os);

private:
  // Output (debug) messages that a timer is starting or stopping
  void msgCreate() const;
  void msgStart() const;
  void msgStop() const;
  void msgRead(double val) const;
};

} // namespace Timers

#endif // CACTUSTIMER_HH
