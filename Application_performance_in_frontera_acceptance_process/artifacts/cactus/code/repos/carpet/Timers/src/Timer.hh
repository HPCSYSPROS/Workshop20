#ifndef TIMER_HH
#define TIMER_HH

#include <iostream>
#include <list>

#include <cctk.h>

#include <TimerTree.hh>

namespace Timers {

/**
   This class allows the user to instrument their code with named
   timers which can then be later queried to determine the amount of
   time spent in the code between "start" and "end" calls. The
   sequence of start and end calls of different timers determines a
   dynamical hierarchical tree structure implemented by the
   TimerNode class.

   To use this class, create a timer object with a particular name:

      Timer timer("MyTimer")

   Now wrap the code to be timed with start() and stop() calls:

      timer.start()

   some code

      timer.stop()

   You can start and stop a timer multiple times. The timer will be
   created as a child of whatever timer is current (i.e. has been
   started and not stopped) at the time of the first start() call.
   Any timers which are started between the start() and stop(),
   whether or not they are in the same file, will be stored as
   children of this timer.

   Timer objects must be started and stopped in a non-overlapping
   manner. Specifically, a timer cannot be stopped if it is not the
   most recently started timer. Doing so will generate an error.

   Timer objects can be allocated as "static" or not - it does not
   matter.
*/

class Timer {

public:
  Timer(std::string name, int tree = 0);
  ~Timer();

  void instantiate();
  void start();
  void stop();
  std::string name() const;
  double getTime();

  static void outputTree(std::string name);
  static void outputTreeXML();

private:
  std::string d_name;
  TimerTree *d_tree;
};

} // namespace Timers

// Macros for using timers in a convenient manner

#define TIMING_BEGIN(name)                                                     \
  do {                                                                         \
    static Timers::timer(name);                                                \
    timer.start();                                                             \
    {

#define TIMING_END                                                             \
  }                                                                            \
  timer.stop();                                                                \
  }                                                                            \
  while (0)

#endif // TIMER_HH
