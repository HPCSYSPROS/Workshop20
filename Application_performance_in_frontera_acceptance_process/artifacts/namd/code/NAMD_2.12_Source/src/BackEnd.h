/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef BACKEND_H
#define BACKEND_H

/*  Base class for providing an API to a front end interface.  */

class BackEnd {
public:

  static void init(int argc, char **argv);  // Must call at program startup
  static void exit(void);  // Must call at program shutdown

  static void suspend(void);  // Let message code run...
  static void awaken(void);  // ...then return
  static void barrier(void);  // Let message code run, then return
  static void ExitSchedOn(int pe); // Tell a remote processor to stop
                                   // his scheduler

};

#endif

