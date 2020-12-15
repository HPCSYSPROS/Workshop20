/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/* This class is just used to test out the GlobalMaster computing
   objects. */

#ifndef GLOBALMASTERTEST_H
#define GLOBALMASTERTEST_H

class GlobalMasterTest : public GlobalMaster {
 protected:
  virtual void calculate();
};

#endif
