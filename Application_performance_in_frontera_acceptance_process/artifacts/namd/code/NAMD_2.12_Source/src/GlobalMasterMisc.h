/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#ifndef GLOBALMASTERMISC_H
#define GLOBALMASTERMISC_H

#include "ComputeHomePatches.h"
#include "GlobalMaster.h"
#include "GlobalMasterEasy.h"
#include "NamdTypes.h"

class GlobalMasterMisc : public GlobalMasterEasy {
public:
  GlobalMasterMisc();
  virtual ~GlobalMasterMisc();
protected:

  virtual void easy_init(const char *);
  virtual void easy_calc(void);

  Vector originalPosition;
  bool firstTime;
};

#endif
