/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "charm++.h"

#ifndef _BCASTCLI_H
#define _BCASTCLI_H

class BroadcastClient {
public:
  int id;
  BroadcastClient(int id);
  ~BroadcastClient();
  void awaken(int id, int tag);

protected: 
  void suspendFor(int tag);

  int suspended;
  int waitForTag;
  CthThread thread;
};

#endif

