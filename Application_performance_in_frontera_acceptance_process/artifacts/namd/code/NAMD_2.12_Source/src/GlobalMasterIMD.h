/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef GLOBALMASTERIMD_H
#define GLOBALMASTERIMD_H

#include "GlobalMaster.h"
#include "imd.h"
#include "ResizeArray.h"

class FloatVector;

class GlobalMasterIMD : public GlobalMaster {
 public: 
  /* initializes this according to the simulation parameters */
  GlobalMasterIMD();
  ~GlobalMasterIMD();

  void send_energies(IMDEnergies *);
  void send_fcoords(int, FloatVector *);

 protected:

  friend class IMDOutput;

  virtual void calculate();

  // Simple function for getting MDComm-style forces from VMD
  void get_vmd_forces();

  // flag for whether to proceed with simulation when there are no connections
  int IMDwait;

  // flag for whether to ignore forces
  int IMDignore;

  // My server socket handle
  void *sock;

  // Connected sockets
  ResizeArray<void *>clients;

  // temporaries in case 3*sizeof(float) != sizeof(FloatVector)
  float *coordtmp;
  int coordtmpsize;
};

#endif

