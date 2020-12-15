/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef IMD_OUTPUT_H__
#define IMD_OUTPUT_H__

#include "imd.h"
class FloatVector;
class GlobalMasterIMD;

// IMDOutput
// This object's only reason for existence is to forward energies and 
// coordinates to GlobalMasterIMD for export to connected IMD connections.
// If Controller could access GlobalMasterIMD directly we wouldn't need
// need this.

class IMDOutput {

public:
  IMDOutput();   
  ~IMDOutput();

  // The GlobalMasterIMD instance passes itself to this object so it can
  // receive energies and coordinates.
  void use_imd(GlobalMasterIMD *);

  // gather_* are called by Controller with the current timesteps and
  // energies.
  void gather_energies(IMDEnergies *energies); 
  void gather_coordinates(int timestep, int N, FloatVector *coords);

  // called by GlobalMasterIMD to set the transfer rate.  Should probably
  // be handled internally by GlobalMasterIMD instead.
  void set_transrate(int newrate) {transrate = newrate; }

private:
  GlobalMasterIMD *imd;
  int transrate;
  int ignore;
};

#endif

