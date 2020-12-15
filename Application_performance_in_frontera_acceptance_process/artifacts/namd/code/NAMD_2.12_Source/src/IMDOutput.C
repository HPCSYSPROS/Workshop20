/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "IMDOutput.h"
#include "GlobalMasterIMD.h"

IMDOutput::IMDOutput() {
  imd = NULL;
  transrate = 1;
  ignore = 0;
}

IMDOutput::~IMDOutput() {
}

void IMDOutput::use_imd(GlobalMasterIMD *g) {
  imd = g;
  ignore = g->IMDignore;
}

void IMDOutput::gather_energies(IMDEnergies *energies) { 
  if (!imd || energies->tstep % transrate) return;
  imd->send_energies(energies);
}

void IMDOutput::gather_coordinates(int timestep, int N, FloatVector *coords) {
  if ( ignore ) {
    imd->step = timestep;
    imd->calculate();
  }
  if (!imd || timestep % transrate) return;
  imd->send_fcoords(N, coords);
}

