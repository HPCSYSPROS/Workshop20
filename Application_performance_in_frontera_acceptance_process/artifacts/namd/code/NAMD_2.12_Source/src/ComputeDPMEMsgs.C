/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeDPMEMsgs.h"
#include "packmsg.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

// DATA MESSAGE

ComputeDPMEDataMsg::ComputeDPMEDataMsg(void) { 
  numParticles = 0;
  particles = 0;
}

ComputeDPMEDataMsg::~ComputeDPMEDataMsg(void) { 
  delete [] particles;
}

PACK_MSG(ComputeDPMEDataMsg,
  PACK(node);
  PACK_AND_NEW_ARRAY(particles,numParticles);
)


// RESULTS MESSAGE

ComputeDPMEResultsMsg::ComputeDPMEResultsMsg(void) { 
  numParticles = 0;
  forces = 0;
}

ComputeDPMEResultsMsg::~ComputeDPMEResultsMsg(void) { 
  delete [] forces;
}

PACK_MSG(ComputeDPMEResultsMsg,
  PACK(node);
  PACK_AND_NEW_ARRAY(forces,numParticles);
)

