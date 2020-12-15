/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeGlobalMsgs.h"
#include "packmsg.h"

//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#include "ComputeMgr.decl.h"


// CONFIG MESSAGE

#if 0
ComputeGlobalConfigMsg::ComputeGlobalConfigMsg(void) { 
}

ComputeGlobalConfigMsg::~ComputeGlobalConfigMsg(void) { 
}

PACK_MSG(ComputeGlobalConfigMsg,
  PACK_RESIZE(aid);
  PACK_RESIZE(gdef);
)
#endif


// DATA MESSAGE

ComputeGlobalDataMsg::ComputeGlobalDataMsg(void) { 
}

ComputeGlobalDataMsg::~ComputeGlobalDataMsg(void) { 
}

PACK_MSG(ComputeGlobalDataMsg,
  PACK(step);
  PACK(count);
  PACK_RESIZE(aid);
  PACK_RESIZE(p);
  PACK_RESIZE(gcom);
  PACK_RESIZE(gmass);
  PACK_RESIZE(fid);
  PACK_RESIZE(tf);
  PACK_RESIZE(gtf);
  PACK_RESIZE(lat);
)


// RESULTS MESSAGE

ComputeGlobalResultsMsg::ComputeGlobalResultsMsg(void) { 
  reconfig = 0;
  resendCoordinates = 0;
}

ComputeGlobalResultsMsg::~ComputeGlobalResultsMsg(void) { 
}

PACK_MSG(ComputeGlobalResultsMsg,
  PACK_RESIZE(aid);
  PACK_RESIZE(f);
  PACK_RESIZE(gforce);
  PACK(seq);
  PACK(totalforces);
  PACK(reconfig);
  PACK(resendCoordinates);
  if ( packmsg_msg->reconfig ) {
    PACK_RESIZE(newaid);
    PACK_RESIZE(newgdef);
  }
)
