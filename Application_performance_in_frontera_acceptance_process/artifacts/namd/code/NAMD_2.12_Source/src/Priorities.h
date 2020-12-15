/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PRIORITIES_H
#define PRIORITIES_H

// pass PRIORITY_SIZE as last argument to new when allocating message
// e.g.,  MyMsg *msg = new (len1, len2, PRIORITY_SIZE) MyMsg;

#define PRIORITY_SIZE ((int) sizeof(int)*8)

// total priority is sequence + type + patch
// always use lowest (most urgent) patch priority applicable

#define SET_PRIORITY(MSG,SEQ,PRIO) { \
  CkSetQueueing(MSG, CK_QUEUEING_IFIFO); \
  *((int*) CkPriorityPtr(MSG)) = (((SEQ)&0xffff)<<15) + (PRIO); }
// sequence priority is 16 bits shifted by 15 to leave sign bit 0

// patch priorities in range [1,252]
// reserve 0 for tuples, use prime to randomize neighbors
#define PATCH_PRIORITY(PID) (((PID)%251)+1)

// the following are in order of decreasing urgency

#define PME_PRIORITY (2<<8)
#define PME_GRID_PRIORITY (PME_PRIORITY+1)
#define PME_TRANS_PRIORITY (PME_PRIORITY+2)
#define PME_TRANS2_PRIORITY (PME_PRIORITY+3)
#define PME_UNTRANS_PRIORITY (PME_PRIORITY+4)
#define PME_UNTRANS2_PRIORITY (PME_PRIORITY+5)

#define MSM_PRIORITY PME_PRIORITY
 
#if defined(NAMD_CUDA) || defined(NAMD_MIC)
// higher priority so offloaded work can overlap
#define PROXY_DATA_PRIORITY (1<<8)
#define PME_OFFLOAD_PRIORITY 0
#define PME_OFFLOAD_UNGRID_PRIORITY (3<<8)
#else
#define PROXY_DATA_PRIORITY (3<<8)
#endif

//used in HomePatch::positionsReady
//gbis positions distributed with normal PROXY_DATA_PRIORITY
//compute priorities are inited in ComputePatch* and added in Compute.C

//use in Compute::patchReady              DONE
#define GB1_COMPUTE_PROXY_PRIORITY (4<<8)
//use in ProxyPatch::boxClosed            DONE
#define GB1_PROXY_RESULTS_PRIORITY (5<<8)
//use in Compute::patchReady              DONE
#define GB1_COMPUTE_HOME_PRIORITY (6<<8)
//used in HomePatch::gbisP2Ready          DONE
#define GB2_PROXY_DATA_PRIORITY (7<<8)
//use in Compute::gbisP2PatchReady        DONE
#define GB2_COMPUTE_PROXY_PRIORITY (8<<8)
//use in ProxyPatch::boxClosed            DONE
#define GB2_PROXY_RESULTS_PRIORITY (9<<8)
//use in Compute::patchReady              DONE
#define GB2_COMPUTE_HOME_PRIORITY (10<<8)
//used in HomePatch::gbisP3Ready          DONE
#define GB3_PROXY_DATA_PRIORITY (11<<8)


// from here on GB computes use normal compute priorities
//use in Compute::gbisP3PatchReady        DONE
#define COMPUTE_PROXY_PRIORITY (12<<8)
//use in ProxyPatch::send
#define PROXY_RESULTS_PRIORITY (13<<8) // DONE
#define PME_UNGRID_PRIORITY (14<<8)
//use in Compute::patchReady              DONE
#define COMPUTE_HOME_PRIORITY (15<<8)
//end gbis

#endif // PRIORITIES_H

