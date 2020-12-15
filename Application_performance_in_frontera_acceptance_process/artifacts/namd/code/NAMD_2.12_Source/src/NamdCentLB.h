/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/NamdCentLB.h,v $
 * $Author: jim $
 * $Date: 2013/08/22 15:17:18 $
 * $Revision: 1.19 $
 *****************************************************************************/

/** \file NamdCentLB.h
 *  This class decides which load balancers to invoke.
 *  By default if nothing is specified in the config file,
 *  AlgSeven and RefineOnly are used. This is when each gets
 *  invoked:
 *
 *  If we choose firstLdbStep = 20 and ldbPeriod = 100, we have the
 *  following timeline (for these particular numbers): 
 *           
 *  Tracing	    :  <---- off ----><------------ on -----------><-- off
 *  ldbCycleNum     :  1    2    3    4    5        6     7        8     9
 *  Iteration Steps : 00===20===40===60===80======160===180======260===280
 *  Ldb Step() No   :            1    2    3        4     5        6     7
 *  Ldb Strategy    :         Alg7  Ref  Ref     Inst   Ref     Inst   Ref
 *                            Alg7					  
 *									  
 *  Alg7 = AlgSeven							  
 *  Ref  = Refine (NamdCentLB.C, Rebalancer.C)				  
 *  Inst = Instrumentation Phase (no real load balancing)		  
 * 
 *  This class also fills in the NAMD load balancing data structures
 *  from the Charm LDBStats
 */


#ifndef _NAMDCENTLB_H_
#define _NAMDCENTLB_H_

#include <CentralLB.h>
#include "NamdCentLB.decl.h"

#include "Node.h"
#include "PatchMap.h"
#include "SimParameters.h"
#include "RefineOnly.h"
#include "Alg7.h"
#include "AlgRecBisection.h"
#include "InfoStream.h"
#include "TorusLB.h"
#include "RefineTorusLB.h"

void CreateNamdCentLB();
NamdCentLB *AllocateNamdCentLB();

class NamdCentLB : public CentralLB {

public:
  NamdCentLB();
  NamdCentLB(CkMigrateMessage *);
  CLBMigrateMsg* Strategy(LDStats* stats);

private:
  bool QueryBalanceNow(int step);
  bool QueryDumpData();
  int buildData(LDStats* stats);
  int requiredProxies(PatchID id, int neighborNodes[]);
#if USE_TOPOMAP 
  int requiredProxiesOnProcGrid(PatchID id, int neighborNodes[]);
#endif
  void dumpDataASCII(char *file, int numProcessors, int numPatches,
		int numComputes);
  void loadDataASCII(char *file, int &numProcessors, int &numPatches,
		int &numComputes);

  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;
};

#endif /* _NAMDCENTLB_H_ */
