/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/NamdHybridLB.h,v $
 * $Author: jim $
 * $Date: 2013/08/29 02:20:58 $
 * $Revision: 1.15 $
 *****************************************************************************/

#ifndef _NAMDHYBRIDLB_H_
#define _NAMDHYBRIDLB_H_

#include <HybridBaseLB.h>
#include "NamdHybridLB.decl.h"

#include "Node.h"
#include "PatchMap.h"
#include "SimParameters.h"
#include "RefineOnly.h"
#include "Alg7.h"
#include "AlgRecBisection.h"
#include "InfoStream.h"
#include "NamdCentLB.h"
#include "NamdDummyLB.h"
#include "TorusLB.h"
#include "RefineTorusLB.h"

void CreateNamdHybridLB();

/**
 * Class for message containing CPU loads.
 */
class LocalLBInfoMsg: public CMessage_LocalLBInfoMsg{
public:
	int n_moves;
	int startPE;
	int endPE;
	MigrateInfo *moves;
	double *cpuloads;

	// Constructor
	LocalLBInfoMsg(): n_moves(0), startPE(0), endPE(0){}

	// Pup method
	void pup(PUP::er &p) {
		int i;
		p | n_moves;
		p | startPE;
		p | endPE;
		for (i=0; i<n_moves; ++i) p | moves[i];
		for (i=0; i<endPE-startPE+1; ++i) p | cpuloads[i];
	}

};

class SplitComputesMsg;

class NamdHybridLB : public HybridBaseLB {

public:
  NamdHybridLB();
  NamdHybridLB(CkMigrateMessage *m):HybridBaseLB(m) {}
  void UpdateLocalLBInfo(LocalLBInfoMsg *msg);
  //void CollectInfo(Location *loc, int n, int fromlevel);
  void splitComputes(SplitComputesMsg*);

private:
  CProxy_NamdHybridLB thisProxy;
  int updateCount;
  int splitCount;
  SplitComputesMsg **splitComputesMsgs;
  bool collectFlag;
  bool updateFlag;
  int parent_backup;
  Location *loc_backup;
  int n_backup;
  int fromlevel_backup;

  int *from_procs;
  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;

	double *peLoads;
	int startPE;
	int endPE;

  bool QueryBalanceNow(int step);
  bool QueryDumpData();
  // LBVectorMigrateMsg* VectorStrategy(LDStats* stats);

  CLBMigrateMsg* Strategy(LDStats* stats);
  LBMigrateMsg* GrpLevelStrategy(LDStats* stats);
  
  int buildData(LDStats* stats);
  int requiredProxies(PatchID id, int neighborNodes[]);
  void dumpDataASCII(char *file, int numProcessors, int numPatches,
                int numComputes);

  // centralized load balancer for load balancing all the children processors
  NamdCentLB *centralLB;
  NamdDummyLB *dummyLB;
};

#endif /* _NAMDHYBRIDLB_H_ */
