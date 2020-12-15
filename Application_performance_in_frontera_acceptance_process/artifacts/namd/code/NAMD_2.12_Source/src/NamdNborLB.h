#ifndef _NAMDNBORLB_H_
#define _NAMDNBORLB_H_

#include "NeighborLB.h"
#include "NamdNborLB.decl.h"

#include "Node.h"
#include "PatchMap.h"
#include "SimParameters.h"
#include "AlgNbor.h"
#include "InfoStream.h"

void CreateNamdNborLB();

class NamdNborLB : public NeighborLB {

public:
  NamdNborLB();
private:
  int act;
  int numNbors;
private:
  int max_neighbors();
  int num_neighbors();
  void neighbors(int* _n);
  bool QueryBalanceNow(int step);
  bool QueryMigrateStep(int _step);
  NLBMigrateMsg* Strategy(NborBaseLB::LDStats* stats, int count);
  int buildData(NborBaseLB::LDStats* stats, int count);
  int requiredProxies(PatchID id, int neighborNodes[]);

  int ldbNum;
  computeInfo *computeArray;
  patchInfo *patchArray;
  processorInfo *processorArray;
};

#endif /* _NAMDCENTLB_H_ */
