/** \file TorusLB.h 
 *  Author: Abhinav S Bhatele
 *  Date Created: June 05th, 2007
 *
 *  Replacement for Alg7.h
 */

#ifndef _TORUSLB_H_
#define _TORUSLB_H_

#include "Rebalancer.h"
#include "RefineTorusLB.h"

class TorusLB : public RefineTorusLB
{
  private:
    processorInfo *bestPe[6];
    processorInfo *goodPe[6];
    processorInfo *badPe[6];

    void strategy();
    void selectPes(processorInfo *p, computeInfo *c);

  public:
    TorusLB(computeInfo *cs, patchInfo *pas, processorInfo *pes, int ncs, 
int npas, int npes);
    ~TorusLB();

}; 

#endif
