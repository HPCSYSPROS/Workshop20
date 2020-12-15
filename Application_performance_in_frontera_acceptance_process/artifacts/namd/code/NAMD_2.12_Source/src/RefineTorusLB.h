/** \file RefineTorusLB.h 
 *  Author: Abhinav S Bhatele
 *  Date Created: June 12th, 2007
 *
 *  Replacement for RefineOnly.h
 */

#ifndef _REFINE_TORUSLB_H_
#define _REFINE_TORUSLB_H_

#include "Rebalancer.h"

class RefineTorusLB : public Rebalancer
{
  private:
    pcpair *bestPe[6];
    pcpair *goodPe[6];

    void strategy();
    void selectPes(processorInfo *p, computeInfo *c);

  public:
    RefineTorusLB(computeInfo *cs, patchInfo *pas, processorInfo *pes, int ncs, 
int npas, int npes, int flag);
    ~RefineTorusLB();
    void binaryRefine();
    int newRefine();

}; 

#endif
