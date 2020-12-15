#ifndef COMPUTECONSFORCE_H
#define COMPUTECONSFORCE_H

#include "ComputeHomePatch.h"
#include "ReductionMgr.h"

class ComputeConsForce : public ComputeHomePatch
{
public:
  ComputeConsForce(ComputeID, PatchID);
  virtual ~ComputeConsForce();
  virtual void doForce(FullAtom*, Results*);
private:
  SubmitReduction *reduction;
};

class ComputeConsTorque : public ComputeHomePatch
{
public:
  ComputeConsTorque(ComputeID, PatchID);
  virtual ~ComputeConsTorque();
  virtual void doForce(FullAtom*, Results*);
private:
  SubmitReduction *reduction;
};

#endif
