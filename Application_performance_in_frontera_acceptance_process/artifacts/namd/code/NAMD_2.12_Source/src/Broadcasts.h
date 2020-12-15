/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef BROADCASTS_H
#define BROADCASTS_H

#include "NamdTypes.h"
#include "Lattice.h"
#include "BroadcastObject.h"

enum {
  SCRIPT_END,
  SCRIPT_RUN,
  SCRIPT_CONTINUE,
  SCRIPT_OUTPUT,
  SCRIPT_FORCEOUTPUT,
  SCRIPT_MEASURE,
  SCRIPT_REINITVELS,
  SCRIPT_RESCALEVELS,
  SCRIPT_RELOADCHARGES,
  SCRIPT_CHECKPOINT,
  SCRIPT_REVERT,
  SCRIPT_CHECKPOINT_STORE,
  SCRIPT_CHECKPOINT_LOAD,
  SCRIPT_CHECKPOINT_SWAP,
  SCRIPT_CHECKPOINT_FREE,
  SCRIPT_ATOMSENDRECV,
  SCRIPT_ATOMSEND,
  SCRIPT_ATOMRECV,
  SCRIPT_MINIMIZE,
  SCRIPT_DUMMY
};

// Tags used in common by all users of broadcast system.
enum {
  velocityRescaleFactorTag,
  positionRescaleFactorTag,
  // For multigrator
  velocityRescaleTensorTag,
  velocityRescaleTensor2Tag,
  velocityRescaleFactor2Tag,
  positionRescaleFactor2Tag,
  // End multigrator
  tcoupleCoefficientTag,
  minimizeCoefficientTag,
  momentumCorrectionTag,
#if USE_BARRIER
  cycleBarrierTag,
#endif
  scriptBarrierTag,
  traceBarrierTag,
  accelMDRescaleFactorTag,
  adaptTemperatureTag, //Tag for adaptive tempering temperature updates to Sequencer
#ifdef MEASURE_NAMD_WITH_PAPI
  papiMeasureTag,
#endif
  dummyTag
};

// Broadcasts used by Contoller <-> Sequencer communication.
struct ControllerBroadcasts
{
  SimpleBroadcastObject<BigReal> velocityRescaleFactor;
  SimpleBroadcastObject<Tensor> positionRescaleFactor;
  // For multigrator
  SimpleBroadcastObject<Tensor> velocityRescaleTensor;
  SimpleBroadcastObject<Tensor> velocityRescaleTensor2;
  SimpleBroadcastObject<BigReal> velocityRescaleFactor2;
  SimpleBroadcastObject<Tensor> positionRescaleFactor2;
  // End multigrator
  SimpleBroadcastObject<BigReal> tcoupleCoefficient;
  SimpleBroadcastObject<BigReal> minimizeCoefficient;
  SimpleBroadcastObject<Vector> momentumCorrection;
#if USE_BARRIER
  SimpleBroadcastObject<int> cycleBarrier;
#endif
  SimpleBroadcastObject<int> scriptBarrier;
  SimpleBroadcastObject<int> traceBarrier;
  SimpleBroadcastObject<Vector> accelMDRescaleFactor;
  SimpleBroadcastObject<BigReal> adaptTemperature; 
#ifdef MEASURE_NAMD_WITH_PAPI
  SimpleBroadcastObject<int> papiMeasureBarrier;
#endif

  ControllerBroadcasts(const LDObjHandle *ldObjPtr = 0) : 
    velocityRescaleFactor(velocityRescaleFactorTag, ldObjPtr),
    positionRescaleFactor(positionRescaleFactorTag, ldObjPtr),
    // For multigrator
    velocityRescaleTensor(velocityRescaleTensorTag, ldObjPtr),
    velocityRescaleTensor2(velocityRescaleTensor2Tag, ldObjPtr),
    velocityRescaleFactor2(velocityRescaleFactor2Tag, ldObjPtr),
    positionRescaleFactor2(positionRescaleFactor2Tag, ldObjPtr),
    // End multigrator
    tcoupleCoefficient(tcoupleCoefficientTag, ldObjPtr),
    minimizeCoefficient(minimizeCoefficientTag, ldObjPtr),
    momentumCorrection(momentumCorrectionTag, ldObjPtr),
#if USE_BARRIER
    cycleBarrier(cycleBarrierTag, ldObjPtr),
#endif
    accelMDRescaleFactor(accelMDRescaleFactorTag, ldObjPtr),
    adaptTemperature(adaptTemperatureTag, ldObjPtr), 
    scriptBarrier(scriptBarrierTag, ldObjPtr),
#ifdef MEASURE_NAMD_WITH_PAPI
	papiMeasureBarrier(papiMeasureTag, ldObjPtr),
#endif
	traceBarrier(traceBarrierTag, ldObjPtr)
  { ; }
};

#endif // BROADCASTS_H

