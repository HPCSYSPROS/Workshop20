/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef SEQUENCER_H
#define SEQUENCER_H

#include "converse.h"
#include "Priorities.h"
#include "PatchTypes.h"

class HomePatch;
class SimParameters;
class SubmitReduction;
class CollectionMgr;
class ControllerBroadcasts;
class LdbCoordinator;
class Random;

class Sequencer
{
    friend class HomePatch;
public:
    Sequencer(HomePatch *p);
    virtual ~Sequencer(void);
    void run(void);             // spawn thread, etc.
    void awaken(void) {
      CthAwakenPrio(thread, CK_QUEUEING_IFIFO, PRIORITY_SIZE, &priority);
    }
    void suspend(void);

protected:
    virtual void algorithm(void);	// subclasses redefine this method

    void integrate(int); // Verlet integrator
    void minimize(); // CG minimizer
      SubmitReduction *min_reduction;

    void runComputeObjects(int migration = 1, int pairlists = 0, int pressureStep = 0);
    int pairlistsAreValid;
    int pairlistsAge;

    void calcFixVirial(Tensor& fixVirialNormal, Tensor& fixVirialNbond, Tensor& fixVirialSlow,
      Vector& fixForceNormal, Vector& fixForceNbond, Vector& fixForceSlow);

    void submitReductions(int);
    void submitHalfstep(int);
    void submitMinimizeReductions(int, BigReal fmax2);
    void submitCollections(int step, int zeroVel = 0);

    void submitMomentum(int step);
    void correctMomentum(int step, BigReal drifttime);

    void saveForce(const int ftag = Results::normal);
    void addForceToMomentum(BigReal, const int ftag = Results::normal, const int useSaved = 0);
    void addForceToMomentum3(const BigReal timestep1, const int ftag1, const int useSaved1,
        const BigReal timestep2, const int ftag2, const int useSaved2,
        const BigReal timestep3, const int ftag3, const int useSaved3);
    void addVelocityToPosition(BigReal);
    
    void addRotDragToPosition(BigReal);
    void addMovDragToPosition(BigReal);

    void minimizeMoveDownhill(BigReal fmax2);
    void newMinimizeDirection(BigReal);
    void newMinimizePosition(BigReal);
    void quenchVelocities();

    void hardWallDrude(BigReal,int);

    void rattle1(BigReal,int);
    // void rattle2(BigReal,int);

    void maximumMove(BigReal);
    void minimizationQuenchVelocity(void);

    void reloadCharges();

    BigReal adaptTempT;         // adaptive tempering temperature
    void adaptTempUpdate(int); // adaptive tempering temperature update

    void rescaleVelocities(int);
    void rescaleaccelMD(int, int, int); // for accelMD
    int rescaleVelocities_numTemps;
    void reassignVelocities(BigReal,int);
    void reinitVelocities(void);
    void rescaleVelocitiesByFactor(BigReal);
    void tcoupleVelocities(BigReal,int);
    void berendsenPressure(int);
      int berendsenPressure_count;
      int checkpoint_berendsenPressure_count;
    void langevinPiston(int);
      int slowFreq;
    void newtonianVelocities(BigReal, const BigReal, const BigReal, 
                             const BigReal, const int, const int, const int);
    void langevinVelocities(BigReal);
    void langevinVelocitiesBBK1(BigReal);
    void langevinVelocitiesBBK2(BigReal);
    // Multigrator
    void scalePositionsVelocities(const Tensor& posScale, const Tensor& velScale);
    void multigratorPressure(int step, int callNumber);
    void scaleVelocities(const BigReal velScale);
    BigReal calcKineticEnergy();
    void multigratorTemperature(int step, int callNumber);
    SubmitReduction *multigratorReduction;
    int doKineticEnergy;
    int doMomenta;
    // End of Multigrator
    
    void cycleBarrier(int,int);
	void traceBarrier(int);
#ifdef MEASURE_NAMD_WITH_PAPI
	void papiMeasureBarrier(int);
#endif
    void terminate(void);

    Random *random;
    SimParameters *const simParams;	// for convenience
    HomePatch *const patch;		// access methods in patch
    SubmitReduction *reduction;
    SubmitReduction *pressureProfileReduction;

    CollectionMgr *const collection;
    ControllerBroadcasts * broadcast;

    int ldbSteps;
    void rebalanceLoad(int timestep);


private:
    CthThread thread;
    unsigned int priority;
    static void threadRun(Sequencer*);

    LdbCoordinator *ldbCoordinator;
};

#endif
