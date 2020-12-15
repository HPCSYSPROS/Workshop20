/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef CONTROLLER_H
#define CONTROLLER_H

#include "converse.h"
#include "Node.h"
#include "common.h"
#include "fstream_namd.h"
#include <string>
#include <map>

class ControllerBroadcasts;
class NamdState;
class SimParameters;
class RequireReduction;
class SubmitReduction;

#ifdef MEM_OPT_VERSION
class CollectionMasterHandler;
#else
class CollectionMaster;
#endif

class Random;
class PressureProfileReduction;

struct ControllerState {
    Tensor langevinPiston_strainRate;
    Tensor berendsenPressure_avg;
    int berendsenPressure_count;
    BigReal smooth2_avg;
};

class Controller : protected ControllerState
{
public:
    Controller(NamdState *s);
    virtual ~Controller(void);
    void run(void);             // spawn thread, etc.
    void awaken(void) { CthAwaken(thread); };
    void resumeAfterTraceBarrier(int);
#ifdef MEASURE_NAMD_WITH_PAPI
	void resumeAfterPapiMeasureBarrier(int step);
#endif

protected:
    friend class ScriptTcl;
    friend class Node;
    friend class CheckpointMsg;
    virtual void algorithm(void);	// subclasses redefine this method

    void integrate(int); // Verlet integrator
    void minimize(); // CG minimizer
      RequireReduction *min_reduction;

    void receivePressure(int step, int minimize = 0);
    void calcPressure(int step, int minimize,
      const Tensor& virial_normal_in, const Tensor& virial_nbond_in, const Tensor& virial_slow_in,
      const Tensor& intVirial_normal, const Tensor& intVirial_nbond, const Tensor& intVirial_slow,
      const Vector& extForce_normal, const Vector& extForce_nbond, const Vector& extForce_slow);

      Tensor pressure_normal;
      Tensor pressure_nbond;
      Tensor pressure_slow;
      Tensor pressure_amd;
      Tensor virial_amd;
      Tensor groupPressure_normal;
      Tensor groupPressure_nbond;
      Tensor groupPressure_slow;
      Tensor controlPressure_normal;
      Tensor controlPressure_nbond;
      Tensor controlPressure_slow;
      int nbondFreq;
      int slowFreq;
      BigReal temp_avg;
      BigReal pressure_avg;
      BigReal groupPressure_avg;
      int avg_count;
      Tensor pressure_tavg;
      Tensor groupPressure_tavg;
      int tavg_count;
    void compareChecksums(int,int=0);
      int computeChecksum;
      int marginViolations;
      int pairlistWarnings;
    void printTiming(int);
    void printMinimizeEnergies(int);
      BigReal min_energy;
      BigReal min_f_dot_f;
      BigReal min_f_dot_v;
      BigReal min_v_dot_v;
      int min_huge_count;
    void printDynamicsEnergies(int);
    void printEnergies(int step, int minimize);
      int numDegFreedom;
      int stepInFullRun;
      BigReal totalEnergy;
      BigReal electEnergy;
      BigReal electEnergySlow;
      BigReal ljEnergy;
      BigReal groLJEnergy;
      BigReal groGaussEnergy;
      BigReal goNativeEnergy;
      BigReal goNonnativeEnergy;
      BigReal goTotalEnergy;
//fepb
      BigReal bondedEnergyDiff_f;
      BigReal electEnergy_f;
      BigReal electEnergySlow_f;
      BigReal ljEnergy_f;
      BigReal ljEnergy_f_left;  // used by WCA repulsive, [s1,s2]
      BigReal exp_dE_ByRT;
      BigReal dE;
      BigReal net_dE;
      BigReal dG;
      int FepNo;
      void printFepMessage(int);
      BigReal fepSum;
//fepe
      BigReal bondedEnergy_ti_1;
      BigReal bondedEnergy_ti_2;
      BigReal electEnergy_ti_1;
      BigReal electEnergySlow_ti_1;
      BigReal ljEnergy_ti_1;
      BigReal electEnergy_ti_2;
      BigReal electEnergySlow_ti_2;
      BigReal ljEnergy_ti_2;
      BigReal net_dEdl_bond_1;
      BigReal net_dEdl_bond_2;
      BigReal net_dEdl_elec_1;
      BigReal net_dEdl_elec_2;
      BigReal net_dEdl_lj_1;
      BigReal net_dEdl_lj_2;
      BigReal cumAlchWork;
      BigReal electEnergyPME_ti_1;
      BigReal electEnergyPME_ti_2;
      int TiNo;
      BigReal recent_dEdl_bond_1;
      BigReal recent_dEdl_bond_2;
      BigReal recent_dEdl_elec_1;
      BigReal recent_dEdl_elec_2;
      BigReal recent_dEdl_lj_1;
      BigReal recent_dEdl_lj_2;
      BigReal recent_alchWork;
      BigReal alchWork;
      int recent_TiNo;
      void printTiMessage(int);

      BigReal drudeBondTemp; // temperature of Drude bonds
      BigReal drudeBondTempAvg;

      BigReal kineticEnergy;
      BigReal kineticEnergyHalfstep;
      BigReal kineticEnergyCentered;
      BigReal temperature;
      // BigReal smooth2_avg;
      BigReal smooth2_avg2;  // avoid internal compiler error
      Tensor pressure;
      Tensor groupPressure;
      int controlNumDegFreedom;
      Tensor controlPressure;
    void enqueueCollections(int);
    void correctMomentum(int step);
    void rescaleVelocities(int);
      BigReal rescaleVelocities_sumTemps;
      int rescaleVelocities_numTemps;
    void reassignVelocities(int);
    void tcoupleVelocities(int);
    void berendsenPressure(int);
      // Tensor berendsenPressure_avg;
      // int berendsenPressure_count;
    void langevinPiston1(int);
    void langevinPiston2(int);
      Tensor langevinPiston_origStrainRate;
      Tensor strainRate_old;  // for langevinPistonBarrier no
      Tensor positionRescaleFactor;  // for langevinPistonBarrier no

    void multigratorPressure(int step, int callNumber);
    BigReal multigratorXi;
    BigReal multigratorXiT;
    Tensor momentumSqrSum;
    void multigratorTemperature(int step, int callNumber);
    std::vector<BigReal> multigratorNu;
    std::vector<BigReal> multigratorNuT;
    std::vector<BigReal> multigratorOmega;
    std::vector<BigReal> multigratorZeta;
    RequireReduction *multigratorReduction;
    BigReal multigatorCalcEnthalpy(BigReal potentialEnergy, int step, int minimize);

    int ldbSteps;
    void rebalanceLoad(int);
      int fflush_count;
    void cycleBarrier(int,int);	
	
	void traceBarrier(int, int);

#ifdef MEASURE_NAMD_WITH_PAPI
	void papiMeasureBarrier(int, int);
#endif

    // void suspend(void) { CthSuspend(); };
    void terminate(void);

    Random *random;
    SimParameters *const simParams;	// for convenience
    NamdState *const state;		// access data in state
    RequireReduction *reduction;
    RequireReduction *amd_reduction;
    SubmitReduction *submit_reduction;

    // data for pressure profile reductions and output
    PressureProfileReduction *ppbonded;
    PressureProfileReduction *ppnonbonded;
    PressureProfileReduction *ppint;
    int pressureProfileSlabs;
    int pressureProfileCount;
    BigReal *pressureProfileAverage;

    CollectionMaster *const collection;
    
    ControllerBroadcasts * broadcast;
    ofstream_namd xstFile;
    void outputExtendedSystem(int step);
    void writeExtendedSystemLabels(ofstream_namd &file);
    void writeExtendedSystemData(int step, ofstream_namd &file);

//fepb
    ofstream_namd fepFile;
    void outputFepEnergy(int step);
    void writeFepEnergyData(int step, ofstream_namd &file);
//fepe
    ofstream_namd tiFile;
    void outputTiEnergy(int step);
    BigReal computeAlchWork(const int step);
    void writeTiEnergyData(int step, ofstream_namd &file);

    // for checkpoint/revert
    int checkpoint_stored;
    Lattice checkpoint_lattice;
    ControllerState checkpoint_state;

    struct checkpoint {
      Lattice lattice;
      ControllerState state;
    };
    std::map<std::string,checkpoint*> checkpoints;
    int checkpoint_task;
    void recvCheckpointReq(const char *key, int task, checkpoint &cp);
    void recvCheckpointAck(checkpoint &cp);

    Lattice origLattice;

//for accelMD
   void rescaleaccelMD (int step, int minimize = 0);
   BigReal accelMDdVAverage;

//JS for adaptive temperature sampling
   void adaptTempInit(int step);
   void adaptTempUpdate(int step, int minimize = 0);
   void adaptTempWriteRestart(int step);
   BigReal *adaptTempPotEnergyAveNum;
   BigReal *adaptTempPotEnergyAveDen;
   BigReal *adaptTempPotEnergyVarNum;
   BigReal *adaptTempPotEnergyAve;
   BigReal *adaptTempPotEnergyVar;
   int     *adaptTempPotEnergySamples;
   BigReal *adaptTempBetaN;
   BigReal adaptTempT;
   BigReal adaptTempDTave;
   BigReal adaptTempDTavenum;
   BigReal adaptTempBetaMin;
   BigReal adaptTempBetaMax;
   int     adaptTempBin;
   int     adaptTempBins;
   BigReal adaptTempDBeta;
   BigReal adaptTempCg;
   BigReal adaptTempDt;
   Bool    adaptTempAutoDt;
   BigReal adaptTempDtMin;
   BigReal adaptTempDtMax;
   ofstream_namd adaptTempRestartFile;
  
private:
    CthThread thread;
    static void threadRun(Controller*);

    double startCTime;
    double startWTime;
    double firstCTime;
    double firstWTime;
    double startBenchTime;

    int computesPartitioned;
};

//Modifications for alchemical fep
static char *FEPTITLE(int X)
{
  static char tmp_string[21];
  sprintf(tmp_string, "FepEnergy: %6d ",X);
  return tmp_string;
}

static char *FEPTITLE2(int X)
{
  static char tmp_string[21];
  sprintf(tmp_string, "FEP:    %7d",X);
  return tmp_string;
}

static char *TITITLE(int X)
{
  static char tmp_string[21];
  sprintf(tmp_string, "TI:     %7d",X);
  return tmp_string;
}
//fepe

#endif // CONTROLLER_H

