/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*****************************************************************************
 * $Source: /home/cvs/namd/cvsroot/namd2/src/Controller.C,v $
 * $Author: jim $
 * $Date: 2016/10/14 16:12:12 $
 * $Revision: 1.1324 $
 *****************************************************************************/

#include "InfoStream.h"
#include "memusage.h"
#include "Node.h"
#include "Molecule.h"
#include "SimParameters.h"
#include "Controller.h"
#include "ReductionMgr.h"
#include "CollectionMaster.h"
#include "Output.h"
#include "strlib.h"
#include "BroadcastObject.h"
#include "NamdState.h"
#include "ScriptTcl.h"
#include "Broadcasts.h"
#include "LdbCoordinator.h"
#include "Thread.h"
#include <math.h>
#include <signal.h>
#include "NamdOneTools.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "Random.h"
#include "imd.h"
#include "IMDOutput.h"
#include "BackEnd.h"
#include <fstream>
#include <iomanip>
#include <errno.h>
#include "qd.h"

#include "ComputeNonbondedMICKernel.h"

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
extern "C" void CApplicationDepositNode0Data(char *);
#endif

#ifndef cbrt
  // cbrt() not in math.h on goneril
  #define cbrt(x)  pow(x,(double)(1.0/3.0))
#endif

//#define DEBUG_PRESSURE
#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

#define XXXBIGREAL 1.0e32

class PressureProfileReduction {
private:
  RequireReduction *reduction;
  const int nslabs;
  const int freq;
  int nelements;
  int count;
  char *name;

  BigReal *average;

public:
  PressureProfileReduction(int rtag, int numslabs, int numpartitions,
      const char *myname, int outputfreq)
  : nslabs(numslabs), freq(outputfreq) {
    name = strdup(myname);
    nelements = 3*nslabs * numpartitions;
    reduction = ReductionMgr::Object()->willRequire(rtag,nelements);

    average = new BigReal[nelements];
    count = 0;
  }
  ~PressureProfileReduction() {
    delete [] average;
    delete reduction;
    free(name);
  }
  // 
  void getData(int firsttimestep, int step, const Lattice &lattice, 
      BigReal *total) {
    reduction->require();

    int i;
    double inv_volume = 1.0 / lattice.volume();
    // accumulate the pressure profile computed for this step into the average.
    int arraysize = 3*nslabs;
    for (i=0; i<nelements; i++) {
      BigReal val = reduction->item(i) * inv_volume;
      average[i] += val;
      total[i % arraysize] += val;
    }
    count++;
    if (!(step % freq)) {
      // convert NAMD internal virial to pressure in units of bar 
      BigReal scalefac = PRESSUREFACTOR * nslabs;
 
      iout << "PPROFILE" << name << ": " << step << " ";
      if (step == firsttimestep) {
        // output pressure profile for this step
        for (i=0; i<nelements; i++) 
          iout << reduction->item(i)*scalefac*inv_volume << " ";
      } else {
        // output pressure profile averaged over the last Freq steps.
        scalefac /= count;
        for (i=0; i<nelements; i++) 
          iout << average[i]*scalefac << " ";
      }
      iout << "\n" << endi; 
      // Clear the average for the next block
      memset(average, 0, nelements*sizeof(BigReal));
      count = 0;
    }
  }
};


Controller::Controller(NamdState *s) :
	computeChecksum(0), marginViolations(0), pairlistWarnings(0),
	simParams(Node::Object()->simParameters),
	state(s),
	collection(CollectionMaster::Object()),
        startCTime(0),
        firstCTime(CmiTimer()),
        startWTime(0),
        firstWTime(CmiWallTimer()),
        startBenchTime(0),
        stepInFullRun(0),
	ldbSteps(0),
	fflush_count(3)
{
    broadcast = new ControllerBroadcasts;
    min_reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_MINIMIZER,1);
    // for accelMD
    if (simParams->accelMDOn) {
       // REDUCTIONS_BASIC wil contain all potential energies and arrive first
       amd_reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_BASIC);
       // contents of amd_reduction be submitted to REDUCTIONS_AMD
       submit_reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_AMD);
       // REDUCTIONS_AMD will contain Sequencer contributions and arrive second
       reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_AMD);
    } else {
       amd_reduction = NULL;
       submit_reduction = NULL;
       reduction = ReductionMgr::Object()->willRequire(REDUCTIONS_BASIC);
    }
    // pressure profile reductions
    pressureProfileSlabs = 0;
    pressureProfileCount = 0;
    ppbonded = ppnonbonded = ppint = NULL;
    if (simParams->pressureProfileOn || simParams->pressureProfileEwaldOn) {
      int ntypes = simParams->pressureProfileAtomTypes;
      int nslabs = pressureProfileSlabs = simParams->pressureProfileSlabs;
      // number of partitions for pairwise interactions
      int npairs = (ntypes * (ntypes+1))/2;
      pressureProfileAverage = new BigReal[3*nslabs];
      memset(pressureProfileAverage, 0, 3*nslabs*sizeof(BigReal));
      int freq = simParams->pressureProfileFreq;
      if (simParams->pressureProfileOn) {
        ppbonded = new PressureProfileReduction(REDUCTIONS_PPROF_BONDED, 
            nslabs, npairs, "BONDED", freq);
        ppnonbonded = new PressureProfileReduction(REDUCTIONS_PPROF_NONBONDED, 
            nslabs, npairs, "NONBONDED", freq);
        // internal partitions by atom type, but not in a pairwise fashion
        ppint = new PressureProfileReduction(REDUCTIONS_PPROF_INTERNAL, 
            nslabs, ntypes, "INTERNAL", freq);
      } else {
        // just doing Ewald, so only need nonbonded
        ppnonbonded = new PressureProfileReduction(REDUCTIONS_PPROF_NONBONDED, 
            nslabs, npairs, "NONBONDED", freq);
      }
    }
    random = new Random(simParams->randomSeed);
    random->split(0,PatchMap::Object()->numPatches()+1);

    rescaleVelocities_sumTemps = 0;  rescaleVelocities_numTemps = 0;
    berendsenPressure_avg = 0; berendsenPressure_count = 0;
    // strainRate tensor is symmetric to avoid rotation
    langevinPiston_strainRate =
	Tensor::symmetric(simParams->strainRate,simParams->strainRate2);
    if ( ! simParams->useFlexibleCell ) {
      BigReal avg = trace(langevinPiston_strainRate) / 3.;
      langevinPiston_strainRate = Tensor::identity(avg);
    } else if ( simParams->useConstantRatio ) {
#define AVGXY(T) T.xy = T.yx = 0; T.xx = T.yy = 0.5 * ( T.xx + T.yy );\
		 T.xz = T.zx = T.yz = T.zy = 0.5 * ( T.xz + T.yz );
      AVGXY(langevinPiston_strainRate);
#undef AVGXY
    }
    langevinPiston_origStrainRate = langevinPiston_strainRate;
    if (simParams->multigratorOn) {
      multigratorXi = 0.0;
      int n = simParams->multigratorNoseHooverChainLength;
      BigReal tau = simParams->multigratorTemperatureRelaxationTime;
      Node *node = Node::Object();
      Molecule *molecule = node->molecule;
      BigReal Nf = molecule->num_deg_freedom();
      BigReal kT0 = BOLTZMANN * simParams->multigratorTemperatureTarget;
      multigratorNu.resize(n);
      multigratorNuT.resize(n);
      multigratorZeta.resize(n);
      multigratorOmega.resize(n);
      for (int i=0;i < n;i++) {
        multigratorNu[i] = 0.0;
        multigratorZeta[i] = 0.0;
        if (i == 0) {
          multigratorOmega[i] = Nf*kT0*tau*tau;
        } else {
          multigratorOmega[i] = kT0*tau*tau;
        }
      }
      multigratorReduction = ReductionMgr::Object()->willRequire(REDUCTIONS_MULTIGRATOR,MULTIGRATOR_REDUCTION_MAX_RESERVED);
    } else {
      multigratorReduction = NULL;
    }
    origLattice = state->lattice;
    smooth2_avg = XXXBIGREAL;
    temp_avg = 0;
    pressure_avg = 0;
    groupPressure_avg = 0;
    avg_count = 0;
    pressure_tavg = 0;
    groupPressure_tavg = 0;
    tavg_count = 0;
    checkpoint_stored = 0;
    drudeBondTemp = 0;
    drudeBondTempAvg = 0;
    cumAlchWork = 0;
}

Controller::~Controller(void)
{
    delete broadcast;
    delete reduction;
    delete min_reduction;
    delete amd_reduction;
    delete submit_reduction;
    delete ppbonded;
    delete ppnonbonded;
    delete ppint;
    delete [] pressureProfileAverage;
    delete random;
    if (multigratorReduction) delete multigratorReduction;
}

void Controller::threadRun(Controller* arg)
{
    arg->algorithm();
}

void Controller::run(void)
{
    // create a Thread and invoke it
    DebugM(4, "Starting thread in controller on this=" << this << "\n");
    thread = CthCreate((CthVoidFn)&(threadRun),(void*)(this),CTRL_STK_SZ);
    CthSetStrategyDefault(thread);
#if CMK_BLUEGENE_CHARM
    BgAttach(thread);
#endif
    awaken();
}


void Controller::algorithm(void)
{
  int scriptTask;
  int scriptSeq = 0;
  BackEnd::awaken();
  while ( (scriptTask = broadcast->scriptBarrier.get(scriptSeq++)) != SCRIPT_END) {
    switch ( scriptTask ) {
      case SCRIPT_OUTPUT:
        enqueueCollections(FILE_OUTPUT);
        outputExtendedSystem(FILE_OUTPUT);
        break;
      case SCRIPT_FORCEOUTPUT:
        enqueueCollections(FORCE_OUTPUT);
        break;
      case SCRIPT_MEASURE:
        enqueueCollections(EVAL_MEASURE);
        break;
      case SCRIPT_REINITVELS:
        iout << "REINITIALIZING VELOCITIES AT STEP " << simParams->firstTimestep
          << " TO " << simParams->initialTemp << " KELVIN.\n" << endi;
        break;
      case SCRIPT_RESCALEVELS:
        iout << "RESCALING VELOCITIES AT STEP " << simParams->firstTimestep
          << " BY " << simParams->scriptArg1 << "\n" << endi;
        break;
      case SCRIPT_CHECKPOINT:
        iout << "CHECKPOINTING AT STEP " << simParams->firstTimestep
          << "\n" << endi;
        checkpoint_stored = 1;
        checkpoint_lattice = state->lattice;
        checkpoint_state = *(ControllerState*)this;
        break;
      case SCRIPT_REVERT:
        iout << "REVERTING AT STEP " << simParams->firstTimestep
          << "\n" << endi;
        if ( ! checkpoint_stored )
          NAMD_die("Unable to revert, checkpoint was never called!");
        state->lattice = checkpoint_lattice;
        *(ControllerState*)this = checkpoint_state;
        break;
      case SCRIPT_CHECKPOINT_STORE:
        iout << "STORING CHECKPOINT AT STEP " << simParams->firstTimestep
          << " TO KEY " << simParams->scriptStringArg1;
        if ( CmiNumPartitions() > 1 ) iout << " ON REPLICA " << simParams->scriptIntArg1;
        iout << "\n" << endi;
        if ( simParams->scriptIntArg1 == CmiMyPartition() ) {
          if ( ! checkpoints.count(simParams->scriptStringArg1) ) {
            checkpoints[simParams->scriptStringArg1] = new checkpoint;
          }
          checkpoint &cp = *checkpoints[simParams->scriptStringArg1];
          cp.lattice = state->lattice;
          cp.state = *(ControllerState*)this;
        } else {
          Node::Object()->sendCheckpointReq(simParams->scriptIntArg1,simParams->scriptStringArg1,
                                            scriptTask, state->lattice, *(ControllerState*)this);
        }
        break;
      case SCRIPT_CHECKPOINT_LOAD:
        iout << "LOADING CHECKPOINT AT STEP " << simParams->firstTimestep
          << " FROM KEY " << simParams->scriptStringArg1;
        if ( CmiNumPartitions() > 1 ) iout << " ON REPLICA " << simParams->scriptIntArg1;
        iout << "\n" << endi;
        if ( simParams->scriptIntArg1 == CmiMyPartition() ) {
          if ( ! checkpoints.count(simParams->scriptStringArg1) ) {
            NAMD_die("Unable to load checkpoint, requested key was never stored.");
          }
          checkpoint &cp = *checkpoints[simParams->scriptStringArg1];
          state->lattice = cp.lattice;
          *(ControllerState*)this = cp.state;
        } else {
          Node::Object()->sendCheckpointReq(simParams->scriptIntArg1,simParams->scriptStringArg1,
                                            scriptTask, state->lattice, *(ControllerState*)this);
        }
        break;
      case SCRIPT_CHECKPOINT_SWAP:
        iout << "SWAPPING CHECKPOINT AT STEP " << simParams->firstTimestep
          << " FROM KEY " << simParams->scriptStringArg1;
        if ( CmiNumPartitions() > 1 ) iout << " ON REPLICA " << simParams->scriptIntArg1;
        iout << "\n" << endi;
        if ( simParams->scriptIntArg1 == CmiMyPartition() ) {
          if ( ! checkpoints.count(simParams->scriptStringArg1) ) {
            NAMD_die("Unable to swap checkpoint, requested key was never stored.");
          }
          checkpoint &cp = *checkpoints[simParams->scriptStringArg1];
          std::swap(state->lattice,cp.lattice);
          std::swap(*(ControllerState*)this,cp.state);
        } else {
          Node::Object()->sendCheckpointReq(simParams->scriptIntArg1,simParams->scriptStringArg1,
                                            scriptTask, state->lattice, *(ControllerState*)this);
        }
        break;
      case SCRIPT_CHECKPOINT_FREE:
        iout << "FREEING CHECKPOINT AT STEP " << simParams->firstTimestep
          << " FROM KEY " << simParams->scriptStringArg1;
        if ( CmiNumPartitions() > 1 ) iout << " ON REPLICA " << simParams->scriptIntArg1;
        iout << "\n" << endi;
        if ( simParams->scriptIntArg1 == CmiMyPartition() ) {
          if ( ! checkpoints.count(simParams->scriptStringArg1) ) {
            NAMD_die("Unable to free checkpoint, requested key was never stored.");
          }
          delete checkpoints[simParams->scriptStringArg1];
          checkpoints.erase(simParams->scriptStringArg1);
        } else {
          Node::Object()->sendCheckpointReq(simParams->scriptIntArg1,simParams->scriptStringArg1,
                                            scriptTask, state->lattice, *(ControllerState*)this);
        }
        break;
      case SCRIPT_ATOMSENDRECV:
      case SCRIPT_ATOMSEND:
      case SCRIPT_ATOMRECV:
        break;
      case SCRIPT_MINIMIZE:
        minimize();
        break;
      case SCRIPT_RUN:
      case SCRIPT_CONTINUE:
        integrate(scriptTask);
        break;
      default:
        NAMD_bug("Unknown task in Controller::algorithm");
    }
    BackEnd::awaken();
  }
  enqueueCollections(END_OF_RUN);
  outputExtendedSystem(END_OF_RUN);
  terminate();
}


//
// XXX static and global variables are unsafe for shared memory builds.
// The use of global and static vars should be eliminated.
//
extern int eventEndOfTimeStep;

// Handle SIGINT so that restart files get written completely.
static int gotsigint = 0;
static void my_sigint_handler(int sig) {
  if (sig == SIGINT) gotsigint = 1;
}
extern "C" {
  typedef void (*namd_sighandler_t)(int);
}

void Controller::integrate(int scriptTask) {
    char traceNote[24];
  
    int step = simParams->firstTimestep;

    const int numberOfSteps = simParams->N;
    const int stepsPerCycle = simParams->stepsPerCycle;

    const int zeroMomentum = simParams->zeroMomentum;

    nbondFreq = simParams->nonbondedFrequency;
    const int dofull = ( simParams->fullElectFrequency ? 1 : 0 );
    if (dofull)
      slowFreq = simParams->fullElectFrequency;
    else
      slowFreq = simParams->nonbondedFrequency;
    if ( step >= numberOfSteps ) slowFreq = nbondFreq = 1;

  if ( scriptTask == SCRIPT_RUN ) {

    reassignVelocities(step);  // only for full-step velecities
    rescaleaccelMD(step);
    adaptTempUpdate(step); // Init adaptive tempering;

    receivePressure(step);
    if ( zeroMomentum && dofull && ! (step % slowFreq) )
						correctMomentum(step);
    printFepMessage(step);
    printTiMessage(step);
    printDynamicsEnergies(step);
    outputFepEnergy(step);
    outputTiEnergy(step);
    if(traceIsOn()){
        traceUserEvent(eventEndOfTimeStep);
        sprintf(traceNote, "s:%d", step);
        traceUserSuppliedNote(traceNote);
    }
    outputExtendedSystem(step);
    rebalanceLoad(step);

  }

    // Handling SIGINT doesn't seem to be working on Lemieux, and it
    // sometimes causes the net-xxx versions of NAMD to segfault on exit, 
    // so disable it for now.
    // namd_sighandler_t oldhandler = signal(SIGINT, 
    //  (namd_sighandler_t)my_sigint_handler);
    for ( ++step ; step <= numberOfSteps; ++step )
    {

        adaptTempUpdate(step);
        rescaleVelocities(step);
	tcoupleVelocities(step);
	berendsenPressure(step);
	langevinPiston1(step);
        rescaleaccelMD(step);
	enqueueCollections(step);  // after lattice scaling!
	receivePressure(step);
        if ( zeroMomentum && dofull && ! (step % slowFreq) )
						correctMomentum(step);
	langevinPiston2(step);
        reassignVelocities(step);

      multigratorTemperature(step, 1);
      multigratorPressure(step, 1);
      multigratorPressure(step, 2);
      multigratorTemperature(step, 2);

        printDynamicsEnergies(step);
        outputFepEnergy(step);
        outputTiEnergy(step);
        if(traceIsOn()){
            traceUserEvent(eventEndOfTimeStep);
            sprintf(traceNote, "s:%d", step);
            traceUserSuppliedNote(traceNote);
        }
  // if (gotsigint) {
  //   iout << iINFO << "Received SIGINT; shutting down.\n" << endi;
  //   NAMD_quit();
  // }
        outputExtendedSystem(step);
#if CYCLE_BARRIER
        cycleBarrier(!((step+1) % stepsPerCycle),step);
#elif  PME_BARRIER
        cycleBarrier(dofull && !(step%slowFreq),step);
#elif  STEP_BARRIER
        cycleBarrier(1, step);
#endif
		
	if(Node::Object()->specialTracing || simParams->statsOn){		
		 int bstep = simParams->traceStartStep;
		 int estep = bstep + simParams->numTraceSteps;		
		 if(step == bstep){
			 traceBarrier(1, step);
		 }
		 if(step == estep){			
			 traceBarrier(0, step);
		 }
	 }

#ifdef MEASURE_NAMD_WITH_PAPI
	if(simParams->papiMeasure) {
		int bstep = simParams->papiMeasureStartStep;
		int estep = bstep + simParams->numPapiMeasureSteps;
		if(step == bstep) {
			papiMeasureBarrier(1, step);
		}
		if(step == estep) {
			papiMeasureBarrier(0, step);
		}
	}
#endif
	 
        rebalanceLoad(step);

#if  PME_BARRIER
        cycleBarrier(dofull && !((step+1)%slowFreq),step);   // step before PME
#endif
    }
    // signal(SIGINT, oldhandler);
}


#define CALCULATE \
  printMinimizeEnergies(step); \
  outputExtendedSystem(step); \
  rebalanceLoad(step); \
  if ( step == numberOfSteps ) return; \
  else ++step;

#define MOVETO(X) \
  if ( step == numberOfSteps ) { \
    if ( minVerbose ) { iout << "LINE MINIMIZER: RETURNING TO " << mid.x << " FROM " << last.x << "\n" << endi; } \
    if ( newDir || (mid.x-last.x) ) { \
      broadcast->minimizeCoefficient.publish(minSeq++,mid.x-last.x); \
    } else { \
      broadcast->minimizeCoefficient.publish(minSeq++,0.); \
      broadcast->minimizeCoefficient.publish(minSeq++,0.); \
      min_reduction->require(); \
      broadcast->minimizeCoefficient.publish(minSeq++,0.); \
    } \
    enqueueCollections(step); \
    CALCULATE \
  } else if ( (X)-last.x ) { \
    broadcast->minimizeCoefficient.publish(minSeq++,(X)-last.x); \
    newDir = 0; \
    last.x = (X); \
    enqueueCollections(step); \
    CALCULATE \
    last.u = min_energy; \
    last.dudx = -1. * min_f_dot_v; \
    last.noGradient = min_huge_count; \
    if ( minVerbose ) { \
      iout << "LINE MINIMIZER: POSITION " << last.x << " ENERGY " << last.u << " GRADIENT " << last.dudx; \
      if ( last.noGradient ) iout << " HUGECOUNT " << last.noGradient; \
      iout << "\n" << endi; \
    } \
  }

struct minpoint {
  BigReal x, u, dudx, noGradient;
  minpoint() : x(0), u(0), dudx(0), noGradient(1) { ; }
};

void Controller::minimize() {
  // iout << "Controller::minimize() called.\n" << endi;

  const int minVerbose = simParams->minVerbose;
  const int numberOfSteps = simParams->N;
  int step = simParams->firstTimestep;
  slowFreq = nbondFreq = 1;
  BigReal linegoal = simParams->minLineGoal;  // 1.0e-3
  const BigReal goldenRatio = 0.5 * ( sqrt(5.0) - 1.0 );

  CALCULATE

  int minSeq = 0;

  // just move downhill until initial bad contacts go away or 100 steps
  int old_min_huge_count = 2000000000;
  int steps_at_same_min_huge_count = 0;
  for ( int i_slow = 0; min_huge_count && i_slow < 100; ++i_slow ) {
    if ( min_huge_count >= old_min_huge_count ) {
      if ( ++steps_at_same_min_huge_count > 10 ) break;
    } else {
      old_min_huge_count = min_huge_count;
      steps_at_same_min_huge_count = 0;
    }
    iout << "MINIMIZER SLOWLY MOVING " << min_huge_count << " ATOMS WITH BAD CONTACTS DOWNHILL\n" << endi;
    broadcast->minimizeCoefficient.publish(minSeq++,1.);
    enqueueCollections(step);
    CALCULATE
  }
  if ( min_huge_count ) {
    iout << "MINIMIZER GIVING UP ON " << min_huge_count << " ATOMS WITH BAD CONTACTS\n" << endi;
  }
  iout << "MINIMIZER STARTING CONJUGATE GRADIENT ALGORITHM\n" << endi;

  int atStart = 2;
  int errorFactor = 10;
  BigReal old_f_dot_f = min_f_dot_f;
  broadcast->minimizeCoefficient.publish(minSeq++,0.);
  broadcast->minimizeCoefficient.publish(minSeq++,0.); // v = f
  int newDir = 1;
  min_f_dot_v = min_f_dot_f;
  min_v_dot_v = min_f_dot_f;
  while ( 1 ) {
    // line minimization
    // bracket minimum on line
    minpoint lo,hi,mid,last;
    BigReal x = 0;
    lo.x = x;
    lo.u = min_energy;
    lo.dudx = -1. * min_f_dot_v;
    lo.noGradient = min_huge_count;
    mid = lo;
    last = mid;
    if ( minVerbose ) {
      iout << "LINE MINIMIZER: POSITION " << last.x << " ENERGY " << last.u << " GRADIENT " << last.dudx;
      if ( last.noGradient ) iout << " HUGECOUNT " << last.noGradient;
      iout << "\n" << endi;
    }
    BigReal tol = fabs( linegoal * min_f_dot_v );
    iout << "LINE MINIMIZER REDUCING GRADIENT FROM " <<
            fabs(min_f_dot_v) << " TO " << tol << "\n" << endi;
    int start_with_huge = last.noGradient;
    min_reduction->require();
    BigReal maxstep = 0.1 / sqrt(min_reduction->item(0));
    x = maxstep; MOVETO(x);
    // bracket minimum on line
    while ( last.u < mid.u ) {
      if ( last.dudx < mid.dudx * (5.5 - x/maxstep)/5. ) {
        x = 6 * maxstep; break;
      }
      lo = mid; mid = last;
      x += maxstep;
      if ( x > 5.5 * maxstep ) break;
      MOVETO(x)
    }
    if ( x > 5.5 * maxstep ) {
      iout << "MINIMIZER RESTARTING CONJUGATE GRADIENT ALGORITHM DUE TO POOR PROGRESS\n" << endi;
      broadcast->minimizeCoefficient.publish(minSeq++,0.);
      broadcast->minimizeCoefficient.publish(minSeq++,0.); // v = f
      newDir = 1;
      old_f_dot_f = min_f_dot_f;
      min_f_dot_v = min_f_dot_f;
      min_v_dot_v = min_f_dot_f;
      continue;
    }
    hi = last;
#define PRINT_BRACKET \
    iout << "LINE MINIMIZER BRACKET: DX " \
         << (mid.x-lo.x) << " " << (hi.x-mid.x) << \
        " DU " << (mid.u-lo.u) << " " << (hi.u-mid.u) << " DUDX " << \
        lo.dudx << " " << mid.dudx << " " << hi.dudx << " \n" << endi;
    PRINT_BRACKET
    // converge on minimum on line
    int itcnt;
    for ( itcnt = 10; itcnt > 0; --itcnt ) {
      int progress = 1;
      // select new position
      if ( mid.noGradient ) {
       if ( ( mid.x - lo.x ) > ( hi.x - mid.x ) ) {  // subdivide left side
	x = (1.0 - goldenRatio) * lo.x + goldenRatio * mid.x;
	MOVETO(x)
	if ( last.u <= mid.u ) {
	  hi = mid; mid = last;
	} else {
	  lo = last;
	}
       } else {  // subdivide right side
	x = (1.0 - goldenRatio) * hi.x + goldenRatio * mid.x;
	MOVETO(x)
	if ( last.u <= mid.u ) {
	  lo = mid; mid = last;
	} else {
	  hi = last;
	}
       }
      } else {
       if ( mid.dudx > 0. ) {  // subdivide left side
        BigReal altxhi = 0.1 * lo.x + 0.9 * mid.x;
        BigReal altxlo = 0.9 * lo.x + 0.1 * mid.x;
        x = mid.dudx*(mid.x*mid.x-lo.x*lo.x) + 2*mid.x*(lo.u-mid.u);
        BigReal xdiv = 2*(mid.dudx*(mid.x-lo.x)+(lo.u-mid.u));
        if ( xdiv ) x /= xdiv; else x = last.x;
        if ( x > altxhi ) x = altxhi;
        if ( x < altxlo ) x = altxlo;
        if ( x-last.x == 0 ) break;
        MOVETO(x)
        if ( last.u <= mid.u ) {
	  hi = mid; mid = last;
	} else {
          if ( lo.dudx < 0. && last.dudx > 0. ) progress = 0;
	  lo = last;
	}
       } else {  // subdivide right side
        BigReal altxlo = 0.1 * hi.x + 0.9 * mid.x;
        BigReal altxhi = 0.9 * hi.x + 0.1 * mid.x;
        x = mid.dudx*(mid.x*mid.x-hi.x*hi.x) + 2*mid.x*(hi.u-mid.u);
        BigReal xdiv = 2*(mid.dudx*(mid.x-hi.x)+(hi.u-mid.u));
        if ( xdiv ) x /= xdiv; else x = last.x;
        if ( x < altxlo ) x = altxlo;
        if ( x > altxhi ) x = altxhi;
        if ( x-last.x == 0 ) break;
        MOVETO(x)
        if ( last.u <= mid.u ) {
	  lo = mid; mid = last;
	} else {
          if ( hi.dudx > 0. && last.dudx < 0. ) progress = 0;
	  hi = last;
	}
       }
      }
      PRINT_BRACKET
      if ( ! progress ) {
        MOVETO(mid.x);
        break;
      }
      if ( fabs(last.dudx) < tol ) break;
      if ( lo.x != mid.x && (lo.u-mid.u) < tol * (mid.x-lo.x) ) break;
      if ( mid.x != hi.x && (hi.u-mid.u) < tol * (hi.x-mid.x) ) break;
    }
    // new direction
    broadcast->minimizeCoefficient.publish(minSeq++,0.);
    BigReal c = min_f_dot_f / old_f_dot_f;
    c = ( c > 1.5 ? 1.5 : c );
    if ( atStart ) { c = 0; --atStart; }
    if ( c*c*min_v_dot_v > errorFactor*min_f_dot_f ) {
      c = 0;
      if ( errorFactor < 100 ) errorFactor += 10;
    }
    if ( c == 0 ) {
      iout << "MINIMIZER RESTARTING CONJUGATE GRADIENT ALGORITHM\n" << endi;
    }
    broadcast->minimizeCoefficient.publish(minSeq++,c); // v = c*v+f
    newDir = 1;
    old_f_dot_f = min_f_dot_f;
    min_f_dot_v = c * min_f_dot_v + min_f_dot_f;
    min_v_dot_v = c*c*min_v_dot_v + 2*c*min_f_dot_v + min_f_dot_f;
  }

}

#undef MOVETO
#undef CALCULATE

// NOTE: Only isotropic case implemented here!
void Controller::multigratorPressure(int step, int callNumber) {
  if (simParams->multigratorOn && !(step % simParams->multigratorPressureFreq)) {
    BigReal P0 = simParams->multigratorPressureTarget;
    BigReal kT0 = BOLTZMANN * simParams->multigratorTemperatureTarget;
    BigReal NG = controlNumDegFreedom;
    BigReal NGfac = 1.0/sqrt(3.0*NG + 1.0);
    BigReal tau = simParams->multigratorPressureRelaxationTime;
    BigReal s = 0.5*simParams->multigratorPressureFreq*simParams->dt/tau;
    {
      // Compute new scaling factors and send them to Sequencer
      BigReal V = state->lattice.volume();
      BigReal Pinst = trace(controlPressure)/3.0;
      BigReal PGsum = trace(momentumSqrSum);
      //
      multigratorXiT = multigratorXi + 0.5*s*NGfac/kT0*( 3.0*V*(Pinst - P0) + PGsum/NG );
      BigReal scale = exp(s*NGfac*multigratorXiT);
      BigReal velScale = exp(-s*NGfac*(1.0 + 1.0/NG)*multigratorXiT);
      // fprintf(stderr, "%d | T %lf P %lf V %1.3lf\n", step, temperature, Pinst*PRESSUREFACTOR, V);
      Tensor scaleTensor = Tensor::identity(scale);
      Tensor volScaleTensor = Tensor::identity(scale);
      Tensor velScaleTensor = Tensor::identity(velScale);
      state->lattice.rescale(volScaleTensor);
      if (callNumber == 1) {
        broadcast->positionRescaleFactor.publish(step,scaleTensor);
        broadcast->velocityRescaleTensor.publish(step,velScaleTensor);
      } else {
        broadcast->positionRescaleFactor2.publish(step,scaleTensor);
        broadcast->velocityRescaleTensor2.publish(step,velScaleTensor);      
      }
    }

    {
      // Wait here for Sequencer to finish scaling and force computation
      reduction->require();
      Tensor virial_normal;
      Tensor virial_nbond;
      Tensor virial_slow;
      Tensor intVirial_normal;
      Tensor intVirial_nbond;
      Tensor intVirial_slow;
      Vector extForce_normal;
      Vector extForce_nbond;
      Vector extForce_slow;
      GET_TENSOR(momentumSqrSum, reduction, REDUCTION_MOMENTUM_SQUARED);
      GET_TENSOR(virial_normal, reduction, REDUCTION_VIRIAL_NORMAL);
      GET_TENSOR(virial_nbond, reduction, REDUCTION_VIRIAL_NBOND);
      GET_TENSOR(virial_slow, reduction, REDUCTION_VIRIAL_SLOW);
      GET_TENSOR(intVirial_normal, reduction, REDUCTION_INT_VIRIAL_NORMAL);
      GET_TENSOR(intVirial_nbond, reduction, REDUCTION_INT_VIRIAL_NBOND);
      GET_TENSOR(intVirial_slow, reduction, REDUCTION_INT_VIRIAL_SLOW);
      GET_VECTOR(extForce_normal, reduction, REDUCTION_EXT_FORCE_NORMAL);
      GET_VECTOR(extForce_nbond, reduction, REDUCTION_EXT_FORCE_NBOND);
      GET_VECTOR(extForce_slow, reduction, REDUCTION_EXT_FORCE_SLOW);
      calcPressure(step, 0, virial_normal, virial_nbond, virial_slow,
        intVirial_normal, intVirial_nbond, intVirial_slow,
        extForce_normal, extForce_nbond, extForce_slow);
      if (callNumber == 2) {
        // Update temperature for the Temperature Cycle that is coming next
        BigReal kineticEnergy = reduction->item(REDUCTION_CENTERED_KINETIC_ENERGY);
        temperature = 2.0 * kineticEnergy / ( numDegFreedom * BOLTZMANN );
      }
    }

    {
      // Update pressure integrator
      BigReal V = state->lattice.volume();
      BigReal Pinst = trace(controlPressure)/3.0;
      BigReal PGsum = trace(momentumSqrSum);
      //
      multigratorXi = multigratorXiT + 0.5*s*NGfac/kT0*( 3.0*V*(Pinst - P0) + PGsum/NG );
    }

  }
}

void Controller::multigratorTemperature(int step, int callNumber) {
  if (simParams->multigratorOn && !(step % simParams->multigratorTemperatureFreq)) {
    BigReal tau = simParams->multigratorTemperatureRelaxationTime;
    BigReal t = 0.5*simParams->multigratorTemperatureFreq*simParams->dt;
    BigReal kT0 = BOLTZMANN * simParams->multigratorTemperatureTarget;
    BigReal Nf = numDegFreedom;
    int n = simParams->multigratorNoseHooverChainLength;
    BigReal T1, T2, v;
    {
      T1 = temperature;
      BigReal kTinst = BOLTZMANN * temperature;
      for (int i=n-1;i >= 0;i--) {
        if (i == 0) {
          BigReal NuOmega = (n > 1) ? multigratorNuT[1]/multigratorOmega[1] : 0.0;
          multigratorNuT[0] = exp(-0.5*t*NuOmega)*multigratorNu[0] + 0.5*Nf*t*exp(-0.25*t*NuOmega)*(kTinst - kT0);
        } else if (i == n-1) {
          multigratorNuT[n-1] = multigratorNu[n-1] + 0.5*t*(pow(multigratorNu[n-2],2.0)/multigratorOmega[n-2] - kT0);
        } else {
          BigReal NuOmega = multigratorNuT[i+1]/multigratorOmega[i+1];
          multigratorNuT[i] = exp(-0.5*t*NuOmega)*multigratorNu[i] + 
          0.5*t*exp(-0.25*t*NuOmega)*(pow(multigratorNu[i-1],2.0)/multigratorOmega[i-1] - kT0);
        }
      }
      BigReal velScale = exp(-t*multigratorNuT[0]/multigratorOmega[0]);
      v = velScale;
      if (callNumber == 1)
        broadcast->velocityRescaleFactor.publish(step,velScale);
      else
        broadcast->velocityRescaleFactor2.publish(step,velScale);
    }

    {
      // Wait here for Sequencer to finish scaling and re-calculating kinetic energy
      multigratorReduction->require();
      BigReal kineticEnergy = multigratorReduction->item(MULTIGRATOR_REDUCTION_KINETIC_ENERGY);
      temperature = 2.0 * kineticEnergy / ( numDegFreedom * BOLTZMANN );
      T2 = temperature;
      if (callNumber == 1 && !(step % simParams->multigratorPressureFreq)) {
        // If this is pressure cycle, receive new momentum product
        GET_TENSOR(momentumSqrSum, multigratorReduction, MULTIGRATOR_REDUCTION_MOMENTUM_SQUARED);
      }
    }

    // fprintf(stderr, "%d | T %lf scale %lf T' %lf\n", step, T1, v, T2);

    {
      BigReal kTinst = BOLTZMANN * temperature;
      for (int i=0;i < n;i++) {
        if (i == 0) {
          BigReal NuOmega = (n > 1) ? multigratorNuT[1]/multigratorOmega[1] : 0.0;
          multigratorNu[0] = exp(-0.5*t*NuOmega)*multigratorNuT[0] + 0.5*Nf*t*exp(-0.25*t*NuOmega)*(kTinst - kT0);
        } else if (i == n-1) {
          multigratorNu[n-1] = multigratorNuT[n-1] + 0.5*t*(pow(multigratorNu[n-2],2.0)/multigratorOmega[n-2] - kT0);
        } else {
          BigReal NuOmega = multigratorNuT[i+1]/multigratorOmega[i+1];
          multigratorNu[i] = exp(-0.5*t*NuOmega)*multigratorNuT[i] + 
          0.5*t*exp(-0.25*t*NuOmega)*(pow(multigratorNu[i-1],2.0)/multigratorOmega[i-1] - kT0);
        }
        multigratorZeta[i] += t*multigratorNuT[i]/multigratorOmega[i];
      }
    }

  }
}

// Calculates Enthalpy for multigrator
BigReal Controller::multigatorCalcEnthalpy(BigReal potentialEnergy, int step, int minimize) {
  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  SimParameters *simParameters = node->simParameters;

  BigReal V = state->lattice.volume();
  BigReal P0 = simParams->multigratorPressureTarget;
  BigReal kT0 = BOLTZMANN * simParams->multigratorTemperatureTarget;
  BigReal NG = controlNumDegFreedom;
  BigReal Nf = numDegFreedom;
  BigReal tauV = simParams->multigratorPressureRelaxationTime;
  BigReal sumZeta = 0.0;
  for (int i=1;i < simParams->multigratorNoseHooverChainLength;i++) {
    sumZeta += multigratorZeta[i];
  }
  BigReal nuOmegaSum = 0.0;
  for (int i=0;i < simParams->multigratorNoseHooverChainLength;i++) {
    nuOmegaSum += multigratorNu[i]*multigratorNu[i]/(2.0*multigratorOmega[i]);
  }
  BigReal W = (3.0*NG + 1.0)*kT0*tauV*tauV;
  BigReal eta = sqrt(kT0*W)*multigratorXi;

  BigReal enthalpy = kineticEnergy + potentialEnergy + eta*eta/(2.0*W) + P0*V + nuOmegaSum + kT0*(Nf*multigratorZeta[0] + sumZeta);

//  if (!(step % 100))
    // fprintf(stderr, "enthalpy %lf %lf %lf %lf %lf %lf %lf\n", enthalpy,
    //   kineticEnergy, potentialEnergy, eta*eta/(2.0*W), P0*V, nuOmegaSum, kT0*(Nf*multigratorZeta[0] + sumZeta));

  return enthalpy;
}

void Controller::berendsenPressure(int step)
{
  if ( simParams->berendsenPressureOn ) {
   berendsenPressure_count += 1;
   berendsenPressure_avg += controlPressure;
   const int freq = simParams->berendsenPressureFreq;
   if ( ! (berendsenPressure_count % freq) ) {
    Tensor factor = berendsenPressure_avg / berendsenPressure_count;
    berendsenPressure_avg = 0;
    berendsenPressure_count = 0;
    // We only use on-diagonal terms (for now)
    factor = Tensor::diagonal(diagonal(factor));
    factor -= Tensor::identity(simParams->berendsenPressureTarget);
    factor *= ( ( simParams->berendsenPressureCompressibility / 3.0 ) *
       simParams->dt * freq / simParams->berendsenPressureRelaxationTime );
    factor += Tensor::identity(1.0);
#define LIMIT_SCALING(VAR,MIN,MAX,FLAG) {\
         if ( VAR < (MIN) ) { VAR = (MIN); FLAG = 1; } \
         if ( VAR > (MAX) ) { VAR = (MAX); FLAG = 1; } }
    int limited = 0;
    LIMIT_SCALING(factor.xx,1./1.03,1.03,limited)
    LIMIT_SCALING(factor.yy,1./1.03,1.03,limited)
    LIMIT_SCALING(factor.zz,1./1.03,1.03,limited)
#undef LIMIT_SCALING
    if ( limited ) {
      iout << iERROR << "Step " << step <<
	" cell rescaling factor limited.\n" << endi;
    }
    broadcast->positionRescaleFactor.publish(step,factor);
    state->lattice.rescale(factor);
   }
  } else {
    berendsenPressure_avg = 0;
    berendsenPressure_count = 0;
  }
}

void Controller::langevinPiston1(int step)
{
  if ( simParams->langevinPistonOn )
  {
    Tensor &strainRate = langevinPiston_strainRate;
    int cellDims = simParams->useFlexibleCell ? 1 : 3;
    BigReal dt = simParams->dt;
    BigReal dt_long = slowFreq * dt;
    BigReal kT = BOLTZMANN * simParams->langevinPistonTemp;
    BigReal tau = simParams->langevinPistonPeriod;
    BigReal mass = controlNumDegFreedom * kT * tau * tau * cellDims;

#ifdef DEBUG_PRESSURE
    iout << iINFO << "entering langevinPiston1, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( (step-1) % slowFreq ) )
    {
      BigReal gamma = 1 / simParams->langevinPistonDecay;
      BigReal f1 = exp( -0.5 * dt_long * gamma );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kT / mass );
      strainRate *= f1;
      if ( simParams->useFlexibleCell ) {
        // We only use on-diagonal terms (for now)
        if ( simParams->useConstantRatio ) {
	  BigReal r = f2 * random->gaussian();
	  strainRate.xx += r;
	  strainRate.yy += r;
	  strainRate.zz += f2 * random->gaussian();
        } else {
	  strainRate += f2 * Tensor::diagonal(random->gaussian_vector());
        }
      } else {
	strainRate += f2 * Tensor::identity(random->gaussian());
      }

#ifdef DEBUG_PRESSURE
      iout << iINFO << "applying langevin, strain rate: " << strainRate << "\n";
#endif
    }

    // Apply surface tension.  If surfaceTensionTarget is zero, we get
    // the default (isotropic pressure) case.
    
    Tensor ptarget;
    ptarget.xx = ptarget.yy = ptarget.zz = simParams->langevinPistonTarget;
    if ( simParams->surfaceTensionTarget != 0. ) {
      ptarget.xx = ptarget.yy = simParams->langevinPistonTarget - 
        simParams->surfaceTensionTarget / state->lattice.c().z;
    }

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
      ( controlPressure - ptarget );

    if (simParams->fixCellDims) {
      if (simParams->fixCellDimX) strainRate.xx = 0;
      if (simParams->fixCellDimY) strainRate.yy = 0;
      if (simParams->fixCellDimZ) strainRate.zz = 0;
    }


#ifdef DEBUG_PRESSURE
    iout << iINFO << "integrating half step, strain rate: " << strainRate << "\n";
#endif

    if ( simParams->langevinPistonBarrier ) {
    if ( ! ( (step-1-slowFreq/2) % slowFreq ) )
    {
      // We only use on-diagonal terms (for now)
      Tensor factor;
      if ( !simParams->useConstantArea ) {
        factor.xx = exp( dt_long * strainRate.xx );
        factor.yy = exp( dt_long * strainRate.yy );
      } else {
        factor.xx = factor.yy = 1;
      }
      factor.zz = exp( dt_long * strainRate.zz );
      broadcast->positionRescaleFactor.publish(step,factor);
      state->lattice.rescale(factor);
#ifdef DEBUG_PRESSURE
      iout << iINFO << "rescaling by: " << factor << "\n";
#endif
    }
    } else { // langevinPistonBarrier
    if ( ! ( (step-1-slowFreq/2) % slowFreq ) )
    {
      if ( ! positionRescaleFactor.xx ) {  // first pass without MTS
      // We only use on-diagonal terms (for now)
      Tensor factor;
      if ( !simParams->useConstantArea ) {
        factor.xx = exp( dt_long * strainRate.xx );
        factor.yy = exp( dt_long * strainRate.yy );
      } else {
        factor.xx = factor.yy = 1;
      }
      factor.zz = exp( dt_long * strainRate.zz );
      broadcast->positionRescaleFactor.publish(step,factor);
      positionRescaleFactor = factor;
      strainRate_old = strainRate;
      }
      state->lattice.rescale(positionRescaleFactor);
#ifdef DEBUG_PRESSURE
      iout << iINFO << "rescaling by: " << factor << "\n";
#endif
    }
    if ( ! ( (step-slowFreq/2) % slowFreq ) )
    {
      // We only use on-diagonal terms (for now)
      Tensor factor;
      if ( !simParams->useConstantArea ) {
        factor.xx = exp( (dt+dt_long) * strainRate.xx - dt * strainRate_old.xx );
        factor.yy = exp( (dt+dt_long) * strainRate.yy - dt * strainRate_old.yy );
      } else {
        factor.xx = factor.yy = 1;
      }
      factor.zz = exp( (dt+dt_long) * strainRate.zz - dt * strainRate_old.zz );
      broadcast->positionRescaleFactor.publish(step+1,factor);
      positionRescaleFactor = factor;
      strainRate_old = strainRate;
    }
    }

    // corrections to integrator
    if ( ! ( step % nbondFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for nbond, ";
#endif
      strainRate -= ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (nbondFreq - 1) * controlPressure_nbond );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if ( ! ( step % slowFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for slow, ";
#endif
      strainRate -= ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (slowFreq - 1) * controlPressure_slow );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if (simParams->fixCellDims) {
      if (simParams->fixCellDimX) strainRate.xx = 0;
      if (simParams->fixCellDimY) strainRate.yy = 0;
      if (simParams->fixCellDimZ) strainRate.zz = 0;
    }

  }

}

void Controller::langevinPiston2(int step)
{
  if ( simParams->langevinPistonOn )
  {
    Tensor &strainRate = langevinPiston_strainRate;
    int cellDims = simParams->useFlexibleCell ? 1 : 3;
    BigReal dt = simParams->dt;
    BigReal dt_long = slowFreq * dt;
    BigReal kT = BOLTZMANN * simParams->langevinPistonTemp;
    BigReal tau = simParams->langevinPistonPeriod;
    BigReal mass = controlNumDegFreedom * kT * tau * tau * cellDims;

    // corrections to integrator
    if ( ! ( step % nbondFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for nbond, ";
#endif
      strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (nbondFreq - 1) * controlPressure_nbond );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }
    if ( ! ( step % slowFreq ) )
    {
#ifdef DEBUG_PRESSURE
      iout << iINFO << "correcting strain rate for slow, ";
#endif
      strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
		( (slowFreq - 1) * controlPressure_slow );
#ifdef DEBUG_PRESSURE
      iout << "strain rate: " << strainRate << "\n";
#endif
    }

    // Apply surface tension.  If surfaceTensionTarget is zero, we get
    // the default (isotropic pressure) case.
   
    Tensor ptarget;
    ptarget.xx = ptarget.yy = ptarget.zz = simParams->langevinPistonTarget;
    if ( simParams->surfaceTensionTarget != 0. ) {
      ptarget.xx = ptarget.yy = simParams->langevinPistonTarget - 
        simParams->surfaceTensionTarget / state->lattice.c().z;
    }

    strainRate += ( 0.5 * dt * cellDims * state->lattice.volume() / mass ) *
      ( controlPressure - ptarget );

 
#ifdef DEBUG_PRESSURE
    iout << iINFO << "integrating half step, strain rate: " << strainRate << "\n";
#endif

    if ( ! ( step % slowFreq ) )
    {
      BigReal gamma = 1 / simParams->langevinPistonDecay;
      BigReal f1 = exp( -0.5 * dt_long * gamma );
      BigReal f2 = sqrt( ( 1. - f1*f1 ) * kT / mass );
      strainRate *= f1;
      if ( simParams->useFlexibleCell ) {
        // We only use on-diagonal terms (for now)
        if ( simParams->useConstantRatio ) {
	  BigReal r = f2 * random->gaussian();
	  strainRate.xx += r;
	  strainRate.yy += r;
	  strainRate.zz += f2 * random->gaussian();
        } else {
	  strainRate += f2 * Tensor::diagonal(random->gaussian_vector());
        }
      } else {
	strainRate += f2 * Tensor::identity(random->gaussian());
      }
#ifdef DEBUG_PRESSURE
      iout << iINFO << "applying langevin, strain rate: " << strainRate << "\n";
#endif
    }

#ifdef DEBUG_PRESSURE
    iout << iINFO << "exiting langevinPiston2, strain rate: " << strainRate << "\n";
#endif
    if (simParams->fixCellDims) {
      if (simParams->fixCellDimX) strainRate.xx = 0;
      if (simParams->fixCellDimY) strainRate.yy = 0;
      if (simParams->fixCellDimZ) strainRate.zz = 0;
    }
  }


}

void Controller::rescaleVelocities(int step)
{
  const int rescaleFreq = simParams->rescaleFreq;
  if ( rescaleFreq > 0 ) {
    rescaleVelocities_sumTemps += temperature;  ++rescaleVelocities_numTemps;
    if ( rescaleVelocities_numTemps == rescaleFreq ) {
      BigReal avgTemp = rescaleVelocities_sumTemps / rescaleVelocities_numTemps;
      BigReal rescaleTemp = simParams->rescaleTemp;
      if ( simParams->adaptTempOn && simParams->adaptTempRescale && (step > simParams->adaptTempFirstStep ) && 
                (!(simParams->adaptTempLastStep > 0) || step < simParams->adaptTempLastStep )) {
        rescaleTemp = adaptTempT;
      }
      BigReal factor = sqrt(rescaleTemp/avgTemp);
      broadcast->velocityRescaleFactor.publish(step,factor);
      //iout << "RESCALING VELOCITIES AT STEP " << step
      //     << " FROM AVERAGE TEMPERATURE OF " << avgTemp
      //     << " TO " << rescaleTemp << " KELVIN.\n" << endi;
      rescaleVelocities_sumTemps = 0;  rescaleVelocities_numTemps = 0;
    }
  }
}

void Controller::correctMomentum(int step) {

    Vector momentum;
    momentum.x = reduction->item(REDUCTION_HALFSTEP_MOMENTUM_X);
    momentum.y = reduction->item(REDUCTION_HALFSTEP_MOMENTUM_Y);
    momentum.z = reduction->item(REDUCTION_HALFSTEP_MOMENTUM_Z);
    const BigReal mass = reduction->item(REDUCTION_MOMENTUM_MASS);

    if ( momentum.length2() == 0. )
      iout << iERROR << "Momentum is exactly zero; probably a bug.\n" << endi;
    if ( mass == 0. )
      NAMD_die("Total mass is zero in Controller::correctMomentum");

    momentum *= (-1./mass);

    broadcast->momentumCorrection.publish(step+slowFreq,momentum);
}

//Modifications for alchemical fep
void Controller::printFepMessage(int step)
{
  if (simParams->alchOn && simParams->alchFepOn 
      && !simParams->alchLambdaFreq) {
    const BigReal alchLambda = simParams->alchLambda;
    const BigReal alchLambda2 = simParams->alchLambda2;
    const BigReal alchTemp = simParams->alchTemp;
    const int alchEquilSteps = simParams->alchEquilSteps;
    iout << "FEP: RESETTING FOR NEW FEP WINDOW "
         << "LAMBDA SET TO " << alchLambda << " LAMBDA2 " << alchLambda2
         << "\nFEP: WINDOW TO HAVE " << alchEquilSteps
         << " STEPS OF EQUILIBRATION PRIOR TO FEP DATA COLLECTION.\n"
         << "FEP: USING CONSTANT TEMPERATURE OF " << alchTemp 
         << " K FOR FEP CALCULATION\n" << endi;
  }
} 
void Controller::printTiMessage(int step)
{
  if (simParams->alchOn && simParams->alchThermIntOn 
      && !simParams->alchLambdaFreq) {
    const BigReal alchLambda = simParams->alchLambda;
    const int alchEquilSteps = simParams->alchEquilSteps;
    iout << "TI: RESETTING FOR NEW WINDOW "
         << "LAMBDA SET TO " << alchLambda 
         << "\nTI: WINDOW TO HAVE " << alchEquilSteps 
         << " STEPS OF EQUILIBRATION PRIOR TO TI DATA COLLECTION.\n" << endi;
  }
} 
//fepe

void Controller::reassignVelocities(int step)
{
  const int reassignFreq = simParams->reassignFreq;
  if ( ( reassignFreq > 0 ) && ! ( step % reassignFreq ) ) {
    BigReal newTemp = simParams->reassignTemp;
    newTemp += ( step / reassignFreq ) * simParams->reassignIncr;
    if ( simParams->reassignIncr > 0.0 ) {
      if ( newTemp > simParams->reassignHold && simParams->reassignHold > 0.0 )
        newTemp = simParams->reassignHold;
    } else {
      if ( newTemp < simParams->reassignHold )
        newTemp = simParams->reassignHold;
    }
    iout << "REASSIGNING VELOCITIES AT STEP " << step
         << " TO " << newTemp << " KELVIN.\n" << endi;
  }
}

void Controller::tcoupleVelocities(int step)
{
  if ( simParams->tCoupleOn )
  {
    const BigReal tCoupleTemp = simParams->tCoupleTemp;
    BigReal coefficient = 1.;
    if ( temperature > 0. ) coefficient = tCoupleTemp/temperature - 1.;
    broadcast->tcoupleCoefficient.publish(step,coefficient);
  }
}

static char *FORMAT(BigReal X)
{
  static char tmp_string[25];
  const double maxnum = 9999999999.9999;
  if ( X > maxnum ) X = maxnum;
  if ( X < -maxnum ) X = -maxnum;
  sprintf(tmp_string," %14.4f",X); 
  return tmp_string;
}

static char *FORMAT(const char *X)
{
  static char tmp_string[25];
  sprintf(tmp_string," %14s",X); 
  return tmp_string;
}

static char *ETITLE(int X)
{
  static char tmp_string[21];
  sprintf(tmp_string,"ENERGY: %7d",X); 
  return  tmp_string;
}

void Controller::receivePressure(int step, int minimize)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    reduction->require();

    Tensor virial_normal;
    Tensor virial_nbond;
    Tensor virial_slow;
#ifdef ALTVIRIAL
    Tensor altVirial_normal;
    Tensor altVirial_nbond;
    Tensor altVirial_slow;
#endif
    Tensor intVirial_normal;
    Tensor intVirial_nbond;
    Tensor intVirial_slow;
    Vector extForce_normal;
    Vector extForce_nbond;
    Vector extForce_slow;

#if 1
    numDegFreedom = molecule->num_deg_freedom();
    int numGroupDegFreedom = molecule->num_group_deg_freedom();
    int numFixedGroups = molecule->num_fixed_groups();
    int numFixedAtoms = molecule->num_fixed_atoms();
#endif
#if 0
    int numAtoms = molecule->numAtoms;
    numDegFreedom = 3 * numAtoms;
    int numGroupDegFreedom = 3 * molecule->numHydrogenGroups;
    int numFixedAtoms =
	( simParameters->fixedAtomsOn ? molecule->numFixedAtoms : 0 );
    int numLonepairs = molecule->numLonepairs;
    int numFixedGroups = ( numFixedAtoms ? molecule->numFixedGroups : 0 );
    if ( numFixedAtoms ) numDegFreedom -= 3 * numFixedAtoms;
    if (numLonepairs) numDegFreedom -= 3 * numLonepairs;
    if ( numFixedGroups ) numGroupDegFreedom -= 3 * numFixedGroups;
    if ( ! ( numFixedAtoms || molecule->numConstraints
	|| simParameters->comMove || simParameters->langevinOn ) ) {
      numDegFreedom -= 3;
      numGroupDegFreedom -= 3;
    }
    if (simParameters->pairInteractionOn) {
      // this doesn't attempt to deal with fixed atoms or constraints
      numDegFreedom = 3 * molecule->numFepInitial;
    }
    int numRigidBonds = molecule->numRigidBonds;
    int numFixedRigidBonds =
	( simParameters->fixedAtomsOn ? molecule->numFixedRigidBonds : 0 );

    // numLonepairs is subtracted here because all lonepairs have a rigid bond
    // to oxygen, but all of the LP degrees of freedom are dealt with above
    numDegFreedom -= ( numRigidBonds - numFixedRigidBonds - numLonepairs);
#endif

    kineticEnergyHalfstep = reduction->item(REDUCTION_HALFSTEP_KINETIC_ENERGY);
    kineticEnergyCentered = reduction->item(REDUCTION_CENTERED_KINETIC_ENERGY);

    BigReal groupKineticEnergyHalfstep = kineticEnergyHalfstep -
	reduction->item(REDUCTION_INT_HALFSTEP_KINETIC_ENERGY);
    BigReal groupKineticEnergyCentered = kineticEnergyCentered -
	reduction->item(REDUCTION_INT_CENTERED_KINETIC_ENERGY);

    BigReal atomTempHalfstep = 2.0 * kineticEnergyHalfstep
					/ ( numDegFreedom * BOLTZMANN );
    BigReal atomTempCentered = 2.0 * kineticEnergyCentered
					/ ( numDegFreedom * BOLTZMANN );
    BigReal groupTempHalfstep = 2.0 * groupKineticEnergyHalfstep
					/ ( numGroupDegFreedom * BOLTZMANN );
    BigReal groupTempCentered = 2.0 * groupKineticEnergyCentered
					/ ( numGroupDegFreedom * BOLTZMANN );

    /*  test code for comparing different temperatures  
    iout << "TEMPTEST: " << step << " " << 
	atomTempHalfstep << " " <<
	atomTempCentered << " " <<
	groupTempHalfstep << " " <<
	groupTempCentered << "\n" << endi;
  iout << "Number of degrees of freedom: " << numDegFreedom << " " <<
    numGroupDegFreedom << "\n" << endi;
     */

    GET_TENSOR(momentumSqrSum,reduction,REDUCTION_MOMENTUM_SQUARED);

    GET_TENSOR(virial_normal,reduction,REDUCTION_VIRIAL_NORMAL);
    GET_TENSOR(virial_nbond,reduction,REDUCTION_VIRIAL_NBOND);
    GET_TENSOR(virial_slow,reduction,REDUCTION_VIRIAL_SLOW);

#ifdef ALTVIRIAL
    GET_TENSOR(altVirial_normal,reduction,REDUCTION_ALT_VIRIAL_NORMAL);
    GET_TENSOR(altVirial_nbond,reduction,REDUCTION_ALT_VIRIAL_NBOND);
    GET_TENSOR(altVirial_slow,reduction,REDUCTION_ALT_VIRIAL_SLOW);
#endif

    GET_TENSOR(intVirial_normal,reduction,REDUCTION_INT_VIRIAL_NORMAL);
    GET_TENSOR(intVirial_nbond,reduction,REDUCTION_INT_VIRIAL_NBOND);
    GET_TENSOR(intVirial_slow,reduction,REDUCTION_INT_VIRIAL_SLOW);

    GET_VECTOR(extForce_normal,reduction,REDUCTION_EXT_FORCE_NORMAL);
    GET_VECTOR(extForce_nbond,reduction,REDUCTION_EXT_FORCE_NBOND);
    GET_VECTOR(extForce_slow,reduction,REDUCTION_EXT_FORCE_SLOW);
    // APH NOTE: These four lines are now done in calcPressure()
    // Vector extPosition = lattice.origin();
    // virial_normal -= outer(extForce_normal,extPosition);
    // virial_nbond -= outer(extForce_nbond,extPosition);
    // virial_slow -= outer(extForce_slow,extPosition);

    kineticEnergy = kineticEnergyCentered;
    temperature = 2.0 * kineticEnergyCentered / ( numDegFreedom * BOLTZMANN );

    if (simParameters->drudeOn) {
      BigReal drudeComKE
        = reduction->item(REDUCTION_DRUDECOM_CENTERED_KINETIC_ENERGY);
      BigReal drudeBondKE
        = reduction->item(REDUCTION_DRUDEBOND_CENTERED_KINETIC_ENERGY);
      int g_bond = 3 * molecule->numDrudeAtoms;
      int g_com = numDegFreedom - g_bond;
      temperature = 2.0 * drudeComKE / (g_com * BOLTZMANN);
      drudeBondTemp = (g_bond!=0 ? (2.*drudeBondKE/(g_bond*BOLTZMANN)) : 0.);
    }

    // Calculate number of degrees of freedom (controlNumDegFreedom)
    if ( simParameters->useGroupPressure )
    {
      controlNumDegFreedom = molecule->numHydrogenGroups - numFixedGroups;
      if ( ! ( numFixedAtoms || molecule->numConstraints
  || simParameters->comMove || simParameters->langevinOn ) ) {
        controlNumDegFreedom -= 1;
      }
    }
    else
    {
      controlNumDegFreedom = numDegFreedom / 3;
    }
    if (simParameters->fixCellDims) {
      if (simParameters->fixCellDimX) controlNumDegFreedom -= 1;
      if (simParameters->fixCellDimY) controlNumDegFreedom -= 1;
      if (simParameters->fixCellDimZ) controlNumDegFreedom -= 1;
    }

    // Calculate pressure tensors using virials
    calcPressure(step, minimize,
      virial_normal, virial_nbond, virial_slow,
      intVirial_normal, intVirial_nbond, intVirial_slow,
      extForce_normal, extForce_nbond, extForce_slow);

#ifdef DEBUG_PRESSURE
    iout << iINFO << "Control pressure = " << controlPressure <<
      " = " << controlPressure_normal << " + " <<
      controlPressure_nbond << " + " << controlPressure_slow << "\n" << endi;
#endif
    if ( simParams->accelMDOn && simParams->accelMDDebugOn && ! (step % simParameters->accelMDOutFreq) ) {
        iout << " PRESS " << trace(pressure)*PRESSUREFACTOR/3. << "\n"
             << " PNORMAL " << trace(pressure_normal)*PRESSUREFACTOR/3. << "\n"
             << " PNBOND " << trace(pressure_nbond)*PRESSUREFACTOR/3. << "\n"
             << " PSLOW " << trace(pressure_slow)*PRESSUREFACTOR/3. << "\n"
             << " PAMD " << trace(pressure_amd)*PRESSUREFACTOR/3. << "\n"
             << " GPRESS " << trace(groupPressure)*PRESSUREFACTOR/3. << "\n"
             << " GPNORMAL " << trace(groupPressure_normal)*PRESSUREFACTOR/3. << "\n"
             << " GPNBOND " << trace(groupPressure_nbond)*PRESSUREFACTOR/3. << "\n"
             << " GPSLOW " << trace(groupPressure_slow)*PRESSUREFACTOR/3. << "\n"
             << endi;
   }
}

//
// Calculates all pressure tensors using virials
//
// Sets variables:
// pressure, pressure_normal, pressure_nbond, pressure_slow
// groupPressure, groupPressure_normal, groupPressure_nbond, groupPressure_slow
// controlPressure, controlPressure_normal, controlPressure_nbond, controlPressure_slow
// pressure_amd
// 
void Controller::calcPressure(int step, int minimize,
  const Tensor& virial_normal_in, const Tensor& virial_nbond_in, const Tensor& virial_slow_in,
  const Tensor& intVirial_normal, const Tensor& intVirial_nbond, const Tensor& intVirial_slow,
  const Vector& extForce_normal, const Vector& extForce_nbond, const Vector& extForce_slow) {

  Tensor virial_normal = virial_normal_in;
  Tensor virial_nbond = virial_nbond_in;
  Tensor virial_slow = virial_slow_in;

  // Tensor tmp = virial_normal;
  // fprintf(stderr, "%1.2lf %1.2lf %1.2lf %1.2lf %1.2lf %1.2lf %1.2lf %1.2lf %1.2lf\n",
  //   tmp.xx, tmp.xy, tmp.xz, tmp.yx, tmp.yy, tmp.yz, tmp.zx, tmp.zy, tmp.zz);

  Node *node = Node::Object();
  Molecule *molecule = node->molecule;
  SimParameters *simParameters = node->simParameters;
  Lattice &lattice = state->lattice;

  BigReal volume;

  Vector extPosition = lattice.origin();
  virial_normal -= outer(extForce_normal,extPosition);
  virial_nbond -= outer(extForce_nbond,extPosition);
  virial_slow -= outer(extForce_slow,extPosition);

  if ( (volume=lattice.volume()) != 0. )
  {

    if (simParameters->LJcorrection && volume) {
#ifdef MEM_OPT_VERSION
      NAMD_bug("LJcorrection not supported in memory optimized build.");
#else
      // Apply tail correction to pressure
      BigReal alchLambda = simParameters->getCurrentLambda(step);
      virial_normal += Tensor::identity(molecule->getVirialTailCorr(alchLambda) / volume);
#endif
    }

    // kinetic energy component included in virials
    pressure_normal = virial_normal / volume;
    groupPressure_normal = ( virial_normal - intVirial_normal ) / volume;

    if (simParameters->accelMDOn) {
      pressure_amd = virial_amd / volume;
      pressure_normal += pressure_amd;
      groupPressure_normal +=  pressure_amd;
    }

    if ( minimize || ! ( step % nbondFreq ) )
    {
      pressure_nbond = virial_nbond / volume;
      groupPressure_nbond = ( virial_nbond - intVirial_nbond ) / volume;
    }

    if ( minimize || ! ( step % slowFreq ) )
    {
      pressure_slow = virial_slow / volume;
      groupPressure_slow = ( virial_slow - intVirial_slow ) / volume;
    }

    pressure = pressure_normal + pressure_nbond + pressure_slow; 
    groupPressure = groupPressure_normal + groupPressure_nbond +
          groupPressure_slow;
  }
  else
  {
    pressure = Tensor();
    groupPressure = Tensor();
  }

  if ( simParameters->useGroupPressure )
  {
    controlPressure_normal = groupPressure_normal;
    controlPressure_nbond = groupPressure_nbond;
    controlPressure_slow = groupPressure_slow;
    controlPressure = groupPressure;
  }
  else
  {
    controlPressure_normal = pressure_normal;
    controlPressure_nbond = pressure_nbond;
    controlPressure_slow = pressure_slow;
    controlPressure = pressure;
  }

  if ( simParameters->useFlexibleCell ) {
    // use symmetric pressure to control rotation
    // controlPressure_normal = symmetric(controlPressure_normal);
    // controlPressure_nbond = symmetric(controlPressure_nbond);
    // controlPressure_slow = symmetric(controlPressure_slow);
    // controlPressure = symmetric(controlPressure);
    // only use on-diagonal components for now
    controlPressure_normal = Tensor::diagonal(diagonal(controlPressure_normal));
    controlPressure_nbond = Tensor::diagonal(diagonal(controlPressure_nbond));
    controlPressure_slow = Tensor::diagonal(diagonal(controlPressure_slow));
    controlPressure = Tensor::diagonal(diagonal(controlPressure));
    if ( simParameters->useConstantRatio ) {
#define AVGXY(T) T.xy = T.yx = 0; T.xx = T.yy = 0.5 * ( T.xx + T.yy );\
   T.xz = T.zx = T.yz = T.zy = 0.5 * ( T.xz + T.yz );
      AVGXY(controlPressure_normal);
      AVGXY(controlPressure_nbond);
      AVGXY(controlPressure_slow);
      AVGXY(controlPressure);
#undef AVGXY
    }
  } else {
    controlPressure_normal =
  Tensor::identity(trace(controlPressure_normal)/3.);
    controlPressure_nbond =
  Tensor::identity(trace(controlPressure_nbond)/3.);
    controlPressure_slow =
  Tensor::identity(trace(controlPressure_slow)/3.);
    controlPressure =
  Tensor::identity(trace(controlPressure)/3.);
  }
}

void Controller::rescaleaccelMD(int step, int minimize)
{
    if ( !simParams->accelMDOn ) return;

    amd_reduction->require();

    // copy all to submit_reduction
    for ( int i=0; i<REDUCTION_MAX_RESERVED; ++i ) {
      submit_reduction->item(i) += amd_reduction->item(i);
    }
    submit_reduction->submit();

    if (step == simParams->firstTimestep) accelMDdVAverage = 0;
//    if ( minimize || ((step < simParams->accelMDFirstStep ) || (step > simParams->accelMDLastStep ))) return;
    if ( minimize || (step < simParams->accelMDFirstStep ) || ( simParams->accelMDLastStep > 0 && step > simParams->accelMDLastStep )) return;

    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    Lattice &lattice = state->lattice;

    const BigReal accelMDE = simParams->accelMDE;
    const BigReal accelMDalpha = simParams->accelMDalpha;
    const BigReal accelMDTE = simParams->accelMDTE;
    const BigReal accelMDTalpha = simParams->accelMDTalpha;
    const int accelMDOutFreq = simParams->accelMDOutFreq;

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal crosstermEnergy;
    BigReal boundaryEnergy;
    BigReal miscEnergy;
    BigReal amd_electEnergy;
    BigReal amd_ljEnergy;
    BigReal amd_electEnergySlow;
    BigReal amd_groLJEnergy;
    BigReal amd_groGaussEnergy;
    BigReal amd_goTotalEnergy;
    BigReal potentialEnergy;
    BigReal factor_dihe = 1;
    BigReal factor_tot = 1;
    BigReal testV;
    BigReal dV = 0;
    Vector  accelMDfactor;
    Tensor vir;
    Tensor vir_dihe;
    Tensor vir_normal;
    Tensor vir_nbond;
    Tensor vir_slow;

    bondEnergy = amd_reduction->item(REDUCTION_BOND_ENERGY);
    angleEnergy = amd_reduction->item(REDUCTION_ANGLE_ENERGY);
    dihedralEnergy = amd_reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    improperEnergy = amd_reduction->item(REDUCTION_IMPROPER_ENERGY);
    crosstermEnergy = amd_reduction->item(REDUCTION_CROSSTERM_ENERGY);
    boundaryEnergy = amd_reduction->item(REDUCTION_BC_ENERGY);
    miscEnergy = amd_reduction->item(REDUCTION_MISC_ENERGY);

    GET_TENSOR(vir_dihe,amd_reduction,REDUCTION_VIRIAL_AMD_DIHE);
    GET_TENSOR(vir_normal,amd_reduction,REDUCTION_VIRIAL_NORMAL);
    GET_TENSOR(vir_nbond,amd_reduction,REDUCTION_VIRIAL_NBOND);
    GET_TENSOR(vir_slow,amd_reduction,REDUCTION_VIRIAL_SLOW);

    if ( !( step % nbondFreq ) ) {
      amd_electEnergy = amd_reduction->item(REDUCTION_ELECT_ENERGY);
      amd_ljEnergy = amd_reduction->item(REDUCTION_LJ_ENERGY);
      amd_groLJEnergy = amd_reduction->item(REDUCTION_GRO_LJ_ENERGY);
      amd_groGaussEnergy = amd_reduction->item(REDUCTION_GRO_GAUSS_ENERGY);
      BigReal goNativeEnergy = amd_reduction->item(REDUCTION_GO_NATIVE_ENERGY);
      BigReal goNonnativeEnergy = amd_reduction->item(REDUCTION_GO_NONNATIVE_ENERGY);
      amd_goTotalEnergy = goNativeEnergy + goNonnativeEnergy;
    } else {
      amd_electEnergy = electEnergy;
      amd_ljEnergy = ljEnergy;
      amd_groLJEnergy = groLJEnergy;
      amd_groGaussEnergy = groGaussEnergy;
      amd_goTotalEnergy = goTotalEnergy;
    }

    if ( !( step % slowFreq ) ) {
      amd_electEnergySlow = amd_reduction->item(REDUCTION_ELECT_ENERGY_SLOW);
    } else {
      amd_electEnergySlow = electEnergySlow;
    }

    potentialEnergy = bondEnergy + angleEnergy + dihedralEnergy +
        improperEnergy + amd_electEnergy + amd_electEnergySlow + amd_ljEnergy +
        crosstermEnergy + boundaryEnergy + miscEnergy +
        amd_goTotalEnergy + amd_groLJEnergy + amd_groGaussEnergy;

    if (simParams->accelMDdihe) {

        testV = dihedralEnergy + crosstermEnergy;
        if ( testV < accelMDE ) {
           factor_dihe = accelMDalpha/(accelMDalpha + accelMDE - testV);
           factor_dihe *= factor_dihe;
           vir = vir_dihe * (factor_dihe - 1.0);
           dV = (accelMDE - testV)*(accelMDE - testV)/(accelMDalpha + accelMDE - testV);
           accelMDdVAverage += dV;
        }  

    } else if (simParams->accelMDdual) {

        testV = dihedralEnergy + crosstermEnergy;
        if ( testV < accelMDE ) {
           factor_dihe = accelMDalpha/(accelMDalpha + accelMDE - testV);
           factor_dihe *= factor_dihe;
           vir = vir_dihe * (factor_dihe - 1.0) ;
           dV = (accelMDE - testV)*(accelMDE - testV)/(accelMDalpha + accelMDE - testV);
        }
 
        testV = potentialEnergy - dihedralEnergy - crosstermEnergy;
        if ( testV < accelMDTE ) {
           factor_tot = accelMDTalpha/(accelMDTalpha + accelMDTE - testV);
           factor_tot *= factor_tot;
           vir += (vir_normal - vir_dihe + vir_nbond + vir_slow) * (factor_tot - 1.0);
           dV += (accelMDTE - testV)*(accelMDTE - testV)/(accelMDTalpha + accelMDTE - testV);
        }
       accelMDdVAverage += dV;

    } else {

        testV = potentialEnergy;
        if ( testV < accelMDE ) {
           factor_tot = accelMDalpha/(accelMDalpha + accelMDE - testV);
           factor_tot *= factor_tot;
           vir = (vir_normal + vir_nbond + vir_slow) * (factor_tot - 1.0);
           dV = (accelMDE - testV)*(accelMDE - testV)/(accelMDalpha + accelMDE - testV);
           accelMDdVAverage += dV;
        }
    } 
 
    accelMDfactor[0]=factor_dihe;
    accelMDfactor[1]=factor_tot;
    accelMDfactor[2]=1;
    broadcast->accelMDRescaleFactor.publish(step,accelMDfactor);
    virial_amd = vir; 

    if ( factor_tot < 0.001 ) {
       iout << iWARN << "accelMD is using a very high boost potential, simulation may become unstable!"
            << "\n" << endi;
    }

    if ( ! (step % accelMDOutFreq) ) {
        if ( !(step == simParams->firstTimestep) ) {
             accelMDdVAverage = accelMDdVAverage/accelMDOutFreq; 
        }
        iout << "ACCELERATED MD: STEP " << step
             << " dV "   << dV
             << " dVAVG " << accelMDdVAverage 
             << " BOND " << bondEnergy
             << " ANGLE " << angleEnergy
             << " DIHED " << dihedralEnergy+crosstermEnergy
             << " IMPRP " << improperEnergy
             << " ELECT " << amd_electEnergy+amd_electEnergySlow
             << " VDW "  << amd_ljEnergy
             << " POTENTIAL "  << potentialEnergy << "\n"
             << endi;

        accelMDdVAverage = 0;

        if (simParams->accelMDDebugOn) {
           BigReal volume;
           Tensor p_normal;
           Tensor p_nbond;
           Tensor p_slow;
           Tensor p;
           if ( (volume=lattice.volume()) != 0. )  {
                 p_normal = vir_normal/volume;
                 p_nbond  = vir_nbond/volume;
                 p_slow = vir_slow/volume;
                 p = vir/volume;
           }    
           iout << " accelMD Scaling Factor: " << accelMDfactor << "\n" 
                << " accelMD PNORMAL " << trace(p_normal)*PRESSUREFACTOR/3. << "\n"
                << " accelMD PNBOND " << trace(p_nbond)*PRESSUREFACTOR/3. << "\n"
                << " accelMD PSLOW " << trace(p_slow)*PRESSUREFACTOR/3. << "\n"
                << " accelMD PAMD " << trace(p)*PRESSUREFACTOR/3. << "\n" 
                << endi;
        }
   }
}

void Controller::adaptTempInit(int step) {
    if (!simParams->adaptTempOn) return;
    iout << iINFO << "INITIALISING ADAPTIVE TEMPERING\n" << endi;
    adaptTempDtMin = 0;
    adaptTempDtMax = 0;
    adaptTempAutoDt = false;
    if (simParams->adaptTempBins == 0) {
      iout << iINFO << "READING ADAPTIVE TEMPERING RESTART FILE\n" << endi;
      std::ifstream adaptTempRead(simParams->adaptTempInFile);
      if (adaptTempRead) {
        int readInt;
        BigReal readReal;
        bool readBool;
        // step
        adaptTempRead >> readInt;
        // Start with min and max temperatures
        adaptTempRead >> adaptTempT;     // KELVIN
        adaptTempRead >> adaptTempBetaMin;  // KELVIN
        adaptTempRead >> adaptTempBetaMax;  // KELVIN
        adaptTempBetaMin = 1./adaptTempBetaMin; // KELVIN^-1
        adaptTempBetaMax = 1./adaptTempBetaMax; // KELVIN^-1
        // In case file is manually edited
        if (adaptTempBetaMin > adaptTempBetaMax){
            readReal = adaptTempBetaMax;
            adaptTempBetaMax = adaptTempBetaMin;
            adaptTempBetaMin = adaptTempBetaMax;
        }
        adaptTempRead >> adaptTempBins;     
        adaptTempRead >> adaptTempCg;
        adaptTempRead >> adaptTempDt;
        adaptTempPotEnergyAveNum = new BigReal[adaptTempBins];
        adaptTempPotEnergyAveDen = new BigReal[adaptTempBins];
        adaptTempPotEnergySamples = new int[adaptTempBins];
        adaptTempPotEnergyVarNum = new BigReal[adaptTempBins];
        adaptTempPotEnergyVar    = new BigReal[adaptTempBins];
        adaptTempPotEnergyAve    = new BigReal[adaptTempBins];
        adaptTempBetaN           = new BigReal[adaptTempBins];
        adaptTempDBeta = (adaptTempBetaMax - adaptTempBetaMin)/(adaptTempBins);
        for(int j = 0; j < adaptTempBins; ++j) {
          adaptTempRead >> adaptTempPotEnergyAve[j];
          adaptTempRead >> adaptTempPotEnergyVar[j];
          adaptTempRead >> adaptTempPotEnergySamples[j];
          adaptTempRead >> adaptTempPotEnergyAveNum[j];
          adaptTempRead >> adaptTempPotEnergyVarNum[j];
          adaptTempRead >> adaptTempPotEnergyAveDen[j];
          adaptTempBetaN[j] = adaptTempBetaMin + j * adaptTempDBeta;
        } 
        adaptTempRead.close();
      }
      else NAMD_die("Could not open ADAPTIVE TEMPERING restart file.\n");
    } 
    else {
      adaptTempBins = simParams->adaptTempBins;
      adaptTempPotEnergyAveNum = new BigReal[adaptTempBins];
      adaptTempPotEnergyAveDen = new BigReal[adaptTempBins];
      adaptTempPotEnergySamples = new int[adaptTempBins];
      adaptTempPotEnergyVarNum = new BigReal[adaptTempBins];
      adaptTempPotEnergyVar    = new BigReal[adaptTempBins];
      adaptTempPotEnergyAve    = new BigReal[adaptTempBins];
      adaptTempBetaN           = new BigReal[adaptTempBins];
      adaptTempBetaMax = 1./simParams->adaptTempTmin;
      adaptTempBetaMin = 1./simParams->adaptTempTmax;
      adaptTempCg = simParams->adaptTempCgamma;   
      adaptTempDt  = simParams->adaptTempDt;
      adaptTempDBeta = (adaptTempBetaMax - adaptTempBetaMin)/(adaptTempBins);
      adaptTempT = simParams->initialTemp; 
      if (simParams->langevinOn)
        adaptTempT = simParams->langevinTemp;
      else if (simParams->rescaleFreq > 0)
        adaptTempT = simParams->rescaleTemp;
      for(int j = 0; j < adaptTempBins; ++j){
          adaptTempPotEnergyAveNum[j] = 0.;
          adaptTempPotEnergyAveDen[j] = 0.;
          adaptTempPotEnergySamples[j] = 0;
          adaptTempPotEnergyVarNum[j] = 0.;
          adaptTempPotEnergyVar[j] = 0.;
          adaptTempPotEnergyAve[j] = 0.;
          adaptTempBetaN[j] = adaptTempBetaMin + j * adaptTempDBeta;
      }
    }
    if (simParams->adaptTempAutoDt > 0.0) {
       adaptTempAutoDt = true;
       adaptTempDtMin =  simParams->adaptTempAutoDt - 0.01;
       adaptTempDtMax =  simParams->adaptTempAutoDt + 0.01;
    }
    adaptTempDTave = 0;
    adaptTempDTavenum = 0;
    iout << iINFO << "ADAPTIVE TEMPERING: TEMPERATURE RANGE: [" << 1./adaptTempBetaMax << "," << 1./adaptTempBetaMin << "] KELVIN\n";
    iout << iINFO << "ADAPTIVE TEMPERING: NUMBER OF BINS TO STORE POT. ENERGY: " << adaptTempBins << "\n";
    iout << iINFO << "ADAPTIVE TEMPERING: ADAPTIVE BIN AVERAGING PARAMETER: " << adaptTempCg << "\n";
    iout << iINFO << "ADAPTIVE TEMPERING: SCALING CONSTANT FOR LANGEVIN TEMPERATURE UPDATES: " << adaptTempDt << "\n";
    iout << iINFO << "ADAPTIVE TEMPERING: INITIAL TEMPERATURE: " << adaptTempT<< "\n";
    if (simParams->adaptTempRestartFreq > 0) {
      NAMD_backup_file(simParams->adaptTempRestartFile);
      adaptTempRestartFile.open(simParams->adaptTempRestartFile);
       if(!adaptTempRestartFile)
        NAMD_die("Error opening restart file");
    }
}

void Controller::adaptTempWriteRestart(int step) {
    if (simParams->adaptTempOn && !(step%simParams->adaptTempRestartFreq)) {
        adaptTempRestartFile.seekp(std::ios::beg);        
        iout << "ADAPTEMP: WRITING RESTART FILE AT STEP " << step << "\n" << endi;
        adaptTempRestartFile << step << " ";
        // Start with min and max temperatures
        adaptTempRestartFile << adaptTempT<< " ";     // KELVIN
        adaptTempRestartFile << 1./adaptTempBetaMin << " ";  // KELVIN
        adaptTempRestartFile << 1./adaptTempBetaMax << " ";  // KELVIN
        adaptTempRestartFile << adaptTempBins << " ";     
        adaptTempRestartFile << adaptTempCg << " ";
        adaptTempRestartFile << adaptTempDt ;
        adaptTempRestartFile << "\n" ;
        for(int j = 0; j < adaptTempBins; ++j) {
          adaptTempRestartFile << adaptTempPotEnergyAve[j] << " ";
          adaptTempRestartFile << adaptTempPotEnergyVar[j] << " ";
          adaptTempRestartFile << adaptTempPotEnergySamples[j] << " ";
          adaptTempRestartFile << adaptTempPotEnergyAveNum[j] << " ";
          adaptTempRestartFile << adaptTempPotEnergyVarNum[j] << " ";
          adaptTempRestartFile << adaptTempPotEnergyAveDen[j] << " ";
          adaptTempRestartFile << "\n";          
        }
        adaptTempRestartFile.flush(); 
    }
}    

void Controller::adaptTempUpdate(int step, int minimize)
{
    //Beta = 1./T
    if ( !simParams->adaptTempOn ) return;
    int j = 0;
    if (step == simParams->firstTimestep) {
        adaptTempInit(step);
        return;
    }
    if ( minimize || (step < simParams->adaptTempFirstStep ) || 
        ( simParams->adaptTempLastStep > 0 && step > simParams->adaptTempLastStep )) return;
    const int adaptTempOutFreq  = simParams->adaptTempOutFreq;
    const bool adaptTempDebug  = simParams->adaptTempDebug;
    //Calculate Current inverse temperature and bin 
    BigReal adaptTempBeta = 1./adaptTempT;
    adaptTempBin   = (int)floor((adaptTempBeta - adaptTempBetaMin)/adaptTempDBeta);

    if (adaptTempBin < 0 || adaptTempBin > adaptTempBins)
        iout << iWARN << " adaptTempBin out of range: adaptTempBin: " << adaptTempBin  
                               << " adaptTempBeta: " << adaptTempBeta 
                              << " adaptTempDBeta: " << adaptTempDBeta 
                               << " betaMin:" << adaptTempBetaMin 
                               << " betaMax: " << adaptTempBetaMax << "\n";
    adaptTempPotEnergySamples[adaptTempBin] += 1;
    BigReal gammaAve = 1.-adaptTempCg/adaptTempPotEnergySamples[adaptTempBin];

    BigReal potentialEnergy;
    BigReal potEnergyAverage;
    BigReal potEnergyVariance;
    potentialEnergy = totalEnergy-kineticEnergy;

    //calculate new bin average and variance using adaptive averaging
    adaptTempPotEnergyAveNum[adaptTempBin] = adaptTempPotEnergyAveNum[adaptTempBin]*gammaAve + potentialEnergy;
    adaptTempPotEnergyAveDen[adaptTempBin] = adaptTempPotEnergyAveDen[adaptTempBin]*gammaAve + 1;
    adaptTempPotEnergyVarNum[adaptTempBin] = adaptTempPotEnergyVarNum[adaptTempBin]*gammaAve + potentialEnergy*potentialEnergy;
    
    potEnergyAverage = adaptTempPotEnergyAveNum[adaptTempBin]/adaptTempPotEnergyAveDen[adaptTempBin];
    potEnergyVariance = adaptTempPotEnergyVarNum[adaptTempBin]/adaptTempPotEnergyAveDen[adaptTempBin];
    potEnergyVariance -= potEnergyAverage*potEnergyAverage;

    adaptTempPotEnergyAve[adaptTempBin] = potEnergyAverage;
    adaptTempPotEnergyVar[adaptTempBin] = potEnergyVariance;
    
    // Weighted integral of <Delta E^2>_beta dbeta <= Eq 4 of JCP 132 244101
    // Integrals of Eqs 5 and 6 is done as piecewise assuming <Delta E^2>_beta
    // is constant for each bin. This is to estimate <E(beta)> where beta \in
    // (beta_i,beta_{i+1}) using Eq 2 of JCP 132 244101
    if ( ! ( step % simParams->adaptTempFreq ) ) {
      // If adaptTempBin not at the edge, calculate improved average:
      if (adaptTempBin > 0 && adaptTempBin < adaptTempBins-1) {
          // Get Averaging Limits:
          BigReal deltaBeta = 0.04*adaptTempBeta; //0.08 used in paper - make variable
          BigReal betaPlus;
          BigReal betaMinus;
          int     nMinus =0;
          int     nPlus = 0;
          if ( adaptTempBeta-adaptTempBetaMin < deltaBeta ) deltaBeta = adaptTempBeta-adaptTempBetaMin;
          if ( adaptTempBetaMax-adaptTempBeta < deltaBeta ) deltaBeta = adaptTempBetaMax-adaptTempBeta;
          betaMinus = adaptTempBeta - deltaBeta;
          betaPlus =  adaptTempBeta + deltaBeta;
          BigReal betaMinus2 = betaMinus*betaMinus;
          BigReal betaPlus2  = betaPlus*betaPlus;
          nMinus = (int)floor((betaMinus-adaptTempBetaMin)/adaptTempDBeta);
          nPlus  = (int)floor((betaPlus-adaptTempBetaMin)/adaptTempDBeta);
          // Variables for <E(beta)> estimate:
          BigReal potEnergyAve0 = 0.0;
          BigReal potEnergyAve1 = 0.0;
          // Integral terms
          BigReal A0 = 0;
          BigReal A1 = 0;
          BigReal A2 = 0;
          //A0 phi_s integral for beta_minus < beta < beta_{i+1}
          BigReal betaNp1 = adaptTempBetaN[adaptTempBin+1]; 
          BigReal DeltaE2Ave;
          BigReal DeltaE2AveSum = 0;
          BigReal betaj;
          for (j = nMinus+1; j <= adaptTempBin; ++j) {
              potEnergyAve0 += adaptTempPotEnergyAve[j]; 
              DeltaE2Ave = adaptTempPotEnergyVar[j];
              DeltaE2AveSum += DeltaE2Ave;
              A0 += j*DeltaE2Ave;
          }
          A0 *= adaptTempDBeta;
          A0 += (adaptTempBetaMin+0.5*adaptTempDBeta-betaMinus)*DeltaE2AveSum;
          A0 *= adaptTempDBeta;
          betaj = adaptTempBetaN[nMinus+1]-betaMinus; 
          betaj *= betaj;
          A0 += 0.5*betaj*adaptTempPotEnergyVar[nMinus];
          A0 /= (betaNp1-betaMinus);

          //A1 phi_s integral for beta_{i+1} < beta < beta_plus
          DeltaE2AveSum = 0;
          for (j = adaptTempBin+1; j < nPlus; j++) {
              potEnergyAve1 += adaptTempPotEnergyAve[j];
              DeltaE2Ave = adaptTempPotEnergyVar[j];
              DeltaE2AveSum += DeltaE2Ave;
              A1 += j*DeltaE2Ave;
          }
          A1 *= adaptTempDBeta;   
          A1 += (adaptTempBetaMin+0.5*adaptTempDBeta-betaPlus)*DeltaE2AveSum;
          A1 *= adaptTempDBeta;
          if ( nPlus < adaptTempBins ) {
            betaj = betaPlus-adaptTempBetaN[nPlus];
            betaj *= betaj;
            A1 -= 0.5*betaj*adaptTempPotEnergyVar[nPlus];
          }
          A1 /= (betaPlus-betaNp1);

          //A2 phi_t integral for beta_i
          A2 = 0.5*adaptTempDBeta*potEnergyVariance;

          // Now calculate a+ and a-
          DeltaE2Ave = A0-A1;
          BigReal aplus = 0;
          BigReal aminus = 0;
          if (DeltaE2Ave != 0) {
            aplus = (A0-A2)/(A0-A1);
            if (aplus < 0) {
                    aplus = 0;
            }
            if (aplus > 1)  {
                    aplus = 1;
            }
            aminus = 1-aplus;
            potEnergyAve0 *= adaptTempDBeta;
            potEnergyAve0 += adaptTempPotEnergyAve[nMinus]*(adaptTempBetaN[nMinus+1]-betaMinus);
            potEnergyAve0 /= (betaNp1-betaMinus);
            potEnergyAve1 *= adaptTempDBeta;
            if ( nPlus < adaptTempBins ) {
                potEnergyAve1 += adaptTempPotEnergyAve[nPlus]*(betaPlus-adaptTempBetaN[nPlus]);
            }
            potEnergyAve1 /= (betaPlus-betaNp1);
            potEnergyAverage = aminus*potEnergyAve0;
            potEnergyAverage += aplus*potEnergyAve1;
          }
          if (simParams->adaptTempDebug) {
       iout << "ADAPTEMP DEBUG:"  << "\n"
            << "     adaptTempBin:    " << adaptTempBin << "\n"
            << "     Samples:   " << adaptTempPotEnergySamples[adaptTempBin] << "\n"
            << "     adaptTempBeta:   " << adaptTempBeta << "\n" 
            << "     adaptTemp:   " << adaptTempT<< "\n"
            << "     betaMin:   " << adaptTempBetaMin << "\n"
            << "     betaMax:   " << adaptTempBetaMax << "\n"
            << "     gammaAve:  " << gammaAve << "\n"
            << "     deltaBeta: " << deltaBeta << "\n"
            << "     betaMinus: " << betaMinus << "\n"
            << "     betaPlus:  " << betaPlus << "\n"
            << "     nMinus:    " << nMinus << "\n"
            << "     nPlus:     " << nPlus << "\n"
            << "     A0:        " << A0 << "\n"
            << "     A1:        " << A1 << "\n"
            << "     A2:        " << A2 << "\n"
            << "     a+:        " << aplus << "\n"
            << "     a-:        " << aminus << "\n"
            << endi;
          }
      }
      else {
          if (simParams->adaptTempDebug) {
       iout << "ADAPTEMP DEBUG:"  << "\n"
            << "     adaptTempBin:    " << adaptTempBin << "\n"
            << "     Samples:   " << adaptTempPotEnergySamples[adaptTempBin] << "\n"
            << "     adaptTempBeta:   " << adaptTempBeta << "\n" 
            << "     adaptTemp:   " << adaptTempT<< "\n"
            << "     betaMin:   " << adaptTempBetaMin << "\n"
            << "     betaMax:   " << adaptTempBetaMax << "\n"
            << "     gammaAve:  " << gammaAve << "\n"
            << "     deltaBeta:  N/A\n"
            << "     betaMinus:  N/A\n"
            << "     betaPlus:   N/A\n"
            << "     nMinus:     N/A\n"
            << "     nPlus:      N/A\n"
            << "     A0:         N/A\n"
            << "     A1:         N/A\n"
            << "     A2:         N/A\n"
            << "     a+:         N/A\n"
            << "     a-:         N/A\n"
            << endi;
          }
      }
      
      //dT is new temperature
      BigReal dT = ((potentialEnergy-potEnergyAverage)/BOLTZMANN+adaptTempT)*adaptTempDt;
      dT += random->gaussian()*sqrt(2.*adaptTempDt)*adaptTempT;
      dT += adaptTempT;
      
     // Check if dT in [adaptTempTmin,adaptTempTmax]. If not try simpler estimate of mean
     // This helps sampling with poor statistics in the bins surrounding adaptTempBin.
      if ( dT > 1./adaptTempBetaMin || dT  < 1./adaptTempBetaMax ) {
        dT = ((potentialEnergy-adaptTempPotEnergyAve[adaptTempBin])/BOLTZMANN+adaptTempT)*adaptTempDt;
        dT += random->gaussian()*sqrt(2.*adaptTempDt)*adaptTempT;
        dT += adaptTempT;
        // Check again, if not then keep original adaptTempTor assign random.
        if ( dT > 1./adaptTempBetaMin ) {
          if (!simParams->adaptTempRandom) {             
             //iout << iWARN << "ADAPTEMP: " << step << " T= " << dT 
             //     << " K higher than adaptTempTmax."
             //     << " Keeping temperature at " 
             //     << adaptTempT<< "\n"<< endi;             
             dT = adaptTempT;
          }
          else {             
             //iout << iWARN << "ADAPTEMP: " << step << " T= " << dT 
             //     << " K higher than adaptTempTmax."
             //     << " Assigning random temperature in range\n" << endi;
             dT = adaptTempBetaMin +  random->uniform()*(adaptTempBetaMax-adaptTempBetaMin);             
             dT = 1./dT;
          }
        } 
        else if ( dT  < 1./adaptTempBetaMax ) {
          if (!simParams->adaptTempRandom) {            
            //iout << iWARN << "ADAPTEMP: " << step << " T= "<< dT 
            //     << " K lower than adaptTempTmin."
            //     << " Keeping temperature at " << adaptTempT<< "\n" << endi; 
            dT = adaptTempT;
          }
          else {
            //iout << iWARN << "ADAPTEMP: " << step << " T= "<< dT 
            //     << " K lower than adaptTempTmin."
            //     << " Assigning random temperature in range\n" << endi;
            dT = adaptTempBetaMin +  random->uniform()*(adaptTempBetaMax-adaptTempBetaMin);
            dT = 1./dT;
          }
        }
        else if (adaptTempAutoDt) {
          //update temperature step size counter
          //FOR "TRUE" ADAPTIVE TEMPERING 
          BigReal adaptTempTdiff = fabs(dT-adaptTempT);
          if (adaptTempTdiff > 0) {
            adaptTempDTave += adaptTempTdiff; 
            adaptTempDTavenum++;
//            iout << "ADAPTEMP: adapTempTdiff = " << adaptTempTdiff << "\n";
          }
          if(adaptTempDTavenum == 100){
                BigReal Frac;
                adaptTempDTave /= adaptTempDTavenum;
                Frac = 1./adaptTempBetaMin-1./adaptTempBetaMax;
                Frac = adaptTempDTave/Frac;
                //if average temperature jump is > 3% of temperature range,
                //modify jump size to match 3%
                iout << "ADAPTEMP: " << step << " FRAC " << Frac << "\n"; 
                if (Frac > adaptTempDtMax || Frac < adaptTempDtMin) {
                    Frac = adaptTempDtMax/Frac;
                    iout << "ADAPTEMP: Updating adaptTempDt to ";
                    adaptTempDt *= Frac;
                    iout << adaptTempDt << "\n" << endi;
                }
                adaptTempDTave = 0;
                adaptTempDTavenum = 0;
          }
        }
      }
      else if (adaptTempAutoDt) {
          //update temperature step size counter
          // FOR "TRUE" ADAPTIVE TEMPERING
          BigReal adaptTempTdiff = fabs(dT-adaptTempT);
          if (adaptTempTdiff > 0) {
            adaptTempDTave += adaptTempTdiff; 
            adaptTempDTavenum++;
//            iout << "ADAPTEMP: adapTempTdiff = " << adaptTempTdiff << "\n";
          }
          if(adaptTempDTavenum == 100){
                BigReal Frac;
                adaptTempDTave /= adaptTempDTavenum;
                Frac = 1./adaptTempBetaMin-1./adaptTempBetaMax;
                Frac = adaptTempDTave/Frac;
                //if average temperature jump is > 3% of temperature range,
                //modify jump size to match 3%
                iout << "ADAPTEMP: " << step << " FRAC " << Frac << "\n"; 
                if (Frac > adaptTempDtMax || Frac < adaptTempDtMin) {
                    Frac = adaptTempDtMax/Frac;
                    iout << "ADAPTEMP: Updating adaptTempDt to ";
                    adaptTempDt *= Frac;
                    iout << adaptTempDt << "\n" << endi;
                }
                adaptTempDTave = 0;
                adaptTempDTavenum = 0;

          }
          
      }
      
      adaptTempT = dT; 
      broadcast->adaptTemperature.publish(step,adaptTempT);
    }
    adaptTempWriteRestart(step);
    if ( ! (step % adaptTempOutFreq) ) {
        iout << "ADAPTEMP: STEP " << step
             << " TEMP "   << adaptTempT
             << " ENERGY " << std::setprecision(10) << potentialEnergy   
             << " ENERGYAVG " << std::setprecision(10) << potEnergyAverage
             << " ENERGYVAR " << std::setprecision(10) << potEnergyVariance;
        iout << "\n" << endi;
   }
   
}


void Controller::compareChecksums(int step, int forgiving) {
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;

    // Some basic correctness checking
    BigReal checksum, checksum_b;
    char errmsg[256];

    checksum = reduction->item(REDUCTION_ATOM_CHECKSUM);
    if ( ((int)checksum) != molecule->numAtoms ) {
      sprintf(errmsg, "Bad global atom count! (%d vs %d)\n",
              (int)checksum, molecule->numAtoms);
      NAMD_bug(errmsg);
    }

    checksum = reduction->item(REDUCTION_COMPUTE_CHECKSUM);
    if ( ((int)checksum) && ((int)checksum) != computeChecksum ) {
      if ( ! computeChecksum ) {
        computesPartitioned = 0;
      } else if ( (int)checksum > computeChecksum && ! computesPartitioned ) {
        computesPartitioned = 1;
      } else {
        NAMD_bug("Bad global compute count!\n");
      }
      computeChecksum = ((int)checksum);
    }

    checksum_b = checksum = reduction->item(REDUCTION_BOND_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcBonds) ) {
      sprintf(errmsg, "Bad global bond count! (%d vs %d)\n",
              (int)checksum, molecule->numCalcBonds);
      if ( forgiving && (((int)checksum) < molecule->numCalcBonds) )
        iout << iWARN << errmsg << endi;
      else NAMD_bug(errmsg);
    }

    checksum = reduction->item(REDUCTION_ANGLE_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcAngles) ) {
      sprintf(errmsg, "Bad global angle count! (%d vs %d)\n",
              (int)checksum, molecule->numCalcAngles);
      if ( forgiving && (((int)checksum) < molecule->numCalcAngles) )
        iout << iWARN << errmsg << endi;
      else NAMD_bug(errmsg);
    }

    checksum = reduction->item(REDUCTION_DIHEDRAL_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcDihedrals) ) {
      sprintf(errmsg, "Bad global dihedral count! (%d vs %d)\n",
              (int)checksum, molecule->numCalcDihedrals);
      if ( forgiving && (((int)checksum) < molecule->numCalcDihedrals) )
        iout << iWARN << errmsg << endi;
      else NAMD_bug(errmsg);
    }

    checksum = reduction->item(REDUCTION_IMPROPER_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcImpropers) ) {
      sprintf(errmsg, "Bad global improper count! (%d vs %d)\n",
              (int)checksum, molecule->numCalcImpropers);
      if ( forgiving && (((int)checksum) < molecule->numCalcImpropers) )
        iout << iWARN << errmsg << endi;
      else NAMD_bug(errmsg);
    }

    checksum = reduction->item(REDUCTION_THOLE_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcTholes) ) {
      sprintf(errmsg, "Bad global Thole count! (%d vs %d)\n",
              (int)checksum, molecule->numCalcTholes);
      if ( forgiving && (((int)checksum) < molecule->numCalcTholes) )
        iout << iWARN << errmsg << endi;
      else NAMD_bug(errmsg);
    }

    checksum = reduction->item(REDUCTION_ANISO_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcAnisos) ) {
      sprintf(errmsg, "Bad global Aniso count! (%d vs %d)\n",
              (int)checksum, molecule->numCalcAnisos);
      if ( forgiving && (((int)checksum) < molecule->numCalcAnisos) )
        iout << iWARN << errmsg << endi;
      else NAMD_bug(errmsg);
    }

    checksum = reduction->item(REDUCTION_CROSSTERM_CHECKSUM);
    if ( checksum_b && (((int)checksum) != molecule->numCalcCrossterms) ) {
      sprintf(errmsg, "Bad global crossterm count! (%d vs %d)\n",
              (int)checksum, molecule->numCalcCrossterms);
      if ( forgiving && (((int)checksum) < molecule->numCalcCrossterms) )
        iout << iWARN << errmsg << endi;
      else NAMD_bug(errmsg);
    }

    checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM);
    if ( ((int)checksum) > molecule->numCalcExclusions &&
         ( ! simParams->mollyOn || step % slowFreq ) ) {
      if ( forgiving )
        iout << iWARN << "High global exclusion count ("
                      << ((int)checksum) << " vs "
                      << molecule->numCalcExclusions << "), possible error!\n"
          << iWARN << "This warning is not unusual during minimization.\n"
          << iWARN << "Decreasing pairlistdist or cutoff that is too close to periodic cell size may avoid this.\n" << endi;
      else {
        char errmsg[256];
        sprintf(errmsg, "High global exclusion count!  (%d vs %d)  System unstable or pairlistdist or cutoff too close to periodic cell size.\n",
                (int)checksum, molecule->numCalcExclusions);
        NAMD_bug(errmsg);
      }
    }

    if ( ((int)checksum) && ((int)checksum) < molecule->numCalcExclusions &&
        ! simParams->goGroPair && ! simParams->qmForcesOn) {
      if ( forgiving || ! simParams->fullElectFrequency ) {
        iout << iWARN << "Low global exclusion count!  ("
          << ((int)checksum) << " vs " << molecule->numCalcExclusions << ")";
        if ( forgiving ) {
          iout << "\n"
            << iWARN << "This warning is not unusual during minimization.\n"
            << iWARN << "Increasing pairlistdist or cutoff may avoid this.\n";
        } else {
          iout << "  System unstable or pairlistdist or cutoff too small.\n";
        }
        iout << endi;
      } else {
        char errmsg[256];
        sprintf(errmsg, "Low global exclusion count!  (%d vs %d)  System unstable or pairlistdist or cutoff too small.\n",
                (int)checksum, molecule->numCalcExclusions);
        NAMD_bug(errmsg);
      }
    }

#ifdef NAMD_CUDA
    checksum = reduction->item(REDUCTION_EXCLUSION_CHECKSUM_CUDA);
    if ( ((int)checksum) > molecule->numCalcFullExclusions &&
         ( ! simParams->mollyOn || step % slowFreq ) ) {
      if ( forgiving )
        iout << iWARN << "High global CUDA exclusion count ("
                      << ((int)checksum) << " vs "
                      << molecule->numCalcFullExclusions << "), possible error!\n"
          << iWARN << "This warning is not unusual during minimization.\n"
          << iWARN << "Decreasing pairlistdist or cutoff that is too close to periodic cell size may avoid this.\n" << endi;
      else {
        char errmsg[256];
        sprintf(errmsg, "High global CUDA exclusion count!  (%d vs %d)  System unstable or pairlistdist or cutoff too close to periodic cell size.\n",
                (int)checksum, molecule->numCalcFullExclusions);
        NAMD_bug(errmsg);
      }
    }

    if ( ((int)checksum) && ((int)checksum) < molecule->numCalcFullExclusions &&
        ! simParams->goGroPair && ! simParams->qmForcesOn) {
      if ( forgiving || ! simParams->fullElectFrequency ) {
        iout << iWARN << "Low global CUDA exclusion count!  ("
          << ((int)checksum) << " vs " << molecule->numCalcFullExclusions << ")";
        if ( forgiving ) {
          iout << "\n"
            << iWARN << "This warning is not unusual during minimization.\n"
            << iWARN << "Increasing pairlistdist or cutoff may avoid this.\n";
        } else {
          iout << "  System unstable or pairlistdist or cutoff too small.\n";
        }
        iout << endi;
      } else {
        char errmsg[256];
        sprintf(errmsg, "Low global CUDA exclusion count!  (%d vs %d)  System unstable or pairlistdist or cutoff too small.\n",
                (int)checksum, molecule->numCalcFullExclusions);
        NAMD_bug(errmsg);
      }
    }
#endif

    checksum = reduction->item(REDUCTION_MARGIN_VIOLATIONS);
    if ( ((int)checksum) && ! marginViolations ) {
      iout << iERROR << "Margin is too small for " << ((int)checksum) <<
        " atoms during timestep " << step << ".\n" << iERROR <<
      "Incorrect nonbonded forces and energies may be calculated!\n" << endi;
    }
    marginViolations += (int)checksum;

    checksum = reduction->item(REDUCTION_PAIRLIST_WARNINGS);
    if ( simParams->outputPairlists && ((int)checksum) && ! pairlistWarnings ) {
      iout << iINFO <<
        "Pairlistdist is too small for " << ((int)checksum) <<
        " computes during timestep " << step << ".\n" << endi;
    }
    if ( simParams->outputPairlists )  pairlistWarnings += (int)checksum;

    checksum = reduction->item(REDUCTION_STRAY_CHARGE_ERRORS);
    if ( checksum ) {
      if ( forgiving )
        iout << iWARN << "Stray PME grid charges ignored!\n" << endi;
      else NAMD_bug("Stray PME grid charges detected!\n");
    }
}

void Controller::printTiming(int step) {

    if ( simParams->outputTiming && ! ( step % simParams->outputTiming ) )
    {
      const double endWTime = CmiWallTimer() - firstWTime;
      const double endCTime = CmiTimer() - firstCTime;

      // fflush about once per minute
      if ( (((int)endWTime)>>6) != (((int)startWTime)>>6 ) ) fflush_count = 2;

      const double elapsedW = 
	(endWTime - startWTime) / simParams->outputTiming;
      const double elapsedC = 
	(endCTime - startCTime) / simParams->outputTiming;

      const double remainingW = elapsedW * (simParams->N - step);
      const double remainingW_hours = remainingW / 3600;

      startWTime = endWTime;
      startCTime = endCTime;

#ifdef NAMD_CUDA
      if ( simParams->outputEnergies < 60 &&
           step < (simParams->firstTimestep + 10 * simParams->outputTiming) ) {
        iout << iWARN << "Energy evaluation is expensive, increase outputEnergies to improve performance.\n" << endi;
      }
#endif
      if ( step >= (simParams->firstTimestep + simParams->outputTiming) ) {
	CmiPrintf("TIMING: %d  CPU: %g, %g/step  Wall: %g, %g/step"
		  ", %g hours remaining, %f MB of memory in use.\n",
		  step, endCTime, elapsedC, endWTime, elapsedW,
		  remainingW_hours, memusage_MB());
        if ( fflush_count ) { --fflush_count; fflush(stdout); }
      }
    }
}

void Controller::printMinimizeEnergies(int step) {

    rescaleaccelMD(step,1);
    receivePressure(step,1);

    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    compareChecksums(step,1);

    printEnergies(step,1);

    min_energy = totalEnergy;
    min_f_dot_f = reduction->item(REDUCTION_MIN_F_DOT_F);
    min_f_dot_v = reduction->item(REDUCTION_MIN_F_DOT_V);
    min_v_dot_v = reduction->item(REDUCTION_MIN_V_DOT_V);
    min_huge_count = (int) (reduction->item(REDUCTION_MIN_HUGE_COUNT));

}

void Controller::printDynamicsEnergies(int step) {

    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    compareChecksums(step);

    printEnergies(step,0);
}

void Controller::printEnergies(int step, int minimize)
{
    Node *node = Node::Object();
    Molecule *molecule = node->molecule;
    SimParameters *simParameters = node->simParameters;
    Lattice &lattice = state->lattice;

    // Drude model ANISO energy is added into BOND energy
    // and THOLE energy is added into ELECT energy

    BigReal bondEnergy;
    BigReal angleEnergy;
    BigReal dihedralEnergy;
    BigReal improperEnergy;
    BigReal crosstermEnergy;
    BigReal boundaryEnergy;
    BigReal miscEnergy;
    BigReal potentialEnergy;
    BigReal flatEnergy;
    BigReal smoothEnergy;

    Vector momentum;
    Vector angularMomentum;
    BigReal volume = lattice.volume();

    bondEnergy = reduction->item(REDUCTION_BOND_ENERGY);
    angleEnergy = reduction->item(REDUCTION_ANGLE_ENERGY);
    dihedralEnergy = reduction->item(REDUCTION_DIHEDRAL_ENERGY);
    improperEnergy = reduction->item(REDUCTION_IMPROPER_ENERGY);
    crosstermEnergy = reduction->item(REDUCTION_CROSSTERM_ENERGY);
    boundaryEnergy = reduction->item(REDUCTION_BC_ENERGY);
    miscEnergy = reduction->item(REDUCTION_MISC_ENERGY);

    if ( minimize || ! ( step % nbondFreq ) )
    {
      electEnergy = reduction->item(REDUCTION_ELECT_ENERGY);
      ljEnergy = reduction->item(REDUCTION_LJ_ENERGY);

      // JLai
      groLJEnergy = reduction->item(REDUCTION_GRO_LJ_ENERGY);
      groGaussEnergy = reduction->item(REDUCTION_GRO_GAUSS_ENERGY);

      goNativeEnergy = reduction->item(REDUCTION_GO_NATIVE_ENERGY);
      goNonnativeEnergy = reduction->item(REDUCTION_GO_NONNATIVE_ENERGY);
      goTotalEnergy = goNativeEnergy + goNonnativeEnergy;

//fepb
      bondedEnergyDiff_f = reduction->item(REDUCTION_BONDED_ENERGY_F);
      electEnergy_f = reduction->item(REDUCTION_ELECT_ENERGY_F);
      ljEnergy_f = reduction->item(REDUCTION_LJ_ENERGY_F);
      ljEnergy_f_left = reduction->item(REDUCTION_LJ_ENERGY_F_LEFT);

      bondedEnergy_ti_1 = reduction->item(REDUCTION_BONDED_ENERGY_TI_1);
      electEnergy_ti_1 = reduction->item(REDUCTION_ELECT_ENERGY_TI_1);
      ljEnergy_ti_1 = reduction->item(REDUCTION_LJ_ENERGY_TI_1);
      bondedEnergy_ti_2 = reduction->item(REDUCTION_BONDED_ENERGY_TI_2);
      electEnergy_ti_2 = reduction->item(REDUCTION_ELECT_ENERGY_TI_2);
      ljEnergy_ti_2 = reduction->item(REDUCTION_LJ_ENERGY_TI_2);
//fepe
    }

    if ( minimize || ! ( step % slowFreq ) )
    {
      electEnergySlow = reduction->item(REDUCTION_ELECT_ENERGY_SLOW);
//fepb
      electEnergySlow_f = reduction->item(REDUCTION_ELECT_ENERGY_SLOW_F);

      electEnergySlow_ti_1 = reduction->item(REDUCTION_ELECT_ENERGY_SLOW_TI_1);
      electEnergySlow_ti_2 = reduction->item(REDUCTION_ELECT_ENERGY_SLOW_TI_2);
      electEnergyPME_ti_1 = reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_1);
      electEnergyPME_ti_2 = reduction->item(REDUCTION_ELECT_ENERGY_PME_TI_2);
//fepe
    }

    if (simParameters->LJcorrection && volume) {
#ifdef MEM_OPT_VERSION
      NAMD_bug("LJcorrection not supported in memory optimized build.");
#else
      // Apply tail correction to energy.
      BigReal alchLambda = simParameters->getCurrentLambda(step);
      BigReal alchLambda2 = simParameters->alchLambda2;

      ljEnergy += molecule->getEnergyTailCorr(alchLambda) / volume;
      ljEnergy_f += molecule->getEnergyTailCorr(alchLambda2) / volume;
      ljEnergy_f_left += molecule->getEnergyTailCorr(alchLambda2) / volume;
      ljEnergy_ti_1 += molecule->tail_corr_dUdl_1 / volume;
      // NB: Rather than duplicate variables, dUdl_2 is stored as the energy.
      //     Put another way, dUdl_2 _is_ the energy if alchLambda = 0.
      ljEnergy_ti_2 += molecule->tail_corr_ener / volume;
#endif
    }

//fepb BKR - Compute alchemical work if using dynamic lambda.  This is here so
//           that the cumulative work can be given during a callback.
    if (simParameters->alchLambdaFreq > 0) {
      if (step <= 
          simParameters->firstTimestep + simParameters->alchEquilSteps) {
        cumAlchWork = 0.0;
      }
      alchWork = computeAlchWork(step);
      cumAlchWork += alchWork;
    }
//fepe

    momentum.x = reduction->item(REDUCTION_MOMENTUM_X);
    momentum.y = reduction->item(REDUCTION_MOMENTUM_Y);
    momentum.z = reduction->item(REDUCTION_MOMENTUM_Z);
    angularMomentum.x = reduction->item(REDUCTION_ANGULAR_MOMENTUM_X);
    angularMomentum.y = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Y);
    angularMomentum.z = reduction->item(REDUCTION_ANGULAR_MOMENTUM_Z);

    // Ported by JLai
    potentialEnergy = (bondEnergy + angleEnergy + dihedralEnergy
  + improperEnergy + electEnergy + electEnergySlow + ljEnergy
  + crosstermEnergy + boundaryEnergy + miscEnergy + goTotalEnergy 
        + groLJEnergy + groGaussEnergy);
    // End of port
    totalEnergy = potentialEnergy + kineticEnergy;
    flatEnergy = (totalEnergy
        + (1.0/3.0)*(kineticEnergyHalfstep - kineticEnergyCentered));
    if ( !(step%slowFreq) ) {
      // only adjust based on most accurate energies
      BigReal s = (4.0/3.0)*( kineticEnergyHalfstep - kineticEnergyCentered);
      if ( smooth2_avg == XXXBIGREAL ) smooth2_avg = s;
      if ( step != simParams->firstTimestep ) {
        smooth2_avg *= 0.9375;
        smooth2_avg += 0.0625 * s;
      }
    }
    smoothEnergy = (flatEnergy + smooth2_avg
        - (4.0/3.0)*(kineticEnergyHalfstep - kineticEnergyCentered));

    if ( simParameters->outputMomenta && ! minimize &&
         ! ( step % simParameters->outputMomenta ) )
    {
      iout << "MOMENTUM: " << step 
           << " P: " << momentum
           << " L: " << angularMomentum
           << "\n" << endi;
    }

    if ( simParameters->outputPressure ) {
      pressure_tavg += pressure;
      groupPressure_tavg += groupPressure;
      tavg_count += 1;
      if ( minimize || ! ( step % simParameters->outputPressure ) ) {
        iout << "PRESSURE: " << step << " "
           << PRESSUREFACTOR * pressure << "\n"
           << "GPRESSURE: " << step << " "
           << PRESSUREFACTOR * groupPressure << "\n";
        if ( tavg_count > 1 ) iout << "PRESSAVG: " << step << " "
           << (PRESSUREFACTOR/tavg_count) * pressure_tavg << "\n"
           << "GPRESSAVG: " << step << " "
           << (PRESSUREFACTOR/tavg_count) * groupPressure_tavg << "\n";
        iout << endi;
        pressure_tavg = 0;
        groupPressure_tavg = 0;
        tavg_count = 0;
      }
    }

    // pressure profile reductions
    if (pressureProfileSlabs) {
      const int freq = simParams->pressureProfileFreq;
      const int arraysize = 3*pressureProfileSlabs;
      
      BigReal *total = new BigReal[arraysize];
      memset(total, 0, arraysize*sizeof(BigReal));
      const int first = simParams->firstTimestep;

      if (ppbonded)    ppbonded->getData(first, step, lattice, total);
      if (ppnonbonded) ppnonbonded->getData(first, step, lattice, total);
      if (ppint)       ppint->getData(first, step, lattice, total);
      for (int i=0; i<arraysize; i++) pressureProfileAverage[i] += total[i];
      pressureProfileCount++;

      if (!(step % freq)) {
        // convert NAMD internal virial to pressure in units of bar 
        BigReal scalefac = PRESSUREFACTOR * pressureProfileSlabs;
   
        iout << "PRESSUREPROFILE: " << step << " ";
        if (step == first) {
          // output pressure profile for this step
          for (int i=0; i<arraysize; i++) {
            iout << total[i] * scalefac << " ";
          }
        } else {
          // output pressure profile averaged over the last count steps.
          scalefac /= pressureProfileCount;
          for (int i=0; i<arraysize; i++) 
            iout << pressureProfileAverage[i]*scalefac << " ";
        }
        iout << "\n" << endi; 
       
        // Clear the average for the next block
        memset(pressureProfileAverage, 0, arraysize*sizeof(BigReal));
        pressureProfileCount = 0;
      }
      delete [] total;
    }
  
   if ( step != simParams->firstTimestep || stepInFullRun == 0 ) {  // skip repeated first step
    if ( stepInFullRun % simParams->firstLdbStep == 0 ) {
     int benchPhase = stepInFullRun / simParams->firstLdbStep;
     if ( benchPhase > 0 && benchPhase < 10 ) {
#ifdef NAMD_CUDA
      if ( simParams->outputEnergies < 60 ) {
        iout << iWARN << "Energy evaluation is expensive, increase outputEnergies to improve performance.\n";
      }
#endif
      iout << iINFO;
      if ( benchPhase < 4 ) iout << "Initial time: ";
      else iout << "Benchmark time: ";
      iout << CkNumPes() << " CPUs ";
      {
        BigReal wallPerStep =
    (CmiWallTimer() - startBenchTime) / simParams->firstLdbStep;
  BigReal ns = simParams->dt / 1000000.0;
  BigReal days = 1.0 / (24.0 * 60.0 * 60.0);
  BigReal daysPerNano = wallPerStep * days / ns;
  iout << wallPerStep << " s/step " << daysPerNano << " days/ns ";
        iout << memusage_MB() << " MB memory\n" << endi;
      }
     }
     startBenchTime = CmiWallTimer();
    }
    ++stepInFullRun;
   }

    printTiming(step);

    Vector pairVDWForce, pairElectForce;
    if ( simParameters->pairInteractionOn ) {
      GET_VECTOR(pairVDWForce,reduction,REDUCTION_PAIR_VDW_FORCE);
      GET_VECTOR(pairElectForce,reduction,REDUCTION_PAIR_ELECT_FORCE);
    }

    // callback to Tcl with whatever we can
#ifdef NAMD_TCL
#define CALLBACKDATA(LABEL,VALUE) \
		labels << (LABEL) << " "; values << (VALUE) << " ";
#define CALLBACKLIST(LABEL,VALUE) \
		labels << (LABEL) << " "; values << "{" << (VALUE) << "} ";
    if (step == simParams->N && node->getScript() && node->getScript()->doCallback()) {
      std::ostringstream labels, values;
#if CMK_BLUEGENEL
      // the normal version below gives a compiler error
      values.precision(16);
#else
      values << std::setprecision(16);
#endif
      CALLBACKDATA("TS",step);
      CALLBACKDATA("BOND",bondEnergy);
      CALLBACKDATA("ANGLE",angleEnergy);
      CALLBACKDATA("DIHED",dihedralEnergy);
      CALLBACKDATA("CROSS",crosstermEnergy);
      CALLBACKDATA("IMPRP",improperEnergy);
      CALLBACKDATA("ELECT",electEnergy+electEnergySlow);
      CALLBACKDATA("VDW",ljEnergy);
      CALLBACKDATA("BOUNDARY",boundaryEnergy);
      CALLBACKDATA("MISC",miscEnergy);
      CALLBACKDATA("POTENTIAL",potentialEnergy);
      CALLBACKDATA("KINETIC",kineticEnergy);
      CALLBACKDATA("TOTAL",totalEnergy);
      CALLBACKDATA("TEMP",temperature);
      CALLBACKLIST("PRESSURE",pressure*PRESSUREFACTOR);
      CALLBACKLIST("GPRESSURE",groupPressure*PRESSUREFACTOR);
      CALLBACKDATA("VOLUME",lattice.volume());
      CALLBACKLIST("CELL_A",lattice.a());
      CALLBACKLIST("CELL_B",lattice.b());
      CALLBACKLIST("CELL_C",lattice.c());
      CALLBACKLIST("CELL_O",lattice.origin());
      labels << "PERIODIC "; values << "{" << lattice.a_p() << " "
		<< lattice.b_p() << " " << lattice.c_p() << "} ";
      if (simParameters->drudeOn) {
        CALLBACKDATA("DRUDEBOND",drudeBondTemp);
      }
      if ( simParameters->pairInteractionOn ) {
        CALLBACKLIST("VDW_FORCE",pairVDWForce);
        CALLBACKLIST("ELECT_FORCE",pairElectForce);
      }
      if (simParameters->alchOn) {
        if (simParameters->alchThermIntOn) {
          CALLBACKLIST("BOND1", bondedEnergy_ti_1);
          CALLBACKLIST("ELEC1", electEnergy_ti_1 + electEnergySlow_ti_1 +
                                electEnergyPME_ti_1);
          CALLBACKLIST("VDW1", ljEnergy_ti_1);
          CALLBACKLIST("BOND2", bondedEnergy_ti_2);
          CALLBACKLIST("ELEC2", electEnergy_ti_2 + electEnergySlow_ti_2 +
                                electEnergyPME_ti_2);
          CALLBACKLIST("VDW2", ljEnergy_ti_2);
          if (simParameters->alchLambdaFreq > 0) {
            CALLBACKLIST("CUMALCHWORK", cumAlchWork);
          }
        } else if (simParameters->alchFepOn) {
          CALLBACKLIST("BOND2", bondEnergy + angleEnergy + dihedralEnergy +
                                improperEnergy + bondedEnergyDiff_f);
          CALLBACKLIST("ELEC2", electEnergy_f + electEnergySlow_f);
          CALLBACKLIST("VDW2", ljEnergy_f);
        } 
      }

      labels << '\0';  values << '\0';  // insane but makes Linux work
      state->callback_labelstring = labels.str();
      state->callback_valuestring = values.str();
      // node->getScript()->doCallback(labelstring.c_str(),valuestring.c_str());
    }
#undef CALLBACKDATA
#endif

    drudeBondTempAvg += drudeBondTemp;

    temp_avg += temperature;
    pressure_avg += trace(pressure)/3.;
    groupPressure_avg += trace(groupPressure)/3.;
    avg_count += 1;

    if ( simParams->outputPairlists && pairlistWarnings &&
				! (step % simParams->outputPairlists) ) {
      iout << iINFO << pairlistWarnings <<
        " pairlist warnings in past " << simParams->outputPairlists <<
	" steps.\n" << endi;
      pairlistWarnings = 0;
    }

    BigReal enthalpy;
    if (simParameters->multigratorOn && ((step % simParameters->outputEnergies) == 0)) {
      enthalpy = multigatorCalcEnthalpy(potentialEnergy, step, minimize);
    }

    // NO CALCULATIONS OR REDUCTIONS BEYOND THIS POINT!!!
    if ( ! minimize &&  step % simParameters->outputEnergies ) return;
    // ONLY OUTPUT SHOULD OCCUR BELOW THIS LINE!!!

    if (simParameters->IMDon && !(step % simParameters->IMDfreq)) {
      IMDEnergies energies;
      energies.tstep = step;
      energies.T = temp_avg/avg_count;
      energies.Etot = totalEnergy;
      energies.Epot = potentialEnergy;
      energies.Evdw = ljEnergy;
      energies.Eelec = electEnergy + electEnergySlow;
      energies.Ebond = bondEnergy;
      energies.Eangle = angleEnergy;
      energies.Edihe = dihedralEnergy + crosstermEnergy;
      energies.Eimpr = improperEnergy;
      Node::Object()->imd->gather_energies(&energies);
    }

    if ( marginViolations ) {
      iout << iERROR << marginViolations <<
        " margin violations detected since previous energy output.\n" << endi;
    }
    marginViolations = 0;

    if ( (step % (10 * (minimize?1:simParameters->outputEnergies) ) ) == 0 )
    {
	iout << "ETITLE:      TS";
	iout << FORMAT("BOND");
	iout << FORMAT("ANGLE");
	iout << FORMAT("DIHED");
	if ( ! simParameters->mergeCrossterms ) iout << FORMAT("CROSS");
	iout << FORMAT("IMPRP");
        iout << "     ";
	iout << FORMAT("ELECT");
	iout << FORMAT("VDW");
	iout << FORMAT("BOUNDARY");
	iout << FORMAT("MISC");
	iout << FORMAT("KINETIC");
        iout << "     ";
	iout << FORMAT("TOTAL");
	iout << FORMAT("TEMP");
	iout << FORMAT("POTENTIAL");
	// iout << FORMAT("TOTAL2");
	iout << FORMAT("TOTAL3");
	iout << FORMAT("TEMPAVG");
	if ( volume != 0. ) {
          iout << "     ";
	  iout << FORMAT("PRESSURE");
	  iout << FORMAT("GPRESSURE");
	  iout << FORMAT("VOLUME");
	  iout << FORMAT("PRESSAVG");
	  iout << FORMAT("GPRESSAVG");
	}
        if (simParameters->drudeOn) {
          iout << "     ";
	  iout << FORMAT("DRUDEBOND");
	  iout << FORMAT("DRBONDAVG");
        }
	// Ported by JLai
	if (simParameters->goGroPair) {
	  iout << "     ";
	  iout << FORMAT("GRO_PAIR_LJ");
	  iout << FORMAT("GRO_PAIR_GAUSS");
	}

	if (simParameters->goForcesOn) {
          iout << "     ";
	  iout << FORMAT("NATIVE");
	  iout << FORMAT("NONNATIVE");
	  //iout << FORMAT("REL_NATIVE");
	  //iout << FORMAT("REL_NONNATIVE");
	  iout << FORMAT("GOTOTAL");
	  //iout << FORMAT("GOAVG");
	}
	// End of port -- JLai

        if (simParameters->alchOn) {
          if (simParameters->alchThermIntOn) {
            iout << "\nTITITLE:     TS";
            iout << FORMAT("BOND1");
            iout << FORMAT("ELECT1");
            iout << FORMAT("VDW1");
            iout << FORMAT("BOND2");
            iout << "     ";
            iout << FORMAT("ELECT2");
            iout << FORMAT("VDW2");
            if (simParameters->alchLambdaFreq > 0) {
              iout << FORMAT("LAMBDA");
              iout << FORMAT("ALCHWORK");
              iout << FORMAT("CUMALCHWORK");
            }
          } else if (simParameters->alchFepOn) {
            iout << "\nFEPTITLE:    TS";
            iout << FORMAT("BOND2");
            iout << FORMAT("ELECT2");
            iout << FORMAT("VDW2");
            if (simParameters->alchLambdaFreq > 0) {
              iout << FORMAT("LAMBDA");
            }
          }
        }

	iout << "\n\n" << endi;

    }

    // N.B.  HP's aCC compiler merges FORMAT calls in the same expression.
    //       Need separate statements because data returned in static array.
    iout << ETITLE(step);
    iout << FORMAT(bondEnergy);
    iout << FORMAT(angleEnergy);
    if ( simParameters->mergeCrossterms ) {
      iout << FORMAT(dihedralEnergy+crosstermEnergy);
    } else {
      iout << FORMAT(dihedralEnergy);
      iout << FORMAT(crosstermEnergy);
    }
    iout << FORMAT(improperEnergy);
    iout << "     ";
    iout << FORMAT(electEnergy+electEnergySlow);
    iout << FORMAT(ljEnergy);
    iout << FORMAT(boundaryEnergy);
    iout << FORMAT(miscEnergy);
    iout << FORMAT(kineticEnergy);
    iout << "     ";
    iout << FORMAT(totalEnergy);
    iout << FORMAT(temperature);
    iout << FORMAT(potentialEnergy);
    // iout << FORMAT(flatEnergy);
    iout << FORMAT(smoothEnergy);
    iout << FORMAT(temp_avg/avg_count);
    if ( volume != 0. )
    {
        iout << "     ";
	iout << FORMAT(trace(pressure)*PRESSUREFACTOR/3.);
	iout << FORMAT(trace(groupPressure)*PRESSUREFACTOR/3.);
	iout << FORMAT(volume);
	iout << FORMAT(pressure_avg*PRESSUREFACTOR/avg_count);
	iout << FORMAT(groupPressure_avg*PRESSUREFACTOR/avg_count);
    }
    if (simParameters->drudeOn) {
        iout << "     ";
	iout << FORMAT(drudeBondTemp);
	iout << FORMAT(drudeBondTempAvg/avg_count);
    }
    // Ported by JLai
    if (simParameters->goGroPair) {
      iout << "     ";
      iout << FORMAT(groLJEnergy);
      iout << FORMAT(groGaussEnergy);
    }

    if (simParameters->goForcesOn) {
      iout << "     ";
      iout << FORMAT(goNativeEnergy);
      iout << FORMAT(goNonnativeEnergy);
      //iout << FORMAT(relgoNativeEnergy);
      //iout << FORMAT(relgoNonnativeEnergy);
      iout << FORMAT(goTotalEnergy);
      //iout << FORMAT("not implemented");
    } // End of port -- JLai

    if (simParameters->alchOn) { 
      if (simParameters->alchThermIntOn) {
        iout << "\n";
        iout << TITITLE(step);
        iout << FORMAT(bondedEnergy_ti_1);
        iout << FORMAT(electEnergy_ti_1 + electEnergySlow_ti_1 + 
                       electEnergyPME_ti_1);
        iout << FORMAT(ljEnergy_ti_1);
        iout << FORMAT(bondedEnergy_ti_2);
        iout << "     ";
        iout << FORMAT(electEnergy_ti_2 + electEnergySlow_ti_2 +
                       electEnergyPME_ti_2);
        iout << FORMAT(ljEnergy_ti_2);
        if (simParameters->alchLambdaFreq > 0) {
          iout << FORMAT(simParameters->getCurrentLambda(step));
          iout << FORMAT(alchWork);
          iout << FORMAT(cumAlchWork);
        }
      } else if (simParameters->alchFepOn) {
        iout << "\n";
        iout << FEPTITLE2(step);
        iout << FORMAT(bondEnergy + angleEnergy + dihedralEnergy 
                       + improperEnergy + bondedEnergyDiff_f);
        iout << FORMAT(electEnergy_f + electEnergySlow_f);
        iout << FORMAT(ljEnergy_f);
        if (simParameters->alchLambdaFreq > 0) {
          iout << FORMAT(simParameters->getCurrentLambda(step));
        }
      }
    }

    iout << "\n\n" << endi;

#if(CMK_CCS_AVAILABLE && CMK_WEB_MODE)
     char webout[80];
     sprintf(webout,"%d %d %d %d",(int)totalEnergy,
	     (int)(potentialEnergy),
	     (int)kineticEnergy,(int)temperature);
     CApplicationDepositNode0Data(webout);
#endif

    if (simParameters->pairInteractionOn) {
      iout << "PAIR INTERACTION:";
      iout << " STEP: " << step;
      iout << " VDW_FORCE: ";
      iout << FORMAT(pairVDWForce.x);
      iout << FORMAT(pairVDWForce.y);
      iout << FORMAT(pairVDWForce.z);
      iout << " ELECT_FORCE: ";
      iout << FORMAT(pairElectForce.x);
      iout << FORMAT(pairElectForce.y);
      iout << FORMAT(pairElectForce.z);
      iout << "\n" << endi;
    }
    drudeBondTempAvg = 0;
    temp_avg = 0;
    pressure_avg = 0;
    groupPressure_avg = 0;
    avg_count = 0;

    if ( fflush_count ) {
      --fflush_count;
      fflush(stdout);
    }
}

void Controller::writeExtendedSystemLabels(ofstream_namd &file) {
  Lattice &lattice = state->lattice;
  file << "#$LABELS step";
  if ( lattice.a_p() ) file << " a_x a_y a_z";
  if ( lattice.b_p() ) file << " b_x b_y b_z";
  if ( lattice.c_p() ) file << " c_x c_y c_z";
  file << " o_x o_y o_z";
  if ( simParams->langevinPistonOn ) {
    file << " s_x s_y s_z s_u s_v s_w";
  }
  file << std::endl;
}

void Controller::writeExtendedSystemData(int step, ofstream_namd &file) {
  Lattice &lattice = state->lattice;
  file.precision(12);
  file << step;
    if ( lattice.a_p() ) file << " " << lattice.a().x << " " << lattice.a().y << " " << lattice.a().z;
    if ( lattice.b_p() ) file << " " << lattice.b().x << " " << lattice.b().y << " " << lattice.b().z;
    if ( lattice.c_p() ) file << " " << lattice.c().x << " " << lattice.c().y << " " << lattice.c().z;
    file << " " << lattice.origin().x << " " << lattice.origin().y << " " << lattice.origin().z;
  if ( simParams->langevinPistonOn ) {
    Vector strainRate = diagonal(langevinPiston_strainRate);
    Vector strainRate2 = off_diagonal(langevinPiston_strainRate);
    file << " " << strainRate.x;
    file << " " << strainRate.y;
    file << " " << strainRate.z;
    file << " " << strainRate2.x;
    file << " " << strainRate2.y;
    file << " " << strainRate2.z;
  }
  file << std::endl;
}

void Controller::enqueueCollections(int timestep)
{
  if ( Output::coordinateNeeded(timestep) ){
    collection->enqueuePositions(timestep,state->lattice);
  }
  if ( Output::velocityNeeded(timestep) )
    collection->enqueueVelocities(timestep);
  if ( Output::forceNeeded(timestep) )
    collection->enqueueForces(timestep);
}

//Modifications for alchemical fep
void Controller::outputFepEnergy(int step) {
 if (simParams->alchOn && simParams->alchFepOn) {
  const int stepInRun = step - simParams->firstTimestep;
  const int alchEquilSteps = simParams->alchEquilSteps;
  const BigReal alchLambda = simParams->alchLambda;
  const BigReal alchLambda2 = simParams->alchLambda2;
  const bool alchEnsembleAvg = simParams->alchEnsembleAvg; 
  const bool FepWhamOn = simParams->alchFepWhamOn;

  if (alchEnsembleAvg && (stepInRun == 0 || stepInRun == alchEquilSteps)) {
    FepNo = 0;
    exp_dE_ByRT = 0.0;
    net_dE = 0.0;
  }
  dE = bondedEnergyDiff_f + electEnergy_f + electEnergySlow_f + ljEnergy_f -
       electEnergy - electEnergySlow - ljEnergy;
  BigReal RT = BOLTZMANN * simParams->alchTemp;

  if (alchEnsembleAvg){
    FepNo++;
    exp_dE_ByRT += exp(-dE/RT);
    net_dE += dE;
  }

  if (simParams->alchOutFreq) { 
    if (stepInRun == 0) {
      if (!fepFile.is_open()) {
        NAMD_backup_file(simParams->alchOutFile);
        fepFile.open(simParams->alchOutFile);
        iout << "OPENING FEP ENERGY OUTPUT FILE\n" << endi;
        if(alchEnsembleAvg){
          fepSum = 0.0;
          fepFile << "#            STEP                 Elec                            "
                  << "vdW                    dE           dE_avg         Temp             dG\n"
                  << "#                           l             l+dl      "
                  << "       l            l+dl         E(l+dl)-E(l)" << std::endl;
        }
        else{
          if(!FepWhamOn){ 
            fepFile << "#            STEP                 Elec                            "
                    << "vdW                    dE         Temp\n"
                    << "#                           l             l+dl      "
                    << "       l            l+dl         E(l+dl)-E(l)" << std::endl;
          } 
        }
      }
      if(!step){
        if(!FepWhamOn){ 
          fepFile << "#NEW FEP WINDOW: "
                  << "LAMBDA SET TO " << alchLambda << " LAMBDA2 " 
                  << alchLambda2 << std::endl;
        }
      }
    }
    if ((alchEquilSteps) && (stepInRun == alchEquilSteps)) {
      fepFile << "#" << alchEquilSteps << " STEPS OF EQUILIBRATION AT "
              << "LAMBDA " << simParams->alchLambda << " COMPLETED\n"
              << "#STARTING COLLECTION OF ENSEMBLE AVERAGE" << std::endl;
    }
    if ((simParams->N) && ((step%simParams->alchOutFreq) == 0)) {
      writeFepEnergyData(step, fepFile);
      fepFile.flush();
    }
    if (alchEnsembleAvg && (step == simParams->N)) {
      fepSum = fepSum + dG;
      fepFile << "#Free energy change for lambda window [ " << alchLambda
              << " " << alchLambda2 << " ] is " << dG 
              << " ; net change until now is " << fepSum << std::endl;
      fepFile.flush();
    }
  }
 }
}

void Controller::outputTiEnergy(int step) {
 if (simParams->alchOn && simParams->alchThermIntOn) {
  const int stepInRun = step - simParams->firstTimestep;
  const int alchEquilSteps = simParams->alchEquilSteps;
  const int alchLambdaFreq = simParams->alchLambdaFreq;

  if (stepInRun == alchEquilSteps) {
    TiNo = 0;
    net_dEdl_bond_1 = 0;
    net_dEdl_bond_2 = 0;
    net_dEdl_elec_1 = 0;
    net_dEdl_elec_2 = 0;
    net_dEdl_lj_1 = 0;
    net_dEdl_lj_2 = 0;
  }
  TiNo++;
  net_dEdl_bond_1 += bondedEnergy_ti_1;
  net_dEdl_bond_2 += bondedEnergy_ti_2;
  net_dEdl_elec_1 += (electEnergy_ti_1 + electEnergySlow_ti_1
                      + electEnergyPME_ti_1);
  net_dEdl_elec_2 += (electEnergy_ti_2 + electEnergySlow_ti_2
                      + electEnergyPME_ti_2);
  net_dEdl_lj_1 += ljEnergy_ti_1;
  net_dEdl_lj_2 += ljEnergy_ti_2;

  // Don't accumulate block averages (you'll get a divide by zero!) or even 
  // open the TI output if alchOutFreq is zero.
  if (simParams->alchOutFreq) {
    if (stepInRun == 0 || (! ((step - 1) % simParams->alchOutFreq))) {
      // output of instantaneous dU/dl now replaced with running average
      // over last alchOutFreq steps (except for step 0)
      recent_TiNo = 0;
      recent_dEdl_bond_1 = 0;
      recent_dEdl_bond_2 = 0;
      recent_dEdl_elec_1 = 0;
      recent_dEdl_elec_2 = 0;
      recent_dEdl_lj_1 = 0;
      recent_dEdl_lj_2 = 0;
      recent_alchWork = 0;
    }
    recent_TiNo++;
    recent_dEdl_bond_1 += bondedEnergy_ti_1;
    recent_dEdl_bond_2 += bondedEnergy_ti_2;
    recent_dEdl_elec_1 += (electEnergy_ti_1 + electEnergySlow_ti_1 
                           + electEnergyPME_ti_1); 
    recent_dEdl_elec_2 += (electEnergy_ti_2 + electEnergySlow_ti_2 
                           + electEnergyPME_ti_2); 
    recent_dEdl_lj_1 += ljEnergy_ti_1;
    recent_dEdl_lj_2 += ljEnergy_ti_2;
    recent_alchWork += alchWork;

    if (stepInRun == 0) {
      if (!tiFile.is_open()) {
        NAMD_backup_file(simParams->alchOutFile);
        tiFile.open(simParams->alchOutFile);
        /* BKR - This has been rather drastically updated to better match 
           stdout. This was necessary for several reasons:
           1) PME global scaling is obsolete (now removed)
           2) scaling of bonded terms was added
           3) alchemical work is now accumulated when switching is active
         */
        iout << "OPENING TI ENERGY OUTPUT FILE\n" << endi;
        tiFile << "#TITITLE:    TS";
        tiFile << FORMAT("BOND1");
        tiFile << FORMAT("AVGBOND1");
        tiFile << FORMAT("ELECT1");
        tiFile << FORMAT("AVGELECT1");
        tiFile << "     ";
        tiFile << FORMAT("VDW1");
        tiFile << FORMAT("AVGVDW1");
        tiFile << FORMAT("BOND2");
        tiFile << FORMAT("AVGBOND2");
        tiFile << FORMAT("ELECT2");
        tiFile << "     ";
        tiFile << FORMAT("AVGELECT2");
        tiFile << FORMAT("VDW2");
        tiFile << FORMAT("AVGVDW2");
        if (alchLambdaFreq > 0) {
          tiFile << FORMAT("ALCHWORK");
          tiFile << FORMAT("CUMALCHWORK");
        }
        tiFile << std::endl;
      }

      if (alchLambdaFreq > 0) {
        tiFile << "#ALCHEMICAL SWITCHING ACTIVE " 
               << simParams->alchLambda << " --> " << simParams->alchLambda2
               << "\n#LAMBDA SCHEDULE: " 
               << "dL: " << simParams->getLambdaDelta() 
               << " Freq: " << alchLambdaFreq;
      }
      else {
        const BigReal alchLambda = simParams->alchLambda;    
        const BigReal bond_lambda_1 = simParams->getBondLambda(alchLambda);
        const BigReal bond_lambda_2 = simParams->getBondLambda(1-alchLambda);
        const BigReal elec_lambda_1 = simParams->getElecLambda(alchLambda);
        const BigReal elec_lambda_2 = simParams->getElecLambda(1-alchLambda);
        const BigReal vdw_lambda_1 = simParams->getVdwLambda(alchLambda);
        const BigReal vdw_lambda_2 = simParams->getVdwLambda(1-alchLambda);
        tiFile << "#NEW TI WINDOW: "
               << "LAMBDA " << alchLambda 
               << "\n#PARTITION 1 SCALING: BOND " << bond_lambda_1
               << " VDW " << vdw_lambda_1 << " ELEC " << elec_lambda_1
               << "\n#PARTITION 2 SCALING: BOND " << bond_lambda_2
               << " VDW " << vdw_lambda_2 << " ELEC " << elec_lambda_2;
      }
      tiFile << "\n#CONSTANT TEMPERATURE: " << simParams->alchTemp << " K"
             << std::endl;
    }

    if ((alchEquilSteps) && (stepInRun == alchEquilSteps)) {
      tiFile << "#" << alchEquilSteps << " STEPS OF EQUILIBRATION AT "
             << "LAMBDA " << simParams->alchLambda << " COMPLETED\n"
             << "#STARTING COLLECTION OF ENSEMBLE AVERAGE" << std::endl;
    }
    if ((step%simParams->alchOutFreq) == 0) {
      writeTiEnergyData(step, tiFile);
      tiFile.flush();
    }
  }
 }
}

/* 
 Work is accumulated whenever alchLambda changes.  In the current scheme we
 always increment lambda _first_, then integrate in time.  Therefore the work 
 is wrt the "old" lambda before the increment.
*/
BigReal Controller::computeAlchWork(const int step) {
  // alchemical scaling factors for groups 1/2 at the previous lambda
  const BigReal oldLambda = simParams->getCurrentLambda(step-1);
  const BigReal bond_lambda_1 = simParams->getBondLambda(oldLambda);
  const BigReal bond_lambda_2 = simParams->getBondLambda(1-oldLambda);
  const BigReal elec_lambda_1 = simParams->getElecLambda(oldLambda);
  const BigReal elec_lambda_2 = simParams->getElecLambda(1-oldLambda);
  const BigReal vdw_lambda_1 = simParams->getVdwLambda(oldLambda);
  const BigReal vdw_lambda_2 = simParams->getVdwLambda(1-oldLambda);
  // alchemical scaling factors for groups 1/2 at the new lambda
  const BigReal alchLambda = simParams->getCurrentLambda(step);
  const BigReal bond_lambda_12 = simParams->getBondLambda(alchLambda);
  const BigReal bond_lambda_22 = simParams->getBondLambda(1-alchLambda);
  const BigReal elec_lambda_12 = simParams->getElecLambda(alchLambda);
  const BigReal elec_lambda_22 = simParams->getElecLambda(1-alchLambda);
  const BigReal vdw_lambda_12 = simParams->getVdwLambda(alchLambda);
  const BigReal vdw_lambda_22 = simParams->getVdwLambda(1-alchLambda); 

  return ((bond_lambda_12 - bond_lambda_1)*bondedEnergy_ti_1 +
          (elec_lambda_12 - elec_lambda_1)*
          (electEnergy_ti_1 + electEnergySlow_ti_1 +  electEnergyPME_ti_1) +
          (vdw_lambda_12 - vdw_lambda_1)*ljEnergy_ti_1 +
          (bond_lambda_22 - bond_lambda_2)*bondedEnergy_ti_2 +
          (elec_lambda_22 - elec_lambda_2)*
          (electEnergy_ti_2 + electEnergySlow_ti_2 + electEnergyPME_ti_2) + 
          (vdw_lambda_22 - vdw_lambda_2)*ljEnergy_ti_2
  );
}

void Controller::writeFepEnergyData(int step, ofstream_namd &file) {
  BigReal eeng = electEnergy + electEnergySlow;
  BigReal eeng_f = electEnergy_f + electEnergySlow_f;
  BigReal dE_Left = eeng_f + ljEnergy_f_left - eeng - ljEnergy;
  BigReal RT = BOLTZMANN * simParams->alchTemp;

  const bool alchEnsembleAvg = simParams->alchEnsembleAvg;
  const int stepInRun = step - simParams->firstTimestep;
  const bool FepWhamOn = simParams->alchFepWhamOn;
  const bool WCARepuOn = simParams->alchFepWCARepuOn;
  const BigReal WCArcut1 = simParams->alchFepWCArcut1;
  const BigReal WCArcut2 = simParams->alchFepWCArcut2;
  const BigReal WCArcut3 = simParams->alchFepWCArcut3;
  const BigReal alchLambda = simParams->alchLambda;
	
  const BigReal alchRepLambda = simParams->alchRepLambda;
  const BigReal alchDispLambda = simParams->alchDispLambda;
  const BigReal alchElecLambda = simParams->alchElecLambda;

  if(stepInRun){
    if(!FepWhamOn){
      fepFile << FEPTITLE(step);
      fepFile << FORMAT(eeng);
      fepFile << FORMAT(eeng_f);
      fepFile << FORMAT(ljEnergy);
      fepFile << FORMAT(ljEnergy_f);
    }
    else{ // FepWhamOn = ON
      if(WCARepuOn){
        if(WCArcut1<WCArcut2) { // [s1,s2]
          fepFile << "FEP_WCA_REP  ";
          fepFile << FORMAT(WCArcut1);
          fepFile << FORMAT(WCArcut2);
          fepFile << FORMAT(1.0);
          fepFile << FORMAT(dE_Left);
        }
        if(WCArcut2<WCArcut3) { // [s2,s3]
          if(WCArcut1<WCArcut2) fepFile << " BREAK ";
          fepFile << "FEP_WCA_REP  ";
          fepFile << FORMAT(WCArcut2);
          fepFile << FORMAT(WCArcut3);
          fepFile << FORMAT(0.0);
          fepFile << FORMAT(dE);
        }
        fepFile << std::endl;
      }
      else if(simParams->alchFepWCADispOn) {
        fepFile << "FEP_WCA_DISP ";
        fepFile << FORMAT(alchDispLambda);
      }
      else if(simParams->alchFepElecOn)	{
        fepFile << "FEP_ELEC     ";
        fepFile << FORMAT(alchElecLambda);
      }
    }
    if( ! WCARepuOn ) {
      fepFile << FORMAT(dE);
    }
    if(alchEnsembleAvg){
      BigReal dE_avg = net_dE/FepNo;
      fepFile << FORMAT(dE_avg);
    }
    if(!FepWhamOn){
      fepFile << FORMAT(temperature);
    }
    if(alchEnsembleAvg){
      dG = -(RT * log(exp_dE_ByRT/FepNo));
      fepFile << FORMAT(dG);
    } 
    if( ! WCARepuOn ) {
      fepFile << std::endl;
    }
  }
}

void Controller::writeTiEnergyData(int step, ofstream_namd &file) {
  tiFile << TITITLE(step);
  tiFile << FORMAT(recent_dEdl_bond_1 / recent_TiNo);
  tiFile << FORMAT(net_dEdl_bond_1 / TiNo);
  tiFile << FORMAT(recent_dEdl_elec_1 / recent_TiNo);
  tiFile << FORMAT(net_dEdl_elec_1/TiNo);
  tiFile << "     ";
  tiFile << FORMAT(recent_dEdl_lj_1 / recent_TiNo);
  tiFile << FORMAT(net_dEdl_lj_1/TiNo);
  tiFile << FORMAT(recent_dEdl_bond_2 / recent_TiNo);
  tiFile << FORMAT(net_dEdl_bond_2 / TiNo);
  tiFile << FORMAT(recent_dEdl_elec_2 / recent_TiNo);
  tiFile << "     ";
  tiFile << FORMAT(net_dEdl_elec_2/TiNo);
  tiFile << FORMAT(recent_dEdl_lj_2 / recent_TiNo);
  tiFile << FORMAT(net_dEdl_lj_2/TiNo);
  if (simParams->alchLambdaFreq > 0) {
    tiFile << FORMAT(recent_alchWork / recent_TiNo);
    tiFile << FORMAT(cumAlchWork);
  }
  tiFile << std::endl;
}
//fepe

void Controller::outputExtendedSystem(int step)
{

  if ( step >= 0 ) {

    // Write out eXtended System Trajectory (XST) file
    if ( simParams->xstFrequency &&
         ((step % simParams->xstFrequency) == 0) )
    {
      if ( ! xstFile.is_open() )
      {
        iout << "OPENING EXTENDED SYSTEM TRAJECTORY FILE\n" << endi;
        NAMD_backup_file(simParams->xstFilename);
        xstFile.open(simParams->xstFilename);
        while (!xstFile) {
          if ( errno == EINTR ) {
            CkPrintf("Warning: Interrupted system call opening XST trajectory file, retrying.\n");
            xstFile.clear();
            xstFile.open(simParams->xstFilename);
            continue;
          }
          char err_msg[257];
          sprintf(err_msg, "Error opening XST trajectory file %s",simParams->xstFilename);
          NAMD_err(err_msg);
        }
        xstFile << "# NAMD extended system trajectory file" << std::endl;
        writeExtendedSystemLabels(xstFile);
      }
      writeExtendedSystemData(step,xstFile);
      xstFile.flush();
    }

    // Write out eXtended System Configuration (XSC) files
    //  Output a restart file
    if ( simParams->restartFrequency &&
         ((step % simParams->restartFrequency) == 0) &&
         (step != simParams->firstTimestep) )
    {
      iout << "WRITING EXTENDED SYSTEM TO RESTART FILE AT STEP "
		<< step << "\n" << endi;
      char fname[140];
      const char *bsuffix = ".old";
      strcpy(fname, simParams->restartFilename);
      if ( simParams->restartSave ) {
        char timestepstr[20];
        sprintf(timestepstr,".%d",step);
        strcat(fname, timestepstr);
        bsuffix = ".BAK";
      }
      strcat(fname, ".xsc");
      NAMD_backup_file(fname,bsuffix);
      ofstream_namd xscFile(fname);
      while (!xscFile) {
        if ( errno == EINTR ) {
          CkPrintf("Warning: Interrupted system call opening XSC restart file, retrying.\n");
          xscFile.clear();
          xscFile.open(fname);
          continue;
        }
        char err_msg[257];
        sprintf(err_msg, "Error opening XSC restart file %s",fname);
        NAMD_err(err_msg);
      } 
      xscFile << "# NAMD extended system configuration restart file" << std::endl;
      writeExtendedSystemLabels(xscFile);
      writeExtendedSystemData(step,xscFile);
      if (!xscFile) {
        char err_msg[257];
        sprintf(err_msg, "Error writing XSC restart file %s",fname);
        NAMD_err(err_msg);
      } 
    }

  }

  //  Output final coordinates
  if (step == FILE_OUTPUT || step == END_OF_RUN)
  {
    int realstep = ( step == FILE_OUTPUT ?
        simParams->firstTimestep : simParams->N );
    iout << "WRITING EXTENDED SYSTEM TO OUTPUT FILE AT STEP "
		<< realstep << "\n" << endi;
    static char fname[140];
    strcpy(fname, simParams->outputFilename);
    strcat(fname, ".xsc");
    NAMD_backup_file(fname);
    ofstream_namd xscFile(fname);
    while (!xscFile) {
      if ( errno == EINTR ) {
        CkPrintf("Warning: Interrupted system call opening XSC output file, retrying.\n");
        xscFile.clear();
        xscFile.open(fname);
        continue;
      }
      char err_msg[257];
      sprintf(err_msg, "Error opening XSC output file %s",fname);
      NAMD_err(err_msg);
    } 
    xscFile << "# NAMD extended system configuration output file" << std::endl;
    writeExtendedSystemLabels(xscFile);
    writeExtendedSystemData(realstep,xscFile);
    if (!xscFile) {
      char err_msg[257];
      sprintf(err_msg, "Error writing XSC output file %s",fname);
      NAMD_err(err_msg);
    } 
  }

  //  Close trajectory file
  if (step == END_OF_RUN) {
    if ( xstFile.is_open() ) {
      xstFile.close();
      iout << "CLOSING EXTENDED SYSTEM TRAJECTORY FILE\n" << endi;
    }
  }

}


void Controller::recvCheckpointReq(const char *key, int task, checkpoint &cp) {  // responding replica
  switch ( task ) {
    case SCRIPT_CHECKPOINT_STORE:
      if ( ! checkpoints.count(key) ) {
        checkpoints[key] = new checkpoint;
      }
      *checkpoints[key] = cp;
      break;
    case SCRIPT_CHECKPOINT_LOAD:
      if ( ! checkpoints.count(key) ) {
        NAMD_die("Unable to load checkpoint, requested key was never stored.");
      }
      cp = *checkpoints[key];
      break;
    case SCRIPT_CHECKPOINT_SWAP:
      if ( ! checkpoints.count(key) ) {
        NAMD_die("Unable to swap checkpoint, requested key was never stored.");
      }
      std::swap(cp,*checkpoints[key]);
      break;
    case SCRIPT_CHECKPOINT_FREE:
      if ( ! checkpoints.count(key) ) {
        NAMD_die("Unable to free checkpoint, requested key was never stored.");
      }
      delete checkpoints[key];
      checkpoints.erase(key);
      break;
  }
}

void Controller::recvCheckpointAck(checkpoint &cp) {  // initiating replica
  if ( checkpoint_task == SCRIPT_CHECKPOINT_LOAD || checkpoint_task == SCRIPT_CHECKPOINT_SWAP ) {
    state->lattice = cp.lattice;
    *(ControllerState*)this = cp.state;
  }
  CkpvAccess(_qd)->process();
}


void Controller::rebalanceLoad(int step)
{
  if ( ! ldbSteps ) { 
    ldbSteps = LdbCoordinator::Object()->getNumStepsToRun();
  }
  if ( ! --ldbSteps ) {
    startBenchTime -= CmiWallTimer();
	Node::Object()->outputPatchComputeMaps("before_ldb", step);
    LdbCoordinator::Object()->rebalance(this);	
	startBenchTime += CmiWallTimer();
    fflush_count = 3;
  }
}

void Controller::cycleBarrier(int doBarrier, int step) {
#if USE_BARRIER
	if (doBarrier) {
	  broadcast->cycleBarrier.publish(step,1);
	  CkPrintf("Cycle time at sync Wall: %f CPU %f\n",
		  CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);
	}
#endif
}

void Controller::traceBarrier(int turnOnTrace, int step) {
	CkPrintf("Cycle time at trace sync (begin) Wall at step %d: %f CPU %f\n", step, CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);	
	CProxy_Node nd(CkpvAccess(BOCclass_group).node);
	nd.traceBarrier(turnOnTrace, step);
	CthSuspend();
}

void Controller::resumeAfterTraceBarrier(int step){
	broadcast->traceBarrier.publish(step,1);
	CkPrintf("Cycle time at trace sync (end) Wall at step %d: %f CPU %f\n", step, CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);	
	awaken();
}

#ifdef MEASURE_NAMD_WITH_PAPI
void Controller::papiMeasureBarrier(int turnOnMeasure, int step){
	CkPrintf("Cycle time at PAPI measurement sync (begin) Wall at step %d: %f CPU %f\n", step, CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);	
	CProxy_Node nd(CkpvAccess(BOCclass_group).node);
	nd.papiMeasureBarrier(turnOnMeasure, step);
	CthSuspend();
}

void Controller::resumeAfterPapiMeasureBarrier(int step){
	broadcast->papiMeasureBarrier.publish(step,1);
	CkPrintf("Cycle time at PAPI measurement sync (end) Wall at step %d: %f CPU %f\n", step, CmiWallTimer()-firstWTime,CmiTimer()-firstCTime);	
	awaken();
}
#endif

void Controller::terminate(void) {
  BackEnd::awaken();
  CthFree(thread);
  CthSuspend();
}

