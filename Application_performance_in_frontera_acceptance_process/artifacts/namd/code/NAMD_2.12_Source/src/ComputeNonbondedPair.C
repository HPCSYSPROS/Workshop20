/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeNonbondedPair.h"
#include "ReductionMgr.h"
#include "Patch.h"
#include "LdbCoordinator.h"
#include "PatchMap.h"
#include "ComputeMgr.h"
#include "Molecule.h"

#include "Node.h"
#include "SimParameters.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

ComputeNonbondedPair::ComputeNonbondedPair(ComputeID c, PatchID pid[], int trans[],
		ComputeNonbondedWorkArrays* _workArrays,
		int minPartition, int maxPartition, int numPartitions)
  : ComputePatchPair(c,pid,trans), workArrays(_workArrays),
    minPart(minPartition), maxPart(maxPartition), numParts(numPartitions)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  if (pressureProfileOn) {
    int n = pressureProfileAtomTypes; 
    pressureProfileData = new BigReal[3*n*n*pressureProfileSlabs];
    pressureProfileReduction = ReductionMgr::Object()->willSubmit(
	REDUCTIONS_PPROF_NONBONDED, 3*pressureProfileSlabs*((n*(n+1))/2));
  } else {
    pressureProfileReduction = NULL;
    pressureProfileData = NULL;
  }
  pairlistsValid = 0;
  pairlistTolerance = 0.;
  params.simParameters = Node::Object()->simParameters;
  params.parameters = Node::Object()->parameters;
  params.random = Node::Object()->rand;
}

void ComputeNonbondedPair::initialize() {
  ComputePatchPair::initialize();
  for (int i=0; i<2; i++) {
    avgPositionBox[i] = patch[i]->registerAvgPositionPickup(this);
    // BEGIN LA
    velocityBox[i] = patch[i]->registerVelocityPickup(this);
    // END LA

    psiSumBox[i] = patch[i]->registerPsiSumDeposit(this);
    intRadBox[i] = patch[i]->registerIntRadPickup(this);
    bornRadBox[i] = patch[i]->registerBornRadPickup(this);
    dEdaSumBox[i] = patch[i]->registerDEdaSumDeposit(this);
    dHdrPrefixBox[i] = patch[i]->registerDHdrPrefixPickup(this);
  }
#ifdef NAMD_CUDA
  register_cuda_compute_pair(cid, patchID, trans);
#endif
}

ComputeNonbondedPair::~ComputeNonbondedPair()
{
  delete reduction;
  delete pressureProfileReduction;
  delete [] pressureProfileData;
  for (int i=0; i<2; i++) {
    if (avgPositionBox[i] != NULL) {
      patch[i]->unregisterAvgPositionPickup(this,&avgPositionBox[i]);
    }
    // BEGIN LA
    if (velocityBox[i] != NULL) {
      patch[i]->unregisterVelocityPickup(this,&velocityBox[i]);
    }
    // END LA

    if (psiSumBox[i] != NULL) {
      patch[i]->unregisterPsiSumDeposit(this,&psiSumBox[i]);
    }
    if (intRadBox[i] != NULL) {
      patch[i]->unregisterIntRadPickup(this,&intRadBox[i]);
    }
    if (bornRadBox[i] != NULL) {
      patch[i]->unregisterBornRadPickup(this,&bornRadBox[i]);
    }
    if (dEdaSumBox[i] != NULL) {
      patch[i]->unregisterDEdaSumDeposit(this,&dEdaSumBox[i]);
    }
    if (dHdrPrefixBox[i] != NULL) {
      patch[i]->unregisterDHdrPrefixPickup(this,&dHdrPrefixBox[i]);
    }
  }
}

int ComputeNonbondedPair::noWork() {

  if (patch[0]->flags.doGBIS) {
    gbisPhase = 1 + (gbisPhase % 3);//1->2->3->1...
  }

#ifndef NAMD_CUDA
  if ( patch[0]->flags.doNonbonded && (numAtoms[0] && numAtoms[1]) ) {
    return 0;  // work to do, enqueue as usual
  } else {
#else
  {
#endif

    if (patch[0]->flags.doGBIS) {
     if (gbisPhase == 1) {
      for (int i=0; i<2; i++) {
        psiSumBox[i]->skip();
        intRadBox[i]->skip();
      }
      if (patch[0]->flags.doNonbonded) return 1;
      else gbisPhase = 2;
     }
     if (gbisPhase == 2) {
      for (int i=0; i<2; i++) {
        bornRadBox[i]->skip();
        dEdaSumBox[i]->skip();
      }
      if (patch[0]->flags.doNonbonded) return 1;
      else gbisPhase = 3;
     }
     if (gbisPhase == 3) {
      for (int i=0; i<2; i++) {
        dHdrPrefixBox[i]->skip();
      }
     }
    }

    // skip all boxes
    for (int i=0; i<2; i++) {
      positionBox[i]->skip();
      forceBox[i]->skip();
      if ( patch[0]->flags.doMolly ) avgPositionBox[i]->skip();
      // BEGIN LA
      if (patch[0]->flags.doLoweAndersen) velocityBox[i]->skip();
      // END LA
    }

    reduction->item(REDUCTION_COMPUTE_CHECKSUM) += 1.;
    reduction->submit();
    if (pressureProfileOn) 
      pressureProfileReduction->submit();

#ifndef NAMD_CUDA
    // Inform load balancer
    LdbCoordinator::Object()->skipWork(ldObjHandle);
#endif

    return 1;  // no work to do, do not enqueue
  }
}

void ComputeNonbondedPair::doForce(CompAtom* p[2], CompAtomExt* pExt[2], Results* r[2])
{
  // Inform load balancer. 
  // I assume no threads will suspend until endWork is called

  //single phase declarations
  int doEnergy = patch[0]->flags.doEnergy;
  int a = 0;  int b = 1;
  // swap to place more atoms in inner loop (second patch)
  if ( numAtoms[0] > numAtoms[1] ) { a = 1; b = 0; }
  CompAtom* v[2];


/*******************************************************************************
 * Prepare Parameters
*******************************************************************************/
  if (!patch[0]->flags.doGBIS || gbisPhase == 1) {

#ifdef TRACE_COMPUTE_OBJECTS
    double traceObjStartTime = CmiWallTimer();
#endif

  DebugM(2,"doForce() called.\n");
  DebugM(2, numAtoms[0] << " patch #1 atoms and " <<
	numAtoms[1] << " patch #2 atoms\n");


  for ( int i = 0; i < reductionDataSize; ++i )
    reductionData[i] = 0;
  if (pressureProfileOn) {
    int n = pressureProfileAtomTypes; 
    memset(pressureProfileData, 0, 3*n*n*pressureProfileSlabs*sizeof(BigReal));
    // adjust lattice dimensions to allow constant pressure
    const Lattice &lattice = patch[0]->lattice;
    pressureProfileThickness = lattice.c().z / pressureProfileSlabs;
    pressureProfileMin = lattice.origin().z - 0.5*lattice.c().z;
  }

    params.reduction = reductionData;
    params.pressureProfileReduction = pressureProfileData;

    params.minPart = minPart;
    params.maxPart = maxPart;
    params.numParts = numParts;

    params.workArrays = workArrays;

    params.pairlists = &pairlists;
    params.savePairlists = 0;
    params.usePairlists = 0;
    if ( patch[0]->flags.savePairlists ) {
      params.savePairlists = 1;
      params.usePairlists = 1;
    } else if ( patch[0]->flags.usePairlists && patch[1]->flags.usePairlists ) {
      if ( ! pairlistsValid ||
           ( patch[0]->flags.maxAtomMovement +
             patch[1]->flags.maxAtomMovement > pairlistTolerance ) ) {
        reductionData[pairlistWarningIndex] += 1;
      } else {
        params.usePairlists = 1;
      }
    }
    if ( ! params.usePairlists ) {
      pairlistsValid = 0;
    }
    params.plcutoff = cutoff;
    params.groupplcutoff = cutoff +
	patch[0]->flags.maxGroupRadius + patch[1]->flags.maxGroupRadius;
    if ( params.savePairlists ) {
      pairlistsValid = 1;
      pairlistTolerance = patch[0]->flags.pairlistTolerance +
                          patch[1]->flags.pairlistTolerance;
      params.plcutoff += pairlistTolerance;
      params.groupplcutoff += pairlistTolerance;
    }


    const Lattice &lattice = patch[0]->lattice;
    params.offset = lattice.offset(trans[a]) - lattice.offset(trans[b]);

    PatchMap* patchMap = PatchMap::Object();
    params.offset_f = params.offset + lattice.unscale(patchMap->center(patchID[a]))
                                    - lattice.unscale(patchMap->center(patchID[b]));

    // Atom Sorting : If we are sorting the atoms along the line connecting
    //   the patch centers, then calculate a normalized vector pointing from
    //   patch a to patch b (i.e. outer loop patch to inner loop patch).
    #if NAMD_ComputeNonbonded_SortAtoms != 0

      params.projLineVec = params.offset_f * ( -1. / params.offset_f.length() );

    #endif

      params.p[0] = p[a];
      params.p[1] = p[b];
      params.pExt[0] = pExt[a]; 
      params.pExt[1] = pExt[b];
#ifdef NAMD_KNL
      params.pFlt[0] = patch[a]->getCompAtomFlt();
      params.pFlt[1] = patch[b]->getCompAtomFlt();
#endif
      // BEGIN LA
      params.doLoweAndersen = patch[0]->flags.doLoweAndersen;
      if (params.doLoweAndersen) {
	  DebugM(4, "opening velocity boxes\n");
	  v[0] = velocityBox[0]->open();
	  v[1] = velocityBox[1]->open();
	  params.v[0] = v[a];
	  params.v[1] = v[b];
      }
      // END LA
#ifndef NAMD_CUDA
      params.ff[0] = r[a]->f[Results::nbond_virial];
      params.ff[1] = r[b]->f[Results::nbond_virial];
#endif
      params.numAtoms[0] = numAtoms[a];
      params.numAtoms[1] = numAtoms[b];
      params.step = patch[0]->flags.step;

      // DMK - Atom Separation (water vs. non-water)
      #if NAMD_SeparateWaters != 0
        params.numWaterAtoms[0] = numWaterAtoms[a];
        params.numWaterAtoms[1] = numWaterAtoms[b];
      #endif


/*******************************************************************************
 * Call Nonbonded Functions
*******************************************************************************/
      if (numAtoms[0] && numAtoms[1]) {//only do if has atoms since gbis noWork doesn't account for no atoms

      //force calculation calls
      if ( patch[0]->flags.doFullElectrostatics )
      {
#ifndef NAMD_CUDA
	params.fullf[0] = r[a]->f[Results::slow_virial];
	params.fullf[1] = r[b]->f[Results::slow_virial];
#endif
	if ( patch[0]->flags.doMolly ) {
          if ( doEnergy )
            calcPairEnergy(&params);
	  else calcPair(&params);
	  CompAtom *p_avg[2];
	  p_avg[0] = avgPositionBox[0]->open();
	  p_avg[1] = avgPositionBox[1]->open();
	  params.p[0] = p_avg[a];
	  params.p[1] = p_avg[b];
	  if ( doEnergy ) calcSlowPairEnergy(&params);
	  else calcSlowPair(&params);
	  avgPositionBox[0]->close(&p_avg[0]);
	  avgPositionBox[1]->close(&p_avg[1]);
        } else if ( patch[0]->flags.maxForceMerged == Results::slow ) {
          if ( doEnergy ) calcMergePairEnergy(&params);
    else calcMergePair(&params);
  } else {
    if ( doEnergy ) calcFullPairEnergy(&params);
    else calcFullPair(&params);
  }
      }
      else
        if ( doEnergy ) calcPairEnergy(&params);
        else calcPair(&params);

      }//end if has atoms
      
      // BEGIN LA
      if (params.doLoweAndersen) {
	  DebugM(4, "closing velocity boxes\n");
	  velocityBox[0]->close(&v[0]);
	  velocityBox[1]->close(&v[1]);
      }
      // END LA
    }// end not gbis

/*******************************************************************************
 * gbis Loop
*******************************************************************************/
if (patch[0]->flags.doGBIS) {
  SimParameters *simParams = Node::Object()->simParameters;
  gbisParams.sequence = sequence();
  gbisParams.doGBIS = patch[0]->flags.doGBIS;
  gbisParams.numPatches = 2;//pair
  gbisParams.gbisPhase = gbisPhase;
  gbisParams.doFullElectrostatics = patch[0]->flags.doFullElectrostatics;
  gbisParams.epsilon_s = simParams->solvent_dielectric;
  gbisParams.epsilon_p = simParams->dielectric;
  gbisParams.rho_0 = simParams->coulomb_radius_offset;
  gbisParams.kappa = simParams->kappa;
  gbisParams.cutoff = simParams->cutoff;
  gbisParams.doSmoothing = simParams->switchingActive;
  gbisParams.a_cut = simParams->alpha_cutoff;
  gbisParams.delta = simParams->gbis_delta;
  gbisParams.beta = simParams->gbis_beta;
  gbisParams.gamma = simParams->gbis_gamma;
  gbisParams.alpha_max = simParams->alpha_max;
  gbisParams.cid = cid;
  gbisParams.patchID[0] = patch[a]->getPatchID();
  gbisParams.patchID[1] = patch[b]->getPatchID();
  gbisParams.maxGroupRadius = patch[0]->flags.maxGroupRadius;
  if (patch[1]->flags.maxGroupRadius > gbisParams.maxGroupRadius)
    gbisParams.maxGroupRadius = patch[1]->flags.maxGroupRadius;
  gbisParams.doEnergy = doEnergy;
  gbisParams.fsMax = simParams->fsMax;
  for (int i = 0; i < numGBISPairlists; i++)
    gbisParams.gbisStepPairlists[i] = &gbisStepPairlists[i];

  //open boxes
  if (gbisPhase == 1) {
      gbisParams.intRad[0] = intRadBox[a]->open();
      gbisParams.intRad[1] = intRadBox[b]->open();
      gbisParams.psiSum[0] = psiSumBox[a]->open();
      gbisParams.psiSum[1] = psiSumBox[b]->open();
      gbisParams.gbInterEnergy=0;
      gbisParams.gbSelfEnergy=0;

  } else if (gbisPhase == 2) {
      gbisParams.bornRad[0] = bornRadBox[a]->open();
      gbisParams.bornRad[1] = bornRadBox[b]->open();
      gbisParams.dEdaSum[0] = dEdaSumBox[a]->open();
      gbisParams.dEdaSum[1] = dEdaSumBox[b]->open();
  } else if (gbisPhase == 3) {
      gbisParams.dHdrPrefix[0] = dHdrPrefixBox[a]->open();
      gbisParams.dHdrPrefix[1] = dHdrPrefixBox[b]->open();
  }

  //make call to calculate GBIS
  if ( !ComputeNonbondedUtil::commOnly ) {
    calcGBIS(&params,&gbisParams);
  }

  //close boxes
  if (gbisPhase == 1) {
      psiSumBox[0]->close(&(gbisParams.psiSum[a]));
      psiSumBox[1]->close(&(gbisParams.psiSum[b]));
  } else if (gbisPhase == 2) {
      dEdaSumBox[0]->close(&(gbisParams.dEdaSum[a]));
      dEdaSumBox[1]->close(&(gbisParams.dEdaSum[b]));


  } else if (gbisPhase == 3) {
      bornRadBox[0]->close(&(gbisParams.bornRad[a]));
      bornRadBox[1]->close(&(gbisParams.bornRad[b]));
      reduction->item(REDUCTION_ELECT_ENERGY) += gbisParams.gbInterEnergy;
      reduction->item(REDUCTION_ELECT_ENERGY) += gbisParams.gbSelfEnergy;
      intRadBox[0]->close(&(gbisParams.intRad[a]));
      intRadBox[1]->close(&(gbisParams.intRad[b]));
      dHdrPrefixBox[0]->close(&(gbisParams.dHdrPrefix[a]));
      dHdrPrefixBox[1]->close(&(gbisParams.dHdrPrefix[b]));
  }

}//end if doGBIS


  if (!patch[0]->flags.doGBIS || gbisPhase == 3) {
  submitReductionData(reductionData,reduction);
  if (pressureProfileOn)
    submitPressureProfileData(pressureProfileData, pressureProfileReduction);

#ifdef TRACE_COMPUTE_OBJECTS
    traceUserBracketEvent(TRACE_COMPOBJ_IDOFFSET+cid, traceObjStartTime, CmiWallTimer());
#endif


  reduction->submit();
  if (pressureProfileOn)
    pressureProfileReduction->submit();
  }//end gbis end phase

}//end do Force

