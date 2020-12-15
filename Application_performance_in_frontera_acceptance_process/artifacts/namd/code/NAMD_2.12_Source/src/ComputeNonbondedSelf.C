/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ComputeNonbondedSelf.h"
#include "ReductionMgr.h"
#include "Patch.h"
#include "LdbCoordinator.h"
#include "Molecule.h"

#include "Node.h"
#include "SimParameters.h"

#define MIN_DEBUG_LEVEL 4
// #define DEBUGM
#include "Debug.h"

ComputeNonbondedSelf::ComputeNonbondedSelf(ComputeID c, PatchID pid,
		ComputeNonbondedWorkArrays* _workArrays,
		int minPartition, int maxPartition, int numPartitions)
  : ComputePatch(c,pid), workArrays(_workArrays),
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

void ComputeNonbondedSelf::initialize() {
  ComputePatch::initialize();
  avgPositionBox = patch->registerAvgPositionPickup(this);
  // BEGIN LA
  velocityBox = patch->registerVelocityPickup(this);
  // END LA

  psiSumBox = patch->registerPsiSumDeposit(this);
  intRadBox = patch->registerIntRadPickup(this);
  bornRadBox = patch->registerBornRadPickup(this);
  dEdaSumBox = patch->registerDEdaSumDeposit(this);
  dHdrPrefixBox = patch->registerDHdrPrefixPickup(this);

#ifdef NAMD_CUDA
  register_cuda_compute_self(cid, patchID);
#endif
}

ComputeNonbondedSelf::~ComputeNonbondedSelf()
{
  delete reduction;
  delete pressureProfileReduction;
  delete [] pressureProfileData;
  if (avgPositionBox != NULL) {
    patch->unregisterAvgPositionPickup(this,&avgPositionBox);
  }
  // BEGIN LA
  if (velocityBox != NULL) {
      patch->unregisterVelocityPickup(this,&velocityBox);
  }
  // END LA

  if (psiSumBox != NULL)
  patch->unregisterPsiSumDeposit(this, &psiSumBox);
  if (intRadBox != NULL)
  patch->unregisterIntRadPickup(this, &intRadBox);
  if (bornRadBox != NULL)
  patch->unregisterBornRadPickup(this, &bornRadBox);
  if (dEdaSumBox != NULL)
  patch->unregisterDEdaSumDeposit(this, &dEdaSumBox);
  if (dHdrPrefixBox != NULL)
  patch->unregisterDHdrPrefixPickup(this, &dHdrPrefixBox);
}

int ComputeNonbondedSelf::noWork() {

  if (patch->flags.doGBIS) {
    gbisPhase = 1 + (gbisPhase % 3);//1->2->3->1...
  }

#ifndef NAMD_CUDA
  if ( patch->flags.doNonbonded && numAtoms ) {
    return 0;  // work to do, enqueue as usual
  } else {
#else
  {
#endif

    if (patch->flags.doGBIS) {
     if (gbisPhase == 1) {
      psiSumBox->skip();
      intRadBox->skip();
      if (patch->flags.doNonbonded) return 1;
      else gbisPhase = 2;
     }
     if (gbisPhase == 2) {
      bornRadBox->skip();
      dEdaSumBox->skip();
      if (patch->flags.doNonbonded) return 1;
      else gbisPhase = 3;
     }
     if (gbisPhase == 3) {
      dHdrPrefixBox->skip();
     }
    }

    // skip all boxes
    positionBox->skip();
    forceBox->skip();
    if ( patch->flags.doMolly ) avgPositionBox->skip();
    // BEGIN LA
    if (patch->flags.doLoweAndersen) velocityBox->skip();
    // END LA

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

void ComputeNonbondedSelf::doForce(CompAtom* p, CompAtomExt* pExt, Results* r)
{
  // Inform load balancer. 
  // I assume no threads will suspend until endWork is called
  //single phase declarations
  CompAtom* v;
  int doEnergy = patch->flags.doEnergy;


/*******************************************************************************
   * Prepare Parameters
 ******************************************************************************/
  if (!patch->flags.doGBIS || gbisPhase == 1) {

#ifdef TRACE_COMPUTE_OBJECTS
  double traceObjStartTime = CmiWallTimer();
#endif

  DebugM(2,"doForce() called.\n");
  DebugM(1,numAtoms << " patch 1 atoms\n");
  DebugM(3, "NUMATOMSxNUMATOMS = " << numAtoms*numAtoms << "\n");

  for ( int i = 0; i < reductionDataSize; ++i ) reductionData[i] = 0;
  if (pressureProfileOn) {
    int n = pressureProfileAtomTypes;
    memset(pressureProfileData, 0, 3*n*n*pressureProfileSlabs*sizeof(BigReal));
    // adjust lattice dimensions to allow constant pressure
    const Lattice &lattice = patch->lattice;
    pressureProfileThickness = lattice.c().z / pressureProfileSlabs;
    pressureProfileMin = lattice.origin().z - 0.5*lattice.c().z;
  }

    plint maxa = (plint)(-1);
    if ( numAtoms > maxa ) {
      char estr[1024];
      sprintf(estr,"patch has %d atoms, maximum allowed is %d",numAtoms,maxa);
      NAMD_die(estr); 
    }

    params.offset = 0.;
    params.offset_f = 0.;
    params.p[0] = p;
    params.p[1] = p;
    params.pExt[0] = pExt;
    params.pExt[1] = pExt;
#ifdef NAMD_KNL
    CompAtomFlt *pFlt = patch->getCompAtomFlt();
    params.pFlt[0] = pFlt;
    params.pFlt[1] = pFlt;
#endif
    params.step = patch->flags.step;
    // BEGIN LA
    params.doLoweAndersen = patch->flags.doLoweAndersen;
    if (params.doLoweAndersen) {
	DebugM(4, "opening velocity box\n");
	v = velocityBox->open();
	params.v[0] = v;
	params.v[1] = v;
    }
    // END LA
#ifndef NAMD_CUDA
    params.ff[0] = r->f[Results::nbond_virial];
    params.ff[1] = r->f[Results::nbond_virial];
#endif
    params.numAtoms[0] = numAtoms;
    params.numAtoms[1] = numAtoms;

    // DMK - Atom Separation (water vs. non-water)
    #if NAMD_SeparateWaters != 0
      params.numWaterAtoms[0] = numWaterAtoms;
      params.numWaterAtoms[1] = numWaterAtoms;
    #endif

    params.reduction = reductionData;
    params.pressureProfileReduction = pressureProfileData;

    params.minPart = minPart;
    params.maxPart = maxPart;
    params.numParts = numParts;

    params.workArrays = workArrays;

    params.pairlists = &pairlists;
    params.savePairlists = 0;
    params.usePairlists = 0;
    if ( patch->flags.savePairlists ) {
      params.savePairlists = 1;
      params.usePairlists = 1;
    } else if ( patch->flags.usePairlists ) {
      if ( ! pairlistsValid ||
           ( 2. * patch->flags.maxAtomMovement > pairlistTolerance ) ) {
        reductionData[pairlistWarningIndex] += 1;
      } else { 
        params.usePairlists = 1;
      }
    }
    if ( ! params.usePairlists ) {
      pairlistsValid = 0;
    }
    params.plcutoff = cutoff;
    params.groupplcutoff = cutoff + 2. * patch->flags.maxGroupRadius;
    if ( params.savePairlists ) {
      pairlistsValid = 1;
      pairlistTolerance = 2. * patch->flags.pairlistTolerance;
      params.plcutoff += pairlistTolerance;
      params.groupplcutoff += pairlistTolerance;
    }


/*******************************************************************************
 * Call Nonbonded Functions
 ******************************************************************************/

    if (numAtoms) {
    if ( patch->flags.doFullElectrostatics )
    {
#ifndef NAMD_CUDA
      params.fullf[0] = r->f[Results::slow_virial];
      params.fullf[1] = r->f[Results::slow_virial];
#endif
      if ( patch->flags.doMolly ) {
        if ( doEnergy ) calcSelfEnergy(&params);
  else calcSelf(&params);
        CompAtom *p_avg = avgPositionBox->open();
        params.p[0] = p_avg;
        params.p[1] = p_avg;
        if ( doEnergy ) calcSlowSelfEnergy(&params);
  else calcSlowSelf(&params);
        avgPositionBox->close(&p_avg);
      } else if ( patch->flags.maxForceMerged == Results::slow ) {
        if ( doEnergy ) calcMergeSelfEnergy(&params);
  else calcMergeSelf(&params);
      } else {
        if ( doEnergy ) calcFullSelfEnergy(&params);
  else calcFullSelf(&params);
      }
    }
    else
      if ( doEnergy ) calcSelfEnergy(&params);
      else calcSelf(&params);
    }//end if atoms
    
    // BEGIN LA
    if (params.doLoweAndersen) {
	DebugM(4, "closing velocity box\n");
	velocityBox->close(&v);
    }
    // END LA
  }//end if not gbis

/*******************************************************************************
 * gbis Loop
*******************************************************************************/
if (patch->flags.doGBIS) {
  SimParameters *simParams = Node::Object()->simParameters;
  gbisParams.sequence = sequence();
  gbisParams.doGBIS = patch->flags.doGBIS;
  gbisParams.numPatches = 1;//self
  gbisParams.gbisPhase = gbisPhase;
  gbisParams.doFullElectrostatics = patch->flags.doFullElectrostatics;
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
  gbisParams.patchID[0] = patch->getPatchID();
  gbisParams.patchID[1] = patch->getPatchID();
  gbisParams.maxGroupRadius = patch->flags.maxGroupRadius;
  gbisParams.doEnergy = doEnergy;
  gbisParams.fsMax = simParams->fsMax;
  for (int i = 0; i < numGBISPairlists; i++)
    gbisParams.gbisStepPairlists[i] = &gbisStepPairlists[i];

  //open boxes
  if (gbisPhase == 1) {
      gbisParams.intRad[0] = intRadBox->open();
      gbisParams.intRad[1] = gbisParams.intRad[0];
      gbisParams.psiSum[0] = psiSumBox->open();
      gbisParams.psiSum[1] = gbisParams.psiSum[0];
      gbisParams.gbInterEnergy=0;
      gbisParams.gbSelfEnergy=0;
  } else if (gbisPhase == 2) {
      gbisParams.bornRad[0] = bornRadBox->open();
      gbisParams.bornRad[1] = gbisParams.bornRad[0];
      gbisParams.dEdaSum[0] = dEdaSumBox->open();
      gbisParams.dEdaSum[1] = gbisParams.dEdaSum[0];
  } else if (gbisPhase == 3) {
      gbisParams.dHdrPrefix[0] = dHdrPrefixBox->open();
      gbisParams.dHdrPrefix[1] = gbisParams.dHdrPrefix[0];
  }

  //make call to calculate GBIS
  if ( !ComputeNonbondedUtil::commOnly ) {
    calcGBIS(&params,&gbisParams);
  }

  //close boxes
  if (gbisPhase == 1) {
      psiSumBox->close(&(gbisParams.psiSum[0]));
  } else if (gbisPhase == 2) {
      dEdaSumBox->close(&(gbisParams.dEdaSum[0]));
  } else if (gbisPhase == 3) {
      bornRadBox->close(&(gbisParams.bornRad[0]));
      reduction->item(REDUCTION_ELECT_ENERGY) += gbisParams.gbInterEnergy;
      reduction->item(REDUCTION_ELECT_ENERGY) += gbisParams.gbSelfEnergy;
      intRadBox->close(&(gbisParams.intRad[0]));
      dHdrPrefixBox->close(&(gbisParams.dHdrPrefix[0]));
  }

}// end if doGBIS


/*******************************************************************************
 * Reduction
*******************************************************************************/
  if (!patch->flags.doGBIS || gbisPhase == 3) {
  submitReductionData(reductionData,reduction);
  if (pressureProfileOn)
    submitPressureProfileData(pressureProfileData, pressureProfileReduction);

  
#ifdef TRACE_COMPUTE_OBJECTS
    traceUserBracketEvent(TRACE_COMPOBJ_IDOFFSET+cid, traceObjStartTime, CmiWallTimer());
#endif

  reduction->submit();
  if (pressureProfileOn)
    pressureProfileReduction->submit();
  }// end not gbis

}

