/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeFullDirect.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "Communicate.h"
#include "Lattice.h"
//#define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "SimParameters.h"

ComputeFullDirect::ComputeFullDirect(ComputeID c) : ComputeHomePatches(c)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  SimParameters *simParams = Node::Object()->simParameters;
  useAvgPositions = 1;
}

ComputeFullDirect::~ComputeFullDirect()
{
  delete reduction;
}

BigReal calc_fulldirect(BigReal *data1, BigReal *results1, int n1,
                        BigReal *data2, BigReal *results2, int n2,
			int selfmode, Lattice *lattice, Tensor &virial)
{
  if ( lattice->a_p() || lattice->b_p() || lattice->c_p() ) {
    #define FULLDIRECT_PERIODIC
    #include "ComputeFullDirectBase.h"
  } else {
    #undef FULLDIRECT_PERIODIC
    #include "ComputeFullDirectBase.h"
  }
}

void ComputeFullDirect::doWork()
{
  int numLocalAtoms;
  BigReal *localData;
  BigReal *localResults;
  BigReal *newLocalResults;
  register BigReal *local_ptr;
  Lattice *lattice;

  int numWorkingPes = (PatchMap::Object())->numNodesWithPatches();

  if ( numWorkingPes > 1 ) NAMD_die("FullDirect not supported for parallel runs.");

  ResizeArrayIter<PatchElem> ap(patchList);

  // Skip computations if nothing to do.
  if ( ! patchList[0].p->flags.doFullElectrostatics )
  {
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      CompAtom *x = (*ap).positionBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).positionBox->close(&x);
      (*ap).forceBox->close(&r);
    }
    reduction->submit();
    return;
  }

  // allocate storage
  numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  localData = new BigReal[4*numLocalAtoms];	// freed at end of this method
  localResults = new BigReal[3*numLocalAtoms];	// freed at end of this method
  newLocalResults = new BigReal[3*numLocalAtoms];  // freed at end of this method

  lattice = &((*(ap.begin())).p->lattice);

  // get positions and charges
  local_ptr = localData;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    CompAtom *x = (*ap).positionBox->open();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      *(local_ptr++) = x[i].position.x;
      *(local_ptr++) = x[i].position.y;
      *(local_ptr++) = x[i].position.z;
      *(local_ptr++) = x[i].charge;
    }

    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
  } 

  // zero out forces
  local_ptr = localResults;
  for(int j=0; j<numLocalAtoms; ++j)
  {
    *(local_ptr++) = 0.;
    *(local_ptr++) = 0.;
    *(local_ptr++) = 0.;
  }

  // perform calculations
  BigReal electEnergy = 0;
  Tensor virial;

#define PEMOD(N) (((N)+numWorkingPes)%numWorkingPes)

  int numStages = numWorkingPes / 2 + 2;
  int lastStage = numStages - 2;
  int sendDataPE = PEMOD(CkMyPe()+1);
  int recvDataPE = PEMOD(CkMyPe()-1);
  int sendResultsPE = PEMOD(CkMyPe()-1);
  int recvResultsPE = PEMOD(CkMyPe()+1);
  int numRemoteAtoms = numLocalAtoms;
  int oldNumRemoteAtoms = 0;
  BigReal *remoteData = 0;
  BigReal *remoteResults = 0;
  register BigReal *remote_ptr;
  register BigReal *end_ptr;

  MOStream *sendDataMsg=CkpvAccess(comm)->
		newOutputStream(sendDataPE, FULLTAG, BUFSIZE);
  MIStream *recvDataMsg=CkpvAccess(comm)->
		newInputStream(recvDataPE, FULLTAG);

  for ( int stage = 0; stage < numStages; ++stage )
  {
    // send remoteResults to sendResultsPE
    if ( stage > 1 )
    {
      DebugM(4,"send remoteResults to sendResultsPE " << sendResultsPE << "\n");
      MOStream *msg=CkpvAccess(comm)->
		newOutputStream(sendResultsPE, FULLFORCETAG, BUFSIZE);
      msg->put(3*oldNumRemoteAtoms,remoteResults);
      delete [] remoteResults;
      msg->end();
      delete msg;
      sendResultsPE = PEMOD(sendResultsPE-1);
    }

    // send remoteData to sendDataPE
    if ( stage < lastStage )
    {
      DebugM(4,"send remoteData to sendDataPE " << sendDataPE << "\n");
      sendDataMsg->put(numRemoteAtoms);
      sendDataMsg->put(4*numRemoteAtoms,(stage?remoteData:localData));
      sendDataMsg->end();
    }

    // allocate new result storage
    if ( stage > 0 && stage <= lastStage )
    {
      DebugM(4,"allocate new result storage\n");
      remoteResults = new BigReal[3*numRemoteAtoms];
      remote_ptr = remoteResults;
      end_ptr = remoteResults + 3*numRemoteAtoms;
      for ( ; remote_ptr != end_ptr; ++remote_ptr ) *remote_ptr = 0.;
    }

    // do calculation
    if ( stage == 0 )
    {  // self interaction
      DebugM(4,"self interaction\n");
      electEnergy += calc_fulldirect(
        localData,localResults,numLocalAtoms,
        localData,localResults,numLocalAtoms,1,lattice,virial);
    }
    else if ( stage < lastStage ||
            ( stage == lastStage && ( numWorkingPes % 2 ) ) )
    {  // full other interaction
      DebugM(4,"full other interaction\n");
      electEnergy += calc_fulldirect(
        localData,localResults,numLocalAtoms,
        remoteData,remoteResults,numRemoteAtoms,0,lattice,virial);
    }
    else if ( stage == lastStage )
    {  // half other interaction
      DebugM(4,"half other interaction\n");
      if ( CkMyPe() < ( numWorkingPes / 2 ) )
        electEnergy += calc_fulldirect(
          localData,localResults,numLocalAtoms/2,
          remoteData,remoteResults,numRemoteAtoms,0,lattice,virial);
      else
        electEnergy += calc_fulldirect(
          localData,localResults,numLocalAtoms,
          remoteData + 4*(numRemoteAtoms/2),
          remoteResults + 3*(numRemoteAtoms/2),
          numRemoteAtoms - (numRemoteAtoms/2), 0,lattice,virial);
    }

    delete [] remoteData;  remoteData = 0;
    oldNumRemoteAtoms = numRemoteAtoms;

    // receive newLocalResults from recvResultsPE
    if ( stage > 1 )
    {
      DebugM(4,"receive newLocalResults from recvResultsPE "
						<< recvResultsPE << "\n");
      MIStream *msg=CkpvAccess(comm)->
		newInputStream(recvResultsPE, FULLFORCETAG);
      msg->get(3*numLocalAtoms,newLocalResults);
      delete msg;
      recvResultsPE = PEMOD(recvResultsPE+1);
      remote_ptr = newLocalResults;
      local_ptr = localResults;
      end_ptr = localResults + 3*numLocalAtoms;
      for ( ; local_ptr != end_ptr; ++local_ptr, ++remote_ptr )
	*local_ptr += *remote_ptr;
    }

    // receive remoteData from recvDataPE
    if ( stage < lastStage )
    {
      DebugM(4,"receive remoteData from recvDataPE "
						<< recvDataPE << "\n");
      recvDataMsg->get(numRemoteAtoms);
      remoteData = new BigReal[4*numRemoteAtoms];
      recvDataMsg->get(4*numRemoteAtoms,remoteData);
    }

  }

  delete sendDataMsg;
  delete recvDataMsg;

  // send out reductions
  DebugM(4,"Full-electrostatics energy: " << electEnergy << "\n");
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += electEnergy;
  reduction->item(REDUCTION_VIRIAL_SLOW_XX) += virial.xx;
  reduction->item(REDUCTION_VIRIAL_SLOW_XY) += virial.xy;
  reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += virial.xz;
  reduction->item(REDUCTION_VIRIAL_SLOW_YX) += virial.yx;
  reduction->item(REDUCTION_VIRIAL_SLOW_YY) += virial.yy;
  reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += virial.yz;
  reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += virial.zx;
  reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += virial.zy;
  reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += virial.zz;
  reduction->submit();

  // add in forces
  local_ptr = localResults;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      f[i].x += *(local_ptr++);
      f[i].y += *(local_ptr++);
      f[i].z += *(local_ptr++);
    }

    (*ap).forceBox->close(&r);
  }

  // free storage
  delete [] localData;		// allocated at beginning of this method
  delete [] localResults;	// allocated at beginning of this method
  delete [] newLocalResults;	// allocated at beginning of this method
}

