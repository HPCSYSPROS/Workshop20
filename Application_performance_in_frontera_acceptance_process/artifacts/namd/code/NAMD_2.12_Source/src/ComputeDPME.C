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
#include "ComputeDPME.h"
#include "ComputeDPMEMsgs.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "SimParameters.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"

#ifdef DPME
#include "dpme2.h"

class ComputeDPMEMaster {
private:
  friend class ComputeDPME;
  ComputeDPME *host;
  ComputeDPMEMaster(ComputeDPME *);
  ~ComputeDPMEMaster();
  void recvData(ComputeDPMEDataMsg *);
  ResizeArray<int> homeNode;
  ResizeArray<int> endForNode;
  int numWorkingPes;
  int numLocalAtoms;
  Pme2Particle *localData;
  SubmitReduction *reduction;
  int runcount;
};

ComputeDPME::ComputeDPME(ComputeID c, ComputeMgr *m) :
  ComputeHomePatches(c), comm(m)
{
  DebugM(4,"ComputeDPME created.\n");
  useAvgPositions = 1;

  int numWorkingPes = (PatchMap::Object())->numNodesWithPatches();

  masterNode = numWorkingPes - 1;
  if ( CkMyPe() == masterNode ) {
    master = new ComputeDPMEMaster(this);
    master->numWorkingPes = numWorkingPes;
  }
  else master = 0;
}

ComputeDPME::~ComputeDPME()
{
  delete master;
}

// These are needed to fix up argument mismatches in DPME.

extern "C" {
  extern int cfftf(int *, double *, double *);
  extern int cfftb(int *, double *, double *);
  extern int cffti1(int *, double *, int *);
  extern int cfftf1(int *n, double *c, double *ch, double *wa, int *ifac);
  extern int cfftb1(int *n, double *c, double *ch, double *wa, int *ifac);
}

int cfftf(int *n, doublecomplex *c, double *wsave) {
  // Casting (doublecomplex*) to (double*) is probably OK.
  return cfftf(n, (double *)c, wsave);
}

int cfftb(int *n, doublecomplex *c, double *wsave) {
  // Casting (doublecomplex*) to (double*) is probably OK.
  return cfftb(n, (double *)c, wsave);
}

int cffti1(int *n, double *wa, double *ifac) {
  // Casting (double*) to (int*) is dangerous if sizes differ!!!
  return cffti1(n, wa, (int *)ifac);
}

int cfftf1(int *n, double *c, double *ch, double *wa, double *ifac) {
  // Casting (double*) to (int*) is dangerous if sizes differ!!!
  return cfftf1(n, c, ch, wa, (int *)ifac);
}

int cfftb1(int *n, double *c, double *ch, double *wa, double *ifac) {
  // Casting (double*) to (int*) is dangerous if sizes differ!!!
  return cfftb1(n, c, ch, wa, (int *)ifac);
}

void ComputeDPME::doWork()
{
  DebugM(4,"Entering ComputeDPME::doWork().\n");

  Pme2Particle *localData;

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
    if ( master ) {
      master->reduction->submit();
    }
    return;
  }

  // allocate storage
  numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  Lattice lattice = patchList[0].p->flags.lattice;

  localData = new Pme2Particle[numLocalAtoms];  // given to message

  // get positions and charges
  Pme2Particle * data_ptr = localData;
  const BigReal coulomb_sqrt = sqrt( COULOMB * ComputeNonbondedUtil::scaling
				* ComputeNonbondedUtil::dielectric_1 );
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    CompAtom *x = (*ap).positionBox->open();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      Vector tmp = lattice.delta(x[i].position);
      data_ptr->x = tmp.x;
      data_ptr->y = tmp.y;
      data_ptr->z = tmp.z;
      data_ptr->cg = coulomb_sqrt * x[i].charge;
      data_ptr->id = x[i].id;
      ++data_ptr;
    }

    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
  }

  // send data to master
  ComputeDPMEDataMsg *msg = new ComputeDPMEDataMsg;
  msg->node = CkMyPe();
  msg->numParticles = numLocalAtoms;
  msg->particles = localData;
  comm->sendComputeDPMEData(msg);
}

void ComputeDPME::recvData(ComputeDPMEDataMsg *msg)
{
  if ( master ) {
    master->recvData(msg);
  }
  else NAMD_die("ComputeDPME::master is NULL!");
}

ComputeDPMEMaster::ComputeDPMEMaster(ComputeDPME *h) :
  host(h), numLocalAtoms(0), runcount(0)
{
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);
  Molecule * molecule = Node::Object()->molecule;
  localData = new Pme2Particle[molecule->numAtoms];
}

ComputeDPMEMaster::~ComputeDPMEMaster()
{
  delete reduction;
  delete [] localData;
}

void ComputeDPMEMaster::recvData(ComputeDPMEDataMsg *msg)
{ 
  DebugM(4,"ComputeDPMEMaster::recvData() " << msg->numParticles
	<< " particles from node " << msg->node << "\n");

  {
    homeNode.add(msg->node);
    Pme2Particle *data_ptr = localData + numLocalAtoms;
    for ( int j = 0; j < msg->numParticles; ++j, ++data_ptr ) {
      *data_ptr = msg->particles[j];
    }
    numLocalAtoms += msg->numParticles;
    endForNode.add(numLocalAtoms);
    delete msg;
  }

  if ( homeNode.size() < numWorkingPes ) return;  // messages outstanding

  DebugM(4,"ComputeDPMEMaster::recvData() running serial code.\n");

  // single processor version

  Lattice lattice = host->getFlags()->lattice;
  SimParameters * simParams = Node::Object()->simParameters;
  int i;

  AtomInfo atom_info;
  atom_info.numatoms = numLocalAtoms;
  atom_info.atompnt = 0;  // not used
  atom_info.freepnt = 0;  // not used
  atom_info.nlocal = numLocalAtoms;
  atom_info.nother = 0;

  if ( ! lattice.orthogonal() ) {
    NAMD_die("DPME only supports orthogonal PBC's.");
  }

  BoxInfo box_info;
  box_info.box.x = lattice.a().x;
  box_info.box.y = lattice.b().y;
  box_info.box.z = lattice.c().z;
  box_info.box.origin = 0.;  // why only one number?
  box_info.prd.x = box_info.box.x;
  box_info.prd.y = box_info.box.y;
  box_info.prd.z = box_info.box.z;
  box_info.prd2.x = 0.5 * box_info.box.x;
  box_info.prd2.y = 0.5 * box_info.box.y;
  box_info.prd2.z = 0.5 * box_info.box.z;
  box_info.cutoff = simParams->cutoff;
  box_info.rs = simParams->cutoff;
  box_info.mc2.x = 2. * ( box_info.prd.x - box_info.cutoff );
  box_info.mc2.y = 2. * ( box_info.prd.y - box_info.cutoff );
  box_info.mc2.z = 2. * ( box_info.prd.z - box_info.cutoff );
  box_info.ewaldcof = ComputeNonbondedUtil::ewaldcof;
  box_info.dtol = simParams->PMETolerance;
  for (i = 0; i <= 8; i++) {
    box_info.recip[i] = 0.; /* assume orthogonal box */
  }
  box_info.recip[0] = 1.0/box_info.box.x;
  box_info.recip[4] = 1.0/box_info.box.y;
  box_info.recip[8] = 1.0/box_info.box.z;

  GridInfo grid_info;
  grid_info.order = simParams->PMEInterpOrder;
  grid_info.nfftgrd.x = simParams->PMEGridSizeX;
  grid_info.nfftgrd.y = simParams->PMEGridSizeY;
  grid_info.nfftgrd.z = simParams->PMEGridSizeZ;
  grid_info.nfft = 0;
  grid_info.volume = lattice.volume();

  PeInfo pe_info;  // hopefully this isn't used for anything
  pe_info.myproc.node = 0;
  pe_info.myproc.nprocs = 1;
  pe_info.myproc.ncube = 0;
  pe_info.inst_node[0] = 0;
  pe_info.igrid = 0;

  PmeVector *localResults;
  double recip_vir[6];
  int time_count = 0;
  int tsteps = 1;
  double mytime = 0.;

  // perform calculations
  BigReal electEnergy = 0;

  // calculate self energy
  Pme2Particle *data_ptr = localData;
  for(i=0; i<numLocalAtoms; ++i)
  {
    electEnergy += data_ptr->cg * data_ptr->cg;
    ++data_ptr;
  }
  electEnergy *= -1. * box_info.ewaldcof / SQRT_PI;

  DebugM(4,"Ewald self energy: " << electEnergy << "\n");

  DebugM(4,"Calling dpme_eval_recip().\n");

  double pme_start_time = 0;
  if ( runcount == 1 ) pme_start_time = CmiTimer();

  electEnergy += dpme_eval_recip( atom_info, localData - 1, &localResults,
			recip_vir, grid_info, box_info, pe_info,
			time_count, tsteps, &mytime );

  if ( runcount == 1 ) {
    iout << iINFO << "PME reciprocal sum CPU time per evaluation: "
         << (CmiTimer() - pme_start_time) << "\n" << endi;
  }

  DebugM(4,"Returned from dpme_eval_recip().\n");

  // REVERSE SIGN OF VIRIAL RETURNED BY DPME
  for(i=0; i<6; ++i) recip_vir[i] *= -1.;

  // send out reductions
  DebugM(4,"Timestep : " << host->getFlags()->step << "\n");
  DebugM(4,"Reciprocal sum energy: " << electEnergy << "\n");
  DebugM(4,"Reciprocal sum virial: " << recip_vir[0] << " " <<
	recip_vir[1] << " " << recip_vir[2] << " " << recip_vir[3] << " " <<
	recip_vir[4] << " " << recip_vir[5] << "\n");
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += electEnergy;
  reduction->item(REDUCTION_VIRIAL_SLOW_XX) += (BigReal)(recip_vir[0]);
  reduction->item(REDUCTION_VIRIAL_SLOW_XY) += (BigReal)(recip_vir[1]);
  reduction->item(REDUCTION_VIRIAL_SLOW_XZ) += (BigReal)(recip_vir[2]);
  reduction->item(REDUCTION_VIRIAL_SLOW_YX) += (BigReal)(recip_vir[1]);
  reduction->item(REDUCTION_VIRIAL_SLOW_YY) += (BigReal)(recip_vir[3]);
  reduction->item(REDUCTION_VIRIAL_SLOW_YZ) += (BigReal)(recip_vir[4]);
  reduction->item(REDUCTION_VIRIAL_SLOW_ZX) += (BigReal)(recip_vir[2]);
  reduction->item(REDUCTION_VIRIAL_SLOW_ZY) += (BigReal)(recip_vir[4]);
  reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += (BigReal)(recip_vir[5]);
  reduction->submit();

  PmeVector *results_ptr = localResults + 1;

  numLocalAtoms = 0;
  for ( i = 0; i < homeNode.size(); ++i ) {
    ComputeDPMEResultsMsg *msg = new ComputeDPMEResultsMsg;
    msg->node = homeNode[i];
    msg->numParticles = endForNode[i] - numLocalAtoms;
    msg->forces = new PmeVector[msg->numParticles];
    for ( int j = 0; j < msg->numParticles; ++j, ++results_ptr ) {
      msg->forces[j] = *results_ptr;
    }
    numLocalAtoms = endForNode[i];
    host->comm->sendComputeDPMEResults(msg,homeNode[i]);
  }

  // reset
  runcount += 1;
  numLocalAtoms = 0;
  homeNode.resize(0);
  endForNode.resize(0);

}

void ComputeDPME::recvResults(ComputeDPMEResultsMsg *msg)
{
  if ( CkMyPe() != msg->node ) {
    NAMD_die("ComputeDPME results sent to wrong node!\n");
    return;
  }
  if ( numLocalAtoms != msg->numParticles ) {
    NAMD_die("ComputeDPME sent wrong number of results!\n");
    return;
  }

  PmeVector *results_ptr = msg->forces;
  ResizeArrayIter<PatchElem> ap(patchList);

  // add in forces
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i)
    {
      f[i].x += results_ptr->x;
      f[i].y += results_ptr->y;
      f[i].z += results_ptr->z;
      ++results_ptr;
    }

    (*ap).forceBox->close(&r);
  }

  delete msg;

}

#endif  // DPME

