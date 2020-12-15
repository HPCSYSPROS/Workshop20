/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "common.h"
#include "InfoStream.h"
#include "Node.h"
#include "SimParameters.h"
#include "ComputeNonbondedUtil.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeDPMTA.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "Lattice.h"
#include "Communicate.h"
#include "InfoStream.h"
#include "ProcessorPrivate.h"

#define MIN_DEBUG_LEVEL 2
// #define DEBUGM
#include "Debug.h"

#ifdef DPMTA

#include <pvmc.h>

// #define DUMP_DPMTA

void ComputeDPMTA::get_FMA_cube(int resize)
{
  Vector boxSize,boxCenter;	// used to see if things change
  PatchMap *patchMap = PatchMap::Object();

  if (usePBC == FALSE)
  {
    //  From these extremes, figure out how many patches we will
    //  have to have in each direction
    SimParameters *simParams = Node::Object()->simParameters;
    int dim_x = patchMap->gridsize_a();
    int dim_y = patchMap->gridsize_b();
    int dim_z = patchMap->gridsize_c();

    boxSize.x = dim_x*simParams->patchDimension;
    boxSize.y = dim_y*simParams->patchDimension;
    boxSize.z = dim_z*simParams->patchDimension;
    BigReal skirt = 2*simParams->patchDimension;

    boxCenter = patchMap->origin();
    boxCenter.x += boxSize.x/2.0;
    boxCenter.y += boxSize.y/2.0;
    boxCenter.z += boxSize.z/2.0;

    //  add the skirt of empty patches by adding 2 patches in every direction
    boxSize.x += skirt;
    boxSize.y += skirt;
    boxSize.z += skirt;
  }
  else
  {
    DebugM(2,"getting patch info for FMA box\n");

    // determine boxSize from the PBC lattice
    // lattice is the same on all patches, so choose first patch
    ResizeArrayIter<PatchElem> ap(patchList);
    DebugM(2,"getting first patch info for FMA box\n");
    ap = ap.begin();
    DebugM(2,"getting lattice from patch for FMA box\n");
    Lattice lattice = (*ap).p->lattice;
    if ( ! lattice.orthogonal() ) {
      NAMD_die("DPMTA (FMA) only supports orthogonal PBC's.");
    }
    DebugM(2,"getting patch dimension for FMA box\n");
    boxSize.x = lattice.a().x;
    boxSize.y = lattice.b().y;
    boxSize.z = lattice.c().z;
    DebugM(2,"boxSize is " << boxSize << "\n");
    boxCenter = lattice.origin();
  }

  // don't bother checking if the center has moved since it depends on the size.
  if (boxsize != boxSize)
  {
	DebugM(2,"resetting FMA box\n");
	// reset the size and center
	boxsize = boxSize;
	boxcenter = boxCenter;

	// reset DPMTA (only reset it after it has been initialized!)
	if (resize && usePBC)
	{
	  PmtaVector center,v1,v2,v3;
	  center.x = boxcenter.x;
	  center.y = boxcenter.y;
	  center.z = boxcenter.z;
	  v1.x = boxsize.x;
	  v2.y = boxsize.y;
	  v3.z = boxsize.z;
	  iout << iINFO << "DPMTA box resized:\n";
	  iout << iINFO << "BOX DIMENSIONS = (" << v1.x << ","
		<< v2.y << "," << v3.z << ")\n";
	  iout << iINFO << "BOX CENTER = (" << center.x << ","
		<< center.y << "," << center.z << ")\n";
	  iout << endi;
	  DebugM(2,"calling PMTAresize()\n");
	  PMTAresize(&v1,&v2,&v3,&center);
	  DebugM(2,"called PMTAresize()\n");
	}
  }
  DebugM(2,"cube center: " << boxcenter << " size=" << boxsize << "\n");
}

ComputeDPMTA::ComputeDPMTA(ComputeID c) : ComputeHomePatches(c)
{
  useAvgPositions = 1;
}

void ComputeDPMTA::initialize()
{
  ComputeHomePatches::initialize();

  DebugM(2,"ComputeDPMTA creating\n");
  // comm should always be initialized by this point...
  // In the (bug) case that it isn't, then initialize it.
  if (CkpvAccess(comm) == NULL)
  {
    NAMD_die("Communication protocol (Converse, PVM, etc.) not initialized.");
  }

  // **** NOTE: node 0 must initialized before any other nodes register.

  //  Set everything to 0
  totalAtoms = 0;
  fmaResults = NULL;
  ljResults = NULL;
  boxcenter = 1;	// reset the array (no divide by zero)
  boxsize = 1;	// reset the array (no divide by zero)
  usePBC = FALSE;	// assume not...

  // all nodes should init
  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

  // Don't need any more initialization  -JCP
  ResizeArrayIter<PatchElem> ap(patchList);
  DebugM(2,"init() getting first patch info for FMA box\n");
  ap = ap.begin();
  DebugM(2,"init() getting lattice from patch for FMA box\n");
  initLattice.x = (*ap).p->lattice.a().x;
  initLattice.y = (*ap).p->lattice.b().y;
  initLattice.z = (*ap).p->lattice.c().z;
  DebugM(2,"init() initLattice is " << initLattice << "\n");

  //  NOTE that the theta value is hardwired to the value of 0.715
  //  as per the recommendation of the Duke developers

  //  NOTE 2: Theta is now an optional config parameter,
  //  but it defaults to 0.715

  // check for PBC
  usePBC = ( patchMap->periodic_a() ? 1 : 0 )
	 + ( patchMap->periodic_b() ? 1 : 0 )
	 + ( patchMap->periodic_c() ? 1 : 0 );
  if ((usePBC != 0) && (usePBC != 3))
  {
    NAMD_die("DPMTA (FMA) does not support 1 or 2 dimensional PBC's.");
  }
  DebugM(2,"Use PBC = " << usePBC << "\n");
  usePBC = (usePBC == 3);	// either PBC "3D" or no PBC
  if ( usePBC ) {
    iout << iWARN << "DPMTA (FMA) pressure tensor is incorrect.\n"
	<< iWARN << "Do not use DPMTA with anisotropic pressure control.\n"
	<< endi;
  }

  //  Get the size of the FMA cube
  DebugM(2,"DPMTA getting FMA cube\n");
  get_FMA_cube(FALSE);
  DebugM(2,"DPMTA got FMA cube\n");

  if (CkMyPe() != 0)
  {
    DebugM(2,"waiting for Init go-ahead\n");
    MIStream *msg2 = CkpvAccess(comm)->newInputStream(ANY, DPMTATAG);
    int dummy;
    msg2->get(dummy);
    delete msg2;
    slavetids=NULL;
    if (PMTAregister() < 0)
    {
	NAMD_die("PMTARegister failed!!");
    }
    DebugM(2,"DPMTA done PMTAinit.\n");
    return;
  }
  DebugM(2,"DPMTA configuring\n");

  // *****************************************
  // ONLY THE MASTER (NODE 0) NEEDS TO DO THIS:

  int numProcs = (PatchMap::Object())->numNodesWithPatches();

  slavetids = new int[numProcs];
  if (slavetids == NULL)
  {
    NAMD_die("Memory allocation failed in FMAInterface::FMAInterface");
  }

  // pvm_spawn is a dummy function under Converse.  Just the array is required.
  pvm_spawn(NULL,NULL,0,NULL,numProcs,slavetids);
  DebugM(2,"DPMTA slavetids allocated\n");

  // reduce function calling time
  SimParameters *simParams = Node::Object()->simParameters;

  //  initialize DPMTA
  PmtaInitData pmta_data;
  memset(&pmta_data,0,sizeof(pmta_data));
  pmta_data.nprocs = numProcs;
  pmta_data.nlevels = simParams->FMALevels;
  pmta_data.mp = simParams->FMAMp;
  pmta_data.mp_lj = 4;
  pmta_data.fft = simParams->FMAFFTOn;
  pmta_data.fftblock = simParams->FMAFFTBlock;
  pmta_data.pbc = usePBC;	// use Periodic boundary condition
  pmta_data.kterm = 0;
  pmta_data.theta = simParams->fmaTheta;
  pmta_data.v1.x = boxsize.x;
  pmta_data.v1.y = 0.;
  pmta_data.v1.z = 0.;
  pmta_data.v2.x = 0.;
  pmta_data.v2.y = boxsize.y;
  pmta_data.v2.z = 0.;
  pmta_data.v3.x = 0.;
  pmta_data.v3.y = 0.;
  pmta_data.v3.z = boxsize.z;
  pmta_data.cellctr.x = boxcenter.x;
  pmta_data.cellctr.y = boxcenter.y;
  pmta_data.cellctr.z = boxcenter.z;
  pmta_data.calling_num = pmta_data.nprocs;
  pmta_data.calling_tids = slavetids;

  iout << iINFO << "DPMTA parameters are:\n";
  iout << iINFO << "  LEVELS = " << pmta_data.nlevels << "\n";
  iout << iINFO << "  NUMBER OF MULTIPOLE TERMS = " << pmta_data.mp << "\n";
  iout << iINFO << "  FFT FLAG = " << pmta_data.fft << "\n";
  iout << iINFO << "  FFT BLOCKING FACTOR = " << pmta_data.fftblock << "\n";
  if ( usePBC ) iout << iINFO << "  SYSTEM IS PERIODIC\n" << endi;
  iout << iINFO << "  BOX DIMENSIONS = (" << pmta_data.v1.x << ","
	<< pmta_data.v2.y << "," << pmta_data.v3.z << ")\n";
  iout << iINFO << "  BOX CENTER = (" << pmta_data.cellctr.x << ","
	<< pmta_data.cellctr.y << "," << pmta_data.cellctr.z << ")\n";
  iout << endi;

  if ( usePBC )
  {
    pmta_data.cellctr.x = 0.;
    pmta_data.cellctr.y = 0.;
    pmta_data.cellctr.z = 0.;
  }

#ifdef DUMP_DPMTA
  FILE *fp;
  fp = fopen("DUMP_DPMTA.init","w");
  fwrite(&pmta_data,sizeof(PmtaInitData),1,fp);
  fclose(fp);
#endif

  DebugM(2,"DPMTA calling PMTAinit.\n");
  if (PMTAinit(&pmta_data,slavetids) >= 0)
  {
	iout << iINFO << "SUCCESSFULLY STARTED DPMTA\n" << endi;
  }
  else
  {
	NAMD_die("Unable to start DPMTA!");
  }

  // tell all nodes that it is OK to register
  MOStream *msg = CkpvAccess(comm)->newOutputStream(ALL, DPMTATAG, BUFSIZE);
  // don't actually put in data...  Nodes just need it as a flag.
  msg->put(TRUE);
  msg->end();
  delete msg;
  DebugM(2,"Init go-ahead\n");
  MIStream *msg1 = CkpvAccess(comm)->newInputStream(ANY, DPMTATAG);
  int dummy1;
  msg1->get(dummy1);
  delete msg1;
  DebugM(2,"got Init go-ahead\n");

  //  Register this master with the other DPMTA processes
  if (PMTAregister() < 0)
  {
	NAMD_die("PMTARegister failed!!");
  }
  DebugM(2,"DPMTA done PMTAinit.\n");
  DebugM(2,"DPMTA configured\n");
}

ComputeDPMTA::~ComputeDPMTA()
{
  DebugM(2,"DPMTA exiting\n");
  //  If this is the master node, then call PMTAexit()
  if (CkMyPe() == 0)	PMTAexit();

  if (fmaResults)
	{
	free(fmaResults);
	fmaResults = NULL;
	}
  delete [] ljResults;
  delete [] slavetids;
  DebugM(2,"DPMTA exited\n");

  delete reduction;
}


void ComputeDPMTA::doWork()
{
  ResizeArrayIter<PatchElem> ap(patchList);
  PmtaParticle *particle_list = NULL;

  // 0. only run when necessary
  // Skip computations if nothing to do.
  if (!patchList[0].p->flags.doFullElectrostatics)
  {
    for (ap = ap.begin(); ap != ap.end(); ap++) {
      CompAtom *x = (*ap).positionBox->open();
      Results *r = (*ap).forceBox->open();
      (*ap).forceBox->close(&r);
      (*ap).positionBox->close(&x);
    }
    reduction->submit();
    return;
  }

  // setup
  // 1. get totalAtoms
  for (totalAtoms=0, ap = ap.begin(); ap != ap.end(); ap++)
     totalAtoms += (*ap).p->getNumAtoms();

  Vector newLattice;
  Vector rescaleFactor;
  if (usePBC)
    {
    ap = ap.begin();
    Lattice lattice = (*ap).p->lattice;
    if ( ! lattice.orthogonal() ) {
      NAMD_die("DPMTA (FMA) only supports orthogonal PBC's.");
    }
    newLattice.x = lattice.a().x;
    newLattice.y = lattice.b().y;
    newLattice.z = lattice.c().z;
    rescaleFactor.x = initLattice.x / newLattice.x;
    rescaleFactor.y = initLattice.y / newLattice.y;
    rescaleFactor.z = initLattice.z / newLattice.z;
    DebugM(2,"Rescale factor = " << initLattice << "/" << newLattice
		<< " = " << rescaleFactor << "\n");
    DebugM(2,"boxcenter = " << boxcenter << "\n");
    }
  else
    {
    rescaleFactor.x = 1;
    rescaleFactor.y = 1;
    rescaleFactor.z = 1;
    }

  // 2. setup atom list
  int i,j;
  particle_list = (PmtaParticle *)calloc(totalAtoms,sizeof(PmtaParticle));
  fmaResults =    (PmtaPartInfo *)calloc(totalAtoms,sizeof(PmtaPartInfo));
  if (!particle_list || !fmaResults)
	{
	NAMD_die("DPMTA Failed to allocate memory.");
	}

  BigReal unitFactor = sqrt(COULOMB * ComputeNonbondedUtil::scaling
				* ComputeNonbondedUtil::dielectric_1);
  DebugM(2,"Charge unit factor = " << unitFactor << "\n");
  for (i=0, ap = ap.begin(); ap != ap.end(); ap++)
  {
    CompAtom *x = (*ap).positionBox->open();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }

    // store each atom in the particle_list
    Vector pos;
    for(j=0; j<(*ap).p->getNumAtoms(); j++)
    {
      // explicitly copy -- two different data structures
      if (usePBC)
	{
	particle_list[i].p.x = rescaleFactor.x * (x[j].position.x-boxcenter.x);
	particle_list[i].p.y = rescaleFactor.y * (x[j].position.y-boxcenter.y);
	particle_list[i].p.z = rescaleFactor.z * (x[j].position.z-boxcenter.z);
	}
      else
	{
	particle_list[i].p.x = x[j].position.x;
	particle_list[i].p.y = x[j].position.y;
	particle_list[i].p.z = x[j].position.z;
	}
      particle_list[i].q = x[j].charge * unitFactor;
      DebugM(1,"atom[" << i << "]=" << x[j] << " "
	      << x[j].charge*unitFactor << "\n");
      i++;
      if (i > totalAtoms)
	{
	iout << iERRORF << iPE << " totalAtoms=" << totalAtoms
	     << " but " << i << " atoms are seen!\n" << endi;
	NAMD_die("FMA: atom counts unequal!");
	}
    }

    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
  } 

  if (i != totalAtoms)
  {
    iout << iERRORF << iPE << " totalAtoms=" << totalAtoms
         << " but " << i << " atoms are seen!\n" << endi;
    NAMD_die("FMA: atom counts unequal!");
  }

  DebugM(2,"DPMTA doWork() there are " << totalAtoms << " atoms in this node.\n");

#ifdef DUMP_DPMTA
  FILE *fp;
  char dump_file[32];
  sprintf(dump_file,"DUMP_DPMTA.%d",(int)CkMyPe());
  fp = fopen(dump_file,"w");
  int32 n32 = i;
  fwrite(&n32,sizeof(int32),1,fp);
  fwrite(particle_list,sizeof(PmtaParticle),i,fp);
  fclose(fp);
#endif

  // 3. (run DPMTA) compute the forces
  if ( PMTAforce(i, particle_list, fmaResults, NULL) < 0 )
    {
      NAMD_die("PMTAforce failed!!");
    }
  DebugM(2,"DPMTA forces done.  Now depositing.\n");

  // 4. deposit
  BigReal potential=0;
  for (i=0, ap = ap.begin(); ap != ap.end(); ap++)
  {
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::slow];

    // deposit here
    for(j=0; j<(*ap).p->getNumAtoms(); j++)
    {
      if (usePBC)
	{
	f[j].x += fmaResults[i].f.x * rescaleFactor.x * rescaleFactor.x;
	f[j].y += fmaResults[i].f.y * rescaleFactor.y * rescaleFactor.y;
	f[j].z += fmaResults[i].f.z * rescaleFactor.z * rescaleFactor.z;
	potential += fmaResults[i].v * rescaleFactor.x;
	}
      else
	{
	f[j].x += fmaResults[i].f.x;
	f[j].y += fmaResults[i].f.y;
	f[j].z += fmaResults[i].f.z;
	potential += fmaResults[i].v;
	}
      i++;
    }

    (*ap).forceBox->close(&r);
  }

  potential *= 0.5;
  DebugM(4,"Full-electrostatics energy: " << potential << "\n");
  reduction->item(REDUCTION_ELECT_ENERGY_SLOW) += potential;
  // DPMTA won't work correctly if scaled anisotropically anyway.  -JCP
  reduction->item(REDUCTION_VIRIAL_SLOW_XX) += potential / 3.;
  reduction->item(REDUCTION_VIRIAL_SLOW_YY) += potential / 3.;
  reduction->item(REDUCTION_VIRIAL_SLOW_ZZ) += potential / 3.;
  reduction->submit();

  // 5. clean-up
  if (totalAtoms > 0)
  {
    free(particle_list);
    free(fmaResults);
  }

  DebugM(2,"DPMTA doWork() done\n");
}

#endif

