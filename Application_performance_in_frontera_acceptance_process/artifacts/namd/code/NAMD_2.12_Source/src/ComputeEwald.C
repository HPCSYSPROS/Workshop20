/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Forwards atoms to master node for force evaluation.
*/

#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeEwald.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "ComputeNonbondedUtil.h"
#include "SimParameters.h"
#include "PmeBase.h"
#include <stdio.h>

//#define DEBUGM
#define MIN_DEBUG_LEVEL 1
#include "Debug.h"

class EwaldParticle {
public:
  float x, y, z;
  float cg;
};

// table mapping atom types to array of structure factor combinations.
// For example, for ntypes=3, we have grids corresponding to atom types
// 0-0, 0-1, 0-2, 1-1, 1-2, 2-2, so type 0 gets added to grids 0,1,2, type 1
// to grids 1,3,4, and typ 2 to 2, 4, 5.
static int *generateAtomTypeTable(int ntypes) {
  int *table = new int[ntypes*ntypes];
  int ind = 0;
  for (int i=0; i<ntypes; i++) {
    for (int j=i; j<ntypes; j++) {
      table[ntypes*i+j] = table[ntypes*j+i] = ind++;
    }
  }
  return table;
}

ComputeEwald::ComputeEwald(ComputeID c, ComputeMgr *m)
	: ComputeHomePatches(c)
{
  DebugM(3,"Constructing client\n");
  comm = m;
  SimParameters *sp = Node::Object()->simParameters;
  kxmax = sp->pressureProfileEwaldX;
  kymax = sp->pressureProfileEwaldY;
  kzmax = sp->pressureProfileEwaldZ;
  
  ktot = (1+kxmax) * (2*kymax+1) * (2*kzmax+1);
  kappa = ComputeNonbondedUtil::ewaldcof;
  pressureProfileSlabs = sp->pressureProfileSlabs;
  numAtomTypes = sp->pressureProfileAtomTypes;
  int nelements = 3*pressureProfileSlabs * (numAtomTypes*(numAtomTypes+1))/2;
  pressureProfileData = new float[nelements];
  reduction = ReductionMgr::Object()->willSubmit(
				REDUCTIONS_PPROF_NONBONDED, nelements);

  // figure out who da masta be
  numWorkingPes = (PatchMap::Object())->numNodesWithPatches();
  masterNode = numWorkingPes - 1;

  recvCount = 0;
  localAtoms = NULL;
  localPartitions = NULL;

  expx = new float[kxmax+1];
  expy = new float[kymax+1];
  expz = new float[kzmax+1];

  if (CkMyPe() == masterNode) {
    eiktotal = new float[2 * ktot * numAtomTypes]; 
    memset(eiktotal, 0, 2 * ktot * numAtomTypes*sizeof(float));
  } else {
    eiktotal = NULL;
  }
  // space for exp(iky), k=-kymax, ..., kymax
  eiky = new floatcomplex[2*kymax+1];
  // space for exp(ikz), k=-kzmax, ..., kzmax
  eikz = new floatcomplex[2*kzmax+1];
  Qk = new float[3*ktot];

  gridsForAtomType = generateAtomTypeTable(numAtomTypes);
}

ComputeEwald::~ComputeEwald()
{
  delete reduction;
  delete [] expx;
  delete [] expy;
  delete [] expz;
  delete [] eiktotal;
  delete [] eiky;
  delete [] eikz;
  delete [] pressureProfileData;
  delete [] Qk;
  delete [] gridsForAtomType;
  
  if (localAtoms) free(localAtoms);
  if (localPartitions) free(localPartitions);
}

void ComputeEwald::doWork() {
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


  lattice = patchList[0].p->lattice;
  Vector o = lattice.origin();

  // recompute pressure profile cell parameters based on current lattice
  pressureProfileThickness = lattice.c().z / pressureProfileSlabs;
  pressureProfileMin = lattice.origin().z - 0.5*lattice.c().z;

  const BigReal coulomb_sqrt = sqrt( COULOMB * ComputeNonbondedUtil::scaling
				* ComputeNonbondedUtil::dielectric_1 );

  // get coordinates and store them
  numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }
  localAtoms = (EwaldParticle *)realloc(localAtoms, numLocalAtoms*sizeof(EwaldParticle));
  localPartitions = (int *)realloc(localPartitions, numLocalAtoms*sizeof(int));

  EwaldParticle *data_ptr = localAtoms;
  int *part_ptr = localPartitions;

  for (ap = ap.begin(); ap != ap.end(); ap++) {
    CompAtom *x = (*ap).positionBox->open();
    // CompAtomExt *xExt = (*ap).p->getCompAtomExtInfo();
    Results *r = (*ap).forceBox->open();
    if ( patchList[0].p->flags.doMolly ) {
      (*ap).positionBox->close(&x);
      x = (*ap).avgPositionBox->open();
    }
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i) {
      // wrap back to unit cell, centered on origin
      Vector pos = x[i].position;
      pos += lattice.wrap_delta(pos) - o;
      *part_ptr++ = x[i].partition; 
      data_ptr->x = pos.x;
      data_ptr->y = pos.y;
      data_ptr->z = pos.z;
      data_ptr->cg = coulomb_sqrt * x[i].charge;
      ++data_ptr;
    }
    (*ap).positionBox->close(&x);
    (*ap).forceBox->close(&r);
  }

  // compute structure factor contribution from local atoms
  // 2*ktot since charm++ uses float instead of floatcomplex
  int msgsize = 2 * numAtomTypes * ktot;
  ComputeEwaldMsg *msg = new (msgsize,0) ComputeEwaldMsg;
  memset(msg->eik, 0, msgsize*sizeof(float));
  compute_structurefactor(msg->eik);

  // send our partial sum
  comm->sendComputeEwaldData(msg);
}

void ComputeEwald::recvData(ComputeEwaldMsg *msg) {
  // sum the data...
  int nvecs = 2 * ktot * numAtomTypes;
  for (int i=0; i<nvecs; i++) {
    eiktotal[i] += msg->eik[i];
  }
  delete msg;
  if (++recvCount == numWorkingPes) {
    recvCount = 0;
    int msgsize = 2 * ktot * numAtomTypes;
    msg = new (msgsize,0) ComputeEwaldMsg; 
    memcpy(msg->eik, eiktotal, msgsize*sizeof(float));
    memset(eiktotal, 0, msgsize*sizeof(float));
    comm->sendComputeEwaldResults(msg);
  }
}

void ComputeEwald::recvResults(ComputeEwaldMsg *msg) {
  // receive total sum
  computePprofile(msg->eik);
  delete msg;
  float scalefac = 1.0 / (2 * M_PI * lattice.volume());

  int nelements = 3*pressureProfileSlabs * (numAtomTypes*(numAtomTypes+1))/2;
  for (int i=0; i<nelements; i++) {
    reduction->item(i) += pressureProfileData[i] * scalefac;
  }
  reduction->submit();
}

void ComputeEwald::compute_structurefactor(float *eik) {

  float recipx = lattice.a_r().x;
  float recipy = lattice.b_r().y;
  float recipz = lattice.c_r().z;
  int i, j;
  for (i=0; i<numLocalAtoms; i++) {
    // compute exp(2 pi i r/L) 
    float krx = 2 * M_PI * localAtoms[i].x * recipx;
    float kry = 2 * M_PI * localAtoms[i].y * recipy;
    float krz = 2 * M_PI * localAtoms[i].z * recipz;
    float cg = localAtoms[i].cg;
    // sum structure factors for each atom type separately
    const int offset = 2*ktot*localPartitions[i];
    
    floatcomplex eikx1(cos(krx), sin(krx));
    floatcomplex eiky1(cos(kry), sin(kry));
    floatcomplex eikz1(cos(krz), sin(krz));

    // store exp(2 pi i y j/Ly) for j = -kymax, ..., kymax
    floatcomplex *ptr = eiky + kymax;
    // n=0 case
    ptr[0] = 1;
    floatcomplex y(1,0);
    for (j=1; j<= kymax; j++) {
      y *= eiky1;
      ptr[j] = y;
      ptr[-j] = y.star();
    }

    // store exp(2 pi i y j/Ly) for j = -kzmax, ..., kzmax
    ptr = eikz + kzmax;
    // n=0 case
    ptr[0] = 1;
    floatcomplex z(1,0);
    for (j=1; j<= kzmax; j++) {
      z *= eikz1;
      ptr[j] = z;
      ptr[-j] = z.star();
    }

    // now loop over all k-vectors, computing S(k)
    floatcomplex ex(cg); 
    int ind = offset; // index into eik
    for (int kx=0; kx <= kxmax; kx++) {
      for (int ky=0; ky <= 2*kymax; ky++) {
        floatcomplex exy = eiky[ky];
        exy *= ex;
				const int max = 2*kzmax;
				const float exyr = exy.r;
				const float exyi = exy.i;
		    const float *eikzd = (const float *)eikz;
#pragma vector always
        for (int kz=0; kz <= max; kz++) {
					float ezr = *eikzd++;
					float ezi = *eikzd++;
					float eikr = ezr*exyr - ezi*exyi;
					float eiki = ezr*exyi + ezi*exyr;

          // add this contribution to each grid to which the atom belongs
          eik[ind  ] += eikr;
          eik[ind+1] += eiki;

          // next k vector
          ind += 2;
        }
      }
      ex *= eikx1;
    }
  }
}

// compute exp(-k^2/4 kappa^2) for k=2pi*recip*n, n=0...K inclusive
static void init_exp(float *xp, int K, float recip, float kappa) {
  float piob = M_PI / kappa;
  piob *= piob;
  float fac = -piob*recip*recip;
  for (int i=0; i<= K; i++)
    xp[i] = exp(fac*i*i);
}

void ComputeEwald::computePprofile(const float *eik) const { 
  float recipx = lattice.a_r().x;
  float recipy = lattice.b_r().y;
  float recipz = lattice.c_r().z;

  init_exp(expx, kxmax, recipx, kappa);
  init_exp(expy, kymax, recipy, kappa);
  init_exp(expz, kzmax, recipz, kappa);

  //float energy = 0;
  float piob = M_PI / kappa;
  piob *= piob;

  // compute exp(-pi^2 m^2 / B^2)/m^2 
  int ind = 0;
  for (int kx=0; kx <= kxmax; kx++) {
    float m11 = recipx * kx;
    m11 *= m11;
    float xfac = expx[kx] * (kx ? 2 : 1);
    for (int ky=-kymax; ky <= kymax; ky++) {
      float m22 = recipy * ky;
      m22 *= m22;
      float xyfac = expy[abs(ky)] * xfac;
      for (int kz=-kzmax; kz <= kzmax; kz++) {
        float m33 = recipz * kz;
        m33 *= m33;
        float msq = m11 + m22 + m33;
        float imsq = msq ? 1.0 / msq : 0;
        float fac = expz[abs(kz)] * xyfac * imsq;

        float pfac = 2*(imsq + piob);
        Qk[ind++] = fac*(1-pfac*m11);
        Qk[ind++] = fac*(1-pfac*m22);
        Qk[ind++] = fac*(1-pfac*m33);
      }
    }
  }

  const int nslabs = pressureProfileSlabs;

  int nelements = 3*nslabs * (numAtomTypes*(numAtomTypes+1))/2;
  memset(pressureProfileData, 0, nelements*sizeof(float));
  int i, j;
  for (i=0; i<numLocalAtoms; i++) {

    float krx = 2 * M_PI * localAtoms[i].x * recipx;
    float kry = 2 * M_PI * localAtoms[i].y * recipy;
    float krz = 2 * M_PI * localAtoms[i].z * recipz;
    float cg = localAtoms[i].cg;
    int atype = localPartitions[i];
    const int *grids = gridsForAtomType+atype*numAtomTypes;

    // determine the slab where this particle is located
    int slab = (int)floor((localAtoms[i].z - pressureProfileMin)/pressureProfileThickness);
    if (slab < 0) slab += nslabs;
    else if (slab >= nslabs) slab -= nslabs;
    float *pprofptr = pressureProfileData + 3*slab;

    floatcomplex eikx1(cos(krx), sin(krx));
    floatcomplex eiky1(cos(kry), sin(kry));
    floatcomplex eikz1(cos(krz), sin(krz));

    // store exp(2 pi i y j/Ly) for j = -kymax, ..., kymax
    floatcomplex *ptr = eiky + kymax;
    // n=0 case
    ptr[0] = 1;
    floatcomplex y(1,0);
    for (j=1; j<= kymax; j++) {
      y *= eiky1;
      ptr[j] = y;
      ptr[-j] = y.star();
    }

    // store exp(2 pi i y j/Ly) for j = -kzmax, ..., kzmax
    ptr = eikz + kzmax;
    // n=0 case
    ptr[0] = 1;
    floatcomplex z(1,0);
    for (j=1; j<= kzmax; j++) {
      z *= eikz1;
      ptr[j] = z;
      ptr[-j] = z.star();
    }

    int ind = 0;
    const float *Qkptr = Qk;
    floatcomplex ex(cg);
    const float *eikptr = eik;
    for (int kx=0; kx <= kxmax; kx++) {
      for (int ky=0; ky <= 2*kymax; ky++) {
        floatcomplex exy = eiky[ky];
        exy *= ex;
        const int kzmax2 = 2*kzmax+1;
        const float exyr = exy.r;
        const float exyi = exy.i;
        const float *eikzptr = (const float *)eikz;
#pragma vector always
        for (int kz=0; kz < kzmax2; kz++) {
          float ezr = *eikzptr++;
          float ezi = *eikzptr++;
          float exyzr = ezr * exyr - ezi * exyi;
          float exyzi = ezr * exyi + ezi * exyr;

          // exyz holds exp(ikr) for this atom
          // loop over atom types, adding contribution from each
          // Re[exp(-ikr) * S(k)]
          // add this contribution to each grid the atom belongs to
          for (int igrid=0; igrid<numAtomTypes; igrid++) {
            float eikr = eikptr[igrid*2*ktot  ];
            float eiki = eikptr[igrid*2*ktot+1];
            int grid = grids[igrid];
            int offset = 3*nslabs*grid;

            float E = exyzr * eikr + exyzi * eiki;
            float vx = Qkptr[0] * E;
            float vy = Qkptr[1] * E;
            float vz = Qkptr[2] * E;
            pprofptr[offset  ] += vx;
            pprofptr[offset+1] += vy;
            pprofptr[offset+2] += vz;
          }
          // next k vector
          eikptr += 2;
          Qkptr += 3;
        }
      }
      ex *= eikx1;
    }
  }
}

