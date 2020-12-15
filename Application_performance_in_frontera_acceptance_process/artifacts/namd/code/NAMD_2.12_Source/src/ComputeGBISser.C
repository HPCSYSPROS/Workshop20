/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*******************************************************************************
 *******************************************************************************
  This serial version of GBIS is out of date and hasn't been vetted;
  it is to be used as a way of testing new implicit solvent models
  since all atomic coordinates are gathered into a single Compute.
 *******************************************************************************
 ******************************************************************************/

/*
  print debug
  1 - once per step
  2 - once per atom
  3 - once per pair
*/
//#define PRINT_SERFORCES
#define BENCH_PERIOD 1000

//sum all forces
#define GBIS_DEDR_FORCE 1
#define GBIS_DEDA_FORCE 1
#define GBIS_COUL_FORCE 1

#include "Vector.h"
#include <limits>
#include "InfoStream.h"
#include "Node.h"
#include "PatchMap.h"
#include "PatchMap.inl"
#include "AtomMap.h"
#include "ComputeGBISser.h"
#include "ComputeGBISserMgr.decl.h"
#include "PatchMgr.h"
#include "Molecule.h"
#include "ReductionMgr.h"
#include "ComputeMgr.h"
#include "ComputeMgr.decl.h"
#include "ComputeGBIS.inl"
// #define DEBUGM
#define MIN_DEBUG_LEVEL 3
#include "Debug.h"
#include "SimParameters.h"
#include "WorkDistrib.h"
#include "varsizemsg.h"
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>

  struct vect {
    BigReal x, y, z;
  };
#ifndef COUL_CONST
#define COUL_CONST 332.0636 // ke [kcal*Ang/e^2] 
#endif

#ifndef TA
#define TA 0.333333333333333 // 1/3
#define TB 0.4               // 2/5
#define TC 0.428571428571428 // 3/7
#define TD 0.444444444444444 // 4/9
#define TE 0.454545454545454 // 5/11
#define DA 1.333333333333333 // 4* 1/3
#define DB 2.4               // 6* 2/5
#define DC 3.428571428571428 // 8* 3/7
#define DD 4.444444444444444 // 10*4/9
#define DE 5.454545454545454 // 12*5/11
#endif

inline void Phase2_PairSer(//doesn't do self energies

//input values
  BigReal r,
  BigReal r2,
  BigReal r_i,
  BigReal qiqj,
  BigReal ai,
  BigReal aj,
  BigReal epsilon_p_i,
  BigReal epsilon_s_i,
  BigReal kappa,
  int exclij,
  BigReal scale14,
  int stat,

//return values
  BigReal & coulEij,
  BigReal & ddrCoulEij,
  BigReal & gbEij,
  BigReal & ddrGbEij,
  BigReal & dEdai,
  BigReal & dEdaj
);

inline void CalcHPairSer (
  BigReal r,//distance
  BigReal r2,
  BigReal ri,
  BigReal rc,//cutoff
  BigReal ri0,
  BigReal rjs,
  BigReal rj0,
  BigReal ris,
  int & dij,//domain 1
  int & dji,//domain 2
  BigReal & dhij,
  BigReal & dhji);
inline void CalcDHPairSer (
  BigReal r,//distance
  BigReal r2,
  BigReal ri,
  BigReal rc,//cutoff
  BigReal ri0,
  BigReal rjs,
  BigReal rj0,
  BigReal ris,
  int & dij,//domain 1
  int & dji,//domain 2
  BigReal & dhij,
  BigReal & dhji);


struct ComputeGBISAtom {
    Position position;
    Real charge;
    Real rho;//coulomb radius
    Real rho0;//coulomb radius
    Real rhos;//coulomb radius
    Mass mass;
    //int type;//atom type for S scaling table
    int id;
    int vdwType;
};

class GBISCoordMsg : public CMessage_GBISCoordMsg {
public:
  int sourceNode;
  int numAtoms;
  int doSlow;
  ComputeGBISAtom *coord;
  int sequence;
};

class GBISForceMsg : public CMessage_GBISForceMsg {
public:
  BigReal gbInterEnergy;
  BigReal gbSelfEnergy;
  BigReal coulEnergy;
  ExtForce *force;
  ExtForce *slowForce;
};

class ComputeGBISserMgr : public CBase_ComputeGBISserMgr {
public:
  ComputeGBISserMgr();
  ~ComputeGBISserMgr();

  void setCompute(ComputeGBISser *c) { gbisCompute = c; }
  void recvCoord(GBISCoordMsg *);
  void recvForce(GBISForceMsg *);

private:
  CProxy_ComputeGBISserMgr gbisProxy;
  ComputeGBISser *gbisCompute;

  int numSources;
  int numArrived;
  GBISCoordMsg **coordMsgs;
  int numAtoms;
  ComputeGBISAtom *coord;
  ExtForce *force;
  ExtForce *slowForce;
  GBISForceMsg *oldmsg;

  BigReal gbSelfEnergy;
  BigReal gbInterEnergy;
  BigReal coulEnergy;
  int timestep;
  clock_t t_start;
  clock_t t_stop;
  clock_t t1, t2;
  double totalnamdtime;
  double totalgbistime;
  double inittime;
  double psitime;
  double alphatime;
  double dEdrtime;
  double dEdasumtime;
  double dEdaprefixtime;
  double dEdalooptime;
  void calcGBISReg(int stat);
  int blockSize;
  int numBlocks; 
  int loop1Iter, loop2Iter, loop3Iter;
  bool all2all;
  BigReal MinErr;
  BigReal MaxErr;
  int numBins;
  BigReal base;
};

ComputeGBISserMgr::ComputeGBISserMgr() :
  gbisProxy(thisgroup), gbisCompute(0), numSources(0), numArrived(0),
  inittime(0), 
  psitime(0), 
  alphatime(0), 
  dEdrtime(0), 
  dEdasumtime(0), 
  dEdaprefixtime(0), 
  dEdalooptime(0), slowForce(0),
  coordMsgs(0), coord(0), force(0), oldmsg(0), numAtoms(0), timestep(0) {
  CkpvAccess(BOCclass_group).computeGBISserMgr = thisgroup;
  t_start = clock();
  t_stop = clock();
  all2all = false;
}

ComputeGBISserMgr::~ComputeGBISserMgr() {
  for ( int i=0; i<numSources; ++i ) { delete coordMsgs[i]; }
  delete [] coordMsgs;
  delete [] coord;
  delete [] force;
  delete [] slowForce;
  delete oldmsg;
}

ComputeGBISser::ComputeGBISser(ComputeID c) :
  ComputeHomePatches(c)
{
  CProxy_ComputeGBISserMgr::ckLocalBranch(
    CkpvAccess(BOCclass_group).computeGBISserMgr)->setCompute(this);

  reduction = ReductionMgr::Object()->willSubmit(REDUCTIONS_BASIC);

}

ComputeGBISser::~ComputeGBISser()
{
}

/******************************************************************************
 * scale from 1 to 0 between transition and cutoff
 * s = 1 for [0,trans]
 * s = s(r) for [trans, cutoff]
 * s = 0 for [cutoff, infinity]
 * meets the following 4 conditions:
 * (1) s(trans)   = 1
 * (3) s'(trans)  = 0
 * (2) s(cutoff)  = 0
 * (4) s'(cutoff) = 0
 ******************************************************************************/
inline void CalcScaleSer (
  BigReal r,//rij
  BigReal t,//transition distance
  BigReal c,//cutoff distance
  BigReal & s,//output s(r)
  BigReal & d ) {//output s'(r)

  if (r <= t) { //[0,trans]
    s = 1;
    d = 0;
  } else if (r < c) { //[trans,cutoff]
    //precompute
    BigReal ct = (c-t);
    BigReal ct_i = 1/ct;
    BigReal ct_i2 = ct_i*ct_i;
    //BigReal ct_i4 = ct_i2*ct_i2;

    BigReal rt = r - t;
    BigReal rt2 = rt*rt;
    BigReal omrtct2 = 1-rt2*ct_i2;
    s = omrtct2*omrtct2;
    
    BigReal rc = r - c;
    d=s*4*rt/(rc*(rt+ct));

  } else { //[cutoff,infinity]
    s = 0;
    d = 0;
  }
}


/**********************************************************
 *
 *  Regular 3 sets of nested loops
 *
 *  supports 1-2-4 scheme
 *
 **********************************************************/
void ComputeGBISserMgr::calcGBISReg(int stat) {
  t1 = clock();
  //printf("GBIS: ComputeGBIS::serial\n");

  /*********************************************************
   *  Initializations
   ********************************************************/
  Molecule *molecule = Node::Object()->molecule;
  SimParameters *simParams = Node::Object()->simParameters;
  BigReal epsilon_s = simParams->solvent_dielectric;
  BigReal epsilon_p = simParams->dielectric;
  BigReal epsilon_s_i = 1/epsilon_s;
  BigReal epsilon_p_i = 1/epsilon_p;
  BigReal rho_0 = simParams->coulomb_radius_offset;
  BigReal cutoff = simParams->cutoff;// + 1.8*0.96;//+max screen radius
  BigReal trans = simParams->switchingDist;
  BigReal fsMax = 1.728;
  BigReal a_cut = simParams->alpha_cutoff-fsMax;//+max screen radius; partial sphere
  //a_cut = 10;
  BigReal a_cut2 = a_cut*a_cut;
  BigReal a_cut_ps = a_cut + fsMax;//+max screen radius; partial sphere
  BigReal cutoff2 = cutoff*cutoff;
  BigReal a_cut_ps2 = a_cut_ps*a_cut_ps;
  BigReal delta = simParams->gbis_delta;
  BigReal beta = simParams->gbis_beta;
  BigReal gamma = simParams->gbis_gamma;
  BigReal alphaMax = simParams->alpha_max;
  int exclude = simParams->exclude;
  BigReal scale14 = (exclude==SCALED14) ? simParams->scale14 : 1.0;


#if PRINT > 0
  //printf("cut=%e, rgbmax=%e, fsmax=%e, gbalpha=%e, gbbeta=%e, gbgamma=%e\n", cutoff, alphaMax, fsMax, delta, beta, gamma);
#endif

  //calc kappa
  BigReal kappa = simParams->kappa;

  //debugging variables
  gbInterEnergy = 0;
  gbSelfEnergy = 0;
  coulEnergy = 0;
  BigReal fx, fy, fz;
  loop1Iter = 0;
  loop2Iter = 0;
  loop3Iter = 0;
  int loop1flops = 0;
  int loop2flops = 0;
  int loop3flops = 0;
  BigReal loop1time = 0;
  BigReal loop2time = 0;
  BigReal loop3time = 0;

  Position ri, rj;
  BigReal rhoi, rhoi0, rhois;
  BigReal rhoj, rhoj0, rhojs;
  BigReal dx, dy, dz;
  BigReal r, r2, r_i, qiqj;
  BigReal rnx, rny, rnz;
  BigReal bornRad;
  BigReal aiaj, aiaj4, expr2aiaj4, fij, finv, expkappa, Dij;
  BigReal rc, shift, ddrshift;
  BigReal tmp_dEda, dEdai, dEdaj, gbEij;
  BigReal ddrfij, ddrfinv, ddrDij, ddrGbEij;
  BigReal coulEij, ddrCoulEij;
  BigReal forceCoul, forcedEdr, forceAlpha;
  BigReal hij, hji, ddrHij, ddrHji;
  BigReal tanhi, daidr, dajdr;
  BigReal nbetapsi, gammapsi2;
  vect *coulForce = new vect[numAtoms];
  vect *gbForceR = new vect[numAtoms];//dEdr
  vect *gbForceA = new vect[numAtoms];//dEda*dadr
  BigReal *a = new BigReal[numAtoms];
  BigReal *psi = new BigReal[numAtoms];
  BigReal *dEda = new BigReal[numAtoms];
  BigReal *ddrHijPrefix = new BigReal[numAtoms];

  //init arrays
  for (int i = 0; i < numAtoms; i++) {
    force[i].force.x =0.0;
    force[i].force.y =0.0;
    force[i].force.z =0.0;
    slowForce[i].force.x =0.0;
    slowForce[i].force.y =0.0;
    slowForce[i].force.z =0.0;
    dEda[i] = 0;

    psi[i] = 0;
    ddrHijPrefix[i] = 0;
    coulForce[i].x = 0;
    coulForce[i].y = 0;
    coulForce[i].z = 0;
    gbForceR[i].x = 0;//dEdr
    gbForceR[i].y = 0;//dEdr
    gbForceR[i].z = 0;//dEdr
    gbForceA[i].x = 0;//dEda*dadr
    gbForceA[i].y = 0;//dEda*dadr
    gbForceA[i].z = 0;//dEda*dadr
  }

  t2 = clock();
  inittime = (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;

  /**********************************************************
   *
   *  Loop 1 : accumulate Hij into psi array
   *  cutoff = a_cut = 25A
   *  every 2 timesteps
   *
   **********************************************************/
  double missTime = 0;
  double hitTime = 0;
  t1 = clock();
  int dij, dji;//which hij domain 0123456 - for debugging only
  BigReal dhij, dhji;
  //int TESTITER = 1;
  //for (int k = 0; k < TESTITER; k++)
  for (int i = 0; i < numAtoms; i++) {
    ri = coord[i].position;
    for (int j = i+1; j < numAtoms; j++) {
      rj = coord[j].position;
      dx = (ri.x - rj.x);
      dy = (ri.y - rj.y);
      dz = (ri.z - rj.z);
      r2 = dx*dx+dy*dy+dz*dz;
      if (r2 > a_cut_ps2) continue;
      r_i = 1.0/sqrt(r2);;
      r = 1.0/r_i;;
      loop1Iter++;
          BigReal rhoi0 = coord[i].rho0;
          BigReal rhojs = coord[j].rhos;
          BigReal rhoj0 = coord[j].rho0;
          BigReal rhois = coord[i].rhos;
      CalcHPairSer(r,r2,r_i,a_cut, coord[i].rho0, coord[j].rhos,
                coord[j].rho0, coord[i].rhos,dij,dji,hij,hji);

#ifdef PRINT_COMP
      CkPrintf("PSI(%04i)[%04i,%04i] = %i%i % .4e % .4e\n",timestep,coord[i].id,coord[j].id,dij,dji,hij,hji);
      //CkPrintf("S_Hij_%i\n",dij);
      //CkPrintf("S_Hji_%i\n",dji);
#endif
      psi[i] += hij;
      psi[j] += hji;

    }
  }
  t2 = clock();
  psitime += (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;
  loop1time += (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;
  t1 = clock();



  /**********************************************************
   *
   *  Atom 1: Calculate alpha based on phi for each atom
   *
   **********************************************************/
  BigReal totPsi = 0;
  for (int i = 0; i < numAtoms; i++) {
    //CkPrintf("Srho[%i] = %.3e %.3e\n",coord[i].id,coord[i].rho0,coord[i].rhos);

    rhoi0 = coord[i].rho0;
    rhoi = coord[i].rho;
    //CkPrintf("GBIS_SER: psi[%03i] = %.4e\n",i,psi[i]);
    totPsi += psi[i];
    psi[i] *= rhoi0;
    bornRad=1/(1/rhoi0-1/rhoi*tanh(psi[i]*(delta+psi[i]*(-beta+gamma*psi[i]))));
    bornRad = (1.0/bornRad < 1.0/alphaMax) ? alphaMax : bornRad;
    a[i] = bornRad;
#ifdef PRINT_COMP
    CkPrintf("BORNRAD(%04i)[%04i] = % .4e\n",timestep,coord[i].id,bornRad);
#endif
  }
  t2 = clock();
  alphatime += (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;
  t1 = clock();



  /**********************************************************
   *
   *  Loop 2 : dEda
   *
   **********************************************************/
  for (int i = 0; i < numAtoms; i++) {
    ri = coord[i].position;
    for (int j = i+1; j < numAtoms; j++) {
      rj = coord[j].position;
      dx = (ri.x - rj.x);//rptI
      dy = (ri.y - rj.y);//rptI
      dz = (ri.z - rj.z);//rptI
      r2 = dx*dx+dy*dy+dz*dz;//rptI
      if (r2 > cutoff2) continue;

      //calculate distance
      loop2Iter++;
      qiqj = (-COUL_CONST*coord[i].charge);
      qiqj *= coord[j].charge;
      r_i = 1.0/sqrt(r2);
      r = 1.0/r_i;//rptI
      rnx = dx*r_i;//normalized vector
      rny = dy*r_i;
      rnz = dz*r_i;

      //calculate scaling for energy smoothing
      BigReal scale, ddrScale;
      //CalcScaleSer(r, trans, cutoff, scale, ddrScale);
      if ( false ) {
        BigReal ratio = r2 / cutoff2;
        scale = ratio - 1;
        scale *= scale;
        ddrScale = 4.0*r*(r2-cutoff2)/cutoff2/cutoff2;
      } else {
        scale = 1;
        ddrScale = 0;
      }


      //BigReal ai = a[i];
      //BigReal aj = a[j];
      int exclij = 0;//molecule->checkexcl(i,j);
      //CkPrintf("GBIS_SER_excl[%i,%i] %i\n",i,j,exclij);
      //CkPrintf("GBIS_DOFULL = %i\n",stat);
      
      Phase2_PairSer(r,r2,r_i,qiqj,a[i],a[j],epsilon_p_i,epsilon_s_i,kappa,exclij,
          scale14,stat,coulEij,ddrCoulEij,gbEij,ddrGbEij,dEdai,dEdaj);
    int id1 = coord[i].id;
    int id2 = coord[j].id;
    if (id1 > id2 ) {
        int tmp = id2;
        id2 = id1;
        id1 = tmp;
      }//print ids as lower index, higher index

      //accumulate energies
      gbInterEnergy   += gbEij  *scale;
      forcedEdr = -(ddrGbEij)*scale-(gbEij)*ddrScale;

      //add GB force
      fx = rnx*forcedEdr;
      fy = rny*forcedEdr;
      fz = rnz*forcedEdr;
#ifdef PRINT_COMP
      CkPrintf("DEDR(%04i)[%04i,%04i] = % .4e\n",timestep,i,j,forcedEdr);
      CkPrintf("DASM(%04i)[%04i,%04i] = % .4e % .4e\n",timestep,i,j,dEdai*scale,dEdaj*scale);
      CkPrintf("P2RM(%04i)[%04i,%04i] = % .4e % .4e % .4e % .4e % .4e\n",timestep,i,j, r, a[i],a[j],epsilon_p_i,epsilon_s_i,kappa);

#endif
      gbForceR[i].x += fx;
      gbForceR[i].y += fy;
      gbForceR[i].z += fz;
      gbForceR[j].x -= fx;
      gbForceR[j].y -= fy;
      gbForceR[j].z -= fz;

      //add dEda
      if (stat) {
        dEda[i] += dEdai*scale;
        dEda[j] += dEdaj*scale;
      }

    }//end j loop
  }//end i loop
  t2 = clock();
  dEdasumtime += (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;
  loop2time += (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;
  t1 = clock();



/*******************************************************************************
*
*   Atom 2 : Calculate dHij Prefix
*
*******************************************************************************/
  for (int i = 0; i < numAtoms; i++) {
    //add diagonal dEda term
    fij = a[i];//inf
    expkappa = exp(-kappa*fij);//0
    Dij = epsilon_p_i - expkappa*epsilon_s_i;//dielectric term
    //add diagonal GB energy term
    gbEij = -COUL_CONST*coord[i].charge*coord[i].charge*Dij/fij;
    gbSelfEnergy += 0.5*gbEij;//self energy
//CkPrintf("self_energy[%03i] = %e\n",i,0.5*gbEij);
    //CkPrintf("gbSelfEnergy[%03i] = % e\n", coord[i].id,0.5*gbEij);

  //calculate dHij prefix
    if (stat) {
      dEdai = -0.5*COUL_CONST*coord[i].charge*coord[i].charge
                  *(kappa*epsilon_s_i*expkappa-Dij/fij)/a[i];
      //CkPrintf("SER_dEdai[%03i]%i = % e\n",i,timestep,dEdai);
      dEda[i] += dEdai;
      BigReal dedasum = dEda[i];
      //CkPrintf("SER_dEdaSum[%03i]%i = % e\n",i,timestep,dEda[i]);

      rhoi0 = coord[i].rho0;
      rhoi = rhoi0+rho_0;
      BigReal psii = psi[i];
      nbetapsi = -beta*psii;
      gammapsi2 = gamma*psii*psii;
      tanhi = tanh(psii*(delta+nbetapsi+gammapsi2));
      daidr = a[i]*a[i]*rhoi0/rhoi*(1-tanhi*tanhi)
           * (delta+nbetapsi+nbetapsi+gammapsi2+gammapsi2+gammapsi2);
      ddrHijPrefix[i] = daidr*dEda[i];
#ifdef PRINT_COMP
    CkPrintf("DHDR(%04i)[%04i] = % .4e\n",timestep,coord[i].id, ddrHijPrefix[i]);
#endif
    }
  }
  t2 = clock();
  dEdaprefixtime += (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;
  //gbEnergy = gbInterEnergy + gbSelfEnergy;
  t1 = clock();

  /**********************************************************
   *
   *  Loop 3 : Calculate dEda*dadr Forces
   *
   **********************************************************/
  if (stat) {
  BigReal dhij, dhji;
  int dij, dji;
  //for (int k = 0; k < TESTITER; k++)
  for (int i = 0; i < numAtoms; i++) {
    ri = coord[i].position;
    //CkPrintf("SER_dHdrPrefix[%i]%i = % .5e\n",coord[i].id,timestep,ddrHijPrefix[i]);

    for (int j = i+1; j < numAtoms; j++) {
      rj = coord[j].position;
      dx = (ri.x - rj.x);//rptI
      dy = (ri.y - rj.y);//rptI
      dz = (ri.z - rj.z);//rptI
      r2 = dx*dx+dy*dy+dz*dz;//rptI
      //loop3flops += 9;
      if (r2 > a_cut_ps2) continue;
      r_i = 1.0/sqrt(r2);//rptI
      loop3Iter++;
      //calculate dHij for pair
      r = 1.0/r_i;
      CalcDHPairSer(r,r2,r_i,a_cut,
          coord[i].rho0,coord[j].rhos,
          coord[j].rho0,coord[i].rhos,
          dij,dji,dhij,dhji);

      //calculate and add dEijdai,j*dai,jdrij force
      forceAlpha = -r_i*(ddrHijPrefix[i]*dhij+ddrHijPrefix[j]*dhji);
#ifdef PRINT_COMP
      CkPrintf("DEDA(%04i)[%04i,%04i] = %i%i % .4e % .4e % .4e\n",timestep,i,j,dij,dji,dhij,dhji, forceAlpha/r_i);
#endif
      fx = dx * forceAlpha;
      fy = dy * forceAlpha;
      fz = dz * forceAlpha;

      gbForceA[i].x += fx;
      gbForceA[i].y += fy;
      gbForceA[i].z += fz;
      gbForceA[j].x -= fx;
      gbForceA[j].y -= fy;
      gbForceA[j].z -= fz;
    }//end inner summation loop
  }//end outer summation loop
  }//end if stat
  t2 = clock();
  dEdalooptime += (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;
  loop3time += (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;
  t1 = clock();
  
  //add forces
  for (int i = 0; i < numAtoms; i++) {
    force[i].force.x =0.0;
    force[i].force.y =0.0;
    force[i].force.z =0.0;
    slowForce[i].force.x =0.0;
    slowForce[i].force.y =0.0;
    slowForce[i].force.z =0.0;

#if GBIS_COUL_FORCE
    force[i].force.x += coulForce[i].x;
    force[i].force.y += coulForce[i].y;
    force[i].force.z += coulForce[i].z;
#endif

#if GBIS_DEDR_FORCE
    force[i].force.x += gbForceR[i].x;
    force[i].force.y += gbForceR[i].y;
    force[i].force.z += gbForceR[i].z;
#endif

#if GBIS_DEDA_FORCE
    if (stat) {
      slowForce[i].force.x += gbForceA[i].x;
      slowForce[i].force.y += gbForceA[i].y;
      slowForce[i].force.z += gbForceA[i].z;
      //CkPrintf("SERIAL SLOW %e %e %e\n",gbForceA[i].x,gbForceA[i].y,gbForceA[i].z);
    }
    //force[i].force.x += gbForceA[i].x;
    //force[i].force.y += gbForceA[i].y;
    //force[i].force.z += gbForceA[i].z;
#endif


#ifdef PRINT_SERFORCES
  BigReal fC, fR, fA;
      fx = 0;
      fy = 0;
      fz = 0;
#if GBIS_COUL_FORCE
      fx += coulForce[i].x;
      fy += coulForce[i].y;
      fz += coulForce[i].z;
      //fC = sqrt(fx*fx+fy*fy+fz*fz);
      //fprintf(stderr, "%i %e %e %e ",i,fx, fy, fz);
#endif
#if GBIS_DEDR_FORCE
      fx += gbForceR[i].x;
      fy += gbForceR[i].y;
      fz += gbForceR[i].z;
      //fR = sqrt(fx*fx+fy*fy+fz*fz);
      //fprintf(stderr, "%e %e %e",fx, fy, fz);
#endif
#if GBIS_DEDA_FORCE
      fx += gbForceA[i].x;
      fy += gbForceA[i].y;
      fz += gbForceA[i].z;
      //fA = sqrt(fx*fx+fy*fy+fz*fz);
      //fprintf(stderr, "%5i % .5e % .5e % .5e\n",i, fx, fy, fz);
#endif
      //fprintf(stderr, "%5i % .5e % .5e % .5e\n",i, fx, fy, fz);
#endif //if print>1
  }
  t2 = clock();
  inittime += (double)((double)t2-(double)t1)/CLOCKS_PER_SEC;

  timestep++;

  delete[] coulForce;
  delete[] gbForceR;
  delete[] gbForceA;
  delete[] a;
  delete[] psi;
  delete[] dEda;
  delete[] ddrHijPrefix;
  //printf("GBIS calcReg COMPLETE!\n");
}

void ComputeGBISserMgr::recvCoord(GBISCoordMsg *msg) {
  //printf("GBIS recvCoord()\n");
  if ( ! numSources ) {//init receiving coord
    numSources = (PatchMap::Object())->numNodesWithPatches();
    coordMsgs = new GBISCoordMsg*[numSources];
    for ( int i=0; i<numSources; ++i ) { coordMsgs[i] = 0; }
    numArrived = 0;
    numAtoms = Node::Object()->molecule->numAtoms;
    coord = new ComputeGBISAtom[numAtoms];
    force = new ExtForce[numAtoms];
    slowForce = new ExtForce[numAtoms];
    timestep = msg[0].sequence;
  }
  

  //receive coord
  int i;
  for ( i=0; i < msg->numAtoms; ++i ) {
    coord[msg->coord[i].id] = msg->coord[i];
  }

  coordMsgs[numArrived] = msg;
  ++numArrived;

  if ( numArrived < numSources ) return;
  numArrived = 0;


  /**********************************************************
   *
   *  All sources arrived; calculate energy, forces
   *
   **********************************************************/
  t_start = clock();
  //how long did namd take to cycle back here
  double namdtime = (double)((double)t_start-(double)t_stop)/CLOCKS_PER_SEC;
  totalnamdtime += namdtime;
  //choice choose
  calcGBISReg(msg->doSlow);
  
  t_stop = clock();
  double gbistime = (double)((double)t_stop-(double)t_start)/CLOCKS_PER_SEC;
  //printf("GBIS: elapsednamd(%i)=%f\n", timestep, namdtime);
  //printf("GBIS: elapsedgbis(%i)=%f\n", timestep, gbistime);
  totalgbistime += gbistime;
  //printf("GBIS: total(%i)=%f\n", timestep, totaltime);
  if (timestep % BENCH_PERIOD == 0) {
    printf("\n");
    printf("GBIS:       t_GB=%f sec for %i steps\n",totalgbistime,BENCH_PERIOD);
    printf("GBIS:       t_MD=%f sec for %i steps\n",totalnamdtime,BENCH_PERIOD);
    printf("GBIS:       init=%f sec for %i steps\n", inittime, BENCH_PERIOD);
    printf("GBIS:        psi=%f sec for %i steps\n", psitime, BENCH_PERIOD);
    printf("GBIS:      alpha=%f sec for %i steps\n", alphatime, BENCH_PERIOD);
    printf("GBIS:    dEdasum=%f sec for %i steps\n", dEdasumtime, BENCH_PERIOD);
    printf("GBIS: dEdaprefix=%f sec for %i steps\n",dEdaprefixtime,BENCH_PERIOD);
    printf("GBIS:   dEdaloop=%f sec for %i steps\n", dEdalooptime,BENCH_PERIOD);
    printf("GBIS:      loop1=%i iters\n", loop1Iter);
    printf("GBIS:      loop2=%i iters\n", loop2Iter);
    printf("GBIS:      loop3=%i iters\n", loop3Iter);
    printf("\n");
    totalgbistime = 0;
    totalnamdtime = 0;
    inittime = 0;
    psitime = 0;
    alphatime = 0;
    dEdrtime = 0;
    dEdasumtime = 0;
    dEdaprefixtime = 0;
    dEdalooptime = 0;
  }

  // distribute forces
  //printf("GBIS distributing forces\n");
  for ( int j=0; j < numSources; ++j ) {
    GBISCoordMsg *cmsg = coordMsgs[j];
    coordMsgs[j] = 0;
    GBISForceMsg *fmsg = new (cmsg->numAtoms, cmsg->numAtoms, 0) GBISForceMsg;
    //does this init slowForce, force?
    for ( int i=0; i < cmsg->numAtoms; ++i ) {
      fmsg->force[i] = force[cmsg->coord[i].id];
      fmsg->slowForce[i] = slowForce[cmsg->coord[i].id];
      /*
      if (i == 1000 ) {
        printf("%i force: <", i);
        printf("%f, ", fmsg->force[i].force.x);
        printf("%f, ", fmsg->force[i].force.y);
        printf("%f>\n", fmsg->force[i].force.z);
        printf("%i slwfr: <", i);
        printf("%f, ", fmsg->slowForce[i].force.x);
        printf("%f, ", fmsg->slowForce[i].force.y);
        printf("%f>\n", fmsg->slowForce[i].force.z);
      }
      */
    }

    //CkPrintf("GBISENERGY[%i] c = % e, b = % e\n",timestep,coulEnergy, gbEnergy);
   if ( ! j ) {
#if GBIS_DEDR_FORCE
      fmsg->gbSelfEnergy = gbSelfEnergy;
      fmsg->gbInterEnergy = gbInterEnergy;
#else
      fmsg->gbSelfEnergy = 0;
      fmsg->gbInterEnergy = 0;
#endif
#if GBIS_COUL_FORCE
      fmsg->coulEnergy = coulEnergy;
#else
      fmsg->coulEnergy = 0;
#endif
    } else {
      fmsg->gbSelfEnergy = 0;
      fmsg->gbInterEnergy = 0;
      fmsg->coulEnergy = 0;
    }
    gbisProxy[cmsg->sourceNode].recvForce(fmsg);
    delete cmsg;
  }
  //printf("GBIS distributing forces COMPLETE!\n");
}

void ComputeGBISserMgr::recvForce(GBISForceMsg *msg) {
  //printf("GBIS recvForce()\n");
  gbisCompute->saveResults(msg);
  delete oldmsg;
  oldmsg = msg;
}

void ComputeGBISser::saveResults(GBISForceMsg *msg) {
  //printf("GBIS saveResults()\n");
  ResizeArrayIter<PatchElem> ap(patchList);

  ExtForce *results_ptr = msg->force;
  ExtForce *results_ptr_slow = msg->slowForce;

  // add in forces
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    //CkPrintf("GBIS%i: ComputeGBIS(%i)::saveResults() openedForceBox",CkMyPe(),cid);
    Results *r = (*ap).forceBox->open();
    Force *f = r->f[Results::nbond];
    Force *sf = r->f[Results::slow];
    int numAtoms = (*ap).p->getNumAtoms();

    for(int i=0; i<numAtoms; ++i) {
      f[i] += results_ptr->force;
      sf[i] += results_ptr_slow->force;
    //CkPrintf("GBIS%i: slow[%i] = % e\n",CkMyPe(),i,sf[i].x);
      ++results_ptr;
      ++results_ptr_slow;
    }
    //CkPrintf("GBIS%i: ComputeGBIS(%i)::saveResults() closedForceBox",CkMyPe(),cid);
    (*ap).forceBox->close(&r);
  }

    //reduction->item(REDUCTION_ELECT_ENERGY) += msg->coulEnergy;
    reduction->item(REDUCTION_ELECT_ENERGY) += msg->gbInterEnergy;
    reduction->item(REDUCTION_ELECT_ENERGY) += msg->gbSelfEnergy;
    //CkPrintf("energies= % e, % e, % e\n",msg->coulEnergy, msg->gbInterEnergy, msg->gbSelfEnergy);
    reduction->submit();
}

/******************************************************************************
 * Send data to node 0
 ******************************************************************************/
void ComputeGBISser::doWork()
{
  //printf("GBIS doWork()\n");
  ResizeArrayIter<PatchElem> ap(patchList);

#if 1
 // Skip computations if nothing to do.
 if ( ! patchList[0].p->flags.doNonbonded )
 {
   for (ap = ap.begin(); ap != ap.end(); ap++) {
     CompAtom *x = (*ap).positionBox->open();
     //CkPrintf("GBIS%i: ComputeGBIS(%i)::doWork() openedForceBox",CkMyPe(),cid);
     Results *r = (*ap).forceBox->open();
     (*ap).positionBox->close(&x);
     //CkPrintf("GBIS%i: ComputeGBIS(%i)::doWork() closedForceBox",CkMyPe(),cid);
     (*ap).forceBox->close(&r);
   }
   reduction->submit();
   return;
 }
#endif


  // allocate message
  int numLocalAtoms = 0;
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    numLocalAtoms += (*ap).p->getNumAtoms();
  }

  GBISCoordMsg *msg = new (numLocalAtoms, 0) GBISCoordMsg;
  msg->sourceNode = CkMyPe();
  msg->numAtoms = numLocalAtoms;
  ComputeGBISAtom *data_ptr = msg->coord;
  SimParameters *simParams = Node::Object()->simParameters;
  msg->doSlow = patchList[0].p->flags.doFullElectrostatics;
  //CkPrintf("SERIAL SLOW %i\n",msg->doSlow);
  msg->sequence = sequence();

  // get positions
  for (ap = ap.begin(); ap != ap.end(); ap++) {
    FullAtomList &atoms = (*ap).p->getAtomList();
    CompAtom *x = (*ap).positionBox->open();
    CompAtomExt *xExt = (*ap).p->getCompAtomExtInfo();
    int numAtoms = (*ap).p->getNumAtoms();
    for(int i=0; i<numAtoms; ++i)
    {
      data_ptr->position = x[i].position;
      data_ptr->charge = x[i].charge;
      data_ptr->mass = atoms[i].mass;
      data_ptr->id = xExt[i].id;
      data_ptr->rho = MassToRadius(data_ptr->mass);
      SimParameters *simParams = Node::Object()->simParameters;
      data_ptr->rho0 = data_ptr->rho - simParams->coulomb_radius_offset;
      data_ptr->rhos = data_ptr->rho0 * MassToScreen(data_ptr->mass);
      data_ptr->vdwType = x[i].vdwType;
      ++data_ptr;
    }

#if 0
    if ( patchList[0].p->flags.doMolly ) { (*ap).avgPositionBox->close(&x); }
    else { (*ap).positionBox->close(&x); }
#endif
    (*ap).positionBox->close(&x);
  }

  CProxy_ComputeGBISserMgr gbisProxy(CkpvAccess(BOCclass_group).computeGBISserMgr);
  gbisProxy[0].recvCoord(msg);

}
#include "ComputeGBISserMgr.def.h"

/*
  Piecewise screening functions Hij dHij/drij
  r   distance
  r2  square distance
  ri  inverse distance
  rc  cutoff
  r0  descreened atom radius
  rs  descreening atom radius
  h   return value
  dh  return value
*/
inline void h0Ser ( BigReal r, BigReal r2, BigReal ri,// 5.3%
Real rc, BigReal r0, BigReal rs, BigReal & h ) {
  h = 0;
}
inline void dh0Ser ( BigReal r, BigReal r2, BigReal ri,// 5.3%
Real rc, BigReal r0, BigReal rs, BigReal & dh ) {
  dh = 0;
}

inline void h1Ser ( BigReal r, BigReal r2, BigReal ri, //18.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {

  BigReal rci = 1.0/rc;
  BigReal rmrs = r-rs;// 4 times
  BigReal rmrsi = 1.0/rmrs;
  BigReal rmrs2 = rmrs*rmrs;
  BigReal rs2 = rs*rs;
  BigReal logr = log(rmrs*rci);
  BigReal rci2 = rci*rci;
  h = 0.125*ri*(1 + 2*r*rmrsi + rci2*(r2 - 4*rc*r - rs2) + 2*logr);
}
inline void dh1Ser ( BigReal r, BigReal r2, BigReal ri, //18.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {

  BigReal rci = 1.0/rc;
  BigReal rmrs = r-rs;// 4 times
  BigReal rmrsi = 1.0/rmrs;
  BigReal rmrs2 = rmrs*rmrs;
  BigReal rs2 = rs*rs;
  BigReal logr = log(rmrs*rci);
  BigReal rci2 = rci*rci;
  dh = ri*ri*(-0.25*logr - (rc*rc - rmrs2)*(rs2 + r2)*0.125*rci2*rmrsi*rmrsi);
}

inline void h2Ser ( BigReal r, BigReal r2, BigReal ri,// 74.5%
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {

    BigReal k = rs*ri; k*=k;//k=(rs/r)^2
    h = rs*ri*ri*k*(TA+k*(TB+k*(TC+k*(TD+k*TE))));
}
inline void dh2Ser ( BigReal r, BigReal r2, BigReal ri,// 74.5%
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {

    BigReal k = rs*ri; k*=k;//k=(rs/r)^2
    dh = -rs*ri*ri*ri*k*(DA+k*(DB+k*(DC+k*(DD+k*DE))));
}

inline void h3Ser ( BigReal r, BigReal r2, BigReal ri,// 1.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {
    BigReal r2mrs2i = 1.0/(r2-rs*rs);
    h = 0.5 * ( rs*r2mrs2i + 0.5 * log((r-rs)/(r+rs))*ri );
}
inline void dh3Ser ( BigReal r, BigReal r2, BigReal ri,// 1.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {
    BigReal rs2 = rs*rs;
    BigReal r2mrs2i = 1.0/(r2-rs2);
    dh = -0.25*ri*(2*(r2+rs2)*rs*r2mrs2i*r2mrs2i + ri*log((r-rs)/(r+rs)));
}

inline void h4Ser ( BigReal r, BigReal r2, BigReal ri,// 0.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {
    BigReal ri2 = ri*ri;
    BigReal r02 = r0*r0;
    BigReal rs2 = rs*rs;
    BigReal r0i = 1.0/r0;
    BigReal rspri = 1.0/(r+rs);
    BigReal logr = log(r0*rspri);
    BigReal r02mrs2 = r02-rs2;
    BigReal rilogr = ri*logr;
    h = 0.25*( r0i*(2-0.5*(r0i*ri*(r2 + r02 - rs2))) - rspri + rilogr );
}
inline void dh4Ser ( BigReal r, BigReal r2, BigReal ri,// 0.4%
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {
    BigReal ri2 = ri*ri;
    BigReal r02 = r0*r0;
    BigReal rs2 = rs*rs;
    BigReal r0i = 1.0/r0;
    BigReal rspri = 1.0/(r+rs);
    BigReal logr = log(r0*rspri);
    BigReal r02mrs2 = r02-rs2;
    BigReal rilogr = ri*logr;
    dh = 0.25*( (-0.5+(r2*r02mrs2 - 2*r*rs*rs2+rs2*r02mrs2)
        * 0.5*ri2*rspri*rspri)*r0i*r0i - ri*rilogr );
}

inline void h5Ser ( BigReal r, BigReal r2, BigReal ri,// 0%, r<0.7Ang
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {
    BigReal rs2 = rs*rs;
    BigReal r2mrs2i = 1/(r2-rs2);
    BigReal rsr2mrs2i = rs*r2mrs2i;
    BigReal rprs = r+rs;
    BigReal rmrs = r-rs;
    BigReal logr = 0.5*ri*log(-rmrs/rprs);
    h = 0.5*( rsr2mrs2i + 2/r0 + logr );
}
inline void dh5Ser ( BigReal r, BigReal r2, BigReal ri,// 0%, r<0.7Ang
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {
    BigReal rs2 = rs*rs;
    BigReal r2mrs2i = 1/(r2-rs2);
    BigReal rsr2mrs2i = rs*r2mrs2i;
    BigReal rprs = r+rs;
    BigReal rmrs = r-rs;
    BigReal logr = 0.5*ri*log(-rmrs/rprs);
    dh = -0.5*ri*((rs2+r2)*rsr2mrs2i*r2mrs2i+logr );
}

inline void h6Ser ( BigReal r, BigReal r2, BigReal ri,//0%, one atom within other
BigReal rc, BigReal r0, BigReal rs, BigReal & h ) {
  h = 0;
}
inline void dh6Ser ( BigReal r, BigReal r2, BigReal ri,//0%, one atom within other
BigReal rc, BigReal r0, BigReal rs, BigReal & dh ) {
  dh = 0;
}

inline void CalcHSer ( BigReal r, BigReal r2, BigReal ri,
BigReal rc, BigReal r0, BigReal rs, BigReal & h, int & d) {
  if( r <= rc - rs && r > 4*rs ) {//II
    h2Ser(r,r2,ri,rc,r0,rs,h); d = 2;
  } else if (r <= rc + rs && r > rc - rs) {//I
    h1Ser(r,r2,ri,rc,r0,rs,h); d = 1;
  } else if (r > rc + rs) {//0
    h0Ser(r,r2,ri,rc,r0,rs,h); d = 0;
  } else if( r <= 4*rs && r > r0 + rs ) {//III
    h3Ser(r,r2,ri,rc,r0,rs,h); d = 3;
  } else if ( r <= r0 + rs && r > (r0>rs?r0-rs:rs-r0) ) {//IV
    h4Ser(r,r2,ri,rc,r0,rs,h); d = 4;
  } else if (r0 < rs ) {//V
    h5Ser(r,r2,ri,rc,r0,rs,h); d = 5;
  } else {//VI
    h6Ser(r,r2,ri,rc,r0,rs,h); d = 6;
  }
}
inline void CalcDHSer ( BigReal r, BigReal r2, BigReal ri,
BigReal rc, BigReal r0, BigReal rs, BigReal & dh, int & d) {
  if( r <= rc - rs && r > 4*rs ) {//II
    dh2Ser(r,r2,ri,rc,r0,rs,dh); d = 2;
  } else if (r <= rc + rs && r > rc - rs) {//I
    dh1Ser(r,r2,ri,rc,r0,rs,dh); d = 1;
  } else if (r > rc + rs) {//0
    dh0Ser(r,r2,ri,rc,r0,rs,dh); d = 0;
  } else if( r <= 4*rs && r > r0 + rs ) {//III
    dh3Ser(r,r2,ri,rc,r0,rs,dh); d = 3;
  } else if ( r <= r0 + rs && r > (r0>rs?r0-rs:rs-r0) ) {//IV
    dh4Ser(r,r2,ri,rc,r0,rs,dh); d = 4;
  } else if (r0 < rs ) {//V
    dh5Ser(r,r2,ri,rc,r0,rs,dh); d = 5;
  } else {//VI
    dh6Ser(r,r2,ri,rc,r0,rs,dh); d = 6;
  }
}
inline void CalcHPairSer (
  BigReal r,//distance
  BigReal r2,//distance squared
  BigReal ri,//inverse distance
  BigReal rc,//cutoff
  BigReal ri0,
  BigReal rjs,
  BigReal rj0,
  BigReal ris,
  int & dij,//domain 1
  int & dji,//domain 2
  BigReal & hij,//output
  BigReal & hji//output
) {
  CalcHSer(r,r2,ri,rc,ri0,rjs,hij,dij);//hij
  CalcHSer(r,r2,ri,rc,rj0,ris,hji,dji);//hji
}
inline void CalcDHPairSer (
  BigReal r,//distance
  BigReal r2,
  BigReal ri,
  BigReal rc,//cutoff
  BigReal ri0,
  BigReal rjs,
  BigReal rj0,
  BigReal ris,
  int & dij,//domain 1
  int & dji,//domain 2
  BigReal & dhij,
  BigReal & dhji
) {
  CalcDHSer(r,r2,ri,rc,ri0,rjs,dhij,dij);//hij
  CalcDHSer(r,r2,ri,rc,rj0,ris,dhji,dji);//hji
}


/*******************************************************************************
********************************************************************************
***********  Phase 2 Inner Loop   **********************************************
********************************************************************************
*******************************************************************************/

/*
 * Calculate coulomb energy and force for single pair of atoms
 */
inline void Calc_Coul_PairSer(
  BigReal r_i,
  BigReal qiqj,
  BigReal epsilon_p_i,
  int exclij,
  BigReal scale14,
  BigReal & coulE,
  BigReal & ddrCoulE
) {
  if (exclij != EXCHCK_FULL) {//not excluded
    //calculate Coulomb Energy
    coulE = -qiqj*epsilon_p_i*r_i;

    //calculate Coulomb Force
    if (exclij == EXCHCK_MOD)
      coulE *= scale14;
    ddrCoulE = -r_i*coulE;
  } else {
    coulE = 0;
    ddrCoulE = 0;
  }
}

/*
 * Calculate GB Energy, GB dEdr force
 * also output intermediate values used in dEda
 */
inline void Calc_dEdr_PairSer(//no longer does i==j
  BigReal r,
  BigReal r2,
  BigReal qiqj,
  BigReal ai,
  BigReal aj,
  BigReal kappa,
  BigReal epsilon_p_i,
  BigReal epsilon_s_i,
  BigReal & aiaj,
  BigReal & expr2aiaj4,
  BigReal & fij,
  BigReal & f_i,
  BigReal & expkappa,
  BigReal & Dij,
  BigReal & gbE,   //return
  BigReal & ddrGbE //return
) {
  //allocate local variables
  BigReal aiaj4,ddrDij,ddrf_i,ddrfij;

  //calculate GB energy
  aiaj = ai*aj;
  aiaj4 = 4*aiaj;
  expr2aiaj4 = exp(-r2/aiaj4);
  fij = sqrt(r2+aiaj*expr2aiaj4);
  f_i = 1/fij;
  if (kappa > 0)
    expkappa = exp(-kappa*fij);
  else
    expkappa = 1.0;
  Dij = epsilon_p_i - expkappa*epsilon_s_i;//dielectric term
  gbE = qiqj*Dij*f_i;

  //calculate energy derivatives
  ddrfij = r*f_i*(1 - 0.25*expr2aiaj4);
  ddrf_i = -ddrfij*f_i*f_i;
  ddrDij = kappa*expkappa*ddrfij*epsilon_s_i;
  ddrGbE = qiqj*(ddrDij*f_i+Dij*ddrf_i);
}

/*
 * Calculate summation element of dEda array
 * must calculate dEdr previously to retreive intermediate values
 */
inline void Calc_dEda_PairSer(
  BigReal r2,
  BigReal ai,
  BigReal aj,
  BigReal qiqj,
  BigReal kappa,
  BigReal aiaj,
  BigReal expkappa,
  BigReal expr2aiaj4,
  BigReal fij,
  BigReal f_i,
  BigReal Dij,
  BigReal epsilon_s_i,
  BigReal & dEdai,//return
  BigReal & dEdaj //return
) {

  BigReal tmp_dEda = 0.5*qiqj*f_i*f_i
                      *(kappa*epsilon_s_i*expkappa-Dij*f_i)
                      *(aiaj+0.25*r2)*expr2aiaj4;//0
  dEdai = tmp_dEda/ai;
  dEdaj = tmp_dEda/aj;
}

/*
 * Calculate Coulomb and GB interaction and dEda element
 * for a pair of atoms
 */
inline void Phase2_PairSer(//doesn't do self energies

//input values
  BigReal r,
  BigReal r2,
  BigReal r_i,
  BigReal qiqj,
  BigReal ai,
  BigReal aj,
  BigReal epsilon_p_i,
  BigReal epsilon_s_i,
  BigReal kappa,
  int exclij,
  BigReal scale14,
  int doSlow,

//return values
  BigReal & coulEij,
  BigReal & ddrCoulEij,
  BigReal & gbEij,
  BigReal & ddrGbEij,
  BigReal & dEdai,
  BigReal & dEdaj
) {

  //calculate Coulomb energy and force
  //Calc_Coul_Pair(r_i,qiqj,epsilon_p_i,exclij,scale14,coulEij,ddrCoulEij);
  coulEij = 0;
  ddrCoulEij = 0;

  //calculate GB energy and force
  BigReal aiaj,expr2aiaj4,fij,f_i,expkappa,Dij;
  Calc_dEdr_PairSer(r,r2,qiqj,ai,aj,kappa,epsilon_p_i,epsilon_s_i,
      aiaj,expr2aiaj4,fij,f_i,expkappa,Dij,gbEij,ddrGbEij);

  //calculate dEda
  if (doSlow) {
    Calc_dEda_PairSer(r2,ai,aj,qiqj,kappa,aiaj,expkappa,expr2aiaj4,
             fij,f_i,Dij,epsilon_s_i,dEdai,dEdaj);
  }
}

