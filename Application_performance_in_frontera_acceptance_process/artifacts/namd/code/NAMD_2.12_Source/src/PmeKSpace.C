/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "PmeKSpace.h"
#include <math.h>
#include <stdlib.h>

#ifndef _NO_ALLOCA_H
#include <alloca.h>
#endif

#include "PmeBase.inl"
#include "SimParameters.h" 
#include "Node.h"  
#include "ComputeMoaMgr.decl.h"  

#define ALLOCA(TYPE,NAME,SIZE) TYPE *NAME = (TYPE *) alloca((SIZE)*sizeof(TYPE))

static void dftmod(double *bsp_mod, double *bsp_arr, int nfft) {
  int j, k;
  double twopi, arg, sum1, sum2;
  double infft = 1.0/nfft;
/* Computes the modulus of the discrete fourier transform of bsp_arr, */
/*  storing it into bsp_mod */
  twopi =  2.0 * M_PI;

  for (k = 0; k <nfft; ++k) {
    sum1 = 0.;
    sum2 = 0.;
    for (j = 0; j < nfft; ++j) {
      arg = twopi * k * j * infft;
      sum1 += bsp_arr[j] * cos(arg);
      sum2 += bsp_arr[j] * sin(arg);
    }
    bsp_mod[k] = sum1*sum1 + sum2*sum2;
  }
}

void compute_b_moduli(double *bm, int K, int order) {
  int i;
  double fr[3];

  double *M = new double[3*order];
  double *dM = new double[3*order];
  double *scratch = new double[K];

  fr[0]=fr[1]=fr[2]=0.0;
  compute_b_spline(fr,M,dM,order);  
  for (i=0; i<order; i++) bm[i] = M[i];
  for (i=order; i<K; i++) bm[i] = 0.0;
  dftmod(scratch, bm, K);
  for (i=0; i<K; i++) bm[i] = 1.0/scratch[i];


  delete [] scratch;
  delete [] dM;
  delete [] M;
}

PmeKSpace::PmeKSpace(PmeGrid grid, int K2_start, int K2_end, int K3_start, int K3_end) 
  : myGrid(grid), k2_start(K2_start), k2_end(K2_end),
	k3_start(K3_start), k3_end(K3_end) {
  int K1, K2, K3, order;
  K1=myGrid.K1; K2=myGrid.K2, K3=myGrid.K3; order=myGrid.order;

  bm1 = new double[K1];
  bm2 = new double[K2];
  bm3 = new double[K3];

  exp1 = new double[K1/2 + 1];
  exp2 = new double[K2/2 + 1];
  exp3 = new double[K3/2 + 1];

  compute_b_moduli(bm1, K1, order);
  compute_b_moduli(bm2, K2, order);
  compute_b_moduli(bm3, K3, order);

}

#ifdef OPENATOM_VERSION
PmeKSpace::PmeKSpace(PmeGrid grid, int K2_start, int K2_end, int K3_start, int K3_end, CProxy_ComputeMoaMgr moaProxy) 
  : myGrid(grid), k2_start(K2_start), k2_end(K2_end),
	k3_start(K3_start), k3_end(K3_end) {
  int K1, K2, K3, order;
  K1=myGrid.K1; K2=myGrid.K2, K3=myGrid.K3; order=myGrid.order;

  SimParameters *simParams = Node::Object()->simParameters;

  bm1 = new double[K1];
  bm2 = new double[K2];
  bm3 = new double[K3];

  exp1 = new double[K1/2 + 1];
  exp2 = new double[K2/2 + 1];
  exp3 = new double[K3/2 + 1];

  compute_b_moduli(bm1, K1, order);
  compute_b_moduli(bm2, K2, order);
  compute_b_moduli(bm3, K3, order);

  if ( simParams->openatomOn ) {
    CkPrintf("######################################################\n");
    CkPrintf("Entering recvBmX loop on processor %d \n", CkMyPe() );
    int i;
    for ( i=0; i<=K1; ++i) {
    CkPrintf("bm1 value pre-send %d = %d \n", i, bm1[i] );
    }
    for ( i=0; i<=K1; ++i) {
    CkPrintf("bm1 reference pre-send %d = %d \n", i, &bm1[i] );
    }
    CkPrintf("bm1: %d \n", *bm1);    
    moaProxy[CkMyPe()].recvB(K2_start, K2_end, K3_start, K3_end, K1, K2, K3, bm1, bm2, bm3, order);
//    qmmmcp->recvBm1(bm1, K1, order);
//    qmmmcp->recvBm2(bm2, K2, order);
//    qmmmcp->recvBm3(bm3, K3, order);
  }
    

}
#endif // OPENATOM_VERSION

PmeKSpace::~PmeKSpace() {
  delete [] bm1;
  delete [] bm2;
  delete [] bm3;
  
  delete [] exp1;
  delete [] exp2;
  delete [] exp3;
}

void PmeKSpace::compute_energy_orthogonal_subset(float q_arr[], double *recips, double *partialVirial, double *partialEnergy, int k1from, int k1to){

    double energy = 0.0;
    double v0 = 0.;
    double v1 = 0.;
    double v2 = 0.;
    double v3 = 0.;
    double v4 = 0.;
    double v5 = 0.;
    
    int k1, k2, k3;
    int K1, K2, K3;
    K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3;
    
    double recipx = recips[0];
    double recipy = recips[1];
    double recipz = recips[2];

    int ind = k1from*(k2_end-k2_start)*(k3_end-k3_start)*2;
        
    for ( k1=k1from; k1<=k1to; ++k1 ) {        
        double m1, m11, b1, xp1;
        b1 = bm1[k1];
        int k1_s = k1<=K1/2 ? k1 : k1-K1;
        m1 = k1_s*recipx;
        m11 = m1*m1;
        xp1 = i_pi_volume*exp1[abs(k1_s)];
        for ( k2=k2_start; k2<k2_end; ++k2 ) {
            double m2, m22, b1b2, xp2;
            b1b2 = b1*bm2[k2];
            int k2_s = k2<=K2/2 ? k2 : k2-K2;
            m2 = k2_s*recipy;
            m22 = m2*m2;
            xp2 = exp2[abs(k2_s)]*xp1;

            k3 = k3_start;
            if (k1==0 && k2==0 && k3==0) {
              q_arr[ind++] = 0.0;
              q_arr[ind++] = 0.0;
              ++k3;
            }
            for (; k3<k3_end; ++k3 ) {
              double m3, m33, xp3, msq, imsq, vir, fac;
              double theta3, theta, q2, qr, qc, C;
              theta3 = bm3[k3] *b1b2;
              m3 = k3*recipz;
              m33 = m3*m3;
              xp3 = exp3[k3];
              qr = q_arr[ind]; qc=q_arr[ind+1];
              q2 = 2*(qr*qr + qc*qc)*theta3;
              if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
              msq = m11 + m22 + m33;
              imsq = 1.0/msq;
              C = xp2*xp3*imsq;
              theta = theta3*C;
              q_arr[ind] *= theta;
              q_arr[ind+1] *= theta;
              vir = -2*(piob+imsq);
              fac = q2*C;
              energy += fac;
              v0 += fac*(1.0+vir*m11);
              v1 += fac*vir*m1*m2;
              v2 += fac*vir*m1*m3;
              v3 += fac*(1.0+vir*m22);
              v4 += fac*vir*m2*m3;
              v5 += fac*(1.0+vir*m33);
              ind += 2;
            }
        }
    }
    
    *partialEnergy = 0.5*energy;
    partialVirial[0] = 0.5*v0;
    partialVirial[1] = 0.5*v1;
    partialVirial[2] = 0.5*v2;
    partialVirial[3] = 0.5*v3;
    partialVirial[4] = 0.5*v4;
    partialVirial[5] = 0.5*v5;
}
static inline void compute_energy_orthogonal_ckloop(int first, int last, void *result, int paraNum, void *param){
  for ( int i = first; i <= last; ++i ) {
    void **params = (void **)param;
    PmeKSpace *kspace = (PmeKSpace *)params[0];
    float *q_arr = (float *)params[1];
    double *recips = (double *)params[2];
    double *partialEnergy = (double *)params[3];
    double *partialVirial = (double *)params[4];
    int *unitDist = (int *)params[5];
    
    int unit = unitDist[0];
    int remains = unitDist[1];
    int k1from, k1to;
    if(i<remains){
        k1from = i*(unit+1);
        k1to = k1from+unit;
    }else{
        k1from = remains*(unit+1)+(i-remains)*unit;
        k1to = k1from+unit-1;
    }
    double *pEnergy = partialEnergy+i;
    double *pVirial = partialVirial+i*6;
    kspace->compute_energy_orthogonal_subset(q_arr, recips, pVirial, pEnergy, k1from, k1to);
  }
}

double PmeKSpace::compute_energy_orthogonal_helper(float *q_arr, const Lattice &lattice, double ewald, double *virial) {
  double energy = 0.0;
  double v0 = 0.;
  double v1 = 0.;
  double v2 = 0.;
  double v3 = 0.;
  double v4 = 0.;
  double v5 = 0.;

  int n;
  int k1, k2, k3, ind;
  int K1, K2, K3;

  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3;

  i_pi_volume = 1.0/(M_PI * lattice.volume());
  piob = M_PI/ewald;
  piob *= piob;


    double recipx = lattice.a_r().x;
    double recipy = lattice.b_r().y;
    double recipz = lattice.c_r().z;
        
    init_exp(exp1, K1, 0, K1, recipx);
    init_exp(exp2, K2, k2_start, k2_end, recipy);
    init_exp(exp3, K3, k3_start, k3_end, recipz);

    double recips[] = {recipx, recipy, recipz};
    int NPARTS=CmiMyNodeSize(); //this controls the granularity of loop parallelism
    int maxParts = ( K1 * ( k2_end - k2_start ) * ( k3_end - k3_start ) + 127 ) / 128;
    if ( NPARTS >  maxParts ) NPARTS = maxParts;
    if ( NPARTS >  K1 ) NPARTS = K1; 
    ALLOCA(double, partialEnergy, NPARTS);
    ALLOCA(double, partialVirial, 6*NPARTS);
    int unitDist[] = {K1/NPARTS, K1%NPARTS};
    
    //parallelize the following loop using CkLoop
    void *params[] = {this, q_arr, recips, partialEnergy, partialVirial, unitDist};

#if     CMK_SMP && USE_CKLOOP
    CkLoop_Parallelize(compute_energy_orthogonal_ckloop, 6, (void *)params, NPARTS, 0, NPARTS-1);
#endif
/*  
    //The transformed loop used to compute energy
    int unit = K1/NPARTS;
    int remains = K1%NPARTS;  
    for(int i=0; i<NPARTS; i++){
        int k1from, k1to;
        if(i<remains){
            k1from = i*(unit+1);
            k1to = k1from+unit;
        }else{
            k1from = remains*(unit+1)+(i-remains)*unit;
            k1to = k1from+unit-1;
        }
        double *pEnergy = partialEnergy+i;
        double *pVirial = partialVirial+i*6;
        compute_energy_orthogonal_subset(q_arr, recips, pVirial, pEnergy, k1from, k1to);
    }
*/    
    
    for(int i=0; i<NPARTS; i++){
        v0 += partialVirial[i*6+0];
        v1 += partialVirial[i*6+1];
        v2 += partialVirial[i*6+2];
        v3 += partialVirial[i*6+3];
        v4 += partialVirial[i*6+4];
        v5 += partialVirial[i*6+5];
        energy += partialEnergy[i];
    }
    
    virial[0] = v0;
    virial[1] = v1;
    virial[2] = v2;
    virial[3] = v3;
    virial[4] = v4;
    virial[5] = v5;
    return energy;
}

double PmeKSpace::compute_energy(float *q_arr, const Lattice &lattice, double ewald, double *virial, int useCkLoop) {
  double energy = 0.0;
  double v0 = 0.;
  double v1 = 0.;
  double v2 = 0.;
  double v3 = 0.;
  double v4 = 0.;
  double v5 = 0.;

  int n;
  int k1, k2, k3, ind;
  int K1, K2, K3;

  K1=myGrid.K1; K2=myGrid.K2; K3=myGrid.K3;

  i_pi_volume = 1.0/(M_PI * lattice.volume());
  piob = M_PI/ewald;
  piob *= piob;

  if ( lattice.orthogonal() ) {
  // if ( 0 ) { // JCP FOR TESTING
    //This branch is the usual call path.
#if     CMK_SMP && USE_CKLOOP
    if ( useCkLoop ) {
        return compute_energy_orthogonal_helper(q_arr, lattice, ewald, virial);
    }
#endif    
    double recipx = lattice.a_r().x;
    double recipy = lattice.b_r().y;
    double recipz = lattice.c_r().z;
        
    init_exp(exp1, K1, 0, K1, recipx);
    init_exp(exp2, K2, k2_start, k2_end, recipy);
    init_exp(exp3, K3, k3_start, k3_end, recipz);

    ind = 0;
    for ( k1=0; k1<K1; ++k1 ) {
      double m1, m11, b1, xp1;
      b1 = bm1[k1];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      m1 = k1_s*recipx;
      m11 = m1*m1;
      xp1 = i_pi_volume*exp1[abs(k1_s)];
      for ( k2=k2_start; k2<k2_end; ++k2 ) {
        double m2, m22, b1b2, xp2;
        b1b2 = b1*bm2[k2];
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        m2 = k2_s*recipy;
        m22 = m2*m2;
        xp2 = exp2[abs(k2_s)]*xp1;
        k3 = k3_start;
        if ( k1==0 && k2==0 && k3==0 ) {
          q_arr[ind++] = 0.0;
          q_arr[ind++] = 0.0;
          ++k3;
        }
        for ( ; k3<k3_end; ++k3 ) {
          double m3, m33, xp3, msq, imsq, vir, fac;
          double theta3, theta, q2, qr, qc, C;
          theta3 = bm3[k3] *b1b2;
          m3 = k3*recipz;
          m33 = m3*m3;
          xp3 = exp3[k3];
          qr = q_arr[ind]; qc=q_arr[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          msq = m11 + m22 + m33;
          imsq = 1.0/msq;
          C = xp2*xp3*imsq;
          theta = theta3*C;
          q_arr[ind] *= theta;
          q_arr[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          v0 += fac*(1.0+vir*m11);
          v1 += fac*vir*m1*m2;
          v2 += fac*vir*m1*m3;
          v3 += fac*(1.0+vir*m22);
          v4 += fac*vir*m2*m3;
          v5 += fac*(1.0+vir*m33);
          ind += 2;
        }
      }
    }
    
  } else if ( cross(lattice.a(),lattice.b()).unit() == lattice.c().unit() ) {
  // } else if ( 0 ) { // JCP FOR TESTING
    Vector recip1 = lattice.a_r();
    Vector recip2 = lattice.b_r();
    Vector recip3 = lattice.c_r();
    double recip3_x = recip3.x;
    double recip3_y = recip3.y;
    double recip3_z = recip3.z;
    init_exp(exp3, K3, k3_start, k3_end, recip3.length());

    ind = 0;
    for ( k1=0; k1<K1; ++k1 ) {
      double b1; Vector m1;
      b1 = bm1[k1];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      m1 = k1_s*recip1;
      // xp1 = i_pi_volume*exp1[abs(k1_s)];
      for ( k2=k2_start; k2<k2_end; ++k2 ) {
        double xp2, b1b2, m2_x, m2_y, m2_z;
        b1b2 = b1*bm2[k2];
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        m2_x = m1.x + k2_s*recip2.x;
        m2_y = m1.y + k2_s*recip2.y;
        m2_z = m1.z + k2_s*recip2.z;
        // xp2 = exp2[abs(k2_s)]*xp1;
        xp2 = i_pi_volume*exp(-piob*(m2_x*m2_x+m2_y*m2_y+m2_z*m2_z));
        k3 = k3_start;
        if ( k1==0 && k2==0 && k3==0 ) {
          q_arr[ind++] = 0.0;
          q_arr[ind++] = 0.0;
          ++k3;
        }
        for ( ; k3<k3_end; ++k3 ) {
          double xp3, msq, imsq, vir, fac;
          double theta3, theta, q2, qr, qc, C;
          double m_x, m_y, m_z;
          theta3 = bm3[k3] *b1b2;
          m_x = m2_x + k3*recip3_x;
          m_y = m2_y + k3*recip3_y;
          m_z = m2_z + k3*recip3_z;
          msq = m_x*m_x + m_y*m_y + m_z*m_z;
          xp3 = exp3[k3];
          qr = q_arr[ind]; qc=q_arr[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          imsq = 1.0/msq;
          C = xp2*xp3*imsq;
          theta = theta3*C;
          q_arr[ind] *= theta;
          q_arr[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          v0 += fac*(1.0+vir*m_x*m_x);
          v1 += fac*vir*m_x*m_y;
          v2 += fac*vir*m_x*m_z;
          v3 += fac*(1.0+vir*m_y*m_y);
          v4 += fac*vir*m_y*m_z;
          v5 += fac*(1.0+vir*m_z*m_z);
          ind += 2;
        }
      }
    }

  } else {
    Vector recip1 = lattice.a_r();
    Vector recip2 = lattice.b_r();
    Vector recip3 = lattice.c_r();
    double recip3_x = recip3.x;
    double recip3_y = recip3.y;
    double recip3_z = recip3.z;

    ind = 0;
    for ( k1=0; k1<K1; ++k1 ) {
      double b1; Vector m1;
      b1 = bm1[k1];
      int k1_s = k1<=K1/2 ? k1 : k1-K1;
      m1 = k1_s*recip1;
      // xp1 = i_pi_volume*exp1[abs(k1_s)];
      for ( k2=k2_start; k2<k2_end; ++k2 ) {
        double b1b2, m2_x, m2_y, m2_z;
        b1b2 = b1*bm2[k2];
        int k2_s = k2<=K2/2 ? k2 : k2-K2;
        m2_x = m1.x + k2_s*recip2.x;
        m2_y = m1.y + k2_s*recip2.y;
        m2_z = m1.z + k2_s*recip2.z;
        // xp2 = exp2[abs(k2_s)]*xp1;
        k3 = k3_start;
        if ( k1==0 && k2==0 && k3==0 ) {
          q_arr[ind++] = 0.0;
          q_arr[ind++] = 0.0;
          ++k3;
        }
        for ( ; k3<k3_end; ++k3 ) {
          double xp3, msq, imsq, vir, fac;
          double theta3, theta, q2, qr, qc, C;
          double m_x, m_y, m_z;
          theta3 = bm3[k3] *b1b2;
          m_x = m2_x + k3*recip3_x;
          m_y = m2_y + k3*recip3_y;
          m_z = m2_z + k3*recip3_z;
          msq = m_x*m_x + m_y*m_y + m_z*m_z;
          // xp3 = exp3[k3];
          xp3 = i_pi_volume*exp(-piob*msq);
          qr = q_arr[ind]; qc=q_arr[ind+1];
          q2 = 2*(qr*qr + qc*qc)*theta3;
          if ( (k3 == 0) || ( k3 == K3/2 && ! (K3 & 1) ) ) q2 *= 0.5;
          imsq = 1.0/msq;
          C = xp3*imsq;
          theta = theta3*C;
          q_arr[ind] *= theta;
          q_arr[ind+1] *= theta;
          vir = -2*(piob+imsq);
          fac = q2*C;
          energy += fac;
          v0 += fac*(1.0+vir*m_x*m_x);
          v1 += fac*vir*m_x*m_y;
          v2 += fac*vir*m_x*m_z;
          v3 += fac*(1.0+vir*m_y*m_y);
          v4 += fac*vir*m_y*m_z;
          v5 += fac*(1.0+vir*m_z*m_z);
          ind += 2;
        }
      }
    }

  }

  virial[0] = 0.5 * v0;
  virial[1] = 0.5 * v1;
  virial[2] = 0.5 * v2;
  virial[3] = 0.5 * v3;
  virial[4] = 0.5 * v4;
  virial[5] = 0.5 * v5;
  return 0.5*energy;
}


void PmeKSpace::init_exp(double *xp, int K, int k_start, int k_end, double recip) {
  int i;
  double fac;
  fac = -piob*recip*recip;
  int i_start = k_start;
  int i_end = k_end;
  if ( k_start > K/2 ) {
    i_start = K - k_end + 1;
    i_end = K - k_start + 1;
  } else if ( k_end > K/2 ) {
    i_end = K/2 + 1;
    i_start = K - k_end + 1;
    if ( k_start < i_start ) i_start = k_start;
  }

  for (i=i_start; i < i_end; i++)
    xp[i] = exp(fac*i*i);
} 
