
/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

//#define BENCHMARK
#define CHECK_PRIORITIES 0
//#define PRINT_COMP

//approximate flop equivalence of transcendental functs
#ifdef BENCHMARK
#define SQRT_FLOPS 15
#define LOG_FLOPS 20
#define EXP_FLOPS 20
#define DIV_FLOPS 1
#endif

#include "ComputeNonbondedUtil.h"
#include "ComputeGBIS.inl"
#include "time.h"

/*******************************************************************************
  When GBIS is active, the method calcGBIS is called from a nonbonded
  compute object ::doForce() three times each timestep with
  gbisPhase = 1, 2, 3.

  At the beginning of the first phase, the pairlists are generated which are
  used for all three phases; these pairlists are re-generated every timestep.
  4 mutually exclusive pairlists are used, in an effort to accelerate the
  iteration over pairs.
  pairlist 0: atoms fall within interaction domain 2 (largest group)
  pairlist 1: atoms fall within interaction domain 1 (next largest)
  pairlist 2: atoms are within the phase 1,3 cutoff (all the rest)
  pairlist 3: atoms are within the phase 2 cutoff (generally longer than previous cutoff)

  These calculations generate nonbonded forces (and associated energies) which
  are calculated in addition to Coulomb and van der Waals forces.
*******************************************************************************/

/*Searches all pair of atoms to create
largest pairlist lasting all cycle*/
inline void pairlistFromAll(
  nonbonded *params,
  GBISParamStruct *gbisParams,
  int minIg,
  int strideIg,
  int maxI
) {
  const BigReal offset_x = params->offset.x;
  const BigReal offset_y = params->offset.y;
  const BigReal offset_z = params->offset.z;

  int unique = (gbisParams->numPatches == 1) ? 1 : 0;

  int maxPairs = params->numAtoms[1];
  int seq = gbisParams->sequence;
  int cid = gbisParams->cid;
  int pe = CkMyPe();
  int gbisPhase = gbisParams->gbisPhase;

  //determine long and short cutoff
  float fsMax = FS_MAX; // gbisParams->fsMax;
  float a_cut = gbisParams->a_cut - fsMax;
  float a_cut_ps2 = (a_cut+fsMax)*(a_cut+fsMax);
  float r_cut = gbisParams->cutoff;
  float a_cut2 = a_cut*a_cut;
  float r_cut2 = r_cut*r_cut;
  int s = 0;// 0 is short cutoff list
  int l = 1;// 1 is long-short cutoff list
#ifdef BENCHMARK
  double t1 = 1.0*clock()/CLOCKS_PER_SEC;
  int nops = 0;
#endif
  Position ri, rj;
  Position ngri, ngrj;
  float r, dr, r2;
  float rhoi0, rhois, rhoj0, rhojs, rhois2, rhojs2, rhois2_16;
  float ami2, amj2, api2, apj2;
  int numGBISPairlists = 4;

  float max_gcut2  = r_cut;
  if (a_cut+fsMax > r_cut)
    max_gcut2 = a_cut+fsMax;
  max_gcut2 += 2.0*gbisParams->maxGroupRadius;
  //CkPrintf("max_cut = %f, %f\n",max_gcut2,gbisParams->maxGroupRadius);
  max_gcut2 *= max_gcut2;

  int maxGroupPairs = params->numAtoms[1];
  short *groupPairs = new short[maxGroupPairs];/*delete*/

  //foreach nonbonded group i
  for (int ngi = minIg; ngi < maxI; /*ngi updated at loop bottom*/ ) {
    int numGroupPairs = 0;
    ngri = params->p[0][ngi].position;
    ngri.x += offset_x;
    ngri.y += offset_y;
    ngri.z += offset_z;

    //find close j-groups; include i-group
    for (int ngj = unique*(ngi+params->p[0][ngi].nonbondedGroupSize); ngj < params->numAtoms[1]; ngj+=params->p[1][ngj].nonbondedGroupSize) {
      
      ngrj = params->p[1][ngj].position;
      dr = ngri.x - ngrj.x;
      r2 = dr*dr;
      dr = ngri.y - ngrj.y;
      r2 += dr*dr;
      dr = ngri.z - ngrj.z;
      r2 += dr*dr;
      if (r2 < max_gcut2) {
        groupPairs[numGroupPairs++] = ngj;
      }
    }

    //for all i in i-group
    int iGroupSize = params->p[0][ngi].nonbondedGroupSize;
    for (int i=ngi; i<ngi+iGroupSize; i++) {
      //CkPrintf("\tFORALL i=%05i\n",params->pExt[0][i].id);
      //extend pairlists
      int *size = new int[numGBISPairlists];
      plint **pairs = new plint*[numGBISPairlists];
      for (int k = 0; k < numGBISPairlists; k++) {
        size[k] = 0;
        pairs[k] = gbisParams->gbisStepPairlists[k]->
                             newlist(maxPairs);
      }

      //load i values
      ri = params->p[0][i].position;
      ri.x += offset_x;
      ri.y += offset_y;
      ri.z += offset_z;
      rhois = gbisParams->intRad[0][2*i+1];
      api2 = a_cut + rhois;
      api2 *= api2;
      ami2 = a_cut - rhois;
      ami2 *= ami2;
      rhois2 = rhois*rhois;
      rhois2_16 = 16.0*rhois2;

      //for all other i's in i-group
      if (unique)
      for (int j = i+1; j < ngi+iGroupSize; j++) {

#ifdef BENCHMARK
          nops ++;
#endif
          rj = params->p[1][j].position;
          dr = ri.x - rj.x;
          r2 = dr*dr;
          dr = ri.y - rj.y;
          r2 += dr*dr;
          dr = ri.z - rj.z;
          r2 += dr*dr;

          rhojs = gbisParams->intRad[1][2*j+1];
          rhojs2 = rhojs*rhojs;
          //r = sqrt(r2);
          amj2 = a_cut - rhojs;
          amj2 *= amj2;

          apj2 = a_cut + rhojs;
          apj2 *= apj2;
          //CkPrintf("IPAIR %5i %5i %5i\n",0*gbisParams->cid,params->pExt[0][i].id,params->pExt[1][j].id);

          if (  r2 < ami2 && r2 < amj2 &&
            r2 > rhois2_16 && r2 > 16.0*rhojs2 ) {
            //belongs to 22
            pairs[0][size[0]++] = j;
          } else if ( r2 < api2 && r2 < apj2 &&
            r2 > ami2 && r2 > amj2 ) {
            //belongs to 11
            pairs[1][size[1]++] = j;
          } else if ( r2 < a_cut_ps2 ) {
            //belongs to other a_cut
            pairs[2][size[2]++] = j;
          } else if ( r2 < r_cut2 ) {
            //belongs to r_cut and (r_cut > a_cut)
            pairs[3][size[3]++] = j;
          }
      }

      //foreach j-group
      for (int g = 0; g < numGroupPairs; g++) {
        int ngj = groupPairs[g];
        int jGroupSize = params->p[1][ngj].nonbondedGroupSize;
        //for each j in j-group
        int maxJ = ngj+jGroupSize;
        for (int j = ngj; j < maxJ; j++) {

          //CkPrintf("for j-atom%i\n",j);
#ifdef BENCHMARK
          nops ++;
#endif
          rj = params->p[1][j].position;
          dr = ri.x - rj.x;
          r2 = dr*dr;
          dr = ri.y - rj.y;
          r2 += dr*dr;
          dr = ri.z - rj.z;
          r2 += dr*dr;
          //CkPrintf("PAIR %5i %5i %5i\n",0*gbisParams->cid,params->pExt[0][i].id,params->pExt[1][j].id);

          rhojs = gbisParams->intRad[1][2*j+1];
          rhojs2 = rhojs*rhojs;
          //r = sqrt(r2);
          amj2 = a_cut - rhojs;
          amj2 *= amj2;

          apj2 = a_cut + rhojs;
          apj2 *= apj2;

          if (  r2 < ami2 && r2 < amj2 &&
            r2 > rhois2_16 && r2 > 16.0*rhojs2 ) {
            //belongs to 22
            //CkPrintf("Pair22 %06i %06i\n",params->pExt[0][i].id,params->pExt[1][j].id);
            pairs[0][size[0]++] = j;
          } else if ( r2 < api2 && r2 < apj2 &&
            r2 > ami2 && r2 > amj2 ) {
            //CkPrintf("Pair11 %06i %06i\n",params->pExt[0][i].id,params->pExt[1][j].id);
            //belongs to 11
            pairs[1][size[1]++] = j;
          } else if ( r2 < a_cut_ps2 ) {
            //CkPrintf("PairAA %06i %06i\n",params->pExt[0][i].id,params->pExt[1][j].id);
            //belongs to other a_cut
            pairs[2][size[2]++] = j;
          } else if ( r2 < r_cut2 ) {
            //CkPrintf("PairRR %06i %06i\n",params->pExt[0][i].id,params->pExt[1][j].id);
            //belongs to r_cut and (r_cut > a_cut)
            pairs[3][size[3]++] = j;
          }
        }//end j inner loop
      }//end j-group loop
    for (int k = 0; k < numGBISPairlists; k++) {
      gbisParams->gbisStepPairlists[k]->newsize(size[k]);
    }
    delete[] size;
    delete[] pairs;
    }//end i all atom loop

    //jump to next nbg for round-robin
    //same iteration patter followed in all three gbisPhases below
    for (int s = 0; s < strideIg; s++) {
      ngi+=params->p[0][ngi].nonbondedGroupSize;
    }
  }//end i-group loop
    delete[] groupPairs;
  for (int k = 0; k < numGBISPairlists; k++)
    gbisParams->gbisStepPairlists[k]->reset();

#ifdef BENCHMARK
  double t2 = 1.0*clock()/CLOCKS_PER_SEC;
  CkPrintf("PHASELST: %8.3f ms @ %8.3f ns/iter for %i iter\n",1000.0*(t2-t1), 1000000000.0*(t2-t1)/nops, nops);
#endif
}//end pairlist method

/*******************************************************************************
 * Calculate GBIS 3 Phases
*******************************************************************************/
void ComputeNonbondedUtil::calcGBIS(nonbonded *params, GBISParamStruct *gbisParams) {
//CkPrintf("SEQ%03i, P%i, CID%05i(%02i,%02i) ENTER\n",gbisParams->sequence,gbisParams->gbisPhase,gbisParams->cid,gbisParams->patchID[0],gbisParams->patchID[1]);
#if CHECK_PRIORITIES
CkPrintf("PE%i, S%09i, P%i\n",CkMyPe(),gbisParams->sequence,gbisParams->gbisPhase);
#endif
  if (params->numAtoms[0] == 0 || params->numAtoms[1] == 0) return;

  const BigReal offset_x = params->offset.x;
  const BigReal offset_y = params->offset.y;
  const BigReal offset_z = params->offset.z;

  int partSize = params->numAtoms[0] / params->numParts;
  int minIg = 0;
  for (int s = 0; s < params->minPart; s++) {
    minIg+=params->p[0][minIg].nonbondedGroupSize;
  }
  int maxI = params->numAtoms[0];
  int strideIg = params->numParts;

  int unique = (gbisParams->numPatches == 1) ? 1 : 0;//should inner loop be unique from ourter loop
  int numGBISPairlists = 4;
  float r_cut = gbisParams->cutoff;
  float fsMax = FS_MAX; // gbisParams->fsMax;//FS_MAX;
  float a_cut = gbisParams->a_cut-fsMax;
  float a_cut2 = a_cut*a_cut;
  float a_cut_ps = a_cut + fsMax;//max screen radius
  float r_cut2 = r_cut*r_cut;
  float a_cut_ps2 = a_cut_ps*a_cut_ps;
  //put PL pointer back to beginning of lists
  for (int k = 0; k < numGBISPairlists; k++)
    gbisParams->gbisStepPairlists[k]->reset();
      

/***********************************************************
* GBIS Phase 1
***********************************************************/
if (gbisParams->gbisPhase == 1) {
//CkPrintf("SEQ%03i, P%i, CID%05i(%02i,%02i):[%03i,%03i]\n",gbisParams->sequence,gbisParams->gbisPhase,gbisParams->cid,gbisParams->patchID[0],gbisParams->patchID[1],minI,maxI);

  pairlistFromAll(params,gbisParams,minIg,strideIg,maxI);

#ifdef BENCHMARK
  int nops = 0;
  int domains[] = {0, 0, 0, 0, 0, 0, 0, 0};
  int numDom = 7;
  double t1 = 1.0*clock()/CLOCKS_PER_SEC;
#endif
  register GBReal psiI;

  register float dr;
  register float r2;
  register float r, r_i, r2_i;
  float rhoi0, rhois, rhojs, rhoj0;
  Position ri, rj;
  register int j;
  int numPairs;
  register float ta = TA;
  register float tb = TB;
  register float tc = TC;
  register float td = TD;
  register float te = TE;

  register float hij,hji;
  int dij,dji;
  float k;
  float rhois2, rhojs2;

  //calculate piecewise-22 Pairs
  int c = 0;
  for (int ngi = minIg; ngi < maxI; /**/ ) {
    int iGroupSize = params->p[0][ngi].nonbondedGroupSize;
  for (int i = ngi; i < ngi+iGroupSize; i++) {
    ri = params->p[0][i].position;
    ri.x += offset_x;
    ri.y += offset_y;
    ri.z += offset_z;
    rhoi0 = gbisParams->intRad[0][2*i+0];
    rhois = gbisParams->intRad[0][2*i+1];
    rhois2 = rhois*rhois;
    psiI = ZERO;
    plint *pairs;
    gbisParams->gbisStepPairlists[c]->nextlist(&pairs,&numPairs);
    for (register int jj = 0; jj < numPairs; jj++) {
#ifdef BENCHMARK
      nops++;
#endif
      j = pairs[jj];
      rj = params->p[1][j].position;

      dr = (ri.x - rj.x);
      r2 = dr*dr;
      dr = (ri.y - rj.y);
      r2 += dr*dr;
      dr = (ri.z - rj.z);
      r2 += dr*dr;
      r2_i = 1.0/r2;

      rhoj0 = gbisParams->intRad[1][2*j+0];
      rhojs = gbisParams->intRad[1][2*j+1];
      rhojs2 = rhojs*rhojs;

      k = rhojs2*r2_i;//k=(rs/r)^2
      hij = rhojs*r2_i*k*(ta+k*(tb+k*(tc+k*(td+k*te))));

      k = rhois2*r2_i;//k=(rs/r)^2
      hji = rhois*r2_i*k*(ta+k*(tb+k*(tc+k*(td+k*te))));

//#if 1
#ifdef PRINT_COMP
      int id1 = params->pExt[0][i].id;
      int id2 = params->pExt[1][j].id;
      float h1 = hij;
      float h2 = hji;
      float r10 = rhoi0;
      float r1s = rhojs;
      float r20 = rhoj0;
      float r2s = rhois;
/*printf("PSIPAIR %05i %05i%9.5f%6.3f%6.3f%2i% 13.5e\n",
id1,id2,sqrt(r2),
rhoi0, rhojs,
2,hij);
printf("PSIPAIR %05i %05i%9.5f%6.3f%6.3f%2i% 13.5e\n",
id2,id1,sqrt(r2),
rhoj0, rhois,
2,hji);
*/
      if (id1 > id2) {
        int tmp = id1;
        id1 = id2;
        id2 = tmp;
        h1 = hji;
        h2 = hij;
        r20 = rhoi0;
        r2s = rhojs;
        r10 = rhoj0;
        r1s = rhois;
      }
//      CkPrintf("PSIPAIR%5i%5i%10.5f%5i%14.8f%7.4f%7.4f\n",id1,id2,sqrt(r2),2,h1,r10,r1s);
//      CkPrintf("PSIPAIR%5i%5i%10.5f%5i%14.8f%7.4f%7.4f\n",id2,id1,sqrt(r2),2,h2,r20,r2s);
      //CkPrintf("PPSI(%04i)[%04i,%04i] = 22\n",gbisParams->cid,id1,id2);
#endif

      psiI += hij;
      gbisParams->psiSum[1][j] += hji;
    }//end inner j
    gbisParams->psiSum[0][i] += psiI;
  }//end outer i
  for (int s = 0; s < strideIg; s++) {
    ngi+=params->p[0][ngi].nonbondedGroupSize;
  }
  }
#ifdef BENCHMARK
  double t2 = 1.0*clock()/CLOCKS_PER_SEC;
  //nops *= (9 + 2*DIV_FLOPS+SQRT_FLOPS) + (14 + 0*DIV_FLOPS + 0*LOG_FLOPS);
  //double flops = 1.0 * nops / (t2 - t1);
  //double gflops = flops / 1000000000.0;
  CkPrintf("PHASE1.1: %8.3f ms @ %8.3f ns/iter for %i iter\n",1000.0*(t2-t1), 1000000000.0*(t2-t1)/nops, nops);
nops = 0;
  t1 = 1.0*clock()/CLOCKS_PER_SEC;
  nops = 0;
#endif

  float rmris, rmrjs;
  float rmrsi;
  float rmrs2;
  float rs2;
  float logri, logrj;
  float rci2;
  float a_cut_i = 1.0 / a_cut;
  float a_cut_i2 = a_cut_i*a_cut_i;

  //calculate piecewise-11 pairs
  c = 1;
  for (int ngi = minIg; ngi < maxI; /**/ ) {
    int iGroupSize = params->p[0][ngi].nonbondedGroupSize;
  for (int i = ngi; i < ngi+iGroupSize; i++) {
    ri = params->p[0][i].position;
    ri.x += offset_x;
    ri.y += offset_y;
    ri.z += offset_z;
    rhoi0 = gbisParams->intRad[0][2*i+0];
    rhois = gbisParams->intRad[0][2*i+1];
    psiI = ZERO;
    plint *pairs;
    gbisParams->gbisStepPairlists[c]->nextlist(&pairs,&numPairs);
    for (register int jj = 0; jj < numPairs; jj++) {
      j = pairs[jj];
      rj = params->p[1][j].position;

      dr = (ri.x - rj.x);
      r2 = dr*dr;
      dr = (ri.y - rj.y);
      r2 += dr*dr;
      dr = (ri.z - rj.z);
      r2 += dr*dr;
      r_i = 1.0/sqrt(r2);
      r = r2*r_i;

      rhoj0 = gbisParams->intRad[1][2*j+0];
      rhojs = gbisParams->intRad[1][2*j+1];

      float tmp1 = 0.125*r_i;
      float tmp2 = r2 - 4.0*a_cut*r;
      float rr = 2.0*r;


      rmrjs = r-rhojs;
      rmris = r-rhois;
      logri = log(rmris*a_cut_i);
      logrj = log(rmrjs*a_cut_i);

      rmrsi = 1.0/rmrjs;
      //rmrs2 = rmrjs*rmrjs;
      rs2 = rhojs*rhojs;
      hij = /*0.125*r_i*/tmp1*(1 + rr*rmrsi +
        a_cut_i2*(/*r2 - 4.0*a_cut*r*/tmp2 - rs2) + logrj+logrj);


      rmrsi = 1.0/rmris;
      //rmrs2 = rmris*rmris;
      rs2 = rhois*rhois;
      hji = /*0.125*r_i*/tmp1*(1 + rr*rmrsi +
        a_cut_i2*(/*r2 - 4.0*a_cut*r*/tmp2 - rs2) + 2.0* logri);
//#if 1
#ifdef PRINT_COMP
      int id1 = params->pExt[0][i].id;
      int id2 = params->pExt[1][j].id;
      float h1 = hij;
      float h2 = hji;
      float r10 = rhoi0;
      float r1s = rhojs;
      float r20 = rhoj0;
      float r2s = rhois;
/*printf("PSIPAIR %05i %05i%9.5f%6.3f%6.3f%2i% 13.5e\n",
id1,id2,sqrt(r2),
rhoi0, rhojs,
1,hij);
printf("PSIPAIR %05i %05i%9.5f%6.3f%6.3f%2i% 13.5e\n",
id2,id1,sqrt(r2),
rhoj0, rhois,
1,hji);*/
      if (id1 > id2) {
        int tmp = id1;
        id1 = id2;
        id2 = tmp;
        h1 = hji;
        h2 = hij;
        r20 = rhoi0;
        r2s = rhojs;
        r10 = rhoj0;
        r1s = rhois;
      }
//      CkPrintf("PSIPAIR%5i%5i%10.5f%5i%14.8f%7.4f%7.4f\n",id1,id2,sqrt(r2),1,h1,r10,r1s);
//      CkPrintf("PSIPAIR%5i%5i%10.5f%5i%14.8f%7.4f%7.4f\n",id2,id1,sqrt(r2),1,h2,r20,r2s);
      //CkPrintf("PSI(%04i)[%04i,%04i] = 11 % .4e % .4e\n",gbisParams->sequence,id1,id2,h1, h2);
      //CkPrintf("PPSI(%04i)[%04i,%04i] = 11\n",gbisParams->cid,id1,id2);
#endif

#ifdef BENCHMARK
      nops++;
#endif
          
      psiI += hij;
      gbisParams->psiSum[1][j] += hji;
    }//end inner j
    gbisParams->psiSum[0][i] += psiI;
  }//end outer i
  for (int s = 0; s < strideIg; s++) {
    ngi+=params->p[0][ngi].nonbondedGroupSize;
  }
  }
#ifdef BENCHMARK
  t2 = 1.0*clock()/CLOCKS_PER_SEC;
  CkPrintf("PHASE1.2: %8.3f ms @ %8.3f ns/iter for %i iter\n",1000.0*(t2-t1), 1000000000.0*(t2-t1)/nops, nops);
nops = 0;
  t1 = 1.0*clock()/CLOCKS_PER_SEC;
  nops = 0;
#endif

  //calculate all other piecewise pairs
  c = 2;
  for (int ngi = minIg; ngi < maxI; /**/ ) {
    int iGroupSize = params->p[0][ngi].nonbondedGroupSize;
  for (int i = ngi; i < ngi+iGroupSize; i++) {
    ri = params->p[0][i].position;
    ri.x += offset_x;
    ri.y += offset_y;
    ri.z += offset_z;
    rhoi0 = gbisParams->intRad[0][2*i+0];
    rhois = gbisParams->intRad[0][2*i+1];
    psiI = ZERO;
    plint *pairs;
    gbisParams->gbisStepPairlists[c]->nextlist(&pairs,&numPairs);
    for (register int jj = 0; jj < numPairs; jj++) {
      j = pairs[jj];
      rj = params->p[1][j].position;

      dr = (ri.x - rj.x);
      r2 = dr*dr;
      dr = (ri.y - rj.y);
      r2 += dr*dr;
      dr = (ri.z - rj.z);
      r2 += dr*dr;
      rhojs = gbisParams->intRad[1][2*j+1];
      rhoj0 = gbisParams->intRad[1][2*j+0];
      r_i = 1.0/sqrt(r2);
      r = r2 * r_i;

      CalcHPair(r,r2,r_i,a_cut,
          rhoi0, rhojs,
          rhoj0, rhois,
          dij,dji,hij,hji);
//#if 1
#ifdef PRINT_COMP
      int id1 = params->pExt[0][i].id;
      int id2 = params->pExt[1][j].id;
      float h1 = hij;
      float h2 = hji;
      int d1 = dij;
      int d2 = dji;
      float r10 = rhoi0;
      float r1s = rhojs;
      float r20 = rhoj0;
      float r2s = rhois;
/*  if (dij > 0 ) {
printf("PSIPAIR %05i %05i%9.5f%6.3f%6.3f%2i% 13.5e\n",
id1,id2,sqrt(r2),
rhoi0, rhojs,
dij,hij);
  }
  if (dji > 0 ) {
printf("PSIPAIR %05i %05i%9.5f%6.3f%6.3f%2i% 13.5e\n",
id2,id1,sqrt(r2),
rhoj0, rhois,
dji,hji);
  }*/
//      CkPrintf("PSIPAIR%5i%5i%10.5f%5i%14.8f%7.4f%7.4f\n",id1,id2,sqrt(r2),d1,h1,r10,r1s);
      if (id1 > id2) {
        int tmp = id1;
        id1 = id2;
        id2 = tmp;
        h1 = hji;
        h2 = hij;
        d1 = dji;
        d2 = dij;
        r20 = rhoi0;
        r2s = rhojs;
        r10 = rhoj0;
        r1s = rhois;
      }
//      CkPrintf("PSIPAIR%5i%5i%10.5f%5i%14.8f%7.4f%7.4f\n",id2,id1,sqrt(r2),d2,h2,r20,r2s);
      //CkPrintf("PSI(%04i)[%04i,%04i] = %i%i % .4e % .4e\n",gbisParams->sequence,id1,id2,d1,d2,h1, h2);
      //CkPrintf("PPSI(%04i)[%04i,%04i] = %i%i\n",gbisParams->cid,id1,id2,d1,d2);
#endif

#ifdef BENCHMARK
      nops++;
#endif
      psiI += hij;
      gbisParams->psiSum[1][j] += hji;
    }//end inner j
    gbisParams->psiSum[0][i] += psiI;
  }//end outer i
  for (int s = 0; s < strideIg; s++) {
    ngi+=params->p[0][ngi].nonbondedGroupSize;
  }
  }

#ifdef BENCHMARK
  t2 = 1.0*clock()/CLOCKS_PER_SEC;
  CkPrintf("PHASE1.3: %8.3f ms @ %8.3f ns/iter for %i iter\n",1000.0*(t2-t1), 1000000000.0*(t2-t1)/nops,nops);
#endif

/***********************************************************
* GBIS Phase 2
***********************************************************/
} else if (gbisParams->gbisPhase == 2) {

  float epsilon_s = gbisParams->epsilon_s;
  float epsilon_p = gbisParams->epsilon_p;
  float epsilon_s_i = 1/epsilon_s;
  float epsilon_p_i = 1/epsilon_p;
  float kappa = gbisParams->kappa;

  //values used in loop
  float r_cut_2 = 1.0 / r_cut2;
  float r_cut_4 = 4.0*r_cut_2*r_cut_2;
  float coulEij=0,ddrCoulEij=0,gbEij=0,ddrGbEij=0;
  float dEdai=0,dEdaj=0, qiqj=0;
  float scale=0, ddrScale=0;
  float rnx=0,rny=0,rnz=0;
  float fx=0,fy=0,fz=0,forceCoul=0, forcedEdr=0;

  int nops = 0;
  double t1 = 1.0*clock()/CLOCKS_PER_SEC;
  float r2;
  float dr;
  BigReal dx, dy, dz;
  float r, r_i, ratio;
  float fIx, fIy, fIz;
  GBReal dEdaSumI;
  float bornRadI, bornRadJ;
  float qi;
  Position ri, rj;
  float aiaj,expr2aiaj4,fij,f_i,expkappa,Dij;
  float aiaj4,ddrDij,ddrf_i,ddrfij,tmp_dEda;

  for (int c = 0; c < 4/*dEdrPLs*/; c++) {
  for (int ngi = minIg; ngi < maxI; /**/ ) {
    int iGroupSize = params->p[0][ngi].nonbondedGroupSize;
  for (int i = ngi; i < ngi+iGroupSize; i++) {
    ri = params->p[0][i].position;
    ri.x += offset_x;
    ri.y += offset_y;
    ri.z += offset_z;
    qi = - COULOMB * params->p[0][i].charge * scaling;
  //printf("ATOM(%05i) %.3e %.3e %.3e\n", params->pExt[0][i].id, params->p[0][i].charge, COULOMB, scaling);
    int numPairs;
    plint *pairs;
    gbisParams->gbisStepPairlists[c]->nextlist(&pairs,&numPairs);
    fIx = fIy = fIz = 0.0;
    dEdaSumI = 0.0;
    bornRadI = gbisParams->bornRad[0][i];
    for (int jj = 0; jj < numPairs; jj++) {
      int j = pairs[jj];
      rj = params->p[1][j].position;

      //distance
      dx = (ri.x - rj.x);
      dy = (ri.y - rj.y);
      dz = (ri.z - rj.z);
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 > r_cut2) continue;
      qiqj = qi*params->p[1][j].charge;
      bornRadJ = gbisParams->bornRad[1][j];
      r_i = 1.0/sqrt(r2);
      r = r2 * r_i;

      //pairwise calculation
/*      Phase2_Pair(r, r2, r_i, qiqj,
          bornRadI,
          bornRadJ,
          epsilon_p_i, epsilon_s_i,
          kappa, gbisParams->doFullElectrostatics,
          gbEij, ddrGbEij, dEdai, dEdaj);
*/

  //calculate GB energy
  aiaj = bornRadI*bornRadJ;
  aiaj4 = 4*aiaj;
  expr2aiaj4 = exp(-r2/aiaj4);
  fij = sqrt(r2+aiaj*expr2aiaj4);
  f_i = 1/fij;
  expkappa = (kappa > 0) ? exp(-kappa*fij) : 1.0;
  Dij = epsilon_p_i - expkappa*epsilon_s_i;
  gbEij = qiqj*Dij*f_i;

  //calculate energy derivatives
  ddrfij = r*f_i*(1 - 0.25*expr2aiaj4);
  ddrf_i = -ddrfij*f_i*f_i;
  ddrDij = kappa*expkappa*ddrfij*epsilon_s_i;
  ddrGbEij = qiqj*(ddrDij*f_i+Dij*ddrf_i);

  //calc dEda
      //NAMD smoothing function
      scale = 1;
      ddrScale = 0;
      if (gbisParams->doSmoothing) {
        scale = r2 * r_cut_2 - 1;
        scale *= scale;
        ddrScale = r*(r2-r_cut2)*r_cut_4;
        //CkPrintf("SCALE %f %f\n",scale,ddrScale);
        gbisParams->gbInterEnergy += gbEij * scale;
        forcedEdr = -(ddrGbEij)*scale-(gbEij)*ddrScale;
      } else {
        gbisParams->gbInterEnergy += gbEij;
        forcedEdr = -ddrGbEij;
      }

      //add dEda
      if (gbisParams->doFullElectrostatics) {
        //gbisParams->dEdaSum[0][i] += dEdai*scale;
        tmp_dEda = 0.5*qiqj*f_i*f_i
                      *(kappa*epsilon_s_i*expkappa-Dij*f_i)
                      *(aiaj+0.25*r2)*expr2aiaj4;//0
        dEdai = tmp_dEda/bornRadI;
        dEdaj = tmp_dEda/bornRadJ;
        dEdaSumI += dEdai*scale;
        gbisParams->dEdaSum[1][j] += dEdaj*scale;
      }

#if 1
//#ifdef PRINT_COMP
      int id1 = params->pExt[0][i].id;
      int id2 = params->pExt[1][j].id;
      float deda1 = dEdai;
      float deda2 = dEdaj;
      float bR1 = bornRadI;
      float bR2 = bornRadJ;
      if (id1 > id2) {
        int tmp = id1;
        id1 = id2;
        id2 = tmp;
        deda1 = dEdaj;
        deda2 = dEdai;
        bR1 = bornRadJ;
        bR2 = bornRadI;
      }
      //CkPrintf("DEDR(%04i)[%04i,%04i] = % .4e\n",gbisParams->sequence,id1,id2,forcedEdr);
      //CkPrintf("DASM(%04i)[%04i,%04i] = % .4e % .4e\n",gbisParams->sequence,id1, id2, deda1,deda2);
      //CkPrintf("P2RM(%04i)[%04i,%04i] = % .4e % .4e % .4e % .4e % .4e\n",gbisParams->sequence,id1, id2, r, bR1,bR2,epsilon_p_i,epsilon_s_i,kappa);
/*CkPrintf("P2PAIR %05i %05i%9.5f%6.3f%6.3f% 13.5e% 13.5e% 13.5e% 13.5e% 13.5e\n",
params->pExt[0][i].id, params->pExt[1][j].id,sqrt(r2),
bornRadI, bornRadJ, forcedEdr, dEdai, qiqj, Dij, scale
);
CkPrintf("P2PAIR %05i %05i%9.5f%6.3f%6.3f% 13.5e% 13.5e% 13.5e% 13.5e% 13.5e\n",
params->pExt[1][j].id, params->pExt[0][i].id,sqrt(r2),
bornRadJ, bornRadI, forcedEdr, dEdaj, qiqj, Dij, scale
);*/
#endif

      forcedEdr *= r_i;
      fx = dx*forcedEdr;
      fy = dy*forcedEdr;
      fz = dz*forcedEdr;

      params->ff[1][j].x -= fx;
      params->ff[1][j].y -= fy;
      params->ff[1][j].z -= fz;

      fIx += fx;
      fIy += fy;
      fIz += fz;

#ifdef BENCHMARK
      nops ++;//= (59 + 4*DIV_FLOPS + 2*EXP_FLOPS+SQRT_FLOPS);//56 w/o nops
#endif

    }//end inner j
    gbisParams->dEdaSum[0][i] += dEdaSumI;
    params->ff[0][i].x += fIx;
    params->ff[0][i].y += fIy;
    params->ff[0][i].z += fIz;

    //self energy of each atom
    if (c == 0 && gbisParams->doEnergy && gbisParams->numPatches == 1) {
      float fij = bornRadI;//inf
      float expkappa = exp(-kappa*fij);//0
      float Dij = epsilon_p_i - expkappa*epsilon_s_i;
      float gbEij = qi*params->p[0][i].charge*Dij/fij;
      gbisParams->gbSelfEnergy += 0.5*gbEij;//self energy
    }
  }// end outer i
  for (int s = 0; s < strideIg; s++) {
    ngi+=params->p[0][ngi].nonbondedGroupSize;
  }//end i
  }//end ig
  }//end c
#ifdef BENCHMARK
  double t2 = 1.0*clock()/CLOCKS_PER_SEC;
  //double flops = 1.0 * nops / (t2 - t1);
  //double gflops = flops / 1000000000.0;
  
  CkPrintf("PHASE2.0: %8.3f ms @ %8.3f ns/iter for %i iter\n",1000.0*(t2-t1), 1000000000.0*(t2-t1)/nops,nops);
#endif

/***********************************************************
* GBIS Phase 3
***********************************************************/
} else if (gbisParams->gbisPhase == 3 && gbisParams->doFullElectrostatics) {

#ifdef BENCHMARK
  //int domains[] = {0, 0, 0, 0, 0, 0, 0, 0};
  double t1 = 1.0*clock()/CLOCKS_PER_SEC;
  double t2;
  int nops = 0;
#endif
  //CkPrintf("GBIS(%3i)[%2i]::P3 %3i(%3i) %3i(%3i)\n",gbisParams->sequence,gbisParams->cid, gbisParams->patchID[0],params->numAtoms[0],gbisParams->patchID[1],params->numAtoms[1]);

  register BigReal dx, dy, dz;
  register float  r2;
  register float r, r_i;
  register float rhoi0, rhois;
  float rhojs, rhoj0;
  float fx, fy, fz;
  float forceAlpha;
  register float fIx, fIy, fIz;

  float dhij;
  float dhji;
  int dij;
  int dji;
  register float dHdrPrefixI;
  float dHdrPrefixJ;
  register Position ri;
  register Position rj;
  register int c, numPairs, jj, j;
  register float k;
  register float da = DA;
  register float db = DB;
  register float dc = DC;
  register float dd = DD;
  register float de = DE;
  float r_i3;

  //piecewise 22
  c = 0;
  for (int ngi = minIg; ngi < maxI; /**/ ) {
    int iGroupSize = params->p[0][ngi].nonbondedGroupSize;
  for (int i = ngi; i < ngi+iGroupSize; i++) {
    ri = params->p[0][i].position;
    ri.x += offset_x;
    ri.y += offset_y;
    ri.z += offset_z;
    rhois = gbisParams->intRad[0][2*i+1];
    plint *pairs;
    gbisParams->gbisStepPairlists[c]->nextlist(&pairs,&numPairs);
    fIx = fIy = fIz = 0.0;
    dHdrPrefixI = gbisParams->dHdrPrefix[0][i];
    for (jj = 0; jj < numPairs; jj++) {
      j = pairs[jj];
      rj = params->p[1][j].position;

      dx = (ri.x - rj.x);
      dy = (ri.y - rj.y);
      dz = (ri.z - rj.z);
      r2 = dx*dx + dy*dy + dz*dz;
      dHdrPrefixJ = gbisParams->dHdrPrefix[1][j];

      r_i = 1.0/sqrt(r2);//rptI takes 50% of loop time
      r_i3 = r_i*r_i*r_i;

      rhojs = gbisParams->intRad[1][2*j+1];

      k = rhojs*r_i; k*=k;//k=(rs/r)^2
      dhij = -rhojs*r_i3*k*
              (da+k*(db+k*(dc+k*(dd+k*de))));

      k = rhois*r_i; k*=k;//k=(rs/r)^2
      dhji = -rhois*r_i3*k*
              (da+k*(db+k*(dc+k*(dd+k*de))));

      //add dEda*dadr force
      forceAlpha = -r_i*(dHdrPrefixI*dhij+dHdrPrefixJ*dhji);
      fx = dx * forceAlpha;
      fy = dy * forceAlpha;
      fz = dz * forceAlpha;

      params->fullf[1][j].x -= fx;
      params->fullf[1][j].y -= fy;
      params->fullf[1][j].z -= fz;

      fIx += fx;
      fIy += fy;
      fIz += fz;

#if 1
//#ifdef PRINT_COMP
      int id1 = params->pExt[0][i].id;
      int id2 = params->pExt[1][j].id;
      float h1 = dhij;
      float h2 = dhji;
      if (id1 > id2) {
        int tmp = id1;
        id1 = id2;
        id2 = tmp;
        h1 = dhji;
        h2 = dhij;
      }
/*CkPrintf("P3PAIR %05i %05i%9.5f% 13.5e% 13.5e% 13.5e% 13.5e% 13.5e %i\n",
params->pExt[0][i].id, params->pExt[1][j].id,sqrt(r2),
dHdrPrefixI, dHdrPrefixJ, dhij, dhji, forceAlpha, 2
);
CkPrintf("P3PAIR %05i %05i%9.5f% 13.5e% 13.5e% 13.5e% 13.5e% 13.5e %i\n",
params->pExt[1][j].id, params->pExt[0][i].id,sqrt(r2),
dHdrPrefixJ, dHdrPrefixI, dhji, dhij, forceAlpha, 2
);*/
      //CkPrintf("DEDA(%04i)[%04i,%04i] = 22 % .4e % .4e % .4e\n",gbisParams->sequence,id1,id2,h1,h2, forceAlpha);
#endif

#ifdef BENCHMARK
      nops++;//= (24 + 2*DIV_FLOPS+SQRT_FLOPS);// 8 w/o nops
#endif

    }//end inner j
    params->fullf[0][i].x += fIx;
    params->fullf[0][i].y += fIy;
    params->fullf[0][i].z += fIz;
  }//end outer i
  for (int s = 0; s < strideIg; s++) {
    ngi+=params->p[0][ngi].nonbondedGroupSize;
  }
  }

#ifdef BENCHMARK
  t2 = 1.0*clock()/CLOCKS_PER_SEC;
  CkPrintf("PHASE3.1: %8.3f ms @ %8.3f ns/iter for %i iter\n",1000.0*(t2-t1), 1000000000.0*(t2-t1)/nops,nops);
nops = 0;
  t1 = 1.0*clock()/CLOCKS_PER_SEC;
  nops = 0;
#endif

  float a_cut_i = 1.0/a_cut;
  float a_cut_i2 = a_cut_i*a_cut_i;
  float a_cut2 = a_cut*a_cut;
  float rmrs;
  float rmrsi;
  float rmrs2;
  float rhois2, rhojs2;
  float logri, logrj;
  float r_i2;

  //piecewise 11
  c = 1;
  for (int ngi = minIg; ngi < maxI; /**/ ) {
    int iGroupSize = params->p[0][ngi].nonbondedGroupSize;
  for (int i = ngi; i < ngi+iGroupSize; i++) {
    ri = params->p[0][i].position;
    ri.x += offset_x;
    ri.y += offset_y;
    ri.z += offset_z;
    rhois = gbisParams->intRad[0][2*i+1];
    rhois2 = rhois*rhois;
    plint *pairs;
    gbisParams->gbisStepPairlists[c]->nextlist(&pairs,&numPairs);
    fIx = fIy = fIz = 0.0;
    dHdrPrefixI = gbisParams->dHdrPrefix[0][i];
    for (jj = 0; jj < numPairs; jj++) {
      j = pairs[jj];
      rj = params->p[1][j].position;
      dHdrPrefixJ = gbisParams->dHdrPrefix[1][j];

      dx = (ri.x - rj.x);
      dy = (ri.y - rj.y);
      dz = (ri.z - rj.z);
      r2 = dx*dx + dy*dy + dz*dz;
      r_i = 1.0/sqrt(r2);//rptI
      r = r2* r_i;
      r_i2 = r_i*r_i;

      rhojs = gbisParams->intRad[1][2*j+1];
      rhojs2 = rhojs*rhojs;


      rmrs = r-rhojs;// 4 times
      rmrsi = 1.0/rmrs;
      rmrs2 = rmrs*rmrs;
      logrj = log(rmrs*a_cut_i);
      dhij = r_i2*(-0.25*logrj - (a_cut2 - rmrs2)*(rhojs2 + r2)*0.125*a_cut_i2*rmrsi*rmrsi);


      rmrs = r-rhois;// 4 times
      rmrsi = 1.0/rmrs;
      rmrs2 = rmrs*rmrs;
      logri = log(rmrs*a_cut_i);
      dhji = r_i2*(-0.25*logri - (a_cut2 - rmrs2)*(rhois2 + r2)*0.125*a_cut_i2*rmrsi*rmrsi);

      //add dEda*dadr force
      forceAlpha = -r_i*(dHdrPrefixI*dhij+dHdrPrefixJ*dhji);
      fx = dx * forceAlpha;
      fy = dy * forceAlpha;
      fz = dz * forceAlpha;

      params->fullf[1][j].x -= fx;
      params->fullf[1][j].y -= fy;
      params->fullf[1][j].z -= fz;

      fIx += fx;
      fIy += fy;
      fIz += fz;

#if 1
//#ifdef PRINT_COMP
      int id1 = params->pExt[0][i].id;
      int id2 = params->pExt[1][j].id;
      float h1 = dhij;
      float h2 = dhji;
      if (id1 > id2) {
        int tmp = id1;
        id1 = id2;
        id2 = tmp;
        h1 = dhji;
        h2 = dhij;
      }
      //CkPrintf("DEDA(%04i)[%04i,%04i] = 11 % .4e % .4e % .4e\n",gbisParams->sequence,id1,id2,h1,h2, forceAlpha);
/*CkPrintf("P3PAIR %05i %05i%9.5f% 13.5e% 13.5e% 13.5e% 13.5e% 13.5e %i\n",
params->pExt[0][i].id, params->pExt[1][j].id,sqrt(r2),
dHdrPrefixI, dHdrPrefixJ, dhij, dhji, forceAlpha, 1
);
CkPrintf("P3PAIR %05i %05i%9.5f% 13.5e% 13.5e% 13.5e% 13.5e% 13.5e %i\n",
params->pExt[1][j].id, params->pExt[0][i].id,sqrt(r2),
dHdrPrefixJ, dHdrPrefixI, dhji, dhij, forceAlpha, 1
);*/
#endif

#ifdef BENCHMARK
      nops++;
#endif

    }//end inner j
    params->fullf[0][i].x += fIx;
    params->fullf[0][i].y += fIy;
    params->fullf[0][i].z += fIz;
  }//end outer i
  for (int s = 0; s < strideIg; s++) {
    ngi+=params->p[0][ngi].nonbondedGroupSize;
  }
  }

#ifdef BENCHMARK
  t2 = 1.0*clock()/CLOCKS_PER_SEC;
  CkPrintf("PHASE3.2: %8.3f ms @ %8.3f ns/iter for %i iter\n",1000.0*(t2-t1), 1000000000.0*(t2-t1)/nops,nops);
nops = 0;
  t1 = 1.0*clock()/CLOCKS_PER_SEC;
  nops = 0;
#endif

  //piecewise all others
  c = 2;
  for (int ngi = minIg; ngi < maxI; /**/ ) {
    int iGroupSize = params->p[0][ngi].nonbondedGroupSize;
  for (int i = ngi; i < ngi+iGroupSize; i++) {
    ri = params->p[0][i].position;
    ri.x += offset_x;
    ri.y += offset_y;
    ri.z += offset_z;
    rhoi0 = gbisParams->intRad[0][2*i+0];
    rhois = gbisParams->intRad[0][2*i+1];
    plint *pairs;
    gbisParams->gbisStepPairlists[c]->nextlist(&pairs,&numPairs);
    fIx = fIy = fIz = 0.0;
    dHdrPrefixI = gbisParams->dHdrPrefix[0][i];
    for (jj = 0; jj < numPairs; jj++) {
      j = pairs[jj];
      rj = params->p[1][j].position;
      dHdrPrefixJ = gbisParams->dHdrPrefix[1][j];

      dx = (ri.x - rj.x);
      dy = (ri.y - rj.y);
      dz = (ri.z - rj.z);
      r2 = dx*dx + dy*dy + dz*dz;

      r_i = 1.0/sqrt(r2);//rptI
      r = r2* r_i;

      rhojs = gbisParams->intRad[1][2*j+1];
      rhoj0 = gbisParams->intRad[1][2*j+0];

      CalcDHPair(r,r2,r_i,a_cut,
          rhoi0,rhojs,
          rhoj0,rhois,
          dij,dji,dhij,dhji);

      //add dEda*dadr force
      forceAlpha = -r_i*(dHdrPrefixI*dhij+dHdrPrefixJ*dhji); // *scaling ?
      fx = dx * forceAlpha;
      fy = dy * forceAlpha;
      fz = dz * forceAlpha;

      fIx += fx;
      fIy += fy;
      fIz += fz;

      params->fullf[1][j].x -= fx;
      params->fullf[1][j].y -= fy;
      params->fullf[1][j].z -= fz;

#if 1
//#ifdef PRINT_COMP
      int id1 = params->pExt[0][i].id;
      int id2 = params->pExt[1][j].id;
      int d1 = dij;
      int d2 = dji;
      float h1 = dhij;
      float h2 = dhji;
      if (id1 > id2) {
        int tmp = id1;
        id1 = id2;
        id2 = tmp;
        d1 = dji;
        d2 = dij;
        h1 = dhji;
        h2 = dhij;
      }
//      CkPrintf("DEDA(%04i)[%04i,%04i] = %i%i % .4e % .4e % .4e\n",gbisParams->sequence,id1,id2,d1,d2,h1,h2, forceAlpha);
//      CkPrintf("DEDA(%04i)[%04i,%04i] = %i%i % .4e % .4e % .4e\n",gbisParams->sequence,id1,id2,d1,d2,h1,h2, forceAlpha);
/*if ( dij > 0 ) {
CkPrintf("P3PAIR %05i %05i%9.5f% 13.5e% 13.5e% 13.5e% 13.5e% 13.5e %i\n",
params->pExt[0][i].id, params->pExt[1][j].id,sqrt(r2),
dHdrPrefixI, dHdrPrefixJ, dhij, dhji, forceAlpha, dij
);
}
if ( dji > 0 ) {
CkPrintf("P3PAIR %05i %05i%9.5f% 13.5e% 13.5e% 13.5e% 13.5e% 13.5e %i\n",
params->pExt[1][j].id, params->pExt[0][i].id,sqrt(r2),
dHdrPrefixJ, dHdrPrefixI, dhji, dhij, forceAlpha, dji
);
}*/
#endif

#ifdef BENCHMARK
      nops++;
#endif

    }//end inner j
    params->fullf[0][i].x += fIx;
    params->fullf[0][i].y += fIy;
    params->fullf[0][i].z += fIz;
  }//end outer i
  for (int s = 0; s < strideIg; s++) {
    ngi+=params->p[0][ngi].nonbondedGroupSize;
  }
  }


#ifdef BENCHMARK
  t2 = 1.0*clock()/CLOCKS_PER_SEC;
  CkPrintf("PHASE3.3: %8.3f ms @ %8.3f ns/iter for %i iter\n",1000.0*(t2-t1), 1000000000.0*(t2-t1)/nops, nops);
#endif

}//end if gbisPhase

}//end calcGBIS
