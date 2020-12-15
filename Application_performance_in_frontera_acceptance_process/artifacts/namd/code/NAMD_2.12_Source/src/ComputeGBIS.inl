/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*******************************************************************************
  These are inline functions which define all the atom-pair calculations
  included in ComputeGBIS.C and ComputeGBISCUDAKernal.h
*******************************************************************************/

#ifndef COMPUTEGBIS_INL
#define COMPUTEGBIS_INL

/* Should psiSum and dEdaSum reductions be performed with doubles? */
//this value is now defined in NamdTypes
typedef float GBReal;

//copied from NamdTypes.h
typedef float Mass;


#ifndef FS_MAX
#define FS_MAX 1.728f /*maximum screened radius*/
#endif

#ifndef COULOMB
#define COULOMB 332.0636f /* ke [kcal*Ang/e^2] common.h */
#endif

#ifndef TA
#define TA 0.333333333333333f // 1/3
#define TB 0.4f               // 2/5
#define TC 0.428571428571428f // 3/7
#define TD 0.444444444444444f // 4/9
#define TE 0.454545454545454f // 5/11
#define DA 1.333333333333333f // 4* 1/3
#define DB 2.4f               // 6* 2/5
#define DC 3.428571428571428f // 8* 3/7
#define DD 4.444444444444444f // 10*4/9
#define DE 5.454545454545454f // 12*5/11
#endif

static inline float FastTanH( float x ) {
  float a = 2.f*x+0.02f;;
  a *= (6.f + a*(3.f + a));
  return (a)/(a+12.f);
}

/******************************************************************************
 * use mass to determine element
 * return Bondi radius
 * returns the radius of i, depends on j though
 ******************************************************************************/
static inline float MassToRadius(Mass mi) {//, Mass mj) {
return
  (mi <   2.50f ) ? 1.20f : // AtmNum = 1;  Elem =  H ; Mass =   1.00
  (mi <   5.47f ) ? 1.40f : // AtmNum = 2;  Elem =  He; Mass =   4.00
  (mi <   7.98f ) ? 1.82f : // AtmNum = 3;  Elem =  Li; Mass =   6.94
  (mi <   9.91f ) ? 2.13f : // AtmNum = 4;  Elem =  Be; Mass =   9.01
  (mi <  11.41f ) ? 2.13f : // AtmNum = 5;  Elem =  B ; Mass =  10.81
  (mi <  13.01f ) ? 1.70f : // AtmNum = 6;  Elem =  C ; Mass =  12.01
  (mi <  15.00f ) ? 1.55f : // AtmNum = 7;  Elem =  N ; Mass =  14.00
  (mi <  17.49f ) ? 1.50f : // AtmNum = 8;  Elem =  O ; Mass =  15.99
  (mi <  19.58f ) ? 1.50f : // AtmNum = 9;  Elem =  F ; Mass =  18.99
  (mi <  21.58f ) ? 1.54f : // AtmNum = 10; Elem =  Ne; Mass =  20.17
  (mi <  23.64f ) ? 2.27f : // AtmNum = 11; Elem =  Na; Mass =  22.98
  (mi <  25.64f ) ? 1.73f : // AtmNum = 12; Elem =  Mg; Mass =  24.30
  (mi <  27.53f ) ? 2.51f : // AtmNum = 13; Elem =  Al; Mass =  26.98
  (mi <  29.53f ) ? 2.10f : // AtmNum = 14; Elem =  Si; Mass =  28.08
  (mi <  31.52f ) ? 1.85f : // AtmNum = 15; Elem =  P ; Mass =  30.97
  (mi <  33.76f ) ? 1.80f : // AtmNum = 16; Elem =  S ; Mass =  32.06
  (mi <  37.28f ) ? 1.70f : // AtmNum = 17; Elem =  Cl; Mass =  35.45
  (mi <  39.29f ) ? 2.75f : // AtmNum = 19; Elem =  K ; Mass =  39.10
  (mi <  49.09f ) ? 1.88f : // AtmNum = 18; Elem =  Ar; Mass =  39.48
  (mi <  61.12f ) ? 1.63f : // AtmNum = 28; Elem =  Ni; Mass =  58.69
  (mi <  64.46f ) ? 1.40f : // AtmNum = 29; Elem =  Cu; Mass =  63.54
  (mi <  67.55f ) ? 1.39f : // AtmNum = 30; Elem =  Zn; Mass =  65.38
  (mi <  71.18f ) ? 1.87f : // AtmNum = 31; Elem =  Ga; Mass =  69.72
  (mi <  73.78f ) ? 2.19f : // AtmNum = 32; Elem =  Ge; Mass =  72.64
  (mi <  76.94f ) ? 1.85f : // AtmNum = 33; Elem =  As; Mass =  74.92
  (mi <  79.43f ) ? 1.90f : // AtmNum = 34; Elem =  Se; Mass =  78.96
  (mi <  81.85f ) ? 1.85f : // AtmNum = 35; Elem =  Br; Mass =  79.90
  (mi <  95.11f ) ? 2.02f : // AtmNum = 36; Elem =  Kr; Mass =  83.79
  (mi < 107.14f ) ? 1.63f : // AtmNum = 46; Elem =  Pd; Mass = 106.42
  (mi < 110.14f ) ? 1.72f : // AtmNum = 47; Elem =  Ag; Mass = 107.86
  (mi < 113.61f ) ? 1.58f : // AtmNum = 48; Elem =  Cd; Mass = 112.41
  (mi < 116.76f ) ? 1.93f : // AtmNum = 49; Elem =  In; Mass = 114.81
  (mi < 120.24f ) ? 2.17f : // AtmNum = 50; Elem =  Sn; Mass = 118.71
  (mi < 124.33f ) ? 2.09f : // AtmNum = 51; Elem =  Sb; Mass = 121.76
  (mi < 127.25f ) ? 1.98f : // AtmNum = 53; Elem =  I ; Mass = 126.90
  (mi < 129.45f ) ? 2.06f : // AtmNum = 52; Elem =  Te; Mass = 127.60
  (mi < 163.19f ) ? 2.16f : // AtmNum = 54; Elem =  Xe; Mass = 131.29
  (mi < 196.02f ) ? 1.75f : // AtmNum = 78; Elem =  Pt; Mass = 195.08
  (mi < 198.78f ) ? 1.66f : // AtmNum = 79; Elem =  Au; Mass = 196.96
  (mi < 202.49f ) ? 1.55f : // AtmNum = 80; Elem =  Hg; Mass = 200.59
  (mi < 205.79f ) ? 1.96f : // AtmNum = 81; Elem =  Tl; Mass = 204.38
  (mi < 222.61f ) ? 2.02f : // AtmNum = 82; Elem =  Pb; Mass = 207.20
  (mi < 119.01f ) ? 1.86f : // AtmNum = 92; Elem =  U ; Mass = 238.02
                   1.50f ;   // Unknown
}

/******************************************************************************
 * Screen radii
 * use masses to determine elements
 * use elements to lookup Sij
 * to scale the coulomb radius
 * from Hawkins, Cramer, Truhlar; 1996
 * mi is descreened atom - calculating it's alpha (outer loop index)
 * mj is descreening atom - contributor (inner loop index)
 ******************************************************************************/
static inline float MassToScreen(Mass mi) {//, Mass mj) {
    return
      (mi <   1.500f) ? 0.85f : //H
      (mi <  12.500f) ? 0.72f : //C
      (mi <  14.500f) ? 0.79f : //N
      (mi <  16.500f) ? 0.85f : //O
      (mi <  19.500f) ? 0.88f : //F
      (mi <  31.500f) ? 0.86f : //P
      (mi <  32.500f) ? 0.96f : //S
                        0.8f ; //all others
}


/******************************************************************************
  Piecewise screening functions Hij dHij/drij
  r   distance
  r2  square distance
  ri  inverse distance
  rc  cutoff
  r0  descreened atom radius
  rs  descreening atom radius
  h   return value
  dh  return value
 ******************************************************************************/
#ifdef GBIS_CUDA
__device__ __forceinline__ void h0
#else
static inline void h0
#endif
( float r, float r2, float ri,//(0)*5.3%
    float rc, float r0, float rs, float & h ) {
  h = 0.f;
}
#ifdef GBIS_CUDA
__device__ __forceinline__ void dh0
#else
static inline void dh0
#endif
( float r, float r2, float ri,//(0)*5.3%
    float rc, float r0, float rs, float & dh ) {
  dh = 0.f;
}

#ifdef GBIS_CUDA
__device__ __forceinline__ void h1
#else
static inline void h1
#endif
( float r, float r2, float ri, //(6+ 11* 2/ 1log)*18.4%
    float rc, float r0, float rs, float & h ) {

  float rci = 1.f/rc;
  float rmrs = r-rs;
  float rmrsi = 1.f/rmrs;
  //float rmrs2 = rmrs*rmrs;
  float rs2 = rs*rs;
  float logr = log(rmrs*rci);
  float rci2 = rci*rci;
  h = 0.125f*ri*(1.f + 2.f*r*rmrsi + rci2*(r2 - 4.f*rc*r - rs2) + 2.f*logr);
}
#ifdef GBIS_CUDA
__device__ __forceinline__ void dh1
#else
static inline void dh1
#endif
( float r, float r2, float ri, //(4+ 13* 2/ 1log)*18.4%
    float rc, float r0, float rs, float & dh ) {

  float rci = 1.f/rc;
  float rmrs = r-rs;// 4 times
  float rmrsi = 1.f/rmrs;
  float rmrs2 = rmrs*rmrs;
  float rs2 = rs*rs;
  float logr = log(rmrs*rci);
  float rci2 = rci*rci;
  dh = ri*ri*(-0.25f*logr - (rc*rc - rmrs2)*(rs2 + r2)*0.125f*rci2*rmrsi*rmrsi);
}

#ifdef GBIS_CUDA
__device__ __forceinline__ void h2
#else
static inline void h2
#endif
( float r, float r2, float ri,//(4+ 10* )*74.5%
float rc, float r0, float rs, float & h ) {

    float k = rs*ri; k*=k;//k=(rs/r)^2
    h = rs*ri*ri*k*(TA+k*(TB+k*(TC+k*(TD+k*TE))));
}
#ifdef GBIS_CUDA
__device__ __forceinline__ void dh2
#else
static inline void dh2
#endif
( float r, float r2, float ri,//(4+ 11* )*74.5%
float rc, float r0, float rs, float & dh ) {

    float k = rs*ri; k*=k;//k=(rs/r)^2
    dh = -rs*ri*ri*ri*k*(DA+k*(DB+k*(DC+k*(DD+k*DE))));
}

#ifdef GBIS_CUDA
__device__ __forceinline__ void h3
#else
static inline void h3
#endif
( float r, float r2, float ri,//(3+ 5* 2/ 1log) 1.4%
float rc, float r0, float rs, float & h ) {
    float r2mrs2i = 1.f/(r2-rs*rs);
    h = 0.5f * ( rs*r2mrs2i + 0.5f * log((r-rs)/(r+rs))*ri );
}
#ifdef GBIS_CUDA
__device__ __forceinline__ void dh3
#else
static inline void dh3
#endif
( float r, float r2, float ri,//(5+ 8* 2/ 1log)*1.4%
float rc, float r0, float rs, float & dh ) {
    float rs2 = rs*rs;
    float r2mrs2i = 1.f/(r2-rs2);
    dh = -0.25f*ri*(2.f*(r2+rs2)*rs*r2mrs2i*r2mrs2i + ri*log((r-rs)/(r+rs)));
}

#ifdef GBIS_CUDA
__device__ __forceinline__ void h4
#else
static inline void h4
#endif
( float r, float r2, float ri,//(6+ 9* 2/ 1log)*0.4%
float rc, float r0, float rs, float & h ) {
    //float ri2 = ri*ri;
    float r02 = r0*r0;
    float rs2 = rs*rs;
    float r0i = 1.f/r0;
    float rspri = 1.f/(r+rs);
    float logr = log(r0*rspri);
    //float r02mrs2 = r02-rs2;
    float rilogr = ri*logr;
    h = 0.25f*( r0i*(2.f- 0.5f*(r0i*ri*(r2 + r02 - rs2))) - rspri + rilogr );
}
#ifdef GBIS_CUDA
__device__ __forceinline__ void dh4
#else
static inline void dh4
#endif
( float r, float r2, float ri,//(6+ 18* 2/ 1log)*0.4%
float rc, float r0, float rs, float & dh ) {
    float ri2 = ri*ri;
    float r02 = r0*r0;
    float rs2 = rs*rs;
    float r0i = 1.f/r0;
    float rspri = 1.f/(r+rs);
    float logr = log(r0*rspri);
    float r02mrs2 = r02-rs2;
    float rilogr = ri*logr;
    dh = 0.25f*( (- 0.5f +(r2*r02mrs2 - 2.f*r*rs*rs2+rs2*r02mrs2)
        * 0.5f *ri2*rspri*rspri)*r0i*r0i - ri*rilogr );
}

#ifdef GBIS_CUDA
__device__ __forceinline__ void h5
#else
static inline void h5
#endif
( float r, float r2, float ri,//(6+ 5* 3/ 1log)*0%, r<0.7Ang
float rc, float r0, float rs, float & h ) {
    float rs2 = rs*rs;
    float r2mrs2i = 1.f/(r2-rs2);
    float rsr2mrs2i = rs*r2mrs2i;
    float rprs = r+rs;
    float rmrs = r-rs;
    float logr = 0.5f*ri*log(-rmrs/rprs);
    h = 0.5f*( rsr2mrs2i + 2.f/r0 + logr );
}
#ifdef GBIS_CUDA
__device__ __forceinline__ void dh5
#else
static inline void dh5
#endif
( float r, float r2, float ri,//(5+ 8* 2/ 1log)*0%, r<0.7Ang
float rc, float r0, float rs, float & dh ) {
    float rs2 = rs*rs;
    float r2mrs2i = 1.f/(r2-rs2);
    float rsr2mrs2i = rs*r2mrs2i;
    float rprs = r+rs;
    float rmrs = r-rs;
    float logr = 0.5f*ri*log(-rmrs/rprs);
    dh = -0.5f*ri*((rs2+r2)*rsr2mrs2i*r2mrs2i+logr );
}

#ifdef GBIS_CUDA
__device__ __forceinline__ void h6
#else
static inline void h6
#endif
( float r, float r2, float ri,//0%, one atom within other
float rc, float r0, float rs, float & h ) {
  h = 0;
}
#ifdef GBIS_CUDA
__device__ __forceinline__ void dh6
#else
static inline void dh6
#endif
( float r, float r2, float ri,//0%, one atom within other
float rc, float r0, float rs, float & dh ) {
  dh = 0;
}

#ifdef GBIS_CUDA
__device__ __forceinline__ void CalcH 
#else
static inline void CalcH 
#endif
( float r, float r2, float ri,
float rc, float r0, float rs, float & h, int & d) {

/*
r - distance
rc - alpha cutoff
rs - screened radius
*/

if (r > 4*rs) { //change this to 1/4 r > rs
  if( r < rc - rs) {//II            68%
    h2(r,r2,ri,rc,r0,rs,h); d = 2;
  } else if (r < rc + rs) {//I      23%
    h1(r,r2,ri,rc,r0,rs,h); d = 1;
  } else /*if (r > rc + rs)*/ {//0  7%
    h0(r,r2,ri,rc,r0,rs,h); d = 0;
  } 
} else {
  if( r > r0 + rs ) {//III          1%
    h3(r,r2,ri,rc,r0,rs,h); d = 3;
  } else if ( r > (r0>rs?r0-rs:rs-r0) ) {//IV 0%
    h4(r,r2,ri,rc,r0,rs,h); d = 4;
  } else if (r0 < rs ) {//V         0%
    h5(r,r2,ri,rc,r0,rs,h); d = 5;
  } else {//VI                      0%
    h6(r,r2,ri,rc,r0,rs,h); d = 6;
  }
}
}
#ifdef GBIS_CUDA
__device__ __forceinline__ void CalcDH
#else
static inline void CalcDH
#endif
( float r, float r2, float ri,
float rc, float r0, float rs, float & dh, int & d) {
if (r > 4*rs) {
  if( r < rc - rs) {//II
    dh2(r,r2,ri,rc,r0,rs,dh); d = 2;
  } else if (r < rc + rs) {//I
    dh1(r,r2,ri,rc,r0,rs,dh); d = 1;
  } else /*if (r > rc + rs)*/ {//0
    dh0(r,r2,ri,rc,r0,rs,dh); d = 0;
  }
} else {
  if( r > r0 + rs ) {//III
    dh3(r,r2,ri,rc,r0,rs,dh); d = 3;
  } else if (r > (r0>rs?r0-rs:rs-r0) ) {//IV
    dh4(r,r2,ri,rc,r0,rs,dh); d = 4;
  } else if (r0 < rs ) {//V
    dh5(r,r2,ri,rc,r0,rs,dh); d = 5;
  } else {//VI
    dh6(r,r2,ri,rc,r0,rs,dh); d = 6;
  }
}
}
#ifdef GBIS_CUDA
__device__  __forceinline__ void CalcHPair
#else
static inline void CalcHPair
#endif
(
  float r,//distance
  float r2,//distance squared
  float ri,//inverse distance
  float rc,//cutoff
  float ri0,
  float rjs,
  float rj0,
  float ris,
  int & dij,//domain 1
  int & dji,//domain 2
  float & hij,//output
  float & hji//output
) {
  CalcH(r,r2,ri,rc,ri0,rjs,hij,dij);//hij
  CalcH(r,r2,ri,rc,rj0,ris,hji,dji);//hji
}
#ifdef GBIS_CUDA
__device__ __forceinline__ void CalcDHPair
#else
static inline void CalcDHPair
#endif
( float r,//distance
  float r2,
  float ri,
  float rc,//cutoff
  float ri0,
  float rjs,
  float rj0,
  float ris,
  int & dij,//domain 1
  int & dji,//domain 2
  float & dhij,
  float & dhji
) {
//                  swapped
  CalcDH(r,r2,ri,rc,ri0,rjs,dhij,dij);//hij
  CalcDH(r,r2,ri,rc,rj0,ris,dhji,dji);//hji
}

/*
 * Calculate GB Energy, GB dEdr force
 * also output intermediate values used in dEda
 */
#ifdef GBIS_CUDA
__device__ __forceinline__ void Calc_dEdr_Pair
#else
static inline void Calc_dEdr_Pair
#endif
(//no longer does i==j
  const float & r,
  const float & r2,
  const float & qiqj,
  const float & ai,
  const float & aj,
  const float & kappa,
  const float & epsilon_p_i,
  const float & epsilon_s_i,
  float & aiaj,
  float & expr2aiaj4,
  float & fij,
  float & f_i,
  float & expkappa,
  float & Dij,
  float & gbE,   //return
  float & ddrGbE //return
) {
  //allocate local variables
  float aiaj4,ddrDij,ddrf_i,ddrfij;

  //calculate GB energy
  aiaj = ai*aj;
  aiaj4 = 4.f*aiaj;
  //printf("exp(%e)\n",(-r2/aiaj4));
  expr2aiaj4 = exp(-r2/aiaj4);
  fij = sqrt(r2+aiaj*expr2aiaj4);
  f_i = 1/fij;
  expkappa = (kappa > 0.f) ? exp(-kappa*fij) : 1.f;
  Dij = epsilon_p_i - expkappa*epsilon_s_i;
  //gbE = -COULOMB*qiqj*Dij*f_i;
  gbE = qiqj*Dij*f_i;

  //calculate energy derivatives
  ddrfij = r*f_i*(1.f - 0.25f*expr2aiaj4);
  ddrf_i = -ddrfij*f_i*f_i;
  ddrDij = kappa*expkappa*ddrfij*epsilon_s_i;
  //ddrGbE = -COULOMB*qiqj*(ddrDij*f_i+Dij*ddrf_i);
  ddrGbE = qiqj*(ddrDij*f_i+Dij*ddrf_i);
}

/*
 * Calculate summation element of dEda array
 * must calculate dEdr previously to retreive intermediate values
 */
#ifdef GBIS_CUDA
__device__ __forceinline__ void Calc_dEda_Pair
#else
static inline void Calc_dEda_Pair
#endif
( const float & r2,
  const float & ai,
  const float & aj,
  const float & qiqj,
  const float & kappa,
  const float & aiaj,
  const float & expkappa,
  const float & expr2aiaj4,
  const float & fij,
  const float & f_i,
  const float & Dij,
  const float & epsilon_s_i,
  float & dEdai,//return
  float & dEdaj //return
) {

  //float tmp_dEda = -0.5*COULOMB*qiqj*f_i*f_i
  float tmp_dEda = 0.5f*qiqj*f_i*f_i
                      *(kappa*epsilon_s_i*expkappa-Dij*f_i)
                      *(aiaj+0.25f*r2)*expr2aiaj4;//0
  dEdai = tmp_dEda/ai;
  dEdaj = tmp_dEda/aj;
}

/*
 * Calculate Coulomb and GB interaction and dEda element
 * for a pair of atoms
 */
#ifdef GBIS_CUDA
__device__ __forceinline__ void Phase2_Pair
#else
static inline void Phase2_Pair
#endif
(//doesn't do self energies

//input values
  const float & r,
  const float & r2,
  const float & r_i,
  const float & qiqj,
  const float & ai,
  const float & aj,
  const float & epsilon_p_i,
  const float & epsilon_s_i,
  const float & kappa,
  const int & doFullElect,

//return values
  float & gbEij,
  float & ddrGbEij,
  float & dEdai,
  float & dEdaj
) {

  //calculate GB energy and force
  float aiaj,expr2aiaj4,fij,f_i,expkappa,Dij;
  Calc_dEdr_Pair(r,r2,qiqj,ai,aj,
      kappa,epsilon_p_i,epsilon_s_i,
      aiaj,expr2aiaj4,fij,f_i,expkappa,
      Dij,gbEij,ddrGbEij);

  //calculate dEda
  if (doFullElect) {
    Calc_dEda_Pair(r2,ai,aj,qiqj,kappa,
              aiaj,expkappa,expr2aiaj4,
              fij,f_i,Dij,epsilon_s_i,dEdai,dEdaj);
  } else {
    dEdai = 0.f;
    dEdaj = 0.f;
  }
}

#if 0
static inline void init_gbisTable (
float **tablePtr,
float kappa,
float maxX,
int numEntriesPerX
) {
  float *table = *tablePtr;//dereference
  float minX = 0;
  int numPts = (maxX-minX) * numEntriesPerX;
  int numVals = 3;
  /*
  table = (float*) malloc(numVals*numPts*sizeof(float));
  for (int i = 0; i < numPts; i++) {
    float x = (1.0*i) / numEntriesPerX + minX; 
      bornRadJ = gbisParams->bornRad[1][j];
      aiaj = bornRadI * bornRadJ;
      aiaj4 = 4.0 * aiaj;
      expr2aiaj4 = exp(-r2/aiaj4);
      fij = sqrt(r2+aiaj*expr2aiaj4);
      f_i = 1.0/fij;
      expkappa = (kappa > 0.0) ? exp(-kappa*fij) : 1.0;
  }
*/
}
#endif

//include once
#endif
