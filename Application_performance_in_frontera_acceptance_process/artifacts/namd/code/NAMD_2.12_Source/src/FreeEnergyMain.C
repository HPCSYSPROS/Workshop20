/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <string.h>
#include <iostream.h>
#ifndef WIN32
#include <strstream.h>
#else
#include <strstrea.h>
#endif
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyGroup.h"
#include "Vector.h"
#include "FreeEnergyVector.h"
#include "FreeEnergyRestrain.h"
#include "FreeEnergyRMgr.h"
#include "FreeEnergyLambda.h"
#include "FreeEnergyLambdMgr.h"
#include "GlobalMasterFreeEnergy.h"

void main() {

  GlobalMasterFreeEnergy  CFE;

  char  AStr[] =
  "urestraint { \
     posi (sn, 1) kf=20 ref=(,   0,0,0) \
     dist (sn, 1) (sn, 2) kf= 10, ref = 10 \
     angl (sn, 1) (sn, 2) (sn, 3) kf=5, ref=180 \
     dihe (sn,1)(sn 2)(sn,3)(sn 4),barr=10 ref=90 \
     posi bound group { (sn1, 1) to (sn1, 3) }, kf=10, low = (1,2, 3 4) \
     dist bound (sn1, 1, cb) (sn1, 2, cb) kf=5, hi=10 \
     angle bound group { (sn, 1) (sn, 3) } (sn, 5), (sn, 6), kf = 3,low=90 \
     dihe bound (sn, 1) (sn,2) (sn 3),(sn ,4) gap=5, low=5, hi = 10, delta= 2 \
     posi pmf (sn, 1, ca) kf=10, low=(1,2,3), high = ( 4 5 6) \
     dist pmf (sn, 1, ca) (sn, 2, ca) kf = 10  low = 3  hi=2 \
     angle pmf (sn 1 ca) (sn 2 ca) group {(sn 3 ca) (sn 4 ca)} kf=1 low=90 hi=180 \
     dihe pmf (sn 1) (sn 2) (sn,3,ca) (sn,4,ca) barr= 5, low = 45  hi=135 \
   } \
   urestraint { \
     dist pmf group { (insulin, 10) to (insulin, 15) } \
              group { (insulin, 30) to (insulin, 35) } kf=20, low=20, hi=10 \
   } \
   pmf { \
     task = grow, \
     time = 30 fs \
     print = 8.fs\
   } \
   mcti { \
     task = up, \
     time = 30fs \
     print = 4fs \
   } \
   urestraint { \
     dist (insulin, 10, ca) (insulin, 20, cb) kf = 10 ref=3; \
   } \
   urestraint { \
     angle (insulin,10) (insulin,20) (insulin,30) kf=10 ref=3; \
     dist pmf (sn, 1, ca) (sn, 2, ca) kf = 10  low = 3  hi=2 \
   } \
   pmf { \
     task = down \
   } \
   pmf { \
     task = fade \
   } \
   mcti{\
     task=stepup; \
     equiltime=12 fs \
     accumtime = 12fs \
     numsteps \
     3 \
   } \
   pmf { \
     task = stop \
     time = 200 fs \
     print = 25 fs \
   } \
   pmf { \
     task = nogrow \
     time = 205 fs \
     print = 27 fs \
   } \
   pmf { \
     task = stop \
     lambda = 0.88 \
   } \
   pmf { \
     task = nogrow \
     lambda = 0.55 \
   } \
   pmf { \
     task = nogrow \
     lambdat=.44 \
   } \
   pmf { \
     task = stop \
     lambda=.33 \
   } \
   pmf { \
     task = up \
     time = 3ps \
     lambda = 1.0 \
     lambdat = 0.5 \
     print = 40.fs \
     accumtime = 100. fs \
     equiltime = 4 \
   } \
   ";

  // put everything into the istrstream
  istrstream* pInput;
  pInput = new istrstream(AStr);
  CFE.SetConfig(pInput);
  CFE.user_initialize();
  
  // this is what's executed every time-step
  int  SomeNum = 10000;
  for (int i=0; i<SomeNum; i++) {
    CFE.user_calculate();
  }

}
