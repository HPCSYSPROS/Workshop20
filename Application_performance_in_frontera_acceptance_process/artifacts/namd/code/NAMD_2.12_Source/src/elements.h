/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef ELEMENTS_DEFS_H
#define ELEMENTS_DEFS_H

#include <lbdb.h>
#include "Set.h"

class InfoRecord {
public:
   double load;
   int Id; // should replace other Ids.
};


class computeInfo : public InfoRecord {
public: 
   /*   int computeId; replaced by Id */
   int patch1, patch2;
   int processor; // caller to ReBalancer MAY leave this field -1, 
   int oldProcessor; // stores the current assignment of the compute object.
   LDObjHandle handle;
};

class patchInfo : public InfoRecord {
public:
   int processor;
   int numAtoms;
   IRSet proxiesOn;  // caller to ReBalancer should fill in the forced proxies
};

class processorInfo: public InfoRecord {
public:
   // int processorNum; replaced by inherited "Id".
   double backgroundLoad; // background work pre-assigned to the processor.
   double idleTime;       // idle time
   double computeLoad;    //load due to computes. The total load is computed
                          // by adding these two.		     
   // Added 10/22/01:  indicate if this processor will migrate its objs.   
   bool  available;
   LargeIRSet patchSet;   // caller to ReBalancer should leave this field NULL.
   LargeIRSet proxies;    // caller to ReBalancer should fill in the forced proxies
   LargeIRSet computeSet; // caller to ReBalancer should leave this field NULL.
   
   // Added 4-29-98: Array to keep track of number of computes that are using
   // each proxy on a processor
   // unsigned char *proxyUsage;
public:
   processorInfo(): backgroundLoad(0.), idleTime(0.), computeLoad(0.), available(true) {}
};

#endif
