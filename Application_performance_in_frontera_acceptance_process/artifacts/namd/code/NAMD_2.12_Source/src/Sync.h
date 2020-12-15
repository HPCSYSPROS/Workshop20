
/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Master BOC.  coordinates startup, close down of each PE
   Also owns pointers to common objects needed by system		
   Many utility static methods are owned by Node.
*/

#ifndef _SYNC_H
#define _SYNC_H

#include "charm++.h"

#include "ProcessorPrivate.h"
#include "Sync.decl.h"

class Sync : public CBase_Sync
{
private:
    struct _clist {
    int pid;
    int step;
    Compute **cbegin;
    Compute **cend;
    int doneMigration;
    } *clist;
    const int INCREASE;
    int capacity;

    int useSync, useProxySync;

    int step;
    int counter;
    int cnum;
    int nPatcheReady;
    int numPatches;

    char homeReady;

    void releaseComputes();
    void triggerCompute();
public:
    Sync(void);
    ~Sync(void);
    inline static Sync *Object() { return CkpvAccess(Sync_instance); }
    void openSync(); 
    int holdComputes(PatchID pid, Compute **cbegin, Compute **cend, int doneMigration, int seq);
    void PatchReady(void);
};

#endif
