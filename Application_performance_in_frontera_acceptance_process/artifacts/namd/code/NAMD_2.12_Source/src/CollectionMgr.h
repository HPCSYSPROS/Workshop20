/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COLLECTIONMGR_H
#define COLLECTIONMGR_H

#include "charm++.h"

#include "main.h"
#include "NamdTypes.h"
#include "BOCgroup.h"
#include "PatchMap.h"
#include "ProcessorPrivate.h"
#include "CollectionMgr.decl.h"


class SlaveInitMsg : public CMessage_SlaveInitMsg
{
public:
  CkChareID master;
};

class CollectionMgr : public CBase_CollectionMgr
{
public:

  static CollectionMgr *Object() { 
    return CkpvAccess(CollectionMgr_instance); 
  }
  CollectionMgr(SlaveInitMsg *msg);
  ~CollectionMgr(void);

  void submitPositions(int seq, FullAtomList &a, Lattice l, int prec);
  void submitVelocities(int seq, int zero, FullAtomList &a);
  void submitForces(int seq, FullAtomList &a, int maxForceUsed, ForceList *f);
  void sendDataStream(const char *);
  
  CkChareID getMasterChareID() { return master; }  

  class CollectVectorInstance
  {
  public:

    CollectVectorInstance(void) : seq(-10) { ; }

    void free() { seq = -10; }
    int notfree() { return ( seq != -10 ); }

    void reset(int s, int p) {
      if ( s == -10 ) NAMD_bug("seq == free in CollectionMgr");
      seq = s;
      precisions = p;
      remaining = PatchMap::Object()->numHomePatches();
      aid.resize(0);
#ifdef MEM_OPT_VERSION
      outRank.resize(0);
#endif
      data.resize(0);
      fdata.resize(0);
    }

    // true -> send it and delete it!
#ifdef MEM_OPT_VERSION
    int append(AtomIDList &a, ResizeArray<int> &oRank, ResizeArray<Vector> &d)
#else
    int append(AtomIDList &a, ResizeArray<Vector> &d)
#endif
    {
      int size = a.size();
      for( int i = 0; i < size; ++i )
      {
	aid.add(a[i]);
#ifdef MEM_OPT_VERSION
    outRank.add(oRank[i]);
#endif
	if ( precisions & 2 ) data.add(d[i]);
	if ( precisions & 1 ) fdata.add(d[i]);
      }
      return ( ! --remaining );
    }

    int seq;
    AtomIDList aid;
#ifdef MEM_OPT_VERSION
    //record the dest output proc rank the atom goes to
    ResizeArray<int> outRank;
#endif
    int precisions;
    ResizeArray<Vector> data;
    ResizeArray<FloatVector> fdata;

  private:
    int remaining;

  };

  class CollectVectorSequence
  {
  public:
#ifdef MEM_OPT_VERSION
    CollectVectorInstance* submitData(int seq, AtomIDList &i, ResizeArray<int> &oRank,
                                      ResizeArray<Vector> &d, int prec=2)
#else
    CollectVectorInstance* submitData(int seq, AtomIDList &i, 
                                      ResizeArray<Vector> &d, int prec=2)
#endif
    {
      CollectVectorInstance **c = data.begin();
      CollectVectorInstance **c_e = data.end();
      for( ; c != c_e && (*c)->seq != seq; ++c );
      if ( c == c_e )
      {
       c = data.begin();
       for( ; c != c_e && (*c)->notfree(); ++c );
       if ( c == c_e ) {
	data.add(new CollectVectorInstance);
	c = data.end() - 1;
       }
       (*c)->reset(seq,prec);
      }
#ifdef MEM_OPT_VERSION
      if ( (*c)->append(i, oRank, d) )
#else
      if ( (*c)->append(i,d) )
#endif
      {
        return *c;
      }
      else
      {
        return 0;
      }
    }

    ResizeArray<CollectVectorInstance*> data;

  };
private:

  CkChareID master;

  CollectVectorSequence positions;
  CollectVectorSequence velocities;
  CollectVectorSequence forces;

};

#endif

