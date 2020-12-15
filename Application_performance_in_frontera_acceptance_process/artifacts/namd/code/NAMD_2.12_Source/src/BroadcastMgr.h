/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Coordinates broadcast of a data type from a Controller/Seq
   to all other Controller/Sequencer type objects (they must
   run in a thread!)
*/

#include "charm++.h"
#include "main.h"
#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "ProcessorPrivate.h"
#include "BroadcastMgr.decl.h"

#ifndef _BCASTMGR_H
#define _BCASTMGR_H

#define BCASTMSGSIZE (9*sizeof(double))

class BroadcastMsg : public CMessage_BroadcastMsg {
friend class BroadcastMgr;
public:
  ~BroadcastMsg() { }
  BroadcastMsg() { }

private:
  // Only seen by BroadcastMgr
  char msg[BCASTMSGSIZE];
  int size;
  int id;
  int tag;
  int node;
};

class BroadcastClient;

class BroadcastClientElem {
public:
  BroadcastClientElem() {}
  BroadcastClientElem(BroadcastClient * c) : broadcastClient(c) {}
  ~BroadcastClientElem() {}

  BroadcastClient *broadcastClient;

  size_t hash() const { return (size_t)broadcastClient; }
  int operator==(const BroadcastClientElem &b) const { 
    return broadcastClient == b.broadcastClient; 
  }
};

class TaggedMsg {
public:
  TaggedMsg() {}
  TaggedMsg(int t) : tag(t) {}
  TaggedMsg(int t, int s, int c, void *m) 
    : tag(t), counter(c), msgSize(s) { memcpy((void*)msg,m,s); }
  ~TaggedMsg() {}

  int tag;
  int counter;
  int msgSize;
  char msg[BCASTMSGSIZE];

  int hash() const { return tag; }
  int operator==(const TaggedMsg &tm) const { return(tag == tm.tag); }
};

class BOID {
public:
  BOID() {}
  BOID(int id) { this->id = id; }
  ~BOID() {}

  int hash() const { return id; }
  int operator==(const BOID &b) const { return id == b.id; }
  int id;

  UniqueSet<BroadcastClientElem> *broadcastSet;
  UniqueSet<TaggedMsg> *taggedMsg;
};

class BroadcastMgr : public CBase_BroadcastMgr
{
public:
  BroadcastMgr() { 
    CkpvAccess(BroadcastMgr_instance) = this; 
  }
  ~BroadcastMgr(void);
	  
  // Singleton Access method
  inline static BroadcastMgr *Object() {
    return CkpvAccess(BroadcastMgr_instance);
  }

  int getbuf(BroadcastClient &b, int tag, void* msg);
  void send(BroadcastClient &b, int tag, void *buf, size_t);
  void subscribe(BroadcastClient &bc);
  void unsubscribe(BroadcastClient &bc);
  void recvBroadcast(BroadcastMsg *msg);

private:
  UniqueSet<BOID> boid;
};

#endif /* _BCASTMGR_H */

