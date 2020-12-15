/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "UniqueSet.h"
#include "UniqueSetIter.h"
#include "BroadcastMgr.decl.h"
#include "BroadcastMgr.h"
#include "BroadcastClient.h"
#include "BroadcastObject.h"
#include "ProcessorPrivate.h"
#define MIN_DEBUG_LEVEL 3
//#define DEBUGM
#include "Debug.h"

BroadcastMgr::~BroadcastMgr(void) {
  UniqueSetIter<BOID> boidIter(boid);
  for (boidIter = boidIter.begin(); boidIter != boidIter.end(); boidIter++) {
    delete boidIter->broadcastSet;
    if (boidIter->taggedMsg) {
      delete boidIter->taggedMsg;
    }
  }
}


int
BroadcastMgr::getbuf(BroadcastClient &b, int tag, void *msg) {
  int rval = -1;
  TaggedMsg *tm;
  BOID* boidTmp = boid.find(BOID(b.id));
  if (!boidTmp) {
    return(-2);
  }
  if ( (tm = (boidTmp->taggedMsg->find(TaggedMsg(tag)))) ) {
    rval = tm->msgSize;
    memcpy(msg, tm->msg, tm->msgSize);
    if (!--(tm->counter)) {
      (boid.find(BOID(b.id)))->taggedMsg->del(TaggedMsg(tag));
    }
  }
  return(rval);
}


void 
BroadcastMgr::send(BroadcastClient &b, int tag, void *buf, size_t size) {
  BroadcastMsg* msg = new BroadcastMsg;
  memcpy((void*)(msg->msg),buf,size);
  msg->size = (int)size;
  msg->tag = tag;
  msg->id = b.id;
  msg->node = CkMyPe();
  CProxy_BroadcastMgr(thisgroup).recvBroadcast(msg);
}

void 
BroadcastMgr::subscribe(BroadcastClient &bc) {
  BOID *b;
  if (!(b = boid.find(BOID(bc.id)))) {
    boid.add(BOID(bc.id));
    b = boid.find(BOID(bc.id));
    b->broadcastSet = new UniqueSet<BroadcastClientElem>;
    b->taggedMsg = new UniqueSet<TaggedMsg>;
  }
  b->broadcastSet->add(BroadcastClientElem(&bc));
}

void 
BroadcastMgr::unsubscribe(BroadcastClient &bc) {
  BOID *b;
  if ( (b = boid.find(BOID(bc.id))) ) {
    b->broadcastSet->del(BroadcastClientElem(&bc));
    if (!b->broadcastSet->size()) {
      delete b->broadcastSet;
      b->broadcastSet = 0;
      delete b->taggedMsg;
      b->taggedMsg = 0;
    }
  }
}

void 
BroadcastMgr::recvBroadcast(BroadcastMsg *msg) {
  BOID *b;
  int counter;
  // Check if msg->id has any registrants
  if ( (b = boid.find(BOID(msg->id))) ) {
    // add message to taggedMsg container
    counter = b->broadcastSet->size();
    if (msg->node == CkMyPe()) counter--; // get rid of sender
    if ( counter < 0 ) NAMD_bug("BroadcastMgr::recvBroadcast counter < 0");
    else if ( counter > 0 ) {
      b->taggedMsg->add(TaggedMsg(msg->tag,msg->size,counter,msg->msg));

      // inform all registrants of mew message
      UniqueSetIter<BroadcastClientElem> bcIter(*(b->broadcastSet));
      for (bcIter = bcIter.begin(); bcIter != bcIter.end(); bcIter++) {
        bcIter->broadcastClient->awaken(msg->id, msg->tag);
      }
    }
  }
  delete msg;
}

#include "BroadcastMgr.def.h"

