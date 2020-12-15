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

#ifndef _BCASTOBJ_H
#define _BCASTOBJ_H

#include "BroadcastMgr.h"
#include "BroadcastClient.h"
#include "LdbCoordinator.h"
#include "common.h"

template <class T> class SimpleBroadcastObject : public BroadcastClient {

  public:

    const LDObjHandle *ldObjPtr;

    SimpleBroadcastObject(int id, const LDObjHandle *h = 0) : BroadcastClient(id), ldObjPtr(h) {
      if ( sizeof(T) > BCASTMSGSIZE ) {
        NAMD_bug("SimpleBroadcastObject instantiated on class larger than BCASTMSGSIZE");
      }
    }
    ~SimpleBroadcastObject() { }

    T get(int tag) {
      T tmp;
      if ( ldObjPtr ) LdbCoordinator::Object()->pauseWork(*ldObjPtr);
      while ( BroadcastMgr::Object()->getbuf(*this, tag, (void*)(&tmp)) < 0 ) {
        suspendFor(tag);
      }
      if ( ldObjPtr ) LdbCoordinator::Object()->startWork(*ldObjPtr);
      return tmp;
    }
    
    void publish(int tag,const T &t ) {
      BroadcastMgr::Object()->send(*this, tag, (void*)(&t), sizeof(T));
    }

};

#endif
