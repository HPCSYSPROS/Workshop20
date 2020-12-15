/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef OWNERBOX_H
#define OWNERBOX_H

#include "charm++.h"

template <class Owner, class Data> class Box;

template <class Owner, class Data> class OwnerBox {

  friend class Box<Owner,Data>;

  int bid, type;

  public:

    OwnerBox(Owner *o, void (Owner::*fn)() ) :
      bid(-1), type(-1),
      owner(o), callback(fn), data(0),
      numberUsers(0), openCount(0), closeCount(0) {};

   OwnerBox(Owner *o, void (Owner::*fn)(), int id ,int tp) :
      bid(id), type(tp),
      owner(o), callback(fn), data(0),
      numberUsers(0), openCount(0), closeCount(0) {};


    ~OwnerBox(void) {
      if (numberUsers) {
        CkPrintf("OwnerBox::~OwnerBox() - still have boxes out there!\n");
      }
    }
        
    void open(Data* d) {
      closeCount = openCount = numberUsers;
      data = d;
      if ( ! closeCount ) close();
    }
  
    inline void close(void);

    inline Box<Owner,Data> *checkOut(int id);

    inline void checkIn(Box<Owner,Data> * box);
  
    int isOpen() {
      return (closeCount != numberUsers || openCount != numberUsers);
    }

    // manipulate counters without Box object
    void clientAdd();
    void clientRemove();
    Data* clientOpen(int count=1) {
      openCount -= count;
      return data;
    }
    void clientClose(int count=1) {
      if ( ! (closeCount -= count) ) close();
    }

  private:
    Owner *owner;
    void (Owner::*callback)(void);
    Data* data;
    int numberUsers, openCount, closeCount;
};

template <class Owner, class Data>
inline void OwnerBox<Owner,Data>::clientAdd(void) {
  if (closeCount != numberUsers || openCount != numberUsers) {
    CkPrintf("OwnerBox::clientAdd() while in use\n");
  }
  ++numberUsers; ++closeCount; ++openCount; 
}

template <class Owner, class Data>
inline Box<Owner,Data> *OwnerBox<Owner,Data>::checkOut(int id) {
  clientAdd();
  return (new Box<Owner,Data>(this,id));
}

template <class Owner, class Data>
inline void OwnerBox<Owner,Data>::clientRemove() {
  if (closeCount != numberUsers || openCount != numberUsers) {
    CkPrintf("OwnerBox::clientRemove() while in use\n");
  }
  if ( ! numberUsers-- ) {
    CkPrintf("OwnerBox::clientRemove() - no registrants remaining\n");
    numberUsers = 0;
  } else {
    closeCount--; openCount--;
  }
}

template <class Owner, class Data>
inline void OwnerBox<Owner,Data>::checkIn(Box<Owner,Data> * box) {
  delete box;
  clientRemove();
}

template <class Owner, class Data>
inline void OwnerBox<Owner,Data>::close(void) {
  if (!closeCount && !openCount) {
    data = 0; closeCount = openCount = numberUsers;
    (owner->*callback)();
  } else {
    CkPrintf("OwnerBox::close() - close called, but closeCount %d openCount %d\n", closeCount, openCount);
  }
}

#endif
