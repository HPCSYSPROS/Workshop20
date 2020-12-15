/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   UniqueSet template - (hashtable)
*/

#ifndef UNIQUESET_H
#define UNIQUESET_H

#include "UniqueSetRaw.h"

// Need this juju to use templated friend below
template <class Type> class UniqueSetIter;

template <class Elem> class UniqueSet {

  friend class UniqueSetIter<Elem>;

  private:

    UniqueSetRaw<Elem> *rep;

  public:

    // Various Constructors
    UniqueSet(void) { rep = new UniqueSetRaw<Elem>(); rep->refCount = 1; }

    UniqueSet(int size) { rep=new UniqueSetRaw<Elem>(size); rep->refCount=1; }

    UniqueSet(const UniqueSet<Elem> &us) {
      if (!--rep->refCount) delete rep; rep = us.rep; rep->refCount++; 
    }

    UniqueSet(const UniqueSet<Elem>* us) { 
      if (!--rep->refCount) delete rep; rep = new UniqueSetRaw<Elem>(*us->rep);
      rep->refCount = 1;
    }

    ~UniqueSet(void) { if (!--rep->refCount) delete rep; }

    UniqueSet<Elem>& operator =(const UniqueSet<Elem>& us) {
      if (!--rep->refCount) delete rep; rep = us.rep; rep->refCount++; 
      return (*this);
    }

    void rehash(void) { rep->rehash(); }
  
    int add(const Elem &elem) { return rep->add(elem); }

    int load(const Elem &elem) { return rep->load(elem); }

    int del(const Elem &elem) { return rep->del(elem); }

    int size(void) const { return rep->size(); }

    Elem *find(const Elem &elem) { return rep->find(elem); }

    void clear(void) { rep->clear(); }

#ifdef DEBUG
    void status(void) { rep->status(); }
#endif
};

#endif
