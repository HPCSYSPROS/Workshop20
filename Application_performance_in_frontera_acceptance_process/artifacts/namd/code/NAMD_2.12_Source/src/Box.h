/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef BOX_H
#define BOX_H

// #define BOX_DEBUG

#include "OwnerBox.h"

template <class Owner, class Data> class Box {

  friend class OwnerBox<Owner,Data>;

  private:

    Box(OwnerBox<Owner,Data>* o, int n = -1):
#ifdef BOX_DEBUG
      state(CLOSED), user(n),
#endif
      ownerBox(o) {}

    ~Box(void) {};

#ifdef BOX_DEBUG
    enum box_state {OPEN, CLOSED} state;
#endif

    OwnerBox<Owner,Data> *ownerBox;	

  public:
#ifdef BOX_DEBUG
  int user;
#endif

    Data* open(void) {
      ownerBox->openCount--;
#ifdef BOX_DEBUG
      if (ownerBox->openCount < 0) NAMD_bug("too many boxes opened");
      if (state != CLOSED) NAMD_bug("box re-opened");
      state = OPEN; 
#endif
      return ownerBox->data; 
    }

    void close(Data ** const t) {
      *t = NULL;
      ownerBox->closeCount--;
#ifdef BOX_DEBUG
      if (ownerBox->closeCount < 0) NAMD_bug("too many boxes closed");
      if (state != OPEN) NAMD_bug("box re-closed");
      state = CLOSED;
#endif
      // Trigger callback!
      if ( ! ownerBox->closeCount ) {
	ownerBox->close();
      }
    }

    void skip(void) {
      ownerBox->openCount--;
      ownerBox->closeCount--;
#ifdef BOX_DEBUG
      if (state != CLOSED) NAMD_bug("box skipped while open");
      if (ownerBox->openCount < 0) NAMD_bug("too many boxes opened");
      if (ownerBox->closeCount < 0) NAMD_bug("too many boxes closed");
#endif
      if ( ! ownerBox->closeCount ) {
        ownerBox->close();
      }
    }

};

#endif // BOX_H
