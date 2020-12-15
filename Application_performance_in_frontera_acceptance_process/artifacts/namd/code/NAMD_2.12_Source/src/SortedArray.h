/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef SRTARRAY_H
#define SRTARRAY_H

#include "SortableResizeArray.h"

template <class Elem> class SortedArray: public SortableResizeArray<Elem> {

  protected:

    int isSorted;

  public:

    SortedArray(void) : SortableResizeArray<Elem>() { 
      isSorted = 1;
    }

    SortedArray(int s) : SortableResizeArray<Elem>(s) { 
      isSorted = 1;
    }

    SortedArray(SortedArray<Elem> &sa) : SortableResizeArray<Elem>(sa) { 
      if(!(isSorted = sa.isSorted)) sort();
      isSorted = 1;
    }

    SortedArray(SortableResizeArray<Elem> &ra) : 
        SortableResizeArray<Elem>(ra) {
      sort(); isSorted = 1;
    }

    SortedArray<Elem>& operator =(SortedArray<Elem> & sa) {
      SortableResizeArray<Elem>::operator=(sa);
      isSorted = sa.isSorted;
      return(*this);
    }

    SortedArray<Elem>& operator =(SortableResizeArray<Elem> &ra) {
      SortableResizeArray<Elem>::operator=(ra);
      sort(); isSorted = 1;
      return(*this);
    }

    int load(const Elem& elem) {
      isSorted = 0;
      return(ResizeArray<Elem>::add(elem));
    }

    int add(const Elem& elem) {
      return(insert(elem));
    }

    int del(const Elem & elem) {
      int found = bsearch(elem);
      if (this->size() != 0 && this->rep[found] == elem) {
        return(SortableResizeArray<Elem>::del(found,1));
      } else {
        return(-1);
      }
    }

    void sort(void) { SortableResizeArray<Elem>::sort(); isSorted = 1; }

    int bsearch(const Elem& elem) { 
      if (!isSorted) sort();
      return (SortableResizeArray<Elem>::bsearch(elem));
    }

    inline int insert(const Elem& elem);

    int index(const Elem& elem) { return bsearch(elem); }

    inline Elem *find(const Elem& elem);
};

template <class Elem>
inline int SortedArray<Elem>::insert(const Elem& elem) {
  int found = bsearch(elem);
  if (found == -1) {
    return (ResizeArray<Elem>::insert(elem, 0));
  }
  if (found == (this->size()-1) && this->rep[found] < elem) {
    return (ResizeArray<Elem>::insert(elem, this->size()));
  } else {
    return (ResizeArray<Elem>::insert(elem, found));
  }
}

template <class Elem>
inline Elem * SortedArray<Elem>::find(const Elem& elem) {
  int found = bsearch(elem);
  if ( found < 0 || found == this->size() ) 
    return ((Elem *)NULL);
  if (this->rep[found] == elem) {
    return (&(this->rep[found]));
  } else {
    return ((Elem *)NULL);
  }
}

#endif
