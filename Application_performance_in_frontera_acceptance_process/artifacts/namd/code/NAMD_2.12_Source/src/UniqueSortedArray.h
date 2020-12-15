/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef USRTARRAY_H
#define USRTARRAY_H

#include "SortedArray.h"

template <class Elem> class UniqueSortedArray : public SortedArray<Elem> {

  public:

    UniqueSortedArray(int s=0) : SortedArray<Elem>(s) { }

    UniqueSortedArray(UniqueSortedArray<Elem> &ua) : 
      SortedArray<Elem>(ua) { }

    UniqueSortedArray(SortedArray<Elem> &sa) : SortedArray<Elem>(sa) { 
      this->uniq(); 
    }

/*
    UniqueSortedArray(ResizeArray<Elem> &ra) : SortedArray<Elem>(ra) {
      this->uniq();
    }
*/
  
    UniqueSortedArray<Elem>& operator =(UniqueSortedArray<Elem> & ua) {
      SortedArray<Elem>::operator=(ua);
      return(*this);
    }
  
    UniqueSortedArray<Elem>& operator =(SortedArray<Elem> &sa) {
      SortedArray<Elem>::operator=(sa);
      this->uniq();
      return(*this);
    }

/*
    UniqueSortedArray<Elem>& operator =(ResizeArray<Elem> &ra) {
      SortedArray<Elem>::operator=(ra);
      this->uniq();
      return(*this);
    }
*/
  
    int add(const Elem& elem) { return(insert(elem)); }

    inline int insert(const Elem& elem);

};

template <class Elem>
inline int 
UniqueSortedArray<Elem>::insert(const Elem& elem) {
  int found = this->bsearch(elem);
  if (found == -1) {
    return ResizeArray<Elem>::insert(elem, 0);
  }
  if (found < this->size() && this->rep[found] == elem) {
    return -2;
  }
  if (found == (this->size()-1) && this->rep[found] < elem) {
    return ResizeArray<Elem>::insert(elem, this->size());
  } else {
    return ResizeArray<Elem>::insert(elem, found);
  }
}
#endif
