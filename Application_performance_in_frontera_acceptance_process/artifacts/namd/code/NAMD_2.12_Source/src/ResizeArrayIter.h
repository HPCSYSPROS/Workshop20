/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef RAITER_H
#define RAITER_H

#include "ResizeArray.h"

// Don't use for speed - use iter if later we will probably want
// to use a better container class for better space or algorithm behavior

template <class T> class ResizeArrayIter {
  private:

    ResizeArray<T> *array;
    int currentIndex;

  public:

    T *operator->(void) { return ((array->rep.array)+currentIndex); }

    ResizeArrayIter(void) {
      array = NULL;
      currentIndex = 0;
    }

    ResizeArrayIter(ResizeArray<T>& ra) {
      array = &ra;
      currentIndex = 0;
    }

    ResizeArrayIter(const ResizeArrayIter<T>& iter) {
      array = iter.array;
      currentIndex = iter.currentIndex;
    }

    ResizeArrayIter<T>& operator= (const ResizeArrayIter<T>& iter) {
      array = iter.array;
      currentIndex = iter.currentIndex;
      return (*this);
    }

    ~ResizeArrayIter(void) {}

    ResizeArrayIter<T> begin(void) const {
      ResizeArrayIter<T> iter;
      iter.array = array;
      iter.currentIndex = 0;
      return(iter);
    }

    ResizeArrayIter<T> end(void) const {
      ResizeArrayIter<T> iter;
      iter.array = array;
      iter.currentIndex = array->size();
      return(iter);
    }
        
    int operator!= (const ResizeArrayIter<T> &iter) const {
      return (iter.currentIndex != currentIndex || iter.array != array);
    }

    int operator== (const ResizeArrayIter<T> &iter) const {
      return (!operator!=(iter));
    }

    ResizeArrayIter<T> operator++(void) {
      currentIndex++;
      return (*this);
    }

    ResizeArrayIter<T> operator++(int) {
      ResizeArrayIter<T> tmp(*this);
      currentIndex++;
      return (tmp);
    }

    T& operator* (void) const {
      return array->rep.array[currentIndex];
    }
};

#endif
