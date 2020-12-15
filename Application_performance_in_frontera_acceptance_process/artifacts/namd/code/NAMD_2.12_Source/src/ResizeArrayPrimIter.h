/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef RAPITER_H
#define RAPITER_H

#include "ResizeArray.h"

// Don't use for speed - use iter if later we will probably want
// to use a better container class for better space or algorithm behavior

template <class T> class ResizeArrayPrimIter {
  private:

    ResizeArray<T> *array;
    int currentIndex;

  public:
    ResizeArrayPrimIter(void) {
      array = NULL;
      currentIndex = 0;
    }

    ResizeArrayPrimIter(ResizeArray<T>& ra) {
      array = &ra;
      currentIndex = 0;
    }

    ResizeArrayPrimIter(const ResizeArrayPrimIter<T>& iter) {
      array = iter.array;
      currentIndex = iter.currentIndex;
    }

    ResizeArrayPrimIter<T>& operator= (const ResizeArrayPrimIter<T>& iter) {
      array = iter.array;
      currentIndex = iter.currentIndex;
      return (*this);
    }
  
    ~ResizeArrayPrimIter(void) {}
  
    ResizeArrayPrimIter<T> begin(void) const {
      ResizeArrayPrimIter<T> iter;
      iter.array = array;
      iter.currentIndex = 0;
      return(iter);
    }

    ResizeArrayPrimIter<T> end(void) const {
      ResizeArrayPrimIter<T> iter;
      iter.array = array;
      iter.currentIndex = array->size();
      return(iter);
    }
      
    int operator!= (const ResizeArrayPrimIter<T> &iter) const {
      return (iter.currentIndex != currentIndex || iter.array != array);
    }

    int operator== (const ResizeArrayPrimIter<T> &iter) const {
      return (!operator!=(iter));
    }
  
    ResizeArrayPrimIter<T> operator++(void) {
      currentIndex++;
      return (*this);
    }

    ResizeArrayPrimIter<T> operator++(int) {
      ResizeArrayPrimIter<T> tmp(*this);
      currentIndex++;
      return (tmp);
    }
  
    T& operator* (void) const {
      return array->operator[](currentIndex);
    }
};

#endif
