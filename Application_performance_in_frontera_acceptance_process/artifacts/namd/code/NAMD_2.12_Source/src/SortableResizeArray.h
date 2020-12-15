/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef SRTABARRAY_H
#define SRTABARRAY_H

#include <string.h>
#include "ResizeArray.h"

template <class Elem> class SortableResizeArray : public ResizeArray<Elem> {

  private:

    inline void swap(int offset, int i, int j) {
      char tmp[sizeof(Elem)];
      register Elem *r = (this->rep.array+offset);
      memcpy(tmp, (char *)&r[i], sizeof(Elem));
      memcpy((char *)&(r[i]), (char *)&(r[j]), sizeof(Elem));
      memcpy((char *)&r[j], tmp, sizeof(Elem));
    }


    inline void siftup(int offset, int i, int size) {
      char tmp[sizeof(Elem)];
      register int j;
      register Elem *r = (this->rep.array+offset);
    
      while ((j = 2*i+1) < size) {
        if (j+1 < size) {
          if (r[j] < r[j+1])
            j = j+1;
          }
          if (r[i] < r[j]) {
            memcpy(tmp, (char *)&r[i], sizeof(Elem));
            memcpy((char *)&(r[i]), (char *)&(r[j]), sizeof(Elem));
            memcpy((char *)&r[j], tmp, sizeof(Elem));
            i = j;
          } else {
            break;
          }
      }
    }

  public:

    SortableResizeArray(void) : ResizeArray<Elem>() { init(); }

    SortableResizeArray(int size) : ResizeArray<Elem>(size) { init(); }

    SortableResizeArray(ResizeArray<Elem> &ra) : 
      ResizeArray<Elem>(ra) { init(); }

    SortableResizeArray(SortableResizeArray<Elem> &ra) :
      ResizeArray<Elem>(ra) { init(); }

    SortableResizeArray(const ResizeArray<Elem>* ra) : 
      ResizeArray<Elem>(ra) { init(); }

    SortableResizeArray(const SortableResizeArray<Elem>* ra) : 
      ResizeArray<Elem>(ra) { init(); }

    SortableResizeArray(Elem* * const r, int numElem, int maxElem) : 
      ResizeArray<Elem>(r,numElem,maxElem) { init(); }

    SortableResizeArray<Elem>& operator =(SortableResizeArray<Elem>& sa) {
      ResizeArray<Elem>::operator=(sa);
      return(*this);
    }

    ~SortableResizeArray(void) { }
  
    void init(void) { }
  
    void sort(void) { sort(0, this->rep.size()-1); }

    // Heap Sort - worst case is O(n log n)
    //      bot = bottom element of sort range
    //      top = top element of sort range
    void sort(int bot, int top) {
      int index, size;
      if (top > this->rep.size()) top = this->rep.size();
      size = top - bot + 1;
    
      // Make all sub-heaps
      for ( index = size/2-1; index > 0; index-- )
        siftup(bot, index, size);
    
      // Take top element of overall heap, and put on top
      for ( index = size; index > 1; index-- ) {
        siftup(bot, 0, index);
        swap(bot, 0, index-1);
      }
    }

    inline void uniq(void);

    // Search returns index of position where elem should be inserted.
    // This is equal to the first position equal to elem
    // or the first item greater than elem
    // if elem is larger than any item, it returns
    // the index just beyond the end of the list
    // We stick with the < operator only
    int bsearch(const Elem& elem) const {
      int test;
      int bot = -1;
      int top = this->size();
      if (this->size() == 0) return (-1);
      while (top - bot > 1) {
        if ( this->rep.array[test = (bot+top)/2] < elem )
          bot = test;
        else
          top = test;
      }
      return(top);
    }

};

template <class Elem>
inline void SortableResizeArray<Elem>::uniq(void) {
  if (this->size()) {
    int oldIndex=0;
    int newIndex=0;
    while (++oldIndex < this->size()) {
      if ( ! ( this->rep.array[oldIndex] == this->rep.array[newIndex] ) ) {
        if (++newIndex != oldIndex)
          memcpy((void *)&(this->rep.array[newIndex]),
                 (void *)&(this->rep.array[oldIndex]),
                 sizeof(Elem));
      } else {
        this->rep.array[oldIndex].~Elem();
      }
    }
    this->rep.arraySize = ++newIndex;
  }
}

#endif
