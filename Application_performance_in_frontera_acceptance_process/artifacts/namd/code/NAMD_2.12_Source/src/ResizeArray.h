/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   ResizeArray template
   Uses simple contingous array allocation in a hidden manner
   so that array object can have items added without limit
   Suffers from memory fragmentation during resizing

   Copy construction and assignment are forbidden to avoid
   confusion with old reference-counted implementation.
   Use copy() and swap() to explicitly transfer storage.
*/

#ifndef RESIZEARRAY_H
#define RESIZEARRAY_H

#include "ResizeArrayRaw.h"

// Need this juju to use templated friend below
template <class Type> class ResizeArrayIter;
template <class Type> class ResizeArrayRaw;

template <class Elem> class ResizeArray {
  friend class ResizeArrayIter<Elem>;

  protected:
    ResizeArrayRaw<Elem> rep;

  public:
    // STL style iterators
    typedef Elem* iterator;
    iterator begin(void) { return rep.array; }
    iterator end(void) { return rep.array + rep.arraySize; }
    typedef const Elem* const_iterator;
    const_iterator const_begin(void) const { return rep.array; }
    const_iterator const_end(void) const { return rep.array + rep.arraySize; }

    // Various Constructors
    ResizeArray(void) { }

    // Constructor make ResizeArray of predefined size
    ResizeArray(int s) {
      rep.resize(s);
    }

private:
    inline ResizeArray(ResizeArray<Elem> &ra);

    inline ResizeArray(const ResizeArray<Elem>* ra);

    inline ResizeArray<Elem> & operator= (ResizeArray<Elem> &ra);

public:
    // Make copy of ResizeArrayRaw (for use in messages)
    void copy(ResizeArray<Elem> &ra) {
      rep.copy(ra.rep);
    }

    // Swap ResizeArrayRaw (for use in messages to avoid copy)
    void swap(ResizeArray<Elem> &ra) {
      // uses synthesized copy constructor and assignment operator
      ResizeArrayRaw<Elem> tmp = rep;
      rep = ra.rep;
      ra.rep = tmp;
    }

    // does some other ResizeArray have a handle to our data
    bool shared() const {
      return 0;
    }

    // Constructor to take-in pre-existing array
    ResizeArray(Elem * * array, int arraySize, int allocSize=0) :
      rep(array, arraySize, allocSize) { }

    ~ResizeArray(void) { rep.free(); }

    // If array is expanded - new elements are default constructed
    // if array is reduced, removed elements have ~Elem() run
    void resize(int i) { rep.resize(i); }

    // destruct elements, free storage, set size to 0
    void clear() { rep.clear(); }

    // Set all elements to a given value (like 0).
    void setall(const Elem &elem) {
      iterator i = begin();
      iterator e = end();
      for ( ; i != e; ++i ) *i = elem;
    }
  
    // Add element to end of array
    int add (const Elem &elem) {
      int end=rep.size();
      rep.ins(elem, end);
      return(end);
    }
  
    // delete num elements from current index
    int del(int index, int num=1) {
      return(rep.del(index,num));
    }

    // insert element at index
    int insert (const Elem& elem, int index) {
      rep.ins(elem,index);
      return (index);
    }

    // array member access (can be lvalue) that grows array.
    inline Elem & item(int i) {
      i = ( i < 0 ? 0 : i );
      if ((i+1) > size())
          resize(i+1);
      return rep.array[i];
    }

    // array member access (can be lvalue) no checks.
    inline Elem & operator[](int index) { return rep.array[index]; }
    inline const Elem & operator[](int index) const { return rep.array[index]; }

    // returns size of ResizeArray
    inline int size(void) const { return rep.size(); }

    // DMK - MIC Support - Allow us to see the buffer's size, not just how much of it is used
    #if NAMD_MIC != 0
      inline int bufSize(void) const { return rep.bufSize(); }
    #endif

    // reduce storage size
    // void reduce(void) { rep.reduce(); }

    inline int find(const Elem &e) const { return rep.find(e); }

};

#endif
