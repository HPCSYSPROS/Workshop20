/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef ARRAY_H
#define ARRAY_H

template <class Elem, int Size> class Array {

  public:
    Elem data[Size];

    // constructor
    Array(void) { ; }

    // destructor
    ~Array(void) { ; }

    // copy constructor
    Array(const Array<Elem,Size> &a2) {
      for ( int i = 0; i < Size; ++i ) { data[i] = a2.data[i]; }
    }

    // assignment operator
    Array<Elem,Size> & operator= (const Array<Elem,Size> &a2) {
      for ( int i = 0; i < Size; ++i ) { data[i] = a2.data[i]; }
      return (*this);
    }

    // set all elements to a given value (like 0)
    Array(const Elem &v) {
      for ( int i = 0; i < Size; ++i ) { data[i] = v; }
    }
    Array<Elem,Size> & operator= (const Elem &v) {
      for ( int i = 0; i < Size; ++i ) { data[i] = v; }
      return (*this);
    }

    // element access
    inline Elem & operator[](int index) { return data[index]; }
    inline const Elem & operator[](int index) const { return data[index]; }

    // STL-style interface
    typedef Elem* iterator;
    iterator begin(void) { return &data[0]; }
    iterator end(void) { return &data[0] + Size; }
    typedef const Elem* const_iterator;
    const_iterator const_begin(void) const { return &data[0]; }
    const_iterator const_end(void) const { return &data[0] + Size; }

};

#endif

