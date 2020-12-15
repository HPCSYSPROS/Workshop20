/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   ReserveArray template
   Uses simple contingous array allocation in a hidden manner
   so that array object can have items added without limit
   Suffers from memory fragmentation during resizing
   Fast access, safe and efficient passing of encapsulated array thru
   function arguments.
*/

#if 0
//  No longer in use and probably not the best solution anyway.

#ifndef RESERVEARRAY_H
#define RESERVEARRAY_H

#define RESERVEARRAY(TYPE,NAME,RSIZE,SIZE) \
  TYPE * NAME; \
  ReserveArray<TYPE,RSIZE> NAME ## _reserve(SIZE,&NAME);

template <class Elem, int reservedSize> class ReserveArray {

  Elem *allocatedStorage;
  Elem reservedStorage[reservedSize];

  public:

  ReserveArray(int size, Elem **userStorage) {
    if ( size > reservedSize ) {
      *userStorage = allocatedStorage = new Elem[size];
    } else if ( size > 0 ) {
      allocatedStorage = 0;
      *userStorage = reservedStorage;
    } else {
      allocatedStorage = 0;
      *userStorage = 0;
    }
  }

  ~ReserveArray() {
    delete [] allocatedStorage;
  }

};

#endif

#endif

