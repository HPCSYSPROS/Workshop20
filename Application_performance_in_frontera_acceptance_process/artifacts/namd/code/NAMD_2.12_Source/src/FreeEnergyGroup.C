/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz, March to May 1998.

#include <memory.h>
#include <iostream>
#include <iomanip>
using namespace std;
#include "charm++.h"
#include "FreeEnergyAssert.h"
#include "FreeEnergyGroup.h"


AGroup::AGroup() {
//------------------------------------------------------------------------
// make room for some ints.
//------------------------------------------------------------------------
  m_NumInts = 0;
  m_pInts = new int[kGroupNumToStart];
  m_MaxNum = kGroupNumToStart;
}


AGroup::~AGroup() {
//------------------------------------------------------------------------
// return borrowed memory to the free store.
//------------------------------------------------------------------------
  delete []m_pInts;
}


void AGroup::Clear() {
//------------------------------------------------------------------------
// leave memory allocation alone.
//------------------------------------------------------------------------
  m_NumInts = 0;
}


void AGroup::Add(int AnInt) {
//------------------------------------------------------------------------
// add an int to the list.  if there's not enough room, make room.
//------------------------------------------------------------------------
  int*  pInts;

  // if there's no room for a new int
  if (m_NumInts == m_MaxNum) {
    // create an array with more space
    m_MaxNum *= kGroupMultiplier;
    pInts = new int[m_MaxNum];
    // fast copy from the full array to the new one (memcpy(dest, src, bytes))
    memcpy(pInts, m_pInts, sizeof(int)*m_NumInts);
    // return the space used for the full array
    delete []m_pInts;
    // point to the bigger array
    m_pInts = pInts;
  }

  // add the int to the int array
  m_pInts[m_NumInts] = AnInt;
  m_NumInts++;
}


int AGroup::operator[] (int Index) {
//------------------------------------------------------------------------
// return an int from this group of ints.
// note, by returning int, rather than int&, this function only allows
// retrieval of an int, not setting an int.  to set, use add.
//------------------------------------------------------------------------
  ASSERT((Index>=0) && (Index<m_NumInts));
  return(m_pInts[Index]);
}


AGroup& AGroup::operator= (AGroup& Group) {
//------------------------------------------------------------------------
// make this object identical to the passed one.
//------------------------------------------------------------------------
  // if there's not enough room here for Group's array
  if (m_MaxNum < Group.m_MaxNum) {
    // free this array and allocate space for the bigger one
    delete []m_pInts;
    m_MaxNum =  Group.m_MaxNum;
    m_pInts = new int[m_MaxNum];
  }
  // fast copy from the passed array to this one (memcpy(dest, src, bytes))
  m_NumInts = Group.m_NumInts;
  memcpy(m_pInts, Group.m_pInts, sizeof(int)*m_NumInts);
  return(*this);
}

void AGroup::List(int NumToList) {
//------------------------------------------------------------------------
// list NumToList integers in the group to standard out.
// if NumToList is negative, list them all.
//------------------------------------------------------------------------
  int  i;

  if (NumToList < 0) {
    NumToList = m_NumInts;
  }
  for (i=0; i<NumToList; i++) {
//    cout << setw(10) << i << "   " << setw(10) << (*this)[i] << std::endl;
    CkPrintf("%10d   %10d\n", i, (*this)[i]);
  }
}
