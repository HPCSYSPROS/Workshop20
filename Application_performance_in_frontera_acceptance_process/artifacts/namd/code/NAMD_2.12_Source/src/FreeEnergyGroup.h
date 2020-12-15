/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

//---------------------------------------------------------------
// AGroup contains a (potentially long) list of integers
// written by David Hurwitz, March to May 1998.
//---------------------------------------------------------------
#if !defined(GROUP_HPP)
  #define GROUP_HPP

const int kGroupNumToStart = 16;    // to start, there's room for this num ints.
const int kGroupMultiplier = 4;     // each time array size is exceeded,
                                    // its size is increased by this many times.

class AGroup {
private:
  int*  m_pInts;      // the list of integers
  int   m_NumInts;    // the number of integers in the list
  int   m_MaxNum;     // the maximum number of integers allowed in
                      // the list without allocating more memory

public:
  AGroup();
  ~AGroup();
  void    Clear();
  void    Add(int AnInt);
  AGroup& operator= (AGroup& Group);
  int     operator[] (int Index);
  int     GetNumInGroup() { return(m_NumInts); }
  void    List(int NumToList=-1);
};

#endif

