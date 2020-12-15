/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "FreeEnergyAssert.h"

void my_assert(const char* Condition, const char* FileName, int LineNumber) {
  iout << std::endl << endi;
  iout << "Assertion: " << "(" << Condition << ")," << " failed" << std::endl << endi;
  iout << "   in: " << FileName << ", " << "line: " << LineNumber << std::endl << endi;
  iout << std::endl << endi;
}
