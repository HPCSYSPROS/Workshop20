/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#if !defined(ENUMS_HPP)
  #define ENUMS_HPP

enum Bool_t  {kFalse=0, kAlsoTrue=1, kTrue=!kFalse};
enum Bound_t {kLower, kUpper, kUnknownBound};
enum feptask_t  {kUp, kDown, kStop, kGrow, kFade, kNoGrow,
              kStepUp, kStepDown, kStepGrow, kStepFade, kUnknownTask};
enum Error_t {kNoProblem, kParsingProblem, kNAMD_Problem};

const double kAlmostOne =       0.9999999999;
const double kAlmostMinusOne = -0.9999999999;
const double kALittle =         0.0000000001;
const double kPi =              3.1415926536;

#endif
