#ifndef BBOXSET_HH
#define BBOXSET_HH

#include "defs.hh"

#include "bboxset1.hh"
#ifdef CARPET_ENABLE_BBOXSET2
#include "bboxset2.hh"
#endif

#ifdef CARPET_USE_BBOXSET2
using namespace bboxset2;
#else
using namespace bboxset1;
#endif

#endif // #ifndef BBOXSET_HH
