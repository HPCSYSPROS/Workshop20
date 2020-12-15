/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Helper class for PatchMgr to manager HomePatch(es)
   It is object contained in the container HomePatchList
*/

#ifndef HOMEPATCHLIST_H
#define HOMEPATCHLIST_H
#include "NamdTypes.h"
#include "SortedArray.h"
#include "ResizeArrayIter.h"

class HomePatch;

class HomePatchElem {
public:
  PatchID   pid;
  HomePatch *patch;

  int operator<(HomePatchElem e) { return (pid < e.pid); }
  int operator==(HomePatchElem e) { return (pid == e.pid); }

  HomePatchElem(PatchID id=-1, HomePatch *p=NULL) : pid(id), patch(p) {};
  ~HomePatchElem() { };
  HomePatchElem& operator=(const HomePatchElem &e) { 
    pid = e.pid; 
    patch = e.patch;  // Do not delete patch!  This op only used to shuffle
		      // we delete the patch here only when the HomePatch is 
		      // moved off!
    return(*this);
  };
};

typedef SortedArray<HomePatchElem> HomePatchList;
typedef ResizeArrayIter<HomePatchElem> HomePatchListIter;

#endif /* HOMEPATCHLIST_H */

