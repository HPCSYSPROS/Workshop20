/*
 *  MGridforceParams.h
 *  
 *
 *  Created by Robert Brunner on 12/5/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef MGRIDFORCEPARAMS_H
#define MGRIDFORCEPARAMS_H

#include "strlib.h"
#include "common.h"
#include "Vector.h"
#include "InfoStream.h"
#include "MStream.h"

#define MGRIDFORCEPARAMS_DEFAULTKEY "BaseGridForceParams"

class MGridforceParams {
public:
  MGridforceParams() {
    gridforceKey = 0;
    gridforceVfile = 0;
    gridforceScale = Vector(0);
    gridforceFile = 0;
    gridforceCol = 0;
    gridforceQcol = 0;
    gridforceVOffset = Vector(0);
    gridforceCont[0] = gridforceCont[1] = gridforceCont[2] = FALSE;
    gridforceVolts = FALSE;
    gridforceLite = FALSE;
    gridforceCheckSize = TRUE;
  }
  
  char *gridforceKey;
  char *gridforceVfile;
  zVector gridforceScale;
  char *gridforceFile;
  char *gridforceCol;
  char *gridforceQcol;
  zVector gridforceVOffset;
  Bool gridforceCont[3];
  Bool gridforceVolts;
  Bool gridforceLite;
  Bool gridforceCheckSize;
  MGridforceParams *next;
};

class MGridforceParamsList {
public:
  MGridforceParamsList() {
    clear();
  }
  
  ~MGridforceParamsList() 
  {
    MGFElem* cur;
    while (head != NULL) {
      cur = head;
      head = cur->nxt;
      delete cur;
    }
    clear();
  }
  
  // The SimParameters bit copy overwrites these values with illegal pointers,
  // So thise throws away the garbage and lets everything be reinitialized
  // from scratch
  void clear() {
    head = tail = NULL;
    n_elements = 0;
  }
  
  MGridforceParams* find_key(const char* key);  
  int index_for_key(const char* key);
  MGridforceParams* at_index(int idx);
  MGridforceParams* add(const char* key);
  
  MGridforceParams *get_first() {
    if (head == NULL) {
      return NULL;
    } else return &(head->elem);
  }
  
  void pack_data(MOStream *msg);  
  void unpack_data(MIStream *msg);
  
  // convert from a string to Bool; returns 1(TRUE) 0(FALSE) or -1(if unknown)
  static int atoBool(const char *s)
  {
    if (!strcasecmp(s, "on")) return 1;
    if (!strcasecmp(s, "off")) return 0;
    if (!strcasecmp(s, "true")) return 1;
    if (!strcasecmp(s, "false")) return 0;
    if (!strcasecmp(s, "yes")) return 1;
    if (!strcasecmp(s, "no")) return 0;
    if (!strcasecmp(s, "1")) return 1;
    if (!strcasecmp(s, "0")) return 0;
    return -1;
  }
  
private:
  class MGFElem {
  public:
    MGridforceParams elem;
    MGFElem* nxt;
  };
  MGFElem* head;
  MGFElem* tail;
  int n_elements;
};

#endif
