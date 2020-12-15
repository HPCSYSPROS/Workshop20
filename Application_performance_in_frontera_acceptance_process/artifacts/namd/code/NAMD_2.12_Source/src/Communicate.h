/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef COMMUNICATE_H
#define COMMUNICATE_H

class MIStream;
class MOStream;

#define ALL      -1
#define ALLBUTME -2
#define BUFSIZE  4096
#define ANY      -1

class Communicate {

private:
  int CsmHandlerIndex;
  int CsmAckHandlerIndex;
  int parent;
  int nchildren;
  int children[2];
  char *ackmsg;

public:
  Communicate(void);
  ~Communicate();
  MIStream *newInputStream(int pe, int tag);
  MOStream *newOutputStream(int pe, int tag, unsigned int bufsize);
  void *getMessage(int PE, int tag);
  void sendMessage(int PE, void *msg, int size);
};

#include "MStream.h"

#endif
