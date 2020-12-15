/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <string.h>
#include <stdlib.h>
#include "Communicate.h"
#include "MStream.h"
#include "charm++.h"

CkpvStaticDeclare(CmmTable, CsmMessages);
CkpvStaticDeclare(int, CsmAcks);

static void CsmHandler(void *msg)
{
  if ( CmiMyRank() ) NAMD_bug("Communicate CsmHandler on non-rank-zero pe");
  // get start of user message
  int *m = (int *) ((char *)msg+CmiMsgHeaderSizeBytes);
  // sending node  & tag act as tags
  CmmPut(CkpvAccess(CsmMessages), 2, m, msg);
}

static void CsmAckHandler(void *msg)
{
  if ( CmiMyRank() ) NAMD_bug("Communicate CsmAckHandler on non-rank-zero pe");
  CmiFree(msg);
  CkpvAccess(CsmAcks) += 1;
}

Communicate::Communicate(void) 
{
  CkpvInitialize(CmmTable, CsmMessages);
  CsmHandlerIndex = CmiRegisterHandler((CmiHandler) CsmHandler);
  CsmAckHandlerIndex = CmiRegisterHandler((CmiHandler) CsmAckHandler);
  CkpvAccess(CsmMessages) = CmmNew();

  int parent_node = 0;
  int self = CkMyNode();
  int range_begin = 0;
  int range_end = CkNumNodes();
  while ( self != range_begin ) {
    parent_node = range_begin;
    ++range_begin;
    int split = range_begin + ( range_end - range_begin ) / 2;
    if ( self < split ) { range_end = split; }
    else { range_begin = split; }
  }
  int send_near = self + 1;
  int send_far = send_near + ( range_end - send_near ) / 2;

  parent = CkNodeFirst(parent_node);
  nchildren = 0;
  if ( send_far < range_end ) children[nchildren++] = CkNodeFirst(send_far);
  if ( send_near < send_far ) children[nchildren++] = CkNodeFirst(send_near);

  CkpvInitialize(int, CsmAcks);
  CkpvAccess(CsmAcks) = nchildren;

  ackmsg = (char *) CmiAlloc(CmiMsgHeaderSizeBytes);
  CmiSetHandler(ackmsg, CsmAckHandlerIndex);
}


Communicate::~Communicate(void) 
{
  CmiFree(ackmsg);
}

MIStream *Communicate::newInputStream(int PE, int tag)
{
  MIStream *st = new MIStream(this, PE, tag);
  return st;
}

MOStream *Communicate::newOutputStream(int PE, int tag, unsigned int bufSize)
{
  MOStream *st = new MOStream(this, PE, tag, bufSize);
  return st;
}

void *Communicate::getMessage(int PE, int tag)
{
  if ( CmiMyRank() ) NAMD_bug("Communicate::getMessage called on non-rank-zero Pe\n");

  int itag[2], rtag[2];
  void *msg;

  itag[0] = (PE==(-1)) ? (CmmWildCard) : PE;
  itag[1] = (tag==(-1)) ? (CmmWildCard) : tag;
  while((msg=CmmGet(CkpvAccess(CsmMessages),2,itag,rtag))==0) {
    CmiDeliverMsgs(0);
    // CmiDeliverSpecificMsg(CsmHandlerIndex);
  }

  CmiSyncSend(parent, CmiMsgHeaderSizeBytes, ackmsg);

  while ( CkpvAccess(CsmAcks) < nchildren ) {
    CmiDeliverMsgs(0);
    // CmiDeliverSpecificMsg(CsmAckHandlerIndex);
  }
  CkpvAccess(CsmAcks) = 0;

  int size = SIZEFIELD(msg);
  for ( int i = 0; i < nchildren; ++i ) {
    CmiSyncSend(children[i],size,(char*)msg);
  }

  return msg;
}

void Communicate::sendMessage(int PE, void *msg, int size)
{
  if ( CmiMyPe() ) NAMD_bug("Communicate::sendMessage not from Pe 0");

  while ( CkpvAccess(CsmAcks) < nchildren ) {
    CmiDeliverMsgs(0);
    // CmiDeliverSpecificMsg(CsmAckHandlerIndex);
  }
  CkpvAccess(CsmAcks) = 0;

  CmiSetHandler(msg, CsmHandlerIndex);
  switch(PE) {
    case ALL:
      NAMD_bug("Unexpected Communicate::sendMessage(ALL,...)");
      //CmiSyncBroadcastAll(size, (char *)msg);
      break;
    case ALLBUTME:
      //CmiSyncBroadcast(size, (char *)msg);
      for ( int i = 0; i < nchildren; ++i ) {
        CmiSyncSend(children[i],size,(char*)msg);
      }
      break;
    default:
      NAMD_bug("Unexpected Communicate::sendMessage(PEL,...)");
      //CmiSyncSend(PE, size, (char *)msg);
      break;
  }
}
