/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include <string.h>
#include "Communicate.h"
#include "MStream.h"
#include "converse.h"

#define MIN_DEBUG_LEVEL 2
//#define DEBUGM
#include "Debug.h"

struct StreamMessage {
  char header[CmiMsgHeaderSizeBytes];
  int PE;
  int tag;
  size_t len; // sizeof the data
  unsigned int index; // index of packet in stream
  unsigned int checksum;
  StreamMessage *next; // for linked list of early packets
  char data[1];
};

MIStream::MIStream(Communicate *c, int p, int t)
{
  cobj = c;
  PE = p;
  tag = t;
  msg = (StreamMessage *) 0;
  early = (StreamMessage *) 0;
  currentIndex = 0;
  checksum = 0;
}

MIStream::~MIStream()
{
  if(msg!=0)
    CmiFree(msg);
}

MOStream::MOStream(Communicate *c, int p, int t, size_t size)
{
  cobj = c;
  PE = p;
  tag = t;
  bufLen = size;
  msgBuf = (StreamMessage *)CmiAlloc(sizeof(StreamMessage)+size);
  msgBuf->PE = CmiMyPe();
  msgBuf->tag = tag;
  msgBuf->len = 0;
  msgBuf->index = 0;
  msgBuf->next = (StreamMessage *)0;
  msgBuf->checksum = 0;
}

MOStream::~MOStream()
{
  if(msgBuf != 0)
    CmiFree(msgBuf);
}

static int checkSum(StreamMessage *msg)
{
  int checksum = 0;
  for ( size_t i=0; i < msg->len; i++ ) {
    checksum += (unsigned char) msg->data[i];
  }
  if ( checksum != msg->checksum ) {
    DebugM(5,"Error on " << msg->tag << ":" << msg->index <<
          " of length " << msg->len <<
          " with checksum " << ((int)checksum) <<
          " vs " << ((int)(msg->checksum)) <<"\n");
    NAMD_bug("MStream checksums do not agree!");
  }
  return 1;
}

MIStream *MIStream::Get(char *buf, size_t len)
{
  while(len) {
    if(msg==0) {
      if ( early && (early->index < currentIndex) ) {
          DebugM(3,"Duplicate message " << early->index <<
            " from Pe(" << msg->PE << ")" <<
            " on stack while waiting for " << currentIndex <<
            " from Pe(" << PE << ").\n");
          NAMD_bug("MIStream::Get - duplicate message on stack!");
      }
      if ( early && (early->index == currentIndex) ) {
        DebugM(2,"Popping message " << currentIndex << " from stack.\n");
        msg = early;
        early = early->next;
        msg->next = (StreamMessage *)0;
      } else {
        DebugM(1,"Receiving message.\n");
        msg = (StreamMessage *) cobj->getMessage(PE, tag);
        checkSum(msg);
      }
      while ( msg->index != currentIndex ) {
        if ( msg->index < currentIndex ) {
          DebugM(3,"Duplicate message " << msg->index <<
            " from Pe(" << msg->PE << ")" <<
            " received while waiting for " << currentIndex <<
            " from Pe(" << PE << ").\n");
          NAMD_bug("MIStream::Get - duplicate message received!");
        }
        DebugM(2,"Pushing message " << msg->index << " on stack.\n");
        if ( (! early) || (early->index > msg->index) ) {
          msg->next = early;
          early = msg;
        } else {
          StreamMessage *cur = early;
          while ( cur->next && (cur->next->index < msg->index) ) {
            cur = cur->next;
          }
          msg->next = cur->next;
          cur->next = msg;
        }
        DebugM(1,"Receiving message again.\n");
        msg = (StreamMessage *) cobj->getMessage(PE, tag);
        checkSum(msg);
      } 
      currentPos = 0;
      currentIndex += 1;
    }  // end of if (msg==0)
    if(currentPos+len <= msg->len) {
      memcpy(buf, &(msg->data[currentPos]), len);
      currentPos += len;
      len = 0;
    } else {
      size_t b = msg->len-currentPos;
      memcpy(buf, &(msg->data[currentPos]), b);
      len -= b;
      buf += b;
      currentPos += b;
    }
    if(currentPos == msg->len) {
      CmiFree(msg);
      msg = 0;
    }
  }
  return this;
}

MOStream *MOStream::Put(char *buf, size_t len)
{
  while(len) {
    if(msgBuf->len + len <= bufLen) {
      memcpy(&(msgBuf->data[msgBuf->len]), buf, len);
      msgBuf->len += len;
      len = 0;
    } else {
      size_t b = bufLen - msgBuf->len;
      memcpy(&(msgBuf->data[msgBuf->len]), buf, b);
      msgBuf->len = bufLen;
      if ( msgBuf->index && ! ((msgBuf->index) % 100) ) {
        DebugM(3,"Sending message " << msgBuf->index << ".\n");
      }
      msgBuf->checksum = 0;
      for ( size_t i=0; i < msgBuf->len; i++ ) {
        msgBuf->checksum += (unsigned char) msgBuf->data[i];
      }
      cobj->sendMessage(PE, (void *)msgBuf, bufLen+sizeof(StreamMessage)-1);
      msgBuf->len = 0;
      msgBuf->index += 1;
      len -= b;
      buf += b;
    }
  }
  return this;
}

void MOStream::end(void)
{
  if ( msgBuf->len == 0 ) return; // don't send empty message
  if ( msgBuf->index && ! ((msgBuf->index) % 100) ) {
    DebugM(3,"Sending message " << msgBuf->index << ".\n");
  }
  msgBuf->checksum = 0;
  for ( size_t i=0; i < msgBuf->len; i++ ) {
    msgBuf->checksum += (unsigned char) msgBuf->data[i];
  }
  cobj->sendMessage(PE,(void*)msgBuf,msgBuf->len+sizeof(StreamMessage)-1);
  msgBuf->len = 0;
  msgBuf->index += 1;
}

