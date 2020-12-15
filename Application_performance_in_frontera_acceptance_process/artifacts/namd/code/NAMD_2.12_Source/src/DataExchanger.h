/**
 ***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
 ***  The Board of Trustees of the University of Illinois.
 ***  All rights reserved.
 **/

#ifndef _DATAEXCHANGER_H
#define _DATAEXCHANGER_H

#include "charm++.h"
#include "main.h"
#include "DataExchanger.decl.h"

#define CPROXY_DE(x) ((CProxy_DataExchanger)(x))
CpvExtern(int, breakScheduler);
CpvExtern(int, inEval);

class DataMessage {
  public:
  char core[CmiMsgHeaderSizeBytes];
  int src, srcPart;
  int size, code;
  char data[1];

  void setMessage(const char *_data, int _src, int _srcPart, int _size, int _handler, int _code) {
    src = _src; srcPart = _srcPart;
    size = _size;
    code = _code;
    memcpy(data,_data,size);
    CmiSetHandler(core,_handler);
  }
};

class ReplicaDcdInitMsg {
  public:
  char core[CmiMsgHeaderSizeBytes];
  int srcPart;
  int dcdIndex;
  char data[1];
};

class ReplicaDcdDataMsg {
  public:
  char core[CmiMsgHeaderSizeBytes];
  int srcPart;
  int dcdIndex;
  int NSAVC, NPRIV, NSTEP, NFILE, with_unitcell, numAtoms;
  double DELTA;
  double unitcell[6];
  char data[1];
};

class ReplicaDcdAckMsg {
  public:
  char core[CmiMsgHeaderSizeBytes];
};

class DataExchanger : public CBase_DataExchanger
{
  public:
    DataExchanger_SDAG_CODE;
    DataExchanger();
    ~DataExchanger(void);
  
    int loop, recvred, sendbcast;
    enum{ TREE_WIDTH=2};
    int numChildren, firstChild, parent;

    //message handlers
    int recv_data_idx;
    int recv_ack_idx;
    int recv_bcast_idx;
    int recv_red_idx;
    int recv_eval_command_idx;
    int recv_eval_result_idx;
};

extern "C" {
void packSend(int dest, int partition, const char *data, int size, int handler, int code=0);
void sendReplicaDcdInit(int dstPart, ReplicaDcdInitMsg *msg, int msgsize);
void sendReplicaDcdData(int dstPart, ReplicaDcdDataMsg *msg, int msgsize);
void sendReplicaDcdAck(int dstPart, ReplicaDcdAckMsg *msg);
void recvReplicaDcdInit(ReplicaDcdInitMsg *msg);
void recvReplicaDcdData(ReplicaDcdDataMsg *msg);
void recvReplicaDcdAck(ReplicaDcdAckMsg *msg);
void recvData(DataMessage *dmsg); 
void recvAck(DataMessage *dmsg); 

void replica_send(char *sndbuf, int sendcount, int destPart, int destPE);
void replica_sendRecv(char *sndbuf, int sendcount, int destPart, int destPE, DataMessage **precvMsg, int srcPart, int srcPE);
void replica_recv(DataMessage **precvMsg, int srcPart, int srcPE);
void replica_barrier();

void replica_bcast(char *buf, int count, int root=0);
void replica_min_double(double *dat, int count);

void replica_eval(char *cmdbuf, int targPart, int targPE, DataMessage **precvMsg);
}
#endif
