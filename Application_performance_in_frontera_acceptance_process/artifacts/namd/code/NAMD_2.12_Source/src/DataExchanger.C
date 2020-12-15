/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/
#include <stdio.h>
#include "DataExchanger.h"
#include "ProcessorPrivate.h"
#include "common.h"
#include "Node.h"
#include "CollectionMaster.h"
#include "Output.h"
#include "ScriptTcl.h"
#include "qd.h"

#if CMK_HAS_PARTITION
#ifdef CmiMyPartitionSize
extern "C" {
  void setDefaultPartitionParams() {
    // if(!CmiMyNodeGlobal()) printf("NAMD setDefaultPartitionParams called\n");
    CmiSetPartitionScheme(3);  // recursive bisection
  }
}
#endif
#endif

static int recvRedCalledEarly;

CpvDeclare(int, breakScheduler);
CpvDeclare(int, inEval);

//functions to receive and invoke chare's entry methods
extern "C" {
  void packSend(int dst, int dstPart, const char *data, int size, int handler, int code) {
    int msgsize = sizeof(DataMessage) + size;
    DataMessage *dmsg = (DataMessage *)CmiAlloc(msgsize);
    dmsg->setMessage(data,CkMyPe(),CmiMyPartition(),size,handler,code);
#if CMK_HAS_PARTITION
    CmiInterSyncSendAndFree(dst,dstPart,msgsize,(char*)dmsg);
#else
    CmiSyncSendAndFree(dst,msgsize,(char*)dmsg);
#endif
  }

  void sendReplicaDcdInit(int dstPart, ReplicaDcdInitMsg *msg, int msgsize) {
    CmiSetHandler(msg->core, CkpvAccess(recv_replica_dcd_init_idx));
#if CMK_HAS_PARTITION
    CmiInterSyncSendAndFree(0,dstPart,msgsize,(char*)msg);
#else
    CmiSyncSendAndFree(0,msgsize,(char*)msg);
#endif
    CollectionMaster::Object()->blockPositions();  // ensure ordering
    CkpvAccess(_qd)->create();  // ensure completion
  }

  void sendReplicaDcdData(int dstPart, ReplicaDcdDataMsg *msg, int msgsize) {
    CmiSetHandler(msg->core, CkpvAccess(recv_replica_dcd_data_idx));
#if CMK_HAS_PARTITION
    CmiInterSyncSendAndFree(0,dstPart,msgsize,(char*)msg);
#else
    CmiSyncSendAndFree(0,msgsize,(char*)msg);
#endif
    CollectionMaster::Object()->blockPositions();  // ensure ordering
    CkpvAccess(_qd)->create();  // ensure completion
  }

  void sendReplicaDcdAck(int dstPart, ReplicaDcdAckMsg *msg) {
    CmiSetHandler(msg->core, CkpvAccess(recv_replica_dcd_ack_idx));
    int msgsize = sizeof(ReplicaDcdAckMsg);
#if CMK_HAS_PARTITION
    CmiInterSyncSendAndFree(0,dstPart,msgsize,(char*)msg);
#else
    CmiSyncSendAndFree(0,msgsize,(char*)msg);
#endif
  }

  void recvReplicaDcdInit(ReplicaDcdInitMsg *msg) {
    Node::Object()->output->recvReplicaDcdInit(msg);
    CmiFree(msg);
  }

  void recvReplicaDcdData(ReplicaDcdDataMsg *msg) {
    Node::Object()->output->recvReplicaDcdData(msg);
    CmiFree(msg);
  }

  void recvReplicaDcdAck(ReplicaDcdAckMsg *msg) {
    CmiFree(msg);
    CollectionMaster::Object()->unblockPositions();
    CkpvAccess(_qd)->process();
  }

  void recvData(DataMessage *dmsg) {
    Pointer msg(dmsg);
    if ( CkpvAccess(BOCclass_group).dataExchanger.isZero() ) NAMD_bug("BOCgroup::dataExchanger is zero in recvData!");
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_data(msg);
  }

  void recvAck(DataMessage *dmsg) {
    if ( CkpvAccess(BOCclass_group).dataExchanger.isZero() ) NAMD_bug("BOCgroup::dataExchanger is zero in recvAck!");
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_ack();
    CmiFree(dmsg);
  }

  void recvBcast(DataMessage *dmsg) {
    if ( CkpvAccess(BOCclass_group).dataExchanger.isZero() ) NAMD_bug("BOCgroup::dataExchanger is zero in recvBcast!");
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_bcast();
    CmiFree(dmsg);
  }

  void recvRed(DataMessage *dmsg) {
    if ( CkpvAccess(BOCclass_group).dataExchanger.isZero() ) {
      ++recvRedCalledEarly;
      CmiFree(dmsg);
      return;
    }
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_red();
    CmiFree(dmsg);
  }

  void recvEvalCommand(DataMessage *dmsg) {
    Pointer msg(dmsg);
    if ( CkpvAccess(BOCclass_group).dataExchanger.isZero() ) NAMD_bug("BOCgroup::dataExchanger is zero in recvEvalCommand!");
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_eval_command(msg);
  }

  void recvEvalResult(DataMessage *dmsg) {
    Pointer msg(dmsg);
    if ( CkpvAccess(BOCclass_group).dataExchanger.isZero() ) NAMD_bug("BOCgroup::dataExchanger is zero in recvEvalResult!");
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_eval_result(msg);
  }

  void replica_send(char *sndbuf, int sendcount, int destPart, int destPE) {
    if ( CpvAccess(inEval) ) {
      packSend(destPE,destPart,sndbuf,sendcount,CkpvAccess(recv_data_idx),1);
      return;
    }
    Pointer sendPointer(sndbuf);
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].send(sendPointer,sendcount,destPart,destPE); 
    CpvAccess(breakScheduler) = 0;
    while(!CpvAccess(breakScheduler)) CsdSchedulePoll();
  }

  void replica_recv(DataMessage **precvMsg, int srcPart, int srcPE) {
    Pointer recvPointer((char *) precvMsg);
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv(recvPointer,srcPart,srcPE);
    CpvAccess(breakScheduler) = 0;
    while(!CpvAccess(breakScheduler)) CsdSchedulePoll();
  }

  void replica_sendRecv(char *sndbuf, int sendcount, int destPart, int destPE, DataMessage **precvMsg, int srcPart, int srcPE)  {
    Pointer sendPointer(sndbuf);
    Pointer recvPointer((char *) precvMsg);
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].sendRecv(sendPointer,sendcount,destPart,destPE,recvPointer,srcPart,srcPE);
    CpvAccess(breakScheduler) = 0;
    while(!CpvAccess(breakScheduler)) CsdSchedulePoll();
  }

  void replica_eval(char *cmdbuf, int targPart, int targPE, DataMessage **precvMsg) {
    Pointer sendPointer(cmdbuf);
    Pointer recvPointer((char *) precvMsg);
    int sendcount = strlen(cmdbuf) + 1;
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].eval(sendPointer,sendcount,targPart,targPE,recvPointer);
    CpvAccess(breakScheduler) = 0;
    while(!CpvAccess(breakScheduler)) CsdSchedulePoll();
  }

  void replica_barrier() {
    for ( ; recvRedCalledEarly > 0; --recvRedCalledEarly ) CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].recv_red();
    CPROXY_DE(CkpvAccess(BOCclass_group).dataExchanger)[CkMyPe()].barrier();
    CpvAccess(breakScheduler) = 0;
    while(!CpvAccess(breakScheduler)) CsdSchedulePoll();
  }

  void replica_bcast(char *buf, int count, int root) {
    if ( root != 0 ) NAMD_bug("replica_bcast root must be zero");
    int rank = CmiMyPartition();
    int size = CmiNumPartitions();
    int i;
    for ( i=1; i<size; i*=2 );
    for ( i/=2; i>0; i/=2 ) {
      if ( rank & (i - 1) ) continue;
      if ( rank & i ) {
        int src = rank - i;
        //CkPrintf("rank %d recv from %d\n", rank, src);
        DataMessage *recvMsg = NULL;
        replica_recv(&recvMsg, src, CkMyPe());
        if ( recvMsg == NULL ) NAMD_bug("recvMsg == NULL in replica_bcast");
        if ( recvMsg->size != count ) NAMD_bug("size != count in replica_bcast");
        memcpy(buf, recvMsg->data, count);
        CmiFree(recvMsg);
      } else {
        int dst = rank + i;
        if ( dst < size ) {
          //CkPrintf("rank %d send to %d\n", rank, dst);
          replica_send(buf, count, dst, CkMyPe());
        }
      }
    }
  }

  void replica_min_double(double *dat, int count) {
    int rank = CmiMyPartition();
    int size = CmiNumPartitions();
    for ( int i=1; i<size; i*=2 ) {
      if ( rank & i ) {
        int dst = rank - i;
        //CkPrintf("rank %d send to %d\n", rank, dst);
        replica_send((char*)dat, count * sizeof(double), dst, CkMyPe());
      } else {
        int src = rank + i;
        if ( src < size ) {
          //CkPrintf("rank %d recv from %d\n", rank, src);
          DataMessage *recvMsg = NULL;
          replica_recv(&recvMsg, src, CkMyPe());
          if ( recvMsg == NULL ) NAMD_bug("recvMsg == NULL in replica_bcast");
          if ( recvMsg->size != count * sizeof(double) ) NAMD_bug("size != count in replica_min_double");
          double *rdat = new double[count];
          memcpy(rdat, recvMsg->data, count * sizeof(double));
          CmiFree(recvMsg);
          for ( int j=0; j<count; ++j ) {
            if ( rdat[j] < dat[j] ) dat[j] = rdat[j];
          }
          delete [] rdat;
        }
      }
      if ( rank & (2 * i - 1) ) break;
    }
    replica_bcast((char*)dat, count * sizeof(double), 0);
  }
} //endof extern C

#if CMK_IMMEDIATE_MSG && CMK_SMP && ! ( CMK_MULTICORE || CMK_SMP_NO_COMMTHD )
extern "C" void CmiPushImmediateMsg(void *msg);

class SleepCommthdMsg {
  public:
  char core[CmiMsgHeaderSizeBytes];
};

void recvSleepCommthdMsg(SleepCommthdMsg *msg) {
  if ( CkMyRank() != CkMyNodeSize() ) NAMD_bug("recvSleepCommthdMsg called on PE instead of communication thread");
  usleep(1000);
  CmiDelayImmediate();  // re-enqueue for next cycle
}
#endif

void initializeReplicaConverseHandlers() {
  CkpvInitialize(int, recv_data_idx);
  CkpvInitialize(int, recv_ack_idx);
  CkpvInitialize(int, recv_bcast_idx);
  CkpvInitialize(int, recv_red_idx);
  CkpvInitialize(int, recv_eval_command_idx);
  CkpvInitialize(int, recv_eval_result_idx);
  CkpvInitialize(int, recv_replica_dcd_init_idx);
  CkpvInitialize(int, recv_replica_dcd_data_idx);
  CkpvInitialize(int, recv_replica_dcd_ack_idx);

  CkpvAccess(recv_data_idx) = CmiRegisterHandler((CmiHandler)recvData);                   
  CkpvAccess(recv_ack_idx) = CmiRegisterHandler((CmiHandler)recvAck);                     
  CkpvAccess(recv_red_idx) = CmiRegisterHandler((CmiHandler)recvRed);                     
  CkpvAccess(recv_bcast_idx) = CmiRegisterHandler((CmiHandler)recvBcast);                 
  CkpvAccess(recv_eval_command_idx) = CmiRegisterHandler((CmiHandler)recvEvalCommand);    
  CkpvAccess(recv_eval_result_idx) = CmiRegisterHandler((CmiHandler)recvEvalResult);      
  CkpvAccess(recv_replica_dcd_init_idx) = CmiRegisterHandler((CmiHandler)recvReplicaDcdInit);
  CkpvAccess(recv_replica_dcd_data_idx) = CmiRegisterHandler((CmiHandler)recvReplicaDcdData);
  CkpvAccess(recv_replica_dcd_ack_idx) = CmiRegisterHandler((CmiHandler)recvReplicaDcdAck);

#if CMK_IMMEDIATE_MSG && CMK_SMP && ! ( CMK_MULTICORE || CMK_SMP_NO_COMMTHD )
  int sleep_commthd_idx = CmiRegisterHandler((CmiHandler)recvSleepCommthdMsg);
  if ( CkMyPe() == 0 && CkNumNodes() == 1 ) {
    CkPrintf("Charm++ communication thread will sleep due to single-process run.\n");
    SleepCommthdMsg *msg = (SleepCommthdMsg *)malloc(sizeof(SleepCommthdMsg));
    CmiSetHandler(msg, sleep_commthd_idx);
    CmiBecomeImmediate(msg);
    CmiPushImmediateMsg(msg);
  }
#endif
}

//======================================================================
// Public functions
//----------------------------------------------------------------------
DataExchanger::DataExchanger()
{
  CpvInitialize(int, breakScheduler);
  CpvAccess(breakScheduler) = 1;
  CpvInitialize(int, inEval);
  CpvAccess(inEval) = 0;
  if(CmiMyPartition() == 0) 
    parent = -1;
  else 
    parent = (CmiMyPartition()+1)/TREE_WIDTH - 1;
  firstChild = (CmiMyPartition()+1)*TREE_WIDTH - 1;
  numChildren = CmiNumPartitions() - firstChild;
  if(numChildren > TREE_WIDTH)
    numChildren = TREE_WIDTH;
  
  recv_data_idx = CkpvAccess(recv_data_idx);
  recv_ack_idx = CkpvAccess(recv_ack_idx);
  recv_red_idx = CkpvAccess(recv_red_idx);
  recv_bcast_idx = CkpvAccess(recv_bcast_idx);
  recv_eval_command_idx = CkpvAccess(recv_eval_command_idx);
  recv_eval_result_idx = CkpvAccess(recv_eval_result_idx);
}

//----------------------------------------------------------------------
DataExchanger::~DataExchanger(void)
{ }


#include "DataExchanger.def.h"
