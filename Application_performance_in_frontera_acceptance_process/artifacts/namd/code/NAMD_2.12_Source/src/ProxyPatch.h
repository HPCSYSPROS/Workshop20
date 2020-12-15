/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PROXYPATCH_H
#define PROXYPATCH_H

#include "Patch.h"

class ProxyDataMsg;
//class ProxyAllMsg;

#define PROXYMSGNOTBUFFERED 0
#define PROXYDATAMSGBUFFERED 1
#define PROXYALLMSGBUFFERED 2

class ProxyPatch : public Patch
{
  public:

     ProxyPatch(PatchID pd);
     virtual ~ProxyPatch(void);

     void receiveData(ProxyDataMsg*);
     void receiveAll(ProxyDataMsg*);
     //include gbis phase 1 data with std message
     void receiveData(ProxyGBISP2DataMsg*);//receive P1 results; begin P2
     void receiveData(ProxyGBISP3DataMsg*);//receive P2 results; begin P3

     void setSpanningTree(int, int*, int);
     int  getSpanningTreeParent() { return parent; }
     int  getSpanningTreeChild(int *);
     const int *getSpanningTreeChildPtr() { return child; }
     inline int getSpanningTreeNChild(void) { return nChild; }

    #if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
     void setSTNodeChildren(int numNids, int *nids);
     int *getSTNodeChildPtr() { return nodeChildren; }
     int getSTNNodeChild() { return numNodeChild; }
     #endif

     ProxyCombinedResultMsg *depositCombinedResultMsg(ProxyCombinedResultMsg *);
	 ProxyCombinedResultMsg *depositCombinedResultRawMsg(ProxyCombinedResultRawMsg *);

#if CMK_PERSISTENT_COMM
  private:
     PersistentHandle localphs;
     PersistentHandle *treephs;
     int               ntreephs;
  public:
     PersistentHandle *getSpanningTreePhs(int &n) { n = ntreephs; return treephs; }
#endif
  protected:

     virtual void boxClosed(int);

  private:

     void sendResults(void);

     //"proxyMsgBufferStatus" indicates whether there's a ProxyDataMsg buffered
     // and waiting to be processed, while "curProxyMsg" points to 
     // the actual msg. This msg will be freed at the next step. --Chao Mei  
     int proxyMsgBufferStatus;
     ProxyDataMsg* curProxyMsg;
     ProxyDataMsg* prevProxyMsg;

     // for spanning tree
     ProxyCombinedResultMsg *msgCBuffer;
     int parent;
//#ifdef NODEAWARE_PROXY_SPANNINGTREE
     /* Moved to Patch.h */
     //int *children;
     //int numChild;
//#else
     /* Moved to Patch.h */
     //int *child; // spanning tree for recvResults()
     //int nChild;
//#endif
     int nWait;
     
#if defined(NODEAWARE_PROXY_SPANNINGTREE) && defined(USE_NODEPATCHMGR)
    CmiNodeLock depositLock;
#endif     
};


#endif

