/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "ccsinterface.h"
#include <conv-ccs.h>

#if(CMK_CCS_AVAILABLE)
#include <stdlib.h>

static int shouldReply=0;
static CcsDelayedReply theReply;

extern "C" void CApplicationDepositNode0Data(char *data)
{
  char *reply;
  int len;

  if(shouldReply == 0) {
    return; 
  }

  len = strlen(data) + 8; /* for the 'namdpr ' and '\0' */
  reply = (char *)malloc(len * sizeof(char));
  strcpy(reply, "namdpr ");
  strcat(reply, data);
  
  /* Do the CcsSendReply */
  CcsSendDelayedReply(theReply, strlen(reply) + 1, reply);
  shouldReply=0;
  free(reply);
}

/*Called by clients on node 0 to ask for perf. data*/
void CApplicationRequestData(void)
{
	shouldReply=1;
	theReply=CcsDelayReply();
}

void CApplicationInit(void)
{
  CcsRegisterHandler("perf_app",(CmiHandler)CApplicationRequestData);
}

#endif
