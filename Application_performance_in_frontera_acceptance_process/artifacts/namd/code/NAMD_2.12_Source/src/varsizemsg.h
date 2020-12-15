/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef VARSIZEMSG_H
#define VARSIZEMSG_H

/* Macros for automatically allocating, packing and unpacking varsize messages.

Usage:

class MyMsg : public CMessage_MyMsg {
public:

  int myscalar;
  int *myarray;
  double myscalar2;
  double *myarray2;

  VARSIZE_DECL(MyMsg);
};

VARSIZE_MSG(MyMsg,
  VARSIZE_ARRAY(myarray);  // array size passed in via new
  VARSIZE_ARRAY(myarray2);  // array size passed in via new
)

*/

#define VARSIZE_DECL(MSGTYPE) \
  static void *alloc(int, int, int *, int); \
  static void *pack(MSGTYPE *); \
  static MSGTYPE* unpack(void *)

#define VARSIZE_MSG(MSGTYPE,MSGDATA) \
void *MSGTYPE::alloc(int varsizemsg_msgnum, int varsizemsg_size, \
                   int *varsizemsg_array, int varsizemsg_priobits) { \
  varsizemsg_size = ALIGN8(varsizemsg_size); \
  MSGTYPE *varsizemsg_msg = 0; \
  { \
    const int varsizemsg_pass = 0; \
    int varsizemsg_totalsize = varsizemsg_size; \
    int varsizemsg_arraycount = 0; \
    MSGDATA \
    varsizemsg_msg = (MSGTYPE *) CkAllocMsg( \
	varsizemsg_msgnum, varsizemsg_totalsize, varsizemsg_priobits); \
  } \
  { \
    const int varsizemsg_pass = 1; \
    int varsizemsg_totalsize = varsizemsg_size; \
    int varsizemsg_arraycount = 0; \
    MSGDATA \
  } \
  return (void *) varsizemsg_msg; \
} \
 \
void *MSGTYPE::pack(MSGTYPE *varsizemsg_msg) { \
  int *varsizemsg_array=0, varsizemsg_arraycount=0, varsizemsg_totalsize=0; \
  { \
    const int varsizemsg_pass = 2; \
    MSGDATA \
  } \
  return (void *) varsizemsg_msg; \
} \
 \
MSGTYPE *MSGTYPE::unpack(void *varsizemsg_buf) { \
  int *varsizemsg_array=0, varsizemsg_arraycount=0, varsizemsg_totalsize=0; \
  MSGTYPE *varsizemsg_msg = (MSGTYPE *) varsizemsg_buf; \
  { \
    const int varsizemsg_pass = 3; \
    MSGDATA \
  } \
  return varsizemsg_msg; \
}

template<class T> inline T* cast_array(T*, char *a) { return (T*) a; }
template<class T> inline T* cast_size(T*, size_t a) { return (T*) a; }

#define VARSIZE_ARRAY(ARRAY) { \
  int varsizemsg_arraysize; \
  switch ( varsizemsg_pass ) { \
  case 0: \
    varsizemsg_arraysize = sizeof(*(varsizemsg_msg->ARRAY)) * \
		varsizemsg_array[varsizemsg_arraycount]; \
    varsizemsg_totalsize += ALIGN8(varsizemsg_arraysize); \
    varsizemsg_arraycount++; \
    break; \
  case 1: \
    varsizemsg_msg->ARRAY = cast_array(varsizemsg_msg->ARRAY, \
	(char *) varsizemsg_msg + varsizemsg_totalsize); \
    varsizemsg_arraysize = sizeof(*(varsizemsg_msg->ARRAY)) * \
		varsizemsg_array[varsizemsg_arraycount]; \
    varsizemsg_totalsize += ALIGN8(varsizemsg_arraysize); \
    varsizemsg_arraycount++; \
    break; \
  case 2: \
    varsizemsg_msg->ARRAY = cast_size(varsizemsg_msg->ARRAY, \
      (char *) (varsizemsg_msg->ARRAY) - (char *) &(varsizemsg_msg->ARRAY) ); \
    break; \
  case 3: \
    varsizemsg_msg->ARRAY = cast_array(varsizemsg_msg->ARRAY, \
      (char *) &(varsizemsg_msg->ARRAY) + (size_t) (varsizemsg_msg->ARRAY) ); \
    break; \
  default: \
    break; \
  } \
}

#endif

