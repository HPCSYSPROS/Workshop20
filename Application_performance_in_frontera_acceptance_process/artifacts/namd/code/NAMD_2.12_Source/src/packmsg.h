/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef PACKMSG_H
#define PACKMSG_H

/* Macros for automatically packing and unpacking messages.

Usage:

PACK_MSG(MyMsg,
  PACK_THIS;  // for the lazy, bit-copies the whole object, may be unsafe
  PACK(myint);
  PACK(myfloat);
  PACK_RESIZE(myvector);  // needs size(), resize(), and begin()
  PACK_ARRAY(myarray,n);  // n must be a message field and is also sent
  PACK_AND_NEW_ARRAY(myarray2,n);  // also calls new [] but not delete []
)

*/

#define ALIGN_8(x)   (((unsigned long)x + 7) & (~7))

#define PACKMSG_CHECKSUM(X)

template<class T> class ResizeArray;

template<class T> inline size_t sizeof_element(ResizeArray<T> &) { return sizeof(T); }

template<class T> inline T* new_array(T*, int n) { return new T[n]; }

#define PACK_MSG(MSGTYPE,MSGDATA) \
void *MSGTYPE::pack(MSGTYPE *packmsg_msg) { \
  PACKMSG_CHECKSUM(unsigned int packmsg_checksum = 0;) \
  int packmsg_size = 0; \
  char *packmsg_cur = 0; \
  { \
    const int packmsg_pass = 0; \
    PACKMSG_CHECKSUM( \
      packmsg_size += sizeof(packmsg_checksum); \
      PACK_MEMORY(&packmsg_size,sizeof(packmsg_size)); \
    ) \
    MSGDATA \
  } \
  void *packmsg_buf = CkAllocBuffer(packmsg_msg,packmsg_size); \
  packmsg_cur = (char *)packmsg_buf; \
  { \
    const int packmsg_pass = 1; \
    PACKMSG_CHECKSUM( \
      packmsg_cur += sizeof(packmsg_checksum); \
      PACK_MEMORY(&packmsg_size,sizeof(packmsg_size)); \
    ) \
    MSGDATA \
  } \
  PACKMSG_CHECKSUM( \
    packmsg_cur = (char *)packmsg_buf; \
    for ( int i=sizeof(packmsg_checksum); i < packmsg_size; i++ ) { \
      packmsg_checksum += (unsigned char) packmsg_cur[i]; \
    } \
    CmiMemcpy(packmsg_buf,(void *)&packmsg_checksum,sizeof(packmsg_checksum)); \
  ) \
  delete packmsg_msg; \
  return packmsg_buf; \
} \
 \
MSGTYPE *MSGTYPE::unpack(void *packmsg_buf) { \
  PACKMSG_CHECKSUM( \
    unsigned int packmsg_checksum = 0; \
    unsigned int packmsg_checksum_orig = 0; \
  ) \
  int packmsg_size = 0; \
  void *packmsg_msg_ = CkAllocBuffer(packmsg_buf,sizeof(MSGTYPE)); \
  MSGTYPE *packmsg_msg = new (packmsg_msg_) MSGTYPE; \
  char *packmsg_cur = (char *)packmsg_buf; \
  { \
    const int packmsg_pass = 2; \
    PACKMSG_CHECKSUM( \
      CmiMemcpy((void *)&packmsg_checksum_orig,(void *)packmsg_cur, \
				sizeof(packmsg_checksum)); \
      packmsg_cur += sizeof(packmsg_checksum); \
      PACK_MEMORY(&packmsg_size,sizeof(packmsg_size)); \
      char *packmsg_cur2 = (char *)packmsg_buf; \
      for ( int i=sizeof(packmsg_checksum); i < packmsg_size; i++ ) { \
        packmsg_checksum += (unsigned char) packmsg_cur2[i]; \
      } \
      if ( packmsg_checksum != packmsg_checksum_orig ) { \
        char errmsg[256]; \
        sprintf(errmsg,"PACKMSG checksums do not agree!  %s(%d): %d vs %d", \
	__FILE__, __LINE__, packmsg_checksum, packmsg_checksum_orig); \
        NAMD_bug(errmsg); \
      } \
    ) \
    MSGDATA \
  } \
  CkFreeMsg(packmsg_buf); \
  return packmsg_msg; \
}

#define PACK_MEMORY(BUF,SIZE) { \
  int ASIZE = ALIGN_8(SIZE); \
  switch ( packmsg_pass ) { \
  case 0: \
    packmsg_size += (ASIZE); \
    break; \
  case 1: \
    CmiMemcpy((void *)packmsg_cur,(void *)(BUF),(SIZE)); \
    packmsg_cur += (ASIZE); \
    break; \
  case 2: \
    CmiMemcpy((void *)(BUF),(void *)packmsg_cur,(SIZE)); \
    packmsg_cur += (ASIZE); \
    break; \
  default: \
    break; \
  } \
}

#define PACK_THIS PACK_MEMORY(packmsg_msg,sizeof(*packmsg_msg));

#define PACK(DATA) PACK_MEMORY(&(packmsg_msg->DATA),sizeof(packmsg_msg->DATA))

#define PACK_RESIZE(DATA) { \
  int packmsg_array_len = packmsg_msg->DATA.size(); \
  PACK_MEMORY(&packmsg_array_len,sizeof(packmsg_array_len)); \
  if ( packmsg_pass == 2 ) packmsg_msg->DATA.resize(packmsg_array_len); \
  int packmsg_array_size = \
    packmsg_array_len * sizeof_element(packmsg_msg->DATA); \
  PACK_MEMORY(packmsg_msg->DATA.begin(),packmsg_array_size); \
}

#define PACK_ARRAY(DATA,LEN) { \
  PACK(LEN); \
  PACK_MEMORY(packmsg_msg->DATA,packmsg_msg->LEN*sizeof(*(packmsg_msg->DATA))); \
}

#define PACK_AND_NEW_ARRAY(DATA,LEN) { \
  PACK(LEN)\
  if ( packmsg_pass == 2 ) { \
    packmsg_msg->DATA = new_array(packmsg_msg->DATA,packmsg_msg->LEN); \
  } \
  PACK_MEMORY(packmsg_msg->DATA,packmsg_msg->LEN*sizeof(*(packmsg_msg->DATA))); \
}

#endif

