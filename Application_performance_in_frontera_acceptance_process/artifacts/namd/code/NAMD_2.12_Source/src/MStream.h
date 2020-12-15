/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef MSTREAM_H
#define MSTREAM_H

#include "Vector.h"
#include <string.h>

class StreamMessage;
class Communicate;

class MIStream {
  private:
    int PE, tag;
    StreamMessage *msg;
    size_t currentPos;
    int currentIndex;
    StreamMessage *early;
    Communicate *cobj;
    unsigned int checksum;
    MIStream *Get(char *buf, size_t len);  // get len bytes from message to buf
  public:
    MIStream(Communicate *c, int pe, int tag);
    ~MIStream();
    MIStream *get(char &data) { 
      return Get(&data,sizeof(char)); 
    }
    MIStream *get(unsigned char &data) { 
      return Get((char *)&data,sizeof(unsigned char)); 
    }
    MIStream *get(short &data) { 
      return Get((char *)&data, sizeof(short)); 
    }
    MIStream *get(unsigned short &data) { 
      return Get((char *)&data, sizeof(unsigned short)); 
    }
    MIStream *get(int &data) { 
      return Get((char *)&data, sizeof(int)); 
    }
    MIStream *get(unsigned int &data) { 
      return Get((char *)&data, sizeof(unsigned int)); 
    }
    MIStream *get(long &data) { 
      return Get((char *)&data, sizeof(long)); 
    }
    MIStream *get(unsigned long &data) { 
      return Get((char *)&data, sizeof(unsigned long)); 
    }
    MIStream *get(float &data) { 
      return Get((char *)&data, sizeof(float)); 
    }
    MIStream *get(double &data) { 
      return Get((char *)&data, sizeof(double)); 
    }
    MIStream *get(size_t len, char *data) { 
      return Get(data,len*sizeof(char)); 
    }
    MIStream *get(size_t len, unsigned char *data) { 
      return Get((char *)data,len*sizeof(unsigned char)); 
    }
    MIStream *get(size_t len, short *data) { 
      return Get((char *)data,len*sizeof(short)); 
    }
    MIStream *get(size_t len, unsigned short *data) { 
      return Get((char *)data,len*sizeof(unsigned short)); 
    }
    MIStream *get(size_t len, int *data) { 
      return Get((char *)data,len*sizeof(int)); 
    }
    MIStream *get(size_t len, unsigned int *data) { 
      return Get((char *)data,len*sizeof(unsigned int)); 
    }
    MIStream *get(size_t len, long *data) { 
      return Get((char *)data,len*sizeof(long)); 
    }
    MIStream *get(size_t len, unsigned long *data) { 
      return Get((char *)data,len*sizeof(unsigned long)); 
    }
    MIStream *get(size_t len, float *data) { 
      return Get((char *)data,len*sizeof(float)); 
    }
    MIStream *get(size_t len, double *data) { 
      return Get((char *)data,len*sizeof(double)); 
    }
    MIStream *get(size_t len, Vector *data) {
      return Get((char *)data, len*sizeof(Vector));
    }
    MIStream *get(Vector *data) {
      return Get((char *)data, sizeof(Vector));
    }
};

class MOStream {
  private:
    int PE, tag;
    unsigned int bufLen;
    StreamMessage *msgBuf;
    Communicate *cobj;
    MOStream *Put(char *buf, size_t len); // put len bytes from buf into message
  public:
    MOStream(Communicate *c, int pe, int tag, size_t bufSize);
    ~MOStream();
    void end(void);
    MOStream *put(char data) { 
      return Put(&data,sizeof(char)); 
    }
    MOStream *put(unsigned char data) { 
      return Put((char *)&data,sizeof(unsigned char)); 
    }
    MOStream *put(short data) { 
      return Put((char *)&data, sizeof(short)); 
    }
    MOStream *put(unsigned short data) { 
      return Put((char *)&data, sizeof(unsigned short)); 
    }
    MOStream *put(int data) { 
      return Put((char *)&data, sizeof(int)); 
    }
    MOStream *put(unsigned int data) { 
      return Put((char *)&data, sizeof(unsigned int)); 
    }
    MOStream *put(long data) { 
      return Put((char *)&data, sizeof(long)); 
    }
    MOStream *put(unsigned long data) { 
      return Put((char *)&data, sizeof(unsigned long)); 
    }
    MOStream *put(float data) { 
      return Put((char *)&data, sizeof(float)); 
    }
    MOStream *put(double data) { 
      return Put((char *)&data, sizeof(double)); 
    }
    MOStream *put(size_t len, char *data) { 
      return Put(data,len*sizeof(char)); 
    }
    MOStream *put(size_t len, unsigned char *data) { 
      return Put((char *)data,len*sizeof(unsigned char)); 
    }
    MOStream *put(size_t len, short *data) { 
      return Put((char *)data,len*sizeof(short)); 
    }
    MOStream *put(size_t len, unsigned short *data) { 
      return Put((char *)data,len*sizeof(unsigned short)); 
    }
    MOStream *put(size_t len, int *data) { 
      return Put((char *)data,len*sizeof(int)); 
    }
    MOStream *put(size_t len, unsigned int *data) { 
      return Put((char *)data,len*sizeof(unsigned int)); 
    }
    MOStream *put(size_t len, long *data) { 
      return Put((char *)data,len*sizeof(long)); 
    }
    MOStream *put(size_t len, unsigned long *data) { 
      return Put((char *)data,len*sizeof(unsigned long)); 
    }
    MOStream *put(size_t len, float *data) { 
      return Put((char *)data,len*sizeof(float)); 
    }
    MOStream *put(size_t len, double *data) { 
      return Put((char *)data,len*sizeof(double)); 
    }
    MOStream *put(Vector *data) { 
      return Put((char *)data,sizeof(Vector)); 
    }
    MOStream *put(size_t len, Vector *data) { 
      return Put((char *)data,len*sizeof(Vector)); 
    }
    MOStream *put(char *data) {
      size_t length = strlen(data);
      put((int)length);
      return put(length, data);
    }
};

#endif
