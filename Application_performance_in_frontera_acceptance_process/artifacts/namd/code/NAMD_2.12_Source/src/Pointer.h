#ifndef _POINTER_H
#define _POINTER_H

class Pointer {
  public:

  void *data;

  Pointer() { }
  Pointer(void *_data): data(_data) { }

  void pup(PUP::er &p) {
    pup_bytes(&p,&data,sizeof(data));
  }
};

#endif

