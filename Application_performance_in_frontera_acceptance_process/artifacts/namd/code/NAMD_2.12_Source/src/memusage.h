/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef MEMUSAGE_H
#define MEMUSAGE_H

unsigned long memusage(const char **source = 0);

inline double memusage_kB() { return memusage() / 1024.; }
inline double memusage_MB() { return memusage() / 1048576.; }

class memusageinit {
public:
  memusageinit();
private:
  static int initialized;
  static unsigned long sbrkval;
  static unsigned long memusage_sbrk();
  friend unsigned long memusage(const char **source);
};

static memusageinit memusageinitobject;

#endif

