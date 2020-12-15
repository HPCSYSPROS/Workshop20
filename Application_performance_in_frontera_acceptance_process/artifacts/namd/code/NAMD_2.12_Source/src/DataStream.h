/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
    Defines a new stream: dout, for data logging.
*/

#ifndef DATASTREAM_H
#define DATASTREAM_H

#include "InfoStream.h"

class datastream : public std::ostringstream
{
  public:
  datastream() {}
  ~datastream() {;}

  void endd();

  /* define how to use the remaining << args */
  /** datastream<<ostream (hot to handle inherited modifiers) **/
  datastream& operator<<(std::ostream& (*f)(std::ostream&)) { f(*this); return(*this); }
  /** datastream<<datastream (how to handle class modifiers) **/
  datastream& operator<<(datastream& (*f)(datastream&)) { return f(*this); }

  #define LOCALMOD(type) datastream& operator<<(type x) \
		{ (std::ostream&)(*this) << x; return(*this); }
  /** << characters **/
  LOCALMOD(char);
  LOCALMOD(unsigned char);
  LOCALMOD(const char *);
  /** << integers **/
  LOCALMOD(int);
  LOCALMOD(long);
  LOCALMOD(short);
  LOCALMOD(unsigned int);
  LOCALMOD(unsigned long);
  LOCALMOD(unsigned short);
  /** << floats **/
  LOCALMOD(float);
  LOCALMOD(double);
  /** << pointers **/
  LOCALMOD(void *);
  LOCALMOD(std::streambuf *);
  #undef LOCALMOD
};

datastream& operator<<(datastream& strm, const Vector &v1);

datastream& operator<<(datastream& strm, const Tensor &t1);

/** modifiers **/
inline datastream& endd(datastream& s)  { s.endd(); return s; }

#define iFILE __FILE__<<'('<<__LINE__<<"): "
#define iINFOF  iINFO << iFILE
#define iWARNF  iWARN << iFILE
#define iERRORF  iERROR << iFILE
#define iDEBUGF  iDEBUG << iFILE


extern datastream dout;

#endif /* DATASTREAM_H */

