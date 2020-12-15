/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
    Defines a new stream: iout, for "i"nforming consoles.
*/

#ifndef INFOSTREAM_H
#define INFOSTREAM_H

#include <iostream>
#include <sstream>
class Vector;
class Tensor;

class infostream : public std::ostringstream
{
  public:
  static infostream& access_iout_obj();

  infostream();
  ~infostream();
  void endi();

  /* define how to use the remaining << args */
  /** infostream<<ostream (hot to handle inherited modifiers) **/
  infostream& operator<<(std::ostream& (*f)(std::ostream&));
  /** infostream<<infostream (how to handle class modifiers) **/
  infostream& operator<<(infostream& (*f)(infostream&));

  #define LOCALMOD(type) infostream& operator<<(type x)
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
#ifdef _MSC_VER
  LOCALMOD(__int64);
  LOCALMOD(unsigned __int64);
#endif
  /** << floats **/
  LOCALMOD(float);
  LOCALMOD(double);
  /** << pointers **/
  LOCALMOD(void *);
  LOCALMOD(std::streambuf *);
  #undef LOCALMOD
};

std::ostream& operator<<(std::ostream& strm, const Vector &v1);
infostream& operator<<(infostream& strm, const Vector &v1);

std::ostream& operator<<(std::ostream& strm, const Tensor &t1);
infostream& operator<<(infostream& strm, const Tensor &t1);

/** modifiers **/
infostream& endi(infostream& s);

/** common messages **/
/** iINFO, iWARN, iERROR, iDEBUG provide initial headings. **/
/** iINFOF, iWARNF, iERRORF, iDEBUGF provide initial headings with file name
    and line numbers. **/
std::ostream& iINFO (std::ostream& s);
std::ostream& iWARN (std::ostream& s);
std::ostream& iERROR(std::ostream& s);
std::ostream& iDEBUG(std::ostream& s);
std::ostream& iPE(std::ostream& s);  // outlined because of CkMyPe()

#define iFILE __FILE__<<'('<<__LINE__<<"): "
#define iINFOF  iINFO << iFILE
#define iWARNF  iWARN << iFILE
#define iERRORF  iERROR << iFILE
#define iDEBUGF  iDEBUG << iFILE

#define iout infostream::access_iout_obj()

#endif /* INFOSTREAM_H */

