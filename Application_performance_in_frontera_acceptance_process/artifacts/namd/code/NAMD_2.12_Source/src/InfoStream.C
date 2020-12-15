/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Defines a new stream: iout, for "i"nforming consoles.
*/

#include "InfoStream.h"
#include "charm++.h"
#include "Vector.h"
#include "Tensor.h"
#include <stdio.h>

CkpvExtern(infostream, iout_obj);

infostream& infostream::access_iout_obj() {
  return CkpvAccess(iout_obj);
}

infostream::infostream() {}

infostream::~infostream() {;}

/* output using CkPrintf() (end by inform) */
void infostream::endi() {
  *this << std::ends;
  std::string infostr = str();
  CkPrintf("%s",infostr.c_str());
#ifndef NAMD_NO_STDOUT_FLUSH
  fflush(stdout);  // since CkPrintf doesn't always flush
#endif
  str("");
}

infostream& endi(infostream& s)  { s.endi(); return s; }

std::ostream& iPE(std::ostream& s) {
  return s << "Pe(" << CkMyPe() << ')';
}

std::ostream& operator<<(std::ostream& strm, const Vector &v1) {
       strm << v1.x << " " << v1.y << " " << v1.z;
       return strm;
}

infostream& operator<<(infostream& strm, const Vector &v1) {
       strm << v1.x << " " << v1.y << " " << v1.z;
       return strm;
}


std::ostream& operator<<(std::ostream& strm, const Tensor &t1) {
       strm << t1.xx << " " << t1.xy << " " << t1.xz << " "
            << t1.yx << " " << t1.yy << " " << t1.yz << " "
            << t1.zx << " " << t1.zy << " " << t1.zz;
       return strm;
}

infostream& operator<<(infostream& strm, const Tensor &t1) {
       strm << t1.xx << " " << t1.xy << " " << t1.xz << " "
            << t1.yx << " " << t1.yy << " " << t1.yz << " "
            << t1.zx << " " << t1.zy << " " << t1.zz;
       return strm;
}

/* define how to use the remaining << args */
/** infostream<<ostream (hot to handle inherited modifiers) **/
infostream& infostream::operator<<(std::ostream& (*f)(std::ostream&)) { f(*this); return(*this); }
/** infostream<<infostream (how to handle class modifiers) **/
infostream& infostream::operator<<(infostream& (*f)(infostream&)) { return f(*this); }

#define LOCALMOD(type) infostream& infostream::operator<<(type x) \
		{ (std::ostream&)(*this) << x; return(*this); }
/** << characters **/
LOCALMOD(char)
LOCALMOD(unsigned char)
LOCALMOD(const char *)
/** << integers **/
LOCALMOD(int)
LOCALMOD(long)
LOCALMOD(short)
LOCALMOD(unsigned int)
LOCALMOD(unsigned long)
LOCALMOD(unsigned short)
#ifdef _MSC_VER
LOCALMOD(__int64)
LOCALMOD(unsigned __int64)
#endif
/** << floats **/
LOCALMOD(float)
LOCALMOD(double)
/** << pointers **/
LOCALMOD(void *)
LOCALMOD(std::streambuf *)
#undef LOCALMOD

/** common messages **/
/** iINFO, iWARN, iERROR, iDEBUG provide initial headings. **/
/** iINFOF, iWARNF, iERRORF, iDEBUGF provide initial headings with file name
    and line numbers. **/
std::ostream& iINFO (std::ostream& s)  { return s << "Info: "; }
std::ostream& iWARN (std::ostream& s)  { return s << "Warning: "; }
std::ostream& iERROR(std::ostream& s)  { return s << "ERROR: "; }
std::ostream& iDEBUG(std::ostream& s)  { return s << "DEBUG: "; }


