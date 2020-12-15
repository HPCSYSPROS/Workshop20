/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Defines a new stream: dout, for data logs.
*/

#include "DataStream.h"
#include "CollectionMgr.h"
#include "Vector.h"
#include "Tensor.h"
#include <stdio.h>

/* output using CkPrintf() (end by inform) */
void datastream::endd() {
  *this << std::ends;
  std::string datastr = str();
  CollectionMgr::Object()->sendDataStream(datastr.c_str());
  str("");
}

datastream& operator<<(datastream& strm, const Vector &v1) {
       strm << v1.x << " " << v1.y << " " << v1.z;
       return strm;
}

datastream& operator<<(datastream& strm, const Tensor &t1) {
       strm << t1.xx << " " << t1.xy << " " << t1.xz << " "
            << t1.yx << " " << t1.yy << " " << t1.yz << " "
            << t1.zx << " " << t1.zy << " " << t1.zz;
       return strm;
}

datastream dout;

