/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "InfoStream.h"
#include "FreeEnergyEnums.h"
#include "FreeEnergyAssert.h"
#include "Vector.h"
#include "FreeEnergyVector.h"

AVector::AVector(Vector_t Type) {
//------------------------------------------------
// init this vector to (0,0,0) for kRegular Type
// init this vector to (x,y,z) for kRandom Type
// where x, y, and z are in the range (0:1)
//
// Note:  The assumption is made that srand() has
//        already been called.
//------------------------------------------------
  if (Type == kRegular) {
    m_x = 0.0;
    m_y = 0.0;
    m_z = 0.0;
  }
  else {
    ASSERT(Type == kRandom);
    Random();
  }
}


AVector::AVector(double x, double y, double z) {
//-------------------------------------------
// init this vector to (x,y,z)
//-------------------------------------------
  m_x = x;
  m_y = y;
  m_z = z;
}


AVector::AVector(const AVector& Vector) {
//-------------------------------------------
// init this vector to Vector
//-------------------------------------------
  m_x = Vector.m_x;
  m_y = Vector.m_y;
  m_z = Vector.m_z;
}


AVector& AVector::Random() {
//------------------------------------------------
// replace this vector with a new random vector
// where x, y, and z are in the range (0:1)
//------------------------------------------------
  double RandMax;
   
  RandMax = RAND_MAX;
  m_x = (double)rand() / RandMax;
  m_y = (double)rand() / RandMax;
  m_z = (double)rand() / RandMax;
  ASSERT(m_x <= 1.0);
  ASSERT(m_y <= 1.0);
  ASSERT(m_z <= 1.0);
  return(*this);
}


Bool_t AVector::Set(double x, double y, double z) {
//-------------------------------------------
// int this vector to (x,y,z)
//-------------------------------------------
  m_x = x;
  m_y = y;
  m_z = z;
  return(kTrue);
}


AVector AVector::operator+ (const AVector& vector) {
//----------------------------------------------------------
// create a new vector: the sum of this one and passed one
//----------------------------------------------------------
  AVector Point(vector.m_x+m_x, vector.m_y+m_y, vector.m_z+m_z);
  return(Point);
}


AVector AVector::operator- (const AVector& vector) {
//----------------------------------------------------------
// create a new vector: this one minus the passed one
//----------------------------------------------------------
  AVector Point(m_x-vector.m_x, m_y-vector.m_y, m_z-vector.m_z);
  return(Point);
}


AVector AVector::operator/ (double divisor) {
//---------------------------------------------------
// create a new vector: this one divided by divisor
//---------------------------------------------------
  AVector Point(m_x/divisor, m_y/divisor, m_z/divisor);
  return(Point);
}


AVector& AVector::operator/= (double divisor) {
//------------------------------------------------------
// divide this vector by divisor and return a ref to it
//------------------------------------------------------
  m_x /= divisor;
  m_y /= divisor;
  m_z /= divisor;
  return(*this);
}


AVector& AVector::operator*= (double scalar) {
//------------------------------------------------------
// multiply this vector by scalar and return a ref to it
//------------------------------------------------------
  m_x *= scalar;
  m_y *= scalar;
  m_z *= scalar;
  return(*this);
}


AVector AVector::operator* (double scalar) {
//---------------------------------------------------
// create a new vector: this one divided by divisor
//---------------------------------------------------
  AVector Point(m_x*scalar, m_y*scalar, m_z*scalar);
  return(Point);
}


void SetEqual(AVector& Vec1, const Vector& Vec2) {
//------------------------------------------------------------------------
// used for casting Vector -> AVector
//------------------------------------------------------------------------
  Vec1.m_x = Vec2.x;
  Vec1.m_y = Vec2.y;
  Vec1.m_z = Vec2.z;
}


void SetEqual(Vector& Vec1,  const AVector& Vec2) {
//------------------------------------------------------------------------
// used for casting AVector -> Vector
//------------------------------------------------------------------------
  Vec1.x = Vec2.m_x;
  Vec1.y = Vec2.m_y;
  Vec1.z = Vec2.m_z;
}


AVector& AVector::operator= (const AVector& vector) {
//-------------------------------------------
// set this vector to the passed one
//-------------------------------------------
  m_x = vector.m_x;
  m_y = vector.m_y;
  m_z = vector.m_z;
  return(*this);
}


AVector& AVector::operator+= (const AVector& vector) {
//--------------------------------------------------
// set this vector to this one plus the passed one
//--------------------------------------------------
  m_x += vector.m_x;
  m_y += vector.m_y;
  m_z += vector.m_z;
  return(*this);
}


double& AVector::operator[] (int index) {
//-------------------------------------------------------------------
// return one element of this vector
// note: this op is used to get AND set an element
//       (since it returns double&)
//-------------------------------------------------------------------
  ASSERT( (index>=0) && (index<3) );

  switch(index) {
    case 0:  return(m_x);
    case 1:  return(m_y);
    case 2:  return(m_z);
    default: return(m_x);   // should never get here
  }
}


double AVector::Dist() {
//-------------------------------------------------------------------
// calculate distance from this point to (0, 0, 0)
//-------------------------------------------------------------------
  return( sqrt(m_x*m_x + m_y*m_y + m_z*m_z) );
}


double AVector::DistSqr() {
//-------------------------------------------------------------------
// calculate distance-squared from this point to (0, 0, 0)
//-------------------------------------------------------------------
  return(m_x*m_x + m_y*m_y + m_z*m_z);
}


double AVector::Dist(const AVector& Vector) {
//-------------------------------------------------------------------
// calculate distance between this point and Vector
//-------------------------------------------------------------------
  double d1 = (m_x - Vector.m_x);
  double d2 = (m_y - Vector.m_y);
  double d3 = (m_z - Vector.m_z);
  return( sqrt(d1*d1 + d2*d2 + d3*d3) );
}


double AVector::DistSqr(const AVector& Vector) {
//-------------------------------------------------------------------
// calculate distance-squared between this point and Vector
//-------------------------------------------------------------------
  double d1 = (m_x - Vector.m_x);
  double d2 = (m_y - Vector.m_y);
  double d3 = (m_z - Vector.m_z);
  return(d1*d1 + d2*d2 + d3*d3);
}


AVector AVector::cross(const AVector& Vector) {
//-------------------------------------------------------------------
// calculate this vector crossed with Vector (this x Vector).
// see Mathematical Handbook, p.118.
//-------------------------------------------------------------------
  AVector CrossProduct;

  CrossProduct.m_x = (m_y * Vector.m_z) - (m_z * Vector.m_y);
  CrossProduct.m_y = (m_z * Vector.m_x) - (m_x * Vector.m_z);
  CrossProduct.m_z = (m_x * Vector.m_y) - (m_y * Vector.m_x);
  return(CrossProduct);
}


double AVector::dot(const AVector& Vector) {
//-------------------------------------------------------------------
// calculate dot product of this vector and Vector
// see Mathematical Handbook, p.118.
//-------------------------------------------------------------------
  return(m_x*Vector.m_x + m_y*Vector.m_y + m_z*Vector.m_z);
}


void AVector::Out() {
//-------------------------------------------------------------------
// write it
//-------------------------------------------------------------------
  char  Str1[20], Str2[20], Str3[20];

  sprintf(Str1, "%8.3f", m_x);
  sprintf(Str2, "%8.3f", m_y);
  sprintf(Str3, "%8.3f", m_z);
  iout << "(" << Str1 << "," << Str2 << "," << Str3 << ")";
}


void AVector::Output() {
//-------------------------------------------------------------------
// write it to standard output
//-------------------------------------------------------------------
  char  Word1[20], Word2[20], Word3[20];

  if ( (fabs(m_x)<99999) && (fabs(m_y)<99999) && (fabs(m_z)<99999) ) {
    sprintf(Word1, "%10.3f", m_x);
    sprintf(Word2, "%10.3f", m_y);
    sprintf(Word3, "%10.3f", m_z);
  }
  else {
    sprintf(Word1, "%10.2e", m_x);
    sprintf(Word2, "%10.2e", m_y);
    sprintf(Word3, "%10.2e", m_z);
  }
  iout << "( " << Word1 << " " << Word2 << " " << Word3 << " )";
}


AVector& AVector::Scale(AVector& SmallVec, AVector& BigVec) {
//-------------------------------------------------------------------
// scale this vector, whose (x,y,z) are in the range (0:1),
// to be in the range (SmallVec:BigVec)
//-------------------------------------------------------------------
  m_x = SmallVec.m_x + (BigVec.m_x - SmallVec.m_x) * m_x;
  m_y = SmallVec.m_y + (BigVec.m_y - SmallVec.m_y) * m_y;
  m_z = SmallVec.m_z + (BigVec.m_z - SmallVec.m_z) * m_z;
  return(*this);
}

