/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

// written by David Hurwitz.

#if !defined(VECTOR_HPP)
  #define VECTOR_HPP

enum Vector_t {kRandom, kRegular};

class AVector {

protected:
  double m_x;
  double m_y;
  double m_z;

public:
  AVector(Vector_t Type=kRegular);
  AVector(double x, double y, double z);
  AVector(const AVector& Vector);
  Bool_t   Set(double x, double y, double z);
  AVector  operator+  (const AVector& Vector);
  AVector  operator-  (const AVector& Vector);
  AVector  operator/  (double divisor);
  AVector& operator/= (double divisor);
  AVector& operator*= (double scalar);
  AVector  operator*  (double scalar);
  AVector& operator=  (const AVector& Vector);
  AVector& operator+= (const AVector& Vector);
  AVector  cross(const AVector& Vector);
  double   dot(const AVector& Vector);
  double&  operator[] (int index);
  double   Dist();
  double   DistSqr();
  double   Dist(const AVector& Vector);
  double   DistSqr(const AVector& Vector);
  void     Out();
  void     Output();
  AVector& Scale(AVector& SmallVec, AVector& BigVec);
  AVector& Random();
  // used for casting Vector -> AVector and vice-versa
  // (I had difficulty overloading operator=)
  friend void SetEqual(AVector& Vec1, const Vector& Vec2);
  friend void SetEqual(Vector& Vec1,  const AVector& Vec2);
};

#endif
