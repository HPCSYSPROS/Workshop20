/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef VECTOR_H
#define VECTOR_H

#include <math.h>
#include <stdio.h>
#include "common.h"

class Vector;

class FloatVector {
 public:
  float x,y,z;
  inline FloatVector(void) { ; }
  inline FloatVector(const Vector &v);
};

#ifdef ARCH_POWERPC

#include "builtins.h"

inline double namd_rsqrt(double x)
{
  double r0, r1, r2;

  /*------------------------------------------*/
  /* use reciprocal sqrt estimate instruction */
  /*------------------------------------------*/
  r0 = __frsqrte(x);

  /*----------------------*/
  /* 1st newton iteration */
  /*----------------------*/
  r1 = r0 + 0.5*r0*(1.0 - x*r0*r0);

  /*----------------------*/
  /* 2nd newton iteration */
  /*----------------------*/
  r2 = r1 + 0.5*r1*(1.0 - x*r1*r1);

  return r2;
}

inline double namd_reciprocal (double x) {
  register double rx;
  
  rx = __fres(x);             // 0th estimate (13 bits)
  rx = rx + rx*(1.0 - x*rx);  // 1st Newton iteration (26 bits)
  rx = rx + rx*(1.0 - x*rx);  // 2nd Newton iteration (52 bits = full mantissa for a double)

  return rx;
}

#else
#define namd_rsqrt(x)        (1.0 / sqrt (x))
#define namd_reciprocal(x)   (1.0 / x)
#endif

class Vector {
   public:
     BigReal x,y,z;
     
     inline Vector(void) : x(-99999), y(-99999), z(-99999) { ; }
//     inline Vector(void) { ; }
     
     inline Vector( BigReal newx, BigReal newy, BigReal newz)
       : x(newx), y(newy), z(newz) { ; }

     inline Vector( BigReal newv )  // allow Vector v = 0; etc.
       : x(newv), y(newv), z(newv) { ; }

     inline Vector(const FloatVector &v) : x(v.x), y(v.y), z(v.z) { ; }

     inline BigReal &operator[](int i) {
       return i==0 ? x
             :i==1 ? y
             :i==2 ? z
             :(NAMD_die("Vector reference out of bounds."), x);

     }

     //  v1 = const;
     inline Vector& operator=(const BigReal &v2) {
       x = v2;
       y = v2;
       z = v2;
       return *this;
     }

     //  v1 += v2;
     inline void operator+=(const Vector &v2) {
       x += v2.x;
       y += v2.y;
       z += v2.z;
     }

     // v1 -= v2;
     inline void operator-=(const Vector &v2) {
       x -= v2.x;
       y -= v2.y;
       z -= v2.z;
     }

     // v1 *= const
     inline void operator*=(const BigReal &v2) {
       x *= v2;
       y *= v2;
       z *= v2;
     }

     // v1 /= const
     inline void operator/=(const BigReal& v2) {
       BigReal v2_recip = namd_reciprocal(v2);
       x *= v2_recip;
       y *= v2_recip;
       z *= v2_recip;
     }

     inline friend int operator == (const Vector& v1, const Vector& v2) {
       return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
     }
     inline friend int operator != (const Vector& v1, const Vector& v2) {
       // return !(v1.x == v2.x && v1.y == v2.y && v1.z == v2.z);
       return v1.x != v2.x || v1.y != v2.y || v1.z != v2.z;
     }

     // addition of two vectors
     inline friend Vector operator+(const Vector& v1, const Vector& v2) {
       return Vector( v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
     }

     // negation
     inline friend Vector operator-(const Vector &v1) {
       return Vector( -v1.x, -v1.y, -v1.z);
     }

     // subtraction
     inline friend Vector operator-(const Vector &v1, const Vector &v2) {
       return Vector( v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
     }
     // inner ("dot") product
     inline friend BigReal operator*(const Vector &v1, const Vector &v2) {
       return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
     }
     // scalar product
     inline friend Vector operator*(const BigReal &f, const Vector &v1) {
       return Vector(f*v1.x, f*v1.y, f*v1.z);
     }
     // scalar product
     inline friend Vector operator*(const Vector &v1, const BigReal &f) {
       return Vector(f*v1.x, f*v1.y, f*v1.z);
     }
     // division by a scalar
     inline friend Vector operator/(const Vector &v1, const BigReal &f) {
//       if (!f)
//         NAMD_die("Division by 0 on a vector operation.");
       return Vector(v1.x/f, v1.y/f, v1.z/f);
     }
     
     // return the norm
     inline BigReal length(void) const {
       return sqrt(x*x+y*y+z*z);
     }
     
     inline BigReal length2(void) const {
       return (x*x + y*y + z*z);
     }

     inline BigReal rlength (void) {
       return namd_rsqrt (x*x + y*y + z*z);
     }

     // return the unit vector in the same direction
     inline Vector unit(void) const {
       return Vector(x, y, z)/length();
     }
     
     
     // one cross product  v3 = cross(v1, v2)
     inline friend Vector cross(const Vector &v1, const Vector &v2) {
       return Vector( v1.y*v2.z-v2.y*v1.z,
                     // -v1.x*v2.z+v2.x*v1.z,
                      v2.x*v1.z-v1.x*v2.z,
                      v1.x*v2.y-v2.x*v1.y  );
     }

     // multiplying a cross product by a scalar is very common
     // one cross product  v3 = k*cross(v1, v2)
     inline friend Vector cross(const Real &k, const Vector &v1, const Vector &v2) {
       return Vector( k*(v1.y*v2.z-v2.y*v1.z),
                      // k*(-v1.x*v2.z+v2.x*v1.z),
                      k*(v2.x*v1.z-v1.x*v2.z),
                      k*(v1.x*v2.y-v2.x*v1.y) );
     }

     inline friend Vector cross(const BigReal &k, const Vector &v1, const Vector &v2) {
       return Vector( k*(v1.y*v2.z-v2.y*v1.z),
                      // k*(-v1.x*v2.z+v2.x*v1.z),
                      k*(v2.x*v1.z-v1.x*v2.z),
                      k*(v1.x*v2.y-v2.x*v1.y) );
     }

     // A = A x B  -- why do you want this function, anyway?
     void cross(const Vector &v2) {
       BigReal xx =  y*v2.z-v2.y*z;
       // BigReal yy = -x*v2.z+v2.x*z;
       BigReal yy = v2.x*z-x*v2.z;
       z =  x*v2.y-v2.x*y;
       y=yy;
       x=xx;
     }

     // returns (*this) * V2
     BigReal dot(const Vector &v2) {
       return x*v2.x + y*v2.y + z*v2.z;
     }

     // set the vector based on a string.  If bad, return FALSE
     // the string can be in the form "x y z" or "x, y, z"
     Bool set(const char *s) {
	double a[3];    // read into doubles, since I don't know what
	char tmp[100];  // a "BigReal" is in real life
	// cheap way to get commas, etc.  a poor regex
       int i=sscanf(s, "%lf%99[ \t,]%lf%99[ \t,]%lf%99s",
                    a, tmp, a+1, tmp, a+2, tmp);
       if (i != 5) return FALSE;
       const char *t = s;       // now count commas (for "1,,,,2,  , 3")
       int flg = 0;                 // and check for "1 2,,3"
       i = 0;
       for (;*t;t++) {
          if (*t == ',') { 
             if (flg == 0) {   // expecting non-whitespace
                return FALSE;  //    so error
             }
             flg = 0;          // now expect non-whitespace
             i++;              // and increment comma counter
          }
          else if (*t != ' ' && *t != '\t') {  // got non-whitespace
             flg = 1;          // so the next can be whitespace or commas
          }
       }
       if (i == 0 || i == 2) {  // allow "1 2 3" or "1, 2,3" forms
          x = a[0]; y = a[1]; z = a[2];
          return TRUE;
       }
       return FALSE;
     }
};

class zVector : public Vector {
  public:
     inline zVector(void) : Vector(0,0,0) { ; }
     inline zVector(const Vector &v) : Vector(v) { ; }
};


class AlignVector : public Vector {
 public:
  BigReal pad;
  inline AlignVector(void) : Vector(0,0,0) { pad = 0.0; }
  inline AlignVector(const Vector &v) : Vector(v) { pad = 0.0; }

  inline AlignVector( BigReal newx, BigReal newy, BigReal newz)
    : Vector (newx, newy, newz) { pad = 0.0; }
  
  inline AlignVector( BigReal newv )  // allow Vector v = 0; etc.
    : Vector (newv) { pad = 0.0; }
  
  inline AlignVector(const FloatVector &v) : Vector (v) { pad = 0.0; }
};


inline FloatVector::FloatVector(const Vector &v) : x(v.x), y(v.y), z(v.z) { ; }

//#define TEST_VECTOR_CLASS
#ifdef TEST_VECTOR_CLASS
main()
{
  Vector v1(1.1,2.2, 3.3);
  Vector v2(-1, 55, 32.1);
  Vector v3(v1+2*v2);
  Vector v4;
  std::cout << v1 << "  " << v2 << "  " << v3 << "  " << v4 << '\n';
  std::cout << v1*v2 << "  "  << v3-v1-2*v2 <<"  "<< v2 * v3 <<"  "<< v3*v2 <<'\n';
  v4 = v3*5 - v2/4;
  std::cout << v4 << "  " << v3*5.0 - v2/4.0 << '\n';
  std::cout << v4[0] << "  "  << v4[1] << "  " << v4[2] << '\n';
//  std::cout.flush();
//  std::cout << v4[3];
  std::cout << cross(v1, v2) << '\n';
  std::cout << v1 << '\n';  
  v1 += v2;
  std::cout << v1 << '\n';
  v1 -= v2;
  std::cout << v1 << '\n';
  {
    Vector v1(1.0, 2.0, 3.0);  // some more examples, but I was too lazy to
    Vector v2 = v1.unit();     // fix the names
    std::cout << v2 << '\n';
    std::cout << v2.dot(v1) << '\n';
    std::cout << v1.length() << '\n';
    v1 *= -1;
    v1 += v2*14;
    std::cout << v1 << '\n';
  }
}
#endif

#endif

