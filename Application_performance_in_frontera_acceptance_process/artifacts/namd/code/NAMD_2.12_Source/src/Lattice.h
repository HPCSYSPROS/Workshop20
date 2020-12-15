/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef LATTICE_H
#define LATTICE_H

#include <stdlib.h>
#include "NamdTypes.h"
#include <math.h>
#include "Tensor.h"

#ifdef POWERPC_TANINT
extern "builtin" double __tanint(double); //IEEE round
#define latticenearbyint(x)  __tanint(x)
#else
#define latticenearbyint(x)  floor((x)+0.5)
#endif

// Use latticenearbyint(X) instead of rint(X) because rint() is sensitive
// to the current rounding mode and floor() is not.  It's just safer.

typedef Vector ScaledPosition;

class Lattice
{
public:
  Lattice(void) : a1(0,0,0), a2(0,0,0), a3(0,0,0),
                  b1(0,0,0), b2(0,0,0), b3(0,0,0),
                  o(0,0,0), p1(0), p2(0), p3(0) {};

  // maps a transformation triplet onto a single integer
  static int index(int i=0, int j=0, int k=0)
  {
    return 9 * (k+1) + 3 * (j+1) + (i+1);
  }

  // sets lattice basis vectors but not origin (fixed center)
  void set(Vector A, Vector B, Vector C)
  {
    set(A,B,C,o);
  }

  // sets lattice basis vectors and origin (fixed center)
  void set(Vector A, Vector B, Vector C, Position Origin)
  {
    a1 = A; a2 = B; a3 = C; o = Origin;
    p1 = ( a1.length2() ? 1 : 0 );
    p2 = ( a2.length2() ? 1 : 0 );
    p3 = ( a3.length2() ? 1 : 0 );
    if ( ! p1 ) a1 = Vector(1.0,0.0,0.0);
    if ( ! p2 ) {
      Vector u1 = a1 / a1.length();
      Vector e_z(0.0,0.0,1.0);
      if ( fabs(e_z * u1) < 0.9 ) { a2 = cross(e_z,a1); }
      else { a2 = cross(Vector(1.0,0.0,0.0),a1); }
      a2 /= a2.length();
    }
    if ( ! p3 ) {
      a3 = cross(a1,a2);
      a3 /= a3.length();
    }
    if ( volume() < 0.0 ) a3 *= -1.0;
    recalculate();
  }

  // rescale lattice dimensions by factor, origin doesn't move
  void rescale(Tensor factor)
  {
    a1 = factor * a1;
    a2 = factor * a2;
    a3 = factor * a3;
    recalculate();
  }

  // rescale a position, keeping origin constant, assume 3D
  void rescale(Position &p, Tensor factor) const
  {
    p -= o;
    p = factor * p;
    p += o;
  }

  // transform scaled position to unscaled position
  Position unscale(ScaledPosition s) const
  {
    return (o + a1*s.x + a2*s.y + a3*s.z);
  }

  // transform unscaled position to scaled position
  ScaledPosition scale(Position p) const
  {
    p -= o;
    return Vector(b1*p,b2*p,b3*p);
  }

  // transforms a position nearest to a SCALED reference position
  Position nearest(Position data, ScaledPosition ref) const
  {
    ScaledPosition sn = scale(data);
    if ( p1 ) {
      sn.x -= latticenearbyint(sn.x - ref.x);
    }
    if ( p2 ) {
      sn.y -= latticenearbyint(sn.y - ref.y);
    }
    if ( p3 ) {
      sn.z -= latticenearbyint(sn.z - ref.z);
    }
    return unscale(sn);
  }

  // transforms a position nearest to a SCALED reference position
  // adds transform for later reversal
  Position nearest(Position data, ScaledPosition ref, Transform *t) const
  {
    ScaledPosition sn = scale(data);
    if ( p1 ) {
      BigReal tmp = sn.x - ref.x;
      BigReal rit = latticenearbyint(tmp);
      sn.x -= rit;
      t->i -= (int) rit;
    }
    if ( p2 ) {
      BigReal tmp = sn.y - ref.y;
      BigReal rit = latticenearbyint(tmp);
      sn.y -= rit;
      t->j -= (int) rit;
    }
    if ( p3 ) {
      BigReal tmp = sn.z - ref.z;
      BigReal rit = latticenearbyint(tmp);
      sn.z -= rit;
      t->k -= (int) rit;
    }
    return unscale(sn);
  }

  // applies stored transform to original coordinates
  Position apply_transform(Position data, const Transform &t) const
  {
    return ( data + t.i*a1 + t.j*a2 + t.k*a3 );
  }

  // reverses cumulative transformations for output
  Position reverse_transform(Position data, const Transform &t) const
  {
    return ( data - t.i*a1 - t.j*a2 - t.k*a3 );
  }

  // calculates shortest vector from p2 to p1 (equivalent to p1 - p2)
  Vector delta(const Position &pos1, const Position &pos2) const
  {
    Vector diff = pos1 - pos2;
#ifdef ARCH_POWERPC   //Prevents stack temporaries
    Vector result = diff;
    if ( p1 ) {
      BigReal fval = latticenearbyint(b1*diff); 
      result.x -= a1.x *fval;    
      result.y -= a1.y *fval;    
      result.z -= a1.z *fval;    
    }
    if ( p2 ) {
      BigReal fval = latticenearbyint(b2*diff);
      result.x -= a2.x * fval;
      result.y -= a2.y * fval;
      result.z -= a2.z * fval;
    }
    if ( p3 ) {
      BigReal fval = latticenearbyint(b3*diff);
      result.x -= a3.x * fval;
      result.y -= a3.y * fval;
      result.z -= a3.z * fval;
    }
    return result;
#else
    BigReal f1 = p1 ? latticenearbyint(b1*diff) : 0.;
    BigReal f2 = p2 ? latticenearbyint(b2*diff) : 0.;
    BigReal f3 = p3 ? latticenearbyint(b3*diff) : 0.;
    diff.x -= f1*a1.x + f2*a2.x + f3*a3.x;
    diff.y -= f1*a1.y + f2*a2.y + f3*a3.y;
    diff.z -= f1*a1.z + f2*a2.z + f3*a3.z;
    return diff;
#endif
  }

  // calculates shortest vector from origin to p1 (equivalent to p1 - o)
  Vector delta(const Position &pos1) const
  {
    Vector diff = pos1 - o;
    Vector result = diff;
    if ( p1 ) result -= a1*latticenearbyint(b1*diff);
    if ( p2 ) result -= a2*latticenearbyint(b2*diff);
    if ( p3 ) result -= a3*latticenearbyint(b3*diff);
    return result;
  }

  // calculates vector to bring p1 closest to origin
  Vector wrap_delta(const Position &pos1) const
  {
    Vector diff = pos1 - o;
    Vector result(0.,0.,0.);
    if ( p1 ) result -= a1*latticenearbyint(b1*diff);
    if ( p2 ) result -= a2*latticenearbyint(b2*diff);
    if ( p3 ) result -= a3*latticenearbyint(b3*diff);
    return result;
  }

  // calculates vector to bring p1 closest to origin in unscaled coordinates
  Vector wrap_nearest_delta(Position pos1) const
  {
    Vector diff = pos1 - o;
    Vector result0(0.,0.,0.);
    if ( p1 ) result0 -= a1*latticenearbyint(b1*diff);
    if ( p2 ) result0 -= a2*latticenearbyint(b2*diff);
    if ( p3 ) result0 -= a3*latticenearbyint(b3*diff);
    diff += result0;
    BigReal dist = diff.length2();
    Vector result(0.,0.,0.);
    for ( int i1=-p1; i1<=p1; ++i1 ) {
      for ( int i2 =-p2; i2<=p2; ++i2 ) {
        for ( int i3 =-p3; i3<=p3; ++i3 ) {
          Vector newresult = i1*a1+i2*a2+i3*a3;
          BigReal newdist = (diff+newresult).length2();
          if ( newdist < dist ) {
            dist = newdist;
            result = newresult;
          }
        }
      }
    }
    return result0 + result;
  }

  Vector offset(int i) const
  {
    return ( (i%3-1) * a1 + ((i/3)%3-1) * a2 + (i/9-1) * a3 );
  }

  static int offset_a(int i) { return (i%3-1); }
  static int offset_b(int i) { return ((i/3)%3-1); }
  static int offset_c(int i) { return (i/9-1); }

  // lattice vectors
  Vector a() const { return a1; }
  Vector b() const { return a2; }
  Vector c() const { return a3; }

  // only if along x y z axes
  int orthogonal() const {
    return ( ! ( a1.y || a1.z || a2.x || a2.z || a3.x || a3.y ) );
  }

  // origin (fixed center of cell)
  Vector origin() const
  {
    return o;
  }

  // reciprocal lattice vectors
  Vector a_r() const { return b1; }
  Vector b_r() const { return b2; }
  Vector c_r() const { return b3; }

  // periodic along this direction
  int a_p() const { return p1; }
  int b_p() const { return p2; }
  int c_p() const { return p3; }

  BigReal volume(void) const
  {
    return ( p1 && p2 && p3 ? cross(a1,a2) * a3 : 0.0 );
  }

private:
  Vector a1,a2,a3; // real lattice vectors
  Vector b1,b2,b3; // reciprocal lattice vectors (more or less)
  Vector o; // origin (fixed center of cell)
  int p1, p2, p3; // periodic along this lattice vector?

  // calculate reciprocal lattice vectors
  void recalculate(void) {
    {
      Vector c = cross(a2,a3);
      b1 = c / ( a1 * c );
    }
    {
      Vector c = cross(a3,a1);
      b2 = c / ( a2 * c );
    }
    {
      Vector c = cross(a1,a2);
      b3 = c / ( a3 * c );
    }
  }

};

#include <pup.h>
PUPbytes (Lattice);  

#endif

