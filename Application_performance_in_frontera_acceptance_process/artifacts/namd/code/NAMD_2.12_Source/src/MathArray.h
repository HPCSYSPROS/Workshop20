/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#ifndef MATHARRAY_H
#define MATHARRAY_H

#include "Array.h"

template <class Elem, int Size> class MathArray : public Array<Elem,Size> {

  public:
    // constructor
    MathArray(void) : Array<Elem,Size>(0) { ; }

    // destructor
    ~MathArray(void) { ; }

    // copy constructor
    MathArray(const Array<Elem,Size> &a2) : Array<Elem,Size>(a2) { ; }

    // assignment operator
    MathArray<Elem,Size> & operator= (const Array<Elem,Size> &a2) {
      Array<Elem,Size>::operator=(a2);
      return (*this);
    }

    // set all elements to a given value (like 0)
    MathArray(const Elem &v) : Array<Elem,Size>(v) { ; }
    MathArray<Elem,Size> & operator= (const Elem &v) {
      Array<Elem,Size>::operator=(v);
      return (*this);
    }

    // bulk mathematical operations
    MathArray<Elem,Size> & operator+= (const Array<Elem,Size> &a2) {
      for ( int i = 0; i < Size; ++i ) { this->data[i] += a2.data[i]; }
      return (*this);
    }
    MathArray<Elem,Size> & operator+= (const Elem &v) {
      for ( int i = 0; i < Size; ++i ) { this->data[i] += v; }
      return (*this);
    }
    MathArray<Elem,Size> & operator-= (const Array<Elem,Size> &a2) {
      for ( int i = 0; i < Size; ++i ) { this->data[i] -= a2.data[i]; }
      return (*this);
    }
    MathArray<Elem,Size> & operator-= (const Elem &v) {
      for ( int i = 0; i < Size; ++i ) { this->data[i] -= v; }
      return (*this);
    }
    MathArray<Elem,Size> & operator*= (const Array<Elem,Size> &a2) {
      for ( int i = 0; i < Size; ++i ) { this->data[i] *= a2.data[i]; }
      return (*this);
    }
    MathArray<Elem,Size> & operator*= (const Elem &v) {
      for ( int i = 0; i < Size; ++i ) { this->data[i] *= v; }
      return (*this);
    }
    MathArray<Elem,Size> & operator/= (const Array<Elem,Size> &a2) {
      for ( int i = 0; i < Size; ++i ) { this->data[i] /= a2.data[i]; }
      return (*this);
    }
    MathArray<Elem,Size> & operator/= (const Elem &v) {
      for ( int i = 0; i < Size; ++i ) { this->data[i] /= v; }
      return (*this);
    }
    friend MathArray<Elem,Size> operator+ (
		const Array<Elem,Size> &a1, const Array<Elem,Size> &a2 ) {
      return (MathArray<Elem,Size>(a1) += a2);
    }
    friend MathArray<Elem,Size> operator+ (
		const Array<Elem,Size> &a1, const Elem &v2 ) {
      return (MathArray<Elem,Size>(a1) += v2);
    }
    friend MathArray<Elem,Size> operator+ (
		const Elem & v1, const Array<Elem,Size> &a2 ) {
      return (MathArray<Elem,Size>(v1) += a2);
    }
    friend MathArray<Elem,Size> operator- (
		const Array<Elem,Size> &a1, const Array<Elem,Size> &a2 ) {
      return (MathArray<Elem,Size>(a1) -= a2);
    }
    friend MathArray<Elem,Size> operator- (
		const Array<Elem,Size> &a1, const Elem &v2 ) {
      return (MathArray<Elem,Size>(a1) -= v2);
    }
    friend MathArray<Elem,Size> operator- (
		const Elem & v1, const Array<Elem,Size> &a2 ) {
      return (MathArray<Elem,Size>(v1) -= a2);
    }
    friend MathArray<Elem,Size> operator* (
		const Array<Elem,Size> &a1, const Array<Elem,Size> &a2 ) {
      return (MathArray<Elem,Size>(a1) *= a2);
    }
    friend MathArray<Elem,Size> operator* (
		const Array<Elem,Size> &a1, const Elem &v2 ) {
      return (MathArray<Elem,Size>(a1) *= v2);
    }
    friend MathArray<Elem,Size> operator* (
		const Elem & v1, const Array<Elem,Size> &a2 ) {
      return (MathArray<Elem,Size>(v1) *= a2);
    }
    friend MathArray<Elem,Size> operator/ (
		const Array<Elem,Size> &a1, const Array<Elem,Size> &a2 ) {
      return (MathArray<Elem,Size>(a1) /= a2);
    }
    friend MathArray<Elem,Size> operator/ (
		const Array<Elem,Size> &a1, const Elem &v2 ) {
      return (MathArray<Elem,Size>(a1) /= v2);
    }
    friend MathArray<Elem,Size> operator/ (
		const Elem & v1, const Array<Elem,Size> &a2 ) {
      return (MathArray<Elem,Size>(v1) /= a2);
    }

};

#endif

