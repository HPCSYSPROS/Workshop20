/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2008 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: Matrix4.h,v $
 *	$Author: jim $	$Locker:  $		$State: Exp $
 *	$Revision: 1.2 $	$Date: 2008/09/17 16:19:54 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * 4 x 4 Matrix, used for a transformation matrix.
 *
 ***************************************************************************/
#ifndef MATRIX_FOUR_H
#define MATRIX_FOUR_H

/// 4x4 matrix class with numerous operators, conversions, etc.
class Matrix4 {
public:
  Matrix4(void) { identity(); }                ///< identity constructor
  Matrix4(double f) { constant(f); }            ///< const elements constructor
  Matrix4(const double *m);                     ///< construct from double array
  Matrix4(const Matrix4& m) { loadmatrix(m); } ///< copy constructor 
  ~Matrix4(void) {}                            ///< destructor
  double mat[16];                               ///< the matrix itself

  /// multiplies a 3D point (first arg) by the Matrix, returns in second arg
  void multpoint3d (const double[3], double[3]) const;

  /// multiplies a 3D norm (first arg) by the Matrix, returns in second arg
  void multnorm3d (const double[3], double[3]) const;

  /// multiplies a 4D point (first arg) by the Matrix, returns in second arg
  void multpoint4d (const double[4], double[4]) const;

  /// clears the matrix (resets it to identity)
  void identity(void);
  
  /// sets the matrix so all items are the given constant value
  void constant(double);
  
  /// inverts the matrix, that is, 
  /// the inverse of the rotation, the inverse of the scaling, and 
  /// the opposite of the translation vector.
  /// returns 0 if there were no problems, -1 if the matrix is singular
  int inverse(void);
  
  /// transposes the matrix
  void transpose(void);
  
  /// replaces this matrix with the given one
  void loadmatrix(const Matrix4 &m);
  Matrix4& operator=(const Matrix4& m) {loadmatrix(m); return *this;}

  /// premultiply the matrix by the given matrix, this->other * this
  void multmatrix(const Matrix4 &);

  /// performs a left-handed rotation around an axis (char == 'x', 'y', or 'z')
  void rot(double, char); // angle in degrees

  /// apply a rotation around the given vector; angle in radians.
  void rotate_axis(const double axis[3], double angle);
  
  /// apply a rotation such that 'x' is brought along the given vector.
  void transvec(double x, double y, double z);
 
  /// apply a rotation such that the given vector is brought along 'x'.
  void transvecinv(double x, double y, double z);

  /// performs a translation
  void translate(double, double, double);
  void translate(double d[3]) { translate(d[0], d[1], d[2]); }

  /// performs scaling
  void scale(double, double, double);
  void scale(double f) { scale(f, f, f); }

  /// sets this matrix to represent a window perspective
  void window(double, double, double, double, double, double);

  /// sets this matrix to a 3D orthographic matrix
  void ortho(double, double, double, double, double, double);

  /// sets this matrix to a 2D orthographic matrix
  void ortho2(double, double, double, double);

  /// This subroutine defines a viewing transformation with the eye at point
  /// (vx,vy,vz) looking at the point (px,py,pz).  Twist is the right-hand
  /// rotation about this line.  The resultant matrix is multiplied with
  /// the top of the transformation stack and then replaces it.  Precisely,
  /// lookat does:
  /// lookat=trans(-vx,-vy,-vz)*rotate(theta,y)*rotate(phi,x)*rotate(-twist,z)
  void lookat(double, double, double, double, double, double, short);
};

/// Transform 3x3 into 4x4 matrix:
void trans_from_rotate(const double mat3[9], Matrix4 *mat4);

/// Print formatted matrix
void print_Matrix4(const Matrix4 *mat4);

#endif

