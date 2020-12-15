#ifndef MATRIX4SYMMETRY_H
#define MATRIX4SYMMETRY_H

#include "NamdTypes.h"
class Matrix4Symmetry {
public:
  BigReal mat[16];
  Matrix4Symmetry();
  Matrix4Symmetry(const BigReal *);
  Matrix4Symmetry(BigReal []);

  void multpoint(BigReal point[3]) const;
  void identity();
  void transpose();
  void multmatrix(const Matrix4Symmetry &);
  void translate(BigReal x, BigReal y, BigReal z);
  void translate(BigReal d[3]);
};
#endif
