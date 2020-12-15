#include "Matrix4Symmetry.h"

Matrix4Symmetry::Matrix4Symmetry() { identity(); }
Matrix4Symmetry::Matrix4Symmetry(const BigReal *m)  { memcpy(mat, m, 16*sizeof(BigReal)); }
Matrix4Symmetry::Matrix4Symmetry(BigReal arr []){
    for(int i = 0; i < 16; i++){mat[i] = arr[i];}
  }
void Matrix4Symmetry::multpoint(BigReal point[3]) const {
    BigReal tmp[3];
    BigReal itmp3 = 1.0f / (point[0]*mat[3] + point[1]*mat[7] +
                            point[2]*mat[11] + mat[15]);
    tmp[0] = itmp3*point[0];
    tmp[1] = itmp3*point[1];
    tmp[2] = itmp3*point[2];
    point[0]=tmp[0]*mat[0] + tmp[1]*mat[4] + tmp[2]*mat[ 8] + itmp3*mat[12];
    point[1]=tmp[0]*mat[1] + tmp[1]*mat[5] + tmp[2]*mat[ 9] + itmp3*mat[13];
    point[2]=tmp[0]*mat[2] + tmp[1]*mat[6] + tmp[2]*mat[10] + itmp3*mat[14];
  }

void Matrix4Symmetry::identity() {
    memset(mat, 0, 16*sizeof(BigReal));
    mat[0]=1.0f;
    mat[5]=1.0f;
    mat[10]=1.0f;
    mat[15]=1.0f;
  }
void Matrix4Symmetry::transpose() {
    BigReal tmp[16];
    int i,j;
    for(i=0;i<4;i++) {
      for(j=0;j<4;j++) {
        tmp[4*i+j] = mat[i+4*j];
      }
    }
    for(i=0;i<16;i++) mat[i] = tmp[i];
  }
  /// premultiply the matrix by the given matrix, this->other * this
void Matrix4Symmetry::multmatrix(const Matrix4Symmetry &m) {
    BigReal tmp[4];
    for (int j=0; j<4; j++) {
      tmp[0] = mat[j];
      tmp[1] = mat[4+j];
      tmp[2] = mat[8+j]; 
      tmp[3] = mat[12+j];
      for (int i=0; i<4; i++) {
        mat[4*i+j] = m.mat[4*i]*tmp[0] + m.mat[4*i+1]*tmp[1] +
          m.mat[4*i+2]*tmp[2] + m.mat[4*i+3]*tmp[3]; 
      }
    } 
  }
void Matrix4Symmetry::translate(BigReal x, BigReal y, BigReal z) {
    Matrix4Symmetry m;		
    m.mat[12] = x;
    m.mat[13] = y;
    m.mat[14] = z;
    multmatrix(m);
  }
void Matrix4Symmetry::translate(BigReal d[3]) { translate(d[0], d[1], d[2]); }
