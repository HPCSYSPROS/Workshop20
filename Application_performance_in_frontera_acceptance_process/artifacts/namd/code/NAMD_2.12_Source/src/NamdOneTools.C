/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   NAMD 1.X functions copied without change for use primarily in WorkDistrib
*/

#include "InfoStream.h"
#include "common.h"
#include "NamdTypes.h"
#include "NamdOneTools.h"
#include "Vector.h"
#include "PDB.h"
#include "Molecule.h"
#define MIN_DEBUG_LEVEL 4
//#define DEBUGM
#include "Debug.h"

//#define DEBUGM
#include "Debug.h"

/************************************************************************/
/*   FUNCTION read_binary_coors						*/
/*   INPUTS:								*/
/*	fname - Filename to read coordinates from			*/
/*	pdbobj - PDB object to place coordinates into			*/
/*	This function reads initial coordinates from a binary 		*/
/*	restart file							*/
/************************************************************************/

void read_binary_coors(char *fname, PDB *pdbobj) {
  Vector *newcoords;	//  Array of vectors to hold coordinates from file

  //  Allocate an array to hold the new coordinates
  newcoords = new Vector[pdbobj->num_atoms()];

  //  Read the coordinate from the file
  read_binary_file(fname,newcoords,pdbobj->num_atoms());

  //  Set the coordinates in the PDB object to the new coordinates
  pdbobj->set_all_positions(newcoords);

  //  Clean up
  delete [] newcoords;

} // END OF FUNCTION read_binary_coors()


void read_binary_file(const char *fname, Vector *data, int n)
{
  int32 filen;          //  Number of atoms read from file
  FILE *fp;             //  File descriptor
  int needToFlip = 0;

  iout << iINFO << "Reading from binary file " << fname << "\n" << endi;

  //  Open the file and die if the open fails
  if ( (fp = Fopen(fname, "rb")) == NULL)
  {
    char errmsg[256];
    sprintf(errmsg, "Unable to open binary file %s", fname);
    NAMD_die(errmsg);
  }

  //  read the number of coordinates in this file
  if (fread(&filen, sizeof(int32), 1, fp) != (size_t)1)
  {
    char errmsg[256];
    sprintf(errmsg, "Error reading binary file %s", fname);
    NAMD_die(errmsg);
  }

  //  read the number of coordinates in this file
  //  check for palindromic number of atoms
  char lenbuf[4];
  memcpy(lenbuf, (const char *)&filen, 4);
  char tmpc;
  tmpc = lenbuf[0]; lenbuf[0] = lenbuf[3]; lenbuf[3] = tmpc;
  tmpc = lenbuf[1]; lenbuf[1] = lenbuf[2]; lenbuf[2] = tmpc;
  if ( ! memcmp((const char *)&filen, lenbuf, 4) ) {
    iout << iWARN << "Number of atoms in binary file " << fname <<
		" is palindromic, assuming same endian.\n" << endi;
  }

  //  Die if this doesn't match the number in our system
  if (filen != n)
  {
    needToFlip = 1;
    memcpy((char *)&filen, lenbuf, 4);
  }
  if (filen != n)
  {
    char errmsg[256];
    sprintf(errmsg, "Incorrect atom count in binary file %s", fname);
    NAMD_die(errmsg);
  }

  if (fread(data, sizeof(Vector), n, fp) != (size_t)n)
  {
    char errmsg[256];
    sprintf(errmsg, "Error reading binary file %s", fname);
    NAMD_die(errmsg);
  }

  Fclose(fp);

  if (needToFlip) { 
    iout << iWARN << "Converting binary file " << fname << "\n" << endi;
    int i;
    char *cdata = (char *) data;
    for ( i=0; i<3*n; ++i, cdata+=8 ) {
      char tmp0, tmp1, tmp2, tmp3;
      tmp0 = cdata[0]; tmp1 = cdata[1];
      tmp2 = cdata[2]; tmp3 = cdata[3];
      cdata[0] = cdata[7]; cdata[1] = cdata[6];
      cdata[2] = cdata[5]; cdata[3] = cdata[4];
      cdata[7] = tmp0; cdata[6] = tmp1;
      cdata[5] = tmp2; cdata[4] = tmp3;
    }
  }

}


/* 
 * Generate a 3x3 transformation matrix for rotation about a given
 * vector v by an angle given in degrees
 */
void vec_rotation_matrix( BigReal angle, Vector v, BigReal m[] ) {
   /* This function contributed by Erich Boleyn (erich@uruk.org) */
   BigReal mag, s, c;
   BigReal xx, yy, zz, xy, yz, zx, xs, ys, zs, one_c;

   s = sin(angle * PI/180.0);
   c = cos(angle * PI/180.0);

   mag = v.length();

   if (mag == 0.0) {
      /* generate an identity matrix and return */
      for ( int i = 0; i < 9; ++i ) m[i] = 0.0;
      m[0] = m[4] = m[8] = 1.0;
      return;
   }

   // normalize the vector 
   v /= mag;

   /*
    *     Arbitrary axis rotation matrix.
    *
    *  This is composed of 5 matrices, Rz, Ry, T, Ry', Rz', multiplied
    *  like so:  Rz * Ry * T * Ry' * Rz'.  T is the final rotation
    *  (which is about the X-axis), and the two composite transforms
    *  Ry' * Rz' and Rz * Ry are (respectively) the rotations necessary
    *  from the arbitrary axis to the X-axis then back.  They are
    *  all elementary rotations.
    *
    *  Rz' is a rotation about the Z-axis, to bring the axis vector
    *  into the x-z plane.  Then Ry' is applied, rotating about the
    *  Y-axis to bring the axis vector parallel with the X-axis.  The
    *  rotation about the X-axis is then performed.  Ry and Rz are
    *  simply the respective inverse transforms to bring the arbitrary
    *  axis back to it's original orientation.  The first transforms
    *  Rz' and Ry' are considered inverses, since the data from the
    *  arbitrary axis gives you info on how to get to it, not how
    *  to get away from it, and an inverse must be applied.
    *
    *  The basic calculation used is to recognize that the arbitrary
    *  axis vector (x, y, z), since it is of unit length, actually
    *  represents the sines and cosines of the angles to rotate the
    *  X-axis to the same orientation, with theta being the angle about
    *  Z and phi the angle about Y (in the order described above)
    *  as follows:
    *
    *  cos ( theta ) = x / sqrt ( 1 - z^2 )
    *  sin ( theta ) = y / sqrt ( 1 - z^2 )
    *
    *  cos ( phi ) = sqrt ( 1 - z^2 )
    *  sin ( phi ) = z
    *
    *  Note that cos ( phi ) can further be inserted to the above
    *  formulas:
    *
    *  cos ( theta ) = x / cos ( phi )
    *  sin ( theta ) = y / sin ( phi )
    *
    *  ...etc.  Because of those relations and the standard trigonometric
    *  relations, it is pssible to reduce the transforms down to what
    *  is used below.  It may be that any primary axis chosen will give the
    *  same results (modulo a sign convention) using thie method.
    *
    *  Particularly nice is to notice that all divisions that might
    *  have caused trouble when parallel to certain planes or
    *  axis go away with care paid to reducing the expressions.
    *  After checking, it does perform correctly under all cases, since
    *  in all the cases of division where the denominator would have
    *  been zero, the numerator would have been zero as well, giving
    *  the expected result.
    */

   // store matrix in (row, col) form, i.e., m(row,col) = m[row*3+col]

   xx = v.x * v.x;
   yy = v.y * v.y;
   zz = v.z * v.z;
   xy = v.x * v.y;
   yz = v.y * v.z;
   zx = v.z * v.x;
   xs = v.x * s;
   ys = v.y * s;
   zs = v.z * s;
   one_c = 1.0 - c;

   m[0] = (one_c * xx) + c;
   m[1] = (one_c * xy) - zs;
   m[2] = (one_c * zx) + ys;
   
   m[3] = (one_c * xy) + zs;
   m[4] = (one_c * yy) + c;
   m[5] = (one_c * yz) - xs;
   
   m[6] = (one_c * zx) - ys;
   m[7] = (one_c * yz) + xs;
   m[8] = (one_c * zz) + c;
}

// multiply vector v by a 3x3 matrix m stored so that m(row,col) = m[row*3+col]
Vector mat_multiply_vec(const Vector &v, BigReal m[]) {
      return Vector( m[0]*v.x + m[1]*v.y + m[2]*v.z,
		     m[3]*v.x + m[4]*v.y + m[5]*v.z,
		     m[6]*v.x + m[7]*v.y + m[8]*v.z );
}

