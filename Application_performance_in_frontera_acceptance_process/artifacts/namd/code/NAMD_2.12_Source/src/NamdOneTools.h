/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Extern's defined for NAMD 1.X functions used in odd places
   Mainly found in WorkDistrib.
*/

#include "common.h"
#include "NamdTypes.h"
#include "Vector.h"
#include "PDB.h"
#include "Molecule.h"

void read_binary_coors(char *fname, PDB *pdbobj);
void read_binary_file(const char *fname, Vector *vels, int n);
void vec_rotation_matrix(BigReal angle, Vector v, BigReal m[]);
Vector mat_multiply_vec(const Vector &v, BigReal m[]);

