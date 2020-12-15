/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
   Function body for calc_fulldirect.
*/

{
  const BigReal coulomb = COULOMB * ComputeNonbondedUtil::scaling
				* ComputeNonbondedUtil::dielectric_1;
  BigReal *dp1 = data1;
  BigReal *rp1 = results1;
  int j_begin = 0;
  register BigReal electEnergy = 0.;
  register BigReal virial_xx = 0.;
  register BigReal virial_xy = 0.;
  register BigReal virial_xz = 0.;
  register BigReal virial_yy = 0.;
  register BigReal virial_yz = 0.;
  register BigReal virial_zz = 0.;

#ifdef FULLDIRECT_PERIODIC
  Vector a1 = lattice->a();
  Vector b1(0);  if ( lattice->a_p() ) b1 = lattice->a_r();
  Vector a2 = lattice->b();
  Vector b2(0);  if ( lattice->b_p() ) b2 = lattice->b_r();
  Vector a3 = lattice->c();
  Vector b3(0);  if ( lattice->c_p() ) b3 = lattice->c_r();
#endif

  for(int i=0; i<n1; ++i)
  {
    register BigReal p_i_x = *(dp1++);
    register BigReal p_i_y = *(dp1++);
    register BigReal p_i_z = *(dp1++);
    register BigReal kq_i = coulomb * *(dp1++);
    register BigReal f_i_x = 0.;
    register BigReal f_i_y = 0.;
    register BigReal f_i_z = 0.;
    if ( selfmode )
    {
      ++j_begin; data2 += 4; results2 += 3;
    }
    register BigReal *dp2 = data2;
    register BigReal *rp2 = results2;
    register int n2c = n2;
    register int j;
    for( j = j_begin; j<n2c; ++j)
    {
      register BigReal p_ij_x = p_i_x - *(dp2++);
      register BigReal p_ij_y = p_i_y - *(dp2++);
      register BigReal p_ij_z = p_i_z - *(dp2++);

#ifdef FULLDIRECT_PERIODIC
      Vector p_ij(p_ij_x,p_ij_y,p_ij_z);
      p_ij -= ( a1*floor(0.5+b1*p_ij) + a2*floor(0.5+b2*p_ij) + a3*floor(0.5+b3*p_ij) );
      p_ij_x = p_ij.x;
      p_ij_y = p_ij.y;
      p_ij_z = p_ij.z;
#endif

      register BigReal r_1;
      r_1 = 1./sqrt(p_ij_x * p_ij_x + p_ij_y * p_ij_y + p_ij_z * p_ij_z);
      register BigReal f = *(dp2++) * kq_i * r_1;
      electEnergy += f;
      f *= r_1*r_1;
      virial_xx += f * p_ij_x * p_ij_x;
      virial_xy += f * p_ij_x * p_ij_y;
      virial_xz += f * p_ij_x * p_ij_z;
      p_ij_x *= f;
      virial_yy += f * p_ij_y * p_ij_y;
      virial_yz += f * p_ij_y * p_ij_z;
      p_ij_y *= f;
      virial_zz += f * p_ij_z * p_ij_z;
      p_ij_z *= f;
      f_i_x += p_ij_x;
      f_i_y += p_ij_y;
      f_i_z += p_ij_z;
      *(rp2++) -= p_ij_x;
      *(rp2++) -= p_ij_y;
      *(rp2++) -= p_ij_z;
    }
    *(rp1++) += f_i_x;
    *(rp1++) += f_i_y;
    *(rp1++) += f_i_z;
  }

  virial.xx += virial_xx;
  virial.xy += virial_xy;
  virial.xz += virial_xz;
  virial.yx += virial_xy;
  virial.yy += virial_yy;
  virial.yz += virial_yz;
  virial.zx += virial_xz;
  virial.zy += virial_yz;
  virial.zz += virial_zz;
  return electEnergy;
}

