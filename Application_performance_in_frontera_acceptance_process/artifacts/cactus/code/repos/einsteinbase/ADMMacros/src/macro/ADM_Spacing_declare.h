/*@@
  @header   ADM_Spacing_declare.h
  @date     June 2002
  @author   Denis Pollney
  @desc
            Declare various spacing dependent scalars.
  @enddesc
@@*/

#ifdef FCODE

      CCTK_REAL dt, dx, dy, dz
      CCTK_REAL idx, idy, idz
      CCTK_REAL i2dx, i2dy, i2dz
      CCTK_REAL i12dx, i12dy, i12dz
      CCTK_REAL idxx, idxy, idxz, idyy, idyz, idzz
      CCTK_REAL i12dxx, i12dyy, i12dzz
      CCTK_REAL i36dxy, i36dxz, i36dyz

#endif

#ifdef CCODE

      CCTK_REAL dt, dx, dy, dz;
      CCTK_REAL idx, idy, idz;
      CCTK_REAL i2dx, i2dy, i2dz;
      CCTK_REAL i12dx, i12dy, i12dz;
      CCTK_REAL idxx, idxy, idxz, idyy, idyz, idzz;
      CCTK_REAL i12dxx, i12dyy, i12dzz;
      CCTK_REAL i36dxy, i36dxz, i36dyz;

#endif
