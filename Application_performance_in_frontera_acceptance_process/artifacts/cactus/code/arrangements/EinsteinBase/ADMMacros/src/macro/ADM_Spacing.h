/*@@
  @header   ADM_Spacing.h
  @date     June 2002
  @author   Denis Pollney
  @desc
            Calculate various spacing dependent scalars.
  @enddesc
@@*/

#ifndef ADM_SPACING_H
#define ADM_SPACING_H

#ifdef FCODE

       dt=CCTK_DELTA_TIME
       dx=CCTK_DELTA_SPACE(1)
       dy=CCTK_DELTA_SPACE(2)
       dz=CCTK_DELTA_SPACE(3)

       idx=1.d0/dx
       idy=1.d0/dy
       idz=1.d0/dz

       i2dx=idx/2.d0
       i2dy=idy/2.d0
       i2dz=idz/2.d0

       i12dx=idx/12.d0
       i12dy=idy/12.d0
       i12dz=idz/12.d0

       idxx=idx**2
       idyy=idy**2
       idzz=idz**2

       i12dxx=idxx/12.d0
       i12dyy=idyy/12.d0
       i12dzz=idzz/12.d0

       idxy=i2dx*i2dy
       idxz=i2dx*i2dz
       idyz=i2dy*i2dz

       i36dxy=idxy/36.d0
       i36dxz=idxz/36.d0
       i36dyz=idyz/36.d0

#endif

#ifdef CCODE

       dt=CCTK_DELTA_TIME;
       dx=CCTK_DELTA_SPACE(0);
       dy=CCTK_DELTA_SPACE(1);
       dz=CCTK_DELTA_SPACE(2);

       idx=1.0/dx;
       idy=1.0/dy;
       idz=1.0/dz;

       i2dx=idx/2.0;
       i2dy=idy/2.0;
       i2dz=idz/2.0;

       i12dx=idx/12.0;
       i12dy=idy/12.0;
       i12dz=idz/12.0;

       idxx=idx*idx;
       idyy=idy*idy;
       idzz=idz*idz;

       i12dxx=idxx/12.0;
       i12dyy=idyy/12.0;
       i12dzz=idzz/12.0;

       idxy=i2dx*i2dy;
       idxz=i2dx*i2dz;
       idyz=i2dy*i2dz;

       i36dxy=idxy/36.0;
       i36dxz=idxz/36.0;
       i36dyz=idyz/36.0;

#endif

#endif
