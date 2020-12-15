/*@@
  @header   ADM_Derivative.h
  @date     June 2002
  @author   Denis Pollney
  @desc
            Derivative operators.
  @enddesc
@@*/

#ifndef ADM_DERIVATIVE_H
#define ADM_DERIVATIVE_H

#include "ADM_Spacing.h"

/*
 *   2nd order operators.
 */ 

#define ADM_DX_2(var,i,j,k)   I2DX*(var(i+1,j,k) - var(i-1,j,k))
#define ADM_DY_2(var,i,j,k)   I2DY*(var(i,j+1,k) - var(i,j-1,k))
#define ADM_DZ_2(var,i,j,k)   I2DZ*(var(i,j,k+1) - var(i,j,k-1))

#define ADM_DXX_2(var,i,j,k)  IDXX*(var(i+1,j,k) - 2.d0*var(i,j,k) \
                                    + var(i-1,j,k))
#define ADM_DYY_2(var,i,j,k)  IDYY*(var(i,j+1,k) - 2.d0*var(i,j,k) \
                                    + var(i,j-1,k))
#define ADM_DZZ_2(var,i,j,k)  IDZZ*(var(i,j,k+1) - 2.d0*var(i,j,k) \
                                    + var(i,j,k-1))

#define ADM_DXY_2(var,i,j,k)  IDXY*(var(i+1,j+1,k) - var(i+1,j-1,k) \
                                    - var(i-1,j+1,k) + var(i-1,j-1,k))
#define ADM_DXZ_2(var,i,j,k)  IDXZ*(var(i+1,j,k+1) - var(i+1,j,k-1) \
                                    - var(i-1,j,k+1) + var(i-1,j,k-1))
#define ADM_DYZ_2(var,i,j,k)  IDYZ*(var(i,j+1,k+1) - var(i,j+1,k-1) \
                                    - var(i,j-1,k+1) + var(i,j-1,k-1))
/*
 *   4th order operators.
 */

#define ADM_DX_4(var,i,j,k)   I12DX*(-var(i+2,j,k) + var(i-2,j,k) \
                                     + 8.d0*(var(i+1,j,k) - var(i-1,j,k)))
#define ADM_DY_4(var,i,j,k)   I12DY*(-var(i,j+2,k) + var(i,j-2,k) \
                                     + 8.d0*(var(i,j+1,k) - var(i,j-1,k)))
#define ADM_DZ_4(var,i,j,k)   I12DZ*(-var(i,j,k+2) + var(i,j,k-2) \
                                     + 8.d0*(var(i,j,k+1) - var(i,j,k-1)))

#define ADM_DXX_4(var,i,j,k)  I12DXX*(-var(i+2,j,k) - var(i-2,j,k) \
                                      + 16.d0*(var(i+1,j,k) + var(i-1,j,k)) \
                                      - 30.d0*var(i,j,k))
#define ADM_DYY_4(var,i,j,k)  I12DYY*(-var(i,j+2,k) - var(i,j-2,k) \
                                      + 16.d0*(var(i,j+1,k) + var(i,j-1,k)) \
                                      - 30.d0*var(i,j,k))
#define ADM_DZZ_4(var,i,j,k)  I12DZZ*(-var(i,j,k+2) - var(i,j,k-2) \
                                      + 16.d0*(var(i,j,k+1) + var(i,j,k-1)) \
                                      - 30.d0*var(i,j,k))

#define ADM_DXY_4(var,i,j,k)  I36DXY* \
  (var(i+2,j+2,k) - var(i+2,j-2,k) - var(i-2,j+2,k) + var(i-2,j-2,k) \
 + 8.d0*(-var(i+2,j+1,k) + var(i+2,j-1,k) - var(i+1,j+2,k) + var(i+1,j-2,k) \
         + var(i-2,j+1,k) - var(i-2,j-1,k) - var(i-1,j-2,k) + var(i-1,j+2,k)) \
 + 64.d0*(var(i+1,j+1,k) - var(i+1,j-1,k) - var(i-1,j+1,k) + var(i-1,j-1,k)))
#define ADM_DXZ_4(var,i,j,k)  I36DXZ* \
  (var(i+2,j,k+2) - var(i+2,j,k-2) - var(i-2,j,k+2) + var(i-2,j,k-2) \
 + 8.d0*(-var(i+2,j,k+1) + var(i+2,j,k-1) - var(i+1,j,k+2) + var(i+1,j,k-2) \
         + var(i-2,j,k+1) - var(i-2,j,k-1) - var(i-1,j,k-2) + var(i-1,j,k+2)) \
 + 64.d0*(var(i+1,j,k+1) - var(i+1,j,k-1) - var(i-1,j,k+1) + var(i-1,j,k-1)))
#define ADM_DYZ_4(var,i,j,k)  I36DYZ* \
  (var(i,j+2,k+2) - var(i,j+2,k-2) - var(i,j-2,k+2) + var(i,j-2,k-2) \
 + 8.d0*(-var(i,j+2,k+1) + var(i,j+2,k-1) - var(i,j+1,k+2) + var(i,j+1,k-2) \
         + var(i,j-2,k+1) - var(i,j-2,k-1) - var(i,j-1,k-2) + var(i,j-1,k+2)) \
 + 64.d0*(var(i,j+1,k+1) - var(i,j+1,k-1) - var(i,j-1,k+1) + var(i,j-1,k-1)))


#define ADM_ADV_DX_P1(var,i,j,k) betax(i,j,k)*IDX*(var(i+1,j,k) - var(i,j,k))
#define ADM_ADV_DY_P1(var,i,j,k) betay(i,j,k)*IDY*(var(i,j+1,k) - var(i,j,k))
#define ADM_ADV_DZ_P1(var,i,j,k) betaz(i,j,k)*IDZ*(var(i,j,k+1) - var(i,j,k))

#define ADM_ADV_DX_M1(var,i,j,k) betax(i,j,k)*IDX*(var(i,j,k) - var(i-1,j,k))
#define ADM_ADV_DY_M1(var,i,j,k) betay(i,j,k)*IDY*(var(i,j,k) - var(i,j-1,k))
#define ADM_ADV_DZ_M1(var,i,j,k) betaz(i,j,k)*IDZ*(var(i,j,k) - var(i,j,k-1))

#define ADM_ADV_DX_P2(var,i,j,k) betax(i,j,k)*I2DX*(-3.d0*var(i,j,k) \
                                            + 4.d0*var(i+1,j,k) - var(i+2,j,k))
#define ADM_ADV_DY_P2(var,i,j,k) betay(i,j,k)*I2DY*(-3.d0*var(i,j,k) \
                                            + 4.d0*var(i,j+1,k) - var(i,j+2,k))
#define ADM_ADV_DZ_P2(var,i,j,k) betaz(i,j,k)*I2DZ*(-3.d0*var(i,j,k) \
                                            + 4.d0*var(i,j,k+1) - var(i,j,k+2))

#define ADM_ADV_DX_M2(var,i,j,k) betax(i,j,k)*I2DX*(3.d0*var(i,j,k) \
                                           - 4.d0*var(i-1,j,k) + var(i-2,j,k))
#define ADM_ADV_DY_M2(var,i,j,k) betay(i,j,k)*I2DY*(3.d0*var(i,j,k) \
                                           - 4.d0*var(i,j-1,k) + var(i,j-2,k))
#define ADM_ADV_DZ_M2(var,i,j,k) betaz(i,j,k)*I2DZ*(3.d0*var(i,j,k) \
                                           - 4.d0*var(i,j,k-1) + var(i,j,k-2))

#endif
