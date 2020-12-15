/*@@
   @file      ConfPhys.c
   @date      September 3rd 1999
   @author    Gabrielle Allen
   @desc
              Conversions between physical and conformal metrics.
              Be very careful using these functions, note that
              conformal_state is not changed, this is up to the
              calling routine
   @enddesc
   @version $Header$
 @@*/

#include "cctk.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusEinstein_StaticConformal_ConfPhys_c)

/********************************************************************
 *********************     Local Data Types   ***********************
 ********************************************************************/

/********************************************************************
 ********************* Local Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ***************** Scheduled Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 ********************* Other Routine Prototypes *********************
 ********************************************************************/

/********************************************************************
 *********************     Local Data   *****************************
 ********************************************************************/

/********************************************************************
 *********************     External Routines   **********************
 ********************************************************************/

/*@@
  @routine    StaticConf_ConfToPhysInPlace
  @date       September 3rd 1999
  @author     Gabrielle Allen
  @desc
  Convert the conformal metric to the physical in place.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory
  @var     nx
  @vdesc   number of points in the x direction
  @vtype   int
  @vio     in
  @var     ny
  @vdesc   number of points in the y direction
  @vtype   int
  @vio     in
  @var     nz
  @vdesc   number of points in the z direction
  @vtype   int
  @vio     in
  @var     psi
  @vdesc   conformal factor
  @vtype   const CCTK_REAL *
  @vio     in
  @var     gxx
  @vdesc   xx componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gxy
  @vdesc   xy componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gxz
  @vdesc   xz componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gyy
  @vdesc   yy componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gyz
  @vdesc   yz componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gzz
  @vdesc   zz componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout

@@*/
void StaticConf_ConfToPhysInPlace(
    CCTK_INT const nx, CCTK_INT const ny, CCTK_INT const nz,
    const CCTK_REAL *restrict const psi, CCTK_REAL *restrict const gxx,
    CCTK_REAL *restrict const gxy, CCTK_REAL *restrict const gxz,
    CCTK_REAL *restrict const gyy, CCTK_REAL *restrict const gyz,
    CCTK_REAL *restrict const gzz) {
  CCTK_WARN(4, "Converting metric: conformal -> physical");

  for (CCTK_INT index = 0; index < nx * ny * nz; ++index) {
    /* this should be faster than psi4 = pow (psi[index], 4) */
    CCTK_REAL psi4 = psi[index];
    psi4 = psi4 * psi4;
    psi4 = psi4 * psi4;

    gxx[index] *= psi4;
    gxy[index] *= psi4;
    gxz[index] *= psi4;
    gyy[index] *= psi4;
    gyz[index] *= psi4;
    gzz[index] *= psi4;
  }
}

/*@@
  @routine    StaticConf_PhysToConfInPlace
  @date       September 3rd 1999
  @author     Gabrielle Allen
  @desc
  Convert the physical metric to the conformal in place.
  @enddesc
  @calls
  @calledby
  @history

  @endhistory
  @var     nx
  @vdesc   number of points in the x direction
  @vtype   int
  @vio     in
  @var     ny
  @vdesc   number of points in the y direction
  @vtype   int
  @vio     in
  @var     nz
  @vdesc   number of points in the z direction
  @vtype   int
  @vio     in
  @var     psi
  @vdesc   conformal factor
  @vtype   const CCTK_REAL *
  @vio     in
  @var     gxx
  @vdesc   xx componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gxy
  @vdesc   xy componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gxz
  @vdesc   xz componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gyy
  @vdesc   yy componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gyz
  @vdesc   yz componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout
  @var     gzz
  @vdesc   zz componenent of metric
  @vtype   CCTK_REAL *
  @vio     inout

@@*/
void StaticConf_PhysToConfInPlace(
    CCTK_INT const nx, CCTK_INT const ny, CCTK_INT const nz,
    const CCTK_REAL *restrict const psi, CCTK_REAL *restrict const gxx,
    CCTK_REAL *restrict const gxy, CCTK_REAL *restrict const gxz,
    CCTK_REAL *restrict const gyy, CCTK_REAL *restrict const gyz,
    CCTK_REAL *restrict const gzz) {
  CCTK_WARN(4, "Converting metric: physical -> conformal");

  for (CCTK_INT index = 0; index < nx * ny * nz; ++index) {
    /* this should be faster than psi4 = pow (psi[index], 4) */
    CCTK_REAL psi4 = psi[index];
    psi4 = psi4 * psi4;
    psi4 = psi4 * psi4;

    /* get the reciprocal for turning divisions into multiplications */
    psi4 = 1.0 / psi4;

    gxx[index] *= psi4;
    gxy[index] *= psi4;
    gxz[index] *= psi4;
    gyy[index] *= psi4;
    gyz[index] *= psi4;
    gzz[index] *= psi4;
  }
}

/********************************************************************
 *********************     Local Routines   *************************
 ********************************************************************/
