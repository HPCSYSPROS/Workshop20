/* $Header$ */

#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "util_ErrorCodes.h"
#include "cctk.h"
#include "../InterpLocalUniform.h"
#include "../common/structs.h"
#include "../common/load.h"
#include "../common/evaluate.h"
#include "../common/store.h"

/* function prototype */
#define FUNCTION_NAME			AEILocalInterp_U_Herm_3cube_3
#include "../template.h"

#define N_DIMS				3
#define MOLECULE_MIN_M			-2
#define MOLECULE_MAX_M			3
#define MOLECULE_SIZE			6

/* which derivative ops do we support? */
#define HAVE_OP_I
#define HAVE_OP_DX
#define HAVE_OP_DY
#define HAVE_OP_DZ
#define HAVE_OP_DXX
#define HAVE_OP_DXY
#define HAVE_OP_DXZ
#define HAVE_OP_DYY
#define HAVE_OP_DYZ
#define HAVE_OP_DZZ

#define XYZ				x, y, z
#define FP_XYZ				fp x, fp y, fp z
#define STRIDE_IJK			stride_i, stride_j, stride_k
#define JACOBIAN_MIJK_STRIDE		Jacobian_mi_stride, Jacobian_mj_stride, Jacobian_mk_stride

#define DATA_STRUCT			data_struct_3d_cube_size6
#define COEFFS_STRUCT			coeffs_struct_3d_cube_size6

#define LOAD_DATA_REAL			AEILocalInterp_load_3dcube6_r
#define LOAD_DATA_REAL4			AEILocalInterp_load_3dcube6_r4
#define LOAD_DATA_REAL8			AEILocalInterp_load_3dcube6_r8
#define LOAD_DATA_REAL16		AEILocalInterp_load_3dcube6_r16
#define LOAD_DATA_COMPLEX		AEILocalInterp_load_3dcube6_c
#define LOAD_DATA_COMPLEX8		AEILocalInterp_load_3dcube6_c8
#define LOAD_DATA_COMPLEX16		AEILocalInterp_load_3dcube6_c16
#define LOAD_DATA_COMPLEX32		AEILocalInterp_load_3dcube6_c32

#define EVALUATE_MOLECULE		AEILocalInterp_eval_3dcube6

#define STORE_COEFFS			AEILocalInterp_store_3dcube6

/* note pathnames are all relative to "../template.c" */
#define COEFFS_I_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-I.compute.c"
#define COEFFS_DX_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-dx.compute.c"
#define COEFFS_DY_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-dy.compute.c"
#define COEFFS_DZ_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-dz.compute.c"
#define COEFFS_DXX_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-dxx.compute.c"
#define COEFFS_DXY_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-dxy.compute.c"
#define COEFFS_DXZ_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-dxz.compute.c"
#define COEFFS_DYY_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-dyy.compute.c"
#define COEFFS_DYZ_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-dyz.compute.c"
#define COEFFS_DZZ_COMPUTE_FILE_NAME	"Hermite/3d.coeffs/3d.cube.order3/coeffs-dzz.compute.c"

/* actual code */
#include "../template.c"
