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
#define FUNCTION_NAME			AEILocalInterp_U_LagMD_1cube_60
#include "../template.h"

#define N_DIMS				1
#define MOLECULE_MIN_M			-3
#define MOLECULE_MAX_M			3
#define MOLECULE_SIZE			7

/* which derivative ops do we support? */
#define HAVE_OP_I
#define HAVE_OP_DX
#define HAVE_OP_DXX

#define XYZ				x
#define FP_XYZ				fp x
#define STRIDE_IJK			stride_i
#define JACOBIAN_MIJK_STRIDE		Jacobian_mi_stride

#define DATA_STRUCT			data_struct_1d_cube_size7
#define COEFFS_STRUCT			coeffs_struct_1d_cube_size7

#define LOAD_DATA_REAL			AEILocalInterp_load_1dcube7_r
#define LOAD_DATA_REAL4			AEILocalInterp_load_1dcube7_r4
#define LOAD_DATA_REAL8			AEILocalInterp_load_1dcube7_r8
#define LOAD_DATA_REAL16		AEILocalInterp_load_1dcube7_r16
#define LOAD_DATA_COMPLEX		AEILocalInterp_load_1dcube7_c
#define LOAD_DATA_COMPLEX8		AEILocalInterp_load_1dcube7_c8
#define LOAD_DATA_COMPLEX16		AEILocalInterp_load_1dcube7_c16
#define LOAD_DATA_COMPLEX32		AEILocalInterp_load_1dcube7_c32

#define EVALUATE_MOLECULE		AEILocalInterp_eval_1dcube7

#define STORE_COEFFS			AEILocalInterp_store_1dcube7

/* note pathnames are all relative to "../template.c" */
#define COEFFS_I_COMPUTE_FILE_NAME	"Lagrange-maximum-degree/1d.coeffs/1d.cube.order6.smooth0/coeffs-I.compute.c"
#define COEFFS_DX_COMPUTE_FILE_NAME	"Lagrange-maximum-degree/1d.coeffs/1d.cube.order6.smooth0/coeffs-dx.compute.c"
#define COEFFS_DXX_COMPUTE_FILE_NAME	"Lagrange-maximum-degree/1d.coeffs/1d.cube.order6.smooth0/coeffs-dxx.compute.c"

/* actual code */
#include "../template.c"
