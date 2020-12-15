#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "ExternalVariables.h"



/* Coefficients taken from
 * http://en.wikipedia.org/wiki/Linear_multistep_method, which cites
 * (Hairer, Nørsett & Wanner 1993, §III.1; Butcher 2003, p. 103).
 * 
 * Hairer, Ernst; Nørsett, Syvert Paul; Wanner, Gerhard (1993),
 * Solving ordinary differential equations I: Nonstiff problems (2nd
 * ed.), Berlin: Springer Verlag, ISBN 978-3-540-56670-0.
 *
 * Butcher, John C. (2003), Numerical Methods for Ordinary
 * Differential Equations, John Wiley, ISBN 978-0-471-96758-3.
 */

/* The following Mathematic expression (see
 * http://en.wikipedia.org/wiki/Linear_multistep_method) calculates
 * Adams-Bashforth coefficients for arbitrary orders:
 *
 * b = Table[
 *   Table[(-1)^j/(j! (s - j - 1)!) Integrate[
 *      Product[If[j == i, 1, u + i], {i, 0, s - 1}], {u, 0, 1}], {j, 0, 
 *     s - 1}], {s, 1, 8}]
 *
 * where s is the desired order.  The coefficients up to order 8 are:
 * 
 * {1}
 *
 * {3/2, -(1/2)}
 *
 * {23/12, -(4/3), 5/12}
 *
 * {55/24, -(59/24), 37/24, -(3/8)}
 *
 * {1901/720, -(1387/360), 109/30, -(637/360), 251/720}
 *
 * {4277/1440, -(2641/480), 4991/720, -(3649/720), 959/480, -(95/288)}
 *
 * {198721/60480, -(18637/2520), 235183/20160, -(10754/945),
 * 135713/20160, -(5603/2520), 19087/60480}
 *
 * {16083/4480, -(1152169/120960), 242653/13440, -(296053/13440),
 * 2102243/120960, -(115747/13440), 32863/13440, -(5257/17280)}
 */



static void order1 (CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
#pragma omp parallel for
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] += dt * RHSVar0[index];
  }
}

static void order2 (CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
  CCTK_REAL const* restrict const RHSVar1 = RHSVars[1];
  CCTK_REAL const f0 = + (3.0/2.0) * dt;
  CCTK_REAL const f1 = - (1.0/2.0) * dt;
#pragma omp parallel for
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] += f0 * RHSVar0[index] + f1 * RHSVar1[index];
  }
}

static void order3 (CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
  CCTK_REAL const* restrict const RHSVar1 = RHSVars[1];
  CCTK_REAL const* restrict const RHSVar2 = RHSVars[2];
  CCTK_REAL const f0 = + (23.0/12.0) * dt;
  CCTK_REAL const f1 = - ( 4.0/ 3.0) * dt;
  CCTK_REAL const f2 = + ( 5.0/12.0) * dt;
#pragma omp parallel for
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] +=
      f0 * RHSVar0[index] + f1 * RHSVar1[index] + f2 * RHSVar2[index];
  }
}

static void order4 (CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
  CCTK_REAL const* restrict const RHSVar1 = RHSVars[1];
  CCTK_REAL const* restrict const RHSVar2 = RHSVars[2];
  CCTK_REAL const* restrict const RHSVar3 = RHSVars[3];
  CCTK_REAL const f0 = + (55.0/24.0) * dt;
  CCTK_REAL const f1 = - (59.0/24.0) * dt;
  CCTK_REAL const f2 = + (37.0/24.0) * dt;
  CCTK_REAL const f3 = - ( 3.0/ 8.0) * dt;
#pragma omp parallel for
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] +=
      f0 * RHSVar0[index] + f1 * RHSVar1[index] + f2 * RHSVar2[index] +
      f3 * RHSVar3[index];
  }
}

static void order5 (CCTK_REAL* restrict const UpdateVar,
                    CCTK_REAL const* restrict *restrict const RHSVars,
                    CCTK_REAL const dt,
                    int const totalsize)
{
  CCTK_REAL const* restrict const RHSVar0 = RHSVars[0];
  CCTK_REAL const* restrict const RHSVar1 = RHSVars[1];
  CCTK_REAL const* restrict const RHSVar2 = RHSVars[2];
  CCTK_REAL const* restrict const RHSVar3 = RHSVars[3];
  CCTK_REAL const* restrict const RHSVar4 = RHSVars[4];
  CCTK_REAL const f0 = + (1901.0/720.0) * dt;
  CCTK_REAL const f1 = - (1387.0/360.0) * dt;
  CCTK_REAL const f2 = + ( 109.0/ 30.0) * dt;
  CCTK_REAL const f3 = - ( 637.0/360.0) * dt;
  CCTK_REAL const f4 = + ( 251.0/720.0) * dt;
#pragma omp parallel for
  for (int index = 0; index < totalsize; index++) {
    UpdateVar[index] +=
      f0 * RHSVar0[index] + f1 * RHSVar1[index] + f2 * RHSVar2[index] +
      f3 * RHSVar3[index] + f4 * RHSVar4[index];
  }
}

/* Array of function pointers */
static
void (* const orders[]) (CCTK_REAL* restrict const UpdateVar,
                         CCTK_REAL const* restrict *restrict const RHSVars,
                         CCTK_REAL const dt,
                         int const totalsize)
= { order1, order2, order3, order4, order5 };
static int const max_order = sizeof orders / sizeof *orders;



void MoL_ABAdd(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_REAL const dt = CCTK_DELTA_TIME;
  
  /* Determine the order of accuracy */
  int order;
  if (CCTK_EQUALS(AB_Type,"1")) {
    order = 1 ;
  } else if (CCTK_EQUALS(AB_Type,"2")) {
    order = 2;
  } else if (CCTK_EQUALS(AB_Type,"3")) {
    order = 3;
  } else if (CCTK_EQUALS(AB_Type,"4")) {
    order = 4;
  } else if (CCTK_EQUALS(AB_Type,"5")) {
    order = 5;
  } else {
    abort();
  }
  if (AB_initially_reduce_order) {
    /* Reduce the order for the first time steps */
    /* int const iteration = cctk_iteration; */
    int const iteration = 1 + lrint((cctk_time - cctk_initial_time) / dt);
    if (order > iteration) {
      order = iteration;
      CCTK_VInfo (CCTK_THORNSTRING,
                  "Reducing Adams-Bashforth order to %d", order);
    }
  }
/* printf ("MoL AB: iter=%d, order=%d\n", cctk_iteration, order); */
  assert (order >= 1 && order <= max_order);
  
  int totalsize = 1;
  for (int arraydim = 0; arraydim < cctk_dim; arraydim++) {
    totalsize *= cctk_ash[arraydim];
  }
  
  for (int var = 0; var < MoLNumEvolvedVariables; var++) {
    CCTK_REAL* restrict const UpdateVar =
      CCTK_VarDataPtrI(cctkGH, 0, EvolvedVariableIndex[var]);
    CCTK_REAL const* restrict RHSVars[order];
    for (int tl = 0; tl < order; tl++) {
      RHSVars[tl] = CCTK_VarDataPtrI(cctkGH, tl, RHSVariableIndex[var]);
    }
    
    /* Add RHS */
    (orders[order-1]) (UpdateVar, RHSVars, dt, totalsize);
  } /* var */
  
  for (int var = 0; var < MoLNumEvolvedArrayVariables; var++) {
    CCTK_REAL* restrict const UpdateVar =
      CCTK_VarDataPtrI(cctkGH, 0, EvolvedArrayVariableIndex[var]);
    CCTK_REAL const* restrict RHSVars[order];
    for (int tl = 0; tl < order; tl++) {
      RHSVars[tl] = CCTK_VarDataPtrI(cctkGH, tl, RHSArrayVariableIndex[var]);
    }
    
    int const groupindex =
      CCTK_GroupIndexFromVarI(EvolvedArrayVariableIndex[var]);
    cGroupDynamicData arraydata;
    int const ierr = CCTK_GroupDynamicData(cctkGH, groupindex, &arraydata);
    if (ierr) {
      CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING, 
                  "The driver does not return group information for group '%s'.",
                  CCTK_GroupName(groupindex));
    }
    int arraytotalsize = 1;
    for (int arraydim = 0; arraydim < arraydata.dim; arraydim++) {
      arraytotalsize *= arraydata.ash[arraydim];
    }
    
    /* Add RHS */
    (orders[order-1]) (UpdateVar, RHSVars, dt, arraytotalsize);
  } /* var */

}
