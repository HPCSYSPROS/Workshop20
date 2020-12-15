// GRHydro_RegisterVars.cc
//
// converted from F90 to improve readability and maintainability
//
// Frank Loeffler

#include <cstdio>
#include <string>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

using namespace std;

extern "C" CCTK_INT GRHydro_UseGeneralCoordinates(const cGH * cctkGH);

// Utility functions to register variables with MoL
// Note: We could check for the return value here, but MoL issues a
//       level 0 warning in that case anyway. If that changes in the
//       future, a check can simply be inserted here.
static void register_evolved(string v1, string v2)
{
  DECLARE_CCTK_PARAMETERS;
  if (use_MoL_slow_multirate_sector)
     MoLRegisterEvolvedGroupSlow(CCTK_GroupIndex(v1.c_str()), CCTK_GroupIndex(v2.c_str()));
  else
     MoLRegisterEvolvedGroup(CCTK_GroupIndex(v1.c_str()), CCTK_GroupIndex(v2.c_str()));
}
static void register_constrained(string v1)
{
  MoLRegisterConstrainedGroup(CCTK_GroupIndex(v1.c_str()));
}
static void register_saveandrestore(string v1)
{
  MoLRegisterSaveAndRestoreGroup(CCTK_GroupIndex(v1.c_str()));
}

// Main function called by Cactus to register variables with MoL

extern "C" void GRHydro_Register(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int general_coordinates = GRHydro_UseGeneralCoordinates(cctkGH);

  // We need some aliased functions, so we first check if they are available
  string needed_funs[7] = {"MoLRegisterEvolvedGroup",
                           "MoLRegisterEvolvedGroupSlow",
                           "MoLRegisterConstrainedGroup",
                           "MoLRegisterSaveAndRestoreGroup",
                           "MoLRegisterEvolved",
                           "MoLRegisterEvolvedSlow",
                           "MoLRegisterConstrained"};
  for (int i = 0; i < 7; i++)
    if (!CCTK_IsFunctionAliased(needed_funs[i].c_str()))
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "The function \"%s\" has not been aliased!",
                 needed_funs[i].c_str());

  // Now we can set which variables have to be registered as which type with MoL
  register_constrained("HydroBase::rho");
  register_constrained("HydroBase::press");
  register_constrained("HydroBase::eps");
  register_constrained("HydroBase::w_lorentz");
  register_constrained("HydroBase::vel");
  if (*stress_energy_state)
  {
    register_constrained("TmunuBase::stress_energy_scalar");
    register_constrained("TmunuBase::stress_energy_vector");
    register_constrained("TmunuBase::stress_energy_tensor");
  }
  if (general_coordinates) {
    register_constrained("GRHydro::lvel");

    register_constrained("grhydro::local_shift");
    register_constrained("grhydro::local_metric");
    register_constrained("grhydro::local_extrinsic_curvature");
  }

  if (CCTK_EQUALS(evolution_method, "GRHydro"))
  {
    // dens and scon
    register_evolved("GRHydro::dens", "GRHydro::densrhs");
    register_evolved("GRHydro::scon", "GRHydro::srhs");

    if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) {
      register_constrained("HydroBase::Bvec");
      if (general_coordinates) {
        register_constrained("GRHydro::lBvec");
      }
      register_evolved("GRHydro::Bcons", "GRHydro::Bconsrhs");
      if(clean_divergence) {
        register_evolved("GRHydro::psidc" , "GRHydro::psidcrhs");
      }
    } else if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro_Avec")) {
      register_constrained("HydroBase::Bvec");
      register_evolved("HydroBase::Avec", "GRHydro::Avecrhs");
      if ( CCTK_EQUALS(Avec_gauge, "lorenz")) {
         register_evolved("HydroBase::Aphi", "GRHydro::Aphirhs");
      }
    }
    // entropycons
    if(CCTK_EQUALS(entropy_evolution_method,"GRHydro")){
      register_evolved("GRHydro::entropycons" , "GRHydro::entropyrhs");
    }

    // tau
    if (CCTK_EQUALS(GRHydro_eos_type, "General"))
      register_evolved("GRHydro::tau" , "GRHydro::taurhs");
    else if (CCTK_EQUALS(GRHydro_eos_type, "Polytype"))
      register_constrained("GRHydro::tau");
    else
      CCTK_WARN(0, "Don't recognize the type of EOS!");

    // lapse, metric, curv
    if(GRHydro_MaxNumSandRVars != 0) { // hack to save some memory since we "know" that someone else will register these as constrained
      register_saveandrestore("admbase::lapse");
      register_saveandrestore("admbase::metric");
      register_saveandrestore("admbase::curv");
    }

    // shift
    if (!CCTK_EQUALS(initial_shift, "none"))
    {
      if (CCTK_EQUALS(shift_evolution_method, "Comoving"))
      {
        register_constrained("admbase::shift");
        register_evolved("GRHydro::GRHydro_coords", 
                         "GRHydro::GRHydro_coords_rhs");
      }
      else if(GRHydro_MaxNumSandRVars != 0) // hack to save some memory since we "know" that someone else will register these as constrained
      { 
        register_saveandrestore("admbase::shift");
      }
    }

    // tracer
    if (evolve_tracer != 0)
      register_evolved("GRHydro::GRHydro_cons_tracers", 
                       "GRHydro::GRHydro_tracer_rhs");
 
    // note that we set in pararamcheck
    // evolve_y_e and evolve_temper, but MoLRegister
    // happens before paramcheck...
    if (CCTK_EQUALS(Y_e_evolution_method,"GRHydro")) {
      register_constrained("HydroBase::Y_e");
      register_evolved("GRHydro::Y_e_con", 
                       "GRHydro::Y_e_con_rhs");
      }

    if (CCTK_EQUALS(temperature_evolution_method,"GRHydro")) {
      register_constrained("HydroBase::temperature");
      register_constrained("HydroBase::entropy");
    }
    // particles
    if (number_of_particles > 0)
      register_evolved("GRHydro::particles", "GRHydro::particle_rhs");
  }
  else if (CCTK_EQUALS(evolution_method, "none"))
  {
    register_constrained("GRHydro::dens");
    register_constrained("GRHydro::scon");
    register_constrained("GRHydro::tau");

    if (CCTK_EQUALS(Bvec_evolution_method, "GRHydro")) {
      register_constrained("HydroBase::Bcons");
      if(clean_divergence) {
        register_constrained("GRHydro::psidc");
      }
    }
  }
}



void GRHydro_Set_Execution_Flags(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   
   const int index1 = CCTK_VarIndex("MoL::MoL_SlowStep");
   const int index2 = CCTK_VarIndex("MoL::MoL_SlowPostStep");
   
   if (index1 < 0 || index2 < 0)
   {
      CCTK_WARN(0, "Error: MoL does not provide MoL_SlowStep or MoL_SlowPostStep. Does your MoL thorn support multirate RK?");
   }
   
   *execute_MoL_Step =     *((CCTK_INT *) (CCTK_VarDataPtrI(cctkGH, 0, index1)));
   *execute_MoL_PostStep = *((CCTK_INT *) (CCTK_VarDataPtrI(cctkGH, 0, index2)));
}

void GRHydro_Reset_Execution_Flags(CCTK_ARGUMENTS)
{
   DECLARE_CCTK_ARGUMENTS;
   
   *execute_MoL_Step = 1;
   *execute_MoL_PostStep = 1;
}




