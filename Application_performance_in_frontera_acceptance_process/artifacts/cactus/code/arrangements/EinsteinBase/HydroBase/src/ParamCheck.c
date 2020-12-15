#include <assert.h>
#include <stdlib.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "util_String.h"

static void callback(int idx, const char *optstring, void *callback_arg);

/** Ensure that all HydroBase initial data that are supposed to be read
    from a file are actually scheduled for the file reader.  */
void HydroBase_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const struct {
    const char *paramname;
    const char *paramvalue;
    const char *varname;
  } vars_to_read[] = {
      {"initial_hydro", initial_hydro, "HydroBase::rho"},
      {"initial_hydro", initial_hydro, "HydroBase::vel[0]"},
      {"initial_hydro", initial_hydro, "HydroBase::vel[1]"},
      {"initial_hydro", initial_hydro, "HydroBase::vel[2]"},
      {"initial_hydro", initial_hydro, "HydroBase::eps"},
      {"initial_hydro", initial_hydro, "HydroBase::press"},
      {"initial_Aphi", initial_Aphi, "HydroBase::Aphi"},
      {"initial_Avec", initial_Avec, "HydroBase::Avec[0]"},
      {"initial_Avec", initial_Avec, "HydroBase::Avec[1]"},
      {"initial_Avec", initial_Avec, "HydroBase::Avec[2]"},
      {"initial_Bvec", initial_Bvec, "HydroBase::Bvec[0]"},
      {"initial_Bvec", initial_Bvec, "HydroBase::Bvec[1]"},
      {"initial_Bvec", initial_Bvec, "HydroBase::Bvec[2]"},
      {"initial_Y_e", initial_Y_e, "HydroBase::Y_e"},
      {"initial_Abar", initial_Abar, "HydroBase::Abar"},
      {"initial_temperature", initial_temperature, "HydroBase::temperature"},
      {"initial_entropy", initial_entropy, "HydroBase::entropy"},
  };

  char *variable_is_read = malloc(CCTK_NumVars());
  assert(variable_is_read);
  for (int i = 0; i < CCTK_NumVars(); ++i) {
    variable_is_read[i] = 0;
  }

  const int nvars = CCTK_TraverseString(filereader_ID_vars, callback,
                                        variable_is_read, CCTK_GROUP_OR_VAR);
  assert(nvars >= 0);

  for (size_t i = 0; i < sizeof(vars_to_read) / sizeof(vars_to_read[0]); ++i) {
    if (CCTK_EQUALS(vars_to_read[i].paramvalue, "read from file")) {
      int const ivar = CCTK_VarIndex(vars_to_read[i].varname);
      assert(ivar >= 0);
      if (!variable_is_read[ivar]) {
        char *msg;
        size_t written = Util_asprintf(
            &msg, "'%s' is initialised using the file reader by '%s', but has "
                  "not been scheduled to be read.  Please set the parameter "
                  "\"IO::filereader_ID_vars\" accordingly.",
            vars_to_read[i].varname, vars_to_read[i].paramname);
        assert(written > 0);
        CCTK_PARAMWARN(msg);
        free(msg);
      }
    }
  }

  free(variable_is_read);
}

/** Mark a variable as to be read from the file reader.  */
static void callback(int idx, const char *optstring, void *callback_arg) {
  assert(idx >= 0 && idx < CCTK_NumVars());
  ((char *)callback_arg)[idx] = 1;
}
