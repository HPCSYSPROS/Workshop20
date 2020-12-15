#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "HydroBase.h"

void HydroBase_Zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const np = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];

#pragma omp parallel for
  for (int i = 0; i < np; ++i) {
    rho[i] = 0.0;
    vel[i] = 0.0;
    vel[i + np] = 0.0;
    vel[i + 2 * np] = 0.0;
    w_lorentz[i] = 1.0;
    eps[i] = 0.0;
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Abar") >= 1) {
#pragma omp parallel for
    for (int i = 0; i < np; ++i) {
      Abar[i] = 0.0;
    }
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") >= 1) {
#pragma omp parallel for
    for (int i = 0; i < np; ++i) {
      temperature[i] = 0.0;
    }
  }

  if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::entropy") >= 1) {
#pragma omp parallel for
    for (int i = 0; i < np; ++i) {
      entropy[i] = 0.0;
    }
  }

  if (CCTK_EQUALS(initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS(initial_data_setup_method, "init_single_level")) {
    /* do nothing */
  } else if (CCTK_EQUALS(initial_data_setup_method, "init_all_levels")) {

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho") >= 2) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        rho_p[i] = 0.0;
        vel_p[i] = 0.0;
        vel_p[i + np] = 0.0;
        vel_p[i + 2 * np] = 0.0;
        w_lorentz_p[i] = 1.0;
        eps_p[i] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Abar") >= 2) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Abar_p[i] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") >= 2) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        temperature_p[i] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::entropy") >= 2) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        entropy_p[i] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho") >= 3) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        rho_p_p[i] = 0.0;
        vel_p_p[i] = 0.0;
        vel_p_p[i + np] = 0.0;
        vel_p_p[i + 2 * np] = 0.0;
        w_lorentz_p_p[i] = 1.0;
        eps_p_p[i] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::rho") >= 4) {
      CCTK_WARN(CCTK_WARN_ABORT,
                "Too many active time levels for HydroBase variables");
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Abar") >= 3) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Abar_p_p[i] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::temperature") >= 3) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        temperature_p_p[i] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::entropy") >= 3) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        entropy_p_p[i] = 0.0;
      }
    }

  } else {
    CCTK_WARN(
        CCTK_WARN_ABORT,
        "Unsupported parameter value for InitBase::initial_data_setup_method");
  }
}

void HydroBase_Y_e_one(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const np = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];

#pragma omp parallel for
  for (int i = 0; i < np; ++i) {
    Y_e[i] = 1.0;
  }

  if (CCTK_EQUALS(initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS(initial_data_setup_method, "init_single_level")) {
    /* do nothing */
  } else if (CCTK_EQUALS(initial_data_setup_method, "init_all_levels")) {

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Y_e") >= 2) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Y_e_p[i] = 1.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Y_e") >= 3) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Y_e_p_p[i] = 1.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Y_e") >= 4) {
      CCTK_WARN(CCTK_WARN_ABORT,
                "Too many active time levels for HydroBase::Y_e");
    }

  } else {
    CCTK_WARN(
        CCTK_WARN_ABORT,
        "Unsupported parameter value for InitBase::initial_data_setup_method");
  }
}

void HydroBase_Bvec_zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const np = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];

#pragma omp parallel for
  for (int i = 0; i < np; ++i) {
    Bvec[i] = 0.0;
    Bvec[i + np] = 0.0;
    Bvec[i + 2 * np] = 0.0;
  }

  if (CCTK_EQUALS(initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS(initial_data_setup_method, "init_single_level")) {
    /* do nothing */
  } else if (CCTK_EQUALS(initial_data_setup_method, "init_all_levels")) {

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Bvec") >= 2) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Bvec_p[i] = 0.0;
        Bvec_p[i + np] = 0.0;
        Bvec_p[i + 2 * np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Bvec") >= 3) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Bvec_p_p[i] = 0.0;
        Bvec_p_p[i + np] = 0.0;
        Bvec_p_p[i + 2 * np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Bvec") >= 4) {
      CCTK_WARN(CCTK_WARN_ABORT,
                "Too many active time levels for HydroBase::Bvec");
    }

  } else {
    CCTK_WARN(
        CCTK_WARN_ABORT,
        "Unsupported parameter value for InitBase::initial_data_setup_method");
  }
}

void HydroBase_Avec_zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const np = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];

#pragma omp parallel for
  for (int i = 0; i < np; ++i) {
    Avec[i] = 0.0;
    Avec[i + np] = 0.0;
    Avec[i + 2 * np] = 0.0;
  }

  if (CCTK_EQUALS(initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS(initial_data_setup_method, "init_single_level")) {
    /* do nothing */
  } else if (CCTK_EQUALS(initial_data_setup_method, "init_all_levels")) {

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Avec") >= 2) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Avec_p[i] = 0.0;
        Avec_p[i + np] = 0.0;
        Avec_p[i + 2 * np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Avec") >= 3) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Avec_p_p[i] = 0.0;
        Avec_p_p[i + np] = 0.0;
        Avec_p_p[i + 2 * np] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Avec") >= 4) {
      CCTK_WARN(CCTK_WARN_ABORT,
                "Too many active time levels for HydroBase::Avec");
    }

  } else {
    CCTK_WARN(
        CCTK_WARN_ABORT,
        "Unsupported parameter value for InitBase::initial_data_setup_method");
  }
}

void HydroBase_Aphi_zero(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const np = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];

#pragma omp parallel for
  for (int i = 0; i < np; ++i) {
    Aphi[i] = 0.0;
  }

  if (CCTK_EQUALS(initial_data_setup_method, "init_some_levels") ||
      CCTK_EQUALS(initial_data_setup_method, "init_single_level")) {
    /* do nothing */
  } else if (CCTK_EQUALS(initial_data_setup_method, "init_all_levels")) {

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Aphi") >= 2) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Aphi_p[i] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Aphi") >= 3) {
#pragma omp parallel for
      for (int i = 0; i < np; ++i) {
        Aphi_p_p[i] = 0.0;
      }
    }

    if (CCTK_ActiveTimeLevels(cctkGH, "HydroBase::Aphi") >= 4) {
      CCTK_WARN(CCTK_WARN_ABORT,
                "Too many active time levels for HydroBase::Aphi");
    }

  } else {
    CCTK_WARN(
        CCTK_WARN_ABORT,
        "Unsupported parameter value for InitBase::initial_data_setup_method");
  }
}

void HydroBase_InitExcisionMask(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int const np = cctk_ash[0] * cctk_ash[1] * cctk_ash[2];

#pragma omp parallel for
  for (int i = 0; i < np; ++i) {
    hydro_excision_mask[i] = HYDRO_EXCISION_NORMAL;
  }
}
