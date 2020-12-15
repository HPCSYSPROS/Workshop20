#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"

#include <string>
#include <cstdio>
#include <cassert>

#include "nuc_eos.hh"

#define DIM(v) (sizeof(v)/sizeof((v)[0]))

extern "C"
void EOS_OMNI_dumptable(CCTK_ARGUMENTS)
{
  using namespace nuc_eos;
  using namespace nuc_eos_private;

  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  int irho,itemp,iye,n;
  int strlength1,strlength2;
  std::string fullpath;

  if (CCTK_MyProc(cctkGH) != 0)
     return;

  fullpath = std::string(out_dir)+"/"+std::string(dump_nuceos_table_name);

  CCTK_INFO("*******************************");
  CCTK_INFO("Dumping nuc_eos table file in ASCII:");
  CCTK_INFO(fullpath.c_str());
  CCTK_INFO("*******************************");

  FILE *fh = fopen(fullpath.c_str(), "w");
  assert(fh);

  fprintf(fh, "# % 20s\n% 4d\n", "nrho:",nrho);
  fprintf(fh, "# % 20s\n% 4d\n", "ntemp:",ntemp);
  fprintf(fh, "# % 20s\n% 4d\n", "nye:",nye);

  fprintf(fh, "# % 20s\n%18.9E\n", "energy shift:",energy_shift);

  fprintf(fh, "# % 20s\n%18.9E %18.9E\n", "rho min and max:",eos_rhomin/RHOGF,eos_rhomax/RHOGF);
  fprintf(fh, "# % 20s\n%18.9E %18.9E\n", "ye min and max:",eos_yemin,eos_yemax);
  fprintf(fh, "# % 20s\n%18.9E %18.9E\n", "temp min and max:",eos_tempmin,eos_tempmax);

  fprintf(fh, "# % 20s\n\n", "log rho points:");
  for(irho=0;irho<nrho;++irho) {
     fprintf(fh, "%18.9E\n", (logrho[irho] - log(RHOGF))/log(10.));
  }
  fprintf(fh, "# % 20s\n", "log temp points:");
  for(itemp=0;itemp<ntemp;++itemp) {
     fprintf(fh, "%18.9E\n", logtemp[itemp]/log(10.));
  }
  fprintf(fh, "# % 20s\n", "ye points:");
  for(iye=0;iye<nye;++iye) {
     fprintf(fh, "%18.9E\n", yes[iye]);
  }

  fprintf(fh, "# % 20s\n% 4d\n", "nvars:",NTABLES);
  fprintf(fh, "# % 20s\n", "table mappings:");
  fprintf(fh, "# % 20s\n", " 1 -> logpress");
  fprintf(fh, "# % 20s\n", " 2 -> logenergy");
  fprintf(fh, "# % 20s\n", " 3 -> entropy");
  fprintf(fh, "# % 20s\n", " 4 -> munu");
  fprintf(fh, "# % 20s\n", " 5 -> cs2");
  fprintf(fh, "# % 20s\n", " 6 -> dedT");
  fprintf(fh, "# % 20s\n", " 7 -> dpdrhoe");
  fprintf(fh, "# % 20s\n", " 8 -> dpderho");
  fprintf(fh, "# % 20s\n", " 9 -> muhat");
  fprintf(fh, "# % 20s\n", "10 -> mu_e");
  fprintf(fh, "# % 20s\n", "11 -> mu_p");
  fprintf(fh, "# % 20s\n", "12 -> mu_n");
  fprintf(fh, "# % 20s\n", "13 -> xa");
  fprintf(fh, "# % 20s\n", "14 -> xh");
  fprintf(fh, "# % 20s\n", "15 -> xn");
  fprintf(fh, "# % 20s\n", "16 -> xp");
  fprintf(fh, "# % 20s\n", "17 -> abar");
  fprintf(fh, "# % 20s\n", "18 -> zbar");
  fprintf(fh, "# % 20s\n", "19 -> gamma");
  
  const double ctable[][2] = { // un-convert units
    { log(10.0), log(PRESSGF) }, // pressure
    { log(10.0), log(EPSGF) }, // eps
    { 1. },
    { 1. },
    { LENGTHGF*LENGTHGF/TIMEGF/TIMEGF }, // cs2
    { EPSGF }, // dedT
    { PRESSGF/RHOGF }, // dpdrhoe
    { PRESSGF/EPSGF },// dpderho
  };
  for(irho=0;irho<nrho;++irho) {
     for(itemp=0;itemp<ntemp;++itemp) {
        for(iye=0;iye<nye;++iye) {
           for(n=0;n<NTABLES;++n) {
              // transpose nvar to front as we go along
              const int idx = n + NTABLES*(irho + nrho*(itemp + ntemp*iye));
              double val = alltables[idx];
              if(n < DIM(ctable)) val = (val - ctable[n][1]) / ctable[n][0];
              fprintf(fh, "% 4d % 4d % 4d % 4d %18.9E\n", irho+1,itemp+1,iye+1,n+1,val);
           }
        }
     }
  }

  fclose(fh);
}
