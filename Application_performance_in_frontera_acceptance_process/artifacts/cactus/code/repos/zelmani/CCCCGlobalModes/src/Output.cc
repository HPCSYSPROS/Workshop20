#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <cstdio>
#include <cstring>
#include <cerrno>

//#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"


#include "util_Table.h"


extern "C" { void CCCCGlobalModes_Output(CCTK_ARGUMENTS);
}

// open file and report errors and abort if failed
#define fopen_check(fn, mode) fopen_check_fun(fn, mode, __LINE__)
static FILE *fopen_check_fun(const char *fn, const char *mode, int line);

void CCCCGlobalModes_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  //  using namespace Carpet;

  int do_stuff = ( (((cctk_iteration) % compute_every) == 0) || cctk_iteration == 0) &&
    (cctk_time >= (start_time*2.03e2)) && (do_shibata || do_saijo || do_CoM);

  if(!do_stuff) return;


  // Create the output directory
  ostringstream outdirbuf;

  if (CCTK_MyProc(cctkGH)==0) {
    CCTK_CreateDirectory (0755, outdirbuf.str().c_str());
    
    outdirbuf << out_dir;

    ostringstream shibata_outfilenamebuf;      
    ostringstream saijo_outfilenamebuf;      
    ostringstream CoM_outfilenamebuf;      
    ostringstream P_outfilenamebuf;      
    
    shibata_outfilenamebuf << outdirbuf.str().c_str() << "/" 
		   << "Shibata_Modes.dat";

    saijo_outfilenamebuf << outdirbuf.str().c_str() << "/" 
		   << "Saijo_Modes.dat";

    CoM_outfilenamebuf << outdirbuf.str().c_str() << "/" 
		   << "CoM.dat";

    P_outfilenamebuf << outdirbuf.str().c_str() << "/" 
		   << "TotalS.dat";

     
    if(do_shibata) {

      FILE* shibata_outfile=
	fopen(shibata_outfilenamebuf.str().c_str(),"r");
      if (! shibata_outfile ) {
	shibata_outfile=fopen_check(shibata_outfilenamebuf.str().c_str(), "w");
	fprintf(shibata_outfile,"%10s %15s %15s %15s %15s %15s %15s %15s %15s\n",
		"# it","time","Ixx","Ixy","Iyy","eta_plus","eta_cross","eta","rho_max");
	fclose(shibata_outfile);
      } else {
	fclose(shibata_outfile);
      }
      
      shibata_outfile=fopen_check(shibata_outfilenamebuf.str().c_str(),"a");
      fprintf(shibata_outfile,"%10d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",
	      cctk_iteration,cctk_time*(4.9261e-6),*Ixx,
	      *Ixy,*Iyy,*eta_plus,*eta_cross,*eta,*global_rho_max);

      fclose(shibata_outfile);

    } // do_shibata 

    if(do_saijo) {

      FILE* saijo_outfile=
	fopen(saijo_outfilenamebuf.str().c_str(),"r");
      if (! saijo_outfile ) {
	saijo_outfile=fopen_check(saijo_outfilenamebuf.str().c_str(),"a");
	fprintf(saijo_outfile,"%10s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n",
		"# it","time","di_re","di_im","quad_re","quad_im","sextu_re","sextu_im","rho_max","mass");
	fclose(saijo_outfile);
      } else {
	fclose(saijo_outfile);
      }
      
      saijo_outfile=fopen_check(saijo_outfilenamebuf.str().c_str(),"a");
      fprintf(saijo_outfile,"%10d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",
	      cctk_iteration,cctk_time*(4.9261e-6),
	      *di_re,*di_im,*quad_re,*quad_im,*sextu_re,*sextu_im,*global_rho_max,*total_mass);

      fclose(saijo_outfile);

    } // do_saijo

    if(do_CoM) {

      FILE* CoM_outfile=
	fopen(CoM_outfilenamebuf.str().c_str(),"r");
      if (! CoM_outfile ) {
	CoM_outfile=fopen_check(CoM_outfilenamebuf.str().c_str(),"a");
	fprintf(CoM_outfile,"%10s %15s %15s %15s %15s %15s %15s %15s\n",
		"# it","time","Mx","My","Mz","Mr","global_rho_max","mass");
	fclose(CoM_outfile);
      } else {
	fclose(CoM_outfile);
      }
      
      CoM_outfile=fopen_check(CoM_outfilenamebuf.str().c_str(),"a");
      fprintf(CoM_outfile,"%10d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",
	      cctk_iteration,cctk_time*(4.9261e-6),
	      *Mx,*My,*Mz,*Mr,*global_rho_max,*total_mass);

      fclose(CoM_outfile);

    } // do_CoM


    if(do_P) {

      FILE* P_outfile=
	fopen(P_outfilenamebuf.str().c_str(),"r");
      if (! P_outfile ) {
	P_outfile=fopen_check(P_outfilenamebuf.str().c_str(),"a");
	fprintf(P_outfile,"%10s %15s %15s %15s %15s %15s %15s %15s\n",
		"# it","time","Px","Py","Pz","Pr","rho_max","mass");
	fclose(P_outfile);
      } else {
	fclose(P_outfile);
      }
      
      P_outfile=fopen_check(P_outfilenamebuf.str().c_str(),"a");
      fprintf(P_outfile,"%10d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",
	      cctk_iteration,cctk_time*(4.9261e-6),
	      *Px,*Py,*Pz,*Pr,*global_rho_max,*total_mass);

      fclose(P_outfile);

    } // do_P


  }
  

  return;

}

static FILE *fopen_check_fun(const char *fn, const char *mode, int line)
{
  FILE *fh = fopen(fn, mode);
  if (! fh) {
    CCTK_VWarn(CCTK_WARN_ABORT, line, __FILE__, CCTK_THORNSTRING,
               "Could not open file '%s': %s", fn, strerror(errno));
  }
  return fh;
}
