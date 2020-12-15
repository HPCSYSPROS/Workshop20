#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <fstream>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <cstring>
#include <cerrno>

//#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"


#include "util_Table.h"


extern "C" { void ZelmaniCoMShift_Output(CCTK_ARGUMENTS);
}

// open file and report errors and abort if failed
#define fopen_check(fn, mode) fopen_check_fun(fn, mode, __LINE__)
static FILE *fopen_check_fun(const char *fn, const char *mode, int line);

void ZelmaniCoMShift_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  //  using namespace Carpet;

  if (!(cctk_time >= (start_time))){
    return;
  }


  // Create the output directory
  ostringstream outdirbuf;

  if (CCTK_MyProc(cctkGH)==0) {
    CCTK_CreateDirectory (0755, outdirbuf.str().c_str());
    
    outdirbuf << out_dir;

    ostringstream CoM_outfilenamebuf;      

    CoM_outfilenamebuf << outdirbuf.str().c_str() << "/" 
		   << "CoMShift.dat";

    FILE* CoM_outfile=
      fopen(CoM_outfilenamebuf.str().c_str(),"r");
    if (! CoM_outfile ) {
      CoM_outfile=fopen_check(CoM_outfilenamebuf.str().c_str(),"a");
      fprintf(CoM_outfile,"%10s %15s %15s %15s %15s %15s %15s\n",
	      "# it","time","Mx","My","Mz","Mr","mass");
	fclose(CoM_outfile);
    } else {
      fclose(CoM_outfile);
    }
    
    CoM_outfile=fopen_check(CoM_outfilenamebuf.str().c_str(),"a");
    fprintf(CoM_outfile,"%10d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",
	    cctk_iteration,cctk_time*(4.9261e-6),
	    *Mx,*My,*Mz,*Mr,*Mass);
    
    fclose(CoM_outfile);
    
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
