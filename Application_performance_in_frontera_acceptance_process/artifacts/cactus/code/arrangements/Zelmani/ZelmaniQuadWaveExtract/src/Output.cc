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

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"


#include "util_Table.h"


extern "C" { void ZelmaniQuadWaveExtract_Output(CCTK_ARGUMENTS);
}

// open file and report errors and abort if failed
#define fopen_check(fn, mode) fopen_check_fun(fn, mode, __LINE__)
static FILE *fopen_check_fun(const char *fn, const char *mode, int line);

void ZelmaniQuadWaveExtract_Output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  if(!(*dostuff) && !(cctk_iteration==0)) {
    return;
  }

  // Create the output directory
  ostringstream outdirbuf;

  if (CCTK_MyProc(cctkGH)==0) {

    //    fprintf(stderr,"\ncctk_iteration: %d\n",*dostuff,cctk_iteration);
    CCTK_CreateDirectory (0755, outdirbuf.str().c_str());
    
    outdirbuf << out_dir;

    ostringstream outfilenamebuf;      
    outfilenamebuf << outdirbuf.str().c_str() << "/" 
		   << "ZelmaniQuadWave.dat";

    //cout << outfilenamebuf.str().c_str() << endl;

    FILE* outfile=fopen(outfilenamebuf.str().c_str(),"r");
    if (! outfile ) {
      outfile=fopen_check(outfilenamebuf.str().c_str(),"a");
      fprintf(outfile,"%10s %15s %15s %15s %15s %15s %15s %15s %15s %15s %18s\n",
	      "# it","time","Idotxx","Idotxy","Idotxz","Idotyy","Idotyz","Idotzz",
	      "rhomax","volume","Bary Mass");
      fclose(outfile);
    } else {
      fclose(outfile);
    }

    //    fprintf(stderr,"blah\n");
    
    //    if(verbose) {
    //  CCTK_INFO("Opening output file");
    //}
    outfile=fopen_check(outfilenamebuf.str().c_str(),"a");
    fprintf(outfile,"%10d %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E %15.8E\n",
    	    cctk_iteration,cctk_time*(4.9261e-6),*Idotxx,
    	    *Idotxy,*Idotxz,*Idotyy,*Idotyz,*Idotzz,*rhomax, *volume, *BaryMass);

    fclose(outfile);
    //if(verbose) {
    // CCTK_INFO("Output file closed");
    //}
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

