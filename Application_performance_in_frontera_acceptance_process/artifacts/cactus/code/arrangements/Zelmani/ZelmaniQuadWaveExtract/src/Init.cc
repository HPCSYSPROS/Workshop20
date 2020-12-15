#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "carpet.hh"
#include "defs.hh"
#include "vect.hh"


#include "util_Table.h"


extern "C" { void ZelmaniQuadWaveExtract_Init(CCTK_ARGUMENTS);
}

void ZelmaniQuadWaveExtract_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  using namespace std;
  using namespace Carpet;

  // we should do something, right?

  // Yeah. Let's do at least SOMETHING useful!

  *dostuff = 0; // don't do anything until I tell you to.

  *BaryMass = 0.0e0;
  *Idotxx = 0.0e0;
  *Idotxy = 0.0e0;
  *Idotxz = 0.0e0;
  *Idotyy = 0.0e0;
  *Idotyz = 0.0e0;
  *Idotzz = 0.0e0;
 
  *volume = 0.0e0;
  *rhomax = 0.0e0;


  /*
  int nx = cctk_lsh[0];
  int ny = cctk_lsh[1];
  int nz = cctk_lsh[2];

  for (int k=0;k<nz;k++)
    for (int j=0;j<ny;j++)
      for (int i=0;i<nx;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	dProperMass[index] = 0.0e0;
	dIdotxx[index] = 0.0e0;
	dIdotyy[index] = 0.0e0;
	dIdotzz[index] = 0.0e0;
	dIdotxy[index] = 0.0e0;
	dIdotxz[index] = 0.0e0;
	dIdotyz[index] = 0.0e0;
	
      }

*/

  return;
}

