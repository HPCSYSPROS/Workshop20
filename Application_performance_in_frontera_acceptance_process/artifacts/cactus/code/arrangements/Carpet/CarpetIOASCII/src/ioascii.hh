#ifndef CARPETIOASCII_HH
#define CARPETIOASCII_HH

#include <vector>

#include <cctk.h>
#include <cctk_Arguments.h>

namespace CarpetIOASCII {

using namespace std;

// Scheduled functions
extern "C" {
int CarpetIOASCIIStartup(void);
void CarpetIOASCIIInit(CCTK_ARGUMENTS);
}

// routines which are independent of the output dimension
static ibbox GetOutputBBox(const cGH *cctkGH, int group, int rl, int m, int c,
                           const ibbox &ext);

static void GetCoordinates(const cGH *cctkGH, int m, const cGroup &groupdata,
                           const ibbox &ext, CCTK_REAL &coord_time,
                           rvect &coord_lower, rvect &coord_upper);

static int GetGridOffset(const cGH *cctkGH, int m, int dir, const char *iparam,
                         const char *iglobal, const char *cparam,
                         const char *cglobal, CCTK_REAL cfallback);
static int CoordToOffset(const cGH *cctkGH, int m, int dir, CCTK_REAL coord,
                         int ifallback);

// Everything is a class template, so that it can easily be
// instantiated for all output dimensions

template <int outdim> struct IOASCII {

  // name of the output directory
  static char *my_out_dir;

  // list of variables to output
  static char *my_out_vars;

  // I/O request description list (for all variables)
  static vector<ioRequest *> requests;

  // Scheduled functions
  static int Startup();

  // Registered functions
  static void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cctkGH);

  static int OutputGH(const cGH *cctkGH);
  static int OutputVarAs(const cGH *cctkGH, const char *varname,
                         const char *alias);
  static int TimeToOutput(const cGH *cctkGH, int vindex);
  static int TriggerOutput(const cGH *cctkGH, int vindex);

  // Other functions
  static void CheckSteerableParameters(const cGH *cctkGH);

  static bool DidOutput(const cGH *cctkGH, int vindex, string basefilename,
                        bool &is_new_file, bool &truncate_file);

  static bool DirectionIsRequested(const vect<int, outdim> &dirs);

  static void OutputDirection(const cGH *cctkGH, int vindex, string alias,
                              string basefilename,
                              const vect<int, outdim> &dirs, bool is_new_file,
                              bool truncate_file);

  static void OpenFile(const cGH *cctkGH, int m, int vindex, string alias,
                       string basefilename, const vect<int, outdim> &dirs,
                       bool is_new_file, bool truncate_file, fstream &file);

  static void CloseFile(const cGH *cctkGH, fstream &file);

  static ivect GetOutputOffset(const cGH *cctkGH, int m,
                               const vect<int, outdim> &dirs);

}; // struct IOASCII

template <int outdim>
void WriteASCII(ostream &os, vector<gdata *> const &gfdatas,
                const bbox<int, dim> &gfext, const int vi, const int time,
                const vect<int, dim> &org, const vect<int, outdim> &dirs,
                const int rl, const int ml, const int m, const int c,
                const int tl, const CCTK_REAL coord_time,
                const vect<CCTK_REAL, dim> &coord_lower,
                const vect<CCTK_REAL, dim> &coord_upper,
                vector<gdata *> const &gfcoords);

} // namespace CarpetIOASCII

#endif // ! defined CARPETIOASCII_HH
