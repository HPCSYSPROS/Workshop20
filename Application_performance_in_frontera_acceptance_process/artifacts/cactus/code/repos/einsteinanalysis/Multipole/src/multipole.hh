#include "cctk.h"
#include "cctk_Arguments.h"

// Multipole_Calc
//      This is the main scheduling file.  Because we are completely local here
//      and do not use cactus arrays etc, we schedule only one function and then
//      like program like one would in C, C++ with this function taking the 
//      place of int main(void).
//
//      This function calls functions to accomplish 3 things:
//        1) Interpolate psi4 onto a sphere
//        2) Integrate psi4 with the ylm's over that sphere
//        2) Output the mode decomposed psi4
extern "C" void Multipole_Calc(CCTK_ARGUMENTS);

