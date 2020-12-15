#***********************************************************************************
#                                                                                  *
# Random Acceleration Molecular Dynamics (RAMD)                                    *
# Implementation for NAMD v2.7                                                     *
# September 2009                                                                   *
#                                                                                  *
# Copyright (c) 2009, EML Research gGmbH, Heidelberg, Germany                      * 
# Author:  Vlad Cojocaru                                                           *
# Email: vlad.cojocaru@eml-r.villa-bosch.de                                        *
#                                                                                  *
# This Tcl script includes all the features of the implementation for AMBER 8      *
# (fortran code by  Tim Johann & Ting Wang in the MCM group @ EML Research gGmbH)  *
# See: http://projects.eml.org/mcmsoft/amberpatches/                               *
#                                                                                  *
# The first Tcl script to run RAMD in NAMD (v2.5+) was written by Harish Vashisth  *
# Ref: Vashisth H et al, Biophys J. 2008 Nov 1;95(9):4193-204. Epub 2008 Aug 1     *
#                                                                                  *
# The Vashisth script inspired some lines of this script, but mainly this script   *
#    is based on the fortran code written for AMBER 8                              *
#                                                                                  *
# The structure of this script is inspired by the Adaptive Biasing Force module    *
#    distributed with NAMD 2.6                                                     *
#                                                                                  *
# The original RAMD method is described in:                                        *
#    Ref1: Luedemann,S.K.,Lounnas,V.and R.C.Wade.,                                 *
#          J Mol Biol, 303:797-811 (2000)                                          *
#    Ref2: Schleinkofer,K.,Sudarko,Winn,P.,Luedemann,S.K.and R.C.Wade,             *
#          EMBO Reports, 6, 584-589 (2005)                                         *
#                                                                                  *
# Disclaimer:  This script is for research purposes only. EML Research does not    *
#              assume any responsibility for the software or its use.              * 
#                                                                                  *
# !!! Quantitative reproducibility of the results obtained with AMBER 8            *
#     is not possible due to a numerical errors and such like                      * 
#                                                                                  *
# !!! This script is under development.                                            *
#                                                                                  *
#   The script along with usage examples is available at                           *
#   http://projects.eml.org/mcm/software                                           *
#***********************************************************************************

#*******************************************************
# Startup                                              *
#*******************************************************
package provide ramd 4.1

#*******************************************************
# Parameter definitions
#*******************************************************

namespace eval ::RAMD {
 set version "4.1"
 if {! [info exists RAMDdir]} { set RAMDdir [file dirname [info script]] }
 # If it fails, try the local directory
 if { $RAMDdir == "" } { set RAMDdir "." }
 
 TclForces		on
 TclForcesScript	$RAMDdir/ramd-4.1_script.tcl

 array set defaults {
  ramdSteps               50
  accel                    0.25
  rMinRamd                 0.01
  forceOutFreq             0
  firstProtAtom            1
  ramdSeed             14253
  mdSteps                  0
  mdStart                 no
  maxDist                 50
  debugLevel               0
 }

 set mandatory "firstRamdAtom lastRamdAtom lastProtAtom"
 set silent "rMinMd"

 array set capitals {}
 foreach param [concat $mandatory $silent [array names defaults]] {
  set capitals([string tolower $param]) $param
  # not set yet
  set alreadySet($param) 0
 }
} ;# namespace


proc ramd { keyword value } {
 set ::RAMD::keyword $keyword
 set ::RAMD::value $value

 namespace eval ::RAMD {
  # Build list of all allowed parameter names
  set list [array names capitals]
  set lowercase [string tolower $keyword]
  # Process parameters
  if {[lsearch $list $lowercase] != -1} {
   set keyword $capitals($lowercase)
   if { $alreadySet($keyword) } {
    print "RAMD> WARNING - multiple definitions of parameter $keyword"	
   }
   set $keyword $value
   set alreadySet($keyword) 1
   return
  } else {
   error [format "Unknown RAMD keyword: %s" $keyword]
  }
 } ;# namespace

} ;# proc ramd

# define upper-case synonyms to proc abf 
proc RAMD { keyword value } {
    ramd $keyword $value
}
proc Ramd { keyword value } {
    ramd $keyword $value
}

