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

enabletotalforces

namespace eval ::RAMD {
 print RAMD:
 print RAMD:   -------------------------------------------------------------------  
 print RAMD:   Random Acceleration Molecular Dynamics Simulation version $version
 print RAMD:   -------------------------------------------------------------------  
 print RAMD:
 #***** Assign default values for the parameters not specified in the configuration file
 foreach option [array names defaults] {
  if {! [info exists $option]} {
   set $option $defaults($option)
   print [format "RAMD: %25s = %s" $option [expr $$option]]
  } elseif { [info exists $option] } {
   print [format "RAMD: %25s = %s" $option [expr $$option]]
  }
 }
 #***** Check if mandatory parameters are specified in the configuration file
 foreach var $mandatory {
  if {! [info exists $var]} {
   error "RAMD: Mandatory parameter '$var' is not set -- cannot start RAMD"
  } else {
   print [format "RAMD: %25s = %s" $var [expr $$var]]
  }
 }
 #***** Check if 'forceOutFreq' is equal to 1; exit with error if that's the case
 if { $forceOutFreq == 1 } { error "RAMD: ERROR: 'forceOutFreq' parameter may not be 1" } 

 #***** Check if 'mdSteps' is specified in the configuration file

 #***** Performed pure RAMD if 'mdSteps' = 0
 if { $mdSteps == 0 } { 
  
  #***** Check that the number of ramd steps is a multiple of 'forceOutFreq'; exit with error if not
  set r [expr "$ramdSteps % $forceOutFreq"]
  if { $r != 0 } { error "RAMD: ERROR: The number of RAMD steps is not a multiple of 'forceOutFreq'" } 
 
  print "RAMD: Pure RAMD simulation is performed" 
  
  #***** If 'mdSteps' is 0 and "mdStart" is yes, give a warning
  if { $mdStart == "yes" } { 
   print "RAMD: WARNING: 'mdStart' has no meaning for pure RAMD simulation; it will be ignored" 
  }

  #***** If 'mdSteps' is 0 and "rMinMd" is set, give a warning
  if { [info exists rMinMd] } {
   print "RAMD: WARNING: 'rMinMd' specified while 'mdSteps' is 0"
   print "RAMD: WARNING: For combined RAMD-MD simulation 'mdSteps' must be greater than 0"
   print "RAMD: WARNING: Ignore 'rMinMd' and perform pure RAMD simulation"
  }
  
 }

 #***** Perform combined RAMD with MD simulation if 'mdSteps' is not 0 and 'rMinMd' is specified
 if { $mdSteps != 0  } { 
  
  if { [info exists rMinMd] } {

   #***** Check that the number of ramd and md steps are each a multiple of 'forceOutFreq'
   #***** Exit with error if that's not the case
   set r1 [expr "$ramdSteps % $forceOutFreq"]
   set r2 [expr "$mdSteps % $forceOutFreq"]
   
   if { $r1 != 0 || $r2 != 0 } { 
    error "RAMD: ERROR: The number of RAMD or MD steps must be multiple of 'forceOutFreq'" 
   } 

   foreach svar $silent {
    print [format "RAMD: %25s = %s" $svar [expr $$svar]]
   }
  
   print "RAMD:"
   print "RAMD: Combined RAMD-MD simulation is performed"
  
  } elseif { ! [info exists rMinMd] } {
   
   #***** If 'mdSteps' is not 0, exit with error if 'rMinMd' is not specified
   error "RAMD: ERROR: parameter 'rMinMd' not set: 'rMinMd' is required if 'mdSteps' is greater than 0"
  
  }

 }
  
 print "RAMD:"
   
 #***** Make a list of all the atoms on which the force will be applied
 set ramdAtoms {}
 for { set i $firstRamdAtom } { $i <= $lastRamdAtom } { incr i } { lappend ramdAtoms $i }
 print "RAMD: Atoms subject to the random acceleration are: $ramdAtoms"
 foreach ramdAtom $ramdAtoms { addatom $ramdAtom }
 #***** Define a group of the ligand atoms; the force will be applied on the center of mass of this group
 set ramdGroup [ addgroup $ramdAtoms ]

 #***** Define a group containing all protein atoms
 set protAtoms {}
 for { set i $firstProtAtom } { $i <= $lastProtAtom } { incr i } { lappend protAtoms $i }
 foreach protAtom $protAtoms { addatom $protAtom }
 set protGroup [ addgroup $protAtoms ]

 #***** Some variable initialization 
 set timeStep 0; set exitFlag 0; 
 set prevLigCOM "0.0 0.0 0.0"; set prevProtCOM "0.0 0.0 0.0"; set prevDist 0;

 #***** Initialization of simulation flags
 if { $mdSteps == 0 } {
  set ramdFlag 1; set mdFlag 0; set ramdStep 0; set mdStep 0;
 } elseif { $mdSteps != 0 } { 
  if { $mdStart == "yes" } { set ramdFlag 0; set mdFlag 1; set ramdStep 0; set mdStep 0; }
  if { $mdStart == "no" } { set ramdFlag 1; set mdFlag 0; set ramdStep 0; set mdStep 0; }  
 }
 
} ;# namespace
 
#***** In root namespace (::) for all procedures we have to add the following procedure definition
proc veclength {v} {
 return [expr {sqrt([veclength2 $v])}]
}
#***** Source the vectors and matrices procedures from VMD
source $RAMD::RAMDdir/vectors.tcl

#***********************************************************
# PROCEDURE TO GENERATE RANDOMLY ORIENTED ACCELERATION 
#***********************************************************
proc genRandAccel { timeStep } {
namespace eval ::RAMD {

 set pi [expr "2.0*asin(1.0)"]

 #***** Generate new random orientation of the ramd force
 set randTheta [expr "rand()"]
 set randPsi [expr "rand()"]
 set theta [expr "2*$pi*$randTheta"]
 set psi [expr "$pi*$randPsi"]
 set rx [expr "cos($theta)*sin($psi)"]
 set ry [expr "sin($theta)*sin($psi)"]
 set rz [expr "cos($psi)"]
 set r "$rx $ry $rz"
 set lenr [veclength $r]

 # Acceleration is given in kcal/mol*A*amu  in the NAMD configuration file (multiply with 418.68 to get A/ps^2)
 set vecAccel [vecscale [expr "$accel"] $r ]
  
 return 
 
} ;# namespace
} ;# proc genRandAccel {timestep}


#*****************************************************************************
# PROCEDURE TO EVALUATE THE DISTANCE TRAVELLED BY THE LIGAND IN N RAMD STEPS
#*****************************************************************************
proc evalWalkDist { timeStep prevLigCOM prevProtCOM currLigCOM currProtCOM } {
namespace eval ::RAMD {
 
 #***** Compute the relative position of the ligand com with regard to the protein com
 set prevRelLigCOM [ vecsub $prevLigCOM $prevProtCOM ]
 set currRelLigCOM [ vecsub $currLigCOM $currProtCOM ]
 
 #***** Compute the distance travelled by the ligand com during a ramd or md stint
 set vecWalkDist [vecsub $currRelLigCOM $prevRelLigCOM]
 set walkDist [veclength $vecWalkDist]

 set vecWalkDistX [lindex $vecWalkDist 0]
 set vecWalkDistY [lindex $vecWalkDist 1]
 set vecWalkDistZ [lindex $vecWalkDist 2]
 
 return  

} ;# namespace
} ;# proc evalWalkDist
  

#**************************************************************
# PROCEDURE TO APPLY THE FORCE WHICH IS CALLED EVERY TIME STEP 
#**************************************************************
proc calcforces {} {
namespace eval ::RAMD {
 
 #***** Terminate NAMD if the ligand has exited from the protein
 if { $exitFlag == 1 } {
  print "EXIT: $timeStep  > MAX DISTANCE LIGAND COM - PROTEIN COM REACHED"
  print "EXIT: $timeStep  > LIGAND EXIT EVENT DETECTED: STOP SIMULATION"
  print "EXIT: $timeStep  > EXIT NAMD"
  set process [pid]
  exec kill -9 $process
 } 

 if { [ array exists coords ] } { array unset coords }
 if { [ array exists masses ] } { array unset masses }
 if { [ array exists extForces ] } { array unset extForces }
 if { [ array exists totForces ] } { array unset totForces }

 #***** Load coordinates for all the atoms and groups defined
 loadcoords coords
 #***** Load masses for all the atoms and groups defined
 loadmasses masses
 #***** Load external forces from previous time step for all the atoms and groups defined
 loadforces extForces
 #***** Load total forces from previous time step for all the atoms and groups defined
 loadtotalforces totForces 
 
 #***** Calculate the mass of the ligand
 set ligMass 0
 foreach ramdAtom $ramdAtoms {
  set ligMass [expr $ligMass + $masses($ramdAtom)]
 }

 #***** Calculate the position of protein and ligand COM
 set protCOM "$coords($protGroup)"
 set ligCOM "$coords($ramdGroup)"
  
 #***** Initialize ramd simulation or combined ramd-md simulation that begins with ramd
 if { $timeStep == 0 && $ramdFlag == 1 && $mdFlag == 0 } {
  
  expr "srand($ramdSeed)"
  
  set vMin [ expr "($rMinRamd)/($ramdSteps)" ]
  if { $mdSteps == 0 } { 
   print "RAMD: $timeStep  ***** INITIALIZE RAMD SIMULATION *****" 
  } else { 
   print "RAMD: $timeStep  ***** INITIALIZE COMBINED RAMD-MD SIMULATION *****" 
  }
  print "RAMD: $timeStep     >>> minimum travelled distance (A): $rMinRamd"
  print "RAMD: $timeStep     >>> minimum velocity (A/fs): $vMin"

  #***** Initial com positions
  set currLigCOM $ligCOM; set currProtCOM $protCOM  
  print "RAMD: $timeStep     >>> LIGAND COM IS: $currLigCOM"
  print "RAMD: $timeStep     >>> PROTEIN COM IS: $currProtCOM"

  #***** Evaluate initial distance between ligand com and protein com
  set currDist [veclength [vecsub $currLigCOM $currProtCOM]]
  print "RAMD: $timeStep     >>> DISTANCE LIGAND COM - PPROTEIN COM IS: DIST = $currDist"

  #***** Generate an initial orientation for the acceleration to be applied when ramd is switched on
  genRandAccel $timeStep
  print "RAMD: $timeStep     >>> INITIAL RANDOM DIRECTION: $r :: ||r|| = $lenr"

  print "RAMD: $timeStep  ***** START WITH $ramdSteps STEPS OF RAMD SIMULATION *****"
 
  #***** Reset the positions of the ligand and protein COMs and the distance ligand com - protein com
  set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

  incr timeStep
  return  

 } 

 #***** Initialize combined ramd-md simulation that begins with standard md 
 if { $timeStep == 0 && $ramdFlag == 0 && $mdFlag == 1 } {

  expr "srand($ramdSeed)"
  
  set vMin [expr "($rMinMd)/($mdSteps)"]

  print "MD: $timeStep  ***** INITIALIZE COMBINED RAMD-MD SIMULATION *****"
  print "MD: $timeStep     >>> minimum travelled distance (A): $rMinMd"
  print "MD: $timeStep     >>> minimum velocity (A/fs): $vMin"

  #***** Initial com positions
  set currLigCOM $ligCOM; set currProtCOM $protCOM
  print "MD: $timeStep     >>> LIGAND COM IS: $currLigCOM"
  print "MD: $timeStep     >>> PROTEIN COM IS: $currProtCOM"
  
  #***** Evaluate initial distance between ligand com and protein com
  set currDist [veclength [vecsub $currLigCOM $currProtCOM]]
  print "MD: $timeStep     >>> DISTANCE LIGAND COM - PPROTEIN COM IS: DIST = $currDist"
  
  #***** Generate an initial orientation for the acceleration to be applied when ramd is switched on
  #***** This is ignored for the first MD stint, and no acceleration is applied
  genRandAccel $timeStep
  print "MD: $timeStep     >>> GENERATED RANDOM DIRECTION: $r :: ||r|| = $lenr"
  print "MD: $timeStep     >>> DIRECTION WILL BE IGNORED FOR THE FIRST $mdSteps MD STEPS"

  print "MD: $timeStep  ***** START WITH $mdSteps STEPS OF STANDARD MD SIMULATION *****"

  #***** Reset the positions of the ligand and protein COMs and the distance ligand com - protein com
  set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

  incr timeStep
  return  

 } 

 
 #***** Perform ramd simulation 
 if { $timeStep != 0 && $ramdFlag == 1 && $mdFlag == 0 && $exitFlag == 0  } {

  #***** Count ramd steps
  incr ramdStep
  if { $debugLevel != 0 } {
   print "RAMD DEBUG: TIMESTEP IS: $timeStep; RAMD STEP IS: $ramdStep; MD STEP IS: $mdStep"
  }
    
  #***** Define and apply the force to each atom of the ligand
  foreach ramdAtom $ramdAtoms {
   set atomMass "$masses($ramdAtom)"
   set atomForce [ vecscale $atomMass $vecAccel ]
   set atomForceValue [ veclength "$atomForce" ]
   addforce $ramdAtom $atomForce
   if { $debugLevel != 0 } { print "RAMD DEBUG: ATOM $ramdAtom: MASS $atomMass: ADD FORCE $atomForceValue" }
   unset atomForce
   unset atomMass
   unset atomForceValue
  }

  #***** Define the force vector that is applied to the center of mass of the ligand
  set force [vecscale $ligMass $vecAccel]
  set fx [lindex $force 0]
  set fy [lindex $force 1]
  set fz [lindex $force 2]
 
  #***** Check the magnitude of the force
  set f [expr "sqrt((($fx)*($fx)+($fy)*($fy)+($fz)*($fz)))"]
  
  #***** Set flag for writing force output
  if { $forceOutFreq != 0 } { set outputFlag [expr "$ramdStep % $forceOutFreq"] }
  
  #***** Write force output every 'forceOutFreq' steps 
  if { [ array exists extForces ] && [ array exists totForces ] && $outputFlag == 0 } {
      
   #***** Write the force vector and direction
   print "RAMD FORCE: $timeStep  > LIGAND COM is: $ligCOM"
   print "RAMD FORCE: $timeStep  > PROTEIN COM IS $protCOM"
   print "RAMD FORCE: $timeStep  > EXTERNAL FORCE VECTOR (F): $force; ||F|| = $f"
   print "RAMD FORCE: $timeStep  > EXTERNAL FORCE DIRECTION (r): $r; ||r|| = $lenr"

   #***** Calculate external and total forces acting on the ligand
   set totLigForceX 0; set totLigForceY 0; set totLigForceZ 0; set totLigForceV 0;
   set extLigForceX 0; set extLigForceY 0; set extLigForceZ 0; set extLigForceV 0;
   foreach ramdAtom $ramdAtoms {
    set atomMass "$masses($ramdAtom)"
    set totAtomForce "$totForces($ramdAtom)"
    set totAtomForceX [ lindex $totAtomForce 0 ]
    set totAtomForceY [ lindex $totAtomForce 1 ]
    set totAtomForceZ [ lindex $totAtomForce 2 ]
    set totAtomForceV [ veclength "$totAtomForce" ]
    set totLigForceX [ expr "$totLigForceX + $totAtomForceX" ]
    set totLigForceY [ expr "$totLigForceY + $totAtomForceY" ]
    set totLigForceZ [ expr "$totLigForceZ + $totAtomForceZ" ]
    set totLigForceV [ expr "$totLigForceV + $totAtomForceV" ]
    set extAtomForce "$extForces($ramdAtom)"
    set extAtomForceX [ lindex "$extAtomForce" 0 ]
    set extAtomForceY [ lindex "$extAtomForce" 1 ]
    set extAtomForceZ [ lindex "$extAtomForce" 2 ]
    set extAtomForceV [ veclength "$extAtomForce" ]
    set extLigForceX [ expr "$extLigForceX + $extAtomForceX" ]
    set extLigForceY [ expr "$extLigForceY + $extAtomForceY" ]
    set extLigForceZ [ expr "$extLigForceZ + $extAtomForceZ" ]
    set extLigForceV [ expr "$extLigForceV + $extAtomForceV" ]
    if { $debugLevel != 0 } { 
     print "RAMD DEBUG: ATOM $ramdAtom: MASS $atomMass: EXT FORCE $extAtomForceV"
     print "RAMD DEBUG: ATOM $ramdAtom: MASS $atomMass: TOT FORCE $totAtomForceV"
    }    
    unset atomMass
    unset totAtomForceX; unset totAtomForceY; unset totAtomForceZ; unset totAtomForceV; unset totAtomForce
    unset extAtomForceX; unset extAtomForceY; unset extAtomForceZ; unset extAtomForceV; unset extAtomForce
   }
   set totLigForce "$totLigForceX $totLigForceY $totLigForceZ"
   set extLigForce "$extLigForceX $extLigForceY $extLigForceZ"
 
   #***** Write external forces acting on the ligand com for debugging purposes
   if { $debugLevel !=0 } {
    print "RAMD DEBUG: $timeStep  > EXTERNAL FORCE ON THE LIGAND COM IS: $extLigForce ($extLigForceV)"
   }

   #***** Write total forces acting on the ligand com
   print "RAMD FORCE: $timeStep  > TOTAL FORCE ON THE LIGAND COM IS: $totLigForce ($totLigForceV)"
   
   unset totLigForce; unset totLigForceV; unset extLigForce; unset extLigForceV;

  } elseif { ! [ array exists extForces ] && [ array exists totForces ] && $outputFlag == 0 } {
  
   error "RAMD: $timeStep  > ERROR: EXTERNAL FORCES NOT PRESENT DURING RAMD STEP: EXIT NAMD"
  
  }

  #***** Set flag for evaluating ramd simulation
  set evalRamdFlag [expr "$ramdStep % $ramdSteps"]

  #***** Evaluate ramd stint
  if { $evalRamdFlag == 0 } {
  
   print "RAMD: $timeStep  ***** EVALUATE $ramdSteps RAMD STEPS AT TIMESTEP $timeStep *****"

   #***** com positions
   set currLigCOM $ligCOM; set currProtCOM $protCOM
   if { $debugLevel !=0 } {
    print "RAMD DEBUG: $timeStep  > PREVIOUS LIGAND COM IS: $prevLigCOM"
    print "RAMD DEBUG: $timeStep  > PREVIOUS PROTEIN COM IS: $prevProtCOM"
    print "RAMD DEBUG: $timeStep  > CURRENT LIGAND COM IS: $currLigCOM"
    print "RAMD DEBUG: $timeStep  > CURRENT PROTEIN COM IS: $currProtCOM"
   }
  
   #***** Evaluate distance between ligand com and protein com
   set currDist [veclength [vecsub $currLigCOM $currProtCOM]]   
   #***** Evaluate the change in the distance between the protein and the ligand com during the ramd stint 
   set diffDist [expr "${currDist}-${prevDist}"]
   print "RAMD: $timeStep     >>> DISTANCE LIGAND COM - PPROTEIN COM IS: $currDist; IT CHANGED BY $diffDist"

   #***** Check if the ligand has exited the protein
   if { $currDist >= $maxDist } { set exitFlag 1; return }
   
   #***** Compute the distance travelled by the ligand com during the ramd stint
   evalWalkDist $timeStep $prevLigCOM $prevProtCOM $currLigCOM $currProtCOM
 
   #***** Evaluate whether a new force direction will be generated
   if { $walkDist <= $rMinRamd } {    

    genRandAccel $timeStep

    print "RAMD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (< $rMinRamd)"
    print "RAMD: $timeStep     >>> CONTINUE WITH $ramdSteps STEPS OF RAMD SIMULATION"
    print "RAMD: $timeStep     >>> CHANGE ACCELERATION DIRECTION TO: $r; ||r|| = $lenr"
    
    #***** Reset the ramd step count
    #***** Reset the positions of the ligand and protein COMs
    set ramdStep 0; set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"
 
    #***** Increment time step and go to the next time step right now
    incr timeStep
    return
  
   } elseif { $walkDist > $rMinRamd && $mdSteps == 0 } {
    print "RAMD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (> $rMinRamd)"
    print "RAMD: $timeStep     >>> CONTINUE WITH $ramdSteps STEPS OF RAMD SIMULATION"
    print "RAMD: $timeStep     >>> KEEP PREVIOUS ACCELERATION DIRECTION: $r; ||r|| = $lenr"
   
    #***** Reset the positions of the ligand and protein COMs
    set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"
 
    #***** Increment time step and go to the next time step right now
    incr timeStep
    return
  
   } elseif { $walkDist > $rMinRamd && $mdSteps != 0 } {
  
    print "RAMD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (> $rMinRamd)"
    print "RAMD: $timeStep     >>> SWITCH TO $mdSteps STEPS OF STANDARD MD SIMULATION"
    
    #***** Reset the flag values
    set ramdStep 0; set ramdFlag 0; set mdFlag 1
  
    #***** Reset the positions of the ligand and protein COMs
    set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

    #***** Increment time step and go to the next time step right now
    incr timeStep
    return
  
   }   
  
   #***** Ensure that the positions of the ligand and protein COMs are reset after the evaluation
   set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

   #***** Increment time step and go to the next time step right now
   incr timeStep
   return
  
  }

  #***** Unset the force
  if { [info exists force] } { unset force }
  if { [info exists fx] } { unset fx }
  if { [info exists fy] } { unset fy }
  if { [info exists fz] } { unset fz }
  if { [info exists f] } { unset f }
  if { [info exists totLigForce] } { unset totLigForce }
  if { [info exists extLigForce] } { unset extLigForce }
  if { [info exists totLigForceV] } { unset totLigForceV }
  if { [info exists extLigForceV] } { unset extLigForceV }
  
  #***** Increment time step and go to the next time step right now
  incr timeStep
  return
   
 }


 #***** Perform standard md simulation 
 if { $timeStep != 0 && $ramdFlag == 0 && $mdFlag == 1 && $exitFlag == 0 } {

  #***** Count the md steps
  incr mdStep
  #***** If debug level increased, check for the definition of the force
  if { $debugLevel != 0 } {
   print "MD DEBUG: TIMESTEP IS: $timeStep; RAMD STEP IS: $ramdStep; MD STEP IS: $mdStep"
   if { ! [ info exists force ] } { 
    print "MD DEBUG: EXTERNAL FORCE NOT DEFINED" 
   } else {
    print "MD DEBUG: WARNING: EXTERNAL FORCE VECTOR DEFINED: $force" 
   }
  }
    
  #***** Set flag for writing force output
  if { $forceOutFreq != 0 } { set outputFlag [expr "$mdStep % $forceOutFreq"] }

  #***** Write force output every 'forceOutFreq' steps 
  if { ! [ array exists extForces ] && [ array exists totForces ] && $outputFlag == 0 } {

   #***** Write the positions of the protein and ligand COMs
   print "MD FORCE: $timeStep  > LIGAND COM is: $ligCOM"
   print "MD FORCE: $timeStep  > PROTEIN COM IS $protCOM"

   #***** Calculate total forces acting on the ligand com
   set totLigForceX 0; set totLigForceY 0; set totLigForceZ 0; set totLigForceV 0
   foreach ramdAtom $ramdAtoms {
    set atomMass "$masses($ramdAtom)"
    set totAtomForce "$totForces($ramdAtom)"
    set totAtomForceX [ lindex $totAtomForce 0 ]
    set totAtomForceY [ lindex $totAtomForce 1 ]
    set totAtomForceZ [ lindex $totAtomForce 2 ]
    set totAtomForceV [ veclength "$totAtomForce" ]
    set totLigForceX [ expr "$totLigForceX + $totAtomForceX" ]
    set totLigForceY [ expr "$totLigForceY + $totAtomForceY" ]
    set totLigForceZ [ expr "$totLigForceZ + $totAtomForceZ" ]
    set totLigForceV [ expr "$totLigForceV + $totAtomForceV" ]
    if { $debugLevel != 0 } { 
     print "MD DEBUG: ATOM $ramdAtom: MASS $atomMass: TOT FORCE $totAtomForceV"
    }    
    unset atomMass
    unset totAtomForceX; unset totAtomForceY; unset totAtomForceZ; unset totAtomForceV; unset totAtomForce
   }
   set totLigForce "$totLigForceX $totLigForceY $totLigForceZ"
  
   #***** Write total forces acting on the ligand COM
   print "MD FORCE: $timeStep  > TOTAL FORCE ON THE LIGAND COM IS: $totLigForce ($totLigForceV)"
  
  } elseif { [ array exists extForces ] && [ array exists totForces ] && $outputFlag == 0  } {
  
   error "MD: $timeStep  > ERROR: EXTERNAL FORCES PRESENT DURING MD STEP: EXIT NAMD"
  
  }
   
  #***** Evaluate standard md stint
  set evalMdFlag [expr "$mdStep % $mdSteps"]

  if { $evalMdFlag == 0 } {
  
   print "MD: $timeStep  ***** EVALUATE $mdSteps STANDARD MD STEPS AT TIMESTEP $timeStep *****"
  
   #***** com positions
   set currLigCOM $ligCOM; set currProtCOM $protCOM
   if { $debugLevel !=0 } {
    print "MD DEBUG: $timeStep  > PREVIOUS LIGAND COM IS: $prevLigCOM"
    print "MD DEBUG: $timeStep  > PREVIOUS PROTEIN COM IS: $prevProtCOM"
    print "MD DEBUG: $timeStep  > CURRENT LIGAND COM IS: $currLigCOM"
    print "MD DEBUG: $timeStep  > CURRENT PROTEIN COM IS: $currProtCOM"
   }
   
   #***** Evaluate distance between ligand com and protein com
   set currDist [veclength [vecsub $currLigCOM $currProtCOM]] 
   #***** Evaluate the change in the distance between the protein the ligand com during the md stint 
   set diffDist [expr "${currDist}-${prevDist}"]
   print "MD: $timeStep     >>> DISTANCE LIGAND COM - PPROTEIN COM IS: $currDist; IT CHANGED BY: $diffDist"

   #***** Check if the ligand has exited the protein
   if { $currDist >= $maxDist } { set exitFlag 1; return }

   #***** Compute the distance travelled by the ligand com during the md stint
   evalWalkDist $timeStep $prevLigCOM $prevProtCOM $currLigCOM $currProtCOM
 
   #***** Evaluate whether a switch to ramd is required
   if { $walkDist <= $rMinMd && $diffDist > 0 } {

    print "MD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (< $rMinMd)"
    print "MD: $timeStep     >>> THE DISTANCE BETWEEN THE LIGAND COM AND THE PROTEIN COM INCREASED"
    print "MD: $timeStep     >>> SWITCH TO $ramdSteps STEPS OF RAMD SIMULATION"
    print "MD: $timeStep     >>> KEEP PREVIOUS ACCELERATION DIRECTION: $r; ||r|| = $lenr"
    
    #***** Reset the flags
    set mdStep 0; set ramdFlag 1; set mdFlag 0
    
    #***** Reset the positions of the ligand and protein COMs
    set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

    #***** Increment time step and go to the next time step right now
    incr timeStep
    return
   
   } elseif { $walkDist <= $rMinMd && $diffDist < 0 } {

    genRandAccel $timeStep

    print "MD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (< $rMinMd)"
    print "MD: $timeStep     >>> THE DISTANCE BETWEEN THE LIGAND COM AND THE PROTEIN COM DECREASED"
    print "MD: $timeStep     >>> SWITCH TO $ramdSteps STEPS OF RAMD SIMULATION"
    print "MD: $timeStep     >>> CHANGE ACCELERATION DIRECTION TO: $r; ||r|| = $lenr"
    
    #***** Reset the flags
    set mdStep 0; set ramdFlag 1; set mdFlag 0
   
    #***** Reset the positions of the ligand and protein COMs
    set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

    #***** Increment time step and go to the next time step right now
    incr timeStep
    return
   
   } elseif { $walkDist > $rMinMd && $diffDist > 0 } {
      
    print "MD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (> $rMinMd)"
    print "MD: $timeStep     >>> THE DISTANCE BETWEEN THE LIGAND COM AND THE PROTEIN COM INCREASED"
    print "MD: $timeStep     >>> CONTINUE WITH $mdSteps STEPS OF STANDARD MD SIMULATION"
   
    #***** Reset the positions of the ligand and protein COMs
    set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

    #***** Increment time step and go to the next time step right now
    incr timeStep
    return
   
   } elseif { $walkDist > $rMinMd && $diffDist < 0 } {
      
    genRandAccel $timeStep

    print "MD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (> $rMinMd)"
    print "MD: $timeStep     >>> THE DISTANCE BETWEEN THE LIGAND COM AND THE PROTEIN COM DECREASED"
    print "MD: $timeStep     >>> SWITCH TO $ramdSteps STEPS OF RAMD SIMULATION"
    print "MD: $timeStep     >>> CHANGE ACCELERATION DIRECTION TO: $r; ||r|| = $lenr"

    set mdStep 0; set ramdFlag 1; set mdFlag 0
   
    #***** Reset the positions of the ligand and protein COMs
    set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

    #***** Increment time step and go to the next time step right now
    incr timeStep
    return
   
   }
  
   #***** Ensure the positions of the ligand and protein COMs are reset
   set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

   #***** Increment time step and go to the next time step right now
   incr timeStep
   return
   
  }
 
  #***** Unset the force
  if { [info exists force] } { unset force }
  if { [info exists fx] } { unset fx }
  if { [info exists fy] } { unset fy }
  if { [info exists fz] } { unset fz }
  if { [info exists f] } { unset f }
  if { [info exists totLigForce] } { unset totLigForce }
  if { [info exists extLigForce] } { unset extLigForce }
  if { [info exists totLigForceV] } { unset totLigForceV }
  if { [info exists extLigForceV] } { unset extLigForceV }

  #***** Increment time step and go to the next time step right now
  incr timeStep
  return
    
 }
 
 #***** Check that during all timesteps at least one of the flags 'ramdFlag' and 'mdFlag' is activated
 #***** Exit with error if that's not the case
 if { $timeStep != 0 && $ramdFlag == 0 && $mdFlag == 0 && $exitFlag == 0  } { 
  error "RAMD: $timeStep  > ERROR: NEITHER THE RAMD NOR THE MD FLAG IS ACTIVATED; EXIT NAMD"    
 }

 #***** Check that during all timesteps the flags 'ramdFlag' and 'mdFlag' are not simultaneously activated
 #***** Exit with error if that's not the case
 if { $timeStep != 0 && $ramdFlag == 1 && $mdFlag == 1 && $exitFlag == 0  } { 
  error "RAMD: $timeStep  > ERROR: BOTH THE RAMD AND MD FLAGS ARE ACTIVATED; EXIT NAMD" 
 }

 #***** Increment time step and go to the next time step right now
 incr timeStep
 return

} ;# namespace
} ;# proc calcforces {}
 
#*****************
# END
#*****************
