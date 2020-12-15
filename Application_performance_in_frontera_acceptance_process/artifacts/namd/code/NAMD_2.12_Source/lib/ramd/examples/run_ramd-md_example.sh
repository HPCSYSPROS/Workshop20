#! /bin/bash

#NAMD_HOME=/path/to/your/directory ## please change this to the path of your NAMD installation
NAMD_HOME=/sw/mcm/app/vlad/namd/2009-08-13/namd2_mvapich2-1.4rc1/Linux-x86_64-intel-mvapich2
#NAMD_HOME=/sw/mcm/app/vlad/namd/2009-08-13/namd2_mvapich2-1.4rc1_tcl8.3/Linux-x86_64-intel-mvapich2
#NAMD_HOME=/sw/mcm/app/vlad/namd/2009-08-13/namd2_mvapich2-1.4rc1_tcl8.5/Linux-x86_64-intel-mvapich2
#NAMD_HOME=/sw/mcm/app/vlad/namd/2010-02-23/namd2_mvapich2-1.4rc1/Linux-x86_64-intel-mvapich2

echo "NAMD executable is : $NAMD_HOME/namd2"
echo "Linking information for the NAMD executable:"
ldd $NAMD_HOME/namd2

fname=dhaawt_DCL                      # the RAMD example; in order to run the RAMD-MD example, change this to "example_ramd-md"
outdir=ramd-md_output              # the RAMD output directory; change to example_ramd-md_output
if [ -d "$outdir" ]; then rm -rf $outdir; mkdir $outdir; else mkdir $outdir; fi

cat << EOF > $fname.ramd-md.namdin
#*** Example of the NAMD configuration file for running plain RAMD simulation
#*** Uses the 'Cornell et al. 1995' force field (topology file in AMBER format)
#*** Calculation should take about 15 min with the NAMD v2.7b on a single CPU

#**** AMBER force field ********************************************************
amber                               on
parmfile                            $fname.top 
ambercoor                           $fname.emd2.rst
readexclusions                     yes
exclude                      scaled1-4
1-4scaling                           0.83333333   #=1/1.2
scnb                                 2

#*** Approximations for nonbonded interactions **********************************
cutoff                              12
switching                           off
pairlistdist                        15
outputpairlists                   1000
stepspercycle                       10
nonbondedfreq                        1
fullelectfrequency                   1
margin                               1

#*** Timestep *******************************************************************
minimization                       off
numsteps                        500000
timestep                           2.0

#*** SHAKE use ******************************************************************
rigidbonds                         all
rigidTolerance                   1e-08

#*** LANGEVIN Temperature and pressure control (comment out this part if using Berendsen bath)
#*** Temperature control (Langevin) *********************************************
#temperature                        300
#langevin                            on
#langevintemp                       300.0
#langevindamping                      2.0
#langevinhydrogen                    on

#*** Pressure control (Langevin) ************************************************
#useGroupPressure                   yes
#useFlexibleCell                     no
#useConstantArea                     no
#LangevinPiston                      on
#LangevinPistonTarget                 1.01325
#LangevinPistonPeriod               100
#LangevinPistonDecay                 50
#LangevinPistonTemp                 300

#*** BERENDSEN Temperature and pressure coupling (comment out this part if using Langevin control) 
#*** Temperature coupling (Berendsen) ******************************************
temperature      300.0
tCouple          on
tCoupleTemp      300.0
tCoupleFile      $fname.tCouple.pdb   ### friction coef.=1.0
tCoupleCol       B

#*** Constant pressure (with Berendsen bath coupling) **************************
useGroupPressure                         yes
useFlexibleCell                          no
useConstantArea                          no
BerendsenPressure                        on
BerendsenPressureTarget                  1.0
BerendsenPressureCompressibility         0.0000457
BerendsenPressureRelaxationTime          100
BerendsenPressureFreq                    10

#*** PME and PBC ****************************************************************
PME		                    on
PMETolerance                     1e-06
PMEGridSpacing                       1
cellBasisVector1           66.7213612   0.0000000   0.0000000		
cellBasisVector2            0.0000000  73.1177434   0.0000000
cellBasisVector3            0.0000000   0.0000000  76.7255299  
cellOrigin                 36.3112348  42.2434630  34.9181096
wrapAll                            off

#*** Output *********************************************************************
outputname                    $outdir/$fname.ramd-md
outputenergies                     100
restartname                   $outdir/$fname.ramd-md.rst
restartfreq                        100  
dcdfile                       $outdir/$fname.ramd-md.cdcd
dcdfreq                            100 
veldcdfile                    $outdir/$fname.ramd-md.vdcd
veldcdfreq                        1000
binaryoutput                       off
binaryrestart                       on

#*** Random Acceleration Molecular Dynamics *************************************

source ../scripts/ramd-4.1.tcl                   
#*** sources the wrapper script ramd-4.1.tcl;
#*** please change the directory '../scripts/' to '$dir' ( the correct path );
#*** directory '$dir' should contain the scripts: ramd-4.1.tcl, ramd-4.1_script.tcl, and vectors.tcl

ramd debugLevel                       0   
#*** activates verbose output if set to something else than 0

ramd mdStart                        yes      
#*** specifies whether combined RAMD-MD simulation will start with MD or RAMD; 
#*** defaults to no; can be set to yes if an intial MD stint is desired

ramd ramdSteps                       50
#*** specifies the number of steps in 1 ramd stint; 
#*** defaults to 50
 
ramd mdSteps                        100
#*** specifies the number of steps in 1 standard md stint; 
#*** defaults to 0 (pure RAMD simulation)

ramd accel                            0.5  
#*** specifies the acceleration to be applied; 
#*** defaults to 0.25 kcal/mol*A*amu

ramd rMinRamd                         0.5  
#*** specifies the minimum distance to be travelled by the ligand in 1 ramd stint; 
#*** defaults to 0.01 Angstr

ramd rMinMd                           2.0   
#*** specifies the minimum distance in Angstr to be travelled by the ligand in 1 md stint; 
#*** required if mdStep is not 0; ignored if mdSteps is 0
 
ramd forceOutFreq                    50
#*** every 'forceOutFreq' steps detailed output of forces will be written; 
#*** defaults to 0 (no detailed output)

ramd maxDist                         30
#*** specifies the distance between the COMs of the ligand and the protein when the simulation is stopped
#*** defaults to 50 Angstr
 
ramd firstProtAtom                    1 
#*** specifies the index of the first protein atom
#*** defaults to 1 (assumes first atom in the system corresponds to first protein atom
 
ramd lastProtAtom                  4596 
#*** specifies the index of the last protein atom
#*** required; simulation exits if this parameter is not set

ramd firstRamdAtom                 4597
#*** specifies the index of the first ligand atom
#*** required; simulation exits if this parameter is not set

ramd lastRamdAtom                  4608
#*** specifies the index of the last ligand atom
#*** required; simulation exits if this parameter is not set

ramd ramdSeed                     14253
#*** specifies the seed for the random number generator (for the generation of acceleration directions)
#*** defaults to 14253 
#*** please change if you wish to run different trajectories
####################################################################################
EOF

$NAMD_HOME/namd2 $fname.ramd-md.namdin >> $outdir/$fname.ramd-md.namdout  # Command to run the example


