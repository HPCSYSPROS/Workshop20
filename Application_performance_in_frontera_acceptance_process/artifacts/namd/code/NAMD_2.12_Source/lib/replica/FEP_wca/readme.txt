
FEP replica exchange instructions:

  Search "User_to_set" to see places you may need to change

  ====================================================================

  BEFORE RUNNING SIMULATION: 

        prepare your output directories:
        ../make_output_dirs.sh output_solv <NUMBER_OF_REPLICAS>

  ====================================================================

  TO RUN SIMULATION:

        User_to_set: be consistent with the number of replicas and name of output directory 
        mpiexec -n 36 $bindir/namd2 +replicas 36 fep_solv.conf --source FEP_wca.namd +stdout output_solv/%d/job0.%d.log

    mpirun -np 8 -hostfile hostfile $bindir/namd2 +replicas 8 fold_alanin.conf --source ../replica.namd +stdout output/%d/job0.%d.log
  the number of MPI ranks (-np) must be a multiple of the number of replicas (+replicas)
 
  ====================================================================

  AFTER SIMULATION:

      Post-process output with wham script:
           sh ./wham
           sh ./get_fe.sh

      Note: first time you post-process, the "post-processing" scripts
      will attempt to automatically compile several small binaries 
      in the wham directory needed for analysis. This automatic compile
      process therefore requires access to the GNU compiler,
      g++, in your PATH.
            i.e., on most Linux distros
               for bash shell:  export PATH=/usr/bin/g++:$PATH
               for c-based shells: setenv PATH "/usr/bin/g++ $PATH"

 
           Output files from post-processing includes:
             1. repu_fe.dat   <- repulsive energy
             2. disp_wham_fe  <- dispersive energy
             3. chrg_wham_fe  <- electrostatic energy

  ====================================================================
 
