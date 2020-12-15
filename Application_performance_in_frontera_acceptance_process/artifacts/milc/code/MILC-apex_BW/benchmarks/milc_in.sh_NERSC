#!/bin/bash
#
# The follwoing need to be defined before calling run_milc
# nx, ny, nz, nt - problem dimensions
# warms - # of warmups to peform
# trajecs - # of trajectories to run
# checklat - lattice checkpoint file name
# save_cmd - set to "forget" or "save_serial $checklat"
# reload_cmd - set to "reload_parallel $checklat" or "continue"
#
# choose which command to execute by setting run_type before calling
function run_milc() {
    if [ "$run_type" == "srun" ]; then
      command="srun -n $N --ntasks-per-socket=$S ./su3_rhmd_hisq"
    elif [ "$run_type" == "aprun" ]; then
      command="aprun -n $N -S $S -d $threads_per_rank -cc numa_node ./su3_rhmd_hisq"
    else
      command="mpirun -n $N -env I_MPI_PIN_DOMAIN socket -env OMP_NUM_THREADS=$threads_per_rank ./su3_rhmd_hisq"
    fi
    echo "$0: Running \"${command}\" "
    if [ "$debug" == "true" ]; then
      command="cat ";
      echo "Debugging script with cat"
    fi
    $command <<EOI
    prompt 0
    nx $nx
    ny $ny
    nz $nz
    nt $nt
    iseed 4563421
    n_pseudo 5
    load_rhmc_params rationals_m001m05m5.test1
    beta 6.3
    n_dyn_masses 3
    dyn_mass .001 .05 .5
    dyn_flavors 2 1 1
    u0  0.890
 
    warms $warms
    trajecs 0
    traj_between_meas 40
    microcanonical_time_step .05
    steps_per_trajectory 20
    cgresid_md_fa_gr 1 1 1
    max_multicg_md_fa_gr  5  5  5
    cgprec_md_fa_gr  2 2 2
    cgresid_md_fa_gr 1 1 1
    max_multicg_md_fa_gr  5  5  5
    cgprec_md_fa_gr  2 2 2
    cgresid_md_fa_gr 1 1 1
    max_multicg_md_fa_gr  5  5  5
    cgprec_md_fa_gr  2 2 2
    cgresid_md_fa_gr 1 1 1
    max_multicg_md_fa_gr  5  5  5
    cgprec_md_fa_gr  2 2 2
    cgresid_md_fa_gr 1 1 1
    max_multicg_md_fa_gr  5  5  5
    cgprec_md_fa_gr  2 2 2
    prec_ff 2
    number_of_pbp_masses 0
    fresh
    $save_cmd

    warms 0
    trajecs $trajecs
    traj_between_meas 2
    microcanonical_time_step .0333
    steps_per_trajectory 30
    cgresid_md_fa_gr 1e-5 1e-5 1e-5
    max_multicg_md_fa_gr  2500  2500  2500
    cgprec_md_fa_gr  2 2 2
    cgresid_md_fa_gr 1e-7 1e-7 1e-7
    max_multicg_md_fa_gr  2500  2500  2500
    cgprec_md_fa_gr  2 2 2
    cgresid_md_fa_gr 1e-7 1e-7 1e-7
    max_multicg_md_fa_gr  2500  2500  2500
    cgprec_md_fa_gr  2 2 2
    cgresid_md_fa_gr 1e-7 1e-7 1e-7
    max_multicg_md_fa_gr  2500  2500  2500
    cgprec_md_fa_gr  2 2 2
    cgresid_md_fa_gr 1e-7 1e-7 1e-7
    max_multicg_md_fa_gr  2500  2500  2500
    cgprec_md_fa_gr  2 2 2
    prec_ff 2
    number_of_pbp_masses 3
    max_cg_prop 500
    max_cg_prop_restarts 1
    npbp_reps 1
    prec_pbp 2
    mass 0.01
    naik_term_epsilon 0
    error_for_propagator 1e-8
    rel_error_for_propagator 0
    mass 0.05
    naik_term_epsilon 0
    error_for_propagator 1e-8
    rel_error_for_propagator 0
    mass 0.5
    naik_term_epsilon 0
    error_for_propagator 1e-8
    rel_error_for_propagator 0
    $reload_cmd
    forget
EOI
}

