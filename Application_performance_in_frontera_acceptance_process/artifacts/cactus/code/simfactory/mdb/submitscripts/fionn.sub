#! /bin/bash
#PBS -l nodes=@NODES@:ppn=@PPN@
#PBS -l walltime=@WALLTIME@
#PBS -N @SHORT_SIMULATION_NAME@
#PBS -A @ALLOCATION@
#PBS -r n
#PBS -m bea
#PBS -M @EMAIL@
#PBS @('@CHAINED_JOB_ID@' != '' ? '-W depend=afterany:@CHAINED_JOB_ID@' : '')@
#PBS -o @RUNDIR@/@SIMULATION_NAME@.out
#PBS -e @RUNDIR@/@SIMULATION_NAME@.err

cd @SOURCEDIR@
@SIMFACTORY@ run @SIMULATION_NAME@ --machine=fionn --restart-id=@RESTART_ID@ @FROM_RESTART_COMMAND@
