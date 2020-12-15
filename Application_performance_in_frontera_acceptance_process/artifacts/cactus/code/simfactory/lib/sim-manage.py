# sim-manage -- manage simulations.
#
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
#

import sys, os, time

optionGroups = ['sim-manage']
Purpose = "manage a Cactus process"

############################################

from libutil import *
import simenv,simlib,restartlib,simrestart

## lets go to work

known_commands = ['create', 'submit', 'run', 'interactive', 'stop', 'archive', 'cleanup', 'get-output-dir',
                'show-output', 'create-submit', 'create-run', 'purge', 'run-debug']
usage_strings = {
                'create': 'create a simulation', 
                'submit': 'submit a simulation', 
                'cleanup': 'clean up one or more simulations',
                'create-submit': 'create and submit a simulation',
                'create-run': 'create and run a simulation',
                'run': 'launch a simulation', 
                'interactive': 'initiate an interactive session on a compute node',
                'stop': 'stop an active simulation',
                'show-output': 'show the output of a simulation',
                'archive': 'archive a simulation or an individual restart',
                'run-debug': 'launch simulation using an interactive debugger',
                'get-output-dir': 'get the output directory for a given simulation/restart',
                'purge': 'purge a simulation'}

############################### MAIN ###############################

def command_purge():
    if len(simenv.OptionsManager.args) == 0:
        fatal("Error: no simulation name specified\nUsage: %s purge <simulationname> <simulationname> ..." % sys.argv[0])
        
    for simulationName in simenv.OptionsManager.args:
        simulationName = simenv.OptionsManager.args.pop(0)
    
        display("Simulation name: %s" % simulationName)
        
        restart = simrestart.SimRestart()
        
        ret = restart.load(simulationName)
        
        if ret == -1:
            fatal("unable to load simulation %s for purging" % simulationName)

        restart.trash()
        
        restart.done()

def command_cleanup():
    if len(simenv.OptionsManager.args):
        simulations = simenv.OptionsManager.args
    else:
        simulations = simlib.GetSimulations()

    for simulationName in simulations:
        display("Cleaning up simulation %s" % simulationName)

        restart = simrestart.SimRestart()
        ret = restart.load(simulationName)

        if simenv.OptionsManager.HasOption('restart-id'):
            rid = simenv.OptionsManager.GetOption('restart-id')
        else:
            rid = restartlib.GetActiveRestartId(restart)
        if rid<0:
            info("Simulation is not active; doing nothing")
            continue
        
        ret = restart.load(simulationName, rid)
        if ret == -1:
            fatal("unable to load simulation %s for cleanup" % simulationName)
                
        info("Cleaning up simulation %s, restart id %s" % (simulationName, rid))
        restart.finish()
        restart.done()

def command_interactive():
    restart = simrestart.SimRestart()
    restart.interactive()

simenv.CommandDoc['create'] = """

Usage
-----

sim create [options] *simname* [ *parfile* | :option:`--parfile` *parfile* ]

  Create a simulation named *simname* using a specific parameter file.

sim create [options] *simname* --testsuite [--select-tests *testspec* ]

  Create a simulation named *simname* which runs the Cactus test
  suite.  See the documentation for selecting tests.

"""

simenv.CommandOptions['create'] = ["parfile", "configuration", "config", "testsuite", "select-tests"]
        
def command_create():
    if simenv.CACTUS_PATH == None:
        fatal("cannot proceed with unknown CACTUS_PATH")
    
    simulationName = simlib.ParseSimulationCommandLine()
    
    if simulationName == None:
        if len(simenv.OptionsManager.args) == 0:
            fatal("no simulation name specified")
            
        simulationName = simenv.OptionsManager.args.pop(0)
        display("Simulation name: %s" % simulationName)

    if not simenv.OptionsManager.HasOption('testsuite'):
        parfile = simlib.GetParFile()
        if parfile == None:
            fatal("no parameter file specified")
    
        display("Parameter file: %s" % parfile)
    else:
        parfile = simlib.GetParFile()
        if parfile != None:
            fatal("parameter file specified together with --testsuite option")

        parfile = "" # Properties system writes 'None' as an empty
                     # string to file, so we have to use something
                     # else

    restart = simrestart.SimRestart()
    restart.create(simulationName, parfile)
    restart.done()

simenv.CommandDoc['create-submit'] = """
This command combines the action of create and submit to perform the two
operations together.
"""

def command_create_submit():
    
    simulationName = simlib.ParseSimulationCommandLine()
    
    if simulationName == None:
        simulationName = simenv.OptionsManager.args.pop(0)
        display("Simulation name: %s" % simulationName)
    else:
        simenv.OptionsManager.args.insert(0, simulationName)
        
    # TODO: don't do this; there userSubmit which checks whether a
    # simulation needs to be created, handle this in the same manner
    # to clean this up
    command_create()
    
    #reinsert our simulationname since command_create popped it off.
    # TODO: do this differently; manage the simulation name handling
    # differently
    simenv.OptionsManager.args.insert(0, simulationName)
    
    command_submit()

def command_create_run():
    
    simulationName = simlib.ParseSimulationCommandLine()
    
    if simulationName == None:
        simulationName = simenv.OptionsManager.args[0]
        display("Simulation name: %s" % simulationName)
    else:
        simenv.OptionsManager.args.insert(0, simulationName)
        
    command_create()
    
    #reinsert our simulationname since command_create popped it off.
    simenv.OptionsManager.args.insert(0, simulationName)
    
    command_run()
    
simenv.CommandDoc['run'] = """

Usage
-----

sim run [options] *simname*

  This command runs an existing simulation interactively; i.e. without
  submitting it to a queuing system.  This can be useful with
  interactive jobs on HPC systems and on workstations and laptops.
  Running codes directly on the head node of an HPC system is usually
  forbidden, and will not usually work.

"""

simenv.CommandOptions['run'] = ["recover", "debugger", "restart-id", "from-restart-id"]

def command_run():
    simulationName = simlib.ParseSimulationCommandLine()
    
    if simulationName == None:
        simulationName = simenv.OptionsManager.args.pop(0)
    
    display("Simulation name: %s" % simulationName)
    
    restart = simrestart.SimRestart()
    ret = restart.load(simulationName)
    
    if ret == -1:
        fatal("unable to load simulation %s for execution" % simulationName)
    
    if simenv.OptionsManager.HasOption('restart-id'):
        restart_id = simenv.OptionsManager.GetOption('restart-id')
        restart.submitRun(simulationName, restart_id)
    else:
        restart.userRun(simulationName)
        restart.done()
    
def command_run_debug():
    simulationName = simlib.ParseSimulationCommandLine()
    
    if simulationName == None:
        simulationName = simenv.OptionsManager.args.pop(0)
    
    display("Simulation name: %s" % simulationName)
    
    restart = simrestart.SimRestart()
    
    ret = restart.load(simulationName)
    
    if ret == -1:
        fatal("unable to load simulation %s for execution" % simulationName)
                
    restart.userRun(simulationName, debug=True)
    restart.done()
    
simenv.CommandDoc['submit'] = """

Usage
-----

sim submit [options] *simname* [options]

  Create a new restart in the simulation and submit it to the queuing
  system of the machine.

"""

simenv.CommandOptions['submit'] = ["recover", "debugger", "restart-id", "job-id", "from-restart-id", "hide", "hide-boring", "hide-dangerous"]

def command_submit():
    simulationName = simlib.ParseSimulationCommandLine()
    
    if simulationName == None:
        
        if len(simenv.OptionsManager.args) == 0:
            fatal("No simulation name specified")
            
        simulationName = simenv.OptionsManager.args.pop(0)
        display("Simulation name: %s" % simulationName)

    restart = simrestart.SimRestart()
    restart.userSubmit(simulationName)
    restart.done()
    

def command_get_output_dir():

    simulationName = None

    if simulationName == None:
        simulationName = simenv.OptionsManager.args.pop(0)

    restart = simrestart.SimRestart()

    restart_id = -1

    if simenv.OptionsManager.HasOption('restart-id'):
        restart_id = simenv.OptionsManager.GetOption('restart-id')
    else:
        restart_id = restartlib.GetMaxRestartID(simulationName)

    restart.load(simulationName, restart_id)

    display(restart.RestartDir)
    
    restart.done()

def command_show_output():
    simulationName = simenv.OptionsManager.args.pop(0)
    
    display("Simulation name: %s" % simulationName)
    
    restart_id = -1

    if simenv.OptionsManager.HasOption('restart-id'):
        restart_id = simenv.OptionsManager.GetOption('restart-id')
    else:
        restart_id = restartlib.GetMaxRestartID(simulationName)
        
    restart = simrestart.SimRestart()
    ret = restart.load(simulationName, restart_id)
    
    if ret == -1:
        fatal("unable to load simulation %s for output" % simulationName)
        
    restart.show_output()
    
    restart.done()
    
def command_stop():
    
    for simulationName in simenv.OptionsManager.args:
        simulationName = simenv.OptionsManager.args.pop(0)
        
        display("Simulation name: %s" % simulationName)
    
        if simenv.OptionsManager.HasOption('restart-id'):
            restart_id = simenv.OptionsManager.GetOption('restart-id')
            restart = restartlib.GetRestartByRestartId(simulationName, restart_id)
            restart.stop()
            restart.done()
            return
    
        if simenv.OptionsManager.HasOption('job-id'):
            job_id = simenv.OptionsManager.GetOption('job-id')
            info("Have job_id: %s" % job_id)
            restart = restartlib.GetRestartByJobId(simulationName, job_id)
            restart.stop()
        
            # stop does not call finish automatically, because finish() may call stop().
            restart.finish()
            restart.done()
            return
    
        restartlib.StopAllActiveRestarts(simulationName)

def command_archive():
    simulationName = simenv.OptionsManager.args.pop(0)
    
    display("Simulation name: %s" % simulationName)
    
    if simenv.OptionsManager.HasOption('restart-id'):
        restart_id = simenv.OptionsManager.GetOption('restart-id')
        restart = restartlib.GetRestartByRestartId(simulationName, restart_id)
    else:
        restart = simrestart.SimRestart()
        ret = restart.load(simulationName)
    
        if ret < 0:
            fatal("unable to load simulation %s for archiving" % simulationName)
    
    restart.archive()
    restart.done()
    
def CommandDispatch():
    global known_commands
    
    if simenv.COMMAND == None:
        command = simenv.OptionsManager.args.pop(0)
    else:
        command = simenv.COMMAND
    
    if command not in known_commands:
        display("Error: unknown command %s" % command)
        simenv.OptionsManager.PrintHelp()
        sys.exit(0)
    
    command = command.replace("-", "_")
    
    info("Executing command: %s" % command)
    exec("command_%s()" % command)

def main():

    simlib.RequireMachine()

    cactusDir = simenv.CACTUS_PATH
    info("Cactus Directory: %s" % cactusDir)
    
    if os.getcwd() != cactusDir:
        warning("Current Working directory does not match Cactus sourcetree, changing to %s" % cactusDir)
        # make sure we're in the Cactus source directory, otherwise all of this will blow way up.
        # os.chdir(cactusDir)
    
    ############################################

    info("simenv.COMMAND: %s" % simenv.COMMAND)
    if len(simenv.OptionsManager.args) < 2 and simenv.COMMAND == None:
        simenv.OptionsManager.PrintHelp()
        sys.exit(0)

    CommandDispatch()
