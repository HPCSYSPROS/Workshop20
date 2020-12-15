# sim -- central dispatch for Sim Factory
#
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
#

import sys, os, shutil

############################################
## APP SPECIFIC DEFINES ##
#
# Usage            when (-h, --help) is used, this is the first line showing how to call the app
# optionGroups    which option groups from etc/options to import. common is always imported. 

Usage = "usage: sim [args] command"

# import all option groups.

optionGroups = ['sim-build', 'sim-info', 'sim-manage', 'sim-sync', 'sim-util']


############################################

############################################
## INIT ##

rp = os.path.realpath(__file__)
paths = rp.split(os.sep)[:-2]

global BASE_PATH
BASE_PATH = os.sep.join(paths)

LIB_PATH  = os.path.join(BASE_PATH, "lib")

sys.path.append(BASE_PATH)
sys.path.append(LIB_PATH)

import libutil
import simenv,simlib,simremote,restartlib
from libutil import *

global known_commands
global associated_binaries

known_commands = {'sync':'sim-sync', 'build':'sim-build'}
known_binaries = ['sim-sync', 'sim-build', 'sim-manage', 'sim-info', 'sim-util']
# import all of our known binaries

for n in known_binaries:
    __import__(n)


############################################

def buildKnownCommands():
    global known_commands
    for binary in known_binaries:
        cmds = sys.modules[binary].known_commands
        
        for c in cmds:
            known_commands[c] = binary

buildKnownCommands()

def buildUsageString():
    global Usage
    
    cmdString = "%s\n\n\nAvailable Commands:\n" % Usage
    iString = []
    
    commands = known_commands.keys()
    commands.sort()
    for cmd in commands:
        binary = known_commands[cmd]
        if sys.modules[binary].usage_strings.has_key(cmd):
            usage = sys.modules[binary].usage_strings[cmd]
        else:
            usage = ""
        ss = "  %-30s%s" % (cmd, usage)
        iString.append(ss)
    
    cmdString = "%s%s" % (cmdString, "\n".join(iString))
    return cmdString
        

def CommandDispatch():
    global known_binaries
    
    command = simenv.OptionsManager.args[0]
    
    if command == 'help':
        simenv.OptionsManager.PrintHelp()
        sys.exit(0)
        
    # set our command in the environment
    simenv.COMMAND = command
    simenv.OptionsManager.args.pop(0)

    if command not in known_commands:
        fatal("unknown command %s" % command)
    
    binary = known_commands[command]
    module = sys.modules[binary]
    
    module.main()

def RemoteExecute(machineName):
    simlib.RequireMachine()

    RemoteEnvironment = simremote.RemoteEnvironment()
    RemoteEnvironment.init(machineName)
    RemoteEnvironment.ExecuteSameCommand(parrotArguments=True, stripArguments=['remote', 'remotecactuspath', 'machine'])
        
def main():
    # a small bit of error checking.
    if len(simenv.OptionsManager.args) == 0:
        simenv.OptionsManager.PrintHelp()
        sys.exit(0)

    ############################################
    
    info("The Simulation Factory: Manage Cactus simulations\n")
    
    info("defs: %s" % simenv.cdb)
    
    if simenv.udb != None:
        info("defs.local: %s" % simenv.udb)
        
        # TODO: Cleaning up all simulations here can lead to fatal
        # errors if a problem with a particular restart is detected.
        # This makes Simfactory then completely ususable. Therefore
        # (and for several other reasons discussed on the mailing
        # list), cleaning up simulations here is disabled.
        #normalize restarts.
        #restartlib.CleanupRestarts()
    
    if simenv.OptionsManager.HasOption("remote"):
        remoteMachine = simenv.OptionsManager.GetOption("remote")
        info("Executing command for remote machine: %s" % remoteMachine)
        RemoteExecute(remoteMachine)
        return

    CommandDispatch()

if __name__ == '__main__':
    simenv.init(base=BASE_PATH, callingExecutable=os.path.join(BASE_PATH, "bin/sim"), usageString=buildUsageString(), optionGroups=optionGroups)
    main()
