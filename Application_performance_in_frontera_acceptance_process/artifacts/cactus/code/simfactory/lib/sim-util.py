# sim-util -- various utility commands for Cactus/simfactory
#
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
#

import sys, os

optionGroups = ['sim-build']
Purpose = "Simfactory utility commands"

############################################
import simdt
import simarchive
import simenv,simlib,simsubs

from libutil import *

## FUNCTIONS ##

known_commands = ['checkout', 'login', 'get-archived-simulation', 'setup', 'setup-silent', 'execute', 'remove-submitscript']
usage_strings = {
                'login': 'launch an interactive shell on a remote machine',
                'checkout': 'checkout thorns using GetComponents',
                'get-archived-simulation': 'retrieve an archived simulation by simulationid',
                'setup': 'interactively configure simfactory',
                'setup-silent': 'silently create a set of reasonable defaults',
                'remove-submitscript': 'remove a submit script from a configuration',
                'execute': 'execute a command'}
                
def CommandDispatch():
    global known_commands
    
    if simenv.COMMAND == None:
        command = simenv.OptionsManager.args.pop(0)
    else:
        command = simenv.COMMAND
        
        if len(simenv.OptionsManager.args) > 0:
            if command == simenv.OptionsManager.args[0]:
                simenv.OptionsManager.args.pop()
            
    if command not in known_commands:
        fatal("unknown command %s" % command)
        
    info("Executing command: %s" % command)
    
    command = command.replace("-", "_")
    
    exec("command_%s()" % command)

def command_remove_submitscript():
    simlib.RequireMachine()
    
    if len(simenv.OptionsManager.args) == 0:
        fatal("No configuration specified")
    
    configuration = simenv.OptionsManager.args.pop(0)
    
    config_path = simlib.BuildPath(simenv.CONFIGS_PATH, configuration)
    
    if not(os.path.exists(config_path)):
        print "path: %s" % config_path
        fatal("Configuration \"%s\" does not exist" % configuration)
    
    submit_path = simlib.BuildPath(config_path, 'SubmitScript')
    
    if not(os.path.exists(submit_path)):
        fatal("Configuration \"%s\" does not have a submit script" % configuration)

    display("Removing SubmitScript for configuration: %s" % configuration)
    
    os.unlink(submit_path)
    
def command_execute():
    simlib.RequireMachine()
    
    cmd = simenv.OptionsManager.args.pop(0)
    simlib.ExecuteCommand(cmd)
    
def command_setup():
    # suppress warnings.
    simenv.VERBOSE = False
    
    try:
        dt = simdt.DecisionTree("defaults")
        dt.setupTree()
        dt.begin()
    except KeyboardInterrupt:
        fatal("\nSetup terminated.")

def command_setup_silent():
    # suppress warnings.
    simenv.VERBOSE = False
    
    try:
        dt = simdt.DecisionTree("defaults-silent")
        dt.setupTree()
        dt.begin()
    except KeyboardInterrupt:
        fatal("\nSetup terminated.")

simenv.CommandDoc['login'] = """

Usage
-----

sim login *machine*

  This command uses the SimFactory machine database to start an
  interactive terminal login session on the remote machine, typically
  using ssh.  This is helpful when the machine cannot be accessed
  directly, or requires a special login process.
"""

simenv.CommandOptions['login'] = []
        
def command_login():
    simlib.RequireMachine()
    
    if len(simenv.OptionsManager.args) == 0:
        fatal("Error: missing machine name\nUsage: %s login <machinename>" % sys.argv[0])
        
    DefineDatabase = simsubs.DefineDatabase()
    
    machineName = simenv.OptionsManager.args.pop()
    
    info("Logging into remote machine %s" % machineName)
    
    machineEntry = simenv.ConfigurationDatabase.GetMachine(machineName)
    localMachineEntry = simenv.LocalMachineEntry
    
    simlib.VerifyKeys(localMachineEntry, ['user', 'sourcebasedir'])
    simlib.VerifyKeys(machineEntry, ['user', 'sourcebasedir', 'basedir'])
    
    DefineDatabase.Set('USER', localMachineEntry.GetKey('user'))
    local_sourcebasedir = DefineDatabase.SubAll(localMachineEntry.GetKey('sourcebasedir'))
    source_name = simlib.GetDirSuffix(local_sourcebasedir)
    
    if source_name != None:
        local_sourcedir = simlib.BuildPath(local_sourcebasedir, source_name)
    
        DefineDatabase.Set('USER', machineEntry.GetKey('user'))
        sourcebasedir = DefineDatabase.SubAll(machineEntry.GetKey('sourcebasedir'))
        sourcedir = simlib.BuildPath(sourcebasedir, source_name)
        DefineDatabase.Set('SOURCEDIR', sourcedir)
    
        path = str()
    
        simulation = simenv.OptionsManager.GetOption('simulation')
    
        if simulation == None:
            msourcebasedir = DefineDatabase.SubAll(machineEntry.GetKey('sourcebasedir'))
            path = simlib.BuildPath(msourcebasedir, source_name)
        else:
            basedir = simenv.OptionsManager.GetOption('basedir')
            if basedir == None:
                basedir = DefineDatabase.SubAll(machineEntry.GetKey('basedir'))
        
            path = simlib.BuildPath(basedir, simulation)
    
        cmd = "cd %s || echo \"Could not change to directory %s\"; $SHELL -l" % (path, path)
        DefineDatabase.Set('SOURCEDIR', local_sourcedir)
    else:
        cmd = "$SHELL -l"
        
    cmd = simlib.GetSSHCommand(simenv.LocalMachine, machineName, cmd, "-t")
    # DefineDatabase.Set('USER', localMachineEntry.GetKey('user'))
    DefineDatabase.Set('USER', machineEntry.GetKey('user'))
    cmd = DefineDatabase.SubAll(cmd)
    localsshsetup = DefineDatabase.SubAll(machineEntry.localsshsetup)
    cmd = "{ %s; } && %s" % (localsshsetup, cmd)
    simlib.ExecuteCommand(cmd)
    
def command_get_archived_simulation():

    simlib.RequireMachine()
    
    machineEntry = simenv.LocalMachineEntry
        
    simlib.VerifyKeys(machineEntry, ['archivetype'])
        
    archiveType = machineEntry.archivetype
    
    ArchiveEngine = simarchive.SimArchive(simenv, archiveType)
    ArchiveEngine.authenticate()
    
    if len(simenv.OptionsManager.args) < 1:
        fatal("Error: wrong number of arguments\nUsage: %s get-archived-simulation [--restart-id=<id>] <simulationid> [<dstpath>]" % sys.argv[0])
        
    simulationid = simenv.OptionsManager.args[0]
    
    if len(simenv.OptionsManager.args) == 1:
        dstPath = "."
    else:
        dstPath = simenv.OptionsManager.args[1]

    display("Retrieving simulation: %s" % simulationid)
    display("Destination path: %s" % dstPath)
    
    if simenv.OptionsManager.HasOption('restart-id'):
        restart_id = simenv.OptionsManager.GetOption('restart-id')
    else:
        restart_id = None
    
    ArchiveEngine.get(simulationid, dstPath, restart_id)
    
def command_checkout():
    
    DefineDatabase = simsubs.DefineDatabase()
    
    simlib.RequireMachine()
    
    thornList = simenv.OptionsManager.args
    
    if len(thornList) == 0:
        fatal("Error: no thornlist(s) specified\nUsage: %s checkout <thornlist>" % sys.argv[0])
    
    thfile = "ThornList"
    
    getComponents = simlib.BuildPath(simenv.CACTUS_PATH, 'GetComponents')
    
    info("Using GetComponents: %s" % getComponents)
    
    if not(os.path.exists(getComponents)):
        fatal("could not find GetComponents at %s" % getComponents)
        
    for th in thornList:
        if th.count("://") > 0:
            local_name = th.split("/").pop()
            if os.path.exists(local_name):
                fatal("thornlist file %s already exists; refusing to overwrite" % local_name)
            
            ret = simlib.ExecuteCommand("curl -O %s" % th)
            
            if ret != 0:
                ret = simlib.ExecuteCommand("wget --no-check-certificate -O %s %s" % (local_name, th))
            
                if ret != 0:
                    fatal("Error: could not download thorn list <%s>" % th)
            
            th = local_name
        
        contents = simlib.GetFileContents(th)
        
        DefineDatabase.AddReplacement('!DEFINE ROOT', '.')
        
        contents = DefineDatabase.SubAll(contents)
        
        simlib.WriteContents(thfile, contents)
        
        cmd = "%s -a %s" % (getComponents, thfile)
        
        simlib.ExecuteCommand(cmd)
        os.unlink(thfile)

def main():

    # make sure our configs directory is created
    try:
        os.mkdir(simenv.CONFIGS_PATH)
    except OSError:
        pass

    CommandDispatch()
