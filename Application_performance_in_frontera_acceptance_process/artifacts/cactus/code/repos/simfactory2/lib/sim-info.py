# sim-info -- information about simfactory
#
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
#

import sys, os, time

optionGroups = ['sim-info']
Purpose = "display some convienient information about simfactory"

############################################

import simarchive,simenv,simlib,restartlib,simrestart

from libutil import *

global known_commands
global known_help
known_commands = ['list-machines', 'print-mdb-entry', 'whoami', 'print-mdb', 'list-archived-simulations',
'list-configurations', 'list-conf', 'list-simulations', 'list-sim']
usage_strings = {'list-machines':'list all the machines in the machine database',
        'print-mdb-entry':'list information about a single machine',
        'whoami':'what is the name of the machine simfactory is running on',
        'print-mdb':'print parsable mdb',
        'list-configurations':'list simfactory cactus configurations',
        'list-simulations':'list simfactory cactus simulations',
        'list-archived-simulations': 'list archived simulations'}
        

############################### INFO FUNCTIONS ###############################


    
def list_archived_simulations():

    simlib.RequireMachine()
    
    machineEntry = simenv.LocalMachineEntry
        
    simlib.VerifyKeys(machineEntry, ['archivetype'])
        
    archiveType = machineEntry.archivetype
    
    ArchiveEngine = simarchive.SimArchive(archiveType)
    ArchiveEngine.authenticate()
    
    simulations = ArchiveEngine.listSimulations()

    if simenv.OptionsManager.GetOption('name-only'):
        for sim in simulations:
            print sim
        return
    
    for sim in simulations:
        display("")
        display("%s: " % sim['SimulationId'])
        display("        Archive Path: %s" % sim['StoredPath'])
        display("     Simulation Name: %s" % sim['SimulationName'])
        display("             Machine: %s" % sim['Machine'])
        display("            Hostname: %s" % sim['Hostname'])
        display("                User: %s" % sim['User'])
        display("                Date: %s" % sim['Date'])
        
        restarts = sim['Restarts']
        
        display("            Restarts: found %s restart(s)" % len(restarts))
        
        for r in restarts:
            display("                      %s" % r)
        
    
def list_sim():
    list_simulations()
def list_simulations():
    
    simlib.RequireMachine()
    
    if len(simenv.OptionsManager.args):
        simulations = simenv.OptionsManager.args
    else:
        simulations = simlib.GetSimulations()
    
    if len(simulations) == 0:
        display("There are no simulations")
        return
        
    simulations.sort()
    
    if simenv.VERBOSE:
        display("Simulations:")

    if simenv.OptionsManager.GetOption('name-only'):
        for sim in simulations:
            print sim
        return

    for sim in simulations:
        restart = simrestart.SimRestart()
        restart.load(sim)
        
        restartIds = restartlib.GetRestartIds(restart)
        activeId = restartlib.GetActiveRestartId(restart, warnOnly=True)
        if activeId == -2:
            display("   %-21s   [ERROR, restart %04d]" % (sim, int(activeId)))
            continue

        active = activeId >= 0
        if not active and len(restartIds) > 0:
            # Use the last restart if there is no active restart
            activeId = restartIds[len(restartIds)-1]
        
        if activeId >= 0:
            ret = restart.load(sim, activeId)
            if ret < 0:
                # Something is wrong with the active restart
                display("   %-21s   [ERROR, restart %04d]" % (sim, int(activeId)))
                continue
        job_id = -1
        if activeId >= 0:
            job_id = restart.GetJobId()
        
        # TODO: The number of presubmitted restarts should also be
        # output.

        if active:
            state_map = dict()
            state_map['H'] = 'ACTIVE (PRESUBMITTED)'
            state_map['Q'] = 'ACTIVE (QUEUED)'
            state_map['R'] = 'ACTIVE (RUNNING)'
            state_map['U'] = 'ACTIVE (FINISHED)'
            state_map['E'] = 'ERROR'
            job_status = restartlib.GetJobStatus(job_id)
            state = state_map[job_status]
            if not state:
                error("Illegal job status '%s'" % job_status)
        else:
            state = 'INACTIVE'
            
        if not simenv.OptionsManager.GetOption('long'):
            display("   %-21s   [%-9s, restart %04d, job id %s]" % (sim, state, int(activeId), job_id))
        else:
            display("   %-21s= %s" % (sim, state))

            # Active Restart Id
            display("        %-16s= %04d" % ("restart id", int(activeId)))

            # Job
            display("        %-16s= %s" % ("job id", job_id))
            
            # Disk usage
            fd = simenv.popen("du -sk %s | cut -f1" % restart.SimulationDir)
            size_in_k = int(fd.read().strip())
            fd.close()
            display("        %-16s= %.1f GByte" % ("Disk usage", size_in_k * 1024 / 1.0e+9))
            
            # Properties
            #PrintManyLeadingSpace(restart.Properties.toString(), 8)
        
        if restart:
            restart.done()

    
def list_machines():

    machines = simenv.ConfigurationDatabase.GetMachines()
    
    if simenv.OptionsManager.GetOption('name-only'):
        for machine in machines:
            print machine
        return

    for machine in machines:
        entry = simenv.ConfigurationDatabase.GetMachine(machine)
        display("   %-12s   %6d nodes,  %4d ppn   %-30s   %s" % (machine, int(entry.nodes), int(entry.ppn), entry.hostname, entry.status))

def print_mdb_entry():    
    if len(simenv.OptionsManager.args) == 0:
        fatal("no machine specified")
        
    for i in range(len(simenv.OptionsManager.args)):
        machineName = simenv.OptionsManager.args[i]
        simenv.ConfigurationDatabase.MachineParser.PrintSection(machineName)
        
def whoami():
    
    machine = simlib.GetMachineName()
    display("Current machine: %s" % machine)

def print_mdb():
    

    if len(simenv.OptionsManager.args) == 0:
        simenv.ConfigurationDatabase.MachineParser.PrintIni()
        return
    else:
        print_mdb_entry()
        return
        
def list_conf():
    list_configuration()
def list_configurations():
    
    
    #display("list_configurations")
    
    simlib.RequireMachine()
    
    configs = simlib.GetConfigurations()
    configs.sort()

    if simenv.OptionsManager.GetOption('name-only'):
        for config in configs:
            print config
        return
    
    if len(configs) == 0:
        display("There are no configurations")
        return

    if simenv.VERBOSE:
        display("Configurations:")
    
    for config in configs:
        display("   %-40s" % config, noNewLine=True)
        exe = simlib.BuildPath(simenv.CACTUS_PATH, 'exe', 'cactus_%s' % config)
        if not(os.path.exists(exe)):
            display("   [incomplete]")
        else:
            statinfo = os.stat(exe)
            tinfo = time.localtime(statinfo.st_mtime)
            
            date = "%s-%02d-%02d %2d:%02d:%02d" % (tinfo.tm_year,tinfo.tm_mon,tinfo.tm_mday, tinfo.tm_hour,tinfo.tm_min,tinfo.tm_sec)
            display("   [built %s]" % date)
############################### MAIN ###############################

def main():
    global known_commands
    
    ############################################
    
    ## HEADER ##
    
    if len(simenv.OptionsManager.args) == 0 and simenv.COMMAND == None:
        simenv.OptionsManager.PrintHelp()
        sys.exit(0)

    cactusDir = simenv.CACTUS_PATH
    
    if simenv.VERBOSE:
        info("Cactus Directory: %s" % cactusDir)
    
    #if os.getcwd() != cactusDir:
        #info("Current Working directory does not match Cactus sourcetree, changing to %s" % cactusDir)
        # make sure we're in the Cactus source directory, otherwise all of this will blow way up.
        #os.chdir(cactusDir)
    
    if simenv.COMMAND != None:
        command = simenv.COMMAND
    else:
        command = simenv.OptionsManager.args[0]
    
    if command not in known_commands:
        fatal("unknown sim-info command %s" % command)
    
    # - not a valid function character
    command = command.replace("-", "_")
    
    # execute the reqisite function.
    globals()[command]()
