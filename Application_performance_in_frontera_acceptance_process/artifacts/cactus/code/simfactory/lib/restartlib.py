# restartlib -- all the helper functions to assist the SimRestart class.
# this seperation allows SimRemote to remain a fairly clean abstraction.

import errno
import time
import shutil
import sys, os, re, math
import random
import simrestart
import simenv,simlib,simsubs
import traceback

from libutil import *

class SimulationLog:
    def __init__(self, restart):
        self.Restart = restart
        self.LogPath = simlib.BuildPath(self.Restart.SimulationDir, "log.txt")
        self.LogPointer = None
        
        self.Open()
        
    def Open(self):
        
        if not(os.path.exists(self.Restart.SimulationDir)):
            try:
                os.makedirs(self.Restart.SimulationDir)
            except OSError, e:
                print "Could not create log directory %s, %s" % (self.Restart.SimulationDir, e)
        try:
            self.LogPointer = open(self.LogPath, "a+")
        except:
            fatal("Could not open simulation log file \"%s\" for writing" % self.LogPath)
        
    def Write(self, statement):
        
        stack = traceback.extract_stack()
        
        ptr = -3
        frame = stack[ptr]
        function = frame[3]

        log_date = time.strftime("%Y-%m-%d %H:%M:%S")
          
        if self.LogPointer == None:
            self.Open()
        
        statement = str(statement)
        
        if statement.count("\n") > 0:
            lines = statement.split("\n")
        else:
            lines = [statement]
        
        for line in lines:
            self.LogPointer.write("[LOG:%s] %s::%s\n" % (log_date, function, line))
            self.LogPointer.flush()
    
    def Close(self):
        self.LogPointer.close()
        self.LogPointer = None
        
class WallTime:
    def __init__(self, walltime=None):
        self.Walltime = walltime
        self.parseWalltime()
    
    def parseWalltime(self):
        
        if self.Walltime == None:
            self.Walltime = '876:00:00' # one year
            warning("Walltime not specified, using %s instead" % self.Walltime)
        
        if self.Walltime.count(":") == 0:
            fatal("Wall time has invalid format, expecting HH[:MM[:SS]]")
        
        parts = self.Walltime.split(":")
        
        self.walltime_hh = "%02d" % int(parts[0])
        self.walltime_mm = "00"
        self.walltime_ss = "00"
        
        if len(parts) == 2:
            self.walltime_mm = "%02d" % int(parts[1])

        if len(parts) == 3:
            self.walltime_mm = "%02d" % int(parts[1])
        
        self.Walltime = "%d:%02d:%02d" % (int(self.walltime_hh), int(self.walltime_mm), int(self.walltime_ss))
        
        self.walltime_seconds = (int(self.walltime_hh) * 3600) + (int(self.walltime_mm) * 60) + int(self.walltime_ss)
        self.walltime_minutes = self.walltime_seconds / 60
        self.walltime_hours = self.walltime_seconds / 3600.0
        
def GetNumberOfRestarts(maxwalltime, walltime):
    return int(math.ceil(walltime.walltime_seconds/maxwalltime.walltime_seconds))

def GetConfiguration():

    if not(simenv.OptionsManager.HasOption("configuration")):
        config = simlib.GetDefaultConfiguration()
        info("Configuration name not specified -- using default configuration \"%s\"" % config)
        return config
    
    return simenv.OptionsManager.GetOption("configuration")
    

def GetExecHost(restart):

    DefineDatabase = simsubs.DefineDatabase()
    
    (machine, machineEntry, sourceBaseDir) = simlib.GetLocalEnvironment()
    
    job_id = restart.GetJobId()
    job_status = 'U'

    if job_id != -1:
        job_status = GetJobStatus(job_id)
    
    if job_status != 'R':
        warning("Job is not running, cannot retreive exechost")
        return None
    
    simlib.VerifyKeys(machineEntry, ['exechost', 'exechostpattern'])
    
    DefineDatabase.Set('JOB_ID', job_id)
    DefineDatabase.Set('USER', machineEntry.user)
    DefineDatabase.Set('SIMULATION_NAME', restart.SimulationName)
    
    exechost = DefineDatabase.SubAll(machineEntry.GetKey('exechost'))
    exechostpattern = DefineDatabase.SubAll(machineEntry.GetKey('exechostpattern'))
    
    output,rc = simlib.ExecuteCommand(exechost, True)
    
    rx = re.compile(exechostpattern, re.MULTILINE)
    
    matches = rx.search(output)
    
    if rc or matches == None:
        warning("Unable to retrieve exechost using pattern %s" % exechostpattern)
        return None
    
    return matches.group(1)
    
# TODO: this is a high-level function; it should live in simrestart,
# not in restart lib. (In general, restartlib should probably not call
# simrestart.)
def CleanupRestarts():

    #dprint("CleanupRestarts...");
    
    if simenv.LocalMachine == None:
        return
    
    for sim in simlib.GetSimulations():

        # lets peek in to see if the simulation has an active restart
        if not HasActiveRestart(sim):
            continue

        restart = simrestart.SimRestart()
        ret = restart.load(sim)
        
        
        # get rate stamp.
        current_time = time.time()
        
        # TODO: this logic need to be moved into simrestart.finish
        # make sure the simulation is at least 1 minute old.
        sim_create_time = restart.GetMarkTime("simulation")
        
        if sim_create_time != None and (current_time - sim_create_time) < 60:
            restart.SimulationLog.Write("Skipping simulation %s, simulation create time %s is not 60 seconds from current mark %s" % (sim, sim_create_time, current_time))
            restart.done()
            continue
            
        mark_time = restart.GetMarkTime()

        # only attempt a cleanup every wo minutes apart.
        if mark_time != None and (current_time - mark_time) < 30:
            restart.SimulationLog.Write("Skipping simulation %s, mark time %s is not 30 seconds from current mark %s" % (sim, mark_time, current_time))
            restart.done()
            continue
        
        # less than zero, not greater than zero.
        if ret < 0:
            continue
        
        active_id = GetActiveRestartId(restart)
         
        if active_id != -1:
            ret = restart.load(sim, active_id)
            
            if ret > 0:
                restart.finish()
                restart.done()

def HasActiveRestart(sim):
    DefineDatabase = simsubs.DefineDatabase()
    DefineDatabase.Set("USER", simenv.LocalMachineEntry.user)

    simdir = DefineDatabase.SubAll(simlib.BuildPath(simenv.LocalMachineEntry.basedir, sim))
    
    for elem in os.listdir(simdir):
        if elem.endswith("-active"):
            return True
    
    return False

def GetRestartByJobId(sim, job_id):

    restart = simrestart.SimRestart()
    restart.load(sim)
    
    rids = GetRestartIds(restart)
    
    for rid in rids:
        restart = simrestart.SimRestart()
        restart.load(sim, rid)

        if job_id == restart.GetJobId():
            return restart
        
        restart.done()

def StopAllActiveRestarts(sim):
    
    restart = simrestart.SimRestart()
    restart.load(sim)
    
    rids = GetRestartIds(restart)

    for rid in rids:
        restart = simrestart.SimRestart()
        restart.load(sim, rid)

        job_id = restart.GetJobId()
        
        if job_id == -1:
            restart.done()
            continue
        
        status = GetJobStatus(job_id)
        
        if status != 'U':
            restart.stop()
        
        restart.done()
    
def GetRestartByRestartId(sim, restart_id):
    restart = simrestart.SimRestart()
    ret = restart.load(sim, restart_id)
    
    if ret < 0:
        fatal("Could not load simulation %s with restart id %s" % (sim, restart_id))
        
    return restart

def CheckActive(simulationName):
    
    # TODO: why do we load the restart?
    restart = simrestart.SimRestart()
    restart.load(simulationName)
    active_id = GetActiveRestartId(restart)
    restart.done()
    return active_id >= 0

    # This was wrong; a restart in the U state is still active
    #if active_id == -1:
    #    return False
    #
    #restart = simrestart.SimRestart()
    #restart.load(simulationName, active_id)
    #
    #job_id = restart.GetJobId()
    #status = False
    #
    #if job_id > 0:
    #    job_status = GetJobStatus(job_id)
    #    if job_status != 'U':
    #        status = True
    #
    #restart.done()
    #
    #return status

def GetActiveJobId(simulationName):
    
    restart = simrestart.SimRestart()
    restart.load(simulationName)

    active_id = GetActiveRestartId(restart)

    if active_id == -1:
        return None
    
    restart.load(simulationName, active_id)
    job_id = restart.GetJobId()
    
    restart.done()
    
    return job_id

def GetMaxJobId(simulationName):
    # TODO: rename this function, because it return the job ID of the
    # penultimate (not the ultimate) restart
    
    restart = simrestart.SimRestart()
    restart.load(simulationName)

    # Subtract one, so that the current restart ID is ignored
    max_id = GetMaxRestartID(simulationName)-1

    if max_id < 0:
        return None
    
    restart.load(simulationName, max_id)
    job_id = restart.GetJobId()
    
    restart.done()
    
    return job_id

def GetJobStatus(job_id):
    # TODO: there needs to be documentation stating what the possible
    # return values are
    
    DefineDatabase = simsubs.DefineDatabase()
    
    (machine, machineEntry, sourceBaseDir) = simlib.GetLocalEnvironment()
    simlib.VerifyKeys(machineEntry, ['getstatus', 'queuedpattern', 'runningpattern', 'statuspattern', 'user'])
    
    status_command  = machineEntry.GetKey('getstatus')
    status_pattern  = machineEntry.GetKey('statuspattern')
    queued_pattern  = machineEntry.GetKey('queuedpattern')
    running_pattern = machineEntry.GetKey('runningpattern')
    
    holding_pattern = machineEntry.GetKey('holdingpattern')
    
    user = machineEntry.GetKey('user')
    
    DefineDatabase.Set('USER', user)
    DefineDatabase.Set('JOB_ID', job_id)
    
    status_command = DefineDatabase.SubAll(status_command)
    status_pattern = DefineDatabase.SubAll(status_pattern)
    queued_pattern = DefineDatabase.SubAll(queued_pattern)
    running_pattern = DefineDatabase.SubAll(running_pattern)
    holding_pattern = DefineDatabase.SubAll(holding_pattern)
    
    #capture output.
    output, ret = simlib.ExecuteCommand(status_command, output=True)
    
    # U == unknown?
    status = 'U'
    
    if ret != None:
        lines = output.split("\n")

        #first, lets see if the job is still in the queue, regardless of whether or not we can determine the job status
        
        # TODO: counting lines doesn't help detecting errors
        InQueue = lines.count(job_id) > 0
    
        matched = []
        
        for line in lines:
            matches = re.search(status_pattern, line)
            if matches != None:
                # queued_pattern
                if re.search(queued_pattern, line):
                    status = 'Q'
                    matched.append(queued_pattern)
                    
                # running_pattern
                if re.search(running_pattern, line):
                    status = 'R'
                    matched.append(running_pattern)
                
                # holding_pattern
                if holding_pattern != None:
                    if re.search(holding_pattern, line):
                        status = 'H'
                        matched.append(holding_pattern)

        # TODO: matches is set in the last loop iteration only; this if
        # statement is bogus
        if matches > 1:
            # TODO: output a better error message; the list of patterns
            # alone is not useful because patterns are difficult to read;
            # better would be also what state the patterns are for
            fatal("multiple status patterns matched: %s" % matched)

        if InQueue and status == 'U':
            status = 'E'
        
    if status == 'U':
        warning("job status is U")
        dprint("queue status return code is: %d" % ret)
        dprint("queue status output is:")
        for line in lines:
            dprint("  lines=[%s]" % line)
        dprint("patterns are:")
        dprint("  status_pattern=[%s]" % status_pattern)
        dprint("  queued_pattern=[%s]" % queued_pattern)
        dprint("  running_pattern=[%s]" % running_pattern)
        dprint("  holding_pattern=[%s]" % holding_pattern)

    return status
    
def GetRestartIds(restart):
    outputpat = re.compile("output-([0-9]+)")
    
    ids = []
    try:
        for ff in os.listdir(restart.SimulationDir):
            matches = outputpat.match(ff)
            if matches:
                iid = matches.group(1)
                ids.append(int(iid))
    except OSError, e:
        if e.errno != errno.ENOENT:
            raise
    
    ids.sort()
    return ids

def GetActiveRestartId(restart, warnOnly=False):
    activepat = re.compile("output-([0-9]+)-active")
    
    ids = []
    try:
        for ff in os.listdir(restart.SimulationDir):
            matches = activepat.match(ff)
            if matches:
                iid = matches.group(1)
                ids.append(int(iid))
    except OSError, e:
        if e.errno != errno.ENOENT:
            raise
        return -2
    
    if len(ids) == 0:
        return -1

    if len(ids) > 1:
        if warnOnly:
            warning("more than one active restart id for simulation %s found, %s" % (restart.SimulationName, ids))
            return -2
        else:
            fatal("more than one active restart id for simulation %s found, %s" % (restart.SimulationName, ids))
    
    return ids[0]
    
def GetCheckpointFiles(workdir):
    files = []
    raw_output, rc = simlib.ExecuteCommand("find %s -name *chkpt.it_*" % workdir, output=True)
    
    if rc:
        warning("Searching for checkpoint files failed")
    else:
        raw_files = raw_output.split("\n")
        for ff in raw_files:
            ff = ff.strip()
            if len(ff) > 0:
                files.append(ff)

    return files
        
def GetExecutable():
    configuration = GetConfiguration()
    
    (machine, machineEntry, sourceBaseDir) = simlib.GetLocalEnvironment()
    
    configPath = simlib.BuildPath(simenv.CONFIGS_PATH, configuration)
    
    if not(simlib.FileExists(configPath)):
        fatal("configuration '%s', which has path '%s' does not exist or is not readable" % (configuration, configPath))
    
    exe = simlib.BuildPath(simenv.CACTUS_PATH, 'exe', 'cactus_%s' % configuration)
    
    if not(simlib.FileExists(exe)):
        fatal("Executable %s for configuration %s does not exist or is not readable" % (exe, configuration))
    
    submitScript = simlib.BuildPath(configPath, "SubmitScript")
    
    if not(simlib.FileExists(submitScript)):
        warning("empty submit script for configuration %s" % configuration)
        submitScript = None
    
    runScript = simlib.BuildPath(configPath, "RunScript")
    
    if not(simlib.FileExists(runScript)):
        fatal("Error: empty/missing run script for configuration %s" % configuration)
        
    return (exe, submitScript, runScript)

def SubmitInteractiveRequest(command):
    
    #give up control to the executing terminal
    simenv.system(command)
    
    #nodes = []
    
    #lines = output.split("\n")
    
    #alias_pattern = "^([A-Za-z0-9-]+)$"
    
    #for i in range(len(lines)):
    #   line = lines[i]
        
    #   if line.startswith("PBS has allocated"):
    #       for j in range(i+1, len(lines)):
    #           subline = lines[j]
    #           if subline.startswith("A total of"):
    #               return nodes
    #           
    #           matches = re.search(alias_pattern, subline)
    #           if matches != None:
    #               node = subline.strip()
    #               nodes.append(node)
    #
    #return nodes
            
def CreateSimulationId(simulationName):
    (machine, machineEntry, sourceBaseDir) = simlib.GetLocalEnvironment()
    
    hostname = machineEntry.hostname
    
    # make sure fd is closed.
    fd = simenv.popen('whoami')
    user = fd.read().strip()
    fd.close()
    
    tt = time.localtime()
    
    timestamp = "%4d.%02d.%02d-%02d.%02d.%02d" % (tt.tm_year, tt.tm_mon, tt.tm_mday, tt.tm_hour, tt.tm_min, tt.tm_sec)
    
    pid = os.getpid()
    
    simulation_id = "simulation-%s-%s-%s-%s-%s-%s" % (simulationName, machine, hostname, user, timestamp, pid)
    
    return simulation_id

def CreatePbsSimulationName(SimRestart):
    
    
    simulationName = "%s-%s" % (SimRestart.SimulationName, SimRestart.LongRestartID)
    
    pid = os.getpid()
    
    if simenv.OptionsManager.HasOption('hide') and simenv.OptionsManager.GetOption('hide'):
        shortString = "sim-%06d" % pid
    elif simenv.OptionsManager.HasOption('hide-boring') and simenv.OptionsManager.GetOption('hide-boring'):
        words = ['headon', 'D3.0', 'a0.6', 'mu0.25', 'PN1.5', 'FMR', '1+log', 'nowaves', 'findAH', 'coarse', 'singleBH', 'PUGH', 'movie']
        
        random.seed()
        randomWord = words[random.randint(0, len(words)-1)]
        
        shortString = "sim-%s-%s" % (randomWord, pid)
    elif simenv.OptionsManager.HasOption('hide-dangerous') and simenv.OptionsManager.GetOption('hide-dangerous'):
        words = ['paramesh', 'D25.0', 'a0.999', 'mu0.01', 'PN4.0', 'CCM', 'spec35', 'maximal', 'string', 'FE', 'tail', 'DSS', 'PRL', 'naked']
         
        random.seed()
        randomWord = words[random.randint(0, len(words)-1)]
        
        shortString = "sim-%s-%s" % (randomWord, pid)
    else:
        shortString = simulationName
    
    shortString = re.sub("[^\x20-\x7E]", "", shortString)
    shortString = re.sub("[\s]", "_", shortString)
    shortString = re.sub("^(?![A-Za-z])", "J", shortString)
    
    # limit to 15 characters.
    shortString = shortString[:15]

    return shortString      
        
def CopyFileWithCaching(srcfile, destdir, cachedir):
    """Copy 'srcfile' into 'destdir' using and maintaining 'cacheDir' as a hard link cache."""

    # os.link == create hard link.
    # os.makedirs() == recursively make directories.
    
    if not(os.path.exists(cachedir)):
        try:
            os.makedirs(cachedir)
        except:
            fatal("could not create cache directory: %s" % cachedir)
        
    filename = simlib.BaseName(srcfile)
    
    cachefile = simlib.BuildPath(cachedir, filename)
    dstfile = simlib.BuildPath(destdir, filename)
    
    if not(os.path.exists(cachefile)):
        if os.path.exists(dstfile):
            try:
                os.remove(dstfile)
            except (Exception, EnvironmentError), e:
                fatal("Could not remove existing destination file %s, %s" % (dstfile, e))
        
        try:
            shutil.copyfile(srcfile, dstfile)
            mode = os.stat(srcfile).st_mode
            os.chmod(dstfile, mode)
        except:
            fatal("Could not copy %s to %s" % (srcfile, dstfile))
        
        try:
            os.link(dstfile, cachefile)
            # If the system creates a symlink instead, copy the file. 
            if os.path.islink(cachefile):
               try:
                  os.remove(cachefile)
               except (Exception, EnvironmentError), e:
                  fatal("Could not remove existing symlink cache file %s, %s" % (cachefile, e))
               raise IOError 
        except:
            # If there is an error, copy the file. It could be that
            # the maximum number of hard links per file has been
            # reached.
            shutil.copyfile(dstfile, cachefile)
            mode = os.stat(dstfile).st_mode
            os.chmod(cachefile, mode)
        return
    
    # cachefile exists
    try:
        os.link(cachefile, dstfile)
        # If the system creates a symlink instead, copy the file. 
        if os.path.islink(dstfile):
           try:
              os.remove(dstfile)
           except (Exception, EnvironmentError), e:
              fatal("Could not remove existing destination file %s, %s" % (dstfile, e))
           raise IOError 
    except:
        # If there is an error, copy the file. It could be that the
        # maximum number of hard links per file has been reached.
        shutil.copyfile(cachefile, dstfile)
        mode = os.stat(cachefile).st_mode
        os.chmod(dstfile, mode)
    
    eq = True
    
    srcstat = os.stat(srcfile)
    dststat = os.stat(dstfile)
    
    eq = eq and srcstat.st_mtime <= dststat.st_mtime
    eq = eq and srcstat.st_size == dststat.st_size
    
    if not eq:
        os.remove(dstfile)
        shutil.copyfile(srcfile, dstfile)
        mode = os.stat(srcfile).st_mode
        os.chmod(dstfile, mode)

        if os.path.exists(cachefile):
            try:
                os.unlink(cachefile)
            except:
                warning("could not remove existing cached executable %s" % cachefile)
                
        try:
            os.link(dstfile, cachefile)
            # If the system creates a symlink instead, copy the file. 
            if os.path.islink(cachefile):
               try:
                   os.remove(cachefile)
               except (Exception, EnvironmentError), e:
                   fatal("Could not remove existing symlink cache file %s, %s" % (cachefile, e))
               raise IOError 
        except:
            # If there is an error, copy the file. It could be that
            # the maximum number of hard links per file has been
            # reached.
            shutil.copyfile(dstfile, cachefile)
            mode = os.stat(dstfile).st_mode
            os.chmod(cachefile, mode)

def CreateInternalDirs(internaldir):
    roots = ['exe', 'cfg', 'run', 'par', 'data']
    mdirs = []
    
    for root in roots:
        fullpath = simlib.BuildPath(internaldir, root)
        try:
            os.makedirs(fullpath)
        except OSError, e:
            fatal("Error: could not make %s directory \"%s\", %s" % (root, fullpath, e))
        
        mdirs.append(fullpath)
    
    return mdirs

def GetMaxRestartID(simulationName):
    
    restart = simrestart.SimRestart()
    restart.load(simulationName)
    
    rids = GetRestartIds(restart)
    
    if len(rids) == 0:
        max_restart_id = -1
    else:
        max_restart_id = rids[len(rids)-1]
    
    restart.done()
    
    return max_restart_id
        
def CreateRestartSkeleton(simulationName):
    
    (machine, machineEntry, sourceBaseDir) = simlib.GetLocalEnvironment()
    
    basedir = simlib.GetBaseDir(machineEntry)
    
    try:
        simlib.MakeDirs(basedir)
    except OSError, e:
        fatal("could not access simulation base directory %s for reading and writing: %s" %
              (basedir, e))

    simulationdir = simlib.BuildPath(basedir, simulationName)
    internaldir = simlib.BuildPath(simulationdir, simenv.INTERNALDIRECTORY)
    
    if simlib.FileExists(internaldir):
        fatal("cannot create job skeleton directory: Simulation \"%s\" already exists and has been initialized" % simulationdir)
    
    try:
        simlib.MakeDirs(simulationdir)
    except OSError, e:
        fatal("could not create simulation skeleton directory \"%s\", %s" % (simulationdir, e))

    try:
        simlib.MakeDirs(internaldir)
    except OSError, e:
        fatal("could not create simulation skeleton directory \"%s\", %s" % (internaldir, e))

    cachedir = simlib.BuildPath(basedir, 'CACHE')
    
    try:
        simlib.MakeDirs(cachedir)
    except OSError:
        fatal("Could not create simulation skeleton directory \"%s\", %s" % (cachedir, e))
    
    return (basedir, simulationdir, internaldir, cachedir)
