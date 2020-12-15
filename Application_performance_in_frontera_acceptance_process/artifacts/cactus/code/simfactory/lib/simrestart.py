"""This module contains the SimRestart class, the encapsulation of a
'Restart' which is a submitted simulation created and managed by
simfactory"""

import os, sys, re, stat
import filecmp
import restartlib
import shutil
import time

from libutil import *

import simenv,simlib,simsubs,simarchive,simproperties

class SimRestart:

    """A SimRestart object is an encapsulation of a 'Restart' which is a
    submitted simulation created and managed by simfactory"""

    def __init__(self):
        """Initialise the SimRestart object"""
        self.JobID = -1
        self.RestartID = -1
        self.LongRestartID = -1
        
        self.SimulationName = None
        self.Properties = None
        self.SimulationLog = None
    
    def ClearMark(self, label="timestamp"):
        
        expression = "%s-([0-9\\.]+$)" % label

        rx = re.compile(expression)
            
        # make this as quick as possible
        for fptr in os.listdir(self.InternalDir):
            matches = rx.match(fptr)
            
            if matches != None:
                path = simlib.BuildPath(self.InternalDir, fptr)
                # remove existing mark
                try:
                    os.unlink(path)
                except OSError, e:
                    fatal("unable to delete simulation mark file %s, %s" % (path, e))
                    
    def mark(self, label="timestamp"):
        
        self.ClearMark(label)

        timestamp = time.time()
        
        np = simlib.BuildPath(self.InternalDir, "%s-%s" % (label, timestamp))
        
        try:
            fptr = open(np, 'w+')
            fptr.close()
        except OSError, e:
            fatal("Unable to write %s file %s, %s" % (label, np, e))
        
        
    def GetMarkTime(self, label="timestamp"):
        expression = "%s-([0-9\\.]+$)" % label
        
        rx = re.compile(expression)

        # make this as quick as possible
        for fptr in os.listdir(self.InternalDir):
            matches = rx.match(fptr)
            
            if matches != None:
                rr = matches.group(1)
                
                try: 
                    return float(rr)
                except:
                    fatal("could not coerce %s %s into a float" % (label, rr))

        return None
        
    def load(self, simulationName, restartid=None):
        """Load this object with information from the restart 'restartid' in simulation 'simulationName' stored on disk"""
        # need basedir, simulationdir, and internaldir
        machineEntry = simenv.LocalMachineEntry
        
        self.SimulationName = simulationName
        self.BaseDir = simlib.GetBaseDir(machineEntry)
        self.SimulationDir = simlib.BuildPath(self.BaseDir, self.SimulationName)
        
        if not(os.path.exists(self.SimulationDir)):
            warning("simulation \"%s\" does not exist or is not readable" % self.SimulationName)
            return -1
        
        self.InternalDir = simlib.BuildPath(self.SimulationDir, simenv.INTERNALDIRECTORY)
        self.CacheDir = simlib.BuildPath(self.BaseDir, 'CACHE')
        
        self.SimulationLog = restartlib.SimulationLog(self)

        sim1_marker = simlib.BuildPath(self.SimulationDir, 'SIMULATION_ID')
        
        if os.path.exists(sim1_marker):
            return -2
        
        if restartid == None:
            propertyFile = simlib.BuildPath(self.InternalDir, "properties.ini")
            
            if not(os.path.exists(propertyFile)):
                warning("properties.ini does not exist for simulation %s" % self.SimulationName)
                return -1
                    
            self.Properties = simproperties.SimProperties()
            self.Properties.init(propertyFile)
            
            return 1
            
        # load the specific restart.
        self.RestartID = int(restartid)
        self.LongRestartID = "%04d" % int(restartid)

        self.RestartDir = simlib.BuildPath(self.SimulationDir, "output-%s" % self.LongRestartID)
        self.InternalDir = simlib.BuildPath(self.RestartDir, simenv.INTERNALDIRECTORY)
        
        self.SimulationLog.Write("For simulation %s, loaded restart id %s, long restart id %s" % (self.SimulationName, self.RestartID, self.LongRestartID))

        self.RestartIDs = restartlib.GetRestartIds(self)
        self.MaxRestartID = self.RestartIDs[-1]
            
        propertyFile = simlib.BuildPath(self.InternalDir, "properties.ini")
        
        if not(os.path.exists(propertyFile)):
            warning("properties.ini does not exist for simulation %s" % self.SimulationName)
            return -1
            
        self.Properties = simproperties.SimProperties()
        self.Properties.init(propertyFile)
        
        self.Initialized = True
        
        return 1
            
    # TODO: this function name is wrong, it prepares recovery
    def PrepareCheckpointing(self, max_restart_id):
    
        # make sure it's an int.
        max_restart_id = int(max_restart_id)
        
        self.SimulationLog.Write("PrepareCheckpointing: max_restart_id: %s" % max_restart_id)
            
        # lets see if from-restart-id is sest
        if self.Properties.HasProperty('from-restart-id'):
            ii = self.Properties.GetProperty('from-restart-id')
            
            try:
                restore_restart_id = int(ii)
            except ValueError:
                fatal("Could not coerce provided from-restart-id %s into an integer" % ii)
        else:
            if max_restart_id >= 0:
                restore_restart_id = max_restart_id
            else:
                info("max_restart_id is < 0, returning False")
                return False

        if restore_restart_id == self.RestartID:
            self.SimulationLog.Write("max_restart_id == self.RestartID, returning false")
            return False

        self.SimulationLog.Write("Restoring from restart id %s, 04d is: %04d" % (restore_restart_id, restore_restart_id))
        
        d_restore_restart_id = "%04d" % restore_restart_id
        
        restore_dir =  simlib.BuildPath(self.SimulationDir, "output-%s" % d_restore_restart_id)
        
        if not(os.path.exists(restore_dir)):
            try:
                os.makedirs(restore_dir)
            except OSError, e:
                fatal("could not make restore_dir %s, %s" % (restore_dir, e))
    
        previousIni = simlib.BuildPath(restore_dir, simenv.INTERNALDIRECTORY, "properties.ini")
        
        work_dir = None
        
        if os.path.exists(previousIni):
            ep = simproperties.SimProperties()
            ep.init(previousIni)
            
            parfile = ep.GetProperty("parfile")
            pf = simlib.FileBaseName(parfile)
            work_dir = simlib.BuildPath(restore_dir, pf)
        
        current_parfile = self.Properties.parfile
        cpf = simlib.FileBaseName(current_parfile)
        cwork_dir = simlib.BuildPath(self.RestartDir, cpf)
        
        self.SimulationLog.Write("Current working directory for simulation is: %s" % cwork_dir)
        
        if not(os.path.exists(cwork_dir)):
            # TODO: state in the message what this directory is needed
            # for
            self.SimulationLog.Write("%s does not exist, creating" % cwork_dir)
            os.makedirs(cwork_dir)
        
        chkpoint_files = restartlib.GetCheckpointFiles(restore_dir)
        if len(chkpoint_files) == 0:
            return False
        # TODO: find checkpoint files from the last iteration, and
        # link only those
            
        # TODO: there are too many log messages about checkpoint files
        for ff in chkpoint_files:
            self.SimulationLog.Write("checkpoint file: %s" % ff)

        for ff in chkpoint_files:
            # TODO: skip this check
            if not(os.path.exists(ff)):
                continue

            dfile = ff.replace(restore_dir, self.RestartDir).strip()
            
            self.SimulationLog.Write("linking %s to %s" % (ff, dfile))
                        
            try:
                dfile = ff.replace(restore_dir, self.RestartDir).strip()
                
                if simenv.VERBOSE:
                    self.SimulationLog.Write("linking %s to %s" % (ff, dfile))
                # TODO: tying linking the file first; this would be
                # faster
                try:
                    if not os.access(os.path.dirname(dfile), os.F_OK):
                        os.makedirs(os.path.dirname(dfile))
                    os.link(ff, dfile)
                except:
                    shutil.copyfile(ff, dfile)
            except OSError, e:
                # TODO: message should mention copying as well
                warning("Could not link checkpoint file %s to %s, %s" % (ff, dfile, e))
                # TODO: this disables recovery, not checkpointing
                # TODO: this doesn't delete any of the checkpoint
                # files that it linked or copied... it should abort
                # instead
                # TODO: silently disabling recovery is a bad idea; it
                # should abort instead
                warning("Disabling checkpointing")
                error("It does not make sense to continue without 'checkpointing'")
                return False

        return True
    
    
    def initRestart(self, simulationName):
        
        ret = self.load(simulationName)
        
        if ret < 0:
            if not simenv.OptionsManager.HasOption('testsuite'):
                parfile = simlib.GetParFile()

                if parfile == None:
                    fatal("could not create simulation %s, no parameter file specified" % simulationName)
            
                display("Parameter file: %s" % parfile)
            else:
                parfile = ""

            self.create(simulationName, parfile)
            self.load(simulationName)

        # Try to clean up
        active_id = restartlib.GetActiveRestartId(self)
        if active_id >= 0:
            info("Restart %d is active" % active_id)
            job_id = restartlib.GetActiveJobId(simulationName)
            if job_id<0:
                fatal("job id is negative")
            job_status = restartlib.GetJobStatus(job_id)
            if job_status == 'E':
                fatal("Could not determine restart status")
            if job_status == 'U':
                info("Cleaning up restart %d before submitting" % active_id)
                restart = SimRestart()
                restart.load(simulationName, active_id)
                restart.finish()
                restart.done()

        (rids, max_restart_id, my_restart_id) = self.GetRestartId()
        
        self.RestartIDs = rids
        self.MaxRestartID = max_restart_id

        info("Found the following restart_ids: %s" % rids)

        if max_restart_id != -1:
            info("Maximum restart id determined to be: %04d" % max_restart_id)

        display("Assigned restart id: %s " % my_restart_id)

        self.RestartID = my_restart_id
        self.LongRestartID = "%04d" % self.RestartID
        
        self.RestartDir = simlib.BuildPath(self.SimulationDir, "output-%s" % self.LongRestartID)
        self.InternalDir = simlib.BuildPath(self.RestartDir, simenv.INTERNALDIRECTORY)
        
        for dd in [self.RestartDir, self.InternalDir]:
            try:
                os.makedirs(dd)
            except OSError, e:
                fatal("could not create directory \"%s\", %s" % (dd, e))
            
        pass

        self.Initialized = True
        
    def userSubmit(self, simulationName):
        """Submit simulation 'simulationName'.  Calls submit()."""

        self.initRestart(simulationName)
        
        assert(self.Initialized)

        ss = simenv.OptionsManager.GetOption('submitscript')

        if ss:
            if not os.path.exists(ss):
                fatal("submit script '%s' does not exist" % ss)
            if os.path.isdir(ss):
                fatal("submit script '%s' is a directory" % ss)
            if os.path.getsize(ss) == 0:
                fatal("submit script '%s' is an empty file" % ss)
            submitScript = ss
        else:
            submitScript = self.Properties.submitscript

        if len(submitScript) == 0 or submitScript == "None":
            fatal("no default submission script defined for simulation %s, submission disabled" % self.SimulationName)

        if not(os.path.exists(submitScript)):
            fatal("could not read submit script %s" % submitScript)
            
        if restartlib.CheckActive(simulationName):
            
            # we presubmit
            chainedJobId = restartlib.GetMaxJobId(simulationName)
            self.ActiveRestartID = restartlib.GetActiveRestartId(self)
            # subtract one, so that the restart that we are just
            # creating is ignored
            maxRestartId = restartlib.GetMaxRestartID(simulationName)-1
            warning("Simulation %s already has an active submission. Chaining this submission onto job id %s" % (simulationName, chainedJobId))
            self.submit(submitScript, chainedJobId, maxRestartId)
        else:
            self.submit(submitScript)
    
    def IsActive(self):
        """Return whether or not this restart is active"""

        active_id = restartlib.GetActiveRestartId(self)
        return active_id == self.RestartID
        
    def preSubmit(self, simulationName, chainedJobId=None, chainedRestartId=None):

        self.initRestart(simulationName)

        assert(self.Initialized)

        ss = simenv.OptionsManager.GetOption('submitscript')

        if ss:
            if not os.path.exists(ss):
                fatal("Submit script '%s' does not exist" % ss)
            submitScript = ss
        else:
            submitScript = self.Properties.submitscript

        if len(submitScript) == 0 or submitScript == "None":
            fatal("no default submission script defined for simulation %s, submission disabled" % self.SimulationName)

        if not(os.path.exists(submitScript)):
            fatal("could not read submit script %s" % submitScript)

        self.submit(submitScript, chainedJobId, chainedRestartId)

    def GetPreviousProperties(self):
        
        assert(self.Initialized)
        
        if self.MaxRestartID == -1:
            return None
        
        restore_dir = simlib.BuildPath(self.SimulationDir, "output-%04d" % self.MaxRestartID)
        previousIni = simlib.BuildPath(restore_dir, simenv.INTERNALDIRECTORY, "properties.ini")

        properties = None
        
        if os.path.exists(previousIni):
            properties = simproperties.SimProperties()
            properties.init(previousIni)
        
        return properties
        
    def makeActive(self):
        
        assert self.Initialized
        assert self.RestartID > -1
        
        # TODO: ensure the simulation is not active
        active_id = restartlib.GetActiveRestartId(self)
        if active_id != -1:
            fatal("Internal error: Cannot submit simulation %s because it is already active" % self.SimulationName)
        
        curdir = os.getcwd()
        
        # set our current working directory to the restart dir in case the 
        # submit script writes out any files.
        # TODO: we must not chdir
        os.chdir(self.RestartDir)
        
        try:
            # create a symlink indicating that this job is active.
            # create a relative symlink, not an absolute one
            # TODO: What should we do on architectures that don't
            # support symbolic links?
            os.symlink("output-%s" % self.LongRestartID, "%s-active" % self.RestartDir)
        except:
            fatal("unable to activate simulation %s, with restart id %s" % (self.SimulationName, self.RestartID))
    
        os.chdir(curdir)
        
        self.SimulationLog.Write("Simulation %s with restart-id %s has been made active" % (self.SimulationName, self.RestartID))
        
    def extractJobId(self, output):
        """Parse the queuing system job id from the output of the submit command"""
        assert(self.Initialized)
        
        machineEntry = simenv.LocalMachineEntry
        
        self.SimulationLog.Write("received raw output: %s" % output)

        submitRegex = machineEntry.submitpattern

        self.SimulationLog.Write("using submitRegex: %s" % submitRegex)

        rx = re.compile(submitRegex, re.MULTILINE)

        matches = rx.search(output)
        
        # if we didn't match anything, just use whatever got output.
        if matches == None:
            job_id = "-1"
        else:
            job_id = matches.group(1)
        
        return job_id
    
    def copyTestsuiteData(self):
        """Copy testsuite data into restart directory"""
        if self.Properties.testsuite == "True":
            display("Copying testsuite data")
            testsuitedir = self.RestartDir
            testexe=os.path.join(testsuitedir,"exe")
            os.mkdir(testexe)
            try:
                os.link(self.Properties.executable, os.path.join(testexe,simlib.BaseName(self.Properties.executable)))
            except OSError, e:
                shutil.copyfile(self.Properties.executable, os.path.join(testexe,simlib.BaseName(self.Properties.executable)))
            
            test = self.Properties.select_tests
            if test == "all":
                includes = ["--include='/arrangements/*/*/test/**'"]
            elif test.endswith(".par"):
                testdirname = test[0:-4]
                includes = ["--include='/arrangements/*/*/test/%s'" % test,
                            "--include='/arrangements/*/*/test/%s'" % testdirname,
                            "--include='/arrangements/*/*/test/%s/**'" % testdirname]
            elif test.count("/") > 1:
                fatal("Test specification '%s' should be no, all, <test>.par, <arrangement> or <arrangement>/<thorn>")
            elif test.count("/") == 0:
                # Arrangement
                includes =[ "--include='/arrangements/%s/*/test/**'" % test]
                if not os.path.exists(os.path.join(self.Properties.sourcedir,"arrangements/%s"%test)):
                    fatal("Arrangement %s not found" % test)
            else:
                # Arrangement/Thorn
                includes = ["--include='/arrangements/%s/test/**'" % test]
                if not os.path.exists(os.path.join(self.Properties.sourcedir,"arrangements/%s"%test)):
                    fatal("Thorn %s not found" % test)
            # Would like to have checks here that the argument names
            # something which is actually present in the thorn list

            localMachineEntry = simenv.LocalMachineEntry        
            rsyncInfo = simlib.GetRsyncInfo(localMachineEntry)
            (rsynccmd, rsyncopts) = rsyncInfo
            rsyncversion = simlib.RsyncVersion(rsyncInfo)
            oldRsync = True     # a safe choice
            rsyncoptions = ['--archive', '--copy-links']
            ruleFile = 'filter.rules'
            if rsyncversion[0] < 3 or oldRsync:
                ruleFile = 'filter.prersync3.rules'
            rsyncoptions = (rsyncoptions +
                            ["--filter 'merge %s/etc/%s'" % (simenv.BASE_PATH, ruleFile)])

            # --prune-empty-dirs could be used to avoid listing all the parent directories, but this doesn't seem to work
            rsyncfiles = (["--include='/arrangements/'",
                           "--include='/arrangements/*/'",
                           "--include='/arrangements/*/*/'",
                           "--include='/arrangements/*/*/test'",
                           "--include='/arrangements/*/*/test/test.ccl'",
                           "--include=/Makefile",
                           "--include=/lib",
                           "--include='/lib/**'",
                           "--include='/configs'",
                           "--include='/configs/%s/'" % self.Properties.configuration,
                           "--include='/configs/%s/ThornList'" % self.Properties.configuration,
                           "--include='/configs/%s/config-data/'" % self.Properties.configuration,
                           "--include='/configs/%s/config-data/cctk_Extradefs.h'" % self.Properties.configuration,
                           "--include='/configs/%s/bindings/'" % self.Properties.configuration,
                           "--include='/configs/%s/bindings/Configuration/'" % self.Properties.configuration,
                           "--include='/configs/%s/bindings/Configuration/Capabilities/'" % self.Properties.configuration,
                           "--include='/configs/%s/bindings/Configuration/Capabilities/cctki_MPI.h'" % self.Properties.configuration] +
                          includes +
                          ["--exclude='*'",
                           "%s/" % self.Properties.sourcedir,
                           "%s/" % testsuitedir])

            cmd = "%s %s %s %s" % (rsynccmd, rsyncopts,
                                   ' '.join(rsyncoptions), ' '.join(rsyncfiles))
            if simlib.ExecuteCommand(cmd) != 0:
                fatal("Rsync of test data for simulation failed")

            
    def submit(self, submitScript, chainedJobId=None, chainedRestartId=None):
        """Low-level job submission function."""

        assert(self.Initialized)
        
        # TODO: this never checks whether the simulation is already
        # active! only inactive simulations can have restarts
        # submitted.

        active_id = restartlib.GetActiveRestartId(self)
        if active_id >= 0:
            info("Simulation is active: presubmitting")
        else:
            info("Simulation is inactive: submitting")

        # If chaining is requested, check that the chained job is in the queue.
        if chainedJobId != None:
            chainedJobStatus = restartlib.GetJobStatus(chainedJobId)
            if chainedJobStatus != "R" and chainedJobStatus != "Q" and chainedJobStatus != "H":
                warning("Job chaining requested but job id %s is not in the queue. Its status is %s. Aborting submission." % (chainedJobId, chainedJobStatus) )
                return

        # TODO: Below, use the same job id as here, since the job may
        # finish in the mean time. Also handle the case where we
        # presubmit for a job that has just finished -- in this case,
        # it should become a regular submit.

        UseChaining = False
    
        DefineDatabase = simsubs.DefineDatabase()

        machineEntry = simenv.LocalMachineEntry
        simlib.VerifyKeys(machineEntry, ['num-threads', 'num-smt', 'ppn', 'maxwalltime', 'hostname'])
        
        submitscript_contents = simlib.GetFileContents(submitScript)

        MaxWalltime = restartlib.WallTime(machineEntry.maxwalltime)
        
        if not(simenv.OptionsManager.HasOption('walltime')):
            Walltime = MaxWalltime
        else:
            Walltime = restartlib.WallTime(simenv.OptionsManager.GetOption('walltime'))

        self.SimulationLog.Write("Restart for simulation %s created with restart id %s, long restart id %s" % (self.SimulationName, self.RestartID, self.LongRestartID))
        self.SimulationLog.Write("Prepping for submission")
        
        existingProperties = self.GetPreviousProperties()

        # import walltime if no --walltime is specified.
        if existingProperties != None and not simenv.OptionsManager.HasOption('walltime') and existingProperties.HasProperty('walltime'):
            Walltime = restartlib.WallTime(existingProperties.GetProperty("walltime"))
            self.SimulationLog.Write("Using walltime %s from previous restart %s" % (existingProperties.GetProperty("walltime"), self.MaxRestartID))
        else:
            self.SimulationLog.Write("No previous walltime available to be used, using walltime %s" % Walltime.Walltime)
            
        if simenv.OptionsManager.HasOption("from-restart-id"):
            self.Properties.AddProperty("from-restart-id", simenv.OptionsManager.GetOption("from-restart-id"))
            self.SimulationLog.Write("from-restart-id was specified, value is %s" % simenv.OptionsManager.GetOption("from-restart-id"))
            
        hostname = machineEntry.hostname
        user = simlib.GetUsername()
        memory = simlib.GetMachineOption('memory')
        cpufreq = simlib.GetMachineOption('cpufreq')
        allocation = simlib.GetMachineOption('allocation')
        
        queue = simlib.GetMachineOption('queue')
        
        (nodes, procs_requested, ppn, num_procs, node_procs,
         procs, num_threads, ppn_used, num_smt) = simlib.GetProcs(existingProperties)
        
        parfile = self.Properties.parfile
        pf = simlib.BaseName(parfile)
        newparpath = simlib.BuildPath(self.RestartDir, pf)
        newsspath = simlib.BuildPath(self.InternalDir, "SubmitScript")
        
        pbsSimulationName = restartlib.CreatePbsSimulationName(self)
                
        new_properties = dict()
        new_properties['SOURCEDIR'] = self.Properties.sourcedir
        new_properties['SIMULATION_NAME'] = self.SimulationName
        new_properties['SHORT_SIMULATION_NAME'] = pbsSimulationName
        new_properties['SIMULATION_ID'] = self.Properties.simulationid
        new_properties['RESTART_ID'] = self.RestartID
        new_properties['RUNDIR'] = self.RestartDir
        new_properties['SCRIPTFILE'] = self.Properties.submitscript
        new_properties['EXECUTABLE'] = self.Properties.executable
        new_properties['PARFILE'] = newparpath
        new_properties['HOSTNAME'] = hostname
        new_properties['USER'] = user
        new_properties['NODES'] = nodes
        new_properties['PROCS_REQUESTED'] = procs_requested
        new_properties['PPN'] = ppn
        new_properties['NUM_PROCS'] = num_procs
        new_properties['NODE_PROCS'] = node_procs
        new_properties['PROCS'] = procs
        new_properties['NUM_THREADS'] = num_threads
        new_properties['PPN_USED'] = ppn_used
        new_properties['NUM_SMT'] = num_smt
        new_properties['MEMORY'] = memory
        new_properties['CPUFREQ'] = cpufreq
        new_properties['ALLOCATION'] = allocation
        new_properties['QUEUE'] = queue
        new_properties['EMAIL'] = machineEntry.email
        
        walltt = Walltime
        
        # always restrict our walltime to maxwalltime if requested walltime
        # is too large.
        if MaxWalltime.walltime_seconds < Walltime.walltime_seconds:
            walltt = MaxWalltime

            # okay, our walltime requested was too large
            # find out if we should use automatic job chaining.
            if chainedJobId == None:
                UseChaining = True
                # TODO: i don't understand the job chaining logic. a
                # restart should be presubmitted (instead of
                # submitted) if there is a restart currently running.
                # yet there is no check for this.
        
        new_properties['WALLTIME'] = walltt.Walltime
        new_properties['WALLTIME_HH'] = walltt.walltime_hh
        new_properties['WALLTIME_MM'] = walltt.walltime_mm
        new_properties['WALLTIME_SS'] = walltt.walltime_ss
        new_properties['WALLTIME_SECONDS'] = walltt.walltime_seconds
        new_properties['WALLTIME_MINUTES'] = walltt.walltime_minutes
        new_properties['WALLTIME_HOURS'] = walltt.walltime_hours            
        new_properties['SIMFACTORY'] = simenv.EXECUTABLE
        new_properties['SUBMITSCRIPT'] = newsspath
        new_properties['SCRIPTFILE'] = newsspath
        new_properties['CONFIGURATION'] = self.Properties.configuration

        for key in DefineDatabase.PresetMacros:
            new_properties[key] = DefineDatabase.PresetMacros.get(key)


        if chainedRestartId == None:
            new_properties['FROM_RESTART_COMMAND'] = ""
        else:
            new_properties['FROM_RESTART_COMMAND'] = "--from-restart-id=%s" % chainedRestartId
            
        if chainedJobId == None:
            new_properties['CHAINED_JOB_ID'] = ""
        else:
            new_properties['CHAINED_JOB_ID'] = chainedJobId
    
        for key in new_properties.keys():
            DefineDatabase.Set(key, new_properties[key])
        
        self.SimulationLog.Write("Defined substituion properties for submission")
        self.SimulationLog.Write(new_properties)
        
        # lets prepare our submit script.
        submitscript_contents = DefineDatabase.SubAll(submitscript_contents)
        
        # store the rest of our keys
        rr = {
            'nodes': nodes,
            'procsrequested': procs_requested,
            'ppn': ppn,
            'numprocs': num_procs,
            'nodeprocs': node_procs,
            'procs': procs,
            'numthreads': num_threads,
            'ppnused':  ppn_used,
            'numsmt': num_smt,
            'queue': queue,
            'allocation': allocation,
            'hostname': hostname,
            'user': user,
            'memory': memory,
            'cpufreq': cpufreq,
            'pbsSimulationName': pbsSimulationName,
            'walltime': walltt.Walltime
        }
        
        if chainedJobId == None:
            rr['chainedjobid'] = "-1"
        else:
            rr['chainedjobid'] = str(chainedJobId)
        
        # cast this as dict() to silence a pychecker warning.
        # TODO: really?
        for key in dict(rr).keys():
            self.Properties.AddProperty(key, rr[key])
        
        # write to our new properties directory.
        self.Properties.Filename = simlib.BuildPath(self.InternalDir, "properties.ini")
        
        info("writing to internalDir: %s" % self.InternalDir)
        self.Properties.Save()
        
        self.SimulationLog.Write("self.Properties: %s" % self.Properties.Filename)
        self.SimulationLog.Write(self.Properties.toString())
        
        info("saving substituted submitscript contents to: %s" % newsspath)
        self.SimulationLog.Write("saving substituted submitscript contents to: %s" % newsspath)
        
        simlib.WriteContents(newsspath, submitscript_contents)
        MakeFileExecutable(newsspath)
        
        self.copyTestsuiteData()

        submitCommand = DefineDatabase.SubAll(machineEntry.submit)
        
        if chainedJobId != None:
            if (submitscript_contents.count(chainedJobId) == 0 and
                submitCommand.count(chainedJobId) == 0):
                info ("submitCommand=\"%s\\n" % submitCommand)
                info ("chainedJobId=\"%s\\n" % str(chainedJobId))
                fatal("Machine %s currently does not support job chaining. Please modify the submit command or submission script to support job chaining." % simenv.LocalMachine)
        
        self.SimulationLog.Write("Executing submission command: %s" % submitCommand)
        display("Executing submit command: %s" % submitCommand)
        
        # we only make the simulation active if this isn't a chained job.
        if chainedJobId == None:
            self.makeActive()
        
        output, rc = simlib.ExecuteCommand(submitCommand, output=True)
        
        self.mark("simulation")

        if rc:
            self.SimulationLog.Write("Could not execute submit commad: %d\n%s" % (rc, output))
            job_id = '-1'
        else:
            job_id = self.extractJobId(output)
        
            self.SimulationLog.Write("After searching raw output, it was determined that the job_id is: %s" % job_id)
            if job_id == '-1':
                self.SimulationLog.Write("The regex did NOT match anything. No job_id means no control.")
        
        self.Properties.AddProperty('jobid', job_id)
        self.Properties.Save()
        
        if job_id == "-1":
            warning("submit either failed or could not determine job id, output:")
            warning(output)
            return
        else:
            display("Submit finished, job id is %s" % job_id)
        
        self.JobID = job_id
        
        os.environ['PBS_JOBID'] = self.JobID
        
        # have to make sure chainedJobId is None, otherwise this is submit() is actually handling a preSubmit() call.    
        if UseChaining and chainedJobId == None:
            numChains = restartlib.GetNumberOfRestarts(MaxWalltime, Walltime)
            
            restart = None
            previousRestart = None
            
            # we already submitted the first chain
            
            # TODO: this doesn't set the correct wall time limit for
            # the last presubmitted restart

            while numChains > 0:
                # TODO: this logic looks fishy. since we already
                # submitted a restart, previousRestart should never be
                # None.
                if restart != None:
                    previousRestart = restart
                    
                restart = SimRestart()
                restart.load(self.SimulationName)
                
                if previousRestart == None:
                    restart.preSubmit(self.SimulationName, self.JobID, self.RestartID)
                else:
                    restart.preSubmit(self.SimulationName, previousRestart.JobID, previousRestart.RestartID)
                
                numChains = numChains - 1
                
        self.SimulationLog.Write("Simulation %s, with restart id %s, and job id %s has been submitted" % (self.SimulationName, self.RestartID, self.JobID))
        
    def GetRestartId(self):
    
        rids = restartlib.GetRestartIds(self)
        
        if len(rids) > 0:
            max_restart_id = rids[len(rids)-1]
        else:
            max_restart_id = -1
    
        if simenv.OptionsManager.HasOption("restart-id"):
            ii = simenv.OptionsManager.GetOption("restart-id")

            try:
                my_restart_id = int(ii)
            except ValueError:
                fatal("Could not coerce provided restart-id %s into an integer" % ii)
        else:
            my_restart_id = max_restart_id + 1
        
        if not(simenv.OptionsManager.HasOption("restart-id")):
            if my_restart_id in rids:
                fatal("assigned restart id %s for simulation %s already in use!" % my_restart_id)
            
        info("Assigned restart_id of: %04d" % my_restart_id)
        
        if my_restart_id > 9999:
            fatal("maximum number of restarts reached. Please use sim purge to clear existing restarts")
                    
        if max_restart_id == my_restart_id:
            max_restart_id = max_restart_id - 1
            
        return (rids, max_restart_id, my_restart_id)
    
    def interactive(self):

        DefineDatabase = simsubs.DefineDatabase()
        
        # need basedir, simulationdir, and internaldir
        machineEntry = simenv.LocalMachineEntry
        
        simlib.VerifyKeys(machineEntry, ['interactivecmd'])
        
        # need to prepare replacements.
        
        MaxWalltime = restartlib.WallTime(machineEntry.maxwalltime)
        
        if not(simenv.OptionsManager.HasOption('walltime')):
            Walltime = MaxWalltime
        else:
            Walltime = restartlib.WallTime(simenv.OptionsManager.GetOption('walltime'))
        
        hostname = machineEntry.hostname
        user = simlib.GetUsername()
        memory = simlib.GetMachineOption('memory')
        cpufreq = simlib.GetMachineOption('cpufreq')
        allocation = simlib.GetMachineOption('allocation')
        queue = simlib.GetMachineOption('queue')

        (nodes, procs_requested, ppn, num_procs, node_procs,
         procs, num_threads, ppn_used, num_smt) = simlib.GetProcs(None)
        
        pbsSimulationName = restartlib.CreatePbsSimulationName(self)
        
        new_properties = dict()
        new_properties['HOSTNAME'] = hostname
        new_properties['USER'] = user
        new_properties['SHORT_SIMULATION_NAME'] = pbsSimulationName
        new_properties['NODES'] = nodes
        new_properties['PROCS_REQUESTED'] = procs_requested
        new_properties['PPN'] = ppn
        new_properties['NUM_PROCS'] = num_procs
        new_properties['NODE_PROCS'] = node_procs
        new_properties['PROCS'] = procs
        new_properties['NUM_THREADS'] = num_threads
        new_properties['PPN_USED'] = ppn_used
        new_properties['NUM_SMT'] = num_smt
        new_properties['MEMORY'] = memory
        new_properties['CPUFREQ'] = cpufreq
        new_properties['ALLOCATION'] = allocation
        new_properties['QUEUE'] = queue
        new_properties['EMAIL'] = machineEntry.email
        new_properties['WALLTIME'] = Walltime.Walltime
        new_properties['WALLTIME_HH'] = Walltime.walltime_hh
        new_properties['WALLTIME_MM'] = Walltime.walltime_mm
        new_properties['WALLTIME_SS'] = Walltime.walltime_ss
        new_properties['WALLTIME_SECONDS'] = Walltime.walltime_seconds
        new_properties['WALLTIME_MINUTES'] = Walltime.walltime_minutes
        new_properties['WALLTIME_HOURS'] = Walltime.walltime_hours          
        new_properties['SIMFACTORY'] = simenv.EXECUTABLE
        
        for key in new_properties.keys():
            DefineDatabase.Set(key, new_properties[key])
        
        display("Entering interactive mode")
        interactivecmd = DefineDatabase.SubAll(machineEntry.interactivecmd)
        simlib.ExecuteCommand(interactivecmd)
    
    def userRun(self, simulationName, debug=False):
        """Run this restart.  Calls run()."""

        self.initRestart(simulationName)
        
        assert(self.Initialized)
        
        machineEntry = simenv.LocalMachineEntry

        # do setup that isn't done because this wasn't submitted.
        
        self.SimulationLog.Write("Creating new properties because this is an independant run, not a run following a submit")
        
        hostname = machineEntry.hostname
        user = simlib.GetUsername()
        memory = simlib.GetMachineOption('memory')
        cpufreq = simlib.GetMachineOption('cpufreq')
        
        pbsSimulationName = restartlib.CreatePbsSimulationName(self)
        existingProperties = self.GetPreviousProperties()
        (nodes, procs_requested, ppn, num_procs, node_procs,
         procs, num_threads, ppn_used, num_smt) = simlib.GetProcs(existingProperties)
        
        # store the rest of our keys
        rr = {
            'nodes': nodes,
            'procsrequested': procs_requested,
            'ppn': ppn,
            'numprocs': num_procs,
            'nodeprocs': node_procs,
            'procs': procs,
            'numthreads': num_threads,
            'ppnused':  ppn_used,
            'numsmt': num_smt,
            'hostname': hostname,
            'user': user,
            'memory': memory,
            'cpufreq': cpufreq,
            'pbsSimulationName': pbsSimulationName
        }
        
        # cast this as a dict to silence a pychecker warning.
        for key in dict(rr).keys():
            self.Properties.AddProperty(key, rr[key])
        
        self.Properties.Filename = simlib.BuildPath(self.RestartDir,  simenv.INTERNALDIRECTORY, "properties.ini")
        self.Properties.Save()
        
        self.SimulationLog.Write("Determined the following properties")
        self.SimulationLog.Write(self.Properties.toString())
        
        self.makeActive()

        self.copyTestsuiteData()
        
        self.run(debug)
    
    def submitRun(self, simulationName, restartId):
        
        self.load(simulationName, restartId)

        assert(self.Initialized)
        
        # if we are following a previous restart, attempt to clean that previous restart up
        # and make this new restart active.
        if simenv.OptionsManager.HasOption("from-restart-id"):
            frid = simenv.OptionsManager.GetOption("from-restart-id")
            
            self.SimulationLog.Write("Following restart-id %s, finishing it." % frid)
            restart = SimRestart()
            restart.load(simulationName, frid)
            restart.finish()
        
            self.makeActive()
        
        active_dir = "%s-active" % self.RestartDir
        if not os.path.exists(active_dir):
            fatal("cannot rerun a restart that is inactive.")
        
        self.run()
        

    def run(self, debug=False):
        
        assert(self.Initialized)
        assert(self.IsActive())
        
        DefineDatabase = simsubs.DefineDatabase()
        
        display("Running simulation %s" % self.SimulationName)
        self.SimulationLog.Write("Prepping for execution/run")
        
        recover_id = self.RestartID - 1
        
        if self.Properties.HasProperty("from-restart-id"):
            recover_id = self.Properties.GetProperty("from-restart-id")

        testsuite = self.Properties.testsuite == "True"
        
        # TODO: this is about recovering, not about checkpointing.
        # read the code below with this in mind.
        checkpointing = self.PrepareCheckpointing(recover_id)
            
        # TODO: what does this property do? it doesn't make its way
        # into the parameter file, so it is presumably unused.
        # actually, it is never used anywhere except here.
        if checkpointing:
            self.Properties.AddProperty("checkpointing", "yes")
            self.SimulationLog.Write("Checkpointing from max restart id of %s" % self.MaxRestartID)
        else:
            self.Properties.AddProperty("checkpointing", "no")

        if self.Properties.parfile != "":
            parname = simlib.FileBaseName(self.Properties.parfile)
            my_workdir = simlib.BuildPath(self.RestartDir, parname)

            # TODO: don't look for the directory first, create it right away
            if not(os.path.exists(my_workdir)):
                try:
                    os.mkdir(my_workdir)
                except OSError, e:
                    fatal("could not make working directory path \"%s\", %s" % (my_workdir, e))
        
        # TODO: does this work? is this overwritten by the run script?
        os.environ['CACTUS_STARTTIME'] = str(int(time.time()))
            
        # do parfile substitution
        parfile = self.Properties.parfile
        if parfile != "":
            pf = simlib.BaseName(parfile)
            if pf.endswith(".rpar"):
                pf = re.sub(r'(.*)\.rpar$', r'\1.par', pf)
            newparpath = simlib.BuildPath(self.RestartDir, pf)

        new_properties = dict()
        new_properties['MACHINE'] = self.Properties.machine
        new_properties['SOURCEDIR'] = self.Properties.sourcedir
        new_properties['SIMULATION_NAME'] = self.SimulationName
        new_properties['SHORT_SIMULATION_NAME'] = self.Properties.pbsSimulationName
        new_properties['SIMULATION_ID'] = self.Properties.simulationid
        new_properties['RESTART_ID'] = self.RestartID
        new_properties['SCRIPTFILE'] = self.Properties.submitscript
        new_properties['SUBMITSCRIPT'] = self.Properties.submitscript
        new_properties['CONFIGURATION'] = self.Properties.configuration
        new_properties['EXECUTABLE'] = self.Properties.executable
        if testsuite:
            new_properties['PARFILE'] = "${TESTSUITE_PARFILE}"
            new_properties['RUNDIR'] = "output-0000"
        else:
            new_properties['RUNDIR'] = self.RestartDir
            new_properties['PARFILE'] = newparpath
        new_properties['HOSTNAME'] = self.Properties.hostname
        new_properties['USER'] = self.Properties.user

        new_properties['NODES'] = self.Properties.nodes
        new_properties['PROCS_REQUESTED'] = self.Properties.procsrequested
        new_properties['PPN'] = self.Properties.ppn
        new_properties['NUM_PROCS'] = self.Properties.numprocs
        new_properties['NODE_PROCS'] = self.Properties.nodeprocs
        new_properties['PROCS'] = self.Properties.procs
        new_properties['NUM_THREADS'] = self.Properties.numthreads
        new_properties['PPN_USED'] = self.Properties.ppnused
        new_properties['NUM_SMT'] = self.Properties.numsmt
        new_properties['MEMORY'] = self.Properties.memory
        new_properties['CPUFREQ'] = self.Properties.cpufreq
        
        if debug:
            new_properties['RUNDEBUG'] = 1
            new_properties['DEBUGGER'] = simenv.OptionsManager.GetOption('debugger')
        else:
            new_properties['RUNDEBUG'] = 0
        
        if self.Properties.HasProperty("walltime"):
            walltime_raw = self.Properties.walltime
            Walltime = restartlib.WallTime(walltime_raw)
            
            new_properties['WALLTIME'] = Walltime.Walltime
            new_properties['WALLTIME_HH'] = Walltime.walltime_hh
            new_properties['WALLTIME_MM'] = Walltime.walltime_mm
            new_properties['WALLTIME_SS'] = Walltime.walltime_ss
            new_properties['WALLTIME_SECONDS'] = Walltime.walltime_seconds
            new_properties['WALLTIME_MINUTES'] = Walltime.walltime_minutes
            new_properties['WALLTIME_HOURS'] = Walltime.walltime_hours
            
            # parfile.
            # make sure the parfile has the correct walltime in it
            DefineDatabase.AddReplacement("TerminationTrigger::max_walltime", Walltime.walltime_hours)
            DefineDatabase.AddReplacement("TriggerTerminationManual::max_walltime", Walltime.walltime_hours)
        else:
            DefineDatabase.AddReplacement("TerminationTrigger::max_walltime", "0")
            DefineDatabase.AddReplacement("TriggerTerminationManual::max_walltime", "0")
            
        #new_properties['CHAINED_JOB_ID'] = ''

        for key in DefineDatabase.PresetMacros:
            new_properties[key] = DefineDatabase.PresetMacros.get(key)

 
        self.SimulationLog.Write("Defined substitution properties for execution/run")
        self.SimulationLog.Write(new_properties)
            
        for key in new_properties.keys():
            DefineDatabase.Set(key, new_properties[key])
        
        if parfile != "":
            if parfile.endswith(".rpar"):
                shutil.copy(parfile, self.RestartDir)
                restart_rpar = simlib.BuildPath(self.RestartDir, simlib.BaseName(parfile))
                MakeFileExecutable(restart_rpar)
                ret = simenv.system(restart_rpar)
                if ret != 0:
                    fatal("Error while executing parameter file script %s" % restart_rpar)
                restart_par = re.sub(r'(.*)\.rpar$', r'\1.par', restart_rpar)
                if not os.path.exists(restart_par):
                    fatal("Parameter file script '%s' did not generate a parameter file called '%s'" %(restart_rpar,restart_par))
                par_contents = simlib.GetFileContents(restart_par)
            else:
                par_contents = simlib.GetFileContents(parfile)

            contents = DefineDatabase.SubAll(par_contents)
            simlib.WriteContents(newparpath, contents)

        self.Properties.Save()

        rs = simenv.OptionsManager.GetOption('runscript')
        
        if rs:
            if not os.path.exists(rs):
                fatal("run script '%s' does not exist" % rs)
            runScript = rs
        else:
            if not self.Properties.HasProperty('runscript'):
                fatal("no default runscript defined for simulation %s" % self.SimulationName)
            runScript = self.Properties.runscript
        
        contents = simlib.GetFileContents(runScript)
        contents = DefineDatabase.SubAll(contents)
        
        preparedRunScript = simlib.BuildPath(self.InternalDir, 'RunScript')
        simlib.WriteContents(preparedRunScript, contents)
        MakeFileExecutable(preparedRunScript)

        # Process testsuite script and store it in restart
        if testsuite:
            testsuitefile = os.path.join(self.Properties.sourcedir,"simfactory/bin/RunTestSuite")
            contents = simlib.GetFileContents(testsuitefile)
            contents = DefineDatabase.SubAll(contents)
            testsuitepath = simlib.BuildPath(self.RestartDir, simlib.BaseName(testsuitefile))
            simlib.WriteContents(testsuitepath, contents)
            MakeFileExecutable(testsuitepath)

            self.Properties.AddProperty('testsuitescript', testsuitepath)
        else:
            self.Properties.AddProperty('testsuitescript', "None")

        if testsuite:
            cmd = "%s %s" %(self.Properties.GetProperty("testsuitescript"), preparedRunScript)
        else:
            cmd = preparedRunScript
        
        # change to RestartDir and then run the simulation.
        # TODO -- we're not supposed to do this apparently. 
        os.chdir(self.RestartDir)

        self.SimulationLog.Write("Executing run command: %s" % cmd)

        simlib.ExecuteReplaceCommand(cmd)

        display("Simfactory Done at date: %s" % simenv.system("date")) # This is wrong because system doesn't return the command output
    
    def done(self):
        """This function cleans up the simulation. Specifically, it closes the simulation log file."""

        self.SimulationLog.Close()
        #self.ClearMark()

    def create(self, simulationName, parfile):
        """Create a simulation"""
        # Why is this a part of the SimRestart class?  Why is there no SimSimulation class?

        DefineDatabase = simsubs.DefineDatabase()
        
        # lets start here.
    
        #the calling app (sim/sim-job/sim-run), etc, will determine where simulationName, parfile come from
        self.SimulationName = simulationName
        self.Parfile = parfile
        
        if self.Parfile != "" and not(simlib.FileExists(self.Parfile)):
            fatal("Specified parfile %s does not exist or is not readable" % self.Parfile)
        
        configuration = restartlib.GetConfiguration()
        configPath = simlib.BuildPath(simenv.CONFIGS_PATH, configuration)

        (exe, submitScript, runScript) = restartlib.GetExecutable()
        
        optionlist = simlib.BuildPath(configPath, 'OptionList')
        
        config_id = simlib.GetFileContents(simlib.BuildPath(configPath, 'CONFIG-ID'), 'no-config-id')
        build_id = simlib.GetFileContents(simlib.BuildPath(configPath, 'BUILD-ID'), 'no-build-id')
        
        (self.BaseDir, self.SimulationDir, self.InternalDir, self.CacheDir) = restartlib.CreateRestartSkeleton(self.SimulationName)
        
        #now that we have the simulationdir made, lets attach our log
        self.SimulationLog = restartlib.SimulationLog(self)
        
        self.SimulationLog.Write("Creating simulation %s" % self.SimulationName)
        self.SimulationLog.Write("Simulation directory: %s" % self.SimulationDir)
        
        display("Skeleton Created")
        display("Job directory: \"%s\"" % self.SimulationDir)
        
        machine = simenv.LocalMachine
        machineEntry = simenv.LocalMachineEntry
        
        simlib.VerifyKeys(machineEntry, ['user', 'email'])
        
        user = machineEntry.user
        email = machineEntry.email

        DefineDatabase.Set("USER", user)
        DefineDatabase.Set("EMAIL", email)
        DefineDatabase.Set("MACHINE", simenv.LocalMachine)
        DefineDatabase.Set("CONFIGURATION", configuration)
        DefineDatabase.Set("BASEDIR", self.BaseDir)
        
        localsourcebasedir = simlib.GetLocalSourceBaseDir()
        
        dirsuffix = simlib.GetDirSuffix(localsourcebasedir)
         
        if dirsuffix == None:
            currdir = os.getcwd()
            fatal("The current directory '%s' is not a subdirectory of the sourcebasedir '%s'.\n This error can also be caused if you call Simfactory from outside a Cactus directory, or if this Cactus directory is not a subdirectory of the sourcebasedir." % (currdir, localsourcebasedir))
        
        sourcedir = simlib.BuildPath(localsourcebasedir, dirsuffix)
        
        #dirname == simulationdir
        
        # need to create a simulation id
        simulation_id = restartlib.CreateSimulationId(self.SimulationName)
        
        propertyFile = simlib.BuildPath(self.InternalDir, "properties.ini")
        
        self.Properties = simproperties.SimProperties(baseProperties={})
        self.Properties.init(propertyFile)
        
        self.Properties.AddProperty('machine', machine)
        self.Properties.AddProperty('simulationid', simulation_id)
        self.Properties.AddProperty('sourcedir', sourcedir)
        self.Properties.AddProperty('configuration', configuration)
        self.Properties.AddProperty('configid', config_id)
        self.Properties.AddProperty('buildid', build_id)
        self.Properties.AddProperty('testsuite', simenv.OptionsManager.HasOption('testsuite'))
        if simenv.OptionsManager.HasOption('testsuite'):
            print "Option --testsuite given"
            test = simenv.OptionsManager.GetOption('select-tests')
            if test == None:
                test = "all"
            self.Properties.AddProperty('select_tests', test)
        
        self.Properties.Save()
        
        #roots = ['exe', 'cfg', 'run', 'par']
        
        (exedir, cfgdir, rundir, pardir, datadir) = restartlib.CreateInternalDirs(self.InternalDir)
        
        # runscript
        
        # exe
        ef = simlib.BaseName(exe)
        exefile = simlib.BuildPath(exedir, ef)
        restartlib.CopyFileWithCaching(exe, exedir, simlib.BuildPath(self.CacheDir, 'exe'))
        MakeFileExecutable(exefile)
        display("Executable: \"%s\"" % exe)
        
        self.Properties.AddProperty('executable', exefile)
        
        # config
        cfgfile = simlib.BaseName(optionlist)
        contents = simlib.GetFileContents(optionlist)
        contents = DefineDatabase.SubAll(contents)
        cfgpath = simlib.BuildPath(cfgdir, cfgfile)
        simlib.WriteContents(cfgpath, contents)
        display("Option list: \"%s\"" % cfgpath)
        
        self.Properties.AddProperty('optionlist', simlib.BuildPath(cfgdir, cfgfile))
        
        # submit script
        if submitScript != None:
            queuefile = simlib.BaseName(submitScript)
            contents = simlib.GetFileContents(submitScript)
            contents = DefineDatabase.SubAll(contents)
            submitpath = simlib.BuildPath(rundir, queuefile)
            simlib.WriteContents(submitpath, contents)
            display("Submit script: \"%s\"" % submitpath)
            self.Properties.AddProperty('submitscript', submitpath)
        else:
            self.Properties.AddProperty('submitscript', "None")
            
        # run script
        if runScript != None:
            runfile = simlib.BaseName(runScript)
            contents = simlib.GetFileContents(runScript)
            contents = DefineDatabase.SubAll(contents)
            runpath = simlib.BuildPath(rundir, runfile)
            simlib.WriteContents(runpath, contents)
            display("Run script: \"%s\"" % runpath)
            self.Properties.AddProperty('runscript', runpath)
        else:
            self.Properties.AddProperty('runscript', "None")
            
        # parfile
        if parfile != "":
            par = simlib.BaseName(parfile)
            contents = simlib.GetFileContents(parfile)
            contents = DefineDatabase.SubAll(contents)
            parpath = simlib.BuildPath(pardir, par)
            simlib.WriteContents(parpath, contents)
            display("Parameter file: \"%s\"" % parpath)
        
            self.Properties.AddProperty('parfile', parpath)
        else:
            self.Properties.AddProperty('parfile', "")
        
        self.Properties.Save()
        
        self.SimulationLog.Write("Simulation Properties:")
        self.SimulationLog.Write(self.Properties.toString())
        
        # data directory
        
        if simenv.OptionsManager.HasOption('datadir'):
            datasrc = simenv.OptionsManager.GetOption('datadir')
            if not(os.path.exists(datadir)):
                fatal("could not open data directory \"%s\" for reading" % datasrc)
            
            shutil.copytree(datasrc, datadir, True)
            display("Data Directory: \"%s\"" % datasrc)
        
        self.SimulationLog.Write("Simulation %s created" % self.SimulationName)
        
    def stop(self):
        
        DefineDatabase = simsubs.DefineDatabase()
        
        job_id = self.GetJobId()
        
        if job_id == -1:
            fatal("cannot stop a job without an associated job_id")

        self.SimulationLog.Write("Stopping simulation %s, with restart id %s, with job id %s" % (self.SimulationName, self.RestartID, job_id))

        status = restartlib.GetJobStatus(job_id)
        
        if status == 'U':
            display("job %s already finished or stopped" % job_id)
            self.SimulationLog.Write("Job status determined to be 'U', unable to stop, or already stopped/finished")
            return
        
        machine = simlib.GetMachineName()
        machineEntry = simenv.ConfigurationDatabase.GetMachine(machine)
        
        stopcmd = machineEntry.stop
        
        DefineDatabase.Reset()
        DefineDatabase.Set("JOB_ID", job_id)
        DefineDatabase.Set("USER", machineEntry.user)
        
        # it might be a regex, make sure it's in python format.
        stopcmd = DefineDatabase.SubAll(stopcmd)
        
        force = simenv.OptionsManager.GetOption('force')
        
        term_file = simlib.BuildPath(self.RestartDir, 'TERMINATE')
        
        if os.path.exists(term_file) and not force:
            info("TERMINATE exists, stopping %s gracefully" % job_id)
        
            self.SimulationLog.Write("TERMINATE file exists for job_id %s, terminating gracefully" % job_id)
            # write 1 to term_file
            fptr = open(term_file, 'w+')
            fptr.write('1')
            fptr.close()
            
            self.finish()
            
            return
        
        display("Forcing job %s to stop without using graceful termination" % job_id)
        simlib.ExecuteCommand(stopcmd)

        return          
    
    def show_output(self):
        DefineDatabase = simsubs.DefineDatabase()

        (machine, machineEntry, sourceBaseDir) = simlib.GetLocalEnvironment()

        DefineDatabase.Set("SIMULATION_NAME", self.SimulationName)
        
        job_id = self.GetJobId()
        job_status = 'U'

        if job_id != -1:
            job_status = restartlib.GetJobStatus(job_id)
        
        
        DefineDatabase.Set("JOB_ID", job_id)
        
        parname = simlib.FileBaseName(self.Properties.parfile)
        workdir = simlib.BuildPath(self.RestartDir, parname)
        
        exechost = restartlib.GetExecHost(self)
        
        if exechost == None and job_status == 'R':
            fatal("could not retreive exechost for running simulation %s" % self.SimulationName)
        
        DefineDatabase.Set("EXECHOST", exechost)
        
        output = False
        
        stdout_output = "(file does not exist)"
        stderr_output = "(file does not exist)"
        formaline_output = "(file does not exist)"
        
        formaline_file = simlib.BuildPath(workdir, 'formaline-jar.txt')
        
        if os.path.exists(formaline_file):
            formaline_output = simlib.GetFileContents(formaline_file)
        
        if job_status == 'R':
            
            if simenv.OptionsManager.GetOption('follow'):
                simlib.VerifyKeys(machineEntry, ['stdout-follow'])
                
                followcmd = DefineDatabase.SubAll(machineEntry.GetKey('stdout-follow'))
                followcmd = "cd '%s' && { %s; }" % (self.RestartDir, followcmd)
                                
                simlib.ExecuteReplaceCommand(followcmd)
                
                output = True
            else:
                
                simlib.VerifyKeys(machineEntry, ['stdout', 'stderr'])
                
                stderrcmd = DefineDatabase.SubAll(machineEntry.GetKey('stderr'))
                stdoutcmd = DefineDatabase.SubAll(machineEntry.GetKey('stdout'))
                stdoutcmd = "cd '%s' && { %s; }" % (self.RestartDir, stdoutcmd)
                stderrcmd = "cd '%s' && { %s; }" % (self.RestartDir, stderrcmd)
                
                ff = simenv.popen(stdoutcmd)
                stdout_output = ff.read()
                ret = ff.close()
                
                if ret != None:
                    warning("stdout command \"%s\" returned status %s" % (stdoutcmd, ret))
                    
                ff = simenv.popen(stderrcmd)
                stderr_output = ff.read()
                ret = ff.close()
                
                if ret != None:
                    warning("stderr command \"%s\" returned status %s" % (stderr_output, ret))
        else:
            errfile = simlib.BuildPath(self.RestartDir, "%s.err" % self.SimulationName)
            outfile = simlib.BuildPath(self.RestartDir, "%s.out" % self.SimulationName)
            
            if os.path.exists(outfile):
                stdout_output = simlib.GetFileContents(outfile)
            
            if os.path.exists(errfile):
                stderr_output = simlib.GetFileContents(errfile)
        
        
        if not output:
            # Show stdout, stderr, and Formaline output
            sep = "=" * 80
            display(sep)
            display("The job's Formaline output is:")
            display(sep)
            display(formaline_output)
            display(sep)
            display("The job's stdout is:")
            display(sep)
            display(stdout_output)
            display(sep)
            display("The job's stderr is:")
            display(sep)
            display(stderr_output)
            display(sep)
            
    def cleanup(self):
        # TODO: When is this routine called? It seems to always
        # perform a cleanup.

        # TODO: One could implement the actual cleanup in this
        # routine, and then have "finish" first check whether the
        # simulation should be cleaned up, and if so, call "cleanup".
        # This would separate the test-whether-to-clean-up and
        # clean-up logics, and would remove the need for the
        # "allowForce" flag.
        
        # allow force.
        self.finish(True)
        
    def finish(self, allowForce=False):

        # TODO: this function "cleans up", it does not "finish".
        # terminology is important here, people get easily confused.
        # the messages need to reflect this.
        self.SimulationLog.Write("For simulation %s, Finishing restart %s" % (self.SimulationName, self.LongRestartID))

        # mark the simulation.
        # TODO: what does this mean?
        self.mark()
        
        active_dir = "%s-active" % self.RestartDir

        if not os.path.exists(active_dir):
            info("Simulation %s, restart %s has already been cleaned up -- doing nothing" % (self.SimulationName, self.RestartID))
            return
        
        if self.Properties.machine != simenv.LocalMachine:
            logonly("Error: cannot clean up a simulation created on machine %s from machine %s, skipping" % (self.Properties.machine, simenv.LocalMachine))
            return

        # force is now an option.
        force = simenv.OptionsManager.GetOption('force')
        
        if not allowForce:
            force = False
        
        self.SimulationLog.Write("Force option: %s" % force)
            
        job_id = self.GetJobId()
        job_status = 'U'

        # if the job_id is -1, it means self.GetJobId was unable to determine a valid job id
        # eg, the simulation didn't actually submit itself.
        # if we have a job_id of -1, there's no point in attempting to get the job status
        # since it will always return 'U'.
        
        if job_id != -1:
            job_status = restartlib.GetJobStatus(job_id)
        
        self.SimulationLog.Write("Job ID: %s, Job Status: %s" % (job_id, job_status))
            
        # The check to see if the job is active happens when it polls for the job
        # status using the function restartlib.GetJobStatus. If that returns anything
        # other than 'U', this job won't get cleaned up.

        if job_status == 'E' and not force:
            self.SimulationLog.Write("Job status is E, meaning the job was found in the queue, however, simfactory was unable to determine exactly what its status was.")
            self.SimulationLog.Write("In order to prevent cleaning up a simulation that may still be active, simfactory will leave this job alone.")
            warning("Job status is E: running/queued but unable to determine its exact status, leaving alone.")
            
            return
            
        if job_status != 'U' and not force:
            self.SimulationLog.Write("Job status is not 'U', it is '%s', and force is false, therefor this is a noop, no cleanup required" % job_status)
            return

        # clean up.
        
        self.SimulationLog.Write("Cleaning up simulation %s, restart %s, with job_status %s" % (self.SimulationName, self.RestartID, job_status))
        
        # TODO: Before doing anything to the simulation, a log entry
        # needs to be written.
        
        # step 1. if this simulation is active and Force == False, the first thing we need to do is stop the job.
        if force:
            self.SimulationLog.Write("Forcing cleanup of an active simulation, stopping first")
            self.stop()

        # step 2. remove active flag
        try:
            os.unlink(active_dir)
        except:
            self.SimulationLog.Write("Unable to remove active_dir %s" % active_dir)
            fatal("cannot clean up simulation %s, restart id %s, unable to remove -active symlink %s" % (self.SimulationName, self.RestartID, active_dir))
                
        # step 2. make terminate file not world writable, and fix other permissions
        #if os.path.exists(self.RestartDir):
        #    try:
        #        os.chmod(self.RestartDir, stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH | stat.S_IXOTH)
        #    except:
        #        self.SimulationLog.Write("unable to chmod %s to 755" % self.RestartDir)
        #        pass
                
        term_file = simlib.BuildPath(self.RestartDir, 'TERMINATE')
        #if os.path.exists(term_file):
        #    try:
        #        os.chmod(term_file, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH)
        #    except:
        #        self.SimulationLog.Write("unable to chmod termination file %s, it might not exist" % term_file)
        #        pass
        try:
            perms = os.stat(term_file).st_mode
            perms &= ~stat.S_IWGRP & ~stat.S_IWOTH
            os.chmod(term_file, perms)
        except:
            pass
        
        #outfile = simlib.BuildPath(self.RestartDir, "%s.out" % self.SimulationName)
        #errfile = simlib.BuildPath(self.RestartDir, "%s.err" % self.SimulationName)
        #
        #files = [outfile, errfile]
        #
        #for f in files:
        #    if os.path.exists(f):
        #        try:
        #            os.chmod(f, stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IROTH)
        #        except:
        #            self.SimulationLog.Write("unable to chmod 644 %s, it may not exist" % f)
        #            pass
            
        # [step 4. remove scratch dir]

        # TODO: This looks like a pattern -- calling BuildPath only if
        # a path does not start with a slash. This could be abstracted
        # into a function, e.g. "create-an-absolute-path-name".
            
        # step 5. remove half written checkpoint files
        try:
            simenv.system("find %s -name *.chkpt.tmp.it_*.* -exec rm -rf {} \;" % self.RestartDir)
        except:
            pass
        
        # step 6. Hard-link (Formaline) tarballs between different restarts
        # TODO: move this into subroutine
         
        # Algorithm:
        # 1. Loop over all worthwhile files in the current restart
        # 2.    Loop over all previous restarts (in descending order)
        # 3.       If the equivalent file does not exist: continue
        # 4.       Else, compare files (starting with size!)
        # 5.       If not equal: continue
        # 6.       Else, replace file with hard link
        
        # All restart ids
        rids = restartlib.GetRestartIds(self)
        
        # All files of the current restart
        files = simenv.popen("cd %s && find . -name \*.tar.gz" %
                             self.RestartDir).read().split("\n")
        for f in files:
            path = "%s/%s" % (self.RestartDir, f)
            # Ignore everything but regular files that are not
            # accessed via symbolic links
            if not (os.path.isfile(path) and not os.path.islink(path)):
                continue
            # Only consider large files
            if os.path.getsize(path) < 1000:
                continue
            # Can't do reversed(sorted(rids)) with Python 2.3
            rids.sort()
            rids.reverse()
            for rid in rids:
                oldpath = "%s/output-%04d/%s" % (self.SimulationDir, rid, f)
                # Ignore everything but regular files that are not
                # accessed via symbolic links
                if not (os.path.isfile(oldpath) and
                        not os.path.islink(oldpath)):
                    continue
                # Ignore the old file if it is the same file
                if os.path.samefile(path, oldpath):
                    continue
                # Do nothing unless the files are identical
                if not filecmp.cmp(path, oldpath):
                    continue
                # Find a temporary name that does not yet exist
                tmppath = "%s.tmp" % path
                if os.path.exists(tmppath):
                    continue    # bail out
                # Replace file by hard link
                try:
                    os.link(oldpath, tmppath)
                    os.rename(tmppath, path)
                    break       # next file
                except OSError, e:
                    # Clean up
                    try:
                        os.unlink(tmppath)
                    except OSError, e:
                        pass
        
        self.SimulationLog.Write("Simulation %s, restart %s, with job id %s has been successfully cleaned up" % (self.SimulationName, self.RestartID, job_id))    
        # TODO: The log entry should now state that the simulation has
        # been successfully deactivated.
        
    
    def archive(self):
        
        machineEntry = simenv.LocalMachineEntry
        
        simlib.VerifyKeys(machineEntry, ['archivetype'])
        
        archiveType = machineEntry.archivetype
        
        ArchiveEngine = simarchive.SimArchive(archiveType, self)
        ArchiveEngine.authenticate()
        ArchiveEngine.store()
        
    def GetJobId(self):
    
        # transitioning between having the jobid in two places.
        # look for the new place. if it has it in the properites file,
        # use it. otherwise, look for job.ini
        if self.Properties.HasProperty('jobid'):
            return self.Properties.jobid
        
        job_file = simlib.BuildPath(self.InternalDir, 'job.ini')
        
        if not(os.path.exists(job_file)):
            return -1
            
        pjob = simproperties.SimProperties()
        pjob.init(job_file)

        return pjob.GetProperty("jobid")
        
    def attachLog(self, simulationdir):
    
        # noop
        
        pass
        return
        
    # all of these will be flushed out as i proceed.
    def trash(self):    
        (machine, machineEntry, sourceBaseDir) = simlib.GetLocalEnvironment()
        
        simlib.VerifyKeys(machineEntry, ['basedir'])
        
        basedir = simlib.GetBaseDir(machineEntry)
        
        trashPath = simlib.BuildPath(basedir, 'TRASH')
        
        if not(os.path.exists(trashPath)):
            try:
                simenv.system("mkdir -p %s" % trashPath)
            except OSError, e:
                fatal("could not make trash path \"%s\", %s" % (trashPath, e))
        
        # okay, we have trash path. lets find out if any of our restarts are running.
        
        rids = restartlib.GetRestartIds(self)
        
        for rid in rids:
            rr = restartlib.GetRestartByRestartId(self.SimulationName, rid)    
            
            job_id = rr.GetJobId()
            
            if job_id == -1:
                continue
            
            job_status = restartlib.GetJobStatus(job_id)
            
            if job_status != 'U':
                fatal("Error: Simulation %s, with Restart ID %s is either queued, holding, or running.\nJob ID %s must be stopped before a purge can happen." % (self.SimulationName, rid, job_id))
                
        # okay, we have no running/queued/holding jobs.
        # lets move the simulation to the trash 
        
        trashFolder = simlib.BuildPath(trashPath, self.Properties.simulationid)
        
        if not(os.path.exists(trashFolder)):
            try:
                simenv.system("mkdir -p %s" % trashFolder)
            except OSError, e:
                fatal("could not make trash folder \"%s\", %s" % (trashFolder, e))
        
        shutil.move(self.SimulationDir, trashFolder)
        
        display("Simulation %s has been moved to trash folder %s" % (self.SimulationName, trashFolder))
        
        return
