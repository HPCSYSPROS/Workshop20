# simlib.py -- shared necessary library functions across all sim-* utilties
#
#
#
#

import os, sys, re
import socket
import os.path
import errno
import simsubs
import simenv

from libutil import *

def GetHostNameAlias():
    
    # first check if our hostname was specified via an option
    if simenv.OptionsManager.HasOption('hostname'):
        return simenv.OptionsManager.GetOption('hostname')
    
    # get the users home directory, check for ~/.hostname
    homedir = os.environ['HOME']
    home_hostname = "%s%s%s" % (homedir, os.sep, ".hostname")
    if FileExists(home_hostname):
        hostname = GetFileContents(home_hostname).strip()
        return hostname
    
    # lets get the hostname
    fd = simenv.popen("hostname")
    hostname = fd.read().strip()
    fd.close()
    
    try:
        (name, aliases, _) = socket.gethostbyname_ex(hostname)
    except:
        # warning("could not resolve hostname %s" % hostname)
        return hostname
    
    # out of all our names, find one that has a dot, if possible.
    names = aliases
    names.insert(0, name)
    names.insert(0, hostname)
    names = filter(lambda n: not (re.match("localhost", n)), names)
    
    for n in names:
        if n.count("."):
            return n

    if len(name) > 0:
        return name
    
    return hostname

def BoolToYesNo(val):

    if val:
        return 'yes'
    
    return 'no'
    

def GetLocalEnvironment():
    
    machine = simenv.LocalMachine
    machineEntry = simenv.LocalMachineEntry        
    
    req_keys = ['user', 'sourcebasedir']
    
    VerifyKeys(machineEntry, req_keys)
    
    sourceBaseDir = GetLocalSourceBaseDir()

    return (machine, machineEntry, sourceBaseDir)

def GetMachineName():

    if simenv.OptionsManager.HasOption("machine"):
        return simenv.OptionsManager.GetOption("machine")
        
    machineName = GetHostNameAlias()
    
    machines = simenv.ConfigurationDatabase.GetMachines()
    
    found_machines = []
    
    for m in machines:
        if m == machineName:
            found_machines.append(m)
        else:
            entry = simenv.ConfigurationDatabase.GetMachine(m)

            if entry.HasKey('hostname') and entry.hostname == machineName:
                found_machines.append(m)
            else:
                if entry.HasKey("aliaspattern"):
                    if Matches(entry.aliaspattern, machineName):
                        found_machines.append(m)
    
    error = None
    if len(found_machines) == 0:
        error = "Unknown machine name %s" % machineName

    if len(found_machines) > 1:
        error = "Could not identify machine %s -- possible matches are %s" % (machineName, found_machines)
    
    if error != None:
        display(error)
            
        return None
        
    return found_machines[0]


def MachineExists(machineName):

    machines = simenv.ConfigurationDatabase.GetMachines()
    
    found_machines = []
    
    for m in machines:
        if m == machineName:
            found_machines.append(m)
        else:
            entry = simenv.ConfigurationDatabase.GetMachine(m)

            if entry.HasKey('hostname') and entry.hostname == machineName:
                found_machines.append(m)
            else:
                if entry.HasKey("aliaspattern"):
                    if Matches(entry.aliaspattern, machineName):
                        found_machines.append(m)
    
    if len(found_machines) == 0:
        return False
        
    if len(found_machines) > 1:
        return False

    return True
    
# Get source base directory
def GetLocalSourceBaseDir():

    if simenv.OptionsManager.HasOption('sourcebasedir'):
        sourcebasedir = simenv.OptionsManager.GetOption('sourcebasedir')
    else:
        return GetSourceBaseDir(simenv.LocalMachineEntry)
        
    if not(sourcebasedir[0] == "/"):
        cwd = os.getcwd()
        sourcebasedir = "%s%s%s" % (cwd, os.sep, sourcebasedir)
    
    if not(os.path.exists(sourcebasedir)):
        warning('Cannot access source base directory "%s"' % sourcebasedir)
    
    if sourcebasedir.count("//") > 0:
        return sourcebasedir[1:]
   
    return sourcebasedir

def BaseName(filename):
    (path, ff) = os.path.split(filename)
    return ff

def FileBaseName(filename):

    # simulate rsplit, because rsplit doesn't exist
    # in python 2.3
    parts = BaseName(filename).split(".")
    
    parts.pop()
    return ".".join(parts)
    
    return file

    # TODO: Should probably move this function to restartlib
def ParseSimulationCommandLine():
    """Derive the parfile, procs and walltime options from the command
    line if they have been given as positional arguments, and return
    the simulation name."""

    parfile = None
    simulationName = None
    
    #first, lets count how many arguments we have.
    
    numArgs = len(simenv.OptionsManager.args)
    
    if numArgs < 3:
        if numArgs == 1:
            argument = simenv.OptionsManager.args.pop(0)
            if argument.count(".par") > 0 or argument.count(".rpar") > 0:
                parfile = argument
                info("Parfile: %s" % parfile)
                simenv.OptionsManager.UpdateOption('parfile', parfile)
                return FileBaseName(parfile)
            return argument
        else:
            return None

    # numArgs >= 3
        
    argument = simenv.OptionsManager.args.pop(0)
    
    # 4 args -- simname, parfile, procs, walltime
    if numArgs == 4:
        simulationName = argument
        argument = simenv.OptionsManager.args.pop(0)

    if argument.endswith(".par") or argument.endswith(".rpar"):
        parfile = argument
        info("Parfile: %s" % parfile)
        simenv.OptionsManager.UpdateOption('parfile', parfile)
        
        if simulationName == None:
            simulationName = FileBaseName(parfile)
    else:
        fatal("Error: argument stream in incorrect format.\nExpecting: %s %s [simulationName] parfile.[r]par procs walltime" % (sys.argv[0], simenv.COMMAND))

    display("Simulation Name: %s" % simulationName)
    
    argument = simenv.OptionsManager.args.pop(0)
    simenv.OptionsManager.UpdateOption('procs', int(argument))
    
    info("Procs: %s" % argument)
    
    argument = simenv.OptionsManager.args.pop(0)
    simenv.OptionsManager.UpdateOption('walltime', argument)
    
    info("Walltime: %s" % argument)
    
    return simulationName
    
# Return the local basedir; the argument "machineEntry" is ignored
def GetBaseDir(machineEntry):
    
    DefineDatabase = simsubs.DefineDatabase()
    
    if simenv.OptionsManager.HasOption('basedir'):
        basedir = simenv.OptionsManager.GetOption('basedir')
    else:
        machineEntry = simenv.LocalMachineEntry

        VerifyKeys(machineEntry, ['user', 'basedir'])
        
        user = machineEntry.user
    
        DefineDatabase.Set('USER', user)
        basedir = DefineDatabase.SubAll(machineEntry.basedir)
        
    if not(basedir[0] == "/"):
        cwd = os.getcwd()
        basedir = "%s%s%s" % (cwd, os.sep, basedir)
    
    #if not(os.path.exists(basedir)):
    #   dprint('Cannot access source base directory "%s"' % basedir)
    
    return basedir

# Return the basedir
def GetBaseDir1(machineEntry):
    DefineDatabase = simsubs.DefineDatabase()
    
    if not(machineEntry.HasKey('user')):
        machineName = machineEntry.nickname
        fatal("machine entry %s does not have a user defined" % machineName)
    
    user = machineEntry.user
    
    DefineDatabase.Set('USER', user)
    basedir = DefineDatabase.SubAll(machineEntry.basedir)
    
    if not(basedir[0] == "/"):
        cwd = os.getcwd()
        basedir = "%s%s%s" % (cwd, os.sep, basedir)
    
    if basedir.count("//") > 0:
        return basedir[1:]
    
    return basedir

def GetSourceBaseDir(machineEntry):
    DefineDatabase = simsubs.DefineDatabase()
    
    if not(machineEntry.HasKey('user')):
        machineName = machineEntry.nickname
        fatal("machine entry %s does not have a user defined" % machineName)

    user = machineEntry.user
    
    DefineDatabase.Set('USER', user)
    sourcebasedir = DefineDatabase.SubAll(machineEntry.sourcebasedir)
        
    if not(sourcebasedir[0] == "/"):
        cwd = os.getcwd()
        sourcebasedir = "%s%s%s" % (cwd, os.sep, sourcebasedir)
    
    if sourcebasedir.count("//") > 0:
        return sourcebasedir[1:]

    return sourcebasedir

def GetDirSuffix(prefix):
    
    if not(os.path.exists(prefix)):
        fatal("Configured sourcebasedir \"%s\" does not exist or is not readable" % prefix)

    curdir = os.getcwd()
    os.chdir(prefix)
    
    real_prefix = os.getcwd()
    os.chdir(curdir)
    
    pattern = "^%s/(.*)$" % real_prefix
    
    p = re.compile(pattern)
    matches = p.match(curdir)
    
    
    if matches == None:
    #    error = "Called from the wrong location."
    #    error = "%s\n%s" % (error, "Current directory is '%s'" % curdir)
    #    error = "%s\n%s" % (error, "but expected a subdirectory of '%s'." % prefix)
    #    error = "%s\n%s" % (error, "It is also be possible that you need to correct your 'sourcebasedir' entry")
    #    error = "%s\n%s" % (error, "in the mdb entry for this machine.")
    #    warning(error)
        
        return None
           
    return matches.group(1)

def GetCactusPath():

    # look for CACTUS_PATH
    if os.environ.has_key('CACTUS_PATH'):
        env_path = os.environ['CACTUS_PATH']
        
        if not(os.path.exists(BuildPath(env_path, "src", "interface.ccl"))):
            fatal("CACTUS_PATH %s is not a valid Cactus source location" % env_path)
        else:
            #dprint("Cactus path env: %s" % env_path)
            return env_path
        

    # Let's see if PWD is defined and take this is set. This has the
    # advantage that symlinks are not dereferenced.
    current_path = os.environ.get("PWD")
    
    # lets take a look at os.cwd() if we either could not find PWD or if it
    # does not point to the actual current working directory
    if current_path is None or \
       os.path.abspath(os.path.realpath(current_path)) != os.getcwd():
        current_path = os.getcwd()
    
    #dprint("DEBUG: current_path: %s" % current_path)
    
    if os.path.exists(BuildPath(current_path, "src", "interface.ccl")):
        #dprint("DEBUG: Cactus path 2: %s" % current_path)
        return current_path
    
    # alright, still no luck. Lets see if we're in a subdirectory of the Cactus source tree
    
    
    parts = current_path.split(os.sep)
    
    current_path = ""
    for part in parts:
        if len(current_path) == 0:
            current_path = "%s" % part
        else:
            current_path = "%s%s%s" % (current_path, os.sep, part)
        
        if os.path.exists(BuildPath(current_path, "src", "interface.ccl")):
            return current_path
        
    # try using BASE_PATH to determine whether or not simfactory is 
    
    # installed in the Cactus source tree.
    basepath = simenv.BASE_PATH
    
    #dprint("DEBUG: basepath: %s" % basepath)
    
    parts = basepath.split(os.sep)
    parts.pop()
    
    path = os.sep.join(parts)
    
    #dprint("DEBUG: after pop, path is: %s" % path)
    
    if os.path.exists(BuildPath(path, "src", "interface.ccl")):
        #dprint("DEBUG: Cactus path 1: %s" % path)
        return path
        
    return None
    
def GetUsername():
    
    if "USER" not in os.environ:
        fd = simenv.popen("whoami")
        username = fd.read().strip()
        fd.close()
        
    else:
        username = os.environ["USER"]
    
    return username

def GetSimulations():
    pass

    # need basedir, simulationdir, and internaldir
    machineEntry = simenv.LocalMachineEntry
    
    base_dir = GetBaseDir(machineEntry)
    
    #simulation_dir = BuildPath(base_dir, SimulationName)
    sims = []
    
    if not(os.path.exists(base_dir)):
        return sims
    
    for simulation in os.listdir(base_dir):
        simpath = BuildPath(base_dir, simulation)
        if os.path.isdir(simpath):
            prop_path = BuildPath(base_dir, simulation, "SIMFACTORY", "properties.ini")
            if simulation not in ['CACHE', 'TRASH'] and os.path.exists(prop_path):
                sims.append(simulation)
    
    #dprint("Simulations: %s" % sims)
    
    return sims

def GetConfigurations():

    config_dir = BuildPath(simenv.CACTUS_PATH, "configs")
    
    ll = []

    if not(os.path.exists(config_dir)):
        return ll
    
    for item in os.listdir(config_dir):
        itempath = BuildPath(config_dir, item)
        if os.path.isdir(itempath):
            ll.append(item)
    
    #dprint("Configurations: %s" % ll)
    return ll

# Construct an ssh command, i.e. a command that is executed
# remotely (via ssh) from a given local command.

# The command constructed here must work as prefix (for
# rsync), therefore we cannot call /bin/bash, since this only
# executes its first argument. We cannot use quoting either,
# since there is no way to close the quotes. We therefore
# ignore sshsetup and friends, and have to assume that
# trampolines work without it. (Maybe we should ensure here
# that these entries are indeed empty.)

# If cmd is set, assume this command should be executed.
# If cmd is None, assume that it will be used as rsync prefix,
# and omit the username@hostname part.
def GetSSHCommand(from_machine, machine, cmd, other_opts='', remotePath=None):
    if not machine:
        fatal('machine is None')
        #return cmd
    
    DefineDatabase = simsubs.DefineDatabase()
    
    local_sourcebasedir = GetLocalSourceBaseDir()
    sourcename = GetDirSuffix(local_sourcebasedir)
    if sourcename == None:
        sourcename = ''
    
    # Find out from where we should connect
    trampoline = GetTrampoline(machine)
    if trampoline != None:
        if not simenv.ConfigurationDatabase.HasMachineEntry(trampoline):
            fatal('Unknown trampoline "%s" specified for machine "%s"' % (trampoline, machine))
    else:
        trampoline = from_machine

    from_entry = simenv.ConfigurationDatabase.GetMachine(trampoline)
    VerifyKeys(from_entry, ['user', 'localsshprefix', 'localsshopts', 'sourcebasedir'])
    from_user          = from_entry.user
    from_sshprefix     = from_entry.localsshprefix
    from_sshopts       = from_entry.localsshopts
    from_sourcebasedir = from_entry.sourcebasedir
    
    from_sourcedir = os.path.join(from_sourcebasedir, sourcename)
    DefineDatabase.Set('SOURCEDIR', from_sourcedir)
    DefineDatabase.Set('USER', from_user)
    from_sshprefix = DefineDatabase.SubAll(from_sshprefix)
    from_sshopts = DefineDatabase.SubAll(from_sshopts)
    
    entry = simenv.ConfigurationDatabase.GetMachine(machine)
    VerifyKeys(entry, ['hostname', 'user', 'sshcmd', 'sshopts', 'sourcebasedir', ])
    hostname      = entry.hostname
    user          = entry.user
    sshcmd        = entry.sshcmd
    sshopts       = entry.sshopts
    sourcebasedir = entry.sourcebasedir
    sourcedir = os.path.join(sourcebasedir, sourcename)
    DefineDatabase.Set('SOURCEDIR', sourcedir)
    DefineDatabase.Set('USER', user)
    
    if cmd != None:
        if remotePath != None:
            cmd = "cd %s && { %s; }" % (remotePath, cmd)
        cmd = "%s@%s %s" % (user, hostname, QuoteSafe(cmd))
    else:
        cmd = ''
    cmd = "%s %s %s %s %s %s" % (from_sshprefix, sshcmd, from_sshopts, sshopts, other_opts, cmd)
    cmd = DefineDatabase.SubAll(cmd)
    
    if trampoline != from_machine:
        cmd = GetSSHCommand(from_machine, trampoline, cmd, other_opts)
    
    return cmd

def ExecuteCommand(command, output=False):

    sys.stdout.flush()
    sys.stderr.flush()
    
    if command == None or len(command) == 0:
        return 0

    envsetup = GetMachineOption('envsetup')
    
    DefineDatabase = simsubs.DefineDatabase()
    DefineDatabase.Set("USER", simenv.LocalMachineEntry.user)
    
    envsetup = DefineDatabase.SubAll(envsetup)
    
    command = "{ %s; } && { %s; }" % (envsetup, command)
    
    if output:
        command = "{ %s; } 2>&1" % command
    
    command = "/bin/bash -c %s" % (QuoteSafe(command))
    
    if (not output) or simenv.VERBOSE:
        info("EXECUTING COMMAND: %s\n" % command)
    else:
        logonly("EXECUTING COMMAND: %s\n" % command)
        
    if not output:
        ret = os.system(command)
        retval = ret
    else:
        fd = os.popen(command)
        output = fd.read()
        ret = fd.close()
        if ret == None:
            ret = 0
        retval = (output, ret)

    if ret != 0:
        dprint("Error while executing command \"%s\":" % command)
        if output:
            dprint("Command output: %s" % output)
        signal = ret & 0xff
        status = ret >> 8
        if signal != 0:
            dprint("   Command was killed by signal %d" % signal)
        else:
            dprint("   Command returned exit status %d" % status)
        
    return retval

def ExecuteReplaceCommand(command, output=False):
    
    if (not output) or simenv.VERBOSE:
        info("EXECUTING: %s\n" % command)
    else:
        logonly("EXECUTING: %s\n" % command)
    envsetup = GetMachineOption('envsetup')
    
    DefineDatabase = simsubs.DefineDatabase()
    DefineDatabase.Set("USER", simenv.LocalMachineEntry.user)
    
    envsetup = DefineDatabase.SubAll(envsetup)
    
    command = "{ %s; } && %s" % (envsetup, command)
        
    sys.stdout.flush()
    sys.stderr.flush()

    pid = os.fork()
    
    if pid > 0:
        retinfo = os.wait()
        return retinfo[1]
    else:
        args = ['bash', '-c', command]
        logonly("executed: /bin/bash -c %s" % " ".join(args))
        os.execv('/bin/bash', args)
    
def GetArguments():
    args = []
    if simenv.OptionsManager.HasOption("argument"):
        args.append(simenv.OptionsManager.GetOption('argument'))
    
    return args
    
def QuoteSafe(cmd):
    return "'%s'" % cmd.replace("'", r"'\''")
    
def VerifyKeys(machineEntry, req_keys):
    for key in req_keys:
        if not(machineEntry.HasKey(key)):
            fatal("machine %s is missing a required key: %s" % (machineEntry.EntryName, key))
    
def GetIOMachine(machineName):
    return GetMachineByKey(machineName, 'iomachine')

def GetTrampoline(machineName):
    return GetMachineByKey(machineName, 'trampoline', True)
    
def GetMachineByKey(machineName, key, useNone=False):
    entry = simenv.ConfigurationDatabase.GetMachine(machineName)
    
    if not(entry.HasKey(key)):
        if useNone:
            return None
        else:
            return machineName

    mm = entry.GetKey(key)
    
    if len(mm.strip()) == 0:
        return None
    
    if not(simenv.ConfigurationDatabase.HasMachineEntry(mm)):
        fatal("specified %s %s for machine %s does not exist" % (key, mm, machineName))
    
    return mm

def GetDefaultConfiguration():

    dEntry = simenv.ConfigurationDatabase.GetConfigOption("default-configuration-name")
    if dEntry == None:
        configName = "sim"
    else:
        configName = dEntry
        
    config = configName
    
    #opts = ['debug', 'optimise', 'unsafe', 'profile']
    #
    #for opt in opts:
    #   if simenv.OptionsManager.RawOptionDefined(opt):
    #       config = "%s-%s" % (config, opt)

    debug    = simenv.OptionsManager.GetOption('debug'   )
    #optimise = simenv.OptionsManager.GetOption('optimise')
    optimise = not debug
    unsafe   = simenv.OptionsManager.GetOption('unsafe'  )
    profile  = simenv.OptionsManager.GetOption('profile' )

    if debug:
            config = '%s-debug' % config
    if optimise == debug:
        if optimise:
            config = '%s-optimise' % config
        else:
            config = '%s-nooptimise' % config
    if unsafe:
        config = '%s-unsafe' % config
    if profile:
        config = '%s-profile' % config

    return config
    
def FileExists(filename):
    return os.path.exists(filename)
    
def DirectoryExists(dirname):
    return os.path.isdir(dirname)

    # TOOD: How to make this abort if there is an error?
def GetFileContents(filename, default="", suppressWarnings=True):
    if filename == None:
        return default

    try:
        fptr = open(filename, "r")
        contents = fptr.read()
        fptr.close()
        return contents
    except:
        if not suppressWarnings:
            warning("could not open %s for reading" % filename)
        
        return default

def WriteContents(filename, contents):
    try:
        fptr = open(filename, "w")
    except IOError, e:
        fatal("Could not open file \"%s\" for writing:\n%s" % (filename, e))
    
    fptr.write(contents)
    fptr.close()
    

def RemoveFile(filename):
    try:
        os.unlink(filename)
    except OSError:
        return

def RenameFile(old, new):
    try:
        os.rename(old, new)
    except OSError:
        return

def MakeDirs(path, exist_ok = True):
    """ Creates all directories in path. If exist_ok is True (default) then no
    exception is raised if (part of) the path already exists.
    """
    try:
        os.makedirs(path)
    except OSError, exc:
        if exist_ok and exc.errno == errno.EEXIST and DirectoryExists(path):
            pass
        else:
            raise

def GetVersion(fileContents):
    # TODO: Can we use re.DOTALL, and remove \n from the pattern?
    regex = r'^\s*VERSION\s*=\s*([^\n#;]*)$'
    rr = re.compile(regex, re.MULTILINE)
    m = rr.search(fileContents)
    if m != None:
        return m.group(1)
    else:
        return ""
    
    
def BuildPath(*args):
    
    # support the old way
    if len(args) == 1 and isinstance(args[0], list):
        return os.sep.join(args[0])
            
    return os.sep.join(args)

def ReplaceLeadingPath(path, old, new):
    """ looks for the sub-path old in path and replaces it by new.
Normalizes path and old to handle symbolic links gracefully."""

    abs_path = os.path.abspath(os.path.realpath(path))
    abs_old = os.path.abspath(os.path.realpath(old))
    if abs_path.startswith(abs_old):
      return abs_path.replace(abs_old, new, 1)
    else:
      return abs_path
    
def GetOptionList(required=False):
    return GetCactusFile(simenv.MDB_BASE_PATH, "optionlist", required)

def GetThornList(required=False):
    tl = GetCactusFile(simenv.ETC_PATH, "thornlist", False)
    if required and tl == None:
        fatal("No thornlist supplied.  Use --thornlist <thornlist> on the command line or set thornlist = <thornlist> in the [default] section of defs.local.ini.")
    else:
        return tl
    
def GetThornListContents(filename):
    machineEntry = simenv.LocalMachineEntry

    contents = GetFileContents(filename, None)
    if contents == None:
        fatal("could not open file '%s' for reading" % filename)

    defineDatabase = simsubs.DefineDatabase()
    enabledThorns = (machineEntry.GetKey('enabled-thorns') + "\n" +
                     machineEntry.GetKey('enabled-thorns-default') + "\n" +
                     machineEntry.GetKey('enabled-thorns-local') + "\n")
    for thorn in enabledThorns.split():
        defineDatabase.AddSubstitution(r'^[ \t]*#DISABLED[ \t]+(%s)\b' % thorn,
                                       '@1@')
    contents = defineDatabase.SubAll(contents)

    defineDatabase = simsubs.DefineDatabase()
    disabledThorns = (machineEntry.GetKey('disabled-thorns') + "\n" +
                      machineEntry.GetKey('disabled-thorns-default') + "\n" +
                      machineEntry.GetKey('disabled-thorns-local') + "\n")
    for thorn in disabledThorns.split():
        defineDatabase.AddSubstitution(r'^[ \t]*(%s)\b' % thorn,
                                       '#DISABLED @1@')
    contents = defineDatabase.SubAll(contents)

    return contents
    
def GetParFile(required=False):
    return GetCactusFile(simenv.ETC_PATH, "parfile", required)

def GetSubmitScript(required=False):
    return GetCactusFile(simenv.MDB_BASE_PATH, "submitscript", required)

def GetRunScript(required=False):
    return GetCactusFile(simenv.MDB_BASE_PATH, "runscript", required)

def GetMachineOption(option, defaultValue=None):
    value = simenv.OptionsManager.GetOption(option)
    
    if value == None:
        value = simenv.LocalMachineEntry.GetKey(option)
        
    if value == None:
        value = defaultValue
        info("for machine option \"%s\", using default value \"%s\"" % (option, value))

    return value

def GetProcs(existingProperties):
    
    machineEntry = simenv.LocalMachineEntry
    
    maxnodes = int(machineEntry.GetKey("nodes"))
    minppn   = int(machineEntry.GetKey("min-ppn"))
    maxppn   = int(machineEntry.GetKey("ppn"))
    
    # Total number of threads
    procs = simenv.OptionsManager.GetOption('procs')
    if not procs and existingProperties:
        procs = existingProperties.procs
    if not procs:
        procs = 1
    try:
        procs = int(procs)
    except ValueError:
        fatal("Total number of threads (procs=%s) is not an integer" % procs)
    
    if procs < 1:
        fatal("Illegal total number of threads specified: specified procs=%d" %
              procs)
    
    # Number of requested cores per node
    if existingProperties != None:
        ppn = int(GetMachineOption('ppn', existingProperties.ppn))
    else:
        ppn = int(GetMachineOption('ppn'))
    if ppn == None:
        ppn = maxppn            # default
        if ppn > procs:
            ppn = procs
    if ppn < minppn or ppn > maxppn:
        fatal("Illegal number of requested cores per node specified: "
              "specified ppn=%s (min-ppn is %s, max-ppn is %s)" %
              (ppn, minppn, maxppn))
    
    # Number of used cores per node
    ppn_used = simenv.OptionsManager.GetOption('ppn-used')
    if ppn_used == None:
        if existingProperties != None:
            ppn_used = existingProperties.ppnused
        else:
            ppn_used = ppn
    try:
        ppn_used = int(ppn_used)
    except ValueError:
        fatal("Number of used cores per node (ppn-used=%s) is not an integer" %
              ppn_used)
    if ppn_used < 1:
        fatal("Illegal number of used cores per node specified: "
              "specified ppn-used=%s" % ppn_used)
    if ppn_used > ppn:
        # Oversubscription
        warning("Too many used cores per node specified: "
                "specified ppn-used=%d (ppn is %d)" % (ppn_used, ppn))
    
    # Number of threads per process
    suggested_num_threads = int(GetMachineOption('num-threads'))
    num_threads = simenv.OptionsManager.GetOption('num-threads')
    if num_threads == None:
        if existingProperties != None:
            num_threads = existingProperties.numthreads
        else:
            num_threads = GetMachineOption('num-threads')
    try:
        num_threads = int(num_threads)
    except ValueError:
        warning("Number of threads per process (num-threads=%s) "
                "is not an integer" % num_threads)
    if num_threads < 1:
        fatal("Illegal number of threads per process specified: "
              "specified num-threads=%d" % num_threads)
    if num_threads > procs:
        warning("Too many threads per process specified: "
                "specified num-threads=%d (ppn-used is %d)" %
                (num_threads, ppn_used))
    
    # Number of SMT threads per core
    suggested_num_smt = int(GetMachineOption('num-smt'))
    num_smt = simenv.OptionsManager.GetOption('num-smt')
    if num_smt == None:
        if existingProperties != None:
            try:
                num_smt = existingProperties.numsmt
            except AttributeError:
                num_smt = GetMachineOption('num-smt')
        else:
            num_smt = GetMachineOption('num-smt')
    try:
        num_smt = int(num_smt)
    except ValueError:
        warning("Number of SMT threads per core (num-smt=%s) "
                "is not an integer" % num_smt)
    
    if num_smt < 1:
        fatal("Illegal number of SMT threads per core specified: "
              "specified num-smt=%d" % num_smt)
    if num_smt > num_threads:
        warning("Too many SMT threads per process specified: "
                "specified num-smt=%d (num-threads is %d)" %
                (num_smt, num_threads))
    
    # Check consistency

    if ppn_used * num_smt % num_threads != 0:
        # More threads than used cores
        warning ("Number of used cores per node, number of SMT threads, "
                 "and number of threads per process are inconsistent: "
                 "ppn-used=%d, num-smt=%d, num-threads=%d_threads "
                 "(ppn-used*num-smt must be "
                 "an integer multiple of num-threads)" %
                 (ppn_used, num_smt, num_threads))
    
    if procs % num_threads != 0:
        # "Uneven" number of threads; there will be too many threads
        warning ("Total number of threads and number of threads per process "
                 "are inconsistent: procs=%d, num-threads=%d "
                 "(procs*num-smt must be an integer multiple of num-threads)" %
                 (procs, num_threads))
    
    if procs % ppn_used != 0:
        # "Uneven" number of threads; some cores will be idle
        warning ("Total number of threads and number of cores per node "
                 "are inconsistent: procs=%d, ppn-used=%d "
                 "(procs must be an integer multiple of ppn-used)" %
                 (procs, ppn_used))
    
    # Calculate dependent quantities
    
    # Number of MPI processes
    num_procs = (procs + num_threads - 1) / num_threads
    if num_procs < 1:
        num_procs = 1
    
    # Number of nodes
    nodes = (procs + (ppn_used * num_smt) - 1) / (ppn_used * num_smt)
    if nodes < 1:
        nodes = 1
    if nodes > maxnodes:
        fatal("Too many nodes specified: nodes=%d (maxnodes is %d)" %
              (nodes, maxnodes))
    
    # Total number of requested cores
    procs_requested = nodes * ppn
    # Number of MPI processes per node
    node_procs = ppn_used * num_smt / num_threads
    
    info("Job allocation information:")
    info("   System:       nodes=%d cores/node=%d threads/process=%d "
         "threads/core=%d" %
         (maxnodes, maxppn, suggested_num_threads, suggested_num_smt))
    info("   Requested:    nodes=%d cores=%d cores/node=%d" %
         (nodes, procs_requested, ppn))
    info("   Run:          processes=%d threads=%d threads/process=%d "
         "threads/core=%d" %
         (num_procs, procs, num_threads, num_smt))
    info("   Distribution: processes/node=%d threads/node=%d" %
         (node_procs, ppn_used))
    info("   Ratio:        threads/core=%.3f cores/thread=%.3f" %
         (1.0*procs/procs_requested, 1.0*procs_requested/procs))
    
    return (nodes, procs_requested, ppn, num_procs, node_procs,
            procs, num_threads, ppn_used, num_smt)

def GetCactusFile(base_path, option, required=False):
    if simenv.OptionsManager.HasOption(option):
        ff = simenv.OptionsManager.GetOption(option)

        if len(ff) == 0 or ff == None:
            if required:
                fatal("specified %s is empty or None" % option)
            else:
                return None

        ff = os.path.expanduser(ff)

        if not(os.path.exists(ff)):

            if not(ff.startswith("/")):
                
                # lets see if it exists relative to the current directory
                cwd = os.getcwd()
                eff = BuildPath(cwd, ff)
                if os.path.exists(eff):
                    ff = eff
                else:
                    folder = "%ss" % option
                    sff = BuildPath(base_path,folder,ff)

                    if not(os.path.exists(sff)):
                        fatal("specified %s \"%s\" does not exist or is not readable\n   tried paths \"%s\"\n   and \"%s\",\n   working directory is \"%s\"" % (option, ff, eff, sff, cwd))
                    else:
                        ff = sff

        if not(os.path.exists(ff)):
            fatal("specified %s %s does not exist or is not readable" % (option, ff))

        return os.path.abspath(ff)

    machineEntry = simenv.LocalMachineEntry

    if required:
        # make sure option is set
        VerifyKeys(machineEntry, [option])

    ff = machineEntry.GetKey(option)

    if ff == None:
        if required:
            fatal("specified %s %s does not exist or is not readable" % (option, ff))
        else:
            return None

    if ff[0] != os.sep:
        # lets see if it exists relative to the current directory
        cwd = os.getcwd()
        eff = BuildPath(cwd, ff)
        
        if os.path.exists(eff):
            ff = eff
        else:
            folder = "%ss" % option
            eff = BuildPath(base_path,folder,ff)
            if os.path.exists(eff):
                ff = eff
            
    if not(os.path.exists(ff)):
        if required:
            fatal("specified %s %s does not exist or is not readable" % (option, ff))
        else:
            ff = None

    return ff
    
def RequireMachine():

    if simenv.LocalMachine == None or simenv.LocalMachineEntry == None:
        fatal("Unknown local machine %s. Please use 'sim setup' to create a local machine entry from the generic template." % GetHostNameAlias())



def GetRsyncInfo(machineEntry):
    return GetSubbedInfo(machineEntry, ['rsynccmd', 'rsyncopts'])
    
def GetSSHInfo(machineEntry):
    return GetSubbedInfo(machineEntry, ['sshcmd', 'sshopts'])

def GetSubbedInfo(machineEntry, keys):

    DefineDatabase = simsubs.DefineDatabase()
    
    VerifyKeys(machineEntry, keys)
    
    data = list()

    for k in keys:
        d = DefineDatabase.SubAll(machineEntry.GetKey(k))
        data.append(d)
    
    return data

def RsyncVersion(rsyncInfo):
    (rsynccmd, rsyncopts) = rsyncInfo
    versionText, rc = ExecuteCommand('%s --version' % rsynccmd, True)
    # use search rather than match to skip possible text output by envsetup
    match = re.search('rsync +version +(\d+).(\d+).(\d+)', versionText)
    if not match or rc:
      fatal("Could not identify rsync version. rsync --version returned '%s'" % versionText.rstrip())
    else:
      (major, minor, revision) = match.groups()
    return (int(major), int(minor), int(revision))
