# sim-build -- build the cactus distribution
#
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
#

import sys, os, shutil

import restartlib,simenv,simlib,simsubs, simproperties

from libutil import *

## FUNCTIONS ##

known_commands = ['build']
usage_strings = {'build':'build cactus either locally or remotely'}
                
def PrepConfiguration(configName):
    
    DefineDatabase = simsubs.DefineDatabase()

    machineEntry = simenv.LocalMachineEntry

    req_keys = ['user', 'email', 'make']
    simlib.VerifyKeys(machineEntry, req_keys)
    
    user = machineEntry.GetKey("user")
    email = machineEntry.GetKey("email")
    
    DefineDatabase.Set('USER', user)
    DefineDatabase.Set('EMAIL', email)

    make = DefineDatabase.SubAll(machineEntry.GetKey('make'))

    subfolders = ['bindings', 'build', 'config-data', 'lib', 'scratch']
    
    configBase = simlib.BuildPath(simenv.CONFIGS_PATH, configName)
    
    propertiesFile = simlib.BuildPath(configBase, "properties.ini")

    if not os.path.exists(configBase):
        
        try:
            os.mkdir(configBase)
        except OSError:
            pass
        
        #simlib.ExecuteCommand("mkdir %s" % configBase)
        for sub in subfolders:
            sf = simlib.BuildPath(configBase, sub)
            simlib.ExecuteCommand("mkdir %s" % sf)
    else: # The config dir exists already
        if not os.path.exists(propertiesFile):
            fatal("Configuration %s has no properties.ini file.  Either this is not a SimFactory configuration, or the configuration has become corrupted.  If you think the configuration has been corrupted, you will need to delete it and run the sim build command again." % (configName))

    prop = simproperties.SimProperties()
    prop.init(propertiesFile)
    
    build_reconfig = GetConfigValue(prop, "reconfig")
    build_clean    = GetConfigValue(prop, "clean")

    build_debug    = GetConfigValue(prop, "debug")
    build_optimise = GetConfigValue(prop, "optimise")
    build_unsafe   = GetConfigValue(prop, "unsafe")
    build_profile  = GetConfigValue(prop, "profile")
    
    if build_debug:
        display("Disabling optimisation because debug was selected")
        build_optimise = False

    # update our properties.ini
    
    prop.RemoveProperty("reconfig")
    prop.RemoveProperty("clean")

    prop.AddProperty("debug",    build_debug   )
    prop.AddProperty("optimise", build_optimise)
    prop.AddProperty("unsafe",   build_unsafe  )
    prop.AddProperty("profile",  build_profile )
    
        
    # Remove the executable, so that it is not accidentally used
    # while it is rebuilt. Do this before the configuration is
    # modified in any way, but only after some basic error
    # checking.
    try:
        os.remove("exe/cactus_%s" % configName)
    except OSError:
        pass
    try:
        shutil.rmtree("exe/%s" % configName, ignore_errors=True)
    except OSError:
        pass

    prop.Save()
        
    storedOptions = simlib.BuildPath(configBase, "OptionList")
    
    hasStoredOptions = os.path.exists(storedOptions)
    
    cctkConfigPath = simlib.BuildPath(configBase, "config-data", "cctk_Config.h")
    
    hasCompleteConfig = os.path.exists(cctkConfigPath)
    
    if hasStoredOptions:
        storedOptionSettings = simlib.GetFileContents(storedOptions)
    else:
        storedOptionSettings = None
    
    info("HasStoredOptions: %s" % hasStoredOptions)
    
    removeConfig = False   # default behaviour

    # we are only updating the OptionList if its explicitly specified with --optionlist, really.
    if hasStoredOptions and not simenv.OptionsManager.RawOptionDefined("optionlist"):
        info("Use stored options list: %s" % storedOptions)
        optionlist = storedOptions
    else:
        optionlist = simlib.GetOptionList(False)

    info("optionlist is: %s" % optionlist)
    
    if optionlist == None or not os.path.exists(optionlist):
        display("Warning: no option list specified, using blank option list")
        optionSettings = str()
    else:
        optionSettings = simlib.GetFileContents(optionlist)
    
    DefineDatabase.AddReplacement("DEBUG", simlib.BoolToYesNo(build_debug))
    DefineDatabase.AddReplacement("OPTIMISE", simlib.BoolToYesNo(build_optimise))
    DefineDatabase.AddReplacement("UNSAFE", simlib.BoolToYesNo(build_unsafe))
    DefineDatabase.AddReplacement("PROFILE", simlib.BoolToYesNo(build_profile))
    
    
    # add support for CACHELINE_BYTES and CACHE_SIZE
    CACHELINE_BYTES = simenv.LocalMachineEntry.GetKey("L3linesize")
    
    if CACHELINE_BYTES != None:
        DefineDatabase.AddReplacement("CACHELINE_BYTES", CACHELINE_BYTES)
    
    CACHE_SIZE = simenv.LocalMachineEntry.GetKey("L3size")

    if CACHE_SIZE != None:
        DefineDatabase.AddReplacement("CACHE_SIZE", CACHE_SIZE)
    
    optionSettings = DefineDatabase.SubAll(optionSettings)

    hasOutdatedConfig = (not hasCompleteConfig) or (storedOptionSettings == None) or (storedOptionSettings != optionSettings)
    
    if hasOutdatedConfig:
        info("hasCompleteConfig: %s" % hasCompleteConfig)

    info("hasOutdatedConfig: %s" % hasOutdatedConfig)
    info("build_reconfig: %s" % build_reconfig)
    
    # if build virtual, disable configging.
        
    if hasOutdatedConfig or build_reconfig:
    
        if storedOptionSettings != None:
            oldVersion = simlib.GetVersion(storedOptionSettings)
        else:
            oldVersion = ""
        
        newVersion = simlib.GetVersion(optionSettings)
        info("oldVersion %s, newVersion %s" % (oldVersion, newVersion))
        
        removeConfig = removeConfig or newVersion != oldVersion

        display("Reconfiguring %s" % configName)
        simlib.RemoveFile("%s.old" % storedOptions)
        simlib.RenameFile(storedOptions, "%s.old" % storedOptions)
        
        # TODO: Write new option list only after the make *-realclean below!
        display("Writing configuration to: %s" % storedOptions)
        
        simlib.WriteContents(storedOptions, optionSettings)
        
        #"echo yes | { $make $configuration_name-config options=configs/$configuration_name/OptionList; }";
        
        if not simenv.OptionsManager.GetOption('virtual'):
            cmd = "cd %s && echo yes | { %s %s-config options=%s; } 2>&1" % (simenv.CACTUS_PATH, make, configName, storedOptions)
            ret = simlib.ExecuteCommand(cmd)
            if ret != 0:
                sys.exit(1)
                
            #force a rebuild
            simlib.RemoveFile(simlib.BuildPath(configBase, 'config-data', 'make.thornlist'))
        
            #remove the old config if necessary
            if removeConfig:
                display("Complete rebuild required")
    else:
        if not hasStoredOptions:
            fatal("Configuration %s has no option list, and no option list was specified." % configName)
    
    if not simenv.OptionsManager.GetOption('virtual'):
        info("build_clean: %s" % build_clean)
        info("removeConfig: %s" % removeConfig)
        
        if build_clean or removeConfig:
            display("Cleaning %s" % configName)
            ret = simlib.ExecuteCommand("cd %s && %s %s-realclean 2>&1" % (simenv.CACTUS_PATH, make, configName))
            if ret != 0:
                sys.exit(1)
    
    ## deal with submit script now
        
    storedSubmitScript = simlib.BuildPath(configBase, "SubmitScript")
    hasStoredSubmitScript = os.path.exists(storedSubmitScript)
    
    if simenv.OptionsManager.GetOption('no-submitscript'):
        if hasStoredSubmitScript:
            display("Removing stored submit script for configuration %s" % configName)
            os.unlink(storedSubmitScript)
            warning("empty submit script will disable submission")
    else:
        # the submit script is entirely optional. no submit script
        # just means submission is disabled.
    
        submitscript = simlib.GetSubmitScript()
        info("SubmitScript is: %s" % submitscript)
    
        if submitscript != None:
            sContents = simlib.GetFileContents(submitscript)
            ssContents = simlib.GetFileContents(storedSubmitScript)
        
            submitScriptOutdated = not hasStoredSubmitScript or simenv.OptionsManager.RawOptionDefined("submitscript")
            
            if sContents != ssContents:
                warning("default submit script contents have changed")
            
            if (submitScriptOutdated or build_reconfig):
                display("Updated script file for configuration %s" % configName)
                simlib.RemoveFile("%s.old" % storedSubmitScript)
                simlib.RenameFile(storedSubmitScript, "%s.old" % storedSubmitScript)
                shutil.copy(submitscript, storedSubmitScript)
            else:
                if not hasStoredSubmitScript:
                    warning("empty submit script will disable submission")
    
    ## deal with run script now
    
    storedRunScript = simlib.BuildPath(configBase, "RunScript")
    hasStoredRunScript = os.path.exists(storedRunScript)
    
    runscript = simlib.GetRunScript(False)
    info("RunScript is: %s" % runscript)
    
    if runscript != None:
        sContents = simlib.GetFileContents(runscript)
        ssContents = simlib.GetFileContents(storedRunScript)
        
        runScriptOutdated = not hasStoredRunScript or simenv.OptionsManager.RawOptionDefined("runscript")
        
        if sContents != ssContents:
            warning("default run script contents have changed")
                
        if (runScriptOutdated or build_reconfig):
            display("Updated runscript file for configuration %s" % configName)
            simlib.RemoveFile("%s.old" % storedRunScript)
            simlib.RenameFile(storedRunScript, "%s.old" % storedRunScript)
            shutil.copy(runscript, storedRunScript)
        else:
            if not hasStoredRunScript:
                fatal("Configuration %s has no run script, and no run script file was specified" % configName)

    ## deal with thorn list now 
    storedThornList = simlib.BuildPath(configBase, "ThornList")
    hasStoredThornList = os.path.exists(storedThornList)
    
    needThornList = not hasStoredThornList
    thornlist = simlib.GetThornList(needThornList)
    info("ThornList is: %s" % thornlist)
        
    if thornlist != None:
        
        tContents = simlib.GetThornListContents(thornlist)
        ttContents = simlib.GetFileContents(storedThornList)
                
        thornListOutdated = not hasStoredThornList or simenv.OptionsManager.RawOptionDefined("thornlist")
        
        if tContents != ttContents:
            warning("default thorn list contents have changed")

        if (thornListOutdated or build_reconfig):
            display("Updated thorn list for configuration %s" % configName)
            simlib.RemoveFile("%s.old" % storedThornList)
            simlib.RenameFile(storedThornList, "%s.old" % storedThornList)
            simlib.WriteContents(storedThornList, tContents)
    else:
        if not hasStoredThornList:
            fatal("Configuration %s has no thorn list, and no thorn list file was specified" % configName)

   
def GetConfigValue(prop, option):

    # This seems confused: Either the properties or the options
    # should take precedence.
    if simenv.OptionsManager.IsNegated(option):
        return False

    if simenv.OptionsManager.RawOptionDefined(option):
        return True

    if prop.HasProperty(option):
        return CoerceBool(prop.GetProperty(option))

    return simenv.OptionsManager.GetOption(option)

def BuildConfiguration(configName):
    
    DefineDatabase = simsubs.DefineDatabase()

    machineEntry = simenv.LocalMachineEntry

    req_keys = ['user', 'email', 'make']
    simlib.VerifyKeys(machineEntry, req_keys)
    
    user = machineEntry.GetKey("user")
    email = machineEntry.GetKey("email")
    
    DefineDatabase.Set('USER', user)
    DefineDatabase.Set('EMAIL', email)
    
    make = DefineDatabase.SubAll(machineEntry.GetKey('make'))
    
    display("Building %s" % configName)
    ret = simlib.ExecuteCommand("cd %s && %s %s 2>&1" % (simenv.CACTUS_PATH, make, configName))
    if ret != 0:
        sys.exit(1)
    
    display("Building utilities for %s" % configName)
    ret = simlib.ExecuteCommand("cd %s && %s %s-utils 2>&1" % (simenv.CACTUS_PATH, make, configName))
    if ret != 0:
        sys.exit(1)
    
    # display("Building documentation for %s" % configName)
    # ret = simlib.ExecuteCommand("cd %s && %s %s-ThornGuide 2>&1" % (simenv.CACTUS_PATH, make, configName))
    # if ret != 0:
    #     sys.exit(1)
        
def BuildVirtual(configName):

    display("Building virtual config %s" % configName)
    
    if not(simenv.OptionsManager.HasOption("virtual-executable")):
        fatal("virtual executable is not set")
    
    vexe = simenv.OptionsManager.GetOption("virtual-executable")
    
    if not(os.path.exists(vexe)):
        fatal("virtual executable %s does not exist or is not readable" % vexe)
    
    newvpath = simlib.BuildPath(simenv.EXE_PATH, "cactus_%s" % configName)

    display("Copying virtual executable from %s to %s" % (vexe, newvpath))

    simlib.MakeDirs(simenv.EXE_PATH)
    shutil.copy(vexe, newvpath)

def CommandDispatch():
    global known_commands
    
    if simenv.COMMAND == None:
        command = simenv.OptionsManager.args.pop(0)
    else:
        command = simenv.COMMAND
        if command == simenv.OptionsManager.args[0]:
            simenv.OptionsManager.args.pop(0)
    
    if command not in known_commands:
        simenv.OptionsManager.args.insert(0, command)
        command_build()
        return
        
    info("Executing command: %s" % command)
    exec("command_%s()" % command)
    
############################### MAIN ###############################

simenv.CommandDoc['build'] = """

Usage
-----

sim build [options] [ *config* ] [options]

  Build a Cactus configuration called *config*.  If *config* is
  omitted, the default configuration name specified in defs.local.ini
  is used.  

"""

simenv.CommandOptions['build'] = "sim-build.ini"

def command_build():
    
    DefineDatabase = simsubs.DefineDatabase()
    
    configList = simenv.OptionsManager.args
    
    if len(configList) == 0:
        configList = [simlib.GetDefaultConfiguration()]
    
    for config in configList:
        display("Using configuration: %s" % config)
    
    machineEntry = simenv.LocalMachineEntry
    
    req_keys = ['user', 'email']
    simlib.VerifyKeys(machineEntry, req_keys)
    
    user = machineEntry.GetKey("user")
    email = machineEntry.GetKey("email")
    
    DefineDatabase.Set('USER', user)
    DefineDatabase.Set('EMAIL', email)
        
    cactusDir = simenv.CACTUS_PATH
    info("Cactus Directory: %s" % cactusDir)
    
    for config in configList:
        PrepConfiguration(config)

        if not simenv.OptionsManager.GetOption('virtual'):
            BuildConfiguration(config)
        else:
            BuildVirtual(config)

    display("Done.")
            
def main():
    simlib.RequireMachine()

    ############################################

    # make sure our configs directory is created
    try:
        os.mkdir(simenv.CONFIGS_PATH)
    except OSError:
        pass
    
    if simenv.CACTUS_PATH == None:
        fatal("cannot proceed with an unknown CACTUS_PATH")

    # get our modes of compilation

    if len(simenv.OptionsManager.args) > 0:
        if simenv.OptionsManager.args[0] == 'build':
            simenv.OptionsManager.args.pop(0)
    else:
        simenv.OptionsManager.args.append('build')
        
    CommandDispatch()
