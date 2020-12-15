# simenv -- environment stuff now as a module, instead of a class

import os, sys
import libutil 

#simenv state
INITIALIZED = False

#module paths
BASE_PATH = None
LIB_PATH = None
ETC_PATH = None
LOG_PATH = None
MDB_PATH = None
MDB_BASE_PATH = None
SYNTAX_PATH = None
OPTIONS_PATH = None
CACTUS_PATH = None
CONFIGS_PATH = None
EXE_PATH = None
INTERNALDIRECTORY = "SIMFACTORY"

#executable and state information
EXECUTABLE = None
COMMAND = None
SETUP = False

#machine
LocalMachine = None
LocalMachineEntry = None

#Required singletons
OptionsManager = None
ConfigurationDatabase = None

#Verbose
VERBOSE = True
StoredVerbose = True
CommandDoc = dict()
CommandOptions = dict()

def popen(cmd):
    global_dict = globals()
    if global_dict['VERBOSE']:
        libutil.info('EXECUTING COMMAND: '+cmd)
    return os.popen(cmd)

def system(cmd):
    global_dict = globals()
    if global_dict['VERBOSE']:
        libutil.info('EXECUTING COMMAND: '+cmd)
    return os.system(cmd)

def init(base, callingExecutable, usageString, optionGroups):

    import simopts
    import simdb
    import simlib
    
    global_dict = globals()
    
    if global_dict['INITIALIZED']:
        return
    
    global_dict['BASE_PATH'] = base
    global_dict['LIB_PATH'] = "%s%s%s" % (global_dict['BASE_PATH'], os.sep, "lib")
    global_dict['ETC_PATH'] = "%s%s%s" % (global_dict['BASE_PATH'], os.sep, "etc")
    global_dict['LOG_PATH'] = "%s%s%s" % (global_dict['BASE_PATH'], os.sep, "../log")
    global_dict['MDB_PATH'] = "%s%s%s%s%s" % (global_dict['BASE_PATH'], os.sep, "mdb", os.sep, "machines")
    global_dict['MDB_BASE_PATH'] = "%s%s%s" % (global_dict['BASE_PATH'], os.sep, "mdb")
    global_dict['SYNTAX_PATH'] = "%s%s%s" % (global_dict['ETC_PATH'], os.sep, "syntax")
    global_dict['OPTIONS_PATH'] = "%s%s%s" % (global_dict['ETC_PATH'], os.sep, "options")
    global_dict['CACTUS_PATH'] = None
    global_dict['VERBOSE'] = True
    
    global_dict['VERSION'] = "???"
    
    global_dict['StoredVerbose'] = global_dict['VERBOSE']
    
    # a couple of paths that are in the cactus root.
    global_dict['CONFIGS_PATH'] = None
    global_dict['EXE_PATH'] = None
    
    global_dict['EXECUTABLE'] = callingExecutable
    global_dict['COMMAND'] = None
    global_dict['SETUP'] = False
    
    global_dict['INTERNALDIRECTORY'] = "SIMFACTORY"

    global_dict['VERBOSE'] = True
    
    global_dict['udb'] = None
    global_dict['cdb'] = None
    
    try:
        fd = os.popen("svnversion %s 2>&1" % global_dict['BASE_PATH'])
        ver = fd.read().strip()
        fd.close()
        
        #if ver.count(":") > 0:
        #   global_dict['VERSION'] = "r%s" % ver.split(":").pop(0)
        global_dict['VERSION'] = ver
    except:
        pass
        
    global_dict['OptionsManager'] = simopts.OptionManager(global_dict['OPTIONS_PATH'])

    if optionGroups != None:
        for og in optionGroups:
            global_dict['OptionsManager'].UseGroup(og)
    
    if usageString != None:
        global_dict['OptionsManager'].SetUsage(usageString)
    
    global_dict['OptionsManager'].Build()
    
    global_dict['StoredVerbose'] = global_dict['VERBOSE'] = global_dict['OptionsManager'].GetOption('verbose')
    
    # more initialization
    GetIniFiles()
    
    
    # init configuration database
    global_dict['ConfigurationDatabase'] = simdb.ConfigurationDatabase(mdbDirectory=global_dict['MDB_PATH'], udb=global_dict['udb'], cdb=global_dict['cdb'])
    global_dict['ConfigurationDatabase'].Load()

    global_dict['CACTUS_PATH'] = simlib.GetCactusPath()
    
    if global_dict['CACTUS_PATH'] == None:
        global_dict['CACTUS_PATH'] = global_dict['BASE_PATH']
        print "Warning: Unable to determine CACTUS_PATH, using %s instead" % global_dict['CACTUS_PATH']
        
    global_dict['CONFIGS_PATH'] = simlib.BuildPath([global_dict['CACTUS_PATH'], "configs"])
    if "CACTUS_CONFIGS_DIR" in os.environ:
        global_dict['CONFIGS_PATH'] = os.environ["CACTUS_CONFIGS_DIR"]
    global_dict['EXE_PATH'] = simlib.BuildPath([global_dict['CACTUS_PATH'], "exe"])
        
    # get our local machine, and our local machine entry.
    global_dict['LocalMachine'] = simlib.GetMachineName()

    if global_dict['LocalMachine'] != None:
        global_dict['LocalMachineEntry'] = global_dict['ConfigurationDatabase'].GetMachine(global_dict['LocalMachine'])
    else:
        global_dict['LocalMachineEntry'] = None

    cmd = str()
    
    count = 0
    for arg in sys.argv:
        if count == 0:
            cmd = arg
        else:
            cmd = "%s \"%s\"" % (cmd, arg)
        
        count = count + 1
        
    cmd = cmd.strip()
    
    if global_dict['VERBOSE']:
        print "Info: Simfactory command: %s" % cmd
        print "Info: Version %s" % global_dict['VERSION']
    
    global_dict['INITIALIZED'] = True
        
def GetIniFiles():

    global_dict = globals()

    ini_files = ['defs', 'defs.local']
    ini_name = ['cdb', 'udb']
    ini_req = ['defs']
    
    for i in range(len(ini_files)):
        ini = ini_files[i]
        attr = ini_name[i]
        
        path = "%s%s%s.ini" % (global_dict['ETC_PATH'], os.sep, ini)
        if os.path.exists(path):
            global_dict[attr] = path
        else:
            global_dict[attr] = None
            if ini in ini_req:
                print "No %s.ini found in %s, and --%s was not specified" % (ini, global_dict['ETC_PATH'], ini)
                sys.exit(1)

def ForceVerbose():
    global_dict = globals()
    
    global_dict['VERBOSE'] = True

def ResetVerbose():
    global_dict = globals()
    
    global_dict['VERBOSE'] = global_dict['StoredVerbose']
