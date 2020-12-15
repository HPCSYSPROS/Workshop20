# sim-sync -- sync source trees across machines using information from the machine database
#
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
#

import sys, os, re

optionGroups = ['sim-sync']
Purpose = "sync cactus sourcetree to a remote machine"

############################################
import simenv,simsubs,simlib
from libutil import *

known_commands = []
usage_strings = {'sync': 'sync cactus sourcetree to a remote machine'}

global local_suffix
global local_sourcebasedir
global localMachineName
global localMachineEntry

local_suffix = None
local_sourcebasedir = None
localMachineName = None
localMachineEntry = None

############################### FUNCTIONS ###############################

def CompileCommand(machineName, rsyncInfo, paths, oldRsync=False):
    
    global localMachineName
    global localMachineEntry
    global local_sourcebasedir
    global local_suffix
    
    DefineDatabase = simsubs.DefineDatabase()

    (rsynccmd, rsyncopts) = rsyncInfo

    rsyncversion = simlib.RsyncVersion(rsyncInfo)

    if machineName == localMachineName:
        fatal("cannot sync to local machine")

    machineEntry = simenv.ConfigurationDatabase.GetMachine(machineName)
    
    # get our IO machine -- if iomachine doesn't exist, it returns itself.
    
    ioMachine = simlib.GetIOMachine(machineName)
    
    if ioMachine != machineName:
        machineEntry = simenv.ConfigurationDatabase.GetMachine(ioMachine)
        machineName = ioMachine
    
    # make sure we have all the info we need.
    mreq_keys = ['user', 'sourcebasedir', 'hostname', 'rsynccmd', 'rsyncopts', 'sshcmd', 'sshopts']
    simlib.VerifyKeys(machineEntry, mreq_keys)

    DefineDatabase.Set('USER', machineEntry.user)
    
    sourcebasedir = DefineDatabase.SubAll(machineEntry.sourcebasedir)
    
    # There need to be two define databases, one for the local and
    # another for the remote system. Currently, USER is for the
    # remote and SOURCEDIR is for the local system -- this is
    # inconsistent.

    #local_sourcebasedir = DefineDatabase.SubAll(localMachineEntry.sourcebasedir)

    source_name = local_suffix
    
    #source_name = simlib.GetDirSuffix(local_sourcebasedir)
    local_sourcedir = simlib.BuildPath(local_sourcebasedir, source_name)
    DefineDatabase.Set('SOURCEDIR', local_sourcedir)
    localsshsetup = DefineDatabase.SubAll(machineEntry.localsshsetup)
    
#    sshcmd  = machineEntry.sshcmd
#    sshopts = machineEntry.sshopts
#    sshcmd = "%s %s" % (sshcmd, sshopts)
#    trampoline = simlib.GetTrampoline(ioMachine)
#    if trampoline:
#        tentry = simenv.ConfigurationDatabase.GetMachine(trampoline)
#        trampoline_sourcebasedir = simlib.GetSourceBaseDir(tentry)
#        trampoline_sourcedir = os.path.join(trampoline_sourcebasedir, source_name)
#        DefineDatabase.Set('SOURCEDIR', trampoline_sourcedir)
#        sshcmd = DefineDatabase.SubAll(sshcmd)
#        DefineDatabase.Set('SOURCEDIR', local_sourcedir)
#    sshcmd = simlib.GetSSHCommand(trampoline, sshcmd)
#    sshcmd = DefineDatabase.SubAll(sshcmd)
    
    sshcmd = simlib.GetSSHCommand(localMachineName, ioMachine, None)
    
    # If there is more than one explicit path given only accept
    # top-level directories This could be extended to do the right
    # thing when multiple explicit paths are given but would require
    # running rsync multiple times
    rsyncfiles = []
    if len(paths) == 1:
      rsyncfiles = paths
      [head, tail] = os.path.split(rsyncfiles[0])
      local_suffix = os.path.join(local_suffix, head)
    else:
      for file in paths:
        if file == '':
          continue
        [head, tail] = os.path.split(file)
        if head != "":
          fatal("Only top level paths may be specified when syncing multiple paths but '%s' given." % file)
        elif not os.path.exists(tail):
          warning("Specified sync path '%s' not found." % tail)
        else:
          rsyncfiles.append(tail)

    # Note:
    # - We omit --times, since the target system should recompile when
    #   a file changed there
    # - We omit --group and --owner, since this probably works for
    #   root only, and all files on the targe system should rather
    #   belong to the user
    # - We omit -D, since we don't expect to have devices or other
    #   special files
    # - Since we omit --times, we add --checksum so that differing
    #   modification times don't lead to a retransmission
    # - NOTE: rsync before 3.0.8 has an error that makes --checksum
    #   unreliable under certain circumstances (see
    #   <https://rsync.samba.org/ftp/rsync/src/rsync-3.0.8-NEWS>). We
    #   cannot use --ignore-times instead because this modifies all
    #   files' timestamps (and is also very expensive).
    rsyncoptions = [
        '--checksum',
        '--compress',
        '--delete',
        '--hard-links',
        '--links',
        '--partial',
        '--perms',
        '--progress',
        '--recursive',
        '--sparse',
        '--stats',
        #'--times',
        '--verbose']

    # We need to use a different rules file for pre 3.0.0 versions of rsync
    ruleFile = 'filter.rules'
    if rsyncversion[0] < 3 or oldRsync:
        info('Old rsync version in use (local version is %s.%s.%s). Upgrade to at least 3.0.0 both locally and remotely for best performance.' % rsyncversion)
        ruleFile = 'filter.prersync3.rules'

    userRuleFile = "%s/etc/%s" % (simenv.BASE_PATH, "filter.local.rules")
    if os.path.exists(userRuleFile):
        rsyncoptions.append("--filter 'merge %s'" % (userRuleFile))

    rsyncoptions.append("--filter 'merge %s/etc/%s'" %
                        (simenv.BASE_PATH, ruleFile))

    fullpath = "%s@%s:%s%s%s" % (machineEntry.user, machineEntry.hostname, sourcebasedir, os.sep, local_suffix)
    
    arguments = " ".join(simlib.GetArguments())
    
    cmd = "%s --rsh=%s --rsync-path=%s %s %s %s" % (rsynccmd, simlib.QuoteSafe(sshcmd), simlib.QuoteSafe(machineEntry.rsynccmd), rsyncopts, machineEntry.rsyncopts, arguments)
    cmd = "%s %s" % (cmd, " ".join(rsyncoptions))
    cmd = "%s %s" % (cmd, " ".join(rsyncfiles))
    cmd = "%s %s" % (cmd, fullpath)
        
    mkdirpath = os.path.join(sourcebasedir, local_suffix)
    mkdircmd = "mkdir -p %s" % simlib.QuoteSafe(mkdirpath)
    #mkdircmd = "%s %s@%s %s" % (sshcmd, machineEntry.user, machineEntry.hostname, simlib.QuoteSafe(mkdircmd))
    mkdircmd = simlib.GetSSHCommand(localMachineName, ioMachine, mkdircmd)

    cmd = "cd %s && { %s; } && { %s || { %s && %s; } }" % (simlib.QuoteSafe(local_sourcedir), localsshsetup, cmd, mkdircmd, cmd)
    
    return cmd
    
############################### MAIN ###############################

simenv.CommandDoc['sync'] = """

Usage
-----

sim sync [options] [ *machine* ]

  Synchronise local Cactus source tree to *machine*.  The first time
  this operation is performed, it can take some time as all the source
  is copied.  Subsequent calls will only synchronise the files which
  have changed.

  You can control exactly what is synced by creating a file
  simfactory/etc/filter.local.rules which is used as an Rsync filter
  rules file (see "man rsync" for a description).
"""

simenv.CommandOptions['sync'] = "sim-sync.ini"

def main():
    ## HEADER ##
    
    simlib.RequireMachine()
    
    global local_suffix
    global local_sourcebasedir
    global localMachineName
    global localMachineEntry

    machineList = []
    pathList = []

    # Include sync-parfiles in paths if --sync-parfiles option is given
    if simenv.OptionsManager.GetOption('sync-parfiles'):
        pathList.extend(simenv.ConfigurationDatabase.GetConfigOption('sync-parfiles'))

    # Include sync-sourcetree in paths if --sync-sourcetree option is given
    if simenv.OptionsManager.GetOption('sync-sourcetree'):
        pathList.extend(simenv.ConfigurationDatabase.GetConfigOption('sync-sources'))

    if simenv.OptionsManager.HasOption('sync-path'):
        pathList.extend(simenv.OptionsManager.GetOption('sync-path'))

    # Build up the list of machines and paths
    addPath = False
    for arg in simenv.OptionsManager.args:
        if arg == 'sync':
            continue
        machineList.append(arg)

    # If no paths were specified, we use everything in sync-parfiles and sync-sources
    if len(pathList) == 0:
        pathList.extend(simenv.ConfigurationDatabase.GetConfigOption('sync-parfiles'))
        pathList.extend(simenv.ConfigurationDatabase.GetConfigOption('sync-sources'))

    if len(machineList) == 0:
        fatal("Error: no machines specified")

    info("using machines: %s" % machineList)
    info("current working directory: %s" % os.getcwd())
    info("synchronising paths: %s" % pathList)

    if simenv.CACTUS_PATH == None:
        fatal("cannot proceed with an unknown CACTUS_PATH")
        
    cactusDir = simenv.CACTUS_PATH
    info("Cactus Directory: %s" % cactusDir)
    
    # get our local_sourcebasedir
    
    local_sourcebasedir = simlib.GetLocalSourceBaseDir()
    local_suffix = simlib.ReplaceLeadingPath(simenv.CACTUS_PATH, local_sourcebasedir, "")
    info("CACTUS_PATH: %s" % simenv.CACTUS_PATH)
    info("local_sourcebasedir: %s" % local_sourcebasedir)
    info("local_suffix: %s" % local_suffix)

    localMachineName = simenv.LocalMachine
    localMachineEntry = simenv.LocalMachineEntry

    # returns list (rsynccmd, rsyncopts)
    rsyncInfo = simlib.GetRsyncInfo(localMachineEntry)

    errors = 0
    for m in machineList:
        cmd = CompileCommand(m, rsyncInfo, pathList)
        ret = simlib.ExecuteCommand(cmd)

        # If we get a protocol error, try again with less modern rsync options
        if ret == 512:
            cmd = CompileCommand(m, rsyncInfo, pathList, True)
            ret = simlib.ExecuteCommand(cmd)

        if ret != 0:
            errors = errors + 1
            display("\nError while syncing to %s." % m)

    if ret == 0:
        info("\nSync complete.")
    else:
        display("\nErrors occurred during synchronisation.")
