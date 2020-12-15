import os, sys, re

from libutil import *
import simenv,simlib,simsubs

class RemoteEnvironment:
    def __init__(self):
        self.RemotePath = None
        pass
    
    def init(self, machineName):
        
        self.localMachineName = simenv.LocalMachine
        self.localMachineEntry = simenv.LocalMachineEntry

        if self.localMachineEntry == None:
            fatal("unknown local machine %s" % self.localMachineName)
        
        self.machineName = machineName
        
        if self.machineName == self.localMachineName:
            fatal("remote machine %s is the same as local machine %s" % (self.machineName, self.localMachineName))

        self.machineEntry = simenv.ConfigurationDatabase.GetMachine(self.machineName)

        if self.machineEntry == None:
            fatal("unknown remote machine %s" % self.machineName)
        
        req_keys = ['envsetup', 'sourcebasedir', 'user', 'hostname']
        simlib.VerifyKeys(self.localMachineEntry, req_keys)
                
        self.localEnvSetup = self.localMachineEntry.GetKey("envsetup")
        
        DefineDatabase = simsubs.DefineDatabase()
        DefineDatabase.Set("USER", simenv.LocalMachineEntry.user)

        self.localEnvSetup = DefineDatabase.SubAll(self.localEnvSetup)
        
    
    def GetRemoteExePath(self, local_path, remotes):
        exe = simenv.EXECUTABLE

        if simenv.OptionsManager.HasOption("remotecactuspath"):
            rcdir = simenv.OptionsManager.GetOption("remotecactuspath")
            
            warning("remote cactus path override, using %s" % rcdir)
            
            local_path = simenv.CACTUS_PATH
            remotes = rcdir
            self.RemotePath = rcdir
        else:
            self.RemotePath = simlib.ReplaceLeadingPath(simenv.CACTUS_PATH, simlib.GetLocalSourceBaseDir(), remotes)
        
        # if locals has a double //, trim off the first one.
        if local_path.startswith(r"//"):
            local_path = local_path[1:]

        if exe[0] == "/":
            exe = simlib.ReplaceLeadingPath(exe, local_path, remotes)
        else:
            path = simlib.GetDirSuffix(local_path)
            exe = simlib.BuildPath(remotes, path, exe)

        return exe
    
    def ExecuteSameCommand(self, parrotArguments=False, stripArguments=None, sshargs=''):
        
        localSourceBaseDir = simlib.GetLocalSourceBaseDir()
        sourceBaseDir = simlib.GetSourceBaseDir(self.machineEntry)
        
        remoteExe = self.GetRemoteExePath(localSourceBaseDir, sourceBaseDir)

        if simenv.COMMAND != None:
            remoteExe = "%s %s" % (remoteExe, simenv.COMMAND)
        
        self.ExecuteCommand(remoteExe, parrotArguments, stripArguments, sshargs)
        
    def ExecuteCommand(self, command, parrotArguments=False, stripArguments=None, sshargs=''):
    
        DefineDatabase = simsubs.DefineDatabase()
        
        cmdargs = ''
        if parrotArguments:
            for aa in simenv.OptionsManager.args:
                cmdargs = '%s %s' % (cmdargs, simlib.QuoteSafe(aa))
            argStr = simenv.OptionsManager.BuildOptionString(stripArguments)
            
            command = "%s %s %s" % (command, cmdargs, argStr)
            
            # lets append --machine <remotemachinename> to the arg string to avoid any confusion.
            command = "%s --machine=%s" % (command, self.machineName)
                
        cmd = simlib.GetSSHCommand(simenv.LocalMachine, self.machineName, command, sshargs, self.RemotePath)
            
        machineEntry = simenv.ConfigurationDatabase.GetMachine(self.machineName)
        localSourceBaseDir = simlib.GetLocalSourceBaseDir()
        source_name = simlib.GetDirSuffix(localSourceBaseDir)
        local_sourcedir = simlib.BuildPath(localSourceBaseDir, source_name)
        DefineDatabase.Set('SOURCEDIR', local_sourcedir)
        cmd = DefineDatabase.SubAll(cmd)
        DefineDatabase.Set('USER', machineEntry.user)
        localsshsetup = DefineDatabase.SubAll(machineEntry.localsshsetup)
        
        cmd = "{ %s; } && { %s; } && %s" % (self.localEnvSetup, localsshsetup, cmd)
        #dprint("Executing Command: %s" % cmd)
        ret = simlib.ExecuteCommand(cmd)
        (status, signal) = (ret >> 8, ret & 0xff)
        if status == 0xff:
            fatal("Could not execute command \"%s\"" % cmd)
