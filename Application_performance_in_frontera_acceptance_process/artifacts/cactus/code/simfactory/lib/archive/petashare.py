import os, sys

from libutil import *
import simenv,simlib

class Driver:

    def __init__(self, simArchive):
        self.SimArchive = simArchive
        machineEntry = simenv.LocalMachineEntry
        
        simlib.VerifyKeys(machineEntry, ['archivetoolspath'])
        
        self.ToolsPath = machineEntry.archivetoolspath
        
        if not(os.path.exists(self.ToolsPath)):
            dprint("Error: could not access archive tools path %s" % self.ToolsPath)
            sys.exit(1)
            
    def authenticate(self, username):
        
        authtool = "pchangeuser"
        authcmd = "%s %s" % (authtool, username)
        fullcmd = "%s/%s 2>&1" % (self.ToolsPath, authcmd)
        
        dprint("Executing authcmd: %s\n" % fullcmd)
        
        dprint("Enter Petashare Password: ",)
        sys.stdout.flush()

        fd = simenv.popen(fullcmd)
        
        output = ""
        while True:
            line = fd.readline()
            if not line: break
            output = "%s%s" % (output, line)
        
        fd.close()
        
        dprint("\n")
        sys.stdout.flush()
        
        if output.count("rcAuthResponse failed") > 0:
            lines = output.split("\n")
            
            for line in lines:
                if line.count("rcAuthResponse failed") > 0:
                    dprint("Error: could not authenticate to petashare: %s" % line)
            return -1
        
        return 0
        
    def store(self, srcPath, dstPath):
    
        # make the dst folder.
        self.cmd("pmkdir -p %s" % dstPath)
        
        rsynccmd = "prsync -rV %s i:%s" % (srcPath, dstPath)
        
        #dprint("rsynccmd: %s" % rsynccmd)
        self.cmd(rsynccmd)
        
        dprint("Archiving complete.")
    
    def get(self, srcPath, dstPath):
    
        if not(os.path.exists(dstPath)):
            os.makedirs(dstPath)
        
        rsynccmd = "prsync -rV i:%s %s" % (srcPath, dstPath)
        
        #dprint("rsynccmd: %s" % rsynccmd)
        self.cmd(rsynccmd)
        
        dprint("Archive retrieval complete.")
    
    
    def lineToFile(self, line):

        line = line.replace("C- ", "")
        
        storedPath = line
        
        removePart = "%s/" % srcPath
        line = line.replace(removePart, "")
        
        return line
        
    def list(self, srcPath):
        output = self.cmd("pls %s" % srcPath)
        
        files = list()
        
        for line in output.split("\n"):
            line = line.strip()
            
            if len(line) == 0:
                continue
            
            if not(line.startswith("C-")):
                continue

            line = self.lineToFile(line)
            
            files.append(line)
        
        return files
    
    def cmd(self, cmd):
        
        fullcmd = "%s/%s 2>&1" % (self.ToolsPath, cmd)
        dprint("Executing petashare command: %s" % fullcmd)
        
        fd = simenv.popen(fullcmd)
        
        output = ""
        while True:
            line = fd.readline()
            if not line: break
            output = "%s%s" % (output, line)
        fd.close()
        
        return output
        
