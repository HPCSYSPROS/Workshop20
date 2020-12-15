import os, sys
import libutil

from libutil import *
import simenv,simlib

class Driver:

    def __init__(self, simArchive):
        self.SimArchive = simArchive
        
        machineEntry = simenv.LocalMachineEntry
        
        simlib.VerifyKeys(machineEntry, ['archivehostname'])
        
        self.HostName = machineEntry.archivehostname
        
    def authenticate(self, username):
        
        authtool = "grid-proxy-info -timeleft"
        
        timeleft = int(simenv.popen(authtool).read().strip())
        
        dprint("time left: %s" % timeleft)
        
        if timeleft > 0:
            return 0
        
        authtool = "myproxy-logon"
        authcmd = "%s -l %s" % (authtool, username)
        
        dprint("Executing authcmd: %s\n" % authcmd)
        
        ret = simenv.system(authcmd)

        if ret != 0:
            dprint("Error: could not authenticate using myproxy-logon")
            return -1
        return 0
        
    def store(self, srcPath, dstPath):
    
        # make the dst folder.
        self.pmkdir(dstPath)
        
        cmd = "put -r %s/* %s" % (srcPath, dstPath)
        
        #dprint("rsynccmd: %s" % rsynccmd)
        self.cmd(cmd)
        
        dprint("Archiving complete.")
    
    def get(self, srcPath, dstPath):
    
        if not(os.path.exists(dstPath)):
            os.makedirs(dstPath)
        
        cmd = "get -r %s %s" % (srcPath, dstPath)
        
        self.cmd(cmd)
        
        dprint("Archive retrieval complete.")
    
    def pmkdir(self, dstPath):
        parts = dstPath.split("/")
        
        url = None
        for part in parts:
            if url == None:
                url = part
            else:
                url = "%s/%s" % (url, part)
            
            self.cmd("mkdir %s" % url)
    
    def lineToFile(self, line):
    
        parts = line.split(" ")
        return parts.pop()
    
    def list(self, srcPath):
        
        output = self.cmd("ls %s" % srcPath)
        
        files = list()
        
        for line in output.split("\n"):
            line = line.strip()
            
            if len(line) == 0:
                continue
            
            if not(line.startswith("d")):
                continue
            
            line = self.lineToFile(line)
            
            files.append(line)
        
        return files
        
    def cmd(self, cmd):
        
        fullcmd = "uberftp %s \"quote wait; parallel 8; %s\" 2>&1" % (self.HostName, cmd)
        dprint("Executing uberftp command: %s" % fullcmd)
        
        fd = simenv.popen(fullcmd)
        
        output = ""
        while True:
            line = fd.readline()
            if not line: break
            output = "%s%s" % (output, line)
        fd.close()
        
        return output
        
