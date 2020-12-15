import shutil, os, os.path, sys
import restartlib
import libutil

from archive import *
from libutil import *
import simenv,simlib

class SimArchive:

    def __init__(self, driver, restart=None):
        self.Restart = restart
        
        self.DriverName = driver
        
        if self.Restart != None:
            self.SimulationName = self.Restart.SimulationName
        
        moduleName = "%s.Driver(self)" % self.DriverName
        
        dprint("loading driver: %s" % moduleName)
        
        self.Driver = eval(moduleName)
    
    def authenticate(self):
    
        machineEntry = simenv.LocalMachineEntry
        
        simlib.VerifyKeys(machineEntry, ['archiveuser'])
        
        username = machineEntry.GetKey('archiveuser')
        
        ret = self.Driver.authenticate(username)
        
        if ret != 0:
            dprint("Error: could not authenticate for archiving using driver %s" % self.DriverName)
            sys.exit(1)
        else:
            dprint("Authentication for archiveuser %s successful" % username)
            
    def store(self):
        
        if self.Restart == None:
            dprint("Error: Cannot store null simulation")
            sys.exit(1)
    
        machineEntry = simenv.LocalMachineEntry
        
        metadataFile = self.buildMetadataFile()

        simlib.VerifyKeys(machineEntry, ['archivebasepath'])
        dstPath = "%s/%s" % (machineEntry.GetKey('archivebasepath'), self.Restart.Properties.simulationid)
        
        if self.Restart.RestartID == -1:
            srcPath = self.Restart.SimulationDir
        else:
            srcPath = self.Restart.RestartDir
            dstPath = "%s/output-%s" % (dstPath, self.Restart.LongRestartID)
        
        dprint("Archive source path: %s" % srcPath)
        dprint("Archive destination path: %s" % dstPath)
        
        self.Driver.store(srcPath, dstPath)
    
    def get(self, simulationid, dstPath, restartid=None):
        
        machineEntry = simenv.LocalMachineEntry
        
        simlib.VerifyKeys(machineEntry, ['archivebasepath'])
        basepath = machineEntry.GetKey('archivebasepath')
        
        srcPath = "%s/%s" % (machineEntry.GetKey('archivebasepath'), simulationid)
        
        if restartid != None:
            if len(restartid) == 4:
                srcPath = "%s/output-%s" % (dstPath, restartid)
            else:
                srcPath = "%s/output-%04d" % (dstPath, int(restartid))
                
        dprint("Archive source path: %s" % srcPath)
        dprint("Archive destination path: %s" % dstPath)
        
        self.Driver.get(srcPath, dstPath)
        
    def listSimulations(self):
        machineEntry = simenv.LocalMachineEntry
        
        simlib.VerifyKeys(machineEntry, ['archivebasepath'])
        basepath = machineEntry.GetKey('archivebasepath')

        return self.GetStoredSimulations(basepath)
    
    def GetStoredSimulations(self, srcPath):
        
        files = self.Driver.list(srcPath)
        
        sims = list()
        
        for line in files:

            if not(line.startswith("simulation-")):
                continue
            
            sim = dict()
            
            parts = line.split("-")
            
            storedPath = "%s/%s" % (srcPath, line)
            
            sim['StoredPath'] = storedPath
            sim['SimulationId'] = line
            sim['SimulationName'] = parts[1]
            sim['Machine'] = parts[2]
            sim['Hostname'] = parts[3]
            sim['User'] = parts[4]
            sim['Date'] = "%s-%s" % (parts[5], parts[6])
            sim['Restarts'] = self.GetRestartsForSimulation(srcPath, sim['SimulationId'])
            
            sims.append(sim)
        
        return sims
    
    def GetRestartsForSimulation(self, srcPath, simulationName):
        simpath = simlib.BuildPath(srcPath, simulationName)
        
        files = self.Driver.list(simpath)
        
        restarts = list()
        
        for line in files:

            if not(line.startswith("output-")):
                continue
            
            id = line.split("-")[1]
            
            #dprint("found restart id: %s" % id)
            restarts.append(id)
        
        return restarts
        
    def buildMetadataFile(self):
    
        # for now, lets just use properties.ini from the restart
        
        return simlib.BuildPath(self.Restart.InternalDir, "properties.ini")
        
        
        