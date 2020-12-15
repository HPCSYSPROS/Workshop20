# 
# MDB -- MachineDatabase wrapper around pyini.
# Convienient access to the machine database.
# 
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
# 
#

import pyini
import os
import sys

from libutil import *
import simenv,simlib

class Entry:
    def __init__(self, name, dict):
        keys = dict.keys()
        
        self.InternalDictionary = dict
        self.EntryName = name
        
        for key in keys:
            io = dict[key]
            setattr(self, key, io.Value)
    
    def GetKeys(self):
        return self.InternalDictionary.keys()
    
    def HasKey(self, key):
        return hasattr(self, key)
    
    # can also just be accessed as MachineEntry.attribute
    def GetKey(self, key):
        if self.HasKey(key):
            return getattr(self, key)
        
        return None
        
    def HasKeys(self, keys):
        for k in keys:
            if not(self.HasKey(k)):
                return False
        
        return True
        
class ConfigurationDatabase:

    def __init__(self, mdbDirectory=None, cdb=None, udb=None):
        self.mdbDirectory = mdbDirectory
        self.mdbFilename = None
        self.cdbFilename = cdb
        self.udbFilename = udb
        
        self.SyntaxFile = "%s%s%s" % (simenv.SYNTAX_PATH, os.sep, "mdb-syntax.ini")
        self.DescriptionSyntaxFile = "%s%s%s" % (simenv.SYNTAX_PATH, os.sep, "description-syntax.ini")
        self.MachineParser = None
        self.ConfigParser = None
        self.syntaxChecker = None
        self.MachineCache = dict()
        self.ConfigCache = dict()
        
    def Load(self, mdbDirectory=None, cdb=None, udb=None):
        
        if mdbDirectory !=None:
            self.mdbDirectory = mdbDirectory
        
        if cdb != None:
            self.cdbFilename = cdb
            
        if udb != None:
            self.udbFilename = udb
            
        if (mdbDirectory == None and self.mdbDirectory == None) or (not(os.path.exists(self.mdbDirectory))):
            fatal("initializing machine database: No database file provided, or is not readable")
        
        if cdb == None and self.cdbFilename == None:
            fatal("cannot initialize configuration database, no database file provided")
        
        self.CheckSyntax()
        # if we made it here, our syntax has passed and self.parser will be available.

    def CheckSyntax(self):
        
        # first, verify the correctness of our mdb syntax file.
        sc = pyini.IniSyntaxChecker(self.DescriptionSyntaxFile, self.SyntaxFile)
        sc.SyntaxCheck()
        
        # now verify the correctness of our passed in mdb database
        self.MachineParser = pyini.IniParser()
        
        for ff in os.listdir(self.mdbDirectory):
            
            if not(ff.endswith(".ini")):
                continue
                
            filePath = simlib.BuildPath(simenv.MDB_PATH, ff)

            self.MachineParser.UpdateFromIni(filePath, True)
        
        # load the cdb database, which at the moment has no syntax file and convert any blocks to lists
        self.ConfigParser = pyini.IniParser(self.cdbFilename)
        self.ConvertBlocks()
        
        if self.udbFilename != None:
            # import new sections to the machine database, but not the config database.
            self.MachineParser.UpdateFromIni(self.udbFilename, True)
            self.ConfigParser.UpdateFromIni(self.udbFilename)

            syntaxChecker = pyini.IniSyntaxChecker(self.SyntaxFile, self.MachineParser, True)
            syntaxChecker.SyntaxCheck()
    
        self.MachineParser.UpdateFromDict(simenv.OptionsManager.MDBKeys)
    
    def ConvertBlocks(self):
        keys = self.ConfigParser.GetGlobalKeys()
        
        for k in keys:
            io = self.ConfigParser.GetGlobalOption(k)
            if io.IsBlock:
                io.ConvertToList()

    # --- ACCESSOR METHODS ---
    
    def GetMachine(self, key):
        tt = "Machine"
        
        cache = getattr(self, "%sCache" % tt)
        parser = getattr(self, "%sParser" % tt)
        
        if cache.has_key(key):
            return cache[key]
        
        if parser.HasSection(key):
            sdict = parser.GetSectionAsDict(key)
            entry = Entry(key, sdict)
            cache[key] = entry
            return entry
        else:
            fatal("retrieving %s entry for %s %s: %s %s doesn't exist" % (tt, tt, key, tt, key))
    
    def HasMachineEntry(self, key):
        tt = "Machine"
        cache = getattr(self, "%sCache" % tt)
        parser = getattr(self, "%sParser" % tt)
        
        if cache.has_key(key):
            return True
        
        return parser.HasSection(key)

    def GetMachines(self):
        tt = "Machine"
        parser = getattr(self, "%sParser" % tt)
        return parser.GetSections()
        
    def GetConfigOption(self, option):
        
        ret = self.ConfigParser.GetGlobalOption(option)
        
        if ret != None:
            return self.ConfigParser.GetGlobalOption(option).Value
        else:
            return None
            
    def HasConfigOption(self, option):
        return self.ConfigParser.HasOption(None, option)

    def GetConfigOptions(self):
        return self.ConfigParser.GetGlobalKeys()
