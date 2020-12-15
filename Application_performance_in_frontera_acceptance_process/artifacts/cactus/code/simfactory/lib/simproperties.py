import libutil
import simlib
import simenv
import pyini
import os, sys

from libutil import *

class SimProperties:
    def __init__(self, baseProperties=None):

        self.IniSection = 'properties'
        self.BaseProperties = baseProperties
    
    def InitBlank(self):
        self.parser = pyini.IniParser()

    def init(self, filename):
        importExisting = True

        self.Filename = filename
        
        if not(os.path.exists(filename)):
            importExisting = False
            
        self.parser = pyini.IniParser(self.Filename)

        if importExisting:
            ret = self.ImportProperties()
            
            if ret < 0:
                warning("Importing properties from %s failed." % self.Filename)
            
        if self.BaseProperties != None:
            for key in self.BaseProperties.keys():
                self.AddProperty(key, self.BaseProperties[key])
    
    def ImportProperties(self):
    
        section = self.parser.GetSectionAsDict(self.IniSection)
        
        if section == None:
            return -1

        for key in section.keys():
            io = section[key]
            
            if io.IsBlock:
                io.ConvertToList()
            
            value = io.Value
            
            setattr(self, key, value)
            
        return 0
        
    def HasProperty(self, key):
        return hasattr(self, key)

    def GetProperty(self, key):
        return getattr(self, key, None)
        
    def AddProperty(self, key, value):
        
        if value is list:
            block = True
            bi = "EOT"
        else:
            block = False
            bi = None
            
        op = pyini.IniOption(self.IniSection, key, value, block)
        op.BlockIdentifier = bi
        
        self.parser.parser.EnterSection(self.IniSection)
        
        #WriteKey(CurrentSection, BlockKey, op)
        self.parser.parser.WriteKey(self.IniSection, key, op)
        setattr(self, key, value)
    
    def RemoveProperty(self, key):
        
        self.parser.parser.EnterSection(self.IniSection)
        
        self.parser.parser.RemoveKey(self.IniSection, key)
        
        if hasattr(self, key):
            delattr(self, key)
        
    def toString(self):
        return self.parser.GetIniAsString()
        
    def Save(self):
        if self.Filename == None:
            fatal("could not write to filename, filename is undefined")
        
        simlib.WriteContents(self.Filename, self.parser.GetIniAsString())
