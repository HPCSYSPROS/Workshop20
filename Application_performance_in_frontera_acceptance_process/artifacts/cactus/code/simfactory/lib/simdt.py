# simtree.py -- A library for parsing ini based decision trees
#
#
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
#
import pyini
import os, shutil
import re
import sys

from libutil import *

import simenv,simlib,simsubs

class DecisionTree:
    
    def __init__(self, tree, parentTree=None):
        self.IniTree = tree
        self.DefineDatabase = simsubs.DefineDatabase()
        
        self.ParentTree = parentTree
    
    def setupTree(self, tree=None):
        
        if tree != None:
            self.IniTree = tree
            
        if self.IniTree == None:
            fatal("no decision tree defined.")
        
        fullFile = simlib.BuildPath(simenv.LIB_PATH, "dt", "%s.ini" % self.IniTree)
        
        if not(os.path.exists(fullFile)):
            fatal("could not open decision tree %s for reading." % fullFile)
        
        info("decision tree ini: %s" % fullFile)
        
        self.Parser = pyini.IniParser(fullFile)
        
        self.loadConfiguration()
        
    def loadConfiguration(self):
        
        if not(self.Parser.HasSection("decisiontree")):
            fatal("invalid decision tree, missing [decisiontree] section")
        
        if not(self.Parser.HasSection("section")):
            fatal("invalid decision tree, missing [section] section")
            
        self.TreeConfiguration = self.Parser.GetSectionAsDict("decisiontree")
        
        dest = self.TreeConfiguration['dest']
        
        if dest == None:
            fatal("invalid tree configuration, no destination defined")
        
        dest = dest.value
        
        (self.DestType, value) = dest.split(":")
        
        if self.DestType not in ['ini']:
            fatal("invalid tree configuration, unknown destination type %s" % type)
        
        eval("self.%s_initializeDestination('%s')" % (self.DestType, value))
        
        self.SectionConfiguration = self.Parser.GetSectionAsDict("section")
        
        if self.SectionConfiguration['fixed'].Value == 'yes':
            self.DefineDatabase.Set(self.extractValue(self.SectionConfiguration['dest'].Value, '@'), 
                self.SectionConfiguration['value'].Value)
            
            info("Section: %s" % self.DefineDatabase.Get('section'))
    
        
    def finishTree(self):
        cmd = "self.%s_finish()" % self.DestType
        return eval(cmd)
        
    def updateValue(self, section, key, value):
        cmd = "self.%s_storeDestinationValue('%s', '%s', '%s')" % (self.DestType, section, key, value)
        eval(cmd)
        
    def extractValue(self, key, sep):
        if key.startswith(sep):
            m = re.match("\%s([^%s]+)\%s" % (sep, sep, sep), key)
            if m == None:
                fatal("invalid key format: %s" % key)
            else:
                return m.group(1)
            
    def begin(self, depth=0):

        c = 'y'
        display("\n")
        
        repeated = False
        
        while c == 'y':
            if self.TreeConfiguration['required'].Value == 'no':
                
                if not repeated:
                    c = self.getInput(self.TreeConfiguration['prompt'].Value, 'N')
                    if c == 'n':
                        return
            else:
                display(self.TreeConfiguration['prompt'].Value)
            
            display("\n")
            
            if self.SectionConfiguration['fixed'].Value == 'no':
                sv = raw_input("%s: " % self.SectionConfiguration['prompt'].Value)
                self.DefineDatabase.Set(self.extractValue(self.SectionConfiguration['dest'].Value, '@'), sv)
                
            sections = self.Parser.GetSections()
            
            ss = list()
            
            for s in sections:
                if s.count(":") > 0:
                    (tt, index) = s.split(":")
                    index = int(index)
                    ss.insert(index-1, s)

            sections = ss
            
            for i in range(0, len(sections)):
                section = sections[i]

                dd = self.Parser.GetSectionAsDict(section)
    
                if section.startswith("action"):
                    self.action(dd)
                    
                if section.startswith("keyvalue:"):
                    self.keyvalue(dd)
                    continue
                
            if self.TreeConfiguration.has_key('next'):
                next_tree = self.TreeConfiguration['next'].value
                
                (tt, nextTree) = next_tree.split(":")
                
                if tt not in ['decisiontree']:
                    fatal("unknown tree next destination %s" % tt)
                    
                tree = DecisionTree(nextTree, self)
                tree.setupTree()
                tree.begin(depth+1)
    
            if self.TreeConfiguration['repeat'].Value == 'yes':
                display("\n")
                c = self.getInput(self.TreeConfiguration['repeatstatement'].Value, 'N')
                repeated = True
            else:
                c = 'n'
                
        if depth == 0:
            self.finishTree()
                
    def dictSortFunc(self, x, y):
        
        if x.count(":") == 0 or y.count(":") == 0:
            return 0
            
        xi = int(x.split(":").pop())
        yi = int(y.split(":").pop())
        return xi - yi
    
    def action(self, adict):
        
        prompt = self.DefineDatabase.SubAll(adict['prompt'].Value)
        
        action = self.DefineDatabase.SubAll(self.extractValue(adict['action'].Value, '%'))
        
        if adict.has_key("check"):
            check = self.DefineDatabase.SubAll(self.extractValue(adict['check'].Value, '%'))
        
            results = eval("self.macro_%s" % check)
            
            if not results:
                return
        
        results = eval("self.macro_%s" % action)
        
        if adict['printresults'].Value == 'yes':
            display("%s: %s" % (prompt, results))
        else:
            display(prompt)
        
        if adict.has_key('dest'):
            dest = self.extractValue(adict['dest'].Value, '@')
            self.DefineDatabase.Set(dest, results)
    
    
    def keyvalue(self, vdict):
        c = 'y'
        
        repeat = 'no'
        
        if vdict.has_key('repeat'):
            repeat = vdict['repeat'].Value
        
        while c == 'y':
            dest = self.DefineDatabase.SubAll(vdict['dest'].Value)
            
            self.performOperation(vdict, 'key')
            
            key = self.DefineDatabase.Get('key')
            
            if len(key) == 0:
                warning("length of key is 0")
                return
            
            self.performOperation(vdict, 'value')
            self.updateValue(dest, key, self.DefineDatabase.Get('value'))

            self.DefineDatabase.Set(key, self.DefineDatabase.Get('value'))
            
            if repeat == 'yes':
                rs = self.DefineDatabase.SubAll(vdict['repeatstatement'].Value)
                c = self.getInput(rs, 'N')
            else:
                c = 'n'
    
    def getInput(self, prompt, default='N'):
        
        Y = 'Y'
        N = 'N'
        
        if Y == default:
            Y = "Y*"
        if N == default:
            N = "N*"
        
        defaultString = "[%s/%s]" % (Y, N)
        
        acceptable = ['y', 'n', 'yes', 'no']
        
        c = raw_input("%s %s: " % (prompt, defaultString)).lower().strip()

        if len(c.strip()) == 0:
            return default.lower()
        
        while c not in acceptable:
            display("invalid input, specify Y or N")
            c = raw_input("%s %s: " % (prompt, defaultString)).lower().strip()
        
        return c
        
    def performOperation(self, vdict, type):
        
        default = str()

        if vdict.has_key("%sdefault" % type):
            default = vdict["%sdefault" % type].Value
            
            if default.startswith("option:"):
                option = default.split(":").pop()
                
                if simenv.OptionsManager.HasOption(option):
                    default = simenv.OptionsManager.GetOption(option)
                else:
                    default = str()
                    
            if default.count("%") > 0:
                macro = self.DefineDatabase.SubAll(self.extractValue(default, '%'))
                default = eval("self.macro_%s" % macro)
            
            if default.count("@") > 0:
                default = self.DefineDatabase.SubAll(default)
            
        self.DefineDatabase.Set('default', default)

        if not(vdict.has_key('%sdest' % type)):
            fatal("invalid keyvalue section, missing key %sdest" % type)
            
        vdest = self.extractValue(vdict['%sdest' % type].Value, '@')
        
        if not(vdict.has_key("%sprompt" % type)):
            results = vdict['%svalue' % type].Value
            
            if results.startswith("option:"):
                option = results.split(":").pop()
                
                if simenv.OptionsManager.HasOption(option):
                    results = simenv.OptionsManager.GetOption(option)
                else:
                    results = default
            
        else:
            prompt = self.DefineDatabase.SubAll(vdict['%sprompt' % type].Value)
            results = raw_input("%s: " % prompt).strip()
        
            if len(results) == 0:
                results = default

        self.DefineDatabase.Set(vdest, results)
        
    ############ INI #################
    
    def ini_initializeDestination(self, value):
    
        fullPath = simlib.BuildPath(simenv.BASE_PATH, value)
        
        self.ini_FullPath = fullPath
        
        if self.ParentTree != None:
            if self.ParentTree.DestType == 'ini':
                if self.ParentTree.ini_FullPath == self.ini_FullPath:
                    self.ini_DestParser = self.ParentTree.ini_DestParser
                    return
        
        self.ini_DestParser = pyini.IniParser(self.ini_FullPath)
            
    def ini_storeDestinationValue(self, section, key, value):
        if section == 'None' or len(section) == 0:
            section = None
        
        if section != None:
            self.ini_DestParser.parser.InitSection(section)
        
        #IniOption(self,section, key, value, IsBlock):
        
        io = pyini.IniOption(section, key, value, False)
        
        self.ini_DestParser.parser.WriteKey(section, key, io)
        
    def ini_finish(self):
        display("\n--------------------SUMMARY--------------------:\n")
        display(self.ini_DestParser.GetIniAsString())
        display("\n------------------END SUMMARY------------------:\n")
    
        if self.TreeConfiguration.has_key('save'):
            save = self.TreeConfiguration['save'].value
            
            if save == "noprompt":
                r = "y"
            else:
                r = self.getInput("Save contents", 'Y')
        else:
            r = self.getInput("Save contents", 'Y')

        if r == "y":
            simlib.WriteContents(self.ini_FullPath, self.ini_DestParser.GetIniAsString())
            display("Contents successfully written to %s" % self.ini_FullPath)
    
    # ini specific macro 
    def ini_GET_KEY_DEFAULT(self, section, key):
    
        if section == 'None' or len(section) == 0:
            section = None
            
        default = str()
        
        vv = self.ini_DestParser.GetOption(section, key)
        
        if vv != None:
            default = vv.Value
        
        return default
        
    ############ MACROS ##############
    
    def macro_GET_KEY_DEFAULT(self, key):
        
        section = self.DefineDatabase.Get('section')
        
        cmd = "self.%s_GET_KEY_DEFAULT('%s', '%s')" % (self.DestType, section, key)
        return eval(cmd)
        
    def macro_GET_DEFAULT_USERNAME(self):
    
        if self.DestType == 'ini':
            vv = self.ini_GET_KEY_DEFAULT('default', 'user')
            
            if vv == None or len(vv) == 0:
                key_tests = ['LOGNAME', 'USER']
                for key in key_tests:
                    if os.environ.has_key(key):
                        return os.environ[key]
            else:
                return vv
    
    def macro_GET_BASEDIR(self):
        return simlib.BuildPath(os.environ['HOME'], "simulations")
    
    def macro_GET_SOURCEBASEDIR(self):
        
        if simenv.OptionsManager.HasOption('setup-sourcebasedir'):
            sourcebasedir = simenv.OptionsManager.GetOption('setup-sourcebasedir')
        else:
            cpath = simenv.CACTUS_PATH
            parts = cpath.split(os.sep)
            parts.pop()
            
            sourcebasedir = os.sep.join(parts)
        
        return sourcebasedir
                
    def macro_GET_LOCAL_MACHINE_NAME(self):
        return simlib.GetHostNameAlias()
    
    def macro_CHECK_CREATE_MACHINE(self):
        
        mm = self.DefineDatabase.Get('machine')
        
        if simlib.MachineExists(mm):
            return False
            
        return True
        
    def macro_CREATE_MACHINE(self):
        mm = self.DefineDatabase.Get('machine')
        
        srcPath = simlib.BuildPath(simenv.MDB_PATH, 'generic.ini')
        dstPath = simlib.BuildPath(simenv.MDB_PATH, '%s.ini' % mm)
        
        shutil.copy(srcPath, dstPath)
        
        self.tParser = pyini.IniParser(dstPath)
        self.tParser.parser.RenameSection('generic', mm)

        options = ['name', 'hostname', 'nickname']
        
        for oo in options:
            io = self.tParser.GetOption(mm, oo)
            io.Value = io.value = mm
            self.tParser.parser.WriteKey(mm, oo, io)

        # sourcebasedir
        io = self.tParser.GetOption(mm, "sourcebasedir")

        sourcebasedir = self.macro_GET_SOURCEBASEDIR()
        
        io.Value = sourcebasedir
        self.tParser.parser.WriteKey(mm, "sourcebasedir", io)
        
        # basedir
        io = self.tParser.GetOption(mm, "basedir")
        io.Value = self.macro_GET_BASEDIR()
        self.tParser.parser.WriteKey(mm, "basedir", io)

        # optionlist
        if (simenv.OptionsManager.RawOptionDefined("optionlist")):
            io = self.tParser.GetOption(mm, "optionlist")
            io.Value = simlib.GetOptionList(False)
            self.tParser.parser.WriteKey(mm, "optionlist", io)
        
        # submitscript
        if (simenv.OptionsManager.RawOptionDefined("submitscript")):
            io = self.tParser.GetOption(mm, "submitscript")
            io.Value = simlib.GetSubmitScript(False)
            self.tParser.parser.WriteKey(mm, "submitscript", io)
        
        # runscript
        if (simenv.OptionsManager.RawOptionDefined("runscript")):
            io = self.tParser.GetOption(mm, "runscript")
            io.Value = simlib.GetRunScript(False)
            self.tParser.parser.WriteKey(mm, "runscript", io)
        
        # ppn
        if (simenv.OptionsManager.RawOptionDefined("ppn")):
            io = self.tParser.GetOption(mm, "ppn")
            io.Value = simenv.OptionsManager.GetOption('ppn')
            self.tParser.parser.WriteKey(mm, "ppn", io)

            # also set 'max-num-threads'
            io = self.tParser.GetOption(mm, "max-num-threads")
            io.Value = simenv.OptionsManager.GetOption('ppn')
            self.tParser.parser.WriteKey(mm, "max-num-threads", io)
        
            # also enable parallel make
            io = self.tParser.GetOption(mm, "make")
            io.Value = "nice make -j%d" % simenv.OptionsManager.GetOption('ppn')
            self.tParser.parser.WriteKey(mm, "make", io)

        # num-threads
        if (simenv.OptionsManager.RawOptionDefined("num-threads")):
            io = self.tParser.GetOption(mm, "num-threads")
            io.Value = simenv.OptionsManager.GetOption('num-threads')
            self.tParser.parser.WriteKey(mm, "num-threads", io)
        
        simlib.WriteContents(dstPath, self.tParser.GetIniAsString())
        return "machine %s [%s] created successfully" % (mm, dstPath)
        
