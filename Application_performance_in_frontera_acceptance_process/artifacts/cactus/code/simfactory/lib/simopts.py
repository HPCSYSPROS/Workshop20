# simopts.py -- provide central repository/interface for importing/exposing known commandline options
#               across several utilities

# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University
#
#
import os
import sys
import pyini
import simenv
import optparse

from copy import copy
from libutil import *

from optparse import OptionParser
from optparse import OptionGroup
from optparse import Option, OptionValueError

def check_nonnull(option, opt, value):
    if value == None or len(value.strip()) == 0:
        raise OptionValueError("option %s, opt %s: argument cannot be null" % (option, opt))
    else:
        return value
            
class TypeNonNull (Option):
    TYPES = Option.TYPES + ("nonnull",)
    TYPE_CHECKER = copy(Option.TYPE_CHECKER)
    TYPE_CHECKER["nonnull"] = check_nonnull

class OptionManager:
    def __init__(self, optiondir):
        self.optionDirectory = optiondir
        self.parsers = dict()
        self.usedGroups = ["common"]
        self.usageString = "usage: %prog [options] arg"
        self.optionParser = None
        self.options = None
        self.args = None
        self.optionParser = None
        self.OptionDictionaries = dict()
        self.Macros = dict()
        self.MDBKeys = dict()
        self.Substitutions = list()
        self.Replacements = list()
        self.Appends = list()
        self.RawOptionString = " ".join(sys.argv)
        self.Aliases = dict()
        
        # find and init our option groups.
        
        self.Discovery()
    
    def PrintHelp(self):
        self.optionParser.print_help()
    
    def RawOptionDefined(self, option):
        opts = self.RawOptionString.split()
        #return any(map(lambda s: s.startswith("--%s" % option), opts))
        for opt in opts:
            if opt.startswith("--%s" % option): return True
        return False
    
    def IsNegated(self, option):
        opts = self.RawOptionString.split()
        #return any(map(lambda s: s.startswith("--no%s" % option), opts))
        for opt in opts:
            if opt.startswith("--no%s" % option): return True
        return False
    
    def BuildOptionString(self, stripArguments=None):
        if self.options == None:
            return ""
        
        keys = self.options.__dict__.keys()
        
        ostr = str()
        
        for key in keys:
            if stripArguments != None and key in stripArguments:
                continue
            
            if key.startswith("__"):
                continue
            
            if key.startswith("["):
                continue
            
            if key in ['define', 'append', 'substitute', 'replace', 'mdbkey']:
                continue
            
            sdict = self.OptionDictionaries[key]
            
            value = getattr(self.options, key)
            
            if value == None:
                continue
            
            if self.Aliases.has_key(key):
                alias = self.Aliases[key]
                if self.RawOptionDefined(alias):
                    key = alias
            
            # this check actually fails if value is not explicitly checked against Boolean True.
            if sdict['arg'].value == 'no':
                if value and self.RawOptionDefined(key):
                    ostr = "%s --%s" % (ostr, key)
                
                if not value and self.IsNegated(key):
                    ostr = "%s --no%s" % (ostr, key)
            else:
                if self.RawOptionDefined(key):
                    ostr = "%s --%s=%s" % (ostr, key, value)
        
        for key in self.Macros:
            ostr = "%s --define=%s=%s" % (ostr, key, self.Macros[key])
        
        for item in self.Substitutions:
            ostr = "%s --substitute=%s" % (ostr, item)
        
        for item in self.Replacements:
            ostr = "%s --replace=%s" % (ostr, item)
        
        for item in self.Appends:
            ostr = "%s --append=%s" % (ostr, item)
        
        for key in self.MDBKeys:
            ostr = "%s --mdbkey=%s=%s" % (ostr, key, self.MDBKeys[key])
        
        return ostr.strip()
    
    # -- ACCESSOR METHODS --
    def HasOption(self, option):
        o = getattr(self.options, option, None)
        
        if o == None:
            return False
        
        return True
    
    def GetOption(self, option):
        return getattr(self.options, option, None)
    
    def UpdateOption(self, option, value):
        setattr(self.options, option, value)
    
    # -- BUILD METHODS --
    
    def Discovery(self):
        
        try:
            os.stat(self.optionDirectory)
        except OSError:
            fatal("while attempting to read option dir, could not open %s for reading" % self.optionDirectory)
        
        files = os.listdir(self.optionDirectory)
        
        for f in files:
            if f.endswith(".ini"):
                parts = os.path.splitext(f)
                basename = parts[0]
                full_path = "%s/%s" % (self.optionDirectory, f)
                self.parsers[basename] = pyini.IniParser(full_path)
    
    def SetUsage(self, ss):
        self.usageString = ss
    
    def HasGroup(self, group):
        return self.parsers.has_key(group)
    
    def UseGroup(self, group):
        if self.HasGroup(group):
            self.usedGroups.append(group)
            return
        
        fatal("when specifying option group %s: group %s does not exist" % (group, group))
    
    def BuildBuiltins(self):
        
        # macro with callback
        
        self.optionParser.add_option("--define", action="callback", type="string", dest='define', nargs=2, help="set additional definition",
        callback=self.callback)
        self.optionParser.add_option("--mdbkey", action="callback", type="string", dest='mdbkey', nargs=2, help="override an mdb key",
        callback=self.callback)
        self.optionParser.add_option("--substitute", action="callback", type="string", dest='substitute', nargs=2, help="perform regex substitution",
        callback=self.callback)
        self.optionParser.add_option("--replace", action="callback", type="string", dest='replace', nargs=1, help="perform string replacement",
        callback=self.callback)
        self.optionParser.add_option("--append", action="callback", type="string", dest='append', nargs=1, help="append to a variable",
        callback=self.callback)
    
    def callback(self, option, opt, value, parser):
        
        if opt == "--define":
            if value[0].startswith("@"):
                value[0] = value[0].replace("@", "").strip()
            
            self.Macros[value[0]] = value[1]
        
        if opt == "--substitute":
            self.Substitutions.append(value)
        
        if opt == "--replace":
            self.Replacements.append(value)
        
        if opt == "--append":
            self.Appends.append(value)
        
        if opt == "--mdbkey":
            self.MDBKeys[value[0]] = value[1]
    
    def Build(self):
        self.optionParser = OptionParser(usage=self.usageString, option_class=TypeNonNull)
        
        self.BuildBuiltins()
        
        for group in self.usedGroups:
            if group == "common":
                self.AddOptions(self.optionParser, group)
            else:
                og = OptionGroup(self.optionParser, "Options from %s" % group)
                self.AddOptions(og, group)
                self.optionParser.add_option_group(og)
        
        (self.options, self.args) = self.optionParser.parse_args()
    
    def AddOptions(self, op, group):
        if not(self.HasGroup(group)):
            fatal("Could not retreive ini parser for option group %s" % group)
        
        parser = self.parsers[group]
        
        options = parser.GetSections()
        
        for o in options:
            sdict = parser.GetSectionAsDict(o)
            self.OptionDictionaries[o] = sdict
            self.AddOption(op, o, sdict)
    
    def AddOption(self, op, option, sdict):
        
        if sdict.has_key('alias'):
            alias_option = sdict['alias'].value
            self.Aliases[alias_option] = option
        
        if sdict['negatable'].value == 'yes' and sdict.has_key('short'):
            fatal("Error parsing option %s: cannot negate a short option" % option)
        
        if sdict.has_key("short"):
            optstring = '"%s", "%s"' % (sdict['short'].value, sdict['long'].value)
        else:
            optstring = '"%s"' % sdict['long'].value
        
        if sdict.has_key('alias'):
            dest = sdict['alias'].value
        else:
            dest = option
        
        optstring = self.AddOptionToOptionString(optstring, 'dest', dest)
        
        if sdict['arg'].value == "yes":
            optstring = self.AddOptionToOptionString(optstring, 'type', sdict['argtype'].value)
            optstring = self.AddOptionToOptionString(optstring, 'metavar', sdict['argformat'].value)
            if sdict.has_key('multiarg') and sdict['multiarg'].value == "yes":
                optstring = self.AddOptionToOptionString(optstring, 'action', 'append')
        else:
            optstring = self.AddOptionToOptionString(optstring, 'action', 'store_true')
            
            #dprint("storing_true for argument %s" % option)
            
            if sdict.has_key('default'):
                
                optstring = '%s,%s=%s' % (optstring, 'default', CoerceBool(sdict['default'].value))
        
        help_str = sdict['desc'].value
        
        if sdict.has_key('default'):
            help_str = "%s %s" % (help_str, '[default: %default]')
        
        if sdict['negatable'].value == 'yes':
            optstring = '%s,%s=%s' % (optstring, 'help', 'optparse.SUPPRESS_HELP')
        else:
            optstring = self.AddOptionToOptionString(optstring, 'help', help_str)
        
        cmd = "op.add_option(%s)" % optstring
        
        #dprint(cmd)
        #dprint('cmd: %s' % cmd)
        
        eval(cmd)
        
        if sdict['negatable'].value == 'yes':
            self.AddMetaOption(op, option, sdict)
            self.AddNegatedOption(op, option, sdict, "no")
            self.AddNegatedOption(op, option, sdict, "no-")
        
        return
    
    def AddNegatedOption(self, op, option, sdict, prefix):
        
        optstring = '"--%s%s"' % (prefix, option)
        
        optstring = self.AddOptionToOptionString(optstring, 'dest', option)
        optstring = self.AddOptionToOptionString(optstring, 'action', 'store_false')
        #optstring = self.AddOptionToOptionString(optstring, 'default', False)
        optstring = '%s,%s=%s' % (optstring, 'help', 'optparse.SUPPRESS_HELP')
        
        
        cmd = "op.add_option(%s)" % optstring
        #dprint('cmd: %s' % cmd)
        
        eval(cmd)
    
    def AddMetaOption(self, op, option, sdict):
        
        sdict['long'].value = "--[no]%s" % option
        
        optstring = '"%s"' % sdict['long'].value
        
        if sdict.has_key('default'):
            optstring = '%s,%s=%s' % (optstring, 'default', CoerceBool(sdict['default'].value))
        
        help_str = sdict['desc'].value
        
        if sdict.has_key('default'):
            help_str = "%s %s" % (help_str, '[default: %default]')
        
        optstring = self.AddOptionToOptionString(optstring, 'dest', option)
        optstring = self.AddOptionToOptionString(optstring, 'action', 'store_false')
        optstring = self.AddOptionToOptionString(optstring, 'metavar', option)
        optstring = self.AddOptionToOptionString(optstring, 'help', help_str)
        
        cmd = "op.add_option(%s)" % optstring
        #dprint('cmd: %s' % cmd)
        eval(cmd)
    
    def AddOptionToOptionString(self, os, option, value):
        return '%s,%s="%s"' % (os, option, value)

