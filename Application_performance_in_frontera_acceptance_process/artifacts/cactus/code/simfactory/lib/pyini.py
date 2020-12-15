# PYINI -- Python INI file parser with support for <<(ID)..(ID) block comments
# 
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
# 
# 

import sys
import re
import os

from libutil import *

global parser
parser = None

class IniParser:

    def __init__(self, filename=None):
    
        global parser
        parser = self
        
        self.filename = filename
        
        if self.filename:
            if os.path.exists(self.filename):
                self.ParseFile()
                return
        
        self.InitBlankParser()
    
    def InitBlankParser(self):
        self.parser = InternalParser()
    
    def PerformParse(self):
        self.parser = InternalParser(self.fptr)
        self.parser.Parse()
        
    def OpenFile(self):
        try:
            self.fptr = open(self.filename, 'r')
        except IOError:
            fatal("Could not open %s for reading" % self.filename)
    
    def ParseFile(self):
        self.OpenFile()
        self.PerformParse()
        self.fptr.close()
    
    def UpdateFromDict(self, optdict):
        
        
        keys = optdict.keys()
        if len(keys) == 0:
            return
        
        sections = self.GetSections()
        
        for section in sections:
            for key in keys:
                value = optdict[key]
                #IniOption(self, section, key, value, IsBlock):
                op = IniOption(section, key, value, False)
                #dprint("overrode %s on mdb entry %s with value %s" % (key, section, value))
                self.parser.WriteKey(section, key, op)
                
    def UpdateFromIni(self, iniFile, importSections=False):
        
        uIniParser = self.__class__(iniFile)
            
        sections =  uIniParser.GetSections()
        
        # if they are global options, we merrrrgeeeeee very carefully.
        
        gkeys = uIniParser.parser.KeyOrder['__globalK__']
        
        for gk in gkeys:
            gio = uIniParser.GetGlobalOption(gk)
            
                
            #dprint("getting option %s" % gk)
            io = self.GetGlobalOption(gk)
            
            if io == None:
                continue
            
            if type(io.Value) == list:
                #convert to list
                gio.ConvertToList()
                
                ll = io.Value
                ll.extend(gio.Value)
                io.UpdateValue(ll)
                #dprint("%s: after extending, io is now: %s" % (self, io.Value))
                self.parser.WriteKey(None, gk, io)
            else:
                io.Value = gio.value
                self.parser.WriteKey(None, gk, io)
                
        for s in sections:
            if s == 'default':
                continue
                
            keys = uIniParser.parser.KeyOrder[s]
                
            if importSections:
                # if import is set to true, then create the section
                # this function will not clobber an already existing section, so its safe
                # to use blindly
                self.parser.InitSection(s)
                            
            if self.HasSection(s):
                for k in keys:
                    io = uIniParser.GetOption(s, k)
                    #dprint("updating %s.%s" % (s, k))
                    self.parser.WriteKey(s, k, io)
                continue
        
        if "default" in sections:
            keys = uIniParser.GetSectionKeys("default")
            ess = self.GetSections()
            
            for ss in ess:
                for k in keys:
                    if not(self.HasOption(ss, k)):
                        io = uIniParser.GetOption("default", k)
                        self.parser.WriteKey(ss, k, io)
        
            
    # more "readable" alias for self.CheckExistence
    def HasOption(self, section, key):
        return self.CheckExistence(section, key)
    
    def HasSection(self, section):
        return self.parser.Section.has_key(section)
    
    def GetOptionByPath(self, path):
        # path is in section.key format
        if path.count(".") == 0:
            return self.GetGlobalOption(path)
        else:
            parts = path.split(".")
            return self.GetSectionOption(parts[0], parts[1])
    
    
    def GetSectionAsDict(self, section):
        if self.HasSection(section):
            return self.parser.Section[section]
        return None
        
    def GetOption(self, section, key):
        if section == None:
            if self.CheckExistence(None, key):
                return self.parser.Globals[key]
        else:
            if self.CheckExistence(section, key):
                return self.parser.Section[section][key]
                
    def GetGlobalOption(self, key):
        return self.GetOption(None, key)
        
    def GetSectionOption(self, section, key):
        return self.GetOption(section, key)
        
    def CheckExistence(self, section, key):
        if section == None:
            if self.parser.Globals.has_key(key):
                return True
        else:
            if self.parser.Section.has_key(section):
                if key == None:
                    return True
                
                if self.parser.Section[section].has_key(key):
                    return True
        return False
    
    def GetSections(self):
        return self.parser.Section.keys()
    
    def GetGlobalKeys(self):
        return self.GetKeys(None)
        
    def GetSectionKeys(self, section):
        return self.GetKeys(section)
    
    def GetKeys(self, section):
        if section == None:
            return self.parser.Globals.keys()
        
        if self.CheckExistence(section, None):
            return self.parser.Section[section].keys()
        
        return None
    
    
    def GetIniAsString(self):
        return self.parser.PrintParsedIni()
        
    def PrintIni(self):
        display(self.parser.PrintParsedIni())
    
    def PrintSection(self, section):
        if not(self.HasSection(section)):
            fatal("cannot print section %s, unknown section" % section)
            
        display(self.parser.PrintSection(section))

class IniOption:
    def __init__(self, section, key, value, IsBlock):
        #dprint("making ini option %s.%s with value %s" % (section, key, value))
        global parser
        
        if value == None:
            value = ""
        
        if value is not str:
            value = str(value)
        
        self.Parser = parser    
        self.IsQuoted = False
        self.Key = key
        self.IsBlock = IsBlock
        self.BlockIdentifier = None
        self.Section = section
        
        self.LeadingComment = None
        self.InlineComment = None
        
        # strip any trailing non-quoted space off the end of the string.
        
        self.OriginalValue = value.rstrip()
        
        self.InterpretQuotes()
        
        self.key = self.Key
        self.value = self.Value
    
    def UpdateValue(self, value):
        self.Value = self.value = value
        
    def ConvertToList(self):
    
        if type(self.Value) == list:
            return
            
        if self.Value == None:
            self.UpdateValue([])
            return
        
        if self.IsBlock:
            value = self.Value
            lines = value.split("\n")
        else:
            if len(self.Value) == 0:
                lines = []
            else:
                lines = [self.Value]
                
        self.UpdateValue(lines)
        
    def InterpretQuotes(self):
        if self.OriginalValue == None or len(self.OriginalValue) == 0:
            self.Value = self.OriginalValue
            return
            
        if self.OriginalValue[0] == "'" or self.OriginalValue[0] == '"':
            #dprint("value %s has quotes" % value)
            try:
                vv = eval(self.OriginalValue)
                self.Value = vv
                self.IsQuoted = True
            except Exception:
                error = "Could not interpret quoted string: %s" % self.OriginalValue
                fatal("Syntax error on line %s of %s: %s" % (self.Parser.parser.CurrentLineNumber, self.Parser.filename, error))
        else:
            self.Value = self.OriginalValue
        

class IniSyntaxChecker:
    def __init__(self, syntaxfile, iniParser, importDefaults=False):
    
        if type(iniParser).__name__ == 'str':
            self.IniParser = IniParser(iniParser)
        else:
            self.IniParser = iniParser
            
        self.SyntaxFile = syntaxfile
        self.ImportDefaults = importDefaults

    def SyntaxCheck(self):
        # init our parser
        # self.IniParser = IniParser(self.IniFile)

        self.SyntaxParser = IniParser(self.SyntaxFile)
        
        # lets get the sections from the Syntax Parser
        # each section represents a key that each section inside the 
        # ini being checked should possess.
        
        syntax_keys = self.SyntaxParser.GetSections()
        
        ini_sections = self.IniParser.GetSections()
        
        if len(ini_sections) == 0:
            warning("checked ini %s contains no sections" % self.IniParser.filename)
            return True
        
        for ini_section in ini_sections:
            for skey in syntax_keys:
                io = self.IniParser.GetOption(ini_section, skey)
                #if io == None:
                    #dprint("Warning: None for GetOption on section %s, key %s" % (ini_section, skey))
                    
                # lets gather some information about it.
                entry = self.GetSyntaxEntry(skey)
                
                
                # Self-test

                #die if $necessity eq 'required' and defined $keydesc->{'default'};
                if entry['necessity'] == 'required' and entry.has_key('default'):
                    fatal("necessity for key %s is required, but a default entry was provided" % skey)
                
                #die if $type ne 'string' and defined $keydesc->{'pattern'};
                if entry['type'] != 'string' and entry.has_key('pattern'):
                    fatal("cannot define a pattern on a non-string key %s" % skey)
                
                #die if defined $keydesc->{'default'} and
                #   $type eq 'string' and defined $keydesc->{'pattern'} and
                #   $default !~ $pattern;
                if entry.has_key('default') and entry['type'] == 'string' \
                    and entry.has_key('pattern') and not(self.PatternMatches(entry['pattern'], entry['default'])):
                    
                    fatal("specified default value %s for key %s does not match converted pattern %s" \
                        % (entry['default'], skey, entry['pattern']))

                #die if $type eq 'string' and defined $keydesc->{'pattern'} and
                #   $example !~ $pattern;
                if not(entry.has_key('example')):
                    fatal("specified key %s does not have an example defined")
                
                if entry['type'] == 'string' and entry.has_key('pattern') and not(self.PatternMatches(entry['pattern'], entry['example'])):
                    fatal("specified example value %s for key %s does not match converted pattern %s" \
                        % (entry['example'], skey, entry['pattern']))
                
                
                # moving along -- type check our actual entry

                # check if required
                
                if entry['necessity'] == 'required' and io == None:
                    fatal("required key %s in section %s is missing" % (skey, ini_section))
                
                # type check
                # if the type is a string/any, then we just accept anything.
                if entry.has_key('type') and (entry['type'] != 'string' and entry['type'] != 'any') and io != None:
                    val_type = entry['type']
                    
                    value = io.Value
                    
                    if val_type == "int":
                        try:
                            int(value)
                        except ValueError:
                            fatal("value for key %s in section %s is not of required type %s" % (skey, ini_section, val_type))
                    
                    if val_type == "double" or type == "float":
                        try:
                            float(value)
                        except TypeError:
                            fatal("value for key %s in section %s is not of required type %s" % (skey, ini_section, val_type))
                        except ValueError:
                            fatal("value for key %s in section %s is not of required type %s" % (skey, ini_section, val_type))
                            
                # check if matches pattern if there is a pattern
                if entry.has_key('pattern') and io != None:
                    value = io.Value
                    if not(self.PatternMatches(entry['pattern'], value)):
                        fatal("specified value %s for key %s does not match converted pattern %s" 
                        % (value, skey, entry['pattern']))
                
                # import default
                if self.ImportDefaults:
                    #dprint("importing defaults for section %s, key %s" % (ini_section, skey))
                    if io == None:
                        #dprint("for section %s, key %s, importing a default value.." % (ini_section, skey))
                        if entry.has_key('default'):
                            io = IniOption(ini_section, skey, entry['default'], False)
                            # if it's default 
                            #dprint("Importing default value %s into section %s, key %s" % (entry['default'], ini_section, skey))
                            self.IniParser.parser.WriteKey(ini_section, skey, io)
                    
                
    
    def PatternMatches(self, pattern, value):
        cpattern = pattern
        
        #dprint("matching %s against pattern %s" % (value, cpattern))
        try:
            p = re.compile(cpattern)
        except:
            fatal("pyini.py[%s]: could not interpret pattern %s, converted from orig pattern %s" % (LineNumber(), cpattern, pattern))
        
        return p.search(value) != None
        
    def GetSyntaxEntry(self, skey):
        sdict = self.SyntaxParser.GetSectionAsDict(skey)
        if sdict == None:
            fatal("Non-recoverable error encountered when attempting to retrieve syntax dict for %s" % skey)
        
        sentry = dict()
        
        keys = sdict.keys()
        
        #translate from our IniOption class down to a flat dict.
        #we won't need any of the extraneous information
        for key in keys:
            sentry[key] = sdict[key].Value
        
        return sentry
        
class InternalParser:
    def __init__(self, fptr=None):
        self.SectionOrder = list()
        self.KeyOrder = dict()
        self.KeyOrder['__globalK__'] = list()
        
        self.Section = dict()
        self.Globals = dict()
        self.Keys = dict()
        
        if fptr != None:
            self.Lines = fptr.readlines()
        else:
            self.Lines = []
        
        self.CurrentSection = None
        self.BlockMode = False
        self.BlockKey = None
        self.BlockText = None
        self.BlockModeBeginLine = 0
        self.BlockIdentifier = None
        self.CurrentLineNumber = 0
        
    # -- PRINTING THE INI --
    
    def PrintParsedIni(self):
        
        build_str = ""
        gkeys = self.KeyOrder['__globalK__']
        
        for key in gkeys:
            value = self.Globals[key]
            
            if value.IsBlock:
                build_str = "%s%s\n" % (build_str, self.PrintBlock(key, value))
            else:
                build_str = "%s%s\n" % (build_str, self.PrintKey(key, value))
                        
        for s in self.SectionOrder:
            build_str = "%s%s" % (build_str, self.PrintSection(s))
        
        return build_str
    
    def PrintSection(self, section):
        build_str = "\n[%s]\n" % section
        
        keys = self.KeyOrder[section]
        
        for key in keys:
            value = self.Section[section][key]
            if value.IsBlock:
                build_str = "%s%s\n" % (build_str, self.PrintBlock(key, value))
            else:
                build_str = "%s%s\n" % (build_str, self.PrintKey(key, value))
        
        
        return build_str
        
    def PrintKey(self, key, value):
        
        if value.IsQuoted:
            build_str = "%-16s= %s" % (key, value.OriginalValue)
        else:
            build_str = "%-16s= %s" % (key, value.Value)

        if value.LeadingComment != None:
            build_str = "\n%s\n%s" % (value.LeadingComment, build_str)

        return build_str
    
    def BuildSpaces(self, num_spaces):
    
        ss = str()
        for i in range(num_spaces):
            if i == 0:
                ss = " "
            else:
                ss = "%s " %ss
        
        return ss
        
    def PrintBlock(self, key, value):

        build_str = "%s = <<%s\n%s\n%s\n" % (key, value.BlockIdentifier, value.Value, value.BlockIdentifier)
        
        if value.LeadingComment != None:
            build_str = "\n%s\n%s" % (value.LeadingComment, build_str)
            
        return build_str
        
    # -- PARSER --
    
    def WriteKey(self, section, key, value):

        # if there is no current section, write it to our globals dictionary
        if section == None:
            #dprint("writing key %s to globals" % key)
            self.Globals[key] = value
            
            if not(self.KeyOrder.has_key('__globalK__')):
                self.KeyOrder['__globalK__'] = list()
            
            if key not in self.KeyOrder['__globalK__']:
                self.KeyOrder['__globalK__'].append(key)
                
            #dprint("globals.%s = %s" % (key, value))
        else:
            #dprint("writing key %s to section %s" % (key, self.CurrentSection))
            self.Section[section][key] = value
            
            if key not in self.KeyOrder[section]:
                self.KeyOrder[section].append(key)
                
            #dprint("%s.%s = %s" % (self.CurrentSection, key, self.Section[self.CurrentSection][key]))
    
    def RemoveKey(self, section, key):
        if section == None:
            if self.Globals.has_key(key):
                del self.Globals[key]
            
            if not(self.KeyOrder.has_key('__globalK__')):
                self.KeyOrder['__globalK__'] = list()
            
            if key in self.KeyOrder['__globalK__']:
                self.KeyOrder['__globalK__'].pop(self.KeyOrder['__globalK__'].index(key))
        else:
            
            if self.Section[section].has_key(key):
                del self.Section[section][key]
            
            if key in self.KeyOrder[section]:
                self.KeyOrder[section].pop(self.KeyOrder[section].index(key))
                
    def InitSection(self, section):
        if not(self.Section.has_key(section)):
            self.Section[section] = dict()
            self.KeyOrder[section] = list()
            self.SectionOrder.append(section)
            
    
    def EnterSection(self, section):
        #dprint("entering section %s" % section)
        self.CurrentSection = section
        self.InitSection(section)
        
        if section not in self.SectionOrder:
            self.SectionOrder.append(section)
    
    def RenameSection(self, old, new):
        if not(self.Section.has_key(old)):
            return

        self.Section[new] = self.Section[old]
        self.KeyOrder[new] = self.KeyOrder[old]
        
        del self.Section[old]
        del self.KeyOrder[old]
        
        self.SectionOrder[self.SectionOrder.index(old)] = new
        
    def HandleError(self, line):
        # extend this to interpret each error to try to give better feedback as to the
        # nature of the error encountered
        
        error = "Unknown statement: %s" % line
        
        if line.count("["):
            if line.count("]"):
                error = "no whitespace allowed in section identifier: %s" % line
            else:
                error = "Missing ] on section identifier: %s" % line
        
        if line.count("="):
            error = "no whitespace allowed in key identifier: %s" % line
        
        fatal(self.MakeError(self.CurrentLineNumber, error))

    def MakeError(self, line_number, error):
        global parser
        return "Syntax error on line %s of %s: %s" % (line_number, parser.filename, error)
        
    def StripInlineComment(self, line):
        
        cchars = ['#', ';']
        
        found = None
        for cc in cchars:
            if line.count(cc) > 0:
                found = cc
            
            if found == None:
                return line
        
        parts = line.split(found)
        
        return parts[0].rstrip()
    
    
    
    
    def Parse(self):
        
        # match [<word>]
        section_re = '\[(\S+)\]'
        
        # match key = value up to a comment [#]
        key_re = '(\S+)\s*=\s*([^#]*)[#]?.*$'
        
        # match block statement key = <<EOT ([#;] comment)
        block_re = '(\S+)\s*=\s*<<(\w+).*$'
        
        # match EOT on a line with nothing else but whitespace
        # TODO: match the same word that began the block instead
        block_end_re = '\s*EOT\s*$'
        
        currentComment = None
        
        for i in range(len(self.Lines)):
            self.CurrentLineNumber = i+1
            line = self.Lines[i]
        
            line = line.strip()
            
            # blank line? ignore.
            if len(line) == 0:
                continue
            
            # if it begins with # or ;, it's a comment, ignore.
            
            if line[0] == '#' or line[0] == ';':
                #dprint("matched comment on line[%s]: %s" % (i, line))
                if currentComment == None:
                    currentComment = line
                else:
                    currentComment = "%s\n%s" % (currentComment, line)
                continue
            
            # -- SEARCH FOR BLOCK END --

            if self.BlockMode:
            
                p = re.compile(block_end_re)
                ep = p.match(line)
                
                if ep == None:
                    if self.BlockText != None:
                        self.BlockText = "%s%s\n" % (self.BlockText, line)
                    else:
                        self.BlockText = "%s\n" % self.StripInlineComment(line)
                    continue
                else:
                    #dprint("matched End of Block with line %s" % line)
                    #dprint("block key: %s\nblock_text:\n%s" % (self.BlockKey, self.BlockText))
                    
                    #IniOption(self, section, key, value, IsBlock):
                    op = IniOption(self.CurrentSection, self.BlockKey, self.BlockText, True)
                    op.BlockIdentifier = self.BlockIdentifier
                    
                    if currentComment != None:
                        op.LeadingComment = currentComment
                        currentComment = None
                    
                    if self.KeyOrder.has_key(self.CurrentSection):
                        if self.BlockKey in self.KeyOrder[self.CurrentSection]:
                            fatal(self.MakeError(self.CurrentLineNumber, "duplicate key %s" % self.BlockKey))
                        
                    self.WriteKey(self.CurrentSection, self.BlockKey, op)

                    self.BlockMode = False
                    self.BlockKey = None
                    self.BlockText = None
                    self.BlockModeBeginLine = 0
                    self.BlockIdentifier = None
                    continue
        
            # -- SEARCH FOR SECTION  --
            
            #dprint("trying to match section on line[%s]: %s" % (i, line))
        
            p = re.compile(section_re)

            m = p.match(line)
        
            if m != None:
                self.EnterSection(m.group(1))
                continue
            
            
            # -- SEARCH FOR BLOCK --
            # note: has to be in this order, because the key regular expression is greedy, and allows 
            # for all manner of characters after the equals sign. For simplicity of the regular expression
            # it would successfully match <<EOT and accept it as a valid value. Doing this search first
            # guarantees that the regular expression for key will never be attempted on this matched line.
            
            p = re.compile(block_re)

            #dprint("trying to match block on line[%s]: %s" % (i, line))
            
            m = p.match(line)
        
            if m != None:
                #dprint("matched block")
                self.BlockMode = True
                self.BlockKey = m.group(1)
                self.BlockText = None
                self.BlockModeBeginLine = self.CurrentLineNumber
                self.BlockIdentifier = m.group(2)
                #dprint("Found block identifier %s" % self.BlockIdentifier)
                block_end_re = '\s*%s\s*$' % self.BlockIdentifier
                continue
                
            # -- SEARCH FOR KEY --

            #dprint("trying to match key on line[%s]: %s" % (i, line))
            
            p = re.compile(key_re)
            
            m = p.match(line)
        
            if m != None:
                #dprint("%s = %s" % (m.group(1), m.group(2)))
                #IniOption(self,section, key, value, IsBlock):
                op = IniOption(self.CurrentSection, m.group(1), m.group(2), False)
                
                if currentComment != None:
                    op.LeadingComment = currentComment
                    currentComment = None
                
                if self.KeyOrder.has_key(self.CurrentSection):
                    if m.group(1) in self.KeyOrder[self.CurrentSection]:
                        fatal(self.MakeError(self.CurrentLineNumber, "duplicate key %s" % m.group(1)))

                self.WriteKey(self.CurrentSection, m.group(1), op)

                continue

            
            # ERROR ERROR ERROR ERROR ERROR
            self.HandleError(line)
        
        if self.BlockMode:
            error = "Runaway <<%s block. No ending %s found" % (self.BlockIdentifier, self.BlockIdentifier)
            fatal(self.MakeError(self.BlockModeBeginLine, error))
