# simsubs.py -- A library for setting/unsetting defines and performing substitutions of said defines
# @YAHOO@
#
# Michael Thomas <mthomas@cct.lsu.edu>
# Center for Computation & Technology
# Louisiana State University 
#

import re
import sys
import simenv

from libutil import *

class DefineDatabase:
    def __init__(self):

        self.PresetMacros = simenv.OptionsManager.Macros
        self.Patterns = dict()
        
        self.Reset()
    
    def Reset(self):
        self.defines = dict()
        
        #load in our preset (command-line defined) macros.
        for key in self.PresetMacros.keys():
            #dprint("loading predefined macro %s with value %s" % (key, self.PresetMacros[key]))
            self.Set(key, self.PresetMacros[key])
    
        # TODO: ES 2012-07-31: Don't Defines and Substitutions need to
        # be compiled as well?
        self.CompileReplacements()
        self.CompileAppends()
        
    def Set(self, define, value):
        
        # if we attempt to set something that's preset, just ignore it.
#        if self.PresetMacros.has_key(define):
#            return
        
        self.defines[define] = value
    
    def Unset(self, define):
        try:
            del self.defines[define]
        except KeyError:
            return
    
    def Has(self, define):
        return self.defines.has_key(define)
        
    def Get(self, define):
        if self.Has(define):
            return self.defines[define]
        
        return None

    def Sub(self, define, ss):
        if self.Has(define):
            macro = "@%s@" % define
            return self.PerformSub(ss, macro, self.defines[define])
        
        return ss
    
    def PerformSub(self, ss, macro, value):
        if value is not str:
            value = str(value)

        return ss.replace(macro, value)

    def ParseEnvCommands(self, ss):
        pattern = "@ENV\((.+)\)@[ \t]+"
        
        env = []
        
        rx = re.compile(pattern)
        
        mm = rx.findall(ss)
        
        for m in mm:
            env.append(m)
        
        for m in env:
            f_m = "@ENV(%s)@" % m
            ss = self.PerformSub(ss, f_m, '')
            
        return (ss, env)

    def PerformRegexSubstitutions(self, ss):
        slist = simenv.OptionsManager.Substitutions
        
        if simenv.ConfigurationDatabase.HasConfigOption("substitutions"):
            slist2 = simenv.ConfigurationDatabase.GetConfigOption("substitutions")
            slist.extend(slist2)
        
        if len(slist) == 0:
            return ss
        
        for rx_pair in slist:
            rx_pair[1] = rx_pair[1].replace("@1@", r"\1")

            rx = re.compile(rx_pair[0], re.MULTILINE)
            ss = rx.replace(rx_pair[1], ss)
        
        return ss
    
    # Add additional define
    #    @VAR@ = VALUE
    def AddDefine(self, var, value):
        if var != var.strip():
            fatal("illegal variable name '%s'" % var)
                
        pattern = r"@%s@" % var
                
        if pattern in self.Patterns:
            fatal("Error: pattern '%s' exists already" % pattern)

        self.Patterns[pattern] = value
                
    # Add additional substitution
    #    PATTERN = VALUE
    def AddSubstitution(self, var, value):
        pattern = var
        
        if pattern in self.Patterns:
            fatal("pattern '%s' exists already" % pattern)
    
        self.Patterns[pattern] = str(value).replace("@1@", r"\1")

    # Add additional replacement
    #    KEY = NEWVALUE
    def AddReplacement(self, var, value):
        if var != var.strip():
            fatal("illegal variable name '%s'" % var)
            
        var1 = var.replace("$", r"\$")
        pattern = r"^[ \t]*%s[ \t]*=.*$" % var1
        
        if pattern in self.Patterns:
            fatal("pattern '%s' exists already" % pattern)

        self.Patterns[pattern] = "=SIMFACTORY-PROTECTED=%s = %s" % (var, value)
    
        # Add additional append
        #   KEY = NEWVALUE OLDVALUE
    def AddAppend(self, var, value):
        if var != var.strip():
            fatal("illegal variable name '%s'" % var)
        
        pattern = r'^[ \t]*%s[ \t]*=[ \t]*(.*?)[ \t]*$' % var
        
        if pattern in self.Patterns:
            fatal("pattern '%s' exists already" % pattern)
            
        self.Patterns[pattern] = r"=SIMFACTORY-PROTECTED=%s = %s \1" % (var, value)
        
        #def CompileDefines
        #def CompileSubstitutions
        
    def CompileReplacements(self):
        
        rlist = simenv.OptionsManager.Replacements
        
        if simenv.ConfigurationDatabase.HasConfigOption("replacements"):
            rlist2 = simenv.ConfigurationDatabase.GetConfigOption("replacements")
            rlist.extend(rlist2)
        
        for r in rlist:
            (var, value) = r.split('=', 1)
            self.AddReplacement(var, value)
    
    def CompileAppends(self):
        
        rlist = simenv.OptionsManager.Appends
        
        if simenv.ConfigurationDatabase.HasConfigOption("appends"):
            rlist2 = simenv.ConfigurationDatabase.GetConfigOption("appends")
            rlist.extend(rlist2)
        
        for r in rlist:
            (var, value) = r.split('=', 1)
            self.AddAppend(var, value)
    
    def SubAll(self, ss):
    
        while True:
            oldss = ss
            
            for macro in self.defines.keys():
                f_m = "@%s@" % macro
                ss = self.PerformSub(ss, f_m, self.Get(macro))
        
            pattern = "@\(([^@]+)\)@"
            
            rx = re.compile(pattern)
            
            mm = rx.findall(ss)
            
            for m in mm:
                local_dict = dict()
                local_dict['ifthen'] = self.ifthen
                
                #translate ? into ifthen
                if m.count("?") >= 1:
                    if m.count("?") > 2:
                        fatal("cannot evaluate complex ternary operation \"%s\", only one ? allowed." % m)
                        
                    m_new = self.TransformTernary(m)
                else:
                    m_new = m
                    
                output = eval(m_new, local_dict)
    
                f_m = "@(%s)@" % m
                ss = self.PerformSub(ss, f_m, output)
        
            pattern = "@ENV\(([^)]+)\)@"
            rx = re.compile(pattern)
            mm = rx.findall(ss)
            for m in mm:
                try:
                    val = os.environ[m]
                except KeyError:
                    fatal('Environment variable "%s" not set' % m)
                f_m = "@ENV(%s)@" % m
                ss = self.PerformSub(ss, f_m, val)
            
            ss = self.PerformRegexSubstitutions(ss)
            ss = self.PerformPatternSubstitutions(ss)
        
            if ss == oldss:
                ss = ss.replace("=SIMFACTORY-PROTECTED=", "")
                return ss
    
    def PerformPatternSubstitutions(self, m):
        
        for pattern in self.Patterns.keys():
            #display("have replacement pattern: %s, replacement: %s" % (pattern, self.Patterns[pattern]))
            
            rx = re.compile(pattern, re.MULTILINE)
            m = rx.sub(self.Patterns[pattern], m)
        
        return m
        
    # TODO: remove this function, update submit scripts instead
    def TransformTernary(self, m):
        
        # statement ? value : value 
        
        #dprint("ternary statement: %s" % m)
        sp = m.split("?")
        
        statement = sp[0].strip()
        replacements = dict()
        
        replacements['ne'] = "!="
        replacements['eq'] = "=="
        replacements['gt'] = ">"
        replacements['lt'] = "<"
        replacements['gte'] = ">="
        replacements['lte'] = "<="
        
        for key in replacements.keys():
            statement = statement.replace(key, replacements[key])
        
        
        # now, we need to go character by character here to find the first : not in quotes.
        count = 1
        quoteChar = None
        
        for c in sp[1]:
            if c in ["'", '"']:
                if quoteChar == None:
                    quoteChar = c
                else:
                    if c == quoteChar:
                        quoteChar = None
            
            if c == ":" and quoteChar == None:
                splitPos = count
                break
            
            count = count + 1
                

        tt = sp[1][1:splitPos-1].strip()
        ff = sp[1][splitPos:].strip()
        
        cmd = "ifthen(%s, %s, %s)" % (statement, tt, ff)
        
        #dprint("translated ternary to: %s" % cmd)
        
        return cmd
        
    def ifthen(self, expr, tval, fval):
        try:
            if expr:
                return tval
            else:
                return fval
        except:
            return fval
