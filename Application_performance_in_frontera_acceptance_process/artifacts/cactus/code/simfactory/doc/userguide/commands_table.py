#!/usr/bin/env python
import sim
import simenv
import sys
import pyini
import os

def genTable(t):
    colss = map(len, t)
    maxcols = max(colss)
    rows = len(t)

    tt = zip(*t)
    lengths = map(lambda els : max(map(len,els)), tt)

    s = ""

    sep = ""
    for j in range(0,len(t[0])):
        sep += ("+"+"-"*(lengths[j]+2))
    sep += "+\n"

    s += sep

    for i in range(0,len(t)):
        for j in range(0,len(t[i])):
            s += "|"
            s += " "+t[i][j].ljust(lengths[j])+" "
        s += "|\n"
        s += sep

    return s

ini = pyini.IniParser()

for f in os.listdir("../../etc/options"):
    if f != ".svn":
        ini.UpdateFromIni("../../etc/options/"+f,True)

essential_commands = ["build", "create", "submit", "create-submit", "purge", "setup", "stop", "sync"]
essential_commands.sort()

cmds = filter(lambda c : c not in essential_commands, list(sim.known_commands))
cmds.sort()

def command_table(cmds):
    table = []
    for cmd in cmds:
        binary = sim.known_commands[cmd]
        if sys.modules[binary].usage_strings.has_key(cmd):
            usage = sys.modules[binary].usage_strings[cmd]
        else:
            usage = ""
    
        table.append([":doc:`/_auto/commands/"+cmd+"`",usage])
    
        f = open("_auto/commands/%s.rst" % cmd, "w")
        print >>f, ".. _command_"+cmd+":"
        print >>f
        print >>f, cmd
        print >>f, "================================================"
        print >>f
        print >>f, "*"+usage+"*"
        print >>f
        if simenv.CommandDoc.has_key(cmd):
            print >>f, simenv.CommandDoc[cmd]
    
        if simenv.CommandOptions.has_key(cmd):
            print >>f, "Options"
            print >>f, "-------"
            # print >>f, ".. list-table::"
            # print >>f, ""
            opts_p = simenv.CommandOptions[cmd] 
            if type(opts_p) == list:
                opts = opts_p
            elif type(opts_p) == str:
                opts = pyini.IniParser("../../etc/options/"+opts_p).GetSections()
            else:
                print "Command options for %s of bad type %s" %(cmd,str(type(opts_p)))
                sys.exit(1)
    
            print >>f, ".. program:: sim %s" % cmd
            print >>f

            for opt in opts:
                s = ini.GetSectionAsDict(opt)
                # print >>f, "   * - "+s["long"].value
                # print >>f, "     - "+s["desc"].value
                print >>f, ".. cmdoption:: %s" % (s["long"].value),
                if "argformat" in s:
                    print >>f,s["argformat"].value
                else:
                    print >>f
                print >>f
                print >>f, "  %s" % s["desc"].value
                print >>f
    
        f.close
    return table

o = open("_auto/commands_table_essential.txt","w")
print >>o, genTable(command_table(essential_commands))
o.close()

o = open("_auto/commands_table_rest.txt","w")
print >>o, genTable(command_table(cmds))
o.close()


