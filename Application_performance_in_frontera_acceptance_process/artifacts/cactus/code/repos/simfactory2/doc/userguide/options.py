#!/usr/bin/env python
import sim
import simenv
import sys
import pyini

ini = pyini.IniParser("../../etc/options/common.ini")

essential_options = ["remote"]
essential_options.sort()

options = filter(lambda c : c not in essential_options, list(ini.GetSections()))
options.sort()

o = open("_auto/options/common-essential.txt","w")
#print >>o, ".. list-table::"
#print >>o, ""

for section in essential_options:
    s = ini.GetSectionAsDict(section)
    # print >>o, "   * - "+s["long"].value
    # print >>o, "     - "+s["desc"].value
    print >>o, ".. cmdoption:: %s" % (s["long"].value),
    if "argformat" in s:
        print >>o,s["argformat"].value
    else:
        print >>o
    print >>o
    print >>o, "  %s" % s["desc"].value
    print >>o

o = open("_auto/options/common-rest.txt","w")
#print >>o, ".. list-table::"
#print >>o, ""

for section in options:
    s = ini.GetSectionAsDict(section)
    # print >>o, "   * - "+s["long"].value
    # print >>o, "     - "+s["desc"].value
    print >>o, ".. cmdoption:: %s" % (s["long"].value),
    if "argformat" in s:
        print >>o,s["argformat"].value
    else:
        print >>o
    print >>o
    print >>o, "  %s" % s["desc"].value
    print >>o

o.close()
