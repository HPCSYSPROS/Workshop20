#!/usr/bin/python

import sys, os, re

############################################
## APP SPECIFIC DEFINES ##
#
# Usage			when (-h, --help) is used, this is the first line showing how to call the app
# optionGroups	which option groups from etc/options to import. common is always imported. 

App = "mdb-merge.py"
Usage = "usage: %prog [options] command simulationname"
optionGroups = []
Purpose = "merge scriptFile run into each mdb entry"

############################################

############################################
## INIT ##

if os.path.basename(sys.argv[0]) == App:
	rp = os.path.realpath(__file__)
	
	paths = rp.split(os.sep)
	
	paths.pop()
	paths.pop()
	
	global BASE_PATH
	BASE_PATH = os.sep.join(paths)
	sys.path.append(BASE_PATH)
	
	#print "BASE_PATH: %s" % BASE_PATH
	
from lib import *

def getCmd(scriptFile):
	global BASE_PATH
	
	filename = "%s/etc/scriptfiles/%s" % (BASE_PATH, scriptFile)
	
	full_command = str()
	try:
		fptr = open(filename, "r")
	except:
		#print "could not open %s for reading" % filename
		return full_command
	
	contents = fptr.read()
	
	lines = contents.split("\n")
	
	cmd = str()
	mpichdir = str()
	
	mpilist = ['MPICHDIR', 'MPICH_DIR', 'OPENMPI_DIR']
	
	for mpidef in mpilist:
	
		if contents.count(mpidef) > 0:

			for line in lines:
				line = line.strip()
				
				try:
					if line.index("time") == 0:
						if line.count("psrun") > 0:
							continue
						cmd = line
						break
				except:
					pass
				
				try:
					if line.index("${%s}" % mpidef) == 0:
						if line.count('run') > 0:
							cmd = line
							break
				except:
					pass
			
			if len(mpichdir) > 0:
				full_command = "%s && %s" % (mpichdir, cmd)
			else:
				full_command = cmd
			break
	
	alt_commands = ['/poe', 'srun', 'aprun', 'mpirun', 'ibrun', 'mpiexec']
	
	for alt in alt_commands:
		if contents.count(alt) > 0:
			for line in lines:
				line = line.strip()

				if line.startswith("#"):
					continue
					
				try:
					if line.count(alt) > 0:
						cmd = line
				except:
					pass
			
			if len(mpichdir) > 0:
				full_command = "%s && %s" % (mpichdir, cmd)
			else:
				full_command = cmd
				
	#print "returning cmd: %s" % (cmd)
	return full_command


def fixScriptFiles():
	global BASE_PATH
	
	scriptdir = "%s/etc/scriptfiles" % BASE_PATH
	cmds = dict()
	
	for sfile in os.listdir(scriptdir):
	
		if sfile.count(".sh") == 0:
			continue
			
		cmd = getCmd(sfile).strip()
		
		fptr = open("%s/%s" % (scriptdir, sfile), "r")
		
		contents = fptr.read()
		
		simcmd = "@SIMFACTORY@ run @SIMULATION_NAME@ --restart-id=@RESTART_ID@ @FROM_RESTART_COMMAND@"
		
		print "replacing: %s" % cmd
		contents = contents.replace(cmd, simcmd)
		
		fptr.close()
		
		fptr = open("%s/etc/submitscripts/%s" % (BASE_PATH, sfile), "w")
		
		fptr.write(contents)
		fptr.write("\n")
		
		fptr.close()
			
def fixMdbEntries():

	global BASE_PATH
	
	cmds = dict()
	
	for mdb_file in os.listdir("%s/etc/mdb/" % BASE_PATH):
		
		print "mdb_file: %s" % mdb_file
		
		if mdb_file.count(".ini") == 0:
			continue
			
		parser = pyini.IniParser("%s/etc/mdb/%s" % (BASE_PATH, mdb_file))
		machine = parser.GetSections()[0]
		
		print "machine: %s" % machine
		
		scriptfile = parser.GetOption(machine, 'submitscript').Value
		cmd = getCmd(scriptfile).strip()
		cmds[machine] = cmd
	
		if cmd.startswith("#"):
			continue
			
		fptr = open("%s/etc/mdb/%s" % (BASE_PATH, mdb_file), "r")
	
		lines = fptr.readlines()
	
		nary = list()
		offset = 0
		for i in range(len(lines)):
			line = lines[i].strip()
			cmd = cmds[machine]
			
			#strip existing mpirun.
			if line.startswith("mpirun="):
				nary.append("mpirun          = %s\n" % cmd)
				continue
			else:
				nary.append(line)

		fptr.close()
		
		
		
		fptr = open("%s/etc/mdb/%s" % (BASE_PATH, mdb_file), "w")
		
		fptr.write("\n".join(nary))
		fptr.write("\n")
		fptr.close()
		
def main():
	global BASE_PATH
	
	#fixMdbEntries()
	fixScriptFiles()

main()