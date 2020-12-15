


Simulation Factory 2.0:

Michael Thomas <mthomas@cct.lsu.edu>
Erik Schnetter <schnetter@cct.lsu.edu>

git: https://bitbucket.org/simfactory/simfactory2.git


1. What's new

	1. Simfactory has been ported from Perl to Python, and requires at minimum Python >= 2.3. 
	2. The databases, udb/cdb/mdb.pm have been changed to an INI-style database and are now located
	   in the etc/ folder. The mdb has been split into individual files contained within the etc/mdb folder.
	3. To override any key in an mdb entry, simply provide that key in the [default] section within the udb.ini file.
	   To override a specific key within one mdb entry, define a [machinename] section and place the key underneath it.
	4. Simfactory has been split into four seperate applications, which are located in the bin/ folder. The application
	   "sim" is a wrapper for these four seperate applications and can be used instead of calling them directly.
	5. Instead of 'sim help', it is now 'sim --help', or 'bin/sim-* --help'
	6. the scriptfiles folder is for antique purposes only. Scriptfiles have been split into two seperate scripts now. The first
	   is used for submission, and those are contained in the etc/submitscripts folder. The second is used by the run command
	   to execute the simulation, and those are in etc/runscripts.
	   
	7. There are several new mdb keys: 
		
		environcmd: prepend this command before *any* local command gets executed.
		submitscript: the path to the new submit script
		runscript: the path to the new run script
		precmd: reserved for future use
		postcmd: reserved for future use
		holdingpattern: if the job is on hold because it is dependant on a non-finished job, this pattern
		                is used to detect that state
		interactivecmd: execute this command to enter into an interactive session
		
2. Submitscripts
	
	Simfactory now only uses scripts for submission, not execution, and what is contained in the submit
   	script has changed. An example submit script is inside the etc/submitscripts folder and is called
    'pbs-normal.sh'. It contains PBS commands, and some example directives that the submitscript supports.
    Since python 2.3 does not have the equivelent to a ? operator, a directive had to be added to support
    this construct. The directive is:
   
    @(ifthen('statement', 'true', 'false'))@
    
    However, limited ? support has been implemented. Currently it only supports one level of ? operation, eg
    'statement ? true : false', and not something like statement ? statement ? : true : false : false. There is also
    support for translating the Perl operators ne, gt, eq, lt, gte, lte inside of a ? statement.
   
    At the end of the submit script is the new simfactory run command, and looks like this:
    
    @SIMFACTORY@ run @SIMULATION_NAME@ --restart-id=@RESTART_ID@ @FROM_RESTART_COMMAND@
    
    This line should be identical in *every* submitscript. There should be no reason at all to change it.

3. Runscripts
	
	All execution logic, such as loading required modules and setting environment variables and the like have been moved into
	the run script. This script is responsible for calling the right command dispatch system (such as mpi). The current runscripts
	are the bottom half of the original scriptfile, containing everything except the #PBS/etc comments that are used to interface
	with the queueing system. 

4. Usage
	
	Simfactory provides two ways of executing commands. First, there is the (mostly) compatible with the previous version 'sim' command, 
	which resides in the base directory of simfactory. This command serves as a junction of the four individiual commands listed below.
	To get help/usage information, ussue the command 'sim --help', or 'sim help'. Either works.
	
	Simfactory has been divided into four separate commands:
	
	1. bin/sim-sync.py: equiv to 'sim sync'
	2. bin/sim-build.py: equiv to 'sim build'
	3. bin/sim-manage.py: equiv to 'sim create/submit', but with the addition of interactive/run. 
	4. bin/sim-info.py: equiv to 'sim list-*, sim print-*'
	
	All commands can be executed remotely using the --remote=<machine> argument.

	You can override mdb entries and defines via the commandline. 
	To override the mdb, you use --mdbkey <keyname>
	<value>, and to override/set a new define, you use --define
	<define> <value>. Replacements, substitutions, and appends(attaches) have been
	carried over as well and function the same as they did in the Perl version.

	examples:

	--mdbkey "@EXECUTABLE@ -L 3 @PARFILE@"
	--define EXECUTABLE exe/cactus_somesimname

4. Limitations
	
	When using 'sim run', or 'sim-manage run', and you are still using mpi, a nodelist is generated that just says localhost. It currently
	does not repeat localhost correctly up to the number of processors you are intending to use. I will be fixing this as I learn
	more about what the correct procedure is. 
