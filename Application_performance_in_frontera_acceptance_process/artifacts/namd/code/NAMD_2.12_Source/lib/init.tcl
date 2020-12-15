# init.tcl --
#
# Default system startup file for Tcl-based applications.  Defines
# "unknown" procedure and auto-load facilities.
#
# RCS: @(#) $Id: init.tcl,v 1.1 2005/07/21 20:45:53 jim Exp $
#
# Copyright (c) 1991-1993 The Regents of the University of California.
# Copyright (c) 1994-1996 Sun Microsystems, Inc.
# Copyright (c) 1998-1999 Scriptics Corporation.
#
#
# This software is copyrighted by the Regents of the University of
# California, Sun Microsystems, Inc., Scriptics Corporation,
# and other parties.  The following terms apply to all files associated
# with the software unless explicitly disclaimed in individual files.
#
# The authors hereby grant permission to use, copy, modify, distribute,
# and license this software and its documentation for any purpose, provided
# that existing copyright notices are retained in all copies and that this
# notice is included verbatim in any distributions. No written agreement,
# license, or royalty fee is required for any of the authorized uses.
# Modifications to this software may be copyrighted by their authors
# and need not follow the licensing terms described here, provided that
# the new terms are clearly indicated on the first page of each file where
# they apply.
#
# IN NO EVENT SHALL THE AUTHORS OR DISTRIBUTORS BE LIABLE TO ANY PARTY
# FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
# ARISING OUT OF THE USE OF THIS SOFTWARE, ITS DOCUMENTATION, OR ANY
# DERIVATIVES THEREOF, EVEN IF THE AUTHORS HAVE BEEN ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# THE AUTHORS AND DISTRIBUTORS SPECIFICALLY DISCLAIM ANY WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE, AND NON-INFRINGEMENT.  THIS SOFTWARE
# IS PROVIDED ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE
# NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
# MODIFICATIONS.
#
# GOVERNMENT USE: If you are acquiring this software on behalf of the
# U.S. government, the Government shall have only "Restricted Rights"
# in the software and related documentation as defined in the Federal 
# Acquisition Regulations (FARs) in Clause 52.227.19 (c) (2).  If you
# are acquiring the software on behalf of the Department of Defense, the
# software shall be classified as "Commercial Computer Software" and the
# Government shall have only "Restricted Rights" as defined in Clause
# 252.227-7013 (c) (1) of DFARs.  Notwithstanding the foregoing, the
# authors grant the U.S. Government and others acting in its behalf
# permission to use and distribute the software in accordance with the
# terms specified in this license. 

if {[info commands package] == ""} {
    error "version mismatch: library\nscripts expect Tcl version 7.5b1 or later but the loaded version is\nonly [info patchlevel]"
}
package require -exact Tcl 8.3

# Compute the auto path to use in this interpreter.

if {[file pathtype [info script]] == "absolute"} {
    lappend auto_path [file dirname [info script]]
} {
    lappend auto_path [file dirname [file join [pwd] [info script]]]
}

# Windows specific end of initialization

if {(![interp issafe]) && [string equal $tcl_platform(platform) "windows"]} {
    namespace eval tcl {
	proc envTraceProc {lo n1 n2 op} {
	    set x $::env($n2)
	    set ::env($lo) $x
	    set ::env([string toupper $lo]) $x
	}
    }
    foreach p [array names env] {
	set u [string toupper $p]
	if {[string compare $u $p]} {
	    switch -- $u {
		COMSPEC -
		PATH {
		    if {![info exists env($u)]} {
			set env($u) $env($p)
		    }
		    trace variable env($p) w [list tcl::envTraceProc $p]
		    trace variable env($u) w [list tcl::envTraceProc $p]
		}
	    }
	}
    }
    if {[info exists p]} {
	unset p
    }
    if {[info exists u]} {
	unset u
    }
    if {![info exists env(COMSPEC)]} {
	if {[string equal $tcl_platform(os) "Windows NT"]} {
	    set env(COMSPEC) cmd.exe
	} else {
	    set env(COMSPEC) command.com
	}
    }
}

# Setup the unknown package handler

package unknown tclPkgUnknown

# Conditionalize for presence of exec.

if {[llength [info commands exec]] == 0} {

    # Some machines, such as the Macintosh, do not have exec. Also, on all
    # platforms, safe interpreters do not have exec.

    set auto_noexec 1
}
set errorCode ""
set errorInfo ""

# Define a log command (which can be overwitten to log errors
# differently, specially when stderr is not available)

if {[llength [info commands tclLog]] == 0} {
    proc tclLog {string} {
	catch {puts stderr $string}
    }
}

# auto_load --
# Checks a collection of library directories to see if a procedure
# is defined in one of them.  If so, it sources the appropriate
# library file to create the procedure.  Returns 1 if it successfully
# loaded the procedure, 0 otherwise.
#
# Arguments: 
# cmd -			Name of the command to find and load.
# namespace (optional)  The namespace where the command is being used - must be
#                       a canonical namespace as returned [namespace current]
#                       for instance. If not given, namespace current is used.

proc auto_load {cmd {namespace {}}} {
    global auto_index auto_oldpath auto_path

    if {[string length $namespace] == 0} {
	set namespace [uplevel 1 [list ::namespace current]]
    }
    set nameList [auto_qualify $cmd $namespace]
    # workaround non canonical auto_index entries that might be around
    # from older auto_mkindex versions
    lappend nameList $cmd
    foreach name $nameList {
	if {[info exists auto_index($name)]} {
	    uplevel #0 $auto_index($name)
	    return [expr {[info commands $name] != ""}]
	}
    }
    if {![info exists auto_path]} {
	return 0
    }

    if {![auto_load_index]} {
	return 0
    }
    foreach name $nameList {
	if {[info exists auto_index($name)]} {
	    uplevel #0 $auto_index($name)
	    # There's a couple of ways to look for a command of a given
	    # name.  One is to use
	    #    info commands $name
	    # Unfortunately, if the name has glob-magic chars in it like *
	    # or [], it may not match.  For our purposes here, a better
	    # route is to use 
	    #    namespace which -command $name
	    if { ![string equal [namespace which -command $name] ""] } {
		return 1
	    }
	}
    }
    return 0
}

# auto_load_index --
# Loads the contents of tclIndex files on the auto_path directory
# list.  This is usually invoked within auto_load to load the index
# of available commands.  Returns 1 if the index is loaded, and 0 if
# the index is already loaded and up to date.
#
# Arguments: 
# None.

proc auto_load_index {} {
    global auto_index auto_oldpath auto_path errorInfo errorCode

    if {[info exists auto_oldpath] && \
	    [string equal $auto_oldpath $auto_path]} {
	return 0
    }
    set auto_oldpath $auto_path

    # Check if we are a safe interpreter. In that case, we support only
    # newer format tclIndex files.

    set issafe [interp issafe]
    for {set i [expr {[llength $auto_path] - 1}]} {$i >= 0} {incr i -1} {
	set dir [lindex $auto_path $i]
	set f ""
	if {$issafe} {
	    catch {source [file join $dir tclIndex]}
	} elseif {[catch {set f [open [file join $dir tclIndex]]}]} {
	    continue
	} else {
	    set error [catch {
		set id [gets $f]
		if {[string equal $id \
			"# Tcl autoload index file, version 2.0"]} {
		    eval [read $f]
		} elseif {[string equal $id "# Tcl autoload index file: each line identifies a Tcl"]} {
		    while {[gets $f line] >= 0} {
			if {[string equal [string index $line 0] "#"] \
				|| ([llength $line] != 2)} {
			    continue
			}
			set name [lindex $line 0]
			set auto_index($name) \
				"source [file join $dir [lindex $line 1]]"
		    }
		} else {
		    error "[file join $dir tclIndex] isn't a proper Tcl index file"
		}
	    } msg]
	    if {[string compare $f ""]} {
		close $f
	    }
	    if {$error} {
		error $msg $errorInfo $errorCode
	    }
	}
    }
    return 1
}

# auto_qualify --
#
# Compute a fully qualified names list for use in the auto_index array.
# For historical reasons, commands in the global namespace do not have leading
# :: in the index key. The list has two elements when the command name is
# relative (no leading ::) and the namespace is not the global one. Otherwise
# only one name is returned (and searched in the auto_index).
#
# Arguments -
# cmd		The command name. Can be any name accepted for command
#               invocations (Like "foo::::bar").
# namespace	The namespace where the command is being used - must be
#               a canonical namespace as returned by [namespace current]
#               for instance.

proc auto_qualify {cmd namespace} {

    # count separators and clean them up
    # (making sure that foo:::::bar will be treated as foo::bar)
    set n [regsub -all {::+} $cmd :: cmd]

    # Ignore namespace if the name starts with ::
    # Handle special case of only leading ::

    # Before each return case we give an example of which category it is
    # with the following form :
    # ( inputCmd, inputNameSpace) -> output

    if {[regexp {^::(.*)$} $cmd x tail]} {
	if {$n > 1} {
	    # ( ::foo::bar , * ) -> ::foo::bar
	    return [list $cmd]
	} else {
	    # ( ::global , * ) -> global
	    return [list $tail]
	}
    }
    
    # Potentially returning 2 elements to try  :
    # (if the current namespace is not the global one)

    if {$n == 0} {
	if {[string equal $namespace ::]} {
	    # ( nocolons , :: ) -> nocolons
	    return [list $cmd]
	} else {
	    # ( nocolons , ::sub ) -> ::sub::nocolons nocolons
	    return [list ${namespace}::$cmd $cmd]
	}
    } elseif {[string equal $namespace ::]} {
	#  ( foo::bar , :: ) -> ::foo::bar
	return [list ::$cmd]
    } else {
	# ( foo::bar , ::sub ) -> ::sub::foo::bar ::foo::bar
	return [list ${namespace}::$cmd ::$cmd]
    }
}

# auto_import --
#
# Invoked during "namespace import" to make see if the imported commands
# reside in an autoloaded library.  If so, the commands are loaded so
# that they will be available for the import links.  If not, then this
# procedure does nothing.
#
# Arguments -
# pattern	The pattern of commands being imported (like "foo::*")
#               a canonical namespace as returned by [namespace current]

proc auto_import {pattern} {
    global auto_index

    # If no namespace is specified, this will be an error case

    if {![string match *::* $pattern]} {
	return
    }

    set ns [uplevel 1 [list ::namespace current]]
    set patternList [auto_qualify $pattern $ns]

    auto_load_index

    foreach pattern $patternList {
        foreach name [array names auto_index] {
            if {[string match $pattern $name] && \
		    [string equal "" [info commands $name]]} {
                uplevel #0 $auto_index($name)
            }
        }
    }
}

# auto_execok --
#
# Returns string that indicates name of program to execute if 
# name corresponds to a shell builtin or an executable in the
# Windows search path, or "" otherwise.  Builds an associative 
# array auto_execs that caches information about previous checks, 
# for speed.
#
# Arguments: 
# name -			Name of a command.

if {[string equal windows $tcl_platform(platform)]} {
# Windows version.
#
# Note that info executable doesn't work under Windows, so we have to
# look for files with .exe, .com, or .bat extensions.  Also, the path
# may be in the Path or PATH environment variables, and path
# components are separated with semicolons, not colons as under Unix.
#
proc auto_execok name {
    global auto_execs env tcl_platform

    if {[info exists auto_execs($name)]} {
	return $auto_execs($name)
    }
    set auto_execs($name) ""

    set shellBuiltins [list cls copy date del erase dir echo mkdir \
	    md rename ren rmdir rd time type ver vol]
    if {[string equal $tcl_platform(os) "Windows NT"]} {
	# NT includes the 'start' built-in
	lappend shellBuiltins "start"
    }
    if {[info exists env(PATHEXT)]} {
	# Add an initial ; to have the {} extension check first.
	set execExtensions [split ";$env(PATHEXT)" ";"]
    } else {
	set execExtensions [list {} .com .exe .bat]
    }

    if {[lsearch -exact $shellBuiltins $name] != -1} {
	return [set auto_execs($name) [list $env(COMSPEC) /c $name]]
    }

    if {[llength [file split $name]] != 1} {
	foreach ext $execExtensions {
	    set file ${name}${ext}
	    if {[file exists $file] && ![file isdirectory $file]} {
		return [set auto_execs($name) [list $file]]
	    }
	}
	return ""
    }

    set path "[file dirname [info nameof]];.;"
    if {[info exists env(WINDIR)]} {
	set windir $env(WINDIR) 
    }
    if {[info exists windir]} {
	if {[string equal $tcl_platform(os) "Windows NT"]} {
	    append path "$windir/system32;"
	}
	append path "$windir/system;$windir;"
    }

    foreach var {PATH Path path} {
	if {[info exists env($var)]} {
	    append path ";$env($var)"
	}
    }

    foreach dir [split $path {;}] {
	# Skip already checked directories
	if {[info exists checked($dir)] || [string equal {} $dir]} { continue }
	set checked($dir) {}
	foreach ext $execExtensions {
	    set file [file join $dir ${name}${ext}]
	    if {[file exists $file] && ![file isdirectory $file]} {
		return [set auto_execs($name) [list $file]]
	    }
	}
    }
    return ""
}

} else {
# Unix version.
#
proc auto_execok name {
    global auto_execs env

    if {[info exists auto_execs($name)]} {
	return $auto_execs($name)
    }
    set auto_execs($name) ""
    if {[llength [file split $name]] != 1} {
	if {[file executable $name] && ![file isdirectory $name]} {
	    set auto_execs($name) [list $name]
	}
	return $auto_execs($name)
    }
    foreach dir [split $env(PATH) :] {
	if {[string equal $dir ""]} {
	    set dir .
	}
	set file [file join $dir $name]
	if {[file executable $file] && ![file isdirectory $file]} {
	    set auto_execs($name) [list $file]
	    return $auto_execs($name)
	}
    }
    return ""
}

}

# package.tcl --
#
# utility procs formerly in init.tcl which can be loaded on demand
# for package management.
#
# RCS: @(#) $Id: init.tcl,v 1.1 2005/07/21 20:45:53 jim Exp $
#
# Copyright (c) 1991-1993 The Regents of the University of California.
# Copyright (c) 1994-1998 Sun Microsystems, Inc.
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

# Create the package namespace
namespace eval ::pkg {
}

# pkg_compareExtension --
#
#  Used internally by pkg_mkIndex to compare the extension of a file to
#  a given extension. On Windows, it uses a case-insensitive comparison
#  because the file system can be file insensitive.
#
# Arguments:
#  fileName	name of a file whose extension is compared
#  ext		(optional) The extension to compare against; you must
#		provide the starting dot.
#		Defaults to [info sharedlibextension]
#
# Results:
#  Returns 1 if the extension matches, 0 otherwise

proc pkg_compareExtension { fileName {ext {}} } {
    global tcl_platform
    if {![string length $ext]} {set ext [info sharedlibextension]}
    if {[string equal $tcl_platform(platform) "windows"]} {
        return [string equal -nocase [file extension $fileName] $ext]
    } else {
        # Some unices add trailing numbers after the .so, so
        # we could have something like '.so.1.2'.
        set root $fileName
        while {1} {
            set currExt [file extension $root]
            if {[string equal $currExt $ext]} {
                return 1
            } 

	    # The current extension does not match; if it is not a numeric
	    # value, quit, as we are only looking to ignore version number
	    # extensions.  Otherwise we might return 1 in this case:
	    #		pkg_compareExtension foo.so.bar .so
	    # which should not match.

	    if { ![string is integer -strict [string range $currExt 1 end]] } {
		return 0
	    }
            set root [file rootname $root]
	}
    }
}

# pkg_mkIndex --
# This procedure creates a package index in a given directory.  The
# package index consists of a "pkgIndex.tcl" file whose contents are
# a Tcl script that sets up package information with "package require"
# commands.  The commands describe all of the packages defined by the
# files given as arguments.
#
# Arguments:
# -direct		(optional) If this flag is present, the generated
#			code in pkgMkIndex.tcl will cause the package to be
#			loaded when "package require" is executed, rather
#			than lazily when the first reference to an exported
#			procedure in the package is made.
# -verbose		(optional) Verbose output; the name of each file that
#			was successfully rocessed is printed out. Additionally,
#			if processing of a file failed a message is printed.
# -load pat		(optional) Preload any packages whose names match
#			the pattern.  Used to handle DLLs that depend on
#			other packages during their Init procedure.
# dir -			Name of the directory in which to create the index.
# args -		Any number of additional arguments, each giving
#			a glob pattern that matches the names of one or
#			more shared libraries or Tcl script files in
#			dir.

proc pkg_mkIndex {args} {
    global errorCode errorInfo
    set usage {"pkg_mkIndex ?-direct? ?-verbose? ?-load pattern? ?--? dir ?pattern ...?"};

    set argCount [llength $args]
    if {$argCount < 1} {
	return -code error "wrong # args: should be\n$usage"
    }

    set more ""
    set direct 1
    set doVerbose 0
    set loadPat ""
    for {set idx 0} {$idx < $argCount} {incr idx} {
	set flag [lindex $args $idx]
	switch -glob -- $flag {
	    -- {
		# done with the flags
		incr idx
		break
	    }
	    -verbose {
		set doVerbose 1
	    }
	    -lazy {
		set direct 0
		append more " -lazy"
	    }
	    -direct {
		append more " -direct"
	    }
	    -load {
		incr idx
		set loadPat [lindex $args $idx]
		append more " -load $loadPat"
	    }
	    -* {
		return -code error "unknown flag $flag: should be\n$usage"
	    }
	    default {
		# done with the flags
		break
	    }
	}
    }

    set dir [lindex $args $idx]
    set patternList [lrange $args [expr {$idx + 1}] end]
    if {[llength $patternList] == 0} {
	set patternList [list "*.tcl" "*[info sharedlibextension]"]
    }

    set oldDir [pwd]
    cd $dir

    if {[catch {eval glob $patternList} fileList]} {
	global errorCode errorInfo
	cd $oldDir
	return -code error -errorcode $errorCode -errorinfo $errorInfo $fileList
    }
    foreach file $fileList {
	# For each file, figure out what commands and packages it provides.
	# To do this, create a child interpreter, load the file into the
	# interpreter, and get a list of the new commands and packages
	# that are defined.

	if {[string equal $file "pkgIndex.tcl"]} {
	    continue
	}

	# Changed back to the original directory before initializing the
	# slave in case TCL_LIBRARY is a relative path (e.g. in the test
	# suite). 

	cd $oldDir
	set c [interp create]

	# Load into the child any packages currently loaded in the parent
	# interpreter that match the -load pattern.

	if {[string length $loadPat]} {
	    if {$doVerbose} {
		tclLog "currently loaded packages: '[info loaded]'"
		tclLog "trying to load all packages matching $loadPat"
	    }
	    if {![llength [info loaded]]} {
		tclLog "warning: no packages are currently loaded, nothing"
		tclLog "can possibly match '$loadPat'"
	    }
	}
	foreach pkg [info loaded] {
	    if {! [string match $loadPat [lindex $pkg 1]]} {
		continue
	    }
	    if {$doVerbose} {
		tclLog "package [lindex $pkg 1] matches '$loadPat'"
	    }
	    if {[catch {
		load [lindex $pkg 0] [lindex $pkg 1] $c
	    } err]} {
		if {$doVerbose} {
		    tclLog "warning: load [lindex $pkg 0] [lindex $pkg 1]\nfailed with: $err"
		}
	    } elseif {$doVerbose} {
		tclLog "loaded [lindex $pkg 0] [lindex $pkg 1]"
	    }
	    if {[string equal [lindex $pkg 1] "Tk"]} {
		# Withdraw . if Tk was loaded, to avoid showing a window.
		$c eval [list wm withdraw .]
	    }
	}
	cd $dir

	$c eval {
	    # Stub out the package command so packages can
	    # require other packages.

	    rename package __package_orig
	    proc package {what args} {
		switch -- $what {
		    require { return ; # ignore transitive requires }
		    default { eval __package_orig {$what} $args }
		}
	    }
	    proc tclPkgUnknown args {}
	    package unknown tclPkgUnknown

	    # Stub out the unknown command so package can call
	    # into each other during their initialilzation.

	    proc unknown {args} {}

	    # Stub out the auto_import mechanism

	    proc auto_import {args} {}

	    # reserve the ::tcl namespace for support procs
	    # and temporary variables.  This might make it awkward
	    # to generate a pkgIndex.tcl file for the ::tcl namespace.

	    namespace eval ::tcl {
		variable file		;# Current file being processed
		variable direct		;# -direct flag value
		variable x		;# Loop variable
		variable debug		;# For debugging
		variable type		;# "load" or "source", for -direct
		variable namespaces	;# Existing namespaces (e.g., ::tcl)
		variable packages	;# Existing packages (e.g., Tcl)
		variable origCmds	;# Existing commands
		variable newCmds	;# Newly created commands
		variable newPkgs {}	;# Newly created packages
	    }
	}

	$c eval [list set ::tcl::file $file]
	$c eval [list set ::tcl::direct $direct]

	# Download needed procedures into the slave because we've
	# just deleted the unknown procedure.  This doesn't handle
	# procedures with default arguments.

	foreach p {pkg_compareExtension} {
	    $c eval [list proc $p [info args $p] [info body $p]]
	}

	if {[catch {
	    $c eval {
		set ::tcl::debug "loading or sourcing"

		# we need to track command defined by each package even in
		# the -direct case, because they are needed internally by
		# the "partial pkgIndex.tcl" step above.

		proc ::tcl::GetAllNamespaces {{root ::}} {
		    set list $root
		    foreach ns [namespace children $root] {
			eval lappend list [::tcl::GetAllNamespaces $ns]
		    }
		    return $list
		}

		# init the list of existing namespaces, packages, commands

		foreach ::tcl::x [::tcl::GetAllNamespaces] {
		    set ::tcl::namespaces($::tcl::x) 1
		}
		foreach ::tcl::x [package names] {
		    set ::tcl::packages($::tcl::x) 1
		}
		set ::tcl::origCmds [info commands]

		# Try to load the file if it has the shared library
		# extension, otherwise source it.  It's important not to
		# try to load files that aren't shared libraries, because
		# on some systems (like SunOS) the loader will abort the
		# whole application when it gets an error.

		if {[pkg_compareExtension $::tcl::file [info sharedlibextension]]} {
		    # The "file join ." command below is necessary.
		    # Without it, if the file name has no \'s and we're
		    # on UNIX, the load command will invoke the
		    # LD_LIBRARY_PATH search mechanism, which could cause
		    # the wrong file to be used.

		    set ::tcl::debug loading
		    load [file join . $::tcl::file]
		    set ::tcl::type load
		} else {
		    set ::tcl::debug sourcing
		    source $::tcl::file
		    set ::tcl::type source
		}

		# As a performance optimization, if we are creating 
		# direct load packages, don't bother figuring out the 
		# set of commands created by the new packages.  We 
		# only need that list for setting up the autoloading 
		# used in the non-direct case.
		if { !$::tcl::direct } {
		    # See what new namespaces appeared, and import commands
		    # from them.  Only exported commands go into the index.
		    
		    foreach ::tcl::x [::tcl::GetAllNamespaces] {
			if {! [info exists ::tcl::namespaces($::tcl::x)]} {
			    namespace import -force ${::tcl::x}::*
			}

			# Figure out what commands appeared
			
			foreach ::tcl::x [info commands] {
			    set ::tcl::newCmds($::tcl::x) 1
			}
			foreach ::tcl::x $::tcl::origCmds {
			    catch {unset ::tcl::newCmds($::tcl::x)}
			}
			foreach ::tcl::x [array names ::tcl::newCmds] {
			    # determine which namespace a command comes from
			    
			    set ::tcl::abs [namespace origin $::tcl::x]
			    
			    # special case so that global names have no leading
			    # ::, this is required by the unknown command
			    
			    set ::tcl::abs \
				    [lindex [auto_qualify $::tcl::abs ::] 0]
			    
			    if {[string compare $::tcl::x $::tcl::abs]} {
				# Name changed during qualification
				
				set ::tcl::newCmds($::tcl::abs) 1
				unset ::tcl::newCmds($::tcl::x)
			    }
			}
		    }
		}

		# Look through the packages that appeared, and if there is
		# a version provided, then record it

		foreach ::tcl::x [package names] {
		    if {[string compare [package provide $::tcl::x] ""] \
			    && ![info exists ::tcl::packages($::tcl::x)]} {
			lappend ::tcl::newPkgs \
			    [list $::tcl::x [package provide $::tcl::x]]
		    }
		}
	    }
	} msg] == 1} {
	    set what [$c eval set ::tcl::debug]
	    if {$doVerbose} {
		tclLog "warning: error while $what $file: $msg"
	    }
	} else {
	    set what [$c eval set ::tcl::debug]
	    if {$doVerbose} {
		tclLog "successful $what of $file"
	    }
	    set type [$c eval set ::tcl::type]
	    set cmds [lsort [$c eval array names ::tcl::newCmds]]
	    set pkgs [$c eval set ::tcl::newPkgs]
	    if {$doVerbose} {
		tclLog "commands provided were $cmds"
		tclLog "packages provided were $pkgs"
	    }
	    if {[llength $pkgs] > 1} {
		tclLog "warning: \"$file\" provides more than one package ($pkgs)"
	    }
	    foreach pkg $pkgs {
		# cmds is empty/not used in the direct case
		lappend files($pkg) [list $file $type $cmds]
	    }

	    if {$doVerbose} {
		tclLog "processed $file"
	    }
	    interp delete $c
	}
    }

    append index "# Tcl package index file, version 1.1\n"
    append index "# This file is generated by the \"pkg_mkIndex$more\" command\n"
    append index "# and sourced either when an application starts up or\n"
    append index "# by a \"package unknown\" script.  It invokes the\n"
    append index "# \"package ifneeded\" command to set up package-related\n"
    append index "# information so that packages will be loaded automatically\n"
    append index "# in response to \"package require\" commands.  When this\n"
    append index "# script is sourced, the variable \$dir must contain the\n"
    append index "# full path name of this file's directory.\n"

    foreach pkg [lsort [array names files]] {
	set cmd {}
	foreach {name version} $pkg {
	    break
	}
	lappend cmd ::pkg::create -name $name -version $version
	foreach spec $files($pkg) {
	    foreach {file type procs} $spec {
		if { $direct } {
		    set procs {}
		}
		lappend cmd "-$type" [list $file $procs]
	    }
	}
	append index "\n[eval $cmd]"
    }

    set f [open pkgIndex.tcl w]
    puts $f $index
    close $f
    cd $oldDir
}

# tclPkgSetup --
# This is a utility procedure use by pkgIndex.tcl files.  It is invoked
# as part of a "package ifneeded" script.  It calls "package provide"
# to indicate that a package is available, then sets entries in the
# auto_index array so that the package's files will be auto-loaded when
# the commands are used.
#
# Arguments:
# dir -			Directory containing all the files for this package.
# pkg -			Name of the package (no version number).
# version -		Version number for the package, such as 2.1.3.
# files -		List of files that constitute the package.  Each
#			element is a sub-list with three elements.  The first
#			is the name of a file relative to $dir, the second is
#			"load" or "source", indicating whether the file is a
#			loadable binary or a script to source, and the third
#			is a list of commands defined by this file.

proc tclPkgSetup {dir pkg version files} {
    global auto_index

    package provide $pkg $version
    foreach fileInfo $files {
	set f [lindex $fileInfo 0]
	set type [lindex $fileInfo 1]
	foreach cmd [lindex $fileInfo 2] {
	    if {[string equal $type "load"]} {
		set auto_index($cmd) [list load [file join $dir $f] $pkg]
	    } else {
		set auto_index($cmd) [list source [file join $dir $f]]
	    } 
	}
    }
}

# tclMacPkgSearch --
# The procedure is used on the Macintosh to search a given directory for files
# with a TEXT resource named "pkgIndex".  If it exists it is sourced in to the
# interpreter to setup the package database.

proc tclMacPkgSearch {dir} {
    foreach x [glob -directory $dir -nocomplain *.shlb] {
	if {[file isfile $x]} {
	    set res [resource open $x]
	    foreach y [resource list TEXT $res] {
		if {[string equal $y "pkgIndex"]} {source -rsrc pkgIndex}
	    }
	    catch {resource close $res}
	}
    }
}

# tclPkgUnknown --
# This procedure provides the default for the "package unknown" function.
# It is invoked when a package that's needed can't be found.  It scans
# the auto_path directories and their immediate children looking for
# pkgIndex.tcl files and sources any such files that are found to setup
# the package database.  (On the Macintosh we also search for pkgIndex
# TEXT resources in all files.)  As it searches, it will recognize changes
# to the auto_path and scan any new directories.
#
# Arguments:
# name -		Name of desired package.  Not used.
# version -		Version of desired package.  Not used.
# exact -		Either "-exact" or omitted.  Not used.

proc tclPkgUnknown {name version {exact {}}} {
    global auto_path tcl_platform env

    if {![info exists auto_path]} {
	return
    }
    # Cache the auto_path, because it may change while we run through
    # the first set of pkgIndex.tcl files
    set old_path [set use_path $auto_path]
    while {[llength $use_path]} {
	set dir [lindex $use_path end]
	# we can't use glob in safe interps, so enclose the following
	# in a catch statement, where we get the pkgIndex files out
	# of the subdirectories
	catch {
	    foreach file [glob -directory $dir -join -nocomplain \
		    * pkgIndex.tcl] {
		set dir [file dirname $file]
		if {[file readable $file] && ![info exists procdDirs($dir)]} {
		    if {[catch {source $file} msg]} {
			tclLog "error reading package index file $file: $msg"
		    } else {
			set procdDirs($dir) 1
		    }
		}
	    }
	}
	set dir [lindex $use_path end]
	set file [file join $dir pkgIndex.tcl]
	# safe interps usually don't have "file readable", nor stderr channel
	if {([interp issafe] || [file readable $file]) && \
		![info exists procdDirs($dir)]} {
	    if {[catch {source $file} msg] && ![interp issafe]}  {
		tclLog "error reading package index file $file: $msg"
	    } else {
		set procdDirs($dir) 1
	    }
	}
	# On the Macintosh we also look in the resource fork 
	# of shared libraries
	# We can't use tclMacPkgSearch in safe interps because it uses glob
	if {(![interp issafe]) && \
		[string equal $tcl_platform(platform) "macintosh"]} {
	    set dir [lindex $use_path end]
	    if {![info exists procdDirs($dir)]} {
		tclMacPkgSearch $dir
		set procdDirs($dir) 1
	    }
	    foreach x [glob -directory $dir -nocomplain *] {
		if {[file isdirectory $x] && ![info exists procdDirs($x)]} {
		    set dir $x
		    tclMacPkgSearch $dir
		    set procdDirs($dir) 1
		}
	    }
	}
	set use_path [lrange $use_path 0 end-1]
	if {[string compare $old_path $auto_path]} {
	    foreach dir $auto_path {
		lappend use_path $dir
	    }
	    set old_path $auto_path
	}
    }
}

# ::pkg::create --
#
#	Given a package specification generate a "package ifneeded" statement
#	for the package, suitable for inclusion in a pkgIndex.tcl file.
#
# Arguments:
#	args		arguments used by the create function:
#			-name		packageName
#			-version	packageVersion
#			-load		{filename ?{procs}?}
#			...
#			-source		{filename ?{procs}?}
#			...
#
#			Any number of -load and -source parameters may be
#			specified, so long as there is at least one -load or
#			-source parameter.  If the procs component of a 
#			module specifier is left off, that module will be
#			set up for direct loading; otherwise, it will be
#			set up for lazy loading.  If both -source and -load
#			are specified, the -load'ed files will be loaded 
#			first, followed by the -source'd files.
#
# Results:
#	An appropriate "package ifneeded" statement for the package.

proc ::pkg::create {args} {
    append err(usage) "[lindex [info level 0] 0] "
    append err(usage) "-name packageName -version packageVersion"
    append err(usage) "?-load {filename ?{procs}?}? ... "
    append err(usage) "?-source {filename ?{procs}?}? ..."

    set err(wrongNumArgs) "wrong # args: should be \"$err(usage)\""
    set err(valueMissing) "value for \"%s\" missing: should be \"$err(usage)\""
    set err(unknownOpt)   "unknown option \"%s\": should be \"$err(usage)\""
    set err(noLoadOrSource) "at least one of -load and -source must be given"

    # process arguments
    set len [llength $args]
    if { $len < 6 } {
	error $err(wrongNumArgs)
    }
    
    # Initialize parameters
    set opts(-name)		{}
    set opts(-version)		{}
    set opts(-source)		{}
    set opts(-load)		{}

    # process parameters
    for {set i 0} {$i < $len} {incr i} {
	set flag [lindex $args $i]
	incr i
	switch -glob -- $flag {
	    "-name"		-
	    "-version"		{
		if { $i >= $len } {
		    error [format $err(valueMissing) $flag]
		}
		set opts($flag) [lindex $args $i]
	    }
	    "-source"		-
	    "-load"		{
		if { $i >= $len } {
		    error [format $err(valueMissing) $flag]
		}
		lappend opts($flag) [lindex $args $i]
	    }
	    default {
		error [format $err(unknownOpt) [lindex $args $i]]
	    }
	}
    }

    # Validate the parameters
    if { [llength $opts(-name)] == 0 } {
	error [format $err(valueMissing) "-name"]
    }
    if { [llength $opts(-version)] == 0 } {
	error [format $err(valueMissing) "-version"]
    }
    
    if { [llength $opts(-source)] == 0 && [llength $opts(-load)] == 0 } {
	error $err(noLoadOrSource)
    }

    # OK, now everything is good.  Generate the package ifneeded statment.
    set cmdline "package ifneeded $opts(-name) $opts(-version) "
    
    set cmdList {}
    set lazyFileList {}

    # Handle -load and -source specs
    foreach key {load source} {
	foreach filespec $opts(-$key) {
	    foreach {filename proclist} {{} {}} {
		break
	    }
	    foreach {filename proclist} $filespec {
		break
	    }
	    
	    if { [llength $proclist] == 0 } {
		set cmd "\[list $key \[file join \$dir [list $filename]\]\]"
		lappend cmdList $cmd
	    } else {
		lappend lazyFileList [list $filename $key $proclist]
	    }
	}
    }

    if { [llength $lazyFileList] > 0 } {
	lappend cmdList "\[list tclPkgSetup \$dir $opts(-name)\
		$opts(-version) [list $lazyFileList]\]"
    }
    append cmdline [join $cmdList "\\n"]
    return $cmdline
}


