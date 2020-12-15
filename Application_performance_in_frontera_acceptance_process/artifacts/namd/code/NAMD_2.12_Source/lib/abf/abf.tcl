         #############################################
         #           Generic ABF code                #
         # Jerome Henin <jerome.henin@uhp-nancy.fr>  #
         #############################################

############
# Startup  #
############

package provide abf 1.6.0

#########################
# Parameter definitions #
#########################

namespace eval ::ABF {

set version "1.6"

if {! [info exists ABFdir]} { set ABFdir [file dirname [info script]] }
# If it fails, try the local directory
if { $ABFdir == "" } { set ABFdir "." }

TclForces		on
TclForcesScript		$ABFdir/abf_script.tcl

array set defaults {
inFiles         {}
outFile         abf.dat
historyFile     none
distFile        none
forceConst      10.0
fullSamples     1000
outputFreq      5000
df              1.0
fMax            60.0
dSmooth         0.0
writeXiFreq     0
writeFxiFreq    0
usMode		no
applyBias       yes
moveBoundary    0
temp		300.0
abf2		{}
}

set mandatory "coordinate xiMin xiMax dxi abf1"

# these settings are not displayed at startup
set silent "restraintList usMode SFM applyBias direction"

array set capitals {}
foreach param [concat $silent $mandatory [array names defaults]] {
    set capitals([string tolower $param]) $param
    # not set yet
    set alreadySet($param) 0
}

} ;# namespace

proc abf { keyword value } {
    set ::ABF::keyword $keyword
    set ::ABF::value $value

namespace eval ::ABF {

    # Build list of all allowed parameter names
    set list [array names capitals]

    set lowercase [string tolower $keyword]

    # Process parameters
    if {[lsearch $list $lowercase] != -1} {
        set keyword $capitals($lowercase)

	if { $alreadySet($keyword) } {
	    print "ABF> WARNING - multiple definitions of parameter " $keyword	
	}
	
	set $keyword $value
	set alreadySet($keyword) 1
	
	return
    } else {
        error [format "Unknown ABF keyword: %s" $keyword]
    }


} ;# namespace
} ;# proc abf

# define upper-case synonyms to proc abf 
proc ABF { keyword value } {
    abf $keyword $value
}
proc Abf { keyword value } {
    abf $keyword $value
}

