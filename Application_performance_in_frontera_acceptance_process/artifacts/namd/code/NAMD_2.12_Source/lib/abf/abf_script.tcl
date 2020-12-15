         #############################################
         #           Generic ABF code                #
         # Jerome Henin <jerome.henin@uhp-nancy.fr>  #
         #############################################


#######################################################
# Separate source files for coordinate-specific stuff #
#######################################################

enabletotalforces

namespace eval ::ABF {

    print ABF>   ---------------------------------------------  
    print ABF>   Adaptive Biasing Force protocol version $version
    print ABF>   ---------------------------------------------  
    print ABF>
    print "ABF> Using coordinate type : $coordinate"
 
    set sourceFile [format "%s/%s.tcl" $ABFdir $coordinate]

    if ![file exists $sourceFile] {
	error "ABF> TCL file $sourceFile not found"
    } else { source $sourceFile }

    # Print info about the coordinate-specific code
    if [info exists ABFcoordID] { print "ABF> $ABFcoordID" }
 
array set forceDist {}


#####################################################
# Default values for optional parameters            #
#####################################################


# SFM : Scaled-Force Method / Langevin
# This is in development phase, and not official
# For more info, read the code and E. Darve & A. Pohorille, J. Chem. Phys. (2001)

	
# First, use coordinate-specific defaults
if {[array exists ABFoptions]} {foreach option [array names ABFoptions] {
    if {! [info exists $option]} {
	set $option $ABFoptions($option)
	if {[lsearch $silent $option] == -1} {
	    print [format "ABF> %16s : %-14s \[default\]" $option [expr $$option]]
	}
    } else {
	if {[lsearch $silent $option] == -1} {
		print [format "ABF> %16s : %s" $option [expr $$option]]
	}
    }
}}

# Then, the generic ones
foreach option [array names defaults] {
    if {! [info exists $option]} {
	set $option $defaults($option)
	if {[lsearch $silent $option] == -1} {
	    print [format "ABF> %16s : %-14s \[default\]" $option [expr $$option]]
	}
    } elseif { (![array exists ABFoptions] || [lsearch [array names ABFoptions] $option] == -1)} {
	if {[lsearch $silent $option] == -1} {
	    print [format "ABF> %16s : %s" $option [expr $$option]]
	}
    }
}


# Test mandatory parameters

foreach var $mandatory {
    if {! [info exists $var]} {
	error "ABF> Mandatory parameter $var is not set -- cannot start ABF"
    } elseif { (![array exists ABFoptions] || [lsearch [array names ABFoptions] $var] == -1)\
				&& ([lsearch [array names defaults] $var] == -1) \
				&& ($var != "coordinate") } {
	print [format "ABF> %16s : %s" $var [expr $$var]]
    }
}


set timestep 0


# These are useful when using "moving boundaries"
# xiMin0 and xiMin0 + nMax * dxi are the limits of the data arrays
# xiMin and xiMax are the limits between the
# sampling/biasing domain and the confining biases
set xiMin0 $xiMin

set nMax [expr {int (($xiMax - $xiMin)/ $dxi) - 1}]
# this is a bit awkward, but it repairs truncation errors
if { [expr {$xiMin + ($nMax + 1.0) * $dxi}] < $xiMax } {incr nMax}


# These variables will keep track of moving boundaries
set firstBin 0
set lastBin  $nMax

	
for {set i 0} {$i <= $nMax} {incr i} {
	array set samples "$i 0"
	array set accForce "$i 0"
}


if { [info exist SFM ] } {
	print "ABF>"
	print "ABF> *****************************************************"
	print "ABF> Will be using the Scaled Force Method, not ABF" 
	print "ABF> *****************************************************"
	print "ABF>"
}

if { $applyBias != "yes" } {
	print "ABF>"
	print "ABF> *****************************************************"
	print "ABF> *** WARNING --- Biasing force will NOT be applied ***"
	print "ABF> *****************************************************"
	print "ABF>"
	set applyBias no
}

if { $usMode != "no" } {

for {set i 0} {$i <= $nMax} {incr i} {
    array set Usampling "$i 0"
}
	print "ABF>"
	print "ABF> **********************************************"
	print "ABF> *** WARNING --- entering US mode, not ABF  ***"
	print "ABF> **********************************************"
	print "ABF>"
	print "ABF> Accumulating sampling data in [expr {$nMax + 1}] bins"
} else {
	print "ABF> Accumulating force data in [expr {$nMax + 1}] bins"
}

set dSmooth [expr {$dSmooth / $dxi} ]

set minSamples [expr {($fullSamples > 2) ? ($fullSamples / 2) : 1}]


# Remove existing history file

if { $historyFile != "none" } {
    set history [open $historyFile w]
    close $history
}



##################################################
# Read previous results if available             #
##################################################

foreach name $inFiles {
    set in [open $name r]

    print "ABF> Loading previous force data from $name ..."
    
    # Skip initial comment lines
    set ret [gets  $in  line]
    while { ($ret > -1) && ("[string index $line 0]" == "#") } {
	set ret [gets  $in  line]
    }
    
    set count 0

    while { $ret > -1 } {
	set xi [lindex $line 0]
	set f [lindex $line 2]
	set n [lindex $line 3]	

	if { $n && ($xi >= $xiMin) && ($xi < $xiMax) } {
		set i [expr {int(($xi - $xiMin)/$dxi)}] ;# Current bin
		incr samples($i) $n
		set accForce($i) [expr {$accForce($i) + $f * $n}]
		incr count $n 
	}
	set ret [gets  $in  line]
    }
    print "ABF> $count sampling points retrieved"

    close $in
}

}

# The calls for startup procedures are after procedure definitions


# In root namespace (::) for all procedures
# We HAVE to add the following procedure definition
proc veclength {v} {
   return [expr {sqrt([veclength2 $v])}]
}
# Source a file from VMD
source $ABF::ABFdir/vectors.tcl


###################################################
# Procedure to be called by NAMD at each timestep #
###################################################


proc calcforces {} {
namespace eval ::ABF {

# First timestep : we don't have forces
if { $timestep == 0 } {

	# must not be equal to $timestep - 1
	set timeStored -2
    
	set xi	[ABFcoord]	;# for timestep 1
	set cur [expr {int(($xi - $xiMin0)/$dxi)}] ;# Current bin

	writeData   ;# write initial data, if any

	if { $restraintsOn } { restraints }
	
	if { $writeXiFreq } { print ABF> Xi at timestep 0 : $xi}
	
	incr timestep
	return	 	;# nothing more to do in this timestep
}


# ***** Collect the force along the RC *****

if { ($xi >= $xiMin) && ($xi < $xiMax) } {
# To do only if we're in the right range of xi
# Force data is from previous timestep
# Same thing for coordinates
	
    set fxi [expr {($timeStored == $timestep - 1)? \
	    [ABForce] - $storedForce : [ABForce]}]

    if { $usMode != "no" } {
    	incr Usampling($cur)
    	# don't touch forces : just keep track of sampling
    } else {
      	incr samples($cur)
    	set accForce($cur) [expr {$accForce($cur) + $fxi}]
    }
    if { $distFile != "none" } { forceDistrib } ;# Store force distribution in bin n

} elseif { [info exists SFM] } {
    # we need it anyway if we use SFM
    set fxi [expr {($timeStored == $timestep - 1)? \
	    [ABForce] - $storedForce : [ABForce]}]
}

if { $writeFxiFreq && [expr {$timestep % $writeFxiFreq} ] == 0 } {
	print "ABF> Fxi at timestep $timestep : $fxi"
}

# Get coord data for current timestep
set old_xi $xi	;# for SFM
set xi	[ABFcoord]
set cur [expr {int(($xi - $xiMin0)/$dxi)}] ;# Current bin

if { $restraintsOn } { restraints }

if { $writeXiFreq && [expr {$timestep % $writeXiFreq} ] == 0 } {
	print "ABF> Xi at timestep $timestep : $xi"
}

# If requested, use moving boundaries

if { $moveBoundary } {

    if { ($samples($firstBin) > $moveBoundary) && ($cur > $firstBin) && ($firstBin < [expr {$lastBin + 1} ]) } {
	set xiMin [expr {$xiMin + $dxi}]
	incr firstBin
	print "ABF> Increased lower boundary to $xiMin"
    }

    if { ($samples($lastBin) > $moveBoundary) && ($cur < $lastBin) && ($firstBin < [expr {$lastBin + 1} ]) } {
	set xiMax [expr {$xiMax - $dxi}]
	set lastBin [expr {$lastBin - 1}]
	print "ABF> Decreased upper boundary to $xiMax"
    }
}

if { $xi >= $xiMin && $xi < $xiMax && ($applyBias != "no")} {

  if {![info exists SFM]} {
    # We're inside the reaction coordinate range : try to apply ABF
    set F [fSmoothed]

    set storedForce [ABFapply $F]
    set timeStored $timestep

  } else {
    # now we're doing SFM
    set ts 0.001	;# ps
    set friction 50.0	;# amu / ps
    set RT [ expr {8.314 * 300.0} ] ;# J/mol
   
    # velocity along xi 
    set vxi [expr {($xi - $old_xi) / $ts}]
    # stochastic force RMS (kcal/mol/A)
    set sigma [expr {(1.0/4.18e3) * sqrt( 2 * $RT * $friction * 10.0 / $ts)}]

    # normally distributed pseudo-random number
    set sq 2.0
    while { $sq >= 1.0 || $sq == 0.0 } {
	set v1 [expr {2.0 * rand() - 1.0}]
	set v2 [expr {2.0 * rand() - 1.0}]
	set sq [expr {$v1*$v1 + $v2*$v2}]
    }
    set random [expr {$v1 * sqrt(-2.0 * log($sq) / $sq)}]
    
    set FL [ expr {$sigma * $random - $friction * $vxi} ]

    # Langevin force minus (approximate) force along xi
    set storedForce [ABFapply [expr {$FL - $fxi}] ]
    set timeStored $timestep
  }
  
} else {
# We're out of the reaction coordinate range : consider applying bounding restraints
# (if applyBias is "no" and xiMin <= xi <= xiMax, we end up here and do nothing )

    # special treatment for angles
    if {($coordinate == "dihedral") && ($forceConst > 0.0) && (($xi < $xiMin) || ($xi >= $xiMax))} {

        # compute the "real" differences
        set dm [expr abs(fmod(abs($xi-$xiMin)+180, 360) - 180)]
        set dM [expr abs(fmod(abs($xi-$xiMax)+180, 360) - 180)]

        # if we are closer to xiMax, restrain down towards xiMax
        # else restrain up towards xiMin
        set rest  [expr {($dm > $dM) ? -$forceConst * $dM : $forceConst * $dm}]
        ABFapply $rest

    } else {
        if { ($xi < $xiMin) && ($forceConst > 0.0) } {
                set rest  [expr {$forceConst * ($xiMin - $xi)}]
                ABFapply $rest
        }

        if { ($xi > $xiMax ) && ($forceConst > 0.0) } {
                set rest  [expr {$forceConst * ($xiMax - $xi)}]
                ABFapply $rest
        }
    }
}

if { ($outputFreq > 0) && ([expr {$timestep % $outputFreq}]==0) } { writeData }

incr timestep	;# Count timesteps

return
} ;# namespace
}

#########################################
# This computes a smooth weighted	#
# average force				#
#########################################

proc fSmoothed {} {
namespace eval ::ABF {

    if { $dSmooth < 1.0 } {
    # Use the data in one bin
	set i [expr {int( ($xi - $xiMin0) / $dxi)}]
	set n $samples($i)
	if { $n < $minSamples } { return 0 }
	set F [expr {- $accForce($i)/$n}]
	if { $fullSamples < 2} {
	    set factor 1.0
	} else {
	 set factor [expr {($n-$minSamples)/(.0+$fullSamples-$minSamples)}]
	}
	if { $factor >= 1.0 } { return $F }
	return [expr {$factor * $F}]
    }

    #  ($xi - $xiMin0) / $dxi is between 0 and nMax
    set x [expr {($xi - $xiMin0) / $dxi - 0.5}]
    # Now x is between -0.5 and nBins - 0.5
    
    set iMin [expr {($x-$dSmooth < 0.0) ? 0 : int($x-$dSmooth)+1}]
    set iMax [expr {($x+$dSmooth >= $nMax ) ? $nMax : int($x+$dSmooth)}]

    set norm 0.0
    set n 0.0
    set F 0.0

    for { set i $iMin} { $i <= $iMax } { incr i } {
	set c [expr {$dSmooth - abs($i-$x)}]
	# c is the weight of bin i
	set norm [expr {$norm + $c}]
	set n [expr {$n + $c * $samples($i)}]
	set F [expr {$F + $c * $accForce($i)}]
    }
    
    if { $n == 0.0 } { return 0 }

    set F [expr {- $F / $n}]
    # F is the opposite of our average force
    set n [expr {$n / $norm}]
    # n is the (weighted) average amount of sampling

    if { $n < $minSamples } { return 0 }
    if { $fullSamples < 2 } {
	set factor 1.0
    } else {
	set factor [expr {($n-$minSamples)/(.0+$fullSamples-$minSamples)}]
    }
    
    if { $factor >= 1.0 } { return $F }
    return [expr {$F * $factor}] 
} ;# namespace
}


#########################################
# Called by calcforces at last timestep #
# and every outputFreq timesteps	#
#########################################

proc writeData {} {
namespace eval ::ABF {

if { $historyFile != "none" && [expr {($timestep / $outputFreq) % 10 }] == 0 } {
    set history [open $historyFile a]

    puts $history "# Timestep $timestep"
    puts $history "#       xi          A(xi)          av_force      n_samples"
    set histFlag 1
} else {
    set histFlag 0
}

set out [open $outFile w]
puts $out "#       xi          A(xi)          av_force      n_samples"

set pmf 0.0
set force 0

for {set i 0} {$i < [array size samples]} {incr i} {

    set oldForce	$force
    set n $samples($i)

    if {$n > 0} {
	set force [expr {$accForce($i) / $n}]
    } else {
	set force 0.0
    }

    if { $i > 0 } {
	set pmf [expr {$pmf - 0.5 * $dxi * ($force + $oldForce)}]
    }

    if { $usMode != "no" } {
	puts $out [format "%12.3f   %12.4f    %12.4f   %12d"\
	    [expr {$xiMin0 + ($i + 0.5) * $dxi}] $pmf $force $Usampling($i)]
	if { $histFlag } {
	    puts $history [format "%12.3f   %12.4f    %12.4f   %12d"\
	      [expr {$xiMin0 + ($i + 0.5) * $dxi}] $pmf $force $Usampling($i)]
	}
    } else {
	puts $out [format "%12.3f   %12.4f    %12.4f   %12d"\
	    [expr {$xiMin0 + ($i + 0.5) * $dxi}] $pmf $force $n]
	if { $histFlag } {
	    puts $history [format "%12.3f   %12.4f    %12.4f   %12d"\
	      [expr {$xiMin0 + ($i + 0.5) * $dxi}] $pmf $force $n]
	}
    }
}

close $out

if { $histFlag } {
    puts $history "&"
    close $history
}

if {$distFile != "none"} {
    set dist [open $distFile w]

    for {set i 0} {$i < [array size samples]} {incr i} {
	puts $dist "&"
	puts $dist "\# [expr {$xiMin0 + ($i + 0.5) * $dxi}]"
	for {set f [expr {- $fMax}]} {$f < $fMax} {set f [expr {$f + $df}]} {
	    if {[info exists forceDist($i\ $f)]} {
	    puts $dist [format "%12.4f   %d" $f $forceDist($i\ $f)] }
	}
    }
    close $dist
}

print "ABF> Data written to output files at timestep $timestep"

} ;# namespace
}

#########################################
# Collect force distribution            #
#########################################

proc forceDistrib {} {
namespace eval ::ABF {
    set f [expr {floor($fxi/$df + 0.5) * $df}]

    if {[info exists forceDist($cur\ $f)]} {
	incr forceDist($cur\ $f)
    } else {
	set forceDist($cur\ $f) 1
    }
} ;# namespace
}

#########################################
# More or less self explanatory         #
#########################################

proc restraints_init {} {
    # this is called only once

namespace eval ::ABF {
set PI	3.1415926536

foreach r [array names rArray] {
    set type [ lindex $rArray($r) 0 ]

    set ok 0
    switch $type {
	"dist" {
	    set ok 1
	    print "ABF> Restraint $r is a distance"
	    set atoms($r) {}
	    foreach i { 1 2 } {
		set seg [lindex [lindex $rArray($r) $i] 0]
		set res [lindex [lindex $rArray($r) $i] 1]
		set atom [lindex [lindex $rArray($r) $i] 2]
		lappend atoms($r) [atomid $seg $res $atom]
	    }
	    foreach a $atoms($r) {addatom $a}

	    print [format "ABF> Atoms: %-16s k: %6.1f kcal/mol/A�    Ref: %5.1f A"  \
		"($atoms($r))" [lindex $rArray($r) 3] [lindex $rArray($r) 4] ]
	}
       "distLin" {
	    set ok 1
	    print "ABF> Restraint $r is a distance (linear)"
	    set atoms($r) {}
	    foreach i { 1 2 } {
		set seg [lindex [lindex $rArray($r) $i] 0]
		set res [lindex [lindex $rArray($r) $i] 1]
		set atom [lindex [lindex $rArray($r) $i] 2]
		lappend atoms($r) [atomid $seg $res $atom]
	    }
	    foreach a $atoms($r) {addatom $a}

	 print [format "ABF> Atoms: %-16s F: %6.1f kcal/mol/A     Bounds: %5.1f to %5.1f A"  \
	   "($atoms($r))" [lindex $rArray($r) 3] [lindex $rArray($r) 4] [lindex $rArray($r) 5] ]
	}
	"dihe" {
	    set ok 1
	    print ABF> Restraint $r is a dihedral angle
	    set atoms($r) {}
	    foreach i { 1 2 3 4 } {
		set seg [lindex [lindex $rArray($r) $i] 0]
		set res [lindex [lindex $rArray($r) $i] 1]
		set atom [lindex [lindex $rArray($r) $i] 2]
		lappend atoms($r) [atomid $seg $res $atom]
	    }
	    foreach a $atoms($r) {addatom $a}
	 print [format "ABF> Atoms: %-16s k : %6.1f kcal/mol/rad�  Ref: %5.1f deg"  \
	    "($atoms($r))" [lindex $rArray($r) 5] [lindex $rArray($r) 6] ]
	}
	"angle" {
	    set ok 1
	    print "ABF> Restraint $r is a valence angle"
	    set atoms($r) {}
	    foreach i { 1 2 3 } {
		set seg [lindex [lindex $rArray($r) $i] 0]
		set res [lindex [lindex $rArray($r) $i] 1]
		set atom [lindex [lindex $rArray($r) $i] 2]
		lappend atoms($r) [atomid $seg $res $atom]
	    }
	    foreach a $atoms($r) {addatom $a}
	 print [format "ABF> Atoms: %-16s k: %6.1f kcal/mol/rad�  Ref: %5.1f deg"  \
	    "($atoms($r))" [lindex $rArray($r) 4] [lindex $rArray($r) 5] ]
	}
	"cyl" {
	    set ok 1
	    print "ABF> Restraint $r is a cylindrical restraint"
	    set atoms($r) {}
	    foreach i [lindex $rArray($r) 1] {
		set seg [lindex $i 0]
		set res [lindex $i 1]
		set atom [lindex $i 2]
		lappend atoms($r) [atomid $seg $res $atom]
	    }
	    set dir [lindex $rArray($r) 4] 

	    # normalize direction and add the group ID at end of list
	    set rArray($r) [lreplace $rArray($r) 4 4 [vecnorm $dir] [addgroup $atoms($r)]]

	    print [format "ABF> Atoms: %-16s k: %6.1f kcal/mol/A�  Ref: %s   Dir: %s"  \
		"($atoms($r))" [lindex $rArray($r) 2] [lindex $rArray($r) 3] $dir ]
	}
	"harm" {
	    set ok 1
	    print "ABF> Restraint $r is a generic harmonic restraint"
	    set seg [lindex [lindex $rArray($r) 1] 0]
	    set res [lindex [lindex $rArray($r) 1] 1]
	    set name [lindex [lindex $rArray($r) 1] 2]

	    set atoms($r) [atomid $seg $res $name]
	    addatom $atoms($r)

	    print [format "ABF> Atom: %-17s k: %s kcal/mol/A�  Ref: %s"  \
		$atoms($r) [lindex $rArray($r) 2] [lindex $rArray($r) 3] ]
	}
    }
    if { $ok == 0 } {
	print "ABF> Restraint $r is of unknown type $type !"
    }
} ;# foreach

} ;# namespace
} ;# proc restraints_init


#########################################
# More or less self explanatory         #
#########################################

proc restraints {} {
namespace eval ABF {

    loadcoords coords

    # Loop on requested restraints
    foreach r [ array names atoms ] {

    set type [ lindex $rArray($r) 0 ]
    switch $type {
	"dist" {
# Apply a harmonic restraint to a distance

	    foreach { a1 a2 } $atoms($r) {}
	    set k [lindex $rArray($r) 3]
	    set ref [lindex $rArray($r) 4]

	    set d [getbond $coords($a1) $coords($a2)]
	    set vect [vecscale [expr {1.0 / $d}] [vecsub $coords($a2) $coords($a1)]]
	    set force [expr {-$k * ($d - $ref)}]

	    addforce $a1 [vecscale [expr {- $force}] $vect]
	    addforce $a2 [vecscale  $force $vect]
	}
	"distLin" {
# Apply a linear bias to a distance

	    foreach { a1 a2 } $atoms($r) {}
	    set d1 [lindex $rArray($r) 4]
	    set d2 [lindex $rArray($r) 5]

	    set d [getbond $coords($a1) $coords($a2)]

	    if { $d >= $d1 && $d < $d2 } {
		set vect [vecscale [expr {1.0 / $d}] [vecsub $coords($a2) $coords($a1)]]
		set force [expr {- [lindex $rArray($r) 3] }]

		addforce $a1 [vecscale [expr {- $force}] $vect]
		addforce $a2 [vecscale $force $vect]
	    }
	}
	"dihe" {
# Apply a restraint to a dihedral

	    foreach { a1 a2 a3 a4 } $atoms($r) {}
	    set k [lindex $rArray($r) 5]
	    set ref [lindex $rArray($r) 6]

	    set phi [getdihedral $coords($a1) $coords($a2) $coords($a3) $coords($a4)]

	    # we need phi between $ref - 180.0 and ref + 180.0
	    set phi [ expr {$phi + 360.0 * floor( ($ref + 180.0 - $phi) / 360.0)}] 

	    # (in radians)
	    set diff [expr {($phi - $ref) * $PI/180.0}]
		    
	    set force [expr {-$k * $diff}]

	    foreach {g1 g2 g3 g4} [dihedralgrad $coords($a1) $coords($a2) $coords($a3) $coords($a4)] {}

	    addforce $a1 [vecscale $g1 $force]
	    addforce $a2 [vecscale $g2 $force]
	    addforce $a3 [vecscale $g3 $force]
	    addforce $a4 [vecscale $g4 $force]
	}
	"angle" {
# Apply a restraint to a valence angle

	    foreach { a1 a2 a3 } $atoms($r) {}
	    set k [lindex $rArray($r) 4]
	    set ref [lindex $rArray($r) 5]

	    set phi [getangle $coords($a1) $coords($a2) $coords($a3)]

	    # (in radians)
	    set diff [expr {($phi - $ref) * $PI/180.0}]
		    
	    set force [expr {-$k * $diff}]
	
	    foreach {g1 g2 g3} [anglegrad $coords($a1) $coords($a2) $coords($a3)] {}
	   
	    addforce $a1 [vecscale $g1 $force]
	    addforce $a2 [vecscale $g2 $force]
	    addforce $a3 [vecscale $g3 $force]
	}
	"cyl" {
	    set k     [lindex $rArray($r) 2]
	    set ref   [lindex $rArray($r) 3]
	    set u     [lindex $rArray($r) 4]
	    set group [lindex $rArray($r) 5]
	    set pos   $coords($group)

	    set vecd [vecsub $pos $ref]

	    set v [vecnorm [veccross $vecd $u]]
	    set vech [veccross $v $u]
	    set dot [vecdot $u $vecd]
	    
	    set d [veclength $vecd]
	    set h [expr {sqrt($d*$d - $dot*$dot)} ]

	    addforce $group [vecscale [expr {$k * $h} ] $vech ] 
	}
	"harm" {
	    set a $atoms($r)
	    foreach { x y z } $coords($a) {}
	    foreach { kx ky kz } [lindex $rArray($r) 2] {}
	    foreach { x0 y0 z0 } [lindex $rArray($r) 3] {}

	    foreach { Fx Fy Fz } { 0. 0. 0.} {}
	    if { $kx } { set Fx [expr {$kx * ($x0 - $x)} ] }
	    if { $ky } { set Fy [expr {$ky * ($y0 - $y)} ] }
	    if { $kz } { set Fz [expr {$kz * ($z0 - $z)} ] }

	    addforce $a [list $Fx $Fy $Fz]
	}
  } ;# switch $type

} ;# foreach

} ;# namespace
} ;# proc restraints

#################################################
# End of procedure definitions			#
#################################################

namespace eval ::ABF {

# coordinate-specific startup procedure
ABFstartup

if { [info exists restraintList] } {

    array set rArray $restraintList

	set restraintsOn 1
	restraints_init
} else {set restraintsOn 0}

} ;# namespace

return	;# otherwise return value is 0
