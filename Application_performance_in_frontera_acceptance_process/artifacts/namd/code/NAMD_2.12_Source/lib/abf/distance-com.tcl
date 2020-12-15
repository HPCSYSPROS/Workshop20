#########################################################
# ABF procedures for a distance between two atom groups #
#########################################################


set ABFcoordID "Distance between COM of two atom groups"

# Define coordinate-specific optional parameters with default values
array set ABFoptions {
temp		300.0
dxi		0.1
dSmooth		0.2
}


###############################################################
# ABFstartup : declares atoms whose coordinates are requested #
###############################################################

proc ABFstartup {} {

    namespace eval ABFcoord {
	set abf1 $::ABF::abf1
	set abf2 $::ABF::abf2

	# we need this for 'loadtotalforces'
	foreach a $abf1 {addatom $a}
	foreach a $abf2 {addatom $a}

	# this one is convenient for 'loadcoords' and 'addforce'
	set group1 [ addgroup $abf1 ]
	set group2 [ addgroup $abf2 ]
    }
}

	
################################################################
# ABFcoord : reads coord, returns value of reaction coordinate #
################################################################

proc ABFcoord {} {
    namespace eval ABFcoord {
	loadcoords coords

	set r [veclength [vecsub $coords($group2) $coords($group1)]]
	return $r
    }
}

############################################################
# ABForce : returns force along reaction coordinate        #
############################################################

proc ABForce {} {

    namespace eval ABFcoord {
	set dr  [vecsub $coords($group2) $coords($group1)]
	set nv  [vecnorm $dr] ;# unity vector

	loadtotalforces forces

	set f1 0.0
	foreach a $abf1 { set f1 [expr {$f1 + [vecdot $forces($a) $nv]}]}
	set f2 0.0
	foreach a $abf2 { set f2 [expr {$f2 + [vecdot $forces($a) $nv]}]}

	return [expr {($f2 - $f1) / 2.0}]
    }
}


###############################################################################
# ABFapply : applies the force given as a parameter along reaction coordinate #
###############################################################################

proc ABFapply {force} {

    set ABFcoord::force $force

    namespace eval ABFcoord {

	set dr  [vecsub $coords($group2) $coords($group1)]
	set r	[veclength $dr]
	set nv	[vecnorm $dr] ;# unity vector group1 -> group2

	# compensate for the Jacobian term 2kT/r
	set force [expr {$force - 2.0 * 0.001986 * $::ABF::temp / $r}]

	set F2 [vecscale  $force $nv]

	addforce $group1 [vecinvert $F2]
	addforce $group2 $F2

	return $force
    }
}
