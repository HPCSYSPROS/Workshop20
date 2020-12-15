###################################################
# ABF procedures for a distance between two atoms #
# Jerome Henin <jerome.henin@uhp-nancy.fr>        #
###################################################

set ABFcoordID "Distance between two atoms (beware of constraints!)"

# Define coordinate-specific optional parameters with default values
array set ABFoptions {
temp		300.0
dxi		0.1
dSmooth		0.3
}

###############################################################
# ABFstartup : declares atoms whose coordinates are requested #
###############################################################

proc ABFstartup {} {
    namespace eval ABFcoord {
	set abf1 $::ABF::abf1
	set abf2 $::ABF::abf2

	addatom $abf1
	addatom $abf2
    }
}

	
################################################################
# ABFcoord : reads coord, returns value of reaction coordinate #
# called first						       #
################################################################

proc ABFcoord {} {

    namespace eval ABFcoord {
	loadcoords coords
	
	set r [veclength [vecsub $coords($abf2) $coords($abf1)]]
	return $r
    }
}

############################################################
# ABForce : returns force along reaction coordinate        #
# called third     					   #
############################################################

proc ABForce {} {

    namespace eval ABFcoord {

	set dr  [vecsub $coords($abf2) $coords($abf1)]
	set nv	[vecnorm $dr] ;# unity vector

	loadtotalforces forces

	set df	[vecsub $forces($abf2) $forces($abf1)]

	# Not including Jacobian term (2kT / r)
	return [expr {[vecdot $df $nv]/ 2.0}]
    }
}


###############################################################################
# ABFapply : applies the force given as a parameter along reaction coordinate #
###############################################################################

proc ABFapply {force} {

    set ABFcoord::force $force

    namespace eval ABFcoord {

	set dr  [vecsub $coords($abf2) $coords($abf1)]
	set r	[veclength $dr]
	set nv	[vecnorm $dr] ;# unity vector abf1 -> abf2

	set force [expr {$force - 2.0 * 0.001986 * $::ABF::temp / $r}]
	# compensate for the Jacobian term

	set F2 [vecscale $force $nv]

	addforce $abf1 [vecinvert $F2]
	addforce $abf2 $F2

	return $force
    }
}
