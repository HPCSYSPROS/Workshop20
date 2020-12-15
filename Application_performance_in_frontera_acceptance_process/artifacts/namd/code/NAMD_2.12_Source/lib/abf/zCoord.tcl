##########################################################################
# ABF procedures for the z coordinate of a group wrt another group       #
##########################################################################

# Using the actual com of the reference group, and applying force to the group

# This IS rigorous : the force is applied along the gradient, and
# dx/dxi is chosen with zero contributions on the reference atoms

set ABFcoordID "z coordinate of a group wrt another group"

# Define coordinate-specific optional parameters with default values
array set ABFoptions {
dxi	    0.1
dSmooth	    0.0
}


###############################################################
# ABFstartup : declares atoms whose coordinates are requested #
###############################################################

proc ABFstartup {} {
    namespace eval ABFcoord {
	set abf1 $::ABF::abf1
	set abf2 $::ABF::abf2

	foreach a $abf2 {addatom $a}
	
	set g1 [addgroup $abf1]
	set g2 [addgroup $abf2]
    }
}


################################################################
# ABFcoord : reads coord, returns value of reaction coordinate #
################################################################

proc ABFcoord {} {
    namespace eval ABFcoord {

	loadcoords coords

	return [expr [lindex $coords($g2) 2] - [lindex $coords($g1) 2]]
    }
}

############################################################
# ABForce : returns force along reaction coordinate        #
############################################################

proc ABForce {} {
    namespace eval ABFcoord {
	loadtotalforces forces

	set f2 0
	foreach a $abf2 {
	    set f2 [expr {$f2 + [lindex $forces($a) 2] } ]
	}
    
	return $f2
    }
}


###############################################################################
# ABFapply : applies the force given as a parameter along reaction coordinate #
###############################################################################

proc ABFapply {force} {
    set ABFcoord::force $force

    namespace eval ABFcoord {
	set F2 [vecscale $force "0.0 0.0 1.0"]

	addforce $g1 [vecinvert $F2]
	addforce $g2 $F2
    }

    return $force
}
