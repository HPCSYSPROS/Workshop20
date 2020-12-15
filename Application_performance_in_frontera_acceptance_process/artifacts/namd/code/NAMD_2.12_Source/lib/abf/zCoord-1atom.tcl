##########################################################################
# ABF procedures for the z coordinate of one atom wrt the COM of a group #
##########################################################################

# Using the actual com of the reference group, and applying force to the group

# This IS rigorous : the force is applied along the gradient, and
# dx/dxi is chosen with zero contributions on the reference atoms

set ABFcoordID "z coordinate of one atom wrt the COM of a group"

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

	addatom $abf2

	set g [addgroup $abf1]
    }
}


################################################################
# ABFcoord : reads coord, returns value of reaction coordinate #
################################################################

proc ABFcoord {} {
    namespace eval ABFcoord {

	loadcoords coords

	return [expr [lindex $coords($abf2) 2] - [lindex $coords($g) 2]]
    }
}

############################################################
# ABForce : returns force along reaction coordinate        #
############################################################

proc ABForce {} {
    namespace eval ABFcoord {

	loadtotalforces forces

	return [lindex $forces($abf2) 2]
    }
}


###############################################################################
# ABFapply : applies the force given as a parameter along reaction coordinate #
###############################################################################

proc ABFapply {force} {
    set ABFcoord::force $force
    
    namespace eval ABFcoord {

	set F2 [vecscale $force "0.0 0.0 1.0"]

	addforce $g [vecinvert $F2]
	addforce $abf2 $F2
    }

    return $force
}
