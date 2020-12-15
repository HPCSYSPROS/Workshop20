##########################################################################
# ABF procedures for the distance along a vector between two atom groups #
# Jerome Henin <jerome.henin@uhp-nancy.fr>                               #
##########################################################################

set ABFcoordID "Distance between two atom groups along a vector"

# Define coordinate-specific optional parameters with default values
array set ABFoptions {
dxi		0.2
direction	{0.0 0.0 1.0}
}


###############################################################
# ABFstartup : declares atoms whose coordinates are requested #
###############################################################

proc ABFstartup {} {
    namespace eval ABFcoord {

	print "ABF> Direction for RC: " $::ABF::direction
	set direction [vecnorm $::ABF::direction]

	set abf1 $::ABF::abf1
	set abf2 $::ABF::abf2
	
	foreach a $abf1 {addatom $a}
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

	set x1 [vecdot $coords($g1) $direction]
	set x2 [vecdot $coords($g2) $direction]
	
	return [expr $x2 - $x1]
    }
}

############################################################
# ABForce : returns force along reaction coordinate        #
############################################################

proc ABForce {} {
    namespace eval ABFcoord {

	loadtotalforces forces

	set f1 0.0
	foreach a $abf1 { set f1 [expr $f1 + [vecdot $forces($a) $direction]]}
	set f2 0.0
	foreach a $abf2 { set f2 [expr $f2 + [vecdot $forces($a) $direction]]}

	return [expr ($f2 - $f1) / 2.0]
    }
}


###############################################################################
# ABFapply : applies the force given as a parameter along reaction coordinate #
###############################################################################

proc ABFapply {force} {
    set ABFcoord::force $force

    namespace eval ABFcoord {
	set F2 [vecscale $force $direction]

	addforce $g1 [vecinvert $F2]
	addforce $g2 $F2
    }
    
    return $force
}

