#!/bin/env namd2 +tclsh
# eabf script by Haohao Fu (fhh2626_at_mail.nankai.edu.cn)
# nightly version 2016.12.5

package require Eabf

set eabf_dimension 0
set eabf_biases []
set eabf_params []
set eabf_outputgrad 1
set eabf_outputfreq 5000

set eabf_readconfig 0
set eabf_column []

set eabf_BOLTZMANN 0.0019872041

# get the name of biases, also get the dimensions
proc eabf_get_abf_biases {eabf_dimension eabf_biases} {

	upvar $eabf_dimension edimension
	upvar $eabf_biases ebiases

	set abfconfig [cv bias abf1 getconfig]
	set lines [split $abfconfig "\n"]

	foreach line $lines {
		set _biasline [string match -nocase *colvars* $line]
		if {$_biasline} {
			set biasline [split [lrange [eabf_splitline $line] 1 end]]
			break
		}
	}

	set edimension [llength $biasline]
	lappend ebiases {*}$biasline
}

# get the parameters of a bias
proc eabf_get_bias_parameters {bias eabf_params} {

	upvar $eabf_params eparams
	global eabf_temperature eabf_BOLTZMANN
	
	set biasconfig [cv colvar $bias getconfig]
	set lines [split $biasconfig "\n"]

	foreach line $lines {
		set _lowerboundary [string match -nocase *lowerboundary* $line]
		set _upperboundary [string match -nocase *upperboundary* $line]
		set _width [string match -nocase *width* $line]
		set _extendedFluctuation [string match -nocase *extendedFluctuation* $line]

		if {$_lowerboundary} {
			set lowerboundary [lindex [eabf_splitline $line] 1]
		}
		if {$_upperboundary} {
			set upperboundary [lindex [eabf_splitline $line] 1]
		}
		if {$_width} {
			set width [lindex [eabf_splitline $line] 1]
		}
		if {$_extendedFluctuation} {
			set extendedFluctuation [lindex [eabf_splitline $line] 1]
		}
	}

	set Krestr [expr $eabf_temperature * $eabf_BOLTZMANN / ($extendedFluctuation * $extendedFluctuation)]
	
	lappend eparams $lowerboundary $upperboundary $width $Krestr
}

# split a line when more tha one space between two string
# stupid tcl!
proc eabf_splitline {line} {
	set line2 [string trim $line]
	regsub -all {[[:blank:]]+} $line2 " " line3
	return [split $line3]
}


# read the label in .colvar.traj file
# judge which column is necessary
proc read_label {label column} {
	upvar $column col
	global eabf_biases

	set splitedlabel [eabf_splitline $label]

	foreach i $eabf_biases {
		set flag 0
		foreach j $splitedlabel {
			if {$j == $i} {
				lappend col [expr $flag - 1]
			}
			set flag [expr $flag + 1]
		}
		set flag 0
		foreach j $splitedlabel {
			set str "r_"
			append str $i
			if {$j == $str} {
				lappend col [expr $flag - 1]
			}
			set flag [expr $flag + 1]
		}
	}
}

		
# called by colvar module every MD step
proc calc_colvar_forces {step} {

	global eabf_readconfig

	# read the information of collective variables and send them to cpp library
	if {$eabf_readconfig == 0} {

		global eabf_biases eabf_params eabf_dimension eabf_column

		eabf_get_abf_biases eabf_dimension eabf_biases
		foreach bias $eabf_biases {
			eabf_get_bias_parameters $bias eabf_params
		}

		read_label [cv printframelabels] eabf_column

		# I hate the interface, but too lazy to rewrite it
		global eabf_outputgrad eabf_outputfreq eabf_temperature eabf_outputname eabf_inputname 
		set lowerboundary [lindex $eabf_params 0]
		set upperboundary [lindex $eabf_params 1]
		set width [lindex $eabf_params 2]
		set krestr [lindex $eabf_params 3]
		set outputfile $eabf_outputname
		set outputfreq $eabf_outputfreq
		set outputgrad 1
		set gradfreq $eabf_outputfreq
		set dimension $eabf_dimension
		set inputfile $eabf_inputname
		set temperature $eabf_temperature
		if {$inputfile == 0} {
			set restart 0
		} else {
			set restart 1
		}

		if {$dimension == 1} {
			startrun $lowerboundary $upperboundary $width $krestr $outputfile $outputfreq $restart $inputfile $outputgrad $gradfreq $temperature $dimension
			setcol [lindex $eabf_column 0] [lindex $eabf_column 1]
		} elseif {$dimension ==2} {
			set lowerboundary2 [lindex $eabf_params 4]
			set upperboundary2 [lindex $eabf_params 5]
			set width2 [lindex $eabf_params 6]
			set krestr2 [lindex $eabf_params 7]
			startrun $lowerboundary $upperboundary $width $krestr $lowerboundary2 $upperboundary2 $width2 $krestr2 $outputfile $outputfreq $restart $inputfile $outputgrad $gradfreq $temperature $dimension
			setcol [lindex $eabf_column 0] [lindex $eabf_column 1] [lindex $eabf_column 2] [lindex $eabf_column 3]
		}
		set eabf_readconfig 1
	}

	updaterun [cv printframe]

}


# calculate the RMSD of two list
# usually, a list stores the grad of a PMF
proc eabf_calc_RMSD {list1 list2} {

	if {[llength $list1] != [llength $list2]} {
		return 9999999
	}

	set num [llength $list1]
	set sum 0

	foreach i $list1 j $list2 {
		set sum [expr ($i-$j)*($i-$j)]
	}

	set RMSD [expr sqrt($sum/$num)]

	return $sum
}

# merge file produced by mwabf
proc eabf_merge_mwabf {outputfile args} {

	# read the dimension, lowerboundary, etc
	set flag 0
	foreach i $args {
		set parfile [open $i r]
		gets $parfile line
		set splitedline [eabf_splitline $line]
		close $parfile

		set dimension [expr ([llength $splitedline] - 1) / 4]
		
		set temp_lowerboundary [lindex $splitedline 1]
		set temp_upperboundary [expr [lindex $splitedline 1] + [lindex $splitedline 2] * [lindex $splitedline 3]]

		if {$flag == 0} {
			set lowerboundary $temp_lowerboundary
			set upperboundary $temp_upperboundary
			if {$dimension == 1} {
				set flag 1
			}

		} else {
			if {$temp_lowerboundary < $lowerboundary} {
				set lowerboundary $temp_lowerboundary
			}
			if {$temp_upperboundary > $upperboundary} {
				set upperboundary $temp_upperboundary
			}
		}

		if {$dimension == 2} {

			set temp_lowerboundary2 [lindex $splitedline 5]
			set temp_upperboundary2 [expr [lindex $splitedline 5] + [lindex $splitedline 6] * [lindex $splitedline 7]]

			if {$flag == 0} {
				set lowerboundary2 $temp_lowerboundary2
				set upperboundary2 $temp_upperboundary2
				set flag 1
			} else {
				if {$temp_lowerboundary2 < $lowerboundary2} {
					set lowerboundary2 $temp_lowerboundary2
				}
				if {$temp_upperboundary2 > $upperboundary2} {
					set upperboundary2 $temp_upperboundary2
				}
			}
		}
	}
	
	set width [lindex $splitedline 2]

	if {$dimension ==1} {
		mergefile $dimension $lowerboundary $upperboundary $width $outputfile {*}$args
	}

	if {$dimension == 2} {
		set width2 [lindex $splitedline 6]
		mergefile $dimension $lowerboundary $upperboundary $width $lowerboundary2 $upperboundary2 $width2 $outputfile {*}$args
	}
	
}

# read the dimension, lowerboundary, width and upperboundary of a grad file
proc eabf_read_grad_head {filename dimension params} {

	upvar $params pr
	upvar $dimension di

	set pr []

	# read the dimension, lowerboundary, etc
	set parfile [open $filename r]
	gets $parfile line
	set splitedline [eabf_splitline $line]
	# dimension
	set di [lindex $splitedline 1]

	# lowerboundary width and upperboundary
	for {set i 0} {$i < $di} {incr i} {
		gets $parfile line
		set splitedline [eabf_splitline $line]
		lappend pr [lindex $splitedline 1]
		lappend pr [lindex $splitedline 2]
		set upperboundary [expr [lindex $splitedline 3] * [lindex $splitedline 2] + [lindex $splitedline 1]]
		lappend pr $upperboundary
	}

	close $parfile
}

# get range of grad
proc eabf_grad_range {dimension params args} {
	
	upvar $dimension di
	upvar $params pr

	set temppr []
	# flag, read first file
	set first_time 1

	foreach i $args {
		eabf_read_grad_head $i di temppr
		if {$first_time == 1} {
			lappend pr {*}$temppr
			set first_time 0
		} else {
			if {[lindex $temppr 0] < [lindex $pr 0]} {
				set pr [lreplace $pr 0 0 [lindex $temppr 0]]
			}
			if {[lindex $temppr 2] > [lindex $pr 2]} {
				set pr [lreplace $pr 2 2 [lindex $temppr 2]]
			}

			if {$di == 2} {
				if {[lindex $temppr 3] < [lindex $pr 3]} {
					set pr [lreplace $pr 3 3 [lindex $temppr 3]]
				}
				if {[lindex $temppr 5] > [lindex $pr 5]} {
					set pr [lreplace $pr 5 5 [lindex $temppr 5]]
				}
			}
		}
	}
}

# set grid of a given range
proc eabf_set_grid {dimension params grid} {

	upvar $grid gr
	upvar $params pr
	
	# set grid
	if {$dimension == 1} {
		set gridnum [expr ( [lindex $pr 2] - [lindex $pr 0] ) / [lindex $pr 1]]
		for {set i 0} {$i < $gridnum} {incr i} {
			lappend gr 0
		}
	}

	if {$dimension == 2} {
		set gridnum [expr ( [lindex $pr 2] - [lindex $pr 0] ) / [lindex $pr 1]]
		set gridnum2 [expr ( [lindex $pr 5] - [lindex $pr 3] ) / [lindex $pr 4]]
		for {set i 0} {$i < [expr $gridnum * $gridnum2]} {incr i} {
			# one dimensional grid to reflect multi-dimensional grads
			lappend gr 0 0
		}
	}
}

# read the grad from file and set it into grid
proc eabf_read_grad_splitwindow {filename count_filename dimension params grid count_grid} {
	
	upvar $grid gr
	upvar $count_grid cgr
	upvar $params pr

	set file_param []
	eabf_read_grad_head $filename dimension file_param
	set grad_num [expr ([lindex $file_param 2] - [lindex $file_param 0]) / [lindex $file_param 1]]
	if {$dimension == 2} {
		set grad_num [expr $grad_num * [expr ([lindex $file_param 5] - [lindex $file_param 3]) / [lindex $file_param 4]]]
	}

	set parfile [open $filename r]
	set countfile [open $count_filename r]

	# these lines are useless
	for {set i 0} {$i < 3} {incr i} {
		gets $parfile temp_line
		gets $countfile temp_line
	}
	if {$dimension == 2} {
		gets $parfile temp_line
		gets $countfile temp_line
	}
	
	for {set i 0} {$i < $grad_num} {incr i} {
		set empty [gets $parfile temp_line]
		gets $countfile count_temp_line
		if {$empty == 0} {
			set i [expr $i - 1]
			continue
		}
		set line [eabf_splitline $temp_line]
		set count_line [eabf_splitline $count_temp_line]
		if {$dimension == 1} {
			if {[lindex $count_line 1] != 0} {
				set pos [expr round(([lindex $line 0] - [lindex $pr 1] / 2.0 - [lindex $pr 0]) / [lindex $pr 1])]
				set gr [lreplace $gr $pos $pos [expr ([lindex $line 1] * [lindex $count_line 1] + [lindex $gr $pos] * [lindex $cgr $pos]) / ([lindex $count_line 1] + [lindex $cgr $pos])]]
				set cgr [lreplace $cgr $pos $pos [expr [lindex $count_line 1] + [lindex $cgr $pos]]]
			}
		} elseif {$dimension == 2} {
			if {[lindex $count_line 2] != 0} {
				set pos [expr round(2.0 * ((([lindex $line 0] - [lindex $pr 1] / 2.0 - [lindex $pr 0]) / [lindex $pr 1] * ([lindex $pr 5] - [lindex $pr 3]) / [lindex $pr 4]) + ([lindex $line 1] - [lindex $pr 4] / 2.0 - [lindex $pr 3]) / [lindex $pr 4]))]
				set gr [lreplace $gr $pos $pos [expr ([lindex $line 2] * [lindex $count_line 2] + [lindex $gr $pos] * [lindex $cgr $pos]) / ([lindex $count_line 2] + [lindex $cgr $pos])]]
				set gr [lreplace $gr [expr $pos + 1] [expr $pos + 1] [expr ([lindex $line 3] * [lindex $count_line 2] + [lindex $gr [expr $pos + 1]] * [lindex $cgr $pos]) / ([lindex $count_line 2] + [lindex $cgr $pos])]]
				set cgr [lreplace $cgr $pos $pos [expr [lindex $count_line 2] + [lindex $cgr $pos]]]
			}
		}
	}
	close $parfile
}
		
# merge grad produced by split windows
proc eabf_merge_split_window {outputname args} {

	set dimension 0
	set params []
	set grid []
	set counts_grid []

	set grad_args []
	set counts_args []

	foreach i $args {
		lappend grad_args $i.grad
		lappend counts_args $i.count
	}

	eabf_grad_range dimension params {*}$grad_args
	eabf_set_grid $dimension params grid
	eabf_set_grid $dimension params counts_grid

	foreach i $grad_args j $counts_args {
		eabf_read_grad_splitwindow $i $j $dimension params grid counts_grid
	}

	# write output
	set mergefile [open $outputname.grad w]
	set countfile [open $outputname.count w]
	if {$dimension == 1} {
		puts $mergefile "# $dimension"
		puts $mergefile "# [lindex $params 0] [lindex $params 1] [expr int(0.0000001 + ([lindex $params 2] - [lindex $params 0]) / [lindex $params 1])] 0 \n"
		puts $countfile "# $dimension"
		puts $countfile "# [lindex $params 0] [lindex $params 1] [expr int(0.0000001 + ([lindex $params 2] - [lindex $params 0]) / [lindex $params 1])] 0 \n"

		for {set i 0} {$i < [llength $grid]} {incr i} {
			puts $mergefile "[expr [lindex $params 0] + $i * [lindex $params 1] + [lindex $params 1] / 2.0] [lindex $grid $i]"
			puts $countfile "[expr [lindex $params 0] + $i * [lindex $params 1] + [lindex $params 1] / 2.0] [lindex $counts_grid $i]"
		}
	}

	if {$dimension == 2} {
		puts $mergefile "# $dimension"
		puts $mergefile "# [lindex $params 0] [lindex $params 1] [expr int(0.0000001 + ([lindex $params 2] - [lindex $params 0]) / [lindex $params 1])] 0"
		puts $mergefile "# [lindex $params 3] [lindex $params 4] [expr int(0.0000001 + ([lindex $params 5] - [lindex $params 3]) / [lindex $params 4])] 0 \n"

		puts $countfile "# $dimension"
		puts $countfile "# [lindex $params 0] [lindex $params 1] [expr int(0.0000001 + ([lindex $params 2] - [lindex $params 0]) / [lindex $params 1])] 0"
		puts $countfile "# [lindex $params 3] [lindex $params 4] [expr int(0.0000001 + ([lindex $params 5] - [lindex $params 3]) / [lindex $params 4])] 0 \n"

		set k 0
		for {set i [lindex $params 0]} {[expr $i - [lindex $params 2]] < -0.00001} {set i [expr $i + [lindex $params 1]]} {
			for {set j [lindex $params 3]} {[expr $j - [lindex $params 5]] < -0.00001} {set j [expr $j + [lindex $params 4]]} {
				puts $mergefile "[expr $i + [lindex $params 1] / 2.0] [expr $j + [lindex $params 4] / 2.0] [lindex $grid $k] [lindex $grid [expr $k + 1]]"
				puts $countfile "[expr $i + [lindex $params 1] / 2.0] [expr $j + [lindex $params 4] / 2.0] [lindex $counts_grid $k]"
				incr k 2
			}
			puts $mergefile ""
			puts $countfile ""
		}
	}
	close $mergefile
	close $countfile
}


# this is the terminal interface

if {[info exist argc]} {
	if {$argc > 0} {
		# ./eabf.tcl -mergemwabf merge.eabf eabf.0 eabf.1 eabf.2
		if {[lindex $argv 0] == "-mergemwabf"} {
			eabf_merge_mwabf [lindex $argv 1] {*}[lrange $argv 2 end]
		}
		# ./eabf.tcl -mergesplitwindow merge.grad xxx1.grad xxx2.grad xxx3.grad
		if {[lindex $argv 0] == "-mergesplitwindow"} {
			eabf_merge_split_window [lindex $argv 1] {*}[lrange $argv 2 end]
		}

	}
}

