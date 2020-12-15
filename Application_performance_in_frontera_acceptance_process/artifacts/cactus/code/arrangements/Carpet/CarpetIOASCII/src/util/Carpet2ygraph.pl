#! /usr/bin/perl -sw
#
# Blame Ian Hawke for not checking that Scott had already written a C
# program to do this add so writing a perl version instead
#
# Given an output file in CarpetIOASCII 1d format, strip to ygraph format.
# The arguments should be direction (x=0,y=1,z=2) and filename.
# Output is to a file. Only the base name should be given.
# The base will be appended with _<level>.xg
#
# Example:
#
# Carpet2ygraph.pl 0 alp.xl alp_x
#
# produces alp_x_<d>.xg where <d> is the refinement level.
#
# The headers of the Carpet files are just copied with the addition of
# the correct time for ygraph, and the different comment marker.
#
# This script is a much altered version of Tom Goodale's convergence
# testing scripts.
#

use FileHandle;

if(@ARGV != 3 or $ARGV[0] !~ /^0|1|2$/)
{
  print "Usage: $0 direction <Inputfile> <Outputfile>\n";
  print "       where direction is either 0, 1, 2\n";
  exit;
}

open(CARPETFILE,   "<$ARGV[1]") || die "Unable to find file \"$ARGV[1]\"."; 

my $direction = $ARGV[0]+9; 
my $flag = 0;
my $newflag = 0;
my $refinementlevel = 0;
my %componentflag = ();
$componentflag{$refinementlevel} = 0;
my @outputdata = ("\n");

#
# Open the output file for the base grid. 
#
my $file = $ARGV[2]."_$refinementlevel.xg";
my %outputfilelist = ();
$outputfilelist{$refinementlevel} =
  new FileHandle(">$file") || die "Unable to open file \"$file\"."; 

#
# Find the correct column for the spatial coordinate; requires a magic number
#
while (<CARPETFILE>)
{
    # skip empty lines
    next if (/^(\s)*$/);

    $line = $_;

    if(/^\#/) # The line is a header comment
    {
	if ($flag==1) # It's a new level and there is data to output
        {
	    my $fh = $outputfilelist{$refinementlevel};
	    print $fh @outputdata;
	    @outputdata=("\n");
	    $flag = 0;
	}
	if ($line =~ /refinement level ([0-9]{1,2})/) # Line gives ref. level
	{
	    $refinementlevel = $1;
	    $line =~ /component ([0-9+])/;
	    $componentflag{$refinementlevel} = $1; # Which component?

	    # If no file exists for this refinement level, 
	    # open and add filehandle to list
	    if (not defined $outputfilelist{$refinementlevel})
	    { 
		$file = $ARGV[2] . "_$refinementlevel.xg";
		$outputfilelist{$refinementlevel} = 
		  new FileHandle(">$file") || die "Unable to open file \"$file\"."; 
	    }
	}
	# Only output the headers if this is the zero component
	# FIXME: what happens if component 0 isn't output first?
	if (0 == $componentflag{$refinementlevel})
	{
	    push(@outputdata, ("\"",$line)); # Add ygraph comment marker
	    $newflag = 0;
	}
	else
	{
	    $flag = 1;
	    @outputdata=("");
	}
    }
    else # The line contains real data
    {
	@data = split(/[ \t]+/,$line);
	if (($newflag==0)) # This is the first line of data
	{
	    $newflag = 1;
	    my $timeset = $data[8]; # Magic number gives the Cactus time
	    @outputdata = ("\n\n\#Time = $timeset \n",@outputdata);
	}
        chomp ($data[12]);
        push(@outputdata, "$data[$direction] $data[12]\n");
    }
}

#
# At end of file, so output final data set.
#
my $fh = $outputfilelist{$refinementlevel};
print $fh @outputdata;

foreach $refinementlevel (keys %outputfilelist)
{
    close($outputfilelist{$refinementlevel});
}
