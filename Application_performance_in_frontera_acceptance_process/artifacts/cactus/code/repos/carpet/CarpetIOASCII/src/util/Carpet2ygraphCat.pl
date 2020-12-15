#! /usr/bin/perl -sw
#
# Given an output file in CarpetIOScalar format or
# CarpetIOASCII scalar or 1d format, strip to ygraph format.

# Output is to a file named after the input file, with its extension ".asc"
# substituted by ".xg".
#
# Example:
#
#   Carpet2ygraphCat.pl alp.x.asc
#
# produces alp.x.xg
#
#   Carpet2ygraph.pl alp.minimum.asc
#
# produces alp.minumum.xg
#
# Some headers of the Carpet files are just copied with the addition of
# the correct time for ygraph, and the different comment marker.
#
# This script is an altered version of Ian Hawke's script.
#


use FileHandle;

if(@ARGV != 1 and @ARGV != 2)
{
  print "Usage: $0 <Inputfile> [<Data Column (0..n)>]\n";
  exit;
}

my $direction;
my $infile = shift;
my $varindex = 0;
$varindex = shift if @ARGV == 1;

open(CARPETFILE, "< $infile") || die "Unable to find input file '$infile'\n.";

if ($infile =~ /\.(x|y|z|d)\.asc$/) {
  $direction = ord($1) - ord('x') + 9;
  print "extracting CarpetIOASCII 1D data in direction $1...\n";
} elsif ($infile =~ /\.\.asc$/) {
  $direction = 8;
  print "extracting CarpetIOASCII 0D data...\n";
} elsif ($infile =~ /\.asc$/) {
  $direction = 0;
  print "extracting CarpetIOScalar data...\n";
} else {
  die "Do not recognize input file format for '$infile'.\n";
}


# deal with scalar output
if ($direction == 0) {
  local @dataColumns = ();
  while (<CARPETFILE>) {
    if (/^# data columns: (.+)/) {
      @dataColumns = split(/[ :]/, $1);
      last;
    }
  }
  # seek back to the beginning of the file
  seek(CARPETFILE, 0, 0);

  my $outfile = $infile;
  $outfile =~ s/asc$/xg/;
  if (exists $dataColumns[$varindex*2 + 1]) {
    local @tokens = split(/::/, $outfile);
    if (@tokens) {
      $outfile = $tokens[0] . '::';
      @tokens = split(/\./, $tokens[1]);
    } else {
      @tokens = split('.', $outfile);
      $outfile = '';
    }
    $tokens[0] = $dataColumns[$varindex*2 + 1];
    $outfile .= join('.', @tokens);
  }
  print "  output to file '$outfile'\n\n";

  local $fh = new FileHandle("> $outfile") || die "Unable to open output file '$outfile'.\n";
  while (<CARPETFILE>) {
    # The line is not a header comment
    if(/^[0-9]/) {
      my @data = split(/\s+/);
      if (not defined $data[$varindex + 2]) {
        print STDERR "\nInput file '$infile' has no values in data column " .
                     "$varindex.\n" . "Please specify a valid column index " .
                     "(starting from 0) !\n\n";
        exit -1;
      }
      print $fh $data[1]." ".$data[$varindex + 2]."\n";
    }
  }
  close ($fh);
  exit;
}


my %data;
my $time = -1;
my $new = 0;
my $lastit = -1;

my @datatoprint;
my $nsets = 1;
my $maxlength = 0;
my @lengths;
my @dataColumns = ();

$lengths[0]=0;
while (<CARPETFILE>)
{
  chomp;
  next if (/^$/);

  @dataColumns = split(/[ :]/, $1) if (/^# data columns: (.+)/);

  #Do nothing for headers!
  next if (/^#/);

  my @dataline = split(/[ \t]+/);
  my $currentit = $dataline[0];
  if ($currentit != $lastit) {
    if ($new) {
      # do not print "Time..." for zero-D data
      push(@datatoprint,"\n\n\#Time = $time\n") if ($direction !~ 8);

      my @sortedcoords = sort {$a <=> $b} (keys %data);
      foreach my $localcoord (@sortedcoords) {
        push(@datatoprint, $localcoord." ".$data{$localcoord}."\n");
      }
      $maxlength = $maxlength > (scalar @sortedcoords) ? $maxlength : (scalar @sortedcoords);
      $lengths[$nsets-1]=(scalar @sortedcoords);
      $nsets++;
      $lengths[$nsets-1]=0;
      %data=();
    }
    $new++;
    $time = $dataline[8];
    $lastit = $currentit;
  }
  my $coord = $dataline[$direction];
  if (not defined $dataline[12 + $varindex]) {
    print STDERR "\nInput file '$infile' has no values in data column " .
                 "$varindex.\n" . "Please specify a valid column index " .
                 "(starting from 0) !\n\n";
    exit -1;
  }
  my $val = $dataline[12 + $varindex];
  $data{$coord} = $val;
}

# do not print "Time..." for zero-D data
push(@datatoprint,"\n\n\#Time = $time\n") if ($direction !~ 8);

my @sortedcoords = sort {$a <=> $b} (keys %data);
foreach my $localcoord (@sortedcoords) {
  push(@datatoprint, $localcoord." ".$data{$localcoord}."\n");
}
$maxlength = $maxlength > (scalar @sortedcoords) ? $maxlength : (scalar @sortedcoords);
$lengths[$nsets-1]=(scalar @sortedcoords);
$nsets++;
$lengths[$nsets-1]=0;

my $outfile = $infile;
$outfile =~ s/asc$/xg/;
if (exists $dataColumns[$varindex*2 + 1]) {
  my @tokens = split(/::/, $outfile);
  if (@tokens) {
    $outfile = $tokens[0] . '::';
    @tokens = split(/\./, $tokens[1]);
  } else {
    @tokens = split('.', $outfile);
    $outfile = '';
  }
  $tokens[0] = $dataColumns[$varindex*2 + 1];
  $outfile .= join('.', @tokens);
}
print "  output to file '$outfile'\n\n";

my $fh = new FileHandle("> $outfile") or die "Unable to open output file '$outfile'.\n";

my $oldline="";
$nouts=0;
my $set=0;
foreach $line (@datatoprint) {
  if ($line =~ "Time") {
    if ($oldline) {
      for (my $i=$lengths[$set-1]; $i<$maxlength;$i++) {
        $nouts++;
        print $fh $oldline;
      }
    }
    $set++;
    print $fh $line;
  }
  else {
    $nouts++;
    print $fh $line;
    $oldline=$line
  }
}

# TR: I don't see why this is needed, it seems to print the last data item twice
# in CarpetIOBasic 0D output
#for (my $i=$lengths[$set-1]; $i<$maxlength;$i++) {
#  $nouts++;
#  print $fh $oldline;
#}
