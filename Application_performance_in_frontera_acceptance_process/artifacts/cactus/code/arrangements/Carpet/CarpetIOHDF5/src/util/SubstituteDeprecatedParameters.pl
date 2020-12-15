#! /usr/bin/perl -w
#/*@@
#  @file      SubstituteDeprecatedParameters.pl
#  @date      Fri 3 December 2004
#  @author    Thomas Radke
#  @desc
#             Perl script to automatically substitute deprecated parameter
#             names in a parfile
#  @enddesc
#@@*/

# Parse the command line
$parfile = shift;
if (! $parfile || shift)
{
  printf STDERR "This perl script automatically substitutes deprecated parameter names in a parfile.\n\n";
  printf STDERR "Usage: $0 <parameter file>\n\n";
  exit;
}

# build the hash table for renamed parameters
$renamed_parameters{"\UIOHDF5::out3D_dir\E"} = 'IOHDF5::out_dir';
$renamed_parameters{"\UIOHDF5::out3D_vars\E"} = 'IOHDF5::out_vars';
$renamed_parameters{"\UIOHDF5::out3D_extension\E"} = 'IOHDF5::out_extension';
$renamed_parameters{"\UIOHDF5::out3D_criterion\E"} = 'IOHDF5::out_criterion';
$renamed_parameters{"\UIOHDF5::out3D_every\E"} = 'IOHDF5::out_every';
$renamed_parameters{"\UIOHDF5::out3D_dt\E"} = 'IOHDF5::out_dt';
$renamed_parameters{"\UIOHDF5::in3D_dir\E"} = 'IOHDF5::in_dir';
$renamed_parameters{"\UIOHDF5::in3D_vars\E"} = 'IOHDF5::in_vars';
$renamed_parameters{"\UIOHDF5::in3D_extension\E"} = 'IOHDF5::in_extension';


open (PARFILE, $parfile) || die "Cannot open parameter file '$parfile' !";
print STDERR "Processing parameter file '$parfile'\n";
while (<PARFILE>)
{
  if (($pre, $parameter, $rest) = /^(\s*)(\w+::\w+)(.*)$/)
  {
    # rename old parameters if found
    if (defined $renamed_parameters{"\U$parameter\E"})
    {
      $_ = $pre . $renamed_parameters{"\U$parameter\E"} . $rest . "\n";
      print STDERR "  renaming parameter '$parameter' to '" . $renamed_parameters{"\U$parameter\E"} . "'\n";
    }
  }
  print;
}
close (PARFILE);
