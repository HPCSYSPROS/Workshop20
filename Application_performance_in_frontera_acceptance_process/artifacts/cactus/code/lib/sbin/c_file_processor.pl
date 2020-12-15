#!/bin/perl
#/*@@
#  @file      c_file_processor.pl
#  @date      Fri Jan 22 18:09:47 1999
#  @author    Tom Goodale / Gerd Lanfermann / Thomas Radke
#  @desc
#             Processes certain things within a C source file
#             which can't be dealt with by the normal C preprocessor.
#
#             This script puts everything after a DECLARE_CCTK macro
#             until the end of the routine into a new block.
#             It also fixes the function names for fortran wrappers.
#  @enddesc
#  @version   $Header$
#@@*/

$home = shift(@ARGV);

# Do we want line directives?
$line_directives = $line_directives eq 'yes';

$fortran_name_file = "$home/fortran_name.pl";

if (! -e "$fortran_name_file" )
{
  die "Unable to get fortran name file $fortran_name_file!";
}

require "$fortran_name_file";

$routine  = '';
$n_arg_left_braces = $n_arg_right_braces = 0;
$do_fix_fnames = 0;

# parse the file up to a ";\n"
$/ = ";\n";
#$* = 1; enable multiline support -- no longer supported in perl 5.10

while (<>)
{
  # split in lines... and collect in routine;
  foreach $mline (split ("\n"))
  {
    # skip one-line comments
    # (note that this is still incomplete for multi-line C comments -
    #  it is not checked if some code follows after the closing '*/')
    if ($mline !~ m/^\s*\/\//m && $mline !~ m/^\s*\/\*.*\*\/\s*$/m)
    {
      # Remove a ; from after the DECLARE_CCTK_* macros
      $mline =~ s/(DECLARE_CCTK_(PARAMETERS|ARGUMENTS))(\s*;)?/$1/m;

      # Remove a ; from after the fileversion macro
      # such a semicolon could lead to warning messages.
      $mline =~ s/^\s*(CCTK_FILEVERSION\s*\([^)]*\))(\s*;)?/$1/m;
      $mline =~ s/^\s*((ONE|TWO|THREE|FOUR|FIVE)_FORTSTRING_(CREATE|PTR)\s*\([^)]*\))(\s*;)?/$1/m;

      # start counting braces
      $n_arg_left_braces++  while ($mline =~ m/({)/gm);
      $n_arg_right_braces++ while ($mline =~ m/(})/gm);

      # check if we have to fix names of fortran wrappers
      $do_fix_fnames = 1 if ($mline =~ /(CCTK_FNAME|CCTK_FORTRAN_COMMON_NAME)/m);
    }

    $routine .= $mline . "\n";

    if ($n_arg_left_braces > 0 && $n_arg_left_braces - $n_arg_right_braces == 0)
    {
      $n_arg_left_braces = $n_arg_right_braces = 0;

      # call the fortran namefix routine/reset routine
      if ($do_fix_fnames)
      {
        fixfnames ($routine);
        $do_fix_fnames = 0;
      }
      else
      {
        print $routine;
      }
      $routine = '';
    }
  }
}

fixfnames ($routine);


sub fixfnames
{
  my $myroutine = shift (@_);
  @flines = split /(;)/m,$myroutine;

#  print $myroutine;

  foreach $fline (@flines)
  {
    while ($fline =~ m:CCTK_FNAME\s*\(([^\)]*)\):m)
    {
      $arglist = $1;
      $arglist =~ s:[\s\n\t]+::gm;

      @args = split(",", $arglist );

      $new = &fortran_name($args[$#args]);

      $fline =~ s:CCTK_FNAME\s*\(([^\)]*)\):$new:m;
    }
    while ($fline =~ m:CCTK_FORTRAN_COMMON_NAME\s*\(([^\)]*)\):m)
    {
      $arglist = $1;
      $arglist =~ s:[\s\n\t]+::gm;

      @args = split(",", $arglist );

      $new = &fortran_common_name($args[$#args]);

      $fline =~ s:CCTK_FORTRAN_COMMON_NAME\s*\(([^\)]*)\):$new:m;
    }

    print $fline;
  }
}
