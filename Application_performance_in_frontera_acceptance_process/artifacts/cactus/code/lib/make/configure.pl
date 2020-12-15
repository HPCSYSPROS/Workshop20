#! /usr/bin/perl -s
#/*@@
#  @file      configure.pl
#  @date      Fri Jan  8 15:06:22 1999
#  @author    Tom Goodale
#  @desc
#  Perl configuration script for the CCTK.
#  Does extra work it would be awkward to do from within the normal
#  autoconf stuff (or at least that I'm too lazy to do 8-)
#  @enddesc
#@@*/

$tmphome = shift(@ARGV);

print "Determining number of fortran underscores...\n";

($retcode,$data) = test_fortran_name();
push(@routines, $data);
# Some compilers do something really strange with common blocks.
($retcode,$data) = test_fortran_common_name();
push(@routines, $data);
push(@routines, "1;");

# Create the perl module to map the fortran names.
open(OUT, ">$tmphome/fortran_name.pl") || die "Cannot create fortran_name.pl\n";

foreach $line (@routines)
{
  print OUT "$line\n";
}

close(OUT);

if ($retcode > 0)
{
  print "Fortran compilation failed ...\n";
  print "COMPILATION WILL FAIL WITH ANY FORTRAN SOURCE CODE\n";
  print " ! Apparently a Fortran compiler was specified (F77/F90), but it does not \n".
        " ! seem to be working. Either make sure to specify a working Fortran compiler, \n".
        " ! do not set F77/F90, or set them to 'none'.\n\n"; 
}

sub test_fortran_name
{
  my($data);
  my($use_f77,$use_f90);
  my($retcode, $line, $name, $case, $n_underscores);
  my($underscore_suffix, $normal_suffix, $case_prefix);

  $use_f77 = 0;
  $use_f90 = 0;

  if ($compiler_f90 && $compiler_f90 ne "" && $compiler_f90 ne "none" &&
      $compiler_f90 !~ /no-fortran-compiler/)
  {
    ($retcode,$case, $n_underscores) = &compile_fortran_name($compiler_f90,$opts_f90);
    if ($retcode <= 0)
    {
      $use_f90 = 1;
    }
  }
  elsif($compiler_f77 && $compiler_f77 ne "" && $compiler_f77 ne "none" &&
        $compiler_f77 !~ /no-fortran-compiler/)
  {
    ($retcode,$case, $n_underscores) = &compile_fortran_name($compiler_f77,$opts_f77);
    if ($retcode <= 0)
    {
      $use_f77 = 1;
    }
  }

  if($use_f90 || $use_f77)
  {
    # Determine the case and number of underscores
    ($underscore_suffix, $normal_suffix, $case_prefix) = &determine_transformation($n_underscores, $case);

    $data = "
sub fortran_name
{
  my(\$old_name) = \@_;
  my(\$new_name);

  \$new_name = \"$case_prefix\$old_name\\E\";

  if(\$new_name =~ m:_: )
  {
    \$new_name = \$new_name.\"$underscore_suffix\";
  }
  else
  {
    \$new_name = \$new_name.\"$normal_suffix\";
  }

  return \$new_name;
}
";
  }
  else
  {
    if ($retcode <= 0)
    {
      print "No Fortran compiler - creating null fortran name conversion routine.\n";
    }
    else
    {
      print "Creating null fortran name conversion routine\n";
    }
$data = "
sub fortran_name
{
  my (\$old_name) = \@_;
  return \"\\L\$old_name\";
}
";
  }

  return ($retcode,$data);
}


sub test_fortran_common_name
{
  my($data);
  my($retcode, $line, $name, $case, $n_underscores);
  my($underscore_suffix, $normal_suffix, $case_prefix);

  $use_f77 = 0;
  $use_f90 = 0;

  if ($compiler_f90 && $compiler_f90 ne "" && $compiler_f90 !~ /none/ &&
      $compiler_f90 !~ /no-fortran-compiler/)
  {
    ($retcode,$case, $n_underscores) = &compile_fortran_common_name($compiler_f90,$opts_f90);
    if ($retcode<=0) {$use_f90 = 1};
  }
  elsif($compiler_f77 && $compiler_f77 ne "" && $compiler_f77 !~ /none/ &&
        $compiler_f77 !~ /no-fortran-compiler/)
  {
    ($retcode,$case, $n_underscores) = &compile_fortran_common_name($compiler_f77,$opts_f77);
    if ($retcode <=0) {$use_f77 = 1};
  }

  if($use_f90 || $use_f77)
  {
    # Determine the case and number of underscores
    ($underscore_suffix, $normal_suffix, $case_prefix) = &determine_transformation($n_underscores, $case);

    $data =  "
sub fortran_common_name
{
  my (\$old_name) = \@_;
  my (\$new_name);

  \$new_name = \"$case_prefix\$old_name\\E\";
  if(\$new_name =~ m:_: )
  {
      \$new_name = \$new_name.\"$underscore_suffix\";
  }
  else
  {
      \$new_name = \$new_name.\"$normal_suffix\";
  }

  return \"$prefix\".\$new_name;
}
";

  }
  else
  {
    if ($retcode <= 0)
    {
      print "No Fortran compiler - creating null fortran common name conversion routine.\n";
    }
    else
    {
      print "Creating null fortran common name conversion routine.\n";
    }
    $data = "

sub fortran_common_name
{
  my (\$old_name) = \@_;

  return \"\\L\$old_name\";
}
";
}

  return ($retcode,$data);
}

sub determine_transformation
{
  my ($n_underscores, $case) = @_;

  if($n_underscores == 0)
  {
    $normal_suffix = "";
    $underscore_suffix = "";
  }

  if($n_underscores == 1)
  {
    $normal_suffix = "_";
    $underscore_suffix = "_";
  }

  if($n_underscores == 2)
  {
    $normal_suffix = "_";
    $underscore_suffix = "__";
  }

  if($case == 0)
  {
    $case_prefix = "\\L";
  }
  if($case == 1)
  {
    $case_prefix = "\\U";
  }

  return ($underscore_suffix, $normal_suffix, $case_prefix);
}


sub compile_fortran_common_name
{
  my($compiler,$opts) = @_;
  my($data);
  my($retcode, $line, $name, $case, $n_underscores);
  my($underscore_suffix, $normal_suffix, $case_prefix);

  # Create a test file
  open(OUT, ">fname_test.f") || die "Cannot open fname_test.f\n";

  print OUT <<EOT;
      subroutine test_name
      real b
      common /test_common/ b
      b = 2.0
      end
EOT

  close OUT;

  # Compile the test file
  print "Compiling test file with $compiler $opts ...\n";
  system("$compiler $opts -c fname_test.f");

  $retcode = $? >> 8;

  if($retcode > 0)
  {
    print "Failed to compile fname_test.f\n";
  }
  else
  {

    # Search the object file for the appropriate symbols
    open(IN, "<fname_test.o") || open(IN, "<fname_test.obj") || die "Cannot open fname_test.o\n";

    $n_underscores = -1;
    while(<IN>)
    {
      my $line = $_;
      # This line may contain several matches
      while($line =~ m:(_[\w_]*)?(TEST_COMMON)(_*):i)
      {
	my $prefix = $1;
	my $name = $2;
	my $underscores = $3;
        my $nextline = $';   #'

	# This is a pain.  If all symbols have underscores, need to remove
	# the first one here.

	if($symbols_preceded_by_underscores)
	{
	  if($prefix =~ m:^_(.*):)
	  {
	    $prefix = $1;
	  }
	}

        my $tmp_case, $tmp_n_underscores;
	if($name =~ m:TEST_COMMON:)
	{
	  $tmp_case = 1;
	}
	if($name =~ m:test_common:)
	{
	  $tmp_case = 0;
	}
	if($underscores eq "")
	{
	  $tmp_n_underscores = 0;
	}
	if($underscores eq "_")
	{
	  $tmp_n_underscores = 1;
	}
	if($underscores eq "__")
	{
	  $tmp_n_underscores = 2;
	}

        # Look for the maximum number of underscores;
        # if the name occurs with fewer underscores it may be debug information
        if($tmp_n_underscores > $n_underscores)
        {
          $case = $tmp_case;
          $n_underscores = $tmp_n_underscores;
        }

        # Remove what currently matched; continue with the remaineder
        $line = $nextline;
      }
    }

    close IN;

    if($case == 1)
    {
      print "Uppercase - ";
    }
    elsif($case == 0)
    {
      print "Lowercase - ";
    }
    if($n_underscores == 0)
    {
      print "No trailing underscore\n";
    }
    elsif($n_underscores == 1)
    {
      print "One trailing underscore\n";
    }
    elsif($n_underscores == 2)
    {
      print "Two trailing underscores\n";
    }
  }
  # Delete the temporary files
  unlink <fname_test.*>;

  return ($retcode,$case,$n_underscores);
}


sub compile_fortran_name
{
  my($compiler,$opts) = @_;
  my($data);
  my($retcode, $line, $name, $case, $n_underscores);
  my($underscore_suffix, $normal_suffix, $case_prefix);

  # Create a test file
  open(OUT, ">fname_test.f") || die "Cannot open fname_test.f\n";

  print OUT <<EOT;
      subroutine test(a)
      integer a
      a = 1
      call test_name(a)
      end
EOT

  close OUT;

  # Compile the test file
  print "Compiling test file with $compiler_f77 $opts_f77 ...\n";
  system("$compiler_f77 $opts_f77 -c fname_test.f");

  $retcode = $? >> 8;

  if($retcode > 0)
  {
    print "Failed to compile fname_test.f\n";
  }
  else
  {

    # Search the object file for the appropriate symbols
    open(IN, "<fname_test.o") || open(IN, "<fname_test.obj") || die "Cannot open fname_test.o\n";

    $n_underscores = -1;
    while(<IN>)
    {
      my $line = $_;
      # This line may contain several matches
      while($line =~ m:(TEST_NAME)(_*):i)
      {
	my $name = $1;
	my $underscores = $2;
        my $nextline = $';   #'

	# Extremely quick hack to sort out problems later on with common block
	# names.
	
        my $tmp_symbols_preceded_by_underscores;
	if($_ =~ m:_TEST_NAME:i)
	{
	  $tmp_symbols_preceded_by_underscores=1;
	}
	else
	{
	  $tmp_symbols_preceded_by_underscores=0;
	}

	# Find out suffixes.
        my $tmp_case, $tmp_n_underscores;
	if($name =~ m:TEST_NAME:)
	{
	  $tmp_case = 1;
	}
	if($name =~ m:test_name:)
	{
	  $tmp_case = 0;
	}
	if($underscores eq "")
	{
	  $tmp_n_underscores = 0;
	}
	if($underscores eq "_")
	{
	  $tmp_n_underscores = 1;
	}
	if($underscores eq "__")
	{
	  $tmp_n_underscores = 2;
	}

        # Look for the maximum number of underscores;
        # if the name occurs with fewer underscores it may be debug information
        if($tmp_n_underscores > $n_underscores)
        {
          $symbols_preceded_by_underscores = $tmp_symbols_preceded_by_underscores;
          $case = $tmp_case;
          $n_underscores = $tmp_n_underscores;
        }

        # Remove what currently matched; continue with the remaineder
        $line = $nextline;
      }
    }

    close IN;

    if($case == 1)
    {
      print "Uppercase - ";
    }
    elsif($case == 0)
    {
      print "Lowercase - ";
    }
    if($n_underscores == 0)
    {
      print "No trailing underscore\n";
    }
    elsif($n_underscores == 1)
    {
      print "One trailing underscore\n";
    }
    elsif($n_underscores == 2)
    {
      print "Two trailing underscores\n";
    }

  }

  # Delete the temporary files
  unlink <fname_test.*>;

  return ($retcode,$case,$n_underscores);
}
