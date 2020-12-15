$top = `pwd` if (! $top);
$config_dir = "$top/config-data" if (! $config_dir);

# Set up the CCTK home directory
if(! $cctk_home)
{
  $cctk_home = $ENV{'CCTK_HOME'} || "$ENV{HOME}/CCTK";
  $cctk_home =~ s:/$::g;
}

$sbin_dir = "$cctk_home/lib/sbin";

die "Unable to find CCTK sbin directory - tried $sbin_dir\n"
  if (! -e "$sbin_dir/parameter_parser.pl");

require "$sbin_dir/CSTUtils.pl";


############################################################
=pod 

=over 

=item Configure($config, $home_dir, $prompt)
 Sets $config_data with $home_dir and $prompt.

=back

=cut

############################################################
sub Configure
{
  my($config,$home_dir,$prompt) = @_;
  my($configs_dir,$tests_dir);

  # Cactus home directory
  $config_data->{"CCTK_DIR"} = $home_dir;

  # Interactive of not
  $config_data->{"PROMPT"} = $prompt;

  # Cactus configurations directory
  if($ENV{"CONFIGS_DIR"})
  {
    $configs_dir = $ENV{"CONFIGS_DIR"};
  }
  else
  {
    $configs_dir = "$homedir/configs";
  }
  $config_data->{"CONFIGSDIR"} = $configs_dir;

  # Cactus test directory
  if ($ENV{"TESTS_DIR"})
  {
    $tests_dir = $ENV{"TESTS_DIR"};
  }
  else
  {
    $tests_dir = $home_dir."/TEST";
  }
  $config_data->{"TESTS_DIR"} = $tests_dir;

  $config_data->{"SEPARATOR"} = "/";
  $config_data->{"CONFIG"} = $config;

  # Get the executable
  $config_data = FindExecutionDetails($config_data);

  return $config_data;
}


############################################################
=pod 

=over 

=item MissingThorns($active, $allthorns)
 Counts missing thorns and assembles a list of those
 missing. 

=back

=cut

############################################################
sub MissingThorns
{
  my($active,$allthorns) = @_;
  my(@at,$th,$foundit,$thornpart,$nmissing,$missing);

  @at = split(' ',$active);
  $nmissing = 0;
  $missing = "";
  foreach $th (@at)
  {
    $th = "\U$th";
    $foundit = 0;

    foreach $tthorn  (split(" ",$allthorns))
    {
      $thornpart = "\U$tthorn";
      if ($thornpart eq $th)
      {
        $foundit = 1;
      }
    }
    if (!$foundit)
    {
      $missing .=" $th";
      $nmissing++;
    }
  }

  return($nmissing,$missing);

}

############################################################
=pod 

=over 

=item ParseParFile($thorn, $arrangement, $parfile, $config_data)
 This subroutine parses the parameter file for active
 thorns and returns a list of those thorns and their
 descriptions.

=back

=cut

############################################################
sub ParseParFile
{
  my($thorn,$arrangement,$parfile,$config_data) = @_;
  my($line,$file,$processing_active);
  my($active,$desc);

  $file =  "arrangements/$arrangement/$thorn/test/$parfile";

  open (IN, "<$file") || die "Can not open $file";

  $processing_active = 0;

  # Give a default test name in case none is specified in the parameter file.
  $desc = "$arrangement/$thorn/test/$parfile";

  $active = "";

  while (<IN>)
  {
    $line = $_;
    $line =~ s/#.*//g;

    if($processing_active == 1)
    {
      if($line =~ m/(.*)\"/)
      {
        $active .= " ".$1;
        $processing_active = 0;
      }
      else
      {
        $active .= " ".$line;
      }
    }
    elsif ($line =~ m/^\s*\!\s*DESC(RIPTION)?\s*\"(.*)\"\s*$/i)
    {
      $desc = $2;
    }
    elsif ($line =~ m/^\s*ActiveThorns\s*=\s*\"(.*)\"/i)
    {
      $active .= " ".$1;
    }
    elsif($line =~ m/^\s*ActiveThorns\s*=\s*\"(.*)/i)
    {
      $active .= " ".$1;
      $processing_active = 1;
    }
  }
  close IN;

  return($active,$desc);

}

############################################################
=pod 

=over 

=item ParseTestConfigs($testdata, $config_data, $rundata)
 Parses the test.ccl files for absolute
 and relative tolerance and nprocs used. 

=back

=cut

############################################################
sub ParseTestConfigs
{
  my($testdata,$config_data,$rundata) = @_;
  my($line_number, $line);
  my($test, $ABSTOL, $RELTOL);

  my $arrangement_dir = "$config_data->{'CCTK_DIR'}${sep}arrangements${sep}";
  foreach $thorn (split(" ",$testdata->{"THORNS"}))
  {
    my $testdir = $arrangement_dir . $testdata->{"$thorn ARRANGEMENT"} .
                  "${sep}$thorn${sep}test";
    my $config_file     = "$testdir${sep}test.ccl";
    my $old_config_file = "$testdir${sep}config";

    if (-r $old_config_file)
    {
      print "  Thorn $thorn has a 'config' file in 'test/'.\n".
            "    Config files are deprecated. Use test.ccl instead.\n";
    }
    if (-r $config_file)
    {
      my @config = &read_file($config_file);
      for($line_number = 0; $line_number < @config; $line_number++)
      {
        $line = $config[$line_number];

        my @insorder_thornabstol_keys, @insorder_thornreltol_keys;
        my $insorder_thornabstol_counter=0; 
        my $insorder_thornreltol_counter=0;
        # Parse tokens
        if ($line =~ m/^\s*ABSTOL\s*(\S*)\s*(\S*)\s*$/i)
        {
          my $newtol=$1;
          my $varRegex=$2;
          if ( $varRegex =~ m/^$/i ) {
             # No regular expression given. Setting default tolerance ".*"
             $varRegex=".*";
          }
          $rundata->{"$thorn ABSTOL"}{$varRegex}=$newtol;
          $ABSTOL=$$rundata{"$thorn ABSTOL"};
        }
        elsif ($line =~ m/^\s*RELTOL\s*(\S*)\s*(\S*)\s*$/i)
        {
          my $newtol=$1;
          my $varRegex=$2;
          if ( $varRegex =~ m/^$/i ) {
             # No regular expression given. Setting default tolerance ".*"
             $varRegex=".*";
          }
          $rundata->{"$thorn RELTOL"}{$varRegex}= $newtol;
          $RELTOL=$$rundata{"$thorn RELTOL"};
        }
        elsif ($line =~ m/^\s*NPROCS\s+(\d+)\s*$/i)
        {
          $NPROCS = $rundata->{"$thorn NPROCS"} = $1;
        }
        elsif ($line =~ m/^\s*EXTENSIONS\s*(.*)/i)
        {
          $testdata->{"EXTENSIONS"} .= "$1 ";
        }
        elsif ($line =~ m/^\s*TEST\s*(.*)/i)
        {
          ($test, $ABSTOL, $RELTOL, $NPROCS, $line_number) =
            &ParseTestBlock($line_number, \@config);
          $rundata->{"$thorn $test ABSTOL"} = $ABSTOL;
          $rundata->{"$thorn $test RELTOL"} = $RELTOL;
          $rundata->{"$thorn $test NPROCS"} = $NPROCS;
        }
        else
        {
          print "  Unrecognised token $line in config file for thorn $thorn\n";
        }
      }
    }
  }

  return $testdata;
}

############################################################
=pod 

=over 

=item ParseTestBlock($line_number, $data)
 This subroutine parses for ABSTOL, RELTOL, and NPROCS.

=back

=cut

############################################################
sub ParseTestBlock
{
  my ($line_number, $data) = @_;
  my ($Test, $NPROCS) = ();
  my (%ABSTOL, %RELTOL) = (); 

  $data->[$line_number] =~ m/^\s*PROVIDES\s*(.*)/i;

  $Test = $1;

  $line_number++;

  if($data->[$line_number] !~ m/^\s*\{\s*$/)
  {
    $line_number++ while($data[$line_number] !~ m:\s*\}\s*:);
  }
  else
  {
    while($data->[$line_number] !~ m:\s*\}\s*:)
    {
      $line_number++;
      if ($data->[$line_number] =~ m/^\s*ABSTOL\s*(\S*)\s*(\S*)\s*$/i)
      {
        my $newtol=$1;
        my $varRegex=$2;
        if ( $varRegex =~ m/^$/i ) {
           # No regular expression given. Setting default tolerance ".*"
           $varRegex=".*";
        }
        $ABSTOL{$varRegex} = $newtol;
        next;
      }
      elsif ($data->[$line_number] =~ m/^\s*RELTOL\s*(\S*)\s*(\S*)\s*$/i)
      {
        my $newtol=$1;
        my $varRegex=$2;
        if ( $varRegex =~ m/^$/i ) {
           # No regular expression given. Setting default tolerance ".*"
           $varRegex=".*";
        }
        $RELTOL{$varRegex} = $newtol;
        next;
      }
      elsif ($data->[$line_number] =~ m/^\s*NPROCS\s+(\d+)\s*$/i)
      {
        $NPROCS = $1;
        next;
      }
      elsif($data->[$line_number] =~ m:\s*\}\s*:)
      {
        # do nothing.
      }
      else
      {
        print STDERR "Error parsing test config block line '$data->[$line_number]'\n";
      }
    }
  }
  return ($Test, \%ABSTOL, \%RELTOL, $NPROCS, $line_number);
}

############################################################
=pod 

=over 

=item Configure()
 Finds archive files for a test if there are any and returns
 the location of the test results.

=back

=cut

############################################################
sub FindTestArchiveFiles
{
  my($test,$thorn,$testdata) = @_;

  my $dir = "$testdata->{\"$thorn TESTSDIR\"}/$test";
  if ( !-d "$dir")
  {
    my $tarname      = undef;
    my $unpackoption = undef;
    if (-f "$dir.tar")     { $tarname ="$dir.tar";     $unpackoption = "";}
    if (-f "$dir.tar.gz")  { $tarname ="$dir.tar.gz";  $unpackoption = "--gzip";}
    if (-f "$dir.tgz")     { $tarname ="$dir.tgz";     $unpackoption = "--gzip";}
    if (-f "$dir.tar.bz2") { $tarname ="$dir.tar.bz2"; $unpackoption = "--bzip";}
    if (-f "$dir.tbz")     { $tarname ="$dir.tbz";     $unpackoption = "--bzip";}
    if (-f "$dir.tar.xz")  { $tarname ="$dir.tar.xz";  $unpackoption = "--xz";}
    if (-f "$dir.txz")     { $tarname ="$dir.txz";     $unpackoption = "--xz";}
    if (defined($unpackoption))
    {
      my $tar       = `gtar --help > /dev/null 2> /dev/null && echo -n gtar || echo -n tar`;
      my $supported = `tar --help | grep -- "$unpackoption" > /dev/null 2> /dev/null && echo -n "yes" || echo -n "no"`;
      if ($supported ne "yes")
      {
        print("  $tarname found, but $tar does not support the option $unpackoption.");
        return $testdata;
      }
      print("  Decompressing testsuite data\n");
      my $basename = `basename $dir`; chomp($basename);
      my $unpackdir = $config_data->{"TESTS_DIR"}.$sep.$config_data->{"CONFIG"}.$sep.$thorn;
      system("mkdir -p $unpackdir");
      system("cd $unpackdir && rm -rf $basename $basename.orig");
      system("cd $unpackdir && $tar -x $unpackoption -f $tarname && mv $basename $basename.orig");
      $testdata->{"$thorn $test UNCOMPRESSED"} = "$unpackdir$sep$basename.orig";
      $dir = "$unpackdir$sep$basename.orig";
    }
  }
  ($testdata->{"$thorn $test UNKNOWNFILES"},$testdata->{"$thorn $test DATAFILES"}) = &FindFiles($dir,$testdata);
  $testdata->{"$thorn $test NDATAFILES"} = scalar(split(" ",$testdata->{"$thorn $test DATAFILES"}));
  return $testdata;
}


############################################################
=pod 

=over 

=item FindTestParameterFiles($testdata, $config_data)
 For each thorn, finds the parameter files residing in
 the test directory.

=back

=cut

############################################################
sub FindTestParameterFiles
{
  my($testdata,$config_data) = @_;
  my($config,$config_dir);
  my($thorn);
  my(%found_thorns) = ();

  $config      = $config_data->{"CONFIG"};
  $configs_dir = $config_data->{"CONFIGSDIR"};
  $sep         = $config_data->{"SEPARATOR"};

  open (AT, "< $configs_dir${sep}$config${sep}ThornList") || print "Cannot find ThornList for $config";

  while (<AT>)
  {
    next if (/^\s*(\#.*|\!.*)$/);

    /^\s*([^\s]*)\s*/;

    $fullthorn = $1;
    next if (! $fullthorn);

    $fullthorn =~ m:^\s*([^\s]*)/([^\s]*)\s*:;

    $arrangement = $1;
    $thorn = $2;

    # skip duplicate entries in the ThornList
    next if (defined $found_thorns{"$arrangement/$thorn"});
    $found_thorns{"$arrangement/$thorn"} = 1;

    $testdata->{"FULL"} .= "$fullthorn ";
    $testdata->{"THORNS"} .= "$thorn ";
    $testdata->{"$thorn ARRANGEMENT"} = "$arrangement";

    if ($testdata->{"ARRANGEMENTS"} !~ m:\s$arrangement\s:)
    {
      $testdata->{"ARRANGEMENTS"} .= "$arrangement ";
    }

    $thorntestdir = "$config_data->{\"CCTK_DIR\"}${sep}arrangements${sep}$fullthorn${sep}test";

    if (-d $thorntestdir)
    {
      $testdata->{"$thorn TESTSDIR"} = $thorntestdir;

      chdir $thorntestdir;

      while ($file=<*.par>)
      {
        $file =~ m:^(.*)\.par$:;
        $filedir = $1;
        if (-d $filedir or -f "$filedir.tar" or
            -f "$filedir.tar.gz"  or -f "$filedir.tgz" or
            -f "$filedir.tar.bz2" or -f "$filedir.tbz" or
            -f "$filedir.tar.xz"  or -f "$filedir.txz")
        {
          $testdata->{"$thorn TESTS"} .= "$filedir ";
          $testdata->{"$thorn NTESTS"}++;
        }
        else
        {
          print "Parameter file $filedir in thorn $thorn but no output directory\n";
        }
      }
    }
  }
  chdir $config_data->{"CCTK_DIR"};

  close AT;

  return $testdata;
}

############################################################
=pod 

=over 

=item FindExecutionDetails($config_data)
 Finds the executable to be used. 

=back

=cut

############################################################
sub FindExecutionDetails
{
  my($config_data) = @_;
  my($config,$dir,$sep,$defns,$defexename,$executable);

  $config = $config_data->{"CONFIG"};
  $sep = $config_data->{"SEPARATOR"};

  # Check the name and directory of executable
  $defns = "$config_data->{\"CONFIGSDIR\"}${sep}$config${sep}config_data${sep}make.config.defn";

  $defexename = "cactus_$config";

  if (-e "$defns")
  {
    open(DEFNS,"<$defns");
    while(<DEFNS>)
    {
      if (/EXE\s*=\s*([\w_-]+)/)
      {
        $defexename = $1;
      }
      if (/EXEDIR\s*=\s*([\w_-]+)/)
      {
        $defexedirname = $1;
      }
    }
    close(DEFNS);
  }

  $executable = &defprompt('  Enter executable name ($exe, relative to Cactus home dir)',"exe$sep$defexename");

  if (! (-e "$config_data->{\"CCTK_DIR\"}$sep$executable"))
  {
    if (-e "$config_data->{\"CCTK_DIR\"}$sep$executable.exe")
    {
      $executable .= ".exe";
    }
    else
    {
      die "Cannot locate executable '$executable'\n";
    }
  }

  $config_data->{"EXE"} = "$config_data->{\"CCTK_DIR\"}$sep$executable";

  $config_data = &FindRunCommand($config_data);

  return $config_data;
}

############################################################
=pod 

=over 

=item FindRunCommand($config_data)
 Sets the mpirun command to run the testsuite. 

=over

=item Run Command Form
 mpirun -np $nprocs $exe $parfile

=back

=back

=cut

############################################################
sub FindRunCommand
{
  my($config_data) = @_;

  # Look to see if MPI is dfined
  my $have_mpi = ParseExtras($config_data);

  my $nprocs = (defined ($ENV{'CCTK_TESTSUITE_RUN_PROCESSORS'}) ?
                $ENV{'CCTK_TESTSUITE_RUN_PROCESSORS'} : $have_mpi ?
                                                        2 : 1);

  if ($have_mpi)
  {
    $config_data->{'NPROCS'} =
      &defprompt('  Enter number of processors ($nprocs)', $nprocs);
  }
  else
  {
      print "No MPI available\n";
      if ($nprocs > 1)
      {
          die("Cannot run on $nprocs processes without an MPI implementation\n");
      }
      $config_data->{'NPROCS'} = 1;
  }

  my $command;
  if (defined ($ENV{'CCTK_TESTSUITE_RUN_COMMAND'}))
  {
    $command = $ENV{'CCTK_TESTSUITE_RUN_COMMAND'};
  }
  else
  {
   $command = $have_mpi ? 'mpirun -np $nprocs $exe $parfile' : '$exe $parfile';
  }
  $config_data->{'COMMAND'} =
    &defprompt("  Enter command to run testsuite", $command);

  return $config_data;
}

############################################################
=pod 

=over 

=item defprompt($pr, $de)
 Handles prompts.

=back

=cut

############################################################
sub defprompt
{
  my ($pr, $de) = @_;
  my ($res);

  if ($config_data->{"PROMPT"} eq "no")
  {
    $res = $de;
  }
  else
  {
    print "$pr [$de] \n";
    print "   --> ";

    $res = <STDIN> if ($prompt eq "yes");
    if ($res =~ m/^$/)
    {
      $res = $de;
    }
    elsif ($res =~ m/^ $/)
    {
      $res = "";
    }
    $res =~ s/\n//;
    print "\n";
  }
  return $res;
}

############################################################
=pod 

=over 

=item ParseExtras($config_data)
 Parses extra definitions present in 
 config-data/cctk_Extradefs.h

=back

=cut

############################################################
sub ParseExtras
{
  my($config_data) = @_;
  my($mpi,$dir,$sep,$extradir,$capabilitydir);

  $sep = $config_data->{"SEPARATOR"};
  $config = $config_data->{"CONFIG"};
  if(defined($ENV{CACTUS_CONFIGS_DIR})) {
    $dir = $ENV{CACTUS_CONFIGS_DIR};
  } else {
    $dir = $config_data->{"CCTK_DIR"}.${sep}."configs";
  }

  $extradir = "${dir}${sep}$config{$sep}config-data${sep}cctk_Extradefs.h";
  $capabilitydir = "${dir}${sep}$config${sep}bindings${sep}Configuration${sep}Capabilities${sep}cctki_MPI.h";

  $mpi = 0;

  if (-e "$extradir")
  {
    open(EXTRA,"<$extradir");
    while(<EXTRA>)
    {
      if (/\#define CCTK_MPI/)
      {
        $mpi = 1;
      }
    }
    close(EXTRA);
  }

  if (-e "$capabilitydir")
  {
    open(CAP,"<$capabilitydir");
    while(<CAP>)
    {
      if (/\#define CCTK_MPI/)
      {
        $mpi = 1;
      }
    }
    close(CAP);
  }

  return $mpi;
}

############################################################
=pod 

=over 

=item InitialiseTestData()
 Sets $testdata to its null values.

=over 

=item $testdata
 Consists of FULL, NNODATAFILES, NRUNNABLE, NUNRUNNABLE,
 RUNNABLETHORNS, UNRUNNABLETHORNS, RUNNABLEARRANGEMENTS,
 and UNRUNNABLEARRANGEMENTS.

=back

=back

=cut

############################################################
sub InitialiseTestData
{
  my($testdata);

  # Complete list of thorns: arrangement/thorn
  $testdata->{"FULL"} = "";

  $testdata->{"NNODATAFILES"} = 0;
  $testdata->{"NRUNNABLE"} = 0;
  $testdata->{"NUNRUNNABLE"} = 0;
  $testdata->{"RUNNABLETHORNS"} = "";
  $testdata->{"UNRUNNABLETHORNS"} = "";
  $testdata->{"RUNNABLEARRANGEMENTS"} = "";
  $testdata->{"UNRUNNABLEARRANGEMENTS"} = "";

  return $testdata;
}

############################################################
=pod 

=over 

=item InitialiseRunData()
 This subroutine sets $runconfig data.

=over

=item $runconfig
 $runconfig consists of ABSTOL (default 1e-12) and
 RELTOL (1e-12).

=back

=back

=cut

############################################################
sub InitialiseRunData
{
  my(%runconfig);

  $runconfig{"ABSTOL"} = 1e-12;
  $runconfig{"RELTOL"} = 1e-12;

  return %runconfig;
}

############################################################
=pod 

=over 

=item PrintDataBase($database)
 Prints any $database.

=back

=cut

############################################################
sub PrintDataBase
{
  my($database) = @_;
  my($field);

  foreach $field ( sort keys %$database )
  {
    # If field is a hash table reference, descend into hash
    if ( ref $database->{$field} eq "HASH" )
    {
       my $subhashfield;
       my $reftoSubhash=$$database{$field};
       my $hashSize = scalar( keys %$reftoSubhash );
       if ( $hashSize == 0 )
       {
         print "$field is an empty hash.\n";
         next;
       }
       print "$field is a hash.  Descending into it.\n";
       foreach $subhashfield ( sort keys %$reftoSubhash )
       {
         print "    Field $subhashfield has value \n       $database->{$field}->{$subhashfield}\n";
       }
    }
    else
    {
       print "$field has value\n   $database->{$field}\n";
    }
  }
}

############################################################
=pod 

=over 

=item PrintToleranceTable($test, $thorn, $testdata, $runconfig)
 Prints the ABSTOL and RELTOL for the given thorn.

=back

=cut

############################################################
sub PrintToleranceTable
{
  my($test,$thorn,$testdata,$runconfig) = @_;
  my($fileabstol,$filereltol,$maxfilenamelen);

  # Get default tolerances for the test
  if (defined($runconfig{"$thorn $test ABSTOL"}->{".*"}))
  {
     $testabstol=$runconfig{"$thorn $test ABSTOL"}->{".*"};
  }
  elsif (defined($runconfig{"$thorn ABSTOL"}->{".*"}) )
  {
     $testabstol=$runconfig{"$thorn ABSTOL"}->{".*"};
  }
  else
  {
     $testabstol=$runconfig{"ABSTOL"};
  }

  if (defined($runconfig{"$thorn $test RELTOL"}->{".*"}))
  {
     $testreltol=$runconfig{"$thorn $test RELTOL"}->{".*"};
  }
  elsif (defined($runconfig{"$thorn RELTOL"}->{".*"}) )
  { 
     $testreltol=$runconfig{"$thorn RELTOL"}->{".*"};
  }
  else
  {
     $testreltol=$runconfig{"RELTOL"};
  }

  # longest file name for table alignment
  $maxfilenamelen = length("(.*)");
  foreach $file (split(" ",$testdata->{"$thorn $test DATAFILES"}))
  {
     $maxfilenamelen = length($file) if ($maxfilenamelen < length($file));
  }

  # Print test's default tolerances and any deviations
  print "------------------------------------------------------------------------\n\n";
  print "  Test $thorn: $test \n";
  print "    \"$testdata->{\"$thorn $test DESC\"}\"\n";

  print "    File"," "x($maxfilenamelen-4),"\tAbs Tol\t\tRel Tol\n";
  print "    --------------------------------------------------------------------\n";
  print "    (.*)"," "x($maxfilenamelen-4),"\t$testabstol\t\t$testreltol\n";
  foreach $file (split(" ",$testdata->{"$thorn $test DATAFILES"}))
  {
     ($fileabstol, $filereltol)=&GetFileTolerances($test,$thorn,\%runconfig,$file);
     if ( $fileabstol == $testabstol ) { $fileabstol="--"; }
     if ( $filereltol == $testreltol ) { $filereltol="--"; }
     if ( $fileabstol ne "--" || $filereltol ne "--" )
     {
        print "    $file"," "x($maxfilenamelen-length($file)),"\t$fileabstol\t\t$filereltol\n";
     }
  }
  print "\n";
}

############################################################
=pod 

=over 

=item GetToleranceFiles()
 Gets the relative and absolute tolerances from test.ccl (?).

=back

=cut

############################################################
sub GetFileTolerances
{
  my($test,$thorn,$runconfig,$file) = @_;

  # Initialize to global tolerances
  my $fileabstol = $runconfig->{"ABSTOL"};
  my $filereltol = $runconfig->{"RELTOL"};

  # References to tolerance hashes
  my $refAbsTolThorn=$$runconfig{"$thorn ABSTOL"};
  my $refRelTolThorn=$$runconfig{"$thorn RELTOL"};
  my $refAbsTolTest=$$runconfig{"$thorn $test ABSTOL"};
  my $refRelTolTest=$$runconfig{"$thorn $test RELTOL"};

  # Cycle through thorn-wide tolerance hashes for file-matching expressions
  my $countmatches;
  if ( scalar( keys %$refAbsTolThorn ) )
  {
     $countmatches=0;
     foreach $varRegex ( sort keys %$refAbsTolThorn )
     {
        if ( $file =~ m/$varRegex/i )
        {
           $fileabstol=$refAbsTolThorn->{$varRegex};
           if ( $varRegex ne ".*" ) { $countmatches++; }
        }
        if ( $countmatches > 1 ) {
           warn "ERROR: Data file $file of thorn $thorn matches more than one regular expression for ABSTOL.\n";
           die  "ABORTING: Please adjust test.ccl of $thorn to avoid multiple matches at the thorn-wide tolerance level (ABSTOL).\n"
        }
     }
  }
  if ( scalar( keys %$refRelTolThorn ) )
  {
     $countmatches=0;
     foreach $varRegex ( sort keys %$refRelTolThorn )
     {
        if ( $file =~ m/$varRegex/i )
        { 
           $filereltol=$refRelTolThorn->{$varRegex}; 
           if ( $varRegex ne ".*" ) { $countmatches++; }
        }
        if ( $countmatches > 1 ) {
           warn "ERROR: Data file $file of thorn $thorn matches more than one regular expression for RELTOL.\n";
           die  "ABORTING: Please adjust test.ccl of $thorn to avoid multiple matches at the thorn-wide tolerance level (RELTOL).\n"
        }
     }
  }

  # Cycle through test-specific tolerance hashes for file-matching expressions
  if ( scalar( keys %$refAbsTolTest ) )
  {
     $countmatches=0;
     foreach $varRegex ( sort keys %$refAbsTolTest )
     {
        if ( $file =~ m/$varRegex/i ) 
            {
           $fileabstol=$refAbsTolTest->{$varRegex}; 
           if ( $varRegex ne ".*" ) { $countmatches++; }
        }
        if ( $countmatches > 1 ) {
           warn "ERROR: Data file $file of test $thorn -- $test matches more than one regular expression for ABSTOL.\n";
           die  "ABORTING: Please adjust test.ccl of $thorn to avoid multiple matches at the $test tolerance level (ABSTOL).\n"
        }
     }
  }
  if ( scalar( keys %$refRelTolTest ) )
  {
     $countmatches=0;
     foreach $varRegex ( sort keys %$refRelTolTest )
     {
        if ( $file =~ m/$varRegex/i ) 
        {
           $filereltol=$refRelTolTest->{$varRegex};
           if ( $varRegex ne ".*" ) { $countmatches++; }
        }
        if ( $countmatches > 1 ) {
           warn "ERROR: Data file $file of test $thorn -- $test matches more than one regular expression for RELTOL.\n";
           die  "ABORTING: Please adjust test.ccl of $thorn to avoid multiple matches at the $test tolerance level (RELTOL).\n"
        }
     }
  }

  return ($fileabstol, $filereltol);

}

############################################################
=pod 

=over 

=item CleanDir($dir)
 Deletes all files in $dir.

=back

=cut

############################################################
sub CleanDir
{
  my($dir) = @_;

  opendir (DIR, $dir);
  @list = (grep (/.+\..+/, readdir (DIR)));
  foreach $entry (@list)
  {
    unlink "$dir/$entry";
  }
  closedir (DIR);
}

############################################################
=pod 

=over 

=item RunCactus($output, $testname, $command)
 Runs Cactus and returns $retcode.

=back

=cut

############################################################
sub RunCactus
{
  my($output,$testname,$command) = @_;
  my($retcode);

  printf "\n  Issuing $command\n";

  $retcode = 0;
  open (CMD, "pwd; $command 2>&1 |");
  open (LOG, "> $testname.log");

  while (<CMD>)
  {
    print LOG if ($output =~ /log/);
    print STDOUT if ($output =~ /stdout/);

    if( /Cactus exiting with return code (.*)/)
    {
      $retcode = $1 + 0;
    }
  }
  close LOG;
  close CMD;

  $retcode = $? >> 8 if($retcode==0);

  print STDOUT "\n\n" if ($output =~ /stdout/);

  return $retcode;
}


############################################################
=pod 

=over 

=item fpabs($val)
 Returns absolute value of $val.

=back

=cut

############################################################
sub fpabs
{
  my ($val) = $_[0];
  $val > 0 ? $val:-$val;
}


############################################################
=pod 

=over 

=item PrintHeader()
 Prints the header.

=back

=cut

############################################################
sub PrintHeader
{
  print <<EOT;

  -------------------------------
   Cactus Code Test Suite Tool
  -------------------------------

EOT
}


############################################################
=pod 

=over 

=item FindFiles($dir, $testdata)
 Finds files by extension in $testdata->{"EXTENSIONS"}
 and returns $unrecognizedfiles and $recognizedfiles.

=back

=cut

############################################################
sub FindFiles
{
  my ($dir,$testdata) = @_;
  my ($unrecognizedfiles,$recognizedfiles, @tmp);

  $recognizedfiles="";
  $unrecognizedfiles="";

  opendir (DIR, $dir);
  @tmp = readdir (DIR);
  closedir (DIR);

  foreach $f (@tmp)
  {
    $f =~ m:.*\.([^\s\.]+)\s*$:;
    $extension = $1;

    if ($f !~ /^(\.\#.*|\.|\.\.|.*\.par|CVS|.svn|.*~)$/)
    {
      if ($extension =~ /.+/ && $testdata->{"EXTENSIONS"} =~ /\b$extension\b/)
      {
        $recognizedfiles .= " $f ";
      }
      else
      {
        $unrecognizedfiles .= " $f";
      }
    }
  }

  return ($unrecognizedfiles,$recognizedfiles);
}


############################################################
=pod 

=over 

=item WriteFullResults($rundata, $testdata, $config_data)
 This subroutine writes a summary of the test results.

=back

=cut

############################################################
sub WriteFullResults
{
  my ($rundata,$testdata,$config_data) = @_;
  my @summary = ();
  my ($separator);

  $separator1 = "========================================================================\n\n";
  $separator2 = "------------------------------------------------------------------------\n\n";

  push (@summary, $separator2);
  push (@summary, "  Warnings for configuration $config_data->{\"CONFIG\"}\n  --------\n\n");

  # Missing thorns for tests

  $message = "  Tests missed for lack of thorns:\n";
  $missingtests = 0;
  foreach $thorn (sort split(' ',$testdata->{'THORNS'}))
  {
    foreach $parfile (sort split(' ',$testdata->{"$thorn TESTS"}))
    {
      my $missing = $testdata->{"$thorn $parfile MISSING"};
      next unless ($missing);
      $message .= "\n    ".$parfile." in ". $thorn."\n";
      $message .= "      (". $testdata->{"$thorn $parfile DESC"}.")\n";
      $message .= "      Missing: $missing\n";
      $missingtests++;
    }
  }
  push (@summary, "$message\n") if ($missingtests > 0);

  # Different number of processors required
  $message = "  Tests missed for different number of processors required:\n";
  $missingtests = 0;
  foreach $thorn (sort split(' ',$testdata->{'THORNS'}))
  {
    foreach $parfile (sort split(' ',$testdata->{"$thorn TESTS"}))
    {
      my $nprocs = $testdata->{"$thorn $parfile NPROCS"};
      next unless ($nprocs);
      $message .= "\n    ".$parfile." in ". $thorn."\n";
      $message .= "      (". $testdata->{"$thorn $parfile DESC"}.")\n";
      $message .= "      Requires $nprocs processors\n";
      $missingtests++;
    }
  }
  push (@summary, "$message\n") if ($missingtests > 0);

  # Different numbers of test files

  $message =  "  Tests with different number of test files:\n\n";

  $extratests = 0;
  foreach $thorn (sort split(" ",$testdata->{"RUNNABLETHORNS"}))
  {
    foreach $parfile (sort split(" ",$testdata->{"$thorn RUNNABLE"}))
    {
      if ($rundata->{"$thorn $parfile NFILEEXTRA"}>0)
      {
        $extratests++;
        $message .= "    $thorn ($parfile)\n";
        $message .= "      Test created $rundata->{\"$thorn $parfile NFILEEXTRA\"} extra files: $testdata->{\"$thorn $parfile FILEEXTRA\"}\n";
      }
    }
  }
  push (@summary, "$message\n") if ($extratests > 0);

  push (@summary, $separator2);

  push (@summary, "  Testsuite Summary for configuration $config_data->{\"CONFIG\"}\n");
  push (@summary, "  -----------------\n\n");

  push (@summary, "  Suitable testsuite parameter files found in:\n\n");

  $tested = 0;
  $nottested = "";
  foreach $thorn (sort split(" ",$testdata->{"THORNS"}))
  {
    $num = scalar(split(" ",$testdata->{"$thorn RUNNABLE"}));
    if ($num > 0)
    {
      push (@summary, "    $thorn [$num]\n");
      $tested++;
    }
    else
    {
      $nottested .= "\n    $thorn";
    }
  }

  push (@summary, "\n");
  push (@summary, "  Details:\n\n");
  foreach $thorn (sort split(" ",$testdata->{"THORNS"}))
  {
    $num = scalar(split(" ",$testdata->{"$thorn RUNNABLE"}));
    if ($num > 0)
    {
      push (@summary, "    $thorn:\n");
      foreach $test (sort split(" ",$testdata->{"$thorn RUNNABLE"}))
      {
    push (@summary, "      $test\n");
      }
    }
  }

  push (@summary, "\n");
  if ($nottested)
  {
    push (@summary, "  Thorns with no valid testsuite parameter files:\n");
    push (@summary, "$nottested\n\n");
  }

  $unknown = 0;
  foreach $thorn (sort split(" ",$testdata->{"RUNNABLETHORNS"}))
  {
    if ($testdata->{"$thorn RUNNABLE"} !~ m:^\s*$:)
    {
      foreach $test (sort split(" ",$testdata->{"$thorn RUNNABLE"}))
      {
        $gotthorn = 0;
        if ($testdata->{sort "$thorn $test UNKNOWNFILES"})
        {
          if (!$unknown)
          {
            push (@summary, "  Thorns with unrecognized test output files:\n");
            $unknown = 1;
          }

          if (!$gotthorn)
          {
            push (@summary, "    $thorn\n");
            $gotthorn = 1;
          }
          push (@summary, "       $test: $testdata->{\"$thorn $test UNKNOWNFILES\"}\n");
        }
      }
    }
  }

  push (@summary, $separator2);

  push (@summary, "  Run details for configuration $config_data->{'CONFIG'}");
  push (@summary, '');

  foreach $thorn (sort split(" ",$testdata->{"RUNNABLETHORNS"}))
  {
    if ($testdata->{"$thorn RUNNABLE"} !~ m:^\s*$:)
    {
      foreach $test (sort split(' ',$testdata->{"$thorn RUNNABLE"}))
      {
        push (@summary, "      $thorn: $test");
        push (@summary, "         $rundata->{\"$thorn $test SUMMARY\"}");
      }
    }
    else
    {
      push (@summary, '      No tests available');
    }
  }
  push (@summary, '');

  push (@summary, $separator1);

  push (@summary, "  Summary for configuration $config_data->{'CONFIG'}");
  push (@summary, '');

  $total = $testdata->{"NUNRUNNABLE"}+$testdata->{"NRUNNABLE"};

  my $date     = `date`;     chomp($date);
  my $hostname = `hostname`; chomp($hostname);

  push (@summary, "    Time                     -> $date");
  push (@summary, "    Host                     -> $hostname");
  push (@summary, "    Processes                -> $config_data->{'NPROCS'}");
  push (@summary, "    User                     -> ".getpwuid($<));
  push (@summary, '');

  push (@summary, "    Total available tests    -> $total");
  push (@summary, "    Unrunnable tests         -> $testdata->{'NUNRUNNABLE'}");
  push (@summary, "    Runnable tests           -> $testdata->{'NRUNNABLE'}");
  push (@summary, '    Total number of thorns   -> '.scalar(split(' ',$testdata->{'THORNS'})));
  push (@summary, "    Number of tested thorns  -> $tested");

  push (@summary, "    Number of tests passed   -> $rundata->{'NPASSED'}");
  push (@summary, '    Number passed only to');
  push (@summary, "               set tolerance -> $rundata->{'NPASSEDTOTOL'}");
  push (@summary, "    Number failed            -> $rundata->{'NFAILED'}");

  if ($rundata->{'NPASSED'})
  {
    push (@summary, '');
    push (@summary, '  Tests passed:');
    push (@summary, '');
    foreach $thorn (sort split(' ',$testdata->{'THORNS'}))
    {
      foreach $file (sort split(' ',$rundata->{"$thorn PASSED"}))
      {
        push (@summary, "    $file (from $thorn)");
      }
    }
  }

  if ($rundata->{'NFAILED'})
  {
    push (@summary, '');
    push (@summary, '  Tests failed:');
    push (@summary, '');
    foreach $thorn (sort split(' ',$testdata->{'THORNS'}))
    {
      foreach $file (sort split(' ',$rundata->{"$thorn FAILED"}))
      {
        push (@summary, "    $file (from $thorn)");
      }
    }
  }

  push (@summary, '');
  push (@summary, $separator1);

  # write summary to both stdout and a summary logfile
  my $logdir = $config_data->{'TESTS_DIR'};
  mkdir ($logdir, 0755) if (not -e $logdir);
  $logdir .= "/$config_data->{'CONFIG'}";
  mkdir ($logdir, 0755) if (not -e $logdir);
  open (LOG, "> $logdir/summary.log")
    or die "Cannot open logfile '$logdir/summary.log'";
  print join ("\n", @summary);
  print LOG join ("\n", @summary);
  close (LOG);
}

############################################################
=pod 

=over 

=item ChooseTests($choice, $testdata)
 Allows the user to choose which tests to execute
 per thorn.

=back

=cut

############################################################
sub ChooseTests
{
  my ($choice,$testdata) = @_;
  my ($count,$arrangement,@myarrs,$arrchoice,$thorn,$mythorns,$mytests);
  my ($testcount,$test,$thornchoice);

  if ($choice =~ m:^A:i)
  {
    print "  No runnable testsuites in arrangements: ";
    print "$testdata->{\"UNRUNNABLEARRANGEMENTS\"}\n\n";

    print "  Arrangements with runnable testsuites:\n";
    $count = 1;
    foreach $arrangement (split(' ',$testdata->{"RUNNABLEARRANGEMENTS"}))
    {
      printf ("   [%2d] $arrangement\n",$count);
      $myarrs[$count] = "$arrangement";
      $count++;
    }
    while (!$arrchoice or $arrchoice eq " ")
    {
      $arrchoice = &defprompt("  Choose arrangement by number:"," ");
    }

    print "  No runnable testsuites in thorns: ";
    foreach  $thorn (split(" ",$testdata->{"UNRUNNABLETHORNS"}))
    {
      if ($testdata->{"$thorn ARRANGEMENT"} =~ m:^$myarrs[$arrchoice]:)
      {
        print "$thorn ";
      }
    }
    print "\n\n";

    print "  Thorns in $myarrs[$arrchoice] with runnable testsuites:\n";
    $count = 1;
    foreach $thorn (split(" ",$testdata->{"RUNNABLETHORNS"}))
    {
      if ($testdata->{"$thorn ARRANGEMENT"} =~ m:^$myarrs[$arrchoice]:)
      {
        printf ("  [%2d] $thorn\n",$count);
        $mythorns[$count] = "$thorn";
        $count++;
      }
    }
    while (!$thornchoice or $thornchoice eq " ")
    {
      $thornchoice = &defprompt("  Choose thorn by number:"," ");
    }
    $testcount = 0;
    printf ("  [ 0] All tests\n");
    foreach $test (split(" ",$testdata->{"$mythorns[$thornchoice] RUNNABLE"}))
    {
      $testcount++;
      printf ("  [%2d] $test\n",$testcount);
      print "       $testdata->{\"$mythorns[$thornchoice] $test DESC\"}\n";
      $mytests[$testcount] = "$test";
    }
    $testchoice = &defprompt("  Choose test:","0");

    if ($testchoice == 0)
    {
      $ntests = $testcount;
      for ($i=0;$i<$testcount;$i++)
      {
        $returntests[2*$i]   = $mytests[$i+1];
        $returntests[2*$i+1] = $mythorns[$thornchoice];
      }
    }
    else
    {
      $ntests = 1;
      $returntests[0] = $mytests[$testchoice];
      $returntests[1] = $mythorns[$thornchoice];
    }
  }
  elsif ($choice =~ m:^T:i)
  {
    $count = 1;
    foreach $thorn (sort split(' ',$testdata->{"RUNNABLETHORNS"}))
    {
      printf ("  [%2d] $thorn\n",$count);
      $mythorns[$count] = "$thorn";
      $count++;
    }
    if ($count > 1)
    {
      while (!$thornchoice or $thornchoice eq " ")
      {
        $thornchoice = &defprompt("  Choose thorn:"," ");
      }
      $testcount = 0;
      printf ("  [ 0] All tests\n");
      foreach $test (split(" ",$testdata->{"$mythorns[$thornchoice] RUNNABLE"}))
      {
        $testcount++;
        printf ("  [%2d] $test\n",$testcount);
        print "       $testdata->{\"$mythorns[$thornchoice] $test DESC\"}\n";
        $mytests[$testcount] = "$test";
      }
      $testchoice = &defprompt("  Choose test:","0");
      if ($testchoice == 0)
      {
        $ntests = $testcount;
        for ($i=0;$i<$testcount;$i++)
        {
          $returntests[2*$i]   = $mytests[$i+1];
          $returntests[2*$i+1] = $mythorns[$thornchoice];
        }
      }
      else
      {
        $ntests = 1;
        $returntests[0] = $mytests[$testchoice];
        $returntests[1] = $mythorns[$thornchoice];
      }
    }
  }

  return ($ntests,@returntests);
}


############################################################
=pod 

=over 

=item RunTest($output, $test, $thorn, $config_data, $testdata)
 This subroutine runs $test of $thorn with $config_data
 and returns $testdata.

=back

=cut

############################################################
sub RunTest
{
  my ($output,$test,$thorn,$config_data,$testdata) = @_;
  my ($test_dir,$config);
  my ($retcode);

  $testdata = &FindTestArchiveFiles($test,$thorn,$testdata);

  my $arrangement = $testdata->{"$thorn ARRANGEMENT"};

  $testdata->{"$thorn $test TESTRUNDIR"} = $config_data->{"TESTS_DIR"}.$sep.$config_data->{"CONFIG"}.$sep.$thorn;

  $testdata->{"$thorn $test TESTOUTPUTDIR"} = $testdata->{"$thorn $test TESTRUNDIR"}.$sep.$test;

  # Make any necessary directories
  &MakeTestRunDir($testdata->{"$thorn $test TESTRUNDIR"});

  my $parfile = TransformDirs($testdata->{"$thorn TESTSDIR"}. "/" . $test . ".par");

  # Clean the output directory for this test
  &CleanDir($testdata->{"$thorn $test TESTOUTPUTDIR"});

  # Run the test from the test thorn directory
  chdir ($testdata->{"$thorn $test TESTRUNDIR"}) ;

  # substitute the ($nprocs, $exe, $parfile) templates in the command
  my $cmd = $config_data->{'COMMAND'};
  $cmd =~ s/\$exe/$config_data->{'EXE'}/g;
  $cmd =~ s/\$nprocs/$config_data->{'NPROCS'}/g;
  $cmd =~ s/\$parfile/$parfile/g;

  $retcode = &RunCactus($output,$test,$cmd);
  chdir $config_data->{"CCTK_DIR"};

  # Deal with the error code
  if($retcode != 0)
  {
    print "Cactus exited with error code $retcode\n";
    print "Please check the logfile $testdata->{\"$thorn $test TESTRUNDIR\"}$sep$test.log\n\n";
    $testdata->{"$thorn FAILED"} .= "$parfile ";
    $testdata->{"NFAILED"}++;
  }

  return $testdata;
}

############################################################
=pod 

=over 

=item CompareTestFiles($test, $thorn, $runconfig, $rundata, $config_data, $testdata)
 Compares output from a particular testsuite.

=back

=cut

############################################################
sub CompareTestFiles
{
  my ($test,$thorn,$runconfig,$rundata,$config_data,$testdata) = @_;
  my ($test_dir,$file,$newfile,$oldfile);
  my ($vmaxdiff,$tmaxdiff,$numlines);

  $test_dir = $testdata->{"$thorn $test TESTOUTPUTDIR"};

  # Add new output files to database
  ($rundata->{"$thorn $test UNKNOWNFILES"},$rundata->{"$thorn $test TESTFILES"}) = &FindFiles("$test_dir",$testdata);
  $rundata->{"$thorn $test NUNKNOWNFILES"} = scalar(split(" ",$rundata->{"$thorn $test UNKNOWNFILES"}));
  $rundata->{"$thorn $test NTESTFILES"} = scalar(split(" ",$rundata->{"$thorn $test TESTFILES"}));

  $rundata->{"$thorn $test NFAILWEAK"}=0;
  $rundata->{"$thorn $test NFAILSTRONG"}=0;

  $abstol = $runconfig->{"ABSTOL"};
  $reltol = $runconfig->{"RELTOL"};

  if ($rundata->{"$thorn $test NTESTFILES"})
  {

    # Compare each file in the archived test directory
    foreach $file (split(" ",$testdata->{"$thorn $test DATAFILES"}))
    {
      my (@maxabsdiff, @absdiff, @valmax) = ();
      my ($filereltol, $fileabstol);

      $newfile = "$test_dir$sep$file";
      # This is the standard location of test data files
      $oldfile = "$testdata->{\"$thorn TESTSDIR\"}${sep}${test}${sep}$file";
      # If there is no data in the standard location, see if it got
      # decompressed to some place else
      if (!-d "$testdata->{\"$thorn TESTSDIR\"}${sep}${test}")
      {
        if (-d "$testdata->{\"$thorn $test UNCOMPRESSED\"}")
        {
          $oldfile = "$testdata->{\"$thorn $test UNCOMPRESSED\"}${sep}$file";
        }
      }

      $rundata->{"$thorn $test $file NINF"}=0;
      $rundata->{"$thorn $test $file NNAN"}=0;
      $rundata->{"$thorn $test $file NINFNOTFOUND"}=0;
      $rundata->{"$thorn $test $file NNANNOTFOUND"}=0;
      $rundata->{"$thorn $test $file NFAILSTRONG"}=0;
      $rundata->{"$thorn $test $file NFAILWEAK"}=0;

      ($fileabstol, $filereltol)=&GetFileTolerances($test,$thorn,$runconfig,$file);

      if ( -s $newfile && -s $oldfile)
      {
        open (INORIG, "<$oldfile") || print "Warning: Archive file $oldfile not found";
        open (INNEW,  "<$newfile") || print "Warning: Test file $newfile not found";

        $numlines = 0;

        while ($oline = <INORIG>)
        {
          # ignore comment lines in old file
          next if ($oline =~ /^\s*(["#].*)?$/);

          $nline = "";
          while ($nline = <INNEW>)
          {
            # ignore comment lines in new file
            last unless ($nline =~ /^\s*(["#].*)?$/);
          }

          # Now lets see if they differ.
          $numlines++;
          $nline = "\L$nline";
          $oline = "\L$oline";
          next if ($nline eq $oline);

          # Yes, they do. Check differences
          if (($nline !~ /(nan|inf)/) && ($oline !~ /(nan|inf)/))
          {
            # This is the new comparison (subtract last two numbers)

            # Make sure that floating point numbers have 'e' if exponential.
            $nline =~ tr/d/e/;
            $oline =~ tr/d/e/;

            my @newvals = split(' ',$nline);
            my @oldvals = split(' ',$oline);

            my $nnew = scalar(@newvals);
            $nold = scalar(@oldvals);

            my $allzero = 1;
            for ($count = 0; $count < $nold; $count++)
            {
              $absdiff[$count] = abs($newvals[$count] - $oldvals[$count]);
              $allzero = 0 if ($absdiff[$count]);
            }
            next if ($allzero);

            # They diff. But do they differ strongly?
            $rundata->{"$thorn $test $file NFAILWEAK"}++;

            # store difference for strong failures
            for ($count = 0; $count < $nold; $count++)
            {
              $maxabsdiff[$count] = $absdiff[$count]
                if ($maxabsdiff[$count] < $absdiff[$count]);
              my $absoldval = abs ($oldvals[$count]);
              my $absnewval = abs ($newvals[$count]);
              $valmax[$count] = $absoldval > $absnewval ?
                                $absoldval : $absnewval;
            }

            for ($count = 0; $count < $nold; $count++)
            {
              my $vreltol = $filereltol * $valmax[$count];
              my $vtol = $fileabstol > $vreltol ? $fileabstol : $vreltol;
              if ($absdiff[$count] >= $vtol) {
                $rundata->{"$thorn $test $file NFAILSTRONG"}++;
                last;
              }
            }
          }
          # Check against nans
          elsif ($nline =~ /nan/ && $oline !~ /nan/)
          {
            $rundata->{"$thorn $test $file NNAN"}++;
            $rundata->{"$thorn $test $file NFAILWEAK"}++;
            $rundata->{"$thorn $test $file NFAILSTRONG"}++;
          }
          # Check against inf
          elsif ($nline =~ /inf/ && $oline !~ /inf/)
          {
            $rundata->{"$thorn $test $file NINF"}++;
            $rundata->{"$thorn $test $file NFAILWEAK"}++;
            $rundata->{"$thorn $test $file NFAILSTRONG"}++;
          }
          elsif ($oline =~ /nan/)
          {
            $rundata->{"$thorn $test $file NNANNOTFOUND"}++;
            $rundata->{"$thorn $test $file NFAILWEAK"}++;
            $rundata->{"$thorn $test $file NFAILSTRONG"}++;
          }
          elsif ($oline =~ /inf/)
          {
            $rundata->{"$thorn $test $file NINFNOTFOUND"}++;
            $rundata->{"$thorn $test $file NFAILWEAK"}++;
            $rundata->{"$thorn $test $file NFAILSTRONG"}++;
          }
          else
          {
            print "TESTSUITE ERROR: Didn't catch case in CompareFiles\n";
          }
        } #while
        while ($nline = <INNEW>)
        {
          # ignore trailing comment lines in new file
          last unless ($nline =~ /^\s*(["#].*)?$/);
        }
        if (!eof(INNEW))
        {
          $rundata->{"$thorn $test $file NFAILWEAK"}++;
          $rundata->{"$thorn $test $file NFAILSTRONG"}++;
        }

      }
      elsif (!-e $newfile && -s $oldfile)
      {
        print "     $file in archive but not created in test\n";
        $rundata->{"$thorn $test NFAILWEAK"}++;
        $rundata->{"$thorn $test NFAILSTRONG"}++;
        $rundata->{"$thorn $test $file NFAILSTRONG"}++;
      }
      elsif (!-e $newfile && -z $oldfile)
      {
        print "     $file in archive but not created in test\n";
        print "       ($file empty in archive)\n";
        $rundata->{"$thorn $test NFAILWEAK"}++;
        $rundata->{"$thorn $test NFAILSTRONG"}++;
      }
      elsif (-e $newfile && -s $oldfile && -z $newfile)
      {
        print "     $file is empty in test\n";
        $rundata->{"$thorn $test NFAILWEAK"}++;
        $rundata->{"$thorn $test NFAILSTRONG"}++;
        $rundata->{"$thorn $test $file NFAILSTRONG"}++;
      }
      elsif (-e $newfile && -z $oldfile && -z $newfile)
      {
        print "     $file empty in both test and archive\n";
      }
      elsif (-e $newfile && -z $oldfile && -s $newfile)
      {
        print "     $file is empty in archive but not in test\n";
        $rundata->{"$thorn $test NFAILWEAK"}++;
        $rundata->{"$thorn $test NFAILSTRONG"}++;
        $rundata->{"$thorn $test $file NFAILSTRONG"}++;
      }
      else
      {
        print "     TESTSUITE ERROR: $newfile not compared to $oldfile\n";
      }

      my $havediffs = 0;
      my @maxreldiff = @maxabsdiff;
      for ($count = 0; $count < $nold; $count++)
      {
        next unless ($maxreldiff[$count]);
        $havediffs = 1;
        if ($valmax[$count] > 0)
        {
          $maxreldiff[$count] /= $valmax[$count];
        }
        else
        {
          print "ERROR: How did I get here, maximum difference is $maxabsdiff[$count] and maximum value is $valmax[$count] for $file\n";
        }
      }
      if ($havediffs)
      {
        splice (@maxabsdiff, $nold);
        splice (@maxreldiff, $nold);
        $rundata->{"$thorn $test $file MAXABSDIFF"} = \@maxabsdiff;
        $rundata->{"$thorn $test $file MAXRELDIFF"} = \@maxreldiff;
      }

      $rundata->{"$thorn $test $file NUMLINES"} = $numlines;

    }
  }
  else
  {
    print "  \n  No files created in test directory\n";
    $rundata->{"$thorn $test NFAILWEAK"} = $testdata->{"$thorn $test NDATAFILES"};
    $rundata->{"$thorn $test NFAILSTRONG"} = $testdata->{"$thorn $test NDATAFILES"};
  }

  return $rundata;
}

############################################################
=pod 

=over 

=item ReportOnTest($test, $thorn, $rundata, $testdata)
 Prints a report of a specific test, e.g. NFAILSTRONG,
 NNANNOTFOUND, NINF, NINFNOTFOUND, NFAILSTRONG. 

=back

=cut

############################################################
sub ReportOnTest
{
  my($test,$thorn,$rundata,$testdata) = @_;
  my($file,$tmp,$summary);
  my @log = ();

  # Different lines in files
  push (@log, '') if (@log);
  foreach $file (sort(split(' ',$testdata->{"$thorn $test DATAFILES"})))
  {
    my $key = "$thorn $test $file";
    next unless ($rundata->{"$key NFAILWEAK"} > 0);

    $rundata->{"$thorn $test NFAILWEAK"}++;
    # push (@log, '');
    if ($rundata->{"$key NFAILSTRONG"} == 0)
    {
      push (@log, "   $file: differences below tolerance on $rundata->{\"$key NFAILWEAK\"} lines");
    }
    else
    {
      $rundata->{"$thorn $test NFAILSTRONG"}++;
      push (@log, "   $file: substantial differences");
      my $tmp = $rundata->{"$key NNAN"};
      push (@log, "      caught  $tmp NaNs in new $file") if $tmp;
      $tmp = $rundata->{"$key NNANNOTFOUND"};
      push (@log, "      did not reproduce  $tmp NaNs from old $file") if $tmp;
      $tmp = $rundata->{"$key NINF"};
      push (@log, "      caught  $tmp Infs in new $file") if $tmp;
      $tmp = $rundata->{"$key NINFNOTFOUND"};
      push (@log, "      did not reproduce  $tmp Infs from old $file") if $tmp;
      $tmp = $rundata->{"$key NFAILSTRONG"};
      push (@log, "      significant differences on $tmp (out of $rundata->{\"$key NUMLINES\"}) lines");

      my $maxabsdiff = $rundata->{"$key MAXABSDIFF"};
      my $maxreldiff = $rundata->{"$key MAXRELDIFF"};
      my @tmp = ();
      for (my $c = 0; $c < scalar (@$maxabsdiff); $c++)
      {
        next unless ($$maxabsdiff[$c] > 0);
        my $col = $c + 1;
        push (@log, "      maximum absolute difference in column $col is $$maxabsdiff[$c]");
        push (@tmp, "      maximum relative difference in column $col is $$maxreldiff[$c]");
      }
      push (@log, @tmp);

      $tmp = $rundata->{"$key NFAILWEAK"} - $rundata->{"$key NFAILSTRONG"};
      push (@log, "      (insignificant differences on $tmp lines)") if $tmp;
    }
  }

  # Give a warning if there were different files created

  # Look for files created by test not in archive
  # (Note this is not so bad)
  push (@log, '') if (@log);
  foreach $file (split (" ",$rundata->{"$thorn $test TESTFILES"}))
  {
    $myfile = quotemeta($file);
    if ($testdata->{"$thorn $test DATAFILES"} !~ m:\b$myfile\b:)
    {
      push (@log, "   $file: not in thorn archive");
      $rundata->{"$thorn $test NFILEEXTRA"}++;
      $rundata->{"$thorn $test FILEEXTRA"} .= " $file";
    }
  }

  # Look for files in archive which are not created in test
  # (Note this is bad)
  push (@log, '') if (@log);
  my %filesmissing = ();
  if ($rundata->{"$thorn $test NTESTFILES"})
  {
    foreach $file (split (" ",$testdata->{"$thorn $test DATAFILES"}))
    {
      $myfile = quotemeta($file);
      if ($rundata->{"$thorn $test TESTFILES"} !~ m:\b$myfile\b:)
      {
        push (@log, "   $file: not created in test");
        $filesmissing{$file} = 1;
      }
    }
  }
  else
  {
    foreach $file (split (" ",$testdata->{"$thorn $test DATAFILES"}))
    {
      $filesmissing{$file} = 1;
    }
  }
  $rundata->{"$thorn $test FILEMISSING"} = \%filesmissing;

  # Ensure final newline character
  push (@log, '') if (@log);

  # write diffs to STDOUT and logfile
  if (@log)
  {
    print "\n", join ("\n", @log);

    # print LOG join ("\n", @log);
    my $logfile = $testdata->{"$thorn $test TESTRUNDIR"} . "/$test.diffs";
    open (LOG, "> $logfile") or die "Couldn't open logfile '$logfile'";
    print LOG join ("\n", @log);
    close (LOG);
  }

  if (! $rundata->{"$thorn $test NFAILWEAK"})
  {
    $summary = "Success: $testdata->{\"$thorn $test NDATAFILES\"} files identical";
    printf("\n  $summary\n");
    $rundata->{"$thorn PASSED"} .= "$test ";
    $rundata->{"NPASSED"}++;
  }
  else
  {
    if (! $rundata->{"$thorn $test NFAILSTRONG"})
    {
      $summary = "Success: $testdata->{\"$thorn $test NDATAFILES\"} files compared, $rundata->{\"$thorn $test NFAILWEAK\"} differ in the last digits";
      printf "\n  $summary\n";
      $rundata->{"$thorn PASSED"} .= "$test ";
      $rundata->{"NPASSED"}++;
      $rundata->{"NPASSEDTOTOL"}++;
    }
    else
    {
      my $nfilesmissing = scalar (keys (%filesmissing));
      my $nfilescompared = $testdata->{"$thorn $test NDATAFILES"} - $nfilesmissing;
      my $nfilesfailweak = $rundata->{"$thorn $test NFAILWEAK"} - $nfilesmissing;
      my $nfilesfailstrong = $rundata->{"$thorn $test NFAILSTRONG"} - $nfilesmissing;
#      $summary = "Failure: $testdata->{\"$thorn $test NDATAFILES\"} files compared, $rundata->{\"$thorn $test NFAILWEAK\"} differ, $rundata->{\"$thorn $test NFAILSTRONG\"} differ significantly";
      $summary = 'Failure: ';
      $summary .= "$nfilesmissing files missing, " if $nfilesmissing;
      $summary .= "$nfilescompared files compared, $nfilesfailweak differ";
      $summary .= ", $nfilesfailstrong differ significantly" if $nfilesfailstrong;
      printf "\n  $summary\n";
      $rundata->{"$thorn FAILED"} .= "$test ";
      $rundata->{"NFAILED"}++;
    }
  }
  $rundata->{"$thorn $test SUMMARY"} = $summary;
  printf ("\n");

  return $rundata;
}


############################################################
=pod 

=over 

=item ResetTestStatistics($rundata, $testdata)
 Resets $rundata values NFAILED, NPASSED to null values.

=back

=cut

############################################################
sub ResetTestStatistics
{
  my($rundata,$testdata) = @_;

  $rundata->{"NFAILED"} = 0;
  $rundata->{"NPASSED"} = 0;
  $rundata->{"NPASSEDTOTOL"} = 0;
  foreach $thorn (split(" ",$testdata->{"THORN"}))
  {
    $rundata->{"$thorn TESTED"} = 0;
  }

  return $rundata;
}



############################################################
=pod 

=over 

=item ParseAllParameterFiles($testdata, $config_data, $rundata)
 Given tests to do in $testdata, determines 

=back

=cut

############################################################
sub ParseAllParameterFiles
{
  my($testdata, $config_data, $rundata) = @_;
  my $nprocs_available = $config_data->{'NPROCS'};

  # Collect thorns needed for each testsuite
  foreach $thorn (split(" ",$testdata->{"THORNS"}))
  {
    $arr = $testdata->{"$thorn ARRANGEMENT"};
    my $nprocs = $rundata->{"$thorn NPROCS"};
    $nprocs = $nprocs_available unless ($nprocs);

    my $nrunnable = 0;
    foreach $testbase (split(" ",$testdata->{"$thorn TESTS"}))
    {
      my $nprocs_required = $rundata->{"$thorn $testbase NPROCS"};
      $nprocs_required = $nprocs unless ($nprocs_required);

      $parfile = "$testbase.par";

      # Set ActiveThorns and Description for this Test
      ($active,$desc) = &ParseParFile($thorn,$arr,$parfile,$config_data);
      $testdata->{"$thorn $testbase ACTIVE"} = $active;
      $testdata->{"$thorn $testbase DESC"} = $desc;

      # Find any missing thorns for this test
      ($nmissing,$missing) =
        &MissingThorns($testdata->{"$thorn $testbase ACTIVE"},
                       $testdata->{"THORNS"});

      # Set whether test is runnable or not
      if($nmissing)
      {
        $testdata->{"$thorn UNRUNNABLE"} .= "$testbase ";
        $testdata->{"$thorn $testbase MISSING"} .= $missing;
        $testdata->{'NUNRUNNABLE'}++;
      }
      elsif ($nprocs_required != $nprocs_available)
      {
        $testdata->{"$thorn UNRUNNABLE"} .= "$testbase ";
        $testdata->{"$thorn $testbase NPROCS"} = $nprocs_required;
        $testdata->{'NUNRUNNABLE'}++;
      }
      else
      {
        $testdata->{"$thorn RUNNABLE"} .= "$testbase ";
        $testdata->{"$thorn TESTED"} = 1;
        $testdata->{'NRUNNABLE'}++;
        $nrunnable++;
      }
    }

    if ($nrunnable)
    {
      $testdata->{'RUNNABLETHORNS'} .= "$thorn ";
      if ($testdata->{'RUNNABLEARRANGEMENTS'} !~ m:\b$arr\s:)
      {
        $testdata->{'RUNNABLEARRANGEMENTS'} .= "$arr ";
      }
    }
    else
    {
      $testdata->{'UNRUNNABLETHORNS'} .= "$thorn ";
    }
  }

  # Last look for arrangements with no runnable tests

  foreach $arr (split(" ",$testdata->{"ARRANGEMENTS"}))
  {
    if ($testdata->{'RUNNABLEARRANGEMENTS'} !~ m:\b$arr\s:)
    {
      $testdata->{'UNRUNNABLEARRANGEMENTS'} .= "$arr ";
    }
  }

  return $testdata;
}



############################################################
=pod 

=over 

=item MakeTestRunDir($dir)
 Given the name $dir, creates test run directory.

=back

=cut

############################################################
sub MakeTestRunDir
{
  my($dir) = @_;

  $dir =~  m:^(.*)/([^/]*)/([^/]*)$:;
  mkdir ($1,0755);
  mkdir ("$1/$2",0755);
  mkdir ("$1/$2/$3",0755);

}

# Views results from one particular test
############################################################
=pod 

=over 

=item ViewResults($test, $thorn, $runconfig, $rundata, $testdata)
 View the results of a test in a list or graphical format.

=back

=cut

############################################################
sub ViewResults
{
  my($test,$thorn,$runconfig,$rundata,$testdata) = @_;
  my($count,$choice,$myfile,@myfiles);

  if ($rundata->{"$thorn $test NTESTFILES"} && $rundata->{"$thorn $test NFAILSTRONG"})
  {
    $myfile = 1; # fixes uninitialized variable warning below
    &debug_print("thorn is '$thorn'");
    &debug_print("test is '$test'");
    &debug and $datafiles = $testdata->{"$thorn $test DATAFILES"};
    &debug_print("DATAFILES are '$datafiles'");
    &debug and my $nfailstrong = $rundata->{"$thorn $test NFAILSTRONG"};
    &debug_print("NFAILSTRING is '$nfailstrong'");
    while ($myfile !~ /^c/i)
    {
      $choice = 1;
      $count = 1;
      &debug_indent;
      my $filesmissing = $rundata->{"$thorn $test FILEMISSING"};
      foreach $file (sort split(" ",$testdata->{"$thorn $test DATAFILES"}))
      {
        &debug_print("considering file '$file'");
        if ($rundata->{"$thorn $test $file NFAILSTRONG"}
            and not $filesmissing->{$file})
        {
          $count==1 and print "  Files which differ strongly:\n";
          print "    [$count] $file\n";
          $myfiles[$count] = $file;
          $count++;
        }
      }
      &debug_dedent;

      if ($count>1) {
        $myfile = &defprompt("  Choose file by number or [c]ontinue","c");
      } else {
        $myfile = 'c';
      }

      while ($myfile !~ /^[c]/i && $choice !~ /^[c]/i)
      {
        print "  File $myfiles[$myfile] of test $test for thorn $thorn\n";
        my $oldfile = "$testdata->{\"$thorn TESTSDIR\"}/$test/$myfiles[$myfile]";
        my $newfile = "$testdata->{\"$thorn $test TESTOUTPUTDIR\"}/$myfiles[$myfile]";
        if (-e $oldfile and -e $newfile)
        {
          $choice = &defprompt("  Choose action [l]ist, [d]iff, [x]graph, [y]graph, [g]nuplot, [c]ontinue","c");
        } else
        {
          $choice = &defprompt("  Choose action [l]ist, [x]graph, [y]graph, [g]nuplot, [c]ontinue","c");
        }

        if ($choice =~ /^l/i)
        {
          if (-s $oldfile)
          {
            print "Archived file: $myfiles[$myfile]\n\n";
            open (ARCHIVE, "<$oldfile");
            while (<ARCHIVE>)
            {
              print;
            }
            close (ARCHIVE);
          }

          if (-s $newfile)
          {
            print "\n\nNew file: $myfiles[$myfile]\n\n";

            open (TEST, "<$newfile");
            while (<TEST>)
            {
              print;
            }
            close (TEST);
            print "\n";
          }
        }
        elsif ($choice =~ /^d/i)
        {
          &debug_print("oldfile is '$oldfile'");
          if (-e $oldfile and -e $newfile)
          {
            print "\n  Performing diff on  <archive> <test>\n\n";
            $command = "diff $oldfile $newfile\n";
            print "$command\n\n";
            system($command);
          } else {
            print "\n  Both versions of this file do not exist.  Please choose a different option.\n";
          }
          print "\n";
        }
        elsif ($choice =~ /^x/i)
        {
          print "  xgraph <archive> <test>\n\n"; # Should this be changed?
          $command = "xgraph ";
          -e $oldfile and $command .= "$oldfile ";
          -e $newfile and $command .= "$newfile ";
          $command .= "&\n";
          print "  $command\n";
          system($command);
        }
        elsif ($choice =~ /^y/i)
        {
          print "  ygraph <archive> <test>\n\n"; # Should this be changed?
          $command = "ygraph ";
          -e $oldfile and $command .= "$oldfile ";
          -e $newfile and $command .= "$newfile ";
          $command .= "&\n";
          print "  $command\n";
          system($command);
        }
        elsif ($choice =~ /^g/i)
        {
          print "  gnuplot <archive> <test>\n\n"; # Should this be changed?
          $command = ("gnuplot -persist <<EOF\n"
                      . "set grid\n"
                      . "plot ");
          -e $oldfile and $command .= "\"$oldfile\" w lp";
          if (-e $oldfile and -e $newfile)
          {
            $command .= ", ";
          }
          -e $newfile and $command .= "\"$newfile\" w lp ";
          $command .= "\nEOF";
          print "  $command\n";
          system($command);
        }
      }
    }
  }

  return;

}

############################################################
=pod 

=over 

=item TransformDirs($in)
 Transforms directory for use with Cygwin.

=back

=cut

############################################################
sub TransformDirs
{
  my ($in) = @_;
  my $out;

  $out = `cygpath -wa $in`;

  chomp $out;

  if ( ! $out )
  {
    $out = $in;
  }
  else
  {
    $out = '"'.$out.'"';
  }

  return $out;
}

1;
