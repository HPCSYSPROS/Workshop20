#! /usr/bin/perl -w

#/*@@
#  @file     ConfigurationParser.pl
#  @date     Tue Feb  8 17:36:48 2000
#  @author   Tom Goodale
#  @desc
#            Parser for configuration.ccl files
#  @enddesc
#  @version  $Header$
#@@*/


#/*@@
#  @routine    CreateConfigurationDatabase
#  @date       Tue Feb  8 17:47:26 2000
#  @author     Tom Goodale
#  @desc
#              Parses the information in all the thorns' configuration.ccl files
#              and creates a database from it
#  @enddesc
#@@*/
sub CreateConfigurationDatabase
{
  my($config_dir, %thorns) = @_;
  my(%cfg) = ();
  my(%thorn_dependencies);

  # Loop through each thorn's configuration file.
  foreach my $thorn (sort keys %thorns)
  {
    my $filename = "$thorns{$thorn}/configuration.ccl";
    next if (! -r $filename);

    # Get the configuration data from it
    &ParseConfigurationCCL($config_dir, $thorn, \%cfg, \%thorns, $filename);

#    my $debug = 1;
#    if($debug)
#    {
#      print "   $thorn\n";
#      print "           Provides:          ", $cfg{"\U$thorn\E PROVIDES"}, "\n"
#        if ($cfg{"\U$thorn\E PROVIDES"});
#      foreach my $provides (split (' ', $cfg{"\U$thorn\E PROVIDES"}))
#      {
#        print "           as version:        ", $cfg{"\U$thorn\E PROVIDES \U$provides\E VERSION"}, "\n"
#          if ($cfg{"\U$thorn\E PROVIDES \U$provides\E VERSION"});
#      }
#      print "           Requires:          ", $cfg{"\U$thorn\E REQUIRES"}, "\n"
#        if ($cfg{"\U$thorn\E REQUIRES"});
#      print "           Optional:          ", $cfg{"\U$thorn\E OPTIONAL"}, "\n"
#        if ($cfg{"\U$thorn\E OPTIONAL"});
#      print "           Optional-ifactive: ", $cfg{"\U$thorn\E OPTIONAL_IFACTIVE"}, "\n"
#        if ($cfg{"\U$thorn\E OPTIONAL_IFACTIVE"});
#    }

    $cfg{"\U$thorn\E USES THORNS"} = '';

    # verify that all required thorns are there in the ThornList
    next if (! $cfg{"\U$thorn\E REQUIRES THORNS"});

    my @missing = ();
    foreach my $required (split (' ', $cfg{"\U$thorn\E REQUIRES THORNS"}))
    {
      push (@missing, $required)
        if ((! $thorns{"$required"}) && (! $thorns{"\U$required\E"}));
    }
    if (@missing == 1)
    {
      &CST_error (0, "Thorn '$thorn' requires thorn '@missing'. " .
                     'Please add this thorn to your ThornList or remove ' .
                     "'$thorn' from it !");
    }
    elsif (@missing > 1)
    {
      &CST_error (0, "Thorn '$thorn' requires thorns '@missing'. " .
                     'Please add these thorns to your ThornList or ' .
                     "remove '$thorn' from it !");
    }

    $cfg{"\U$thorn\E USES THORNS"} .= $cfg{"\U$thorn\E REQUIRES THORNS"} . ' ';
  }

  # Turn optional capabilities into required capabilities, if the
  # capability is provided. This way we don't have to treat required
  # and optional requirements differently.
  my %providedcaps;
  foreach my $thorn (sort keys %thorns)
  {
      if ($cfg{"\U$thorn\E PROVIDES"})
      {
          foreach my $providedcap (split (' ', $cfg{"\U$thorn\E PROVIDES"}))
          {
              $providedcaps{$providedcap} = 1;
          }
      }
  }
  foreach my $thorn (sort keys %thorns)
  {
      if ($cfg{"\U$thorn\E REQUIRES"})
      {
          foreach my $requiredcap (split (' ', $cfg{"\U$thorn\E REQUIRES"}))
          {
              $cfg{"\U$thorn\E ACTIVATES"} .= "$requiredcap ";
          }
      }
  }
  foreach my $thorn (sort keys %thorns)
  {
      if ($cfg{"\U$thorn\E OPTIONAL"})
      {
          foreach my $optionalcap (split (' ', $cfg{"\U$thorn\E OPTIONAL"}))
          {
              if ($providedcaps{$optionalcap})
              {
                  $cfg{"\U$thorn\E REQUIRES"} .= "$optionalcap ";
                  $cfg{"\U$thorn\E ACTIVATES"} .= "$optionalcap ";
              }
          }
      }
      if ($cfg{"\U$thorn\E OPTIONAL_IFACTIVE"})
      {
          foreach my $optionalcap (split (' ', $cfg{"\U$thorn\E OPTIONAL_IFACTIVE"}))
          {
              if ($providedcaps{$optionalcap})
              {
                  $cfg{"\U$thorn\E REQUIRES"} .= "$optionalcap ";
                  # nothing is activated
              }
          }
      }
  }

  foreach my $thorn (sort keys %thorns)
  {
    # verify that all required capabilities are there in the ThornList
    next if (! $cfg{"\U$thorn\E REQUIRES"});

    foreach my $requiredcap (split (' ', $cfg{"\U$thorn\E REQUIRES"}))
    {
      my @found = ();
      foreach my $thorncap (sort keys %thorns)
      {
        foreach my $cap (split (' ', $cfg{"\U$thorncap\E PROVIDES"}))
        {
          push (@found, $thorncap)
            if ("\U$cap\E" eq "\U$requiredcap\E");
        }
      }

      # there must be exactly one thorn providing a required capability
      if (@found == 0)
      {
        &CST_error (0, "Thorn '$thorn' requires the capability " .
                       "'$requiredcap'.\n" .
                       "     Please add a thorn that provides '$requiredcap' " .
                       "to your ThornList or remove '$thorn' from it !")
      }
      elsif (@found > 1)
      {
        &CST_error (0, "More than one thorn provides the capability " .
                       "'$requiredcap'.\n" .
                       "     These thorns are: '@found'.\n" .
                       "     Please use only one !\n");
      }
      elsif ( $cfg{"\U$thorn\E REQUIRES \U$requiredcap\E VERSION"} )
      {
        if ( &CheckForCompatibleVersion(
                $cfg{"\U".$found[0]."\E PROVIDES \U$requiredcap\E VERSION"},
                $cfg{"\U$thorn\E REQUIRES \U$requiredcap\E VERSION"}) == 0 )
        {
          &CST_error (0, "Thorn '$thorn' requires the capability " .
                       "'$requiredcap' in version ".
                       $cfg{"\U$thorn\E REQUIRES \U$requiredcap\E VERSION"}.
                       ". Thorn ".$found[0]." provides $requiredcap, but ".
                       "in version ".
                       $cfg{"\U".$found[0]."\E PROVIDES \U$requiredcap\E VERSION"}.
                       ".\n");
        }
        $cfg{"\U$thorn\E USES THORNS"} .= $found[0] . ' ';
      }
      else
      {
        $cfg{"\U$thorn\E USES THORNS"} .= $found[0] . ' ';
      }
    }
  }

  # Translate capability to thorn names
  my %capabilities;
  foreach my $thorn (sort keys %thorns)
  {
      next if ! $cfg{"\U$thorn\E PROVIDES"};
      foreach my $cap (split (' ', $cfg{"\U$thorn\E PROVIDES"}))
      {
          $capabilities{"\U$cap\E"} = $thorn;
      }
  }
  foreach my $thorn (sort keys %thorns)
  {
      my $activates = '';
      foreach my $cap (split (' ', $cfg{"\U$thorn\E ACTIVATES"}))
      {
          my $cap_thorn = $capabilities{"\U$cap\E"};
          $activates .= " $cap_thorn";
      }
      $cfg{"\U$thorn\E ACTIVATES THORNS"} = $activates;
  }

  # Check for cyclic dependencies
  # create a hash with thorn-> used thorns (no prefix)
  foreach my $thorn (sort keys %thorns)
  {
    $thorn_dependencies{uc($thorn)}=$cfg{"\U$thorn\E USES THORNS"};
    $thorn_dependencies{uc($thorn)} =~ s/\b$thorn\b//i;
  }

  my $message = &find_dep_cycles(%thorn_dependencies);
  if ("" ne $message)
  {
    $message =~ s/^\s*//g;
    $message =~ s/\s*$//g;
    $message =~ s/\s+/->/g;
    $message = "Found a cyclic dependency in configuration requirements:$message\n";
    &CST_error(0, $message);
  }

  return \%cfg;
}


#/*@@
#  @routine    CompareVersionStrings
#  @date       Tue Oct 20 23:17:18 2015
#  @author     Frank Loeffler
#  @desc
#  Compares two version strings: first non-numeric prefix lexically, next
#  numeric prefix of remainder numerically, and so on.
#  @enddesc
#@@*/
sub CompareVersionStrings
{ 
  my($v1, $v2) = @_;
  my($nan1, $nan2, $num1, $num2, $ret);
  # the loop body strips recognized parts from the strings
  while($v1 ne "" or $v2 ne "") {
    # compare non-numeric prefix followed by numeric value if they exist
    # remove found sub-string from input
    $v1 =~ s/^([^0-9]*)([0-9]*)(.*)/$3/;
    $nan1 = $1;
    $num1 = $2;
    $v2 =~ s/^([^0-9]*)([0-9]*)(.*)/$3/;
    $nan2 = $1;
    $num2 = $2;
    $ret = ($nan1 cmp $nan2) || ($num1 <=> $num2);
    return $ret if ($ret != 0);
  }
  return 0;
}

#/*@@
#  @routine    CheckForCompatibleVersion
#  @date       Tue Oct 20 23:17:18 2015
#  @author     Frank Loeffler
#  @desc
#  Checks that two versions strings are compatible. The first argument is a raw
#  version string, the second argument has also an operator as prefix, which is
#  used to determine if these two match. Returns 1 for success and 0 for failure.
#  @enddesc
#@@*/
sub CheckForCompatibleVersion
{
  my($v1,$fv2) = @_;
  my($op, $v2, $cmp);
  $fv2 =~ m/(<<|<=|-|>=|>>)(.*)/;
  $op = $1;
  $v2 = $2;
  $cmp = &CompareVersionStrings($v1, $v2);
  return 1 if ($op eq '<<' and $cmp <  0);
  return 1 if ($op eq '<=' and $cmp <= 0);
  return 1 if ($op eq '='  and $cmp == 0);
  return 1 if ($op eq '>=' and $cmp >= 0);
  return 1 if ($op eq '>>' and $cmp >  0);
  return 0;
}

#/*@@
#  @routine    ParseConfigurationCCL
#  @date       Tue Feb  8 19:23:18 2000
#  @author     Tom Goodale
#  @desc
#  Parses a configuration.ccl file and generates a database of the values
#  @enddesc
#@@*/
sub ParseConfigurationCCL
{
  my($config_dir, $thorn, $cfg, $thorns, $filename) = @_;
  my(@data);
  my($line_number, $line);
  my($provides, $script, $lang, $options);
  my($optional, $define);

  # Initialise some stuff to prevent perl -w from complaining.

  $cfg->{"\U$thorn\E PROVIDES"} = '';
  $cfg->{"\U$thorn\E REQUIRES"} = '';
  $cfg->{"\U$thorn\E REQUIRES THORNS"} = '';
  $cfg->{"\U$thorn\E OPTIONAL"} = '';
  $cfg->{"\U$thorn\E OPTIONAL_IFACTIVE"} = '';
  $cfg->{"\U$thorn\E ACTIVATES"} = '';
  $cfg->{"\U$thorn\E OPTIONS"}  = '';

  # Read the data
  @data = &read_file($filename);

  for($line_number = 0; $line_number < @data; $line_number++)
  {
    $line = $data[$line_number];
    # Parse the line
    if($line =~ m/^\s*PROVIDES\s*/i)
    {
      $lang = $script = '';
      ($provides, $script, $lang, $options, $line_number, $version) = &ParseProvidesBlock($line_number, \@data);
      if ($provides !~ m{^[A-Za-z0-9_.]+$}) {
        &CST_error (0, "Illegal capability name '$provides' line '$line' in configure.ccl of thorn '$thorn'");
      }
      if ($lang !~ m{^[A-Za-z0-9_.]*$}) {
        &CST_error (0, "Illegal script language '$lang' line '$line' in configure.ccl of thorn '$thorn'");
      }
      $cfg->{"\U$thorn\E PROVIDES"} .= "$provides ";
      $cfg->{"\U$thorn\E PROVIDES \U$provides\E VERSION"} = "$version";
      if($script)
      {
        $cfg->{"\U$thorn\E PROVIDES \U$provides\E SCRIPT"} = "$thorns->{$thorn}/$script";
      }
      $cfg->{"\U$thorn\E PROVIDES \U$provides\E LANG"} = $lang;
      $cfg->{"\U$thorn\E PROVIDES \U$provides\E OPTIONS"} = $options;

#      if ($script)
#      {
#        print "Running configuration script '$script'\n";
#
#        &ParseConfigScript($config_dir, $provides, $lang, $script,
#                           $thorn, $cfg);
#        print "\n";
#      }

      next;
    }
    elsif($line =~ m/^\s*REQUIRES\s+THORNS\s*:\s*(.*)/i)
    {
      my $newlist = $1;
      $newlist =~ s/\b$thorn\b//i;
      $newlist =~ s/,/ /g;
      my $oldlist = $cfg->{"\U$thorn\E REQUIRES THORNS"};
      my $list = $oldlist . ' ' . $newlist;
      $list = join (' ', sort split (' ', $list));
      $cfg->{"\U$thorn\E REQUIRES THORNS"} = $list;
#      if ($cfg->{"\U$thorn\E REQUIRES THORNS"})
#      {
#        &CST_error (3, '\'Requires Thorns\' will not be supported in release beta-14' .
#        "\n Please adjust thorn \U$thorn\E to use \'Requires\' instead.");
#      }
    }
    elsif($line =~ m/^\s*REQUIRES\s+(.*)/i)
    {
      my $cap = $1;
      if ($cap !~ m{^([A-Za-z0-9_.]+ *(\( *(<<|<=|=|>=|>>) *[0-9a-zA-Z.+-:]+ *\))?)+$}) {
        &CST_error (0, "Illegal required capability '$cap' line '$line' in configure.ccl of thorn '$thorn'");
      }
      while ($cap =~ m/ *([A-Za-z0-9_.]+)( *\((.+)\))?/g)
      {
        my $capability = $1;
        my $version    = $3;
        $version =~ s/ //g;
        $cfg->{"\U$thorn\E REQUIRES"} .= "$capability ";
        if ($version)
        {
          $cfg->{"\U$thorn\E REQUIRES \U$capability\E VERSION"} .= "$version";
        }
      }
    }
    elsif($line =~ m/^\s*OPTIONAL\s+/i)
    {
      ($optional, $define, $line_number) = &ParseOptionalBlock($filename, $line_number, \@data);
      if ($optional !~ m{^[A-Za-z0-9_. ]+$}) {
        &CST_error (0, "Illegal optional capability '$optional' line '$line' in configure.ccl of thorn '$thorn'");
      }
      $cfg->{"\U$thorn\E OPTIONAL"} .= "$optional ";
      $cfg->{"\U$thorn\E OPTIONAL \U$optional\E DEFINE"} = $define;
    }
    elsif($line =~ m/^\s*OPTIONAL_IFACTIVE\+*/i)
    {
      ($optional, $define, $line_number) = &ParseOptionalBlock($filename, $line_number, \@data);
      if ($optional !~ m{^[A-Za-z0-9_. ]+$}) {
        &CST_error (0, "Illegal optional capability '$optional' line '$line' in configure.ccl of thorn '$thorn'");
      }
      $cfg->{"\U$thorn\E OPTIONAL_IFACTIVE"} .= "$optional ";
      $cfg->{"\U$thorn\E OPTIONAL_IFACTIVE \U$optional\E DEFINE"} = $define;
    }
    elsif($line =~ m/^\s*NO_SOURCE\s*$/i)
    {
      $cfg->{"\U$thorn\E OPTIONS"} .= "NO_SOURCE";
    }
    else
    {
      chomp($line);
      &CST_error (0, "Unrecognised line '$line' in configure.ccl of thorn '$thorn'");
    }
  }
}


#/*@@
#  @routine    ParseProvidesBlock
#  @date       Mon May  8 15:52:40 2000
#  @author     Tom Goodale
#  @desc
#  Parses the PROVIDES block in a configuration.ccl file.
#  @enddesc
#@@*/
sub ParseProvidesBlock
{
  my ($line_number, $data) = @_;
  my ($provides, $script, $lang, $options, $version);

  $provides = "";
  $script   = "";
  $lang     = "";
  $version  = "0.0.1";
  $options  = [];

  $data->[$line_number] =~ m/^\s*PROVIDES\s*(.*)/i;

  $provides = $1;

  $line_number++;
  if($data->[$line_number] !~ m/^\s*\{\s*$/)
  {
    &CST_error (0, "Error parsing provides block line '$data->[$line_number]' $file_name:$line_number ".
                   'Missing { at start of block');
    $line_number++ while(defined($data->[$line_number]) and $data->[$line_number] !~ m:\s*\}\s*:);
  }
  else
  {
    while(defined($data->[$line_number]) and $data->[$line_number] !~ m:\s*\}\s*:)
    {
      $line_number++;
      if($data->[$line_number] =~ m/^\s*SCRIPT\s*(.*)$/i)
      {
        $script = $1;
        next;
      }
      elsif($data->[$line_number] =~ m/^\s*LANG[^\s]*\s*(.*)$/i)
      {
        $lang = $1;
        next;
      }
      elsif($data->[$line_number] =~ m/^\s*OPTIONS[^\s]*\s*(.*)$/i)
      {
        push(@$options, split(' ',$1));
        next;
      }
      elsif($data->[$line_number] =~ m/^\s*VERSION\s+(.+)$/i)
      {
        $version = $1;
        if ($1 !~ m/[0-9]([0-9a-z.+-:]*)/i)
        {
          print STDERR "Error in version specification '"+$version+"'. "+
                       "Only alphanumeric characters and . + - : are allowed, "+
                       "and a version has to start with a digit."
          &CST_error (0, 'Unrecognised version');
        }
        next;
      }
      elsif($data->[$line_number] =~ m:\s*\}\s*:)
      {
        # do nothing.
      }
      else
      {
        print STDERR "Error parsing provides block line '$data->[$line_number]'\n";
        &CST_error (0, 'Unrecognised statement');
      }
    }
  }

  return ($provides, $script, $lang, $options, $line_number, $version);
}


#/*@@
#  @routine    ParseOptionalBlock
#  @date       Mon May  8 15:52:40 2000
#  @author     Tom Goodale
#  @desc
#  Parses the OPTIONAL or OPTIONAL_IFACTIVE block in a configuration.ccl file.
#  @enddesc
#@@*/
sub ParseOptionalBlock
{
  my ($file_name, $line_number, $data) = @_;
  my ($optional, $define);

  $data->[$line_number] =~ m/^\s*OPTIONAL(_IFACTIVE)?\s*(.*)/i;

  $optional = $2;

  $define = "";

  $line_number++;

  if($data->[$line_number] !~ m/^\s*\{\s*$/)
  {
    &CST_error (0, "Error parsing optional block line '$data->[$line_number]' $file_name:$line_number".
                ' Missing { at start of block.');
    $line_number++ while(defined($data->[$line_number]) and $data->[$line_number] !~ m:\s*\}\s*:);
  }
  else
  {
    while(defined($data->[$line_number]) and $data->[$line_number] !~ m:\s*\}\s*:)
    {
      $line_number++;
      if($data->[$line_number] =~ m/^\s*DEFINE\s*(.*)$/i)
      {
        if($define eq "")
        {
          $define = $1;
          next;
        }
        else
        {
          &CST_error (0, "Error parsing optional block line '$data->[$line_number]' " . 'Only one define allowed.');
        }
      }
      elsif($data->[$line_number] =~ m:\s*\}\s*:)
      {
        # do nothing.
      }
      else
      {
        &CST_error (0, "Error parsing provides block line '$data->[$line_number]' " . 'Unrecognised statement.');
      }
    }
  }

  return ($optional, $define, $line_number);
}

1;
