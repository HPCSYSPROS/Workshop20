#/*@@
#  @file      ProcessConfiguration.pl
#  @date      Mon May  8 15:52:08 2000
#  @author    Tom Goodale
#  @desc 
#  
#  @enddesc 
#@@*/

#/*@@
#  @routine    SplitThorns
#  @date       Mon May  8 16:04:59 2000
#  @author     Tom Goodale
#  @desc 
#  Splits the thorns hash into those with source and those without source
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#
#@@*/
sub SplitThorns
{
  my ($configuration_database, $thorns, $source_thorns, $nosource_thorns) = @_;

  foreach $thorn (sort keys %$thorns)
  {
    if($configuration_database->{"\U$thorn OPTIONS\E"} =~ m/NO_SOURCE/i)
    {
      $nosource_thorns->{"$thorn"} = $thorns->{"$thorn"};
    }
    else
    {
      $source_thorns->{"$thorn"} = $thorns->{"$thorn"};
    }
  }
}

#/*@@
#  @routine    ProcessConfiguration
#  @date       Thu Aug 26 22:09:26 2004
#  @author     Tom Goodale
#  @desc 
#  Runs all the configuration scripts belonging to thorns.
#  Has to setup the environment correctly, and as a result
#  modifies the config-info file for future reference.
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#
#@@*/
sub ProcessConfiguration
{
  my($config_dir,$config_database, $thorns, $config_file) = @_;

  my $thorn;
  my $provides;
  my @allowed_opts;

  # Find the master list of allowed options
  foreach $thorn (sort keys %thorns)
  {
#    print "DEBUG: Thorn $thorn\n";

    foreach $provides (split(' ',$config_database->{"\U$thorn\E PROVIDES"}))
    {
#      print "DEBUG: Provides $provides\n";

      if(@{$config_database->{"\U$thorn\E PROVIDES \U$provides\E OPTIONS"}} != 0)
      {
#        print "DEBUG: $thorn provides $provides with options\n";
        push(@allowed_opts, @{$config_database->{"\U$thorn\E PROVIDES \U$provides\E OPTIONS"}});
      }
    }
  }

#  print "DEBUG: allowed options are @allowed_opts\n";

  # Now get all the configuration options.
  my $env;
  my $optfile;
  my $configinfo;
  my $headers;

  ($configinfo,$headers) = ParseConfigInfo($config_file);

  if($ENV{"options"})
  {
    $optfile = ParseOptionsFile($ENV{"options"})
  }
  else
  {
    $optfile = {};
  }

  $env = GetOptionsFromEnv(\%ENV, \@allowed_opts);

  my $modified = AmalgamateOptions($env,$optfile,$configinfo,\@allowed_opts);

  # Write a new config-info file if anything has changed
  if($modified)
  {
    WriteNewConfigInfo($config_file,$headers,$configinfo);
  }

  # Now setup the environment
  foreach $option (@allowed_opts)
  {
    if (defined($configinfo->{$option}))
    {
      $ENV{$option} = $configinfo->{$option};
    }
  }

  # Ok, can now run the configuration scripts.

  my %thorns_todo = ();
  map { $thorns_todo{"\U$_\E"} = 1; } keys %thorns;
  my %requirements_done = ();

  my $made_progress = 1;
#  print "DEBUG: Processing thorn provisions:\n";
  while (keys %thorns_todo and $made_progress)
  {
#    print "DEBUG:    Thorns left to do: " . scalar(%thorns_todo) . "\n";
#    map { print "DEBUG:       - $_\n"; } sort keys %thorns_todo;
    $made_progress = 0;
    THORN: foreach $thorn (sort keys %thorns_todo)
    {
#      print "DEBUG:       Checking thorn $thorn\n";
      my @provides_list = split(' ',$config_database->{"\U$thorn\E PROVIDES"});
#      print "DEBUG:          Provides: @provides_list\n";
      my @requires_list = split(' ',$config_database->{"\U$thorn\E REQUIRES"});
#      print "DEBUG:          Requires: @requires_list\n";

      my %need = ();
      map { $need{"\U$_\E"} = 1; } @requires_list;
      map { delete $need{"\U$_\E"}; } @provides_list;
      map { next THORN unless exists $requirements_done{"\U$_\E"}; } keys %need;

#      print "DEBUG:       Processing thorn $thorn\n";
      foreach my $provides (@provides_list)
      {
#        print "DEBUG:          Processing provision $provides\n";
        my $script = $config_database->{"\U$thorn\E PROVIDES \U$provides\E SCRIPT"};
        my $lang   = $config_database->{"\U$thorn\E PROVIDES \U$provides\E LANG"};

        if ($script)
        {
#          print "DEBUG: Running configuration script '$script'\n";
  
          &ParseConfigScript($config_dir, $provides, $lang, $script,
                             $thorn, $config_database);
#          print "DEBUG: \n";
        }

        # Add make definitions to the environment, so that they are
        # available to the following scripts
        my $config = $config_database->{"\U$thorn $provides\E MAKE_DEFINITION"};
        my %options = $config =~ /^\s*(\w+)\s*=(.*)$/mg;
        foreach my $option (keys %options)
        {
          my $value = $options{$option};
          $value =~ s/^\s*//;
          $value =~ s/\s*$//;
#          print "DEBUG: Thorn $thorn, providing $provides, setting \$ENV{$option}=\"$value\"\n";
          $ENV{$option} = $value;
        }

        $requirements_done{"\U$provides\E"} = 1;
      }

      delete $thorns_todo{"\U$thorn\E"};
      $made_progress = 1;
    }
#    if (! $made_progress)
#    {
#      print "DEBUG:    Provided requirements: " . scalar($requirements_done) . "\n";
#      map { print "DEBUG:       - $_\n"; } sort keys %requirements_done;
#    }
#    die unless $made_progress;
  }
}

1;
