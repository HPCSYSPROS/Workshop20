#! /usr/bin/perl -sw
#/*@@
#  @file    setup_configuration.pl
#  @date    Fri Jan  8 13:48:48 1999
#  @author  Tom Goodale
#  @desc
#           Setup file for a new configuration in the CCTK
#           Invocation is
#           setup_configuration [-reconfig] [-config_file=<options>] config_name
#  @enddesc
#  @version $Header$
#@@*/

chop ($top = `pwd`);

$configure = "sh $top/lib/make/configure";

if($#ARGV > -1)
{
  $config = shift(@ARGV);
}
else
{
  chop ($config = `uname`);
}
# Replace slashes with underscores.
$config =~ s:[/\\]:_:g;

# Check for troublesome configuration names and abort in that case
my @invalid = ("build", "clean", "cleandeps", "cleanobjs", "config",
               "configinfo", "delete", "editthorns", "realclean",
               "rebuild", "reconfig", "utils", "testsuite",
               "thornlist", "Thornguide", "update", "examples");
for my $i (@invalid)
{
  if ($config =~ /-$i$/)
  {
    die "The suffix -$i in the configuration name $config is invalid.\n".
        "Please choose another name for the configuration.\n";
  }
}

# Work out if there is a user default file
if($ENV{CACTUSRC_DIR})
{
  if (-e "$ENV{CACTUSRC_DIR}/.cactus/config")
  {
    $default_file = "$ENV{CACTUSRC_DIR}/.cactus/config";
  }
}
elsif (-e "$ENV{HOME}/.cactus/config")
{
  $default_file = "$ENV{HOME}/.cactus/config";
}

# Work out where the config directory is
$configs_dir = $ENV{CONFIGS_DIR};
$configs_dir = 'configs' if (! $configs_dir);

# The configs directory doesn't exist.
if (! -d "$configs_dir" && ! -l "$configs_dir")
{
  print "Completely new cactus build.  Creating config database\n";

  mkdir("$configs_dir", 0755) ||
    die "Internal error - couldn't create '$configs_dir'\n";
}

chdir "$configs_dir" ||
  die "Internal error - couldn't enter '$configs_dir'\n";

# The specified configuration doesn't exist
if (! -d "$config" && ! -l "$config")
{
  print "Creating new configuration $config.\n";

  for $dir ("$config", "$config/build", "$config/lib", "$config/scratch", "$config/config-data")
  {
    mkdir("$dir",0755) || die "Internal error - couldn't create $dir";
  }
}
else
{
  print "Reconfiguring $config.\n";
}


%CONFIGURED = &SetConfigureEnv();

$configure_command = &DetermineConfigureCommand($configure);

chdir "$config/config-data" ||
  die "Internal error - couldn't enter '$configs_dir/$config/config-data'\n";

# remove cached configure options
unlink 'config.cache' if (-f 'config.cache');

system("$configure_command");
$retcode = $? >> 8;

chdir '..';

open(INFO, ">config-info") ||
  die "Internal error - couldn't create '$configs_dir/$config/config-info'\n";

print INFO "# CONFIGURATION  : $config\n";
print INFO "# CONFIG-DATE    : " . gmtime(time()) . " (GMT)\n";
chop ($hostname = `hostname`);
print INFO "# CONFIG-HOST    : $hostname\n";
print INFO "# CONFIG-STATUS  : $retcode\n";
print INFO "# CONFIG-OPTIONS :\n";
foreach $setting (sort keys %CONFIGURED)
{
  print INFO "$setting=$CONFIGURED{$setting}\n";
}

close(INFO);

exit $retcode;


#/*@@
#  @routine
#  @date       Fri Feb 19 19:53:48 1999
#  @author     Tom Goodale
#  @desc
#  Sets the environment for running the configure script.
#  @enddesc
#@@*/
sub SetConfigureEnv
{
  local($line_number, $commandline);
  # Set a default name for the configuration
  $ENV{"EXE"} = "cactus_$config";

  # Set variables from makefile command line first
  $commandline = $ENV{"MAKEFLAGS"};
  $line_number = 0;
#  while ($commandline =~ /^(.*)\s+(\w+)\s*=\s*([_+\-\.\w\\\/\s]*)\s*/)
  while ($commandline =~ /^(.*)\s*\b(\w+)\s*=\s*(.*?)\s*$/)
  {
    if ($2 ne 'options' && $2 ne 'VERBOSE')
    {
      print "Using configuration options from configure line\n"
        if (!$line_number);
      $line_number++;
      # Remember it for writing to config-info
      $option = AddQuotes($3);
      $CONFIGURED{$2} = $option;

      print "  Setting $2 to '$option'\n";
    }
    $commandline=$1;
    #  print "New commandline = <$commandline>\n";
  }
  print "End of options from configure line\n"
    if ($line_number);

  # Add variables from user configuration options file
  if($config_file)
  {
    # Turn path to options file to an absolute path
    $filename = $config_file;
    if($config_file =~ m:^/:)
    {
      # Do nothing
    }
    elsif($config_file =~ m:^~:)
    {
      $filename =~ s/^~/$ENV{"HOME"}/
    }
    else
    {
      $filename = "$top/$config_file";
    }

    # The user has specified a configuration file
    print "Adding configuration options from '$config_file'...\n";
    ParseOptionsFile ($filename, \%ENV, \%CONFIGURED);
    print "End of options from '$config_file'.\n";
  }

  # Add variables from either ${CACTUS_CONFIG_FILES} or a user default file
  if (defined $ENV{CACTUS_CONFIG_FILES})
  {
    foreach $file (split (':', $ENV{CACTUS_CONFIG_FILES}))
    {
      print "Adding configuration options from '$file'...\n";
      ParseOptionsFile ($file, \%ENV, \%CONFIGURED);
      print "End of options from '$file'.\n";
    }
  }
  elsif ($default_file)
  {
    print "Adding configuration options from user defaults in $default_file...\n";
    ParseOptionsFile ($default_file, \%ENV, \%CONFIGURED);
    print "End of options from user defaults.\n";
  }

  return %CONFIGURED;
}

sub ParseOptionsFile
{
  my($file, $env, $options) = @_;
  my($line_number);


  open(INFILE, "< $file") || die "Cannot open configuration file '$file'\n";

  $line_number = 0;
  while(<INFILE>)
  {
    $line_number++;

    #Ignore comments.
    s/\#(.*)$//g;

    #Remove spaces at end of lines
    s/\s*$//;
    s/\n//g;                # Different from chop...

    #Ignore blank lines
    next if (m:^\s*$:);

    # Match lines of the form
    #     keyword value
    # or  keyword = value
    if (/^\s*(\w+)[=\s]+(.*)\s*/)
    {
      # only set it if it wasn't already
      if(! $options->{$1})
      {
        print "  Setting $1 to '$2'\n";
        $env->{$1} = $2;
        # Remember it for writing to config-info
        $options->{$1} = AddQuotes($2);
      }
    }
    else
    {
      print "Could not parse configuration line $line_number...\n'$_'\n";
    }
  }
  close(INFILE);
}

sub DetermineConfigureCommand
{
  my($configure_command) = @_;


  $configure_command .= ' --build='  . $ENV{BUILD_MACHINE}
    if($ENV{BUILD_MACHINE});
  $configure_command .= ' --target=' . $ENV{TARGET_MACHINE}
    if($ENV{TARGET_MACHINE});
  $configure_command .= ' --host='   . $ENV{HOST_MACHINE}
    if($ENV{HOST_MACHINE});

  return $configure_command;
}

sub AddQuotes
{
  local($arg) = @_;

  if ($arg =~ /\\/)
  {
    $arg =~ s:\\::g;
#    $arg = "\'$arg\'";
  }

  # When we grab an arg off the MAKEFLAGS it has $s doubled.
  $arg =~ s/\$\$/\$/g;

  return $arg;
}
