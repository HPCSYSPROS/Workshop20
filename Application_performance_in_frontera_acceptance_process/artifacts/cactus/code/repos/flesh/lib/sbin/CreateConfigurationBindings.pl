#/*@@
#  @file      CreateConfigurationBindings.pl
#  @date      Thu Mar 25 14:25:13 2004
#  @author    Yaakoub Y El-Khamra
#  @desc
#             New Configuration.ccl script processing
#  @enddesc
#  @version   $Header$
#@@*/

require "$sbin_dir/CSTUtils.pl";

#/*@@
#  @routine    CreateConfigurationBindings
#  @date       Thu Mar 25 14:25:13 2004
#  @author     Yaakoub Y El-Khamra
#  @desc
#              Creates the configuration bindings.
#  @enddesc
#@@*/
sub CreateConfigurationBindings
{
  my($bindings_dir, $cfg, $thorns)=@_;
  my($field, $providedcap, $thorn, $temp,$defs,$incs,$deps);
  my(%linker_thorns, %linker_cfg, $linker_list, $linkerdirs, $linkerlibs);

  if(! $build_dir)
  {
    $build_dir = "$bindings_dir/build";
  }

  if(! -d $bindings_dir)
  {
    mkdir("$bindings_dir", 0755) || die "Unable to create $bindings_dir";
  }

  $start_dir = `pwd`;

  chdir $bindings_dir;

  if(! -d 'Configuration')
  {
    mkdir('Configuration', 0755) || die 'Unable to create Configuration directory';
  }

  if(! -d 'include')
  {
    mkdir('include', 0755) || die 'Unable to create include directory';
  }

  chdir 'Configuration';

  if(! -d "Capabilities")
  {
    mkdir("Capabilities", 0755) || die "Unable to create Capabilities directory";
  }
  if(! -d 'Thorns')
  {
    mkdir('Thorns', 0755) || die "Unable to create Thorns directory";
  }

  # Put all the provided capabilities where they belong
  foreach my $thorn (sort keys %thorns)
  {
      # We know that all the requirements have been satisfied, so all
      # we need to do is put the provides where they belong.
      if ($cfg->{"\U$thorn\E PROVIDES"})
      {
          foreach my $providedcap (sort split (' ', $cfg->{"\U$thorn\E PROVIDES"}))
          {
              die if $providedcap !~ m{^[A-Za-z0-9_.]+$};
              
              my $defs = '';
              my $incs = '';
              my $deps = '';

              # Add requirements recursively
              foreach my $requiredcap (sort split (' ', $cfg->{"\U$thorn\E REQUIRES"}))
              {
                  next if $requiredcap eq $providedcap;
                  $defs .= "include $bindings_dir/Configuration/Capabilities/make.\U$requiredcap\E.defn\n";
                  $incs .= "#include \"cctki_\U$requiredcap\E.h\"\n";
                  $deps .= "include $bindings_dir/Configuration/Capabilities/make.\U$requiredcap\E.deps\n";
              }
              
              # Put include_dirs and make.definition in one file:
              # make.capability.defn
              my $incdir = $cfg->{"\U$thorn $providedcap\E INCLUDE_DIRECTORY"};
              if ($incdir)
              {
                  $defs .= "INC_DIRS += $incdir\n";
              }
              my $incdir_f = $cfg->{"\U$thorn $providedcap\E INCLUDE_DIRECTORY_FORTRAN"};
              if ($incdir_f)
              {
                  $defs .= "INC_DIRS_F += $incdir_f\n";
              }
              elsif ($incdir)
              {
                  $defs .= "INC_DIRS_F += $incdir\n";
              }
              
              my $makedef = $cfg->{"\U$thorn $providedcap\E MAKE_DEFINITION"};
              if ($makedef)
              {
                  $defs .= $makedef;
              }
              
              $defs .= "HAVE_CAPABILITY_$providedcap = 1\n";
              
              # Put INCLUDE and DEFINE in one file: capability.h
              my $inc = $cfg->{"\U$thorn $providedcap\E INCLUDE"};
              if ($inc)
              {
                  # Prepend #include
                  $inc =~ s/^(.*)/#include $1/gm;
                  $incs .= $inc;
              }
              
              my $def = $cfg->{"\U$thorn $providedcap\E DEFINE"};
              if ($def)
              {
                  # Prepend #define
                  $def =~ s/^(.*)/#define $1/gm;
                  $incs .= $def;
              }
              
              $incs .= "#define HAVE_CAPABILITY_$providedcap 1\n";
              
              # Put make.capability.deps in one file:
              # make.capabiltiy.deps
              my $makedep = $cfg->{"\U$thorn $providedcap\E MAKE_DEPENDENCY"};
              if ($makedep)
              {
                  $deps .= $makedep;
              }
              
              &WriteFile("Capabilities/make.\U$providedcap\E.defn",\$defs);
              &WriteFile("Capabilities/cctki_\U$providedcap\E.h",\$incs);
              &WriteFile("Capabilities/make.\U$providedcap\E.deps",\$deps);
              
              # Create a list of thorns that provide a capability with
              # a library (or library search paths)
              if ($cfg->{"\U$thorn $providedcap\E LIBRARY"})
              {
                  $linker_thorns{"$thorn"} = $thorn;
                  #$linker_cfg{"\U$thorn\E USES"} = $cfg->{"\U$thorn\E USES THORNS"};
              }
              
              if ($cfg->{"\U$thorn $providedcap\E LIBRARY_DIRECTORY"})
              {
                  $linker_thorns{"$thorn"} = $thorn;
                  #$linker_cfg{"\U$thorn\E USES"} = $cfg->{"\U$thorn\E USES THORNS"};
              }
          }
      }
  }

  # here we add the files to the thorns that require capabilities
  foreach $thorn (sort keys %thorns)
  {
    # we know that all the requirements have been satisfied
    # so all we need to do is make references to the capabilities
    # from the requirements since we can have multiple provides,
    # we make each capability separate

    $defs = '';
    $incs = '';
    $deps = '';

    if ($cfg->{"\U$thorn\E REQUIRES"})
    {
      foreach $requiredcap (sort split (' ', $cfg->{"\U$thorn\E REQUIRES"}))
      {
        # put reference to provided capability
        $defs .= "include $bindings_dir/Configuration/Capabilities/make.\U$requiredcap\E.defn\n";
        $incs .= "#include \"../Capabilities/cctki_\U$requiredcap\E.h\"\n";
        $deps .= "include $bindings_dir/Configuration/Capabilities/make.\U$requiredcap\E.deps\n";
      }
    }

    if ($cfg->{"\U$thorn\E REQUIRES"})
    {
      # write everything to file
      # (write the files even if they are empty)
      &WriteFile("./Thorns/make.$thorn.defn",\$defs);
      &WriteFile("./Thorns/cctki_$thorn.h",\$incs);
      &WriteFile("./Thorns/make.$thorn.deps",\$deps);
    }
    else
    {
      # remove the files
      # (we cannot have old files staying around)
      unlink "./Thorns/make.$thorn.defn";
      unlink "./Thorns/cctki_$thorn.h";
      unlink "./Thorns/make.$thorn.deps";
    }

  }

  # Sort the linker thorns
  $linkerdirs = 'LIBDIRS +=';
  $linkerlibs = 'LIBS +=';

  $linker_list = &TopoSort(\%linker_thorns, \%linker_cfg, $cfg);
  foreach $thorn (split (' ', $linker_list))
  {
    foreach $providedcap (sort split (' ', $cfg->{"\U$thorn\E PROVIDES"}))
    {
      $linkerdirs .= ' ' . $cfg->{"\U$thorn $providedcap\E LIBRARY_DIRECTORY"};
      $linkerlibs .= ' ' . $cfg->{"\U$thorn $providedcap\E LIBRARY"};
    }
  }
  $temp = $linkerdirs . "\n" . $linkerlibs . "\n";
  &WriteFile("make.link",\$temp);

  # write cctk_Capabilities.h file to bindings/include
  # this file adds the if_i_am_thorn stuff
  foreach $thorn (sort keys %thorns)
  {
    $temp = '';
    if ($cfg->{"\U$thorn\E REQUIRES"})
    {
      $temp .= "#include \"../Configuration/Thorns/cctki_$thorn.h\"\n";
    }
    &WriteFile("../include/$thorn/cctk_Capabilities.h",\$temp);
  }
  &WriteFile("../include/CactusBindings/cctk_Capabilities.h",
             "#include \"../Configuration/Thorns/cctki_Cactus.h\"\n");
}

return 1;
