#/*@@
#  @file      CreateParameterBindings.pl
#  @date      Sun Jul 25 00:52:36 1999
#  @author    Tom Goodale
#  @desc
#             Parameter binding stuff
#  @enddesc
#  @version   $Header$
#@@*/

#/*@@
#  @routine    CreateParameterBindings
#  @date       Thu Jan 28 15:27:16 1999
#  @author     Tom Goodale
#  @desc
#  Create the bindings used for the parameters.
#  @enddesc
#@@*/

sub CreateParameterBindings
{
  my($bindings_dir, $rhparameter_db, $rhinterface_db) = @_;
  my($start_dir);
  my($line);
  my(%these_parameters);
  my($implementation, $thorn);
  my($files);
  my(%routines);
  my($structure, %structures);
  my(%header_files);


  if(! -d $bindings_dir)
  {
    mkdir("$bindings_dir", 0755) || die "Unable to create $bindings_dir";
  }
  $start_dir = `pwd`;

  chdir $bindings_dir;

  if(! -d 'Parameters')
  {
    mkdir('Parameters', 0755) || die 'Unable to create Parameters directory';
  }

  if(! -d 'include')
  {
    mkdir('include', 0755) || die 'Unable to create include directory';
  }

  # Generate all global parameters
  %these_parameters = &get_global_parameters($rhparameter_db);

  $dataout = &CreateParameterBindingFile('CCTK_BindingsParametersGlobal', 'GLOBAL_PARAMETER_STRUCT', \%these_parameters, $rhparameter_db);

  &WriteFile('Parameters/Global.c',\$dataout);

  $files = 'Global.c';
  $structures{'GLOBAL_PARAMETER_STRUCT'} = 'cctk_params_global';

  # Generate the global data header file

  $dataout = &CreateCStructureParameterHeader('CCTK_BindingsParametersGlobal', 'GLOBAL_PARAMETER_STRUCT', \%these_parameters, $rhparameter_db);

  &WriteFile('include/ParameterCGlobal.h',\$dataout);

  $header_files{'GLOBAL'} = 'ParameterCGlobal.h';

  # Generate all restricted parameters
  foreach $implementation (split(' ',$rhinterface_db->{'IMPLEMENTATIONS'}))
  {
    $rhinterface_db->{"IMPLEMENTATION \U$implementation\E THORNS"} =~ m:([^ ]+):;

    $thorn = $1;

    %these_parameters = &GetThornParameterList($thorn, 'RESTRICTED', $rhparameter_db);

    if((keys %these_parameters > 0))
    {
#      $dataout = &CreateParameterBindingFile("CCTK_BindingsParameters${implementation}_restricted", "RESTRICTED_\U$implementation\E_STRUCT", \%these_parameters, $rhparameter_db);
#
#      &WriteFile("Parameters/${implementation}_restricted.c",\$dataout);
#
#      $files .= " ${implementation}_restricted.c";
      $routines{"CCTK_BindingsParameters${implementation}_restricted"} = "$implementation";
      $structures{"RESTRICTED_\U$implementation\E_STRUCT"} = "${implementation}rest";

      # Generate the data header file
      $dataout = &CreateCStructureParameterHeader("CCTK_BindingsParameters${implementation}_restricted", "RESTRICTED_\U$implementation\E_STRUCT", \%these_parameters, $rhparameter_db);

      &WriteFile("include/ParameterCRestricted\U$implementation\E.h",\$dataout);

      $header_files{"\U$implementation\E RESTRICTED"} = "ParameterCRestricted\U$implementation\E.h";
    }
  }

  # Generate all private parameters
  foreach $thorn (split(' ',$rhinterface_db->{'THORNS'}))
  {
    %these_parameters = &GetThornParameterList($thorn, 'PRIVATE', $rhparameter_db);

    if((keys %these_parameters > 0))
    {
#      $dataout = &CreateParameterBindingFile("CCTK_BindingsParameters${thorn}_private", "PRIVATE_\U$thorn\E_STRUCT", \%these_parameters, $rhparameter_db);
#
#      &WriteFile("Parameters/${thorn}_private.c",\$dataout);
#
#      $files .= " ${thorn}_private.c";
      $routines{"CCTK_BindingsParameters${thorn}_private"} = "$thorn";

      # Generate the data header file
      $dataout = &CreateCStructureParameterHeader("CCTK_BindingsParameters${thorn}_private", "PRIVATE_\U$thorn\E_STRUCT", \%these_parameters, $rhparameter_db);

      $structures{"PRIVATE_\U$thorn\E_STRUCT"} = "${thorn}priv";

      &WriteFile("include/ParameterCPrivate\U$thorn\E.h",\$dataout);

      $header_files{"\U$thorn\E PRIVATE"} = "ParameterCPrivate\U$thorn\E.h";
    }
  }


  @thorns = split(' ',$rhinterface_db->{'THORNS'});
  @data = ();
  push(@data, map ('extern int CCTKi_BindingsCreate' . $_ . 'Parameters(void);', @thorns));
  push(@data, '');
  push(@data, map ('extern int CCTKi_Bindings' . $_ . 'ParameterExtensions(void);', @thorns));
  push(@data, '');
  push(@data, '');

  push(@data, 'int CCTKi_BindingsParametersInitialise(void);');
  push(@data, 'int CCTKi_BindingsParametersInitialise(void)');
  push(@data, '{');

  push(@data, map ('  CCTKi_BindingsCreate' . $_ . 'Parameters();', @thorns));
  push(@data, '');
  push(@data, map ('  CCTKi_Bindings' . $_ . 'ParameterExtensions();', @thorns));
  push(@data, '');

  push(@data, '  return 0;');
  push(@data, '}');
  push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline

  $dataout = join ("\n", @data);
  &WriteFile('Parameters/BindingsParameters.c',\$dataout);

  $newfilelist = NewParamStuff($rhparameter_db, $rhinterface_db);

  $dataout = "SRCS = BindingsParameters.c $files $newfilelist";
  &WriteFile('Parameters/make.code.defn',\$dataout);

  # Create the appropriate thorn parameter headers

  foreach $thorn (split(' ',$rhinterface_db->{'THORNS'}))
  {
    $dataout = &CreateFortranThornParameterBindings($thorn, $rhparameter_db, $rhinterface_db);
    mkdir("include/$thorn");
    &WriteFile("include/${thorn}/FParameters.h",\$dataout);

    $implementation = $rhinterface_db->{"\U$thorn\E IMPLEMENTS"};

    @data = ();
    push(@data, '/*@@');
    push(@data, "   \@header  ${thorn}/CParameters.h");
    push(@data, '   @author  Automatically generated by CreateParameterBindings.pl');
    push(@data, '   @desc');
    push(@data, "            Declares parameters of thorn $thorn");
    push(@data, '   @enddesc');
    push(@data, ' @@*/');
    push(@data, '');
    push(@data, '');

    push(@data, "#ifndef _\U$thorn\E_PARAMETERS_H_");
    push(@data, "#define _\U$thorn\E_PARAMETERS_H_ 1");
    push(@data, '');
    push(@data, '#include "CParameterStructNames.h"');

    if($header_files{'GLOBAL'})
    {
      push(@data, "#include \"$header_files{'GLOBAL'}\"");
    }

    if($header_files{"\U$implementation\E RESTRICTED"})
    {
      push(@data, "#include \"" . $header_files{"\U$implementation\E RESTRICTED"} . "\"");
    }

    push(@data, "#include \"" . $header_files{"\U$thorn\E PRIVATE"} . "\"")
      if($header_files{"\U$thorn\E PRIVATE"});

    foreach $friend (split(' ',$rhparameter_db->{"\U$thorn\E SHARES implementations"}))
    {
      push(@data, "#include \"ParameterCRestricted\U$friend\E.h\"");
    }
    push(@data, '');

    push(@data, '#define DECLARE_CCTK_PARAMETERS \\');
    push(@data, '  DECLARE_GLOBAL_PARAMETER_STRUCT_PARAMS \\')
      if($header_files{'GLOBAL'});
    push(@data, "  DECLARE_RESTRICTED_\U$implementation\E_STRUCT_PARAMS \\")
      if($header_files{"\U$implementation\E RESTRICTED"});
    push(@data, "  DECLARE_PRIVATE_\U$thorn\E_STRUCT_PARAMS \\")
      if($header_files{"\U$thorn\E PRIVATE"});

    # double this loop, add #defines for RESTRICTED_STRUCT.$realname, set these via #ifdef, use these #defines

    my @data2;
    foreach $friend (split(' ',$rhparameter_db->{"\U$thorn\E SHARES implementations"}))
    {
      $rhinterface_db->{"IMPLEMENTATION \U$friend\E THORNS"} =~ m:([^ ]*):;
      $friend_thorn = $1;

      foreach $parameter (split(' ',$rhparameter_db->{"\U$thorn SHARES $friend\E variables"}))
      {
        my $realname = $rhparameter_db->{"\U$thorn $parameter\E realname"};

        $type = $rhparameter_db->{"\U$friend_thorn $realname\E type"};
        $array_size = $rhparameter_db->{"\U$friend_thorn $realname\E array_size"};
        $type_string = &get_c_type_string($type,$realname);

        my $varprefix = '';

        if($array_size)
        {
          $varprefix = ' const *';
        }

        #push(@data, "  CCTK_DECLARE_INIT ($type_string$varprefix const, $parameter, RESTRICTED_\U$friend\E_STRUCT.$realname); \\");
        push(@data, "  CCTK_DECLARE_INIT ($type_string$varprefix const, $parameter, CCTK_PARAMETER__${friend_thorn}__${realname}); \\");
        push(@data2, 
             "#ifndef CCTK_PARAMETER__${friend_thorn}__${realname}\n" .
             "#  define CCTK_PARAMETER__${friend_thorn}__${realname} RESTRICTED_\U$friend\E_STRUCT.$realname\n" .
             "#endif");
      }
    }

    push(@data, '');
    push(@data, @data2);
    push(@data, '');
    push(@data, "#endif  /* _\U$thorn\E_PARAMETERS_H_ */");
    push(@data, "\n");  # workaround for perl 5.004_04 to add a trailing newline

    $dataout = join ("\n", @data);
    &WriteFile("include/${thorn}/CParameters.h",\$dataout);
  }

# Write this one to a temporary file and read it back in
# Can probably do this better

  open(OUT, "| $ENV{'PERLINTERP'} $cctk_home/lib/sbin/c_file_processor.pl $top/config-data > include/CParameterStructNames_temp.h") || die 'Cannot create CParameterStructNames.h by running c_file_processor.pl';

  foreach $structure (sort keys %structures)
  {
    print OUT "#define $structure CCTK_FORTRAN_COMMON_NAME($structures{$structure})\n";
  }
  print OUT "\n";

  close OUT;

  open(IN,'< include/CParameterStructNames_temp.h');
  $dataout = join ('', <IN>);
  close IN;

  &WriteFile('include/CParameterStructNames.h',\$dataout);

  foreach $thorn (split(' ',$rhinterface_db->{'THORNS'}))
  {
    @data = ();
    push(@data, '/* get the CCTK datatype definitions */');
    push(@data, '#include "cctk_Types.h"');
    push(@data, '');
    push(@data, '#ifdef CCODE');
    push(@data, "#include \"${thorn}/CParameters.h\"");
    push(@data, '#elif FCODE');
    push(@data, "#include \"${thorn}/FParameters.h\"");
    push(@data, '#endif');
    push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline

    $dataout = join ("\n", @data);
    &WriteFile("include/${thorn}/cctk_Parameters.h",\$dataout);
  }


  chdir $start_dir;
}


sub NewParamStuff
{
  my($rhparameter_db, $rhinterface_db) = @_;
  my($line);
  my(%these_parameters);
  my($implementation, $thorn);
  my($files);
  my(%routines);
  my($structure, %structures);
  my(%header_files);
  my($block);
  my($filelist);
  my(@creationdata);
  my(@extensiondata);
  my(@accumulationdata);
  my(@data);

  foreach $thorn (split(' ',$rhinterface_db->{'THORNS'}))
  {
    $imp = $rhinterface_db->{"\U$thorn\E IMPLEMENTS"};

    push(@data, '/*@@');
    push(@data, "   \@file    ${thorn}_Parameters.c");
    push(@data, '   @author  Automatically generated by CreateParameterBindings.pl');
    push(@data, '   @desc');
    push(@data, '            Creates/extends parameters for this thorn');
    push(@data, '   @enddesc');
    push(@data, ' @@*/');
    push(@data, '');
    push(@data, '');

    push(@data, '#include <stdio.h>');
    push(@data, '#include <stdlib.h>');
    push(@data, '#include <stdarg.h>');
    push(@data, '');
    push(@data, '#include "cctk_Config.h"');
    push(@data, '#include "cctk_Constants.h"');
    push(@data, '#include "ParameterBindings.h"');
    push(@data, '#include "CParameterStructNames.h"');

    foreach $block ('GLOBAL', 'RESTRICTED', 'PRIVATE')
    {
      %these_parameters = &GetThornParameterList($thorn, $block, $rhparameter_db);

      if((keys %these_parameters > 0))
      {
        if($block eq 'GLOBAL')
        {
          push(@data, '#include "ParameterCGlobal.h"');
        }
#        elsif($block eq 'RESTRICTED')
#        {
#          push(@data, "#include \"ParameterCRestricted\U$imp\E.h\"");
#        }
#        elsif($block eq 'PRIVATE')
#        {
#          push(@data, "#include \"ParameterCPrivate\U$thorn\E.h\"");
#        }
#        else
#        {
#          die 'Internal error';
#        }

#        print "Generating $block parameters for $thorn, providing $imp\n";
        push(@creationdata,&CreateParameterRegistrationStuff($block, $thorn, $imp, $rhparameter_db, %these_parameters));
        push(@accumulationdata,&CreateParameterAccumulationStuff($block, $thorn, $rhparameter_db, %these_parameters));
      }
    }


    # Now the parameter extensions
    foreach $block (split(' ',$rhparameter_db->{"\U$thorn\E SHARES implementations"}))
    {
      push(@data, "#include \"ParameterCRestricted\U$block\E.h\"");

#      print "Generating $block extension from $thorn\n";
      push(@extensiondata,&CreateParameterExtensionStuff($block, $thorn, $rhparameter_db));
    }

    # definition of private and restricted parameters follows right here
    push(@data, '');
    push(@data, "/* structure containing all private parameters of thorn $thorn */");
    %these_parameters = &GetThornParameterList($thorn, 'PRIVATE', $rhparameter_db);
    $block = &CreateParameterBindingFile('', "PRIVATE_\U$thorn\E_STRUCT", \%these_parameters, $rhparameter_db);
    push(@data, $block);
    push(@data, "/* structure containing all restricted parameters of thorn $thorn */");
    %these_parameters = &GetThornParameterList($thorn, 'RESTRICTED', $rhparameter_db);
    $block = &CreateParameterBindingFile('', "RESTRICTED_\U$imp\E_STRUCT", \%these_parameters, $rhparameter_db);
    push(@data, $block);

    push(@data, "int CCTKi_BindingsCreate${thorn}Parameters(void);");
    push(@data, "int CCTKi_BindingsCreate${thorn}Parameters(void)");
    push(@data, '{');

    push(@data, @creationdata);

    push(@data, '  return 0;');
    push(@data, '}');

    push(@data, '');
    push(@data, "int CCTKi_Bindings${thorn}ParameterExtensions(void);");
    push(@data, "int CCTKi_Bindings${thorn}ParameterExtensions(void)");
    push(@data, '{');

    push(@data, @extensiondata);
    push(@data, @accumulationdata);

    push(@data, '  return 0;');
    push(@data, '}');
    push(@data, "\n");  # workaround for perl 5.004_04 to add a trailing newline

    $dataout = join ("\n", @data);
    &WriteFile("Parameters/${thorn}_Parameters.c",\$dataout);

    @data=();
    @creationdata=();
    @extensiondata=();
    @accumulationdata=();

    $filelist .= " ${thorn}_Parameters.c";
  }

  return $filelist;
}


sub CreateParameterRegistrationStuff
{
  my($block, $thorn, $imp, $rhparameter_db, %these_parameters) = @_;
  my(@data);
  my($line);
  my($structure, $type, $n_ranges);

  if($block eq 'GLOBAL')
  {
    $structure='GLOBAL_PARAMETER_STRUCT';
  }
  elsif($block eq 'RESTRICTED')
  {
    $structure="RESTRICTED_\U$imp\E_STRUCT";
  }
  elsif($block eq 'PRIVATE')
  {
    $structure = "PRIVATE_\U$thorn\E_STRUCT";
  }
  else
  {
    die 'Internal error';
  }

#  print "Thorn is $thorn\n";
#  print "Structure is $structure\n";

  foreach $parameter (sort keys %these_parameters)
  {

#    print "This param is $parameter\n";

    my $type = $rhparameter_db->{"\U$thorn $parameter\E type"};

#    print "Type is $type\n";

    my $n_ranges = $rhparameter_db->{"\U$thorn $parameter\E ranges"};

#    print "N_ranges is $n_ranges\n";

    my $quoted_default = $rhparameter_db->{"\U$thorn $parameter\E default"};

#    $quoted_default =~ s:\"::g;  The database now strips all unescaped quotes.

    # Set steerable details
    my $steerable = $rhparameter_db->{"\U$thorn $parameter\E steerable"};
    if ($steerable =~ /^never$/i || $steerable =~/^$/)
    {
      $steerable_type = 'CCTK_STEERABLE_NEVER';
    }
    elsif ($steerable =~ /^always$/i)
    {
      $steerable_type = 'CCTK_STEERABLE_ALWAYS';
    }
    elsif ($steerable =~ /^recover$/i)
    {
      $steerable_type = 'CCTK_STEERABLE_RECOVER';
    }
    else
    {
      $message = "Illegal steerable type ($steerable) for parameter $parameter in $thorn";
      &CST_error(0,$message,'',__LINE__,__FILE__);
    }

    my $realname = $rhparameter_db->{"\U$thorn $parameter\E realname"};

    my $array_size = $rhparameter_db->{"\U$thorn $parameter\E array_size"};

    my $dereference = '';

    if(! $array_size)
    {
      $array_size = 0;
      $dereference = '&'
    }

    my $accumulator_expression = $rhparameter_db->{"\U$thorn $parameter\E accumulator-expression"};

    if(! $accumulator_expression)
    {
      $accumulator_expression = 'NULL';
    }
    else
    {
      $accumulator_expression = "\"$accumulator_expression\"";
    }

#    my $accumulator_base = $rhparameter_db->{"\U$thorn $parameter\E accumulator-base"};
#
#    if(! $accumulator_expression)
#    {
#      $accumulator_base = 'NULL';
#    }

    $line="  CCTKi_ParameterCreate(\"$parameter\",\n" .
          "                        \"$thorn\",\n" .
          "                        \"$type\",\n" .
          "                        \"$block\",\n" .
          "                        $steerable_type,\n" .
          "                        " . $rhparameter_db->{"\U$thorn $parameter\E description"} . ",\n" .
          "                        \"" . $quoted_default . "\",\n" .
          "                        $dereference($structure.$realname),\n" .
          "                        $array_size,\n" .
          "                        $accumulator_expression,\n" .
#          "                        $accumulator_base,\n" .
          "                        $n_ranges";

    for($range=1; $range <= $n_ranges; $range++)
    {
      $quoted_range = $rhparameter_db->{"\U$thorn $parameter\E range $range range"};
      $range_description = $rhparameter_db->{"\U$thorn $parameter\E range $range description"};

      if($range_description !~ m:\":)
      {
        $range_description = "\"$range_description\"";
      }

      $range_description =~ s:,$::;

      #$quoted_range =~ s:\":\\\":g;
      $quoted_range =~ s:\"::g;
      $quoted_range =~ s:^\s*::;
      $quoted_range =~ s:\s*$::;

      # escape all backslashes so that they aren't treated
      # as the beginning of an escape sequence in C strings
      $quoted_range =~ s:\\:\\\\:g;

      $line .= ",\n                        \"".$quoted_range."\",$range_description";
    }

    $line .= ");\n";

    push(@data, $line);
  }

  return @data;
}


sub CreateParameterExtensionStuff
{
  my($block, $thorn, $rhparameter_db) = @_;
  my(@data);
  my($line);
  my($structure, $type, $n_ranges, $range, $quoted_range, $range_description);

#  print "Extending $block from $thorn\n";

  foreach $parameter (split(' ',$rhparameter_db->{"\U$thorn\E SHARES \U$block\E variables"}))
  {
    my $realname = $rhparameter_db->{"\U$thorn $parameter\E realname"};

    $n_ranges = $rhparameter_db->{"\U$thorn $parameter\E ranges"};

    for($range=1; $range <= $n_ranges; $range++)
    {
      $quoted_range = $rhparameter_db->{"\U$thorn $parameter\E range $range range"};
      $range_description = $rhparameter_db->{"\U$thorn $parameter\E range $range description"};

      if($range_description !~ m:\":)
      {
        $range_description = "\"$range_description\"";
      }

      #$quoted_range =~ s:\":\\\":g;
      $quoted_range =~ s:\"::g;
      $quoted_range =~ s:^\s*::;
      $quoted_range =~ s:\s*$::;

      push(@data, "  CCTKi_ParameterAddRange(\"$block\",");
      push(@data, "                          \"$realname\",");
      push(@data, "                          \"$thorn\",");
      push(@data, "                          \"$quoted_range\",");
      push(@data, "                         $range_description);");
      push(@data, '');

#      print "Adding \"$quoted_range\" to $parameter\n";
    }
  }

  return @data;
}

sub CreateParameterAccumulationStuff
{
  my($block, $thorn, $rhparameter_db, %these_parameters) = @_;
  my(@data);

  @data = ();

  foreach $parameter (sort keys %these_parameters)
  {
    my $accumulator_base = $rhparameter_db->{"\U$thorn $parameter\E accumulator-base"};

    if($accumulator_base)
    {

#      print "accumulator_base = $accumulator_base\n";

      $accumulator_base =~ m/([^:]+)::(.+)/;

      my $importhorn = $1;
      my $accparam   = $2;

      push(@data, "  CCTKi_ParameterAccumulatorBase(\"$thorn\",");
      push(@data, "                          \"$parameter\",");
      push(@data, "                          \"$importhorn\",");
      push(@data, "                          \"$accparam\");");
      push(@data, "");
    }
  }

  return @data;
}


1;
