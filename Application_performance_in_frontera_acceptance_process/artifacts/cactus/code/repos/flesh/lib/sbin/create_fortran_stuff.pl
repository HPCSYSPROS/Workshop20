#! /usr/bin/perl
#/*@@
#  @file      create_fortran_stuff.pl
#  @date      Tue Jan 12 09:52:35 1999
#  @author    Tom Goodale
#  @desc
#  Create the Fortran parameter stuff
#  @enddesc
#version $Header$
#@@*/

sub CreateFortranThornParameterBindings
{
  my($thorn, $rhparameter_db, $rhinterface_db) = @_;
  my($line);
  my(%these_parameters);
  my($implementation);
  my(@data);
  my(@file);
  my(%alias_names);
  my(%num_aliases);

  push(@file, '#define DECLARE_CCTK_PARAMETERS \\');

  # Generate all global parameters
  %these_parameters = &get_global_parameters($rhparameter_db);

  if((keys %these_parameters) > 0)
  {
    @data = &CreateFortranCommonDeclaration('cctk_params_global', \%these_parameters, $rhparameter_db);

    if (@data)
    {
      push(@file, join ("&&\\\n", @data) . "&&\\");
    }
  }

  # Generate all restricted parameters of this thorn
  %these_parameters = &GetThornParameterList($thorn, 'RESTRICTED', $rhparameter_db);

  if((keys %these_parameters > 0))
  {
    $implementation = $rhinterface_db->{"\U$thorn\E IMPLEMENTS"};

    @data = &CreateFortranCommonDeclaration("${implementation}rest", \%these_parameters, $rhparameter_db);

    if (@data)
    {
      push(@file, join ("&&\\\n", @data) . "&&\\");
    }
  }

  # Generate all private parameters of this thorn
  %these_parameters = &GetThornParameterList($thorn, 'PRIVATE', $rhparameter_db);

  if((keys %these_parameters > 0))
  {
    @data = &CreateFortranCommonDeclaration("${thorn}priv", \%these_parameters, $rhparameter_db);

    if (@data)
    {
      push(@file, join ("&&\\\n", @data) . "&&\\");
    }
  }

  # Parameters from friends

  # This number can be local to each thorn - it doesn't matter if
  # members of a common block get different names in different
  # thorns, especially if the variable isn't being used !
  $num_aliases = 0;

#  print "DEBUG ********************************************\n";
#  print "DEBUG thorn is $thorn\n";
  foreach $friend (split(' ',$rhparameter_db->{"\U$thorn\E SHARES implementations"}))
  {
#    print "DEBUG friend is $friend\n";

    # Determine which thorn provides this friend implementation
    $rhinterface_db->{"IMPLEMENTATION \U$friend\E THORNS"} =~ m:([^ ]*):;

    $friend_thorn = $1;

    %these_parameters = &GetThornParameterList($friend_thorn, 'RESTRICTED', $rhparameter_db);

    %alias_names = ();

    foreach $parameter (sort keys %these_parameters)
    {
#      print "DEBUG parameter is $parameter\n";
      my $foundit = 0;
      my $thornparam;
      my $name = "";

      foreach $thornparam (split(/\s+/,$rhparameter_db->{"\U$thorn SHARES $friend\E variables"}))
      {
#        print "DEBUG thorn parameter is $thornparam\n";
#        print "DEBUG       realname is " . $rhparameter_db->{"\U$thorn $thornparam\E realname"} ."\n";
          
        if($rhparameter_db->{"\U$thorn $thornparam\E realname"} =~ m/^$parameter$/i)
        {
#          print "DEBUG ... " . $rhparameter_db->{"\U$thorn $thornparam\E realname"} . "\n";
          $name = $thornparam; #$rhparameter_db->{"\U$thorn $thornparam\E realname"};
          last;
        }
      }

      if($name eq "")
      {
        $name = "CCTKH$num_aliases";
        $num_aliases++;
      }
      $alias_names{$parameter} = "$name";
    }

    @data = &CreateFortranCommonDeclaration("${friend}rest", \%these_parameters, $rhparameter_db, \%alias_names);

    if (@data)
    {
      push(@file, join ("&&\\\n", @data) . "&&\\");
    }
  }

  push(@file, '');
  push(@file, '');

  return join ("\n", @file);
}


sub CreateFortranCommonDeclaration
{
  my($common_block, $rhparameters, $rhparameter_db, $rhaliases) = @_;
  my($line,@data);
  my(%parameters);
  my($type, $type_string);
  my($definition);
  my($aliases);

  if(defined $rhaliases)
  {
    $aliases = scalar(keys %$rhaliases);
  }
  else
  {
    $aliases = 0;
  }

  # Create the data

  $definition = "COMMON /$common_block/";

  $sepchar = '';

  foreach $parameter (&order_params($rhparameters,$rhparameter_db))
  {
    my $type = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E type"};

    my $type_string = &get_fortran_type_string($type);

    my $array_size = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E array_size"};

    my $suffix = '';

    if($array_size)
    {
      $suffix = "($array_size)";
    }

    my $name;

    if($aliases)
    {
      $name = $rhaliases->{$parameter};
    }
    else
    {
      $name = $rhparameter_db->{"\U$rhparameters->{$parameter} $parameter\E realname"};
    }
      
    $line = "CCTK_DECLARE($type_string,$name,$suffix)";
    $definition .= "$sepchar$name";

    push(@data, $line);

    $sepchar = ', ';
  }

  push(@data, $definition);

  return @data;
}


sub get_fortran_type_string
{
  my($type) = @_;
  my($type_string);


  if($type eq 'KEYWORD' ||
     $type eq 'STRING'  ||
     $type eq 'SENTENCE')
  {
    $type_string = 'CCTK_STRING';
  }
  elsif($type eq 'BOOLEAN' ||
        $type eq 'INT')
  {
    $type_string = 'CCTK_INT';
  }
  elsif($type eq 'REAL')
  {
    $type_string = 'CCTK_REAL';
  }
  else
  {
    $message = "Unknown parameter type '$type'";
    &CST_error(0,$message,'',__LINE__,__FILE__);
  }

  return $type_string;

}

1;
