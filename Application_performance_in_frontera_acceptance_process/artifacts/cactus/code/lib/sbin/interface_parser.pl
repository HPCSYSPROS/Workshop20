#! /usr/bin/perl -w

#/*@@
#  @file    interface_parser.pl
#  @date    Wed Sep 16 15:07:11 1998
#  @author  Tom Goodale
#  @desc
#           Parses interface.ccl files
#  @enddesc
#  @version $Header$
#@@*/

#/*@@
#  @routine    create_interface_database
#  @date       Wed Sep 16 15:07:11 1998
#  @author     Tom Goodale
#  @desc
#  Creates a database of all the interfaces
#  @enddesc
#@@*/

sub create_interface_database
{
  my($n_system,@inargs) = @_;
  my(%system_database);
  my(%thorns, @thorns);
  my(%interface_data);

  %system_database = @inargs[0..2*$n_system-1];
  %thorns = @inargs[2*$n_system..$#inargs];
  @thorns = sort keys %thorns;

  #  Loop through each  thorn's interface file.
  foreach my $thorn (@thorns)
  {
    print "   $thorn\n";
    #       Get the arrangement name for the thorn
    $thorns{$thorn} =~ m:.*/arrangements/([^/]*)/[^/]*:;
    my $arrangement = $1;

    #       Read the data
    my @indata = &read_file("$thorns{$thorn}/interface.ccl");

    #       Get the interface data from it
    &parse_interface_ccl($arrangement, $thorn, \@indata, \%interface_data);

    &PrintInterfaceStatistics($thorn, \%interface_data);
  }

  &cross_index_interface_data (\@thorns, \%interface_data);

  return %interface_data;
}



sub cross_index_interface_data
{
  my($thorns_ref, $interface_data_ref) = @_;
  my(%implementations);
  my($implementation);
  my(%ancestors);
  my(%friends);
  my($thorn,$thorn_implements,$ancestor_imp,$thorn_ancestor,$message,$hint);

  foreach $thorn (@$thorns_ref)
  {
    $implementation = $interface_data_ref->{"\U$thorn\E IMPLEMENTS"};
    if($implementation =~ m:^\s*$:)
    {
      $message = "Thorn $thorn doesn't specify an implementation";
      $hint = "All compiled thorns must specify an implementation in their interface.ccl file with the format IMPLEMENTS: <implementation>";
      &CST_error(0,$message,$hint,__LINE__,__FILE__);
      next;
    }

    # Put if statement around this to prevent perl -w from complaining.
    if($interface_data_ref->{"IMPLEMENTATION \U$implementation\E THORNS"})
    {
      $interface_data_ref->{"IMPLEMENTATION \U$implementation\E THORNS"} .= "$thorn ";
    }
    else
    {
      $interface_data_ref->{"IMPLEMENTATION \U$implementation\E THORNS"} = "$thorn ";
    }

    $implementations{"\U$implementation\E"} = "$implementation";
  }

  $interface_data_ref->{"THORNS"} = join(" ", @$thorns_ref);

  foreach $implementation (sort keys %implementations)
  {

    # Put if statement around this to prevent perl -w from complaining.
    if($interface_data_ref->{"IMPLEMENTATIONS"})
    {
      $interface_data_ref->{"IMPLEMENTATIONS"} .= $implementations{"\U$implementation\E"} . " ";
    }
    else
    {
      $interface_data_ref->{"IMPLEMENTATIONS"} = $implementations{"\U$implementation\E"} . " ";
    }

    &check_implementation_consistency($implementation, $interface_data_ref);

    my %ancestors = ();
    &get_implementation_ancestors($implementation, $interface_data_ref, \%ancestors);

    $interface_data_ref->{"IMPLEMENTATION \U$implementation\E ANCESTORS"} = join(" ",(sort keys %ancestors));

    $interface_data_ref->{"IMPLEMENTATION \U$implementation\E FRIENDS"} = &get_friends_of_me($implementation, \%implementations, $interface_data_ref);

  }

  # Create Hash table with thorns as ancestors
  foreach $thorn (@$thorns_ref)
  {
    $thorn_implements = $interface_data_ref->{"\U$thorn\E IMPLEMENTS"};
    foreach $ancestor_imp ( split(' ', $interface_data_ref->{"\U$thorn INHERITS\E"}))
    {
      next if($ancestor_imp eq '');
      $thorn_ancestor{uc($thorn)} .= $interface_data_ref->{"IMPLEMENTATION \U$ancestor_imp\E THORNS"}. ' ';
    }
  }

  # Call find_dep_cycles to find and report cycles
  $message = &find_dep_cycles(%thorn_ancestor);
  if ("" ne  $message)
  {
   $message  =~ s/^\s*//g;
   $message  =~ s/\s*$//g;
   $message =~ s/\s+/->/g;
   $message = "Found a cyclic dependency in implementation inheritance: ".$message."\n";
   &CST_error(0,$message,$hint,__LINE__,__FILE__);
  }

  foreach $thorn (@$thorns_ref)
  {
    &check_interface_consistency($thorn, $interface_data_ref);
  }

  foreach $implementation (sort keys %implementations)
  {
    my %friends = ();
    &get_implementation_friends($implementation, $interface_data_ref, \%friends);
    $interface_data_ref->{"IMPLEMENTATION \U$implementation\E FRIENDS"} = join(" ",(sort keys %friends));
  }
}


sub get_friends_of_me
{
  my ($implementation, $implementations_ref, $interface_data_ref) = @_;
  my @friends = ();

  foreach my $other_implementation (sort keys %$implementations_ref)
  {

    $interface_data_ref->{"IMPLEMENTATION \U$other_implementation\E THORNS"} =~ m:(\w+):;

    my $thorn = $1;

    foreach my $friend (split(" ", $interface_data_ref->{"\U$thorn\E FRIEND"}))
    {
      push (@friends, $other_implementation) if ($friend =~ m:$implementation:i);
    }
  }

  return join (' ', @friends);
}


sub get_implementation_friends
{
  my($implementation, $interface_data_ref, $friends_ref) = @_;

  $interface_data_ref->{"IMPLEMENTATION \U$implementation\E THORNS"} =~ m:(\w+):;

  my $thorn = $1;

  # Recurse
  foreach my $friend (split(" ", $interface_data_ref->{"\U$thorn\E FRIEND"}),
                      split(" ", $interface_data_ref->{"IMPLEMENTATION \U$implementation\E FRIENDS"}))
  {
    if(! $friends_ref->{"\U$friend\E"})
    {
      $friends_ref->{"\U$friend\E"} = 1;
      if(! $interface_data_ref->{"IMPLEMENTATION \U$friend\E THORNS"})
      {
        my $message = "$implementation is friends with $friend - non-existent implementation";
        &CST_error(0,$message,"",__LINE__,__FILE__);
        next;
      }
      &get_implementation_friends($friend, $interface_data_ref, $friends_ref);
    }
  }
}

sub get_implementation_ancestors
{
  my($implementation, $interface_data_ref, $ancestors_ref) = @_;

  $interface_data_ref->{"IMPLEMENTATION \U$implementation\E THORNS"} =~ m:(\w+):;

  my $thorn = $1;

  # Recurse.
  foreach my $ancestor (split(" ", $interface_data_ref->{"\U$thorn\E INHERITS"}))
  {
    if(! $ancestors_ref->{"\U$ancestor\E"})
    {
      $ancestors_ref->{"\U$ancestor\E"} = 1;
      if(! $interface_data_ref->{"IMPLEMENTATION \U$ancestor\E THORNS"})
      {
        # Implementation not found; give extensive information
        %info = &buildthorns("$cctk_home/arrangements","thorns");
        $suggest_thorns = "";
        foreach $thorninfo (sort keys %info)
        {
         $info{"$thorninfo"} =~ /^([^\s]+)/;
         $testimp = $1;
         if ($testimp =~ m:^$ancestor$:i)
         {
           $suggest_thorns .= "\n        $thorninfo";
         }
        }
        $message = "$implementation (thorn $thorn) inherits from $ancestor\n";
        $message .= "     No thorn in your current ThornList implements $ancestor\n";
        $message .= "     Either remove $thorn, or add a thorn to your\n";
        $message .= "      ThornList implementing $ancestor\n";
        if ($suggest_thorns !~ m:^$:)
        {
          $message .= "     Available thorns in arrangements directory implementing $ancestor:";
          $message .= "$suggest_thorns";
        }
        else
        {
          $message .= "     No thorns in arrangements directory implement $ancestor";
        }
        &CST_error(0,$message,"",__LINE__,__FILE__);

        next;
      }

      &get_implementation_ancestors($ancestor, $interface_data_ref, $ancestors_ref);
    }
  }
}

sub check_implementation_consistency
{
  my($implementation, $interface_data_ref) = @_;
  my(@thorns);
  my($thorn);
  my($thing);
  my(%inherits);
  my(%friend);
  my(%public_groups);
  my(%private_groups);
  my(%variables);
  my($n_errors);
  my($group);
  my(%attributes);

  # Find out which thorns provide this implementation.
  @thorns = split(" ", $interface_data_ref->{"IMPLEMENTATION \U$implementation\E THORNS"});

  if(scalar(@thorns) > 1)
  {
    foreach $thorn (@thorns)
    {
      # Record the inheritance
      foreach $thing (split(" ", $interface_data_ref->{"\U$thorn\E INHERITS"}))
      {
        if($thing =~ m:\w:)
        {
          # Put if statement around this to prevent perl -w from complaining.
          if($inherits{"\U$thing\E"})
          {
            $inherits{"\U$thing\E"} .= "$thorn ";
          }
          else
          {
            $inherits{"\U$thing\E"} = "$thorn ";
          }
        }
      }

      # Record the friends
      foreach $thing (split(" ", $interface_data_ref->{"\U$thorn\E FRIEND"}))
      {
        if($thing =~ m:\w:)
        {
          # Put if statement around this to prevent perl -w from complaining.
          if($friend{"\U$thing\E"})
          {
            $friend{"\U$thing\E"} .= "$thorn ";
          }
          else
          {
            $friend{"\U$thing\E"} = "$thorn ";
          }
        }
      }

      # Record the public groups
      foreach $thing (split(" ", $interface_data_ref->{"\U$thorn\E PUBLIC GROUPS"}))
      {
        if($thing =~ m:\w:)
        {
          # Put if statement around this to prevent perl -w from complaining.
          if($public_groups{"\U$thing\E"})
          {
            $public_groups{"\U$thing\E"} .= "$thorn ";
          }
          else
          {
            $public_groups{"\U$thing\E"} = "$thorn ";
          }
        }
      }

      # Record the protected groups
      foreach $thing (split(" ", $interface_data_ref->{"\U$thorn\E PROTECTED GROUPS"}))
      {
        if($thing =~ m:\w:)
        {
          # Put if statement around this to prevent perl -w from complaining.
          if($protected_groups{"\U$thing\E"})
          {
            $protected_groups{"\U$thing\E"} .= "$thorn ";
          }
          else
          {
            $protected_groups{"\U$thing\E"} = "$thorn ";
          }
        }
      }
    }

    $n_thorns = @thorns;

    # Check the consistency of the inheritance
    foreach $thing (sort keys %inherits)
    {
      if(split(' ', $inherits{$thing}) != $n_thorns)
      {
        &CST_error(0,
                   "Inconsistent implementation of '$implementation' " .
                   "provided by thorns '@thorns': not all inherit '$thing'",
                   '', __LINE__, __FILE__);
        $n_errors++;
      }
    }

    # Check the consistency of the friendships
    foreach $thing (sort keys %friend)
    {
      if(split(" ", $friend{$thing}) != $n_thorns)
      {
        $message  = "Inconsistent implementations of $implementation\n";
        $message .= "Implemented by thorns " . join(" ", @thorns) . "\n";
        $message .= "Not all are friends of: $thing";
        &CST_error(0,$message,"",__LINE__,__FILE__);
        $n_errors++;
      }
    }

    # Check the consistency of the public groups
    foreach $thing (sort keys %public_groups)
    {
      if(split(" ", $public_groups{$thing}) != $n_thorns)
      {
          $message  = "Inconsistent implementations of $implementation\n";
          $message .= "Implemented by thorns " . join(" ", @thorns) . "\n";
          $message .= "Not all declare public group: $thing";
          &CST_error(0,$message,"",__LINE__,__FILE__);
          $n_errors++;
      }
    }

    # Check the consistency of the protected groups
    foreach $thing (sort keys %protected_groups)
    {
      if(split(" ", $protected_groups{$thing}) != $n_thorns)
      {
        $message  = "Inconsistent implementations of $implementation\n";
        $message .= "Implemented by thorns " . join(" ", @thorns) . "\n";
        $message .= "Not all declare protected group: $thing";
        &CST_error(0,$message,"",__LINE__,__FILE__);
        $n_errors++;
      }
    }

    # Check consistancy of group definitions
    foreach $group ((sort keys %public_groups), (sort keys %protected_groups))
    {
      %variables = ();
      %attributes = ();

      foreach $thorn (@thorns)
      {
        # Remember which variables are defined in this group.
        foreach $thing (split(" ",$interface_data_ref->{"\U$thorn GROUP $group\E"}))
        {
          # Put if statement around this to prevent perl -w from complaining.
          if($variables{"\U$thing\E"})
          {
            $variables{"\U$thing\E"} .= "$thorn ";
          }
          else
          {
            $variables{"\U$thing\E"} = "$thorn ";
          }
        }

        # Check variable type definition.
        if($attributes{"VTYPE"})
        {
          if($attributes{"VTYPE"} ne $interface_data_ref->{"\U$thorn GROUP $group\E VTYPE"})
          {
            $message  = "Inconsistent implementations of $implementation";
            $message .= " in thorns " . join(" ", @thorns) . ". ";
            $message .= "Group $group has inconsistent variable type ($attributes{\"VTYPE\"} and $interface_data_ref->{\"\\U$thorn GROUP $group\\E VTYPE\"}). ";
            $hint = "All public and protected groups implementing $implementation must have groups with consistent properties";
            &CST_error(0,$message,$hint,__LINE__,__FILE__);
            $n_errors++;
          }
        }
        else
        {
          $attributes{"VTYPE"} = $interface_data_ref->{"\U$thorn GROUP $group\E VTYPE"};
        }

        # Check group type definition.
        if($attributes{"GTYPE"})
        {
          if($attributes{"GTYPE"} ne $interface_data_ref->{"\U$thorn GROUP $group\E GTYPE"})
          {
            $message  = "Inconsistent implementations of $implementation";
            $message .= " in thorns " . join(" ", @thorns) . ". ";
            $message .= "Group $group has inconsistent group type ($attributes{\"GTYPE\"} and $interface_data_ref->{\"\U$thorn GROUP $group\E GTYPE\"}). ";
            $hint = "All public and protected groups implementing $implementation must have groups with consistent properties";
            &CST_error(0,$message,$hint,__LINE__,__FILE__);
            $n_errors++;
          }
        }
        else
        {
          $attributes{"GTYPE"} = $interface_data_ref->{"\U$thorn GROUP $group\E GTYPE"};
        }

        # Check the number of time levels is consistent.
        if($attributes{"TIMELEVELS"})
        {
          if($attributes{"TIMELEVELS"} ne $interface_data_ref->{"\U$thorn GROUP $group\E TIMELEVELS"})
          {
            $message  = "Inconsistent implementations of $implementation\n";
            $message .= "Implemented by thorns " . join(" ", @thorns) . "\n";
            $message .= "Group $group has inconsistent time levels";
            &CST_error(0,$message,"",__LINE__,__FILE__);
            $n_errors++;
          }
        }
        else
        {
          $attributes{"TIMELEVELS"} = $interface_data_ref->{"\U$thorn GROUP $group\E TIMELEVELS"};
        }

        # Check the size array sizes are consistent.
        if($attributes{"SIZE"})
        {
          if($attributes{"SIZE"} ne $interface_data_ref->{"\U$thorn GROUP $group\E SIZE"})
          {
            $message  = "Inconsistent implementations of $implementation\n";
            $message .= "Implemented by thorns " . join(" ", @thorns) . "\n";
            $message .= "Group $group has inconsistent size";
            &CST_error(0,$message,"",__LINE__,__FILE__);
            $n_errors++;
          }
        }
        else
        {
          $attributes{"SIZE"} = $interface_data_ref->{"\U$thorn GROUP $group\E SIZE"};
        }

        # Check the ghostsize array sizes are consistent.
        if($attributes{"GHOSTSIZE"})
        {
          if($attributes{"GHOSTSIZE"} ne $interface_data_ref->{"\U$thorn GROUP $group\E GHOSTSIZE"})
          {
            $message  = "Inconsistent implementations of $implementation\n";
            $message .= "Implemented by thorns " . join(" ", @thorns) . "\n";
            $message .= "Group $group has inconsistent ghostsize";
            &CST_error(0,$message,"",__LINE__,__FILE__);
            $n_errors++;
          }
        }
        else
        {
          $attributes{"GHOSTSIZE"} = $interface_data_ref->{"\U$thorn GROUP $group\E GHOSTSIZE"};
        }

        # Check the distribution of arrays are consistent.
        if($attributes{"DISTRIB"})
        {
          if($attributes{"DISTRIB"} ne $interface_data_ref->{"\U$thorn GROUP $group\E DISTRIB"})
          {
            $message  = "Inconsistent implementations of $implementation\n";
            $message .= "Implemented by thorns " . join(" ", @thorns) . "\n";
            $message .= "      Group $group has inconsistent distribution";
            &CST_error(0,$message,"",__LINE__,__FILE__);
            $n_errors++;
          }
        }
        else
        {
          $attributes{"DISTRIB"} = $interface_data_ref->{"\U$thorn GROUP $group\E DISTRIB"};
        }

        # Check the dimensions are consistant
        if($attributes{"DIM"} && $attributes{"GTYPE"} ne "SCALAR")
        {
          if($attributes{"DIM"} ne $interface_data_ref->{"\U$thorn GROUP $group\E DIM"})
          {
            $message  = "Inconsistent implementations of $implementation\n";
            $message .= "Implemented by thorns " . join(" ", @thorns) . "\n";
            $message .= "Group $group has inconsistent dimension";
            &CST_error(0,$message,"",__LINE__,__FILE__);
            $n_errors++;
          }
        }
        else
        {
          $attributes{"DIM"} = $interface_data_ref->{"\U$thorn GROUP $group\E DIM"};
        }
      }
    }
  }
  else
  {
    # No need to do a consistency check if only one thorn
    # provides this implementation.

  }

}


#/*@@
#  @routine    check_interface_consistency
#  @date       Sun Jun 3 2001
#  @author     Gabrielle Allen
#  @desc
#  Check consistency of the interfaces files
#  @enddesc
#@@*/

sub check_interface_consistency
{
  my($thorn, $interface_data_ref) = @_;
  my($implementation);
  my($group,$var1,$var2,$group1,$group2);
  my($ancestor_imp,$ancestor_thorn,$ancestor2_imp,$ancestor2);
  my($message);

  # Find implementation
  $implementation =  $interface_data_ref->{"\U$thorn\E IMPLEMENTS"};

  # Loop over ancestors
  foreach $ancestor_imp (split " ",$interface_data_ref->{"IMPLEMENTATION \U$implementation\E ANCESTORS"})
  {
    # Need one thorn which implements this ancestor (we already have checked consistency)
    $ancestor_thorn = $interface_data_ref->{"IMPLEMENTATION \U$ancestor_imp\E THORNS"};
    if ($ancestor_thorn =~ m:(\w+)[^\w]*:)
    {
      $ancestor_thorn = $1;
    }
    foreach $group1 (split ' ', $interface_data_ref->{"\U$ancestor_thorn\E PUBLIC GROUPS"})
    {
      foreach $var1 (split ' ', $interface_data_ref->{"\U$ancestor_thorn\E GROUP \U$group1\E"})
      {
        foreach $ancestor2_imp (split " ",$interface_data_ref->{"IMPLEMENTATION \U$implementation\E ANCESTORS"})
        {
          $ancestor2 = $interface_data_ref->{"IMPLEMENTATION \U$ancestor2_imp\E THORNS"};
          if ($ancestor2 =~ m:(\w+)[^\w]*:)
          {
            $ancestor2 = $1;
          }
          # skip the second ancestor if it is the first one
          next if (uc($ancestor2) eq uc($ancestor_thorn));

          foreach $group2 (split ' ', $interface_data_ref->{"\U$ancestor2\E PUBLIC GROUPS"})
          {
            if (uc($group1) eq uc($group2))
            {
              $message = "Group $group1 from ancestor implementation $ancestor_imp in thorn $thorn has same name as \n     a public group: $group2 in ancestor implementation $ancestor2_imp (e.g. thorn $ancestor2)";
              &CST_error(1,$message,"",__LINE__,__FILE__);
            }
            if (uc($var1) eq uc($group2))
            {
              $message = "Variable $var1 in group $group1 from ancestor implementation $ancestor_imp in thorn $thorn has same name as \n     a public group: $group2 in ancestor implementation $ancestor2_imp (e.g. thorn $ancestor2)";
              &CST_error(1,$message,"",__LINE__,__FILE__);
            }
            foreach $var2 (split ' ', $interface_data_ref->{"\U$ancestor2\E GROUP \U$group2\E"})
            {
              if (uc($var2) eq uc($var1))
              {
                $message = "Variable $var1 in group $group1 from ancestor $ancestor_imp in thorn $thorn has same name as \n     variable $var2 in public group: $group2 in ancestor implementation $ancestor2_imp (e.g. thorn $ancestor2)";
                &CST_error(0,$message,"",__LINE__,__FILE__);
              }
            }
          }
        }
      }
    }

    foreach $group (split " ",$interface_data_ref->{"\U$thorn\E PRIVATE GROUPS"} . ' '. $interface_data_ref->{"\U$thorn\E PUBLIC GROUPS"} )
    {
      if ($interface_data_ref->{"\U$ancestor_thorn\E PUBLIC GROUPS"} =~ m:(\b$group\b):)
      {
        $message = "Group $group in thorn $thorn has same name as \n     public group in ancestor implementation $ancestor_imp (e.g. thorn $ancestor_thorn)";
        &CST_error(0,$message,"",__LINE__,__FILE__);
      }
      foreach $var (split " ", $interface_data_ref->{"\U$thorn\E GROUP \U$group\E"})
      {
        foreach $pub_anc(split " ", $interface_data_ref->{"\U$ancestor_thorn\E PUBLIC GROUPS"})
        {
          if ($interface_data_ref->{"\U$ancestor_thorn\E GROUP \U$pub_anc\E"} =~  m/\b$var\b/i)
          {
            $message = "Variable $var in group $group in thorn $thorn has same name as \n     a variable in public group: $pub_anc in ancestor implementation $ancestor_imp (e.g. thorn $ancestor_thorn)";
            &CST_error(0,$message,"",__LINE__,__FILE__);
          }

        }


      }

    }
  }
}



#/*@@
#  @routine    parse_interface_ccl
#  @date       Wed Sep 16 15:07:11 1998
#  @author     Tom Goodale
#  @desc
#  Parses an interface.ccl file and generates a database of the values.
#  @enddesc
#@@*/

sub parse_interface_ccl
{
  my($arrangement, $thorn, $data_ref, $interface_data_ref) = @_;
  my($line_number, $line, $block, $type, $variable, $description);
  my($data);
  my($implementation);
  my($option,%options);
  my(%known_groups);
  my(%known_variables);


  # Initialise some stuff to prevent perl -w from complaining.

  $interface_data_ref->{"\U$thorn INHERITS\E"} = "";
  $interface_data_ref->{"\U$thorn FRIEND\E"} = "";
  $interface_data_ref->{"\U$thorn PUBLIC GROUPS\E"} = "";
  $interface_data_ref->{"\U$thorn PROTECTED GROUPS\E"} = "";
  $interface_data_ref->{"\U$thorn PRIVATE GROUPS\E"} = "";
  $interface_data_ref->{"\U$thorn USES HEADER\E"} = "";
  $interface_data_ref->{"\U$thorn FUNCTIONS\E"} = "";
  $interface_data_ref->{"\U$thorn PROVIDES FUNCTION\E"} = " ";
  $interface_data_ref->{"\U$thorn REQUIRES FUNCTION\E"} = " ";
  $interface_data_ref->{"\U$thorn USES FUNCTION\E"} = " ";
  $interface_data_ref->{"\U$thorn ARRANGEMENT\E"} = "$arrangement";

  #   The default block is private.
  $block = "PRIVATE";

  for($line_number = 0; $line_number < @$data_ref; $line_number++)
  {
    $line = $data_ref->[$line_number];

    #       Parse the line
    if($line =~ m/^\s*(PUBLIC|PROTECTED|PRIVATE)\s*$/i)
    {
      #           It's a new block.
      $block = "\U$1\E";
    }
    elsif ($line =~ m/^\s*IMPLEMENTS\s*:/i)
    {
      if ($line =~ m/^\s*IMPLEMENTS\s*:\s*([a-z]+[a-z_0-9]*)\s*$/i)
      {
        if(!$implementation)
        {
          $implementation = $1;
          $interface_data_ref->{"\U$thorn\E IMPLEMENTS"} = $implementation;
        }
        else
        {
          $message = "Multiple implementations specified in $thorn";
          $hint = "A thorn can only specify one implementation in its interface.ccl file, with the format implements:<implementation>";
          &CST_error(0,$message,$hint,__LINE__,__FILE__);
        }
      }
      else
      {
        $message = "Implementation line has wrong format in $thorn";
        $hint = "A thorn must specify one implementation in its interface.ccl file with the format IMPLEMENTS: <implementation>";
        &CST_error(0,$message,$hint,__LINE__,__FILE__);
      }
    }
    # implementation names can be separated by ,\s, where , are stripped out below
    elsif ($line =~ m/^\s*(INHERITS|FRIEND)\s*:(([,\s]*[a-zA-Z]+[a-zA-Z_0-9]*)*[,\s]*)$/i)
    {
      $interface_data_ref->{"\U$thorn $1\E"} .= $2;
      $interface_data_ref->{"\U$thorn $1\E"}=~s/,/ /g;
    }
    elsif ($line =~ m/^\s*(PUBLIC|PROTECTED|PRIVATE)\s*:\s*$/i)
    {
      $block = "\U$1\E";
    }
    elsif ($line =~ m/^\s*PROVIDES\s*FUNCTION\s*([a-zA-Z_0-9]+)\s*WITH\s*(.+)\s*$/i)
    {
      $funcname = $1;
      $provided_by = $2;

      if($provided_by =~ m/^(.*)\s+LANGUAGE\s+(.*\S)\s*$/i)
      {
        $provided_by          = $1;
        $provided_by_language = "\U$2";
        if ($provided_by_language eq 'FORTRAN')
        {
          $provided_by_language = 'Fortran';
        }
        elsif ($provided_by_language ne 'C')
        {
          my $message = "The providing function $provided_by in thorn $thorn " .
                        "has an invalid language specification.";
          my $hint = "Language must be either C or Fortran.";
          &CST_error(0, $message, $hint, __LINE__, __FILE__);
        }
      }
      else
      {
#        $provided_by_language = "Fortran";
#        $provided_by_language = "C";
        $message = "The providing function $provided_by in thorn $thorn does not have a specified language. Please add, e.g., \"LANGUAGE C\"";
        &CST_error(0,$message,"",__LINE__,__FILE__);

      }

      if($funcname eq $provided_by) {
        my $message = "The providing function $provided_by in thorn $thorn " .
                      "has a name that is identical to the name of the provided " .
                      "function $funcname. The names must be different.";
        my $hint = "Rename the providing function by prefixing its name with ".
                   "'${thorn}_'.";
        &CST_error(0, $message, $hint, __LINE__, __FILE__);
      }

      $interface_data_ref->{"\U$thorn PROVIDES FUNCTION\E"} .= "$funcname ";
      $interface_data_ref->{"\U$thorn PROVIDES FUNCTION\E $funcname WITH"} .= "$provided_by ";
      $interface_data_ref->{"\U$thorn PROVIDES FUNCTION\E $funcname LANG"} .= "$provided_by_language ";

    }
    elsif ($line =~ m/^\s*REQUIRES\s*FUNCTION\s*([a-zA-Z_0-9]+)\s*$/i)
    {
      $funcname = $1;
      $interface_data_ref->{"\U$thorn REQUIRES FUNCTION\E"} .= "$funcname ";
    }
    elsif ($line =~ m/^\s*USES\s*FUNCTION\s*([a-zA-Z_0-9]+)\s*$/i)
    {
      $funcname = $1;
      $interface_data_ref->{"\U$thorn USES FUNCTION\E"} .= "$funcname ";
    }
    elsif ($line =~ m/^\s*([a-zA-Z][a-zA-Z_0-9:]+)\s*FUNCTION\s*([a-zA-Z_0-9]+)\s*\((.*)\)\s*$/i)
    {
      $rettype  = $1;
      $funcname = $2;
      $rest     = $3;

      $funcargs = $rest;

      $interface_data_ref->{"\U$thorn FUNCTIONS\E"} .= "${funcname} ";
      $interface_data_ref->{"\U$thorn FUNCTION\E $funcname ARGS"} .= "${funcargs} ";
      $interface_data_ref->{"\U$thorn FUNCTION\E $funcname RET"} .= "${rettype} ";
    }
    elsif ($line =~ m/^\s*SUBROUTINE\s*([a-zA-Z_0-9]+)\s*\((.*)\)\s*$/i)
    {
      $rettype  = "void";
      $funcname = $1;
      $rest     = $2;

      $funcargs = $rest;

      $interface_data_ref->{"\U$thorn FUNCTIONS\E"} .= "${funcname} ";
      $interface_data_ref->{"\U$thorn FUNCTION\E $funcname ARGS"} .= "${funcargs} ";
      $interface_data_ref->{"\U$thorn FUNCTION\E $funcname RET"} .= "${rettype} ";
    }
    elsif ($line =~ m/^\s*(CCTK_)?(CHAR|BYTE|INT|INT1|INT2|INT4|INT8|INT16|REAL|REAL4|REAL8|REAL16|COMPLEX|COMPLEX8|COMPLEX16|COMPLEX32)\s*(([a-zA-Z][a-zA-Z_0-9]*)\s*(\[([^]]+)\])?)\s*(.*)\s*$/i)
    {
#      for($i = 1; $i < 10; $i++)
#      {
#        print "$i is ${$i}\n";
#      }
      my $vtype = $2;
      my $current_group = "$4";
      my $isgrouparray = $5;
      my $grouparray_size = $6;
      my $options_list = $7;

#      print "line is [$line]\n";
#      print "group name is [$current_group]\n";
#      print "options list is [$options_list]\n";

      if($known_groups{"\U$current_group\E"})
      {
        &CST_error(0,"Duplicate group $current_group in thorn $thorn",'',
                   __LINE__,__FILE__);
        if($data_ref->[$line_number+1] =~ m:\{:)
        {
            &CST_error(1,'Skipping interface block','',__LINE__,__FILE__);
            $line_number++ until ($data_ref->[$line_number] =~ m:\}:);
        }
        next;
      }
      else
      {
        $known_groups{"\U$current_group\E"} = 1;

        # Initialise some stuff to prevent perl -w from complaining.
        $interface_data_ref->{"\U$thorn GROUP $current_group\E"} = "";
      }

      $interface_data_ref->{"\U$thorn $block GROUPS\E"} .= " $current_group";
      $interface_data_ref->{"\U$thorn GROUP $current_group\E VTYPE"} = "\U$vtype\E";

      # Grab optional group description from end of $options_list
      if ($options_list =~ /(=?)\s*"([^"]*)"\s*$/)
      {
        if (!$1)
        {
          if ($data_ref->[$line_number+1] =~ m/^\s*\{\s*$/)
          {
            &CST_error(1,"Group description for $current_group in thorn " .
                       "$thorn must be placed at end of variable block " .
                       "when variable block present",'',
                       __LINE__,__FILE__);
          } else
          {
            $description = $2;
            $quoted_description = quotemeta ($description);
            $options_list =~ s/\s*"$quoted_description"//;
          }
        }
      }

      # split(/\s*=\s*|\s+/, $options_list);
      %options = SplitWithStrings($options_list,$thorn);

      # Parse the options
      foreach $option (sort keys %options)
      {

#        print "DEBUG $option is $options{$option}\n";

        if($option =~ m:DIM|DIMENSION:i)
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E DIM"} = $options{$option};
        }
        elsif($option =~ m:TYPE:i)
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E GTYPE"} = "\U$options{$option}\E";
        }
        elsif($option =~ m:TIMELEVELS:i)
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E TIMELEVELS"} = "\U$options{$option}\E";
        }
        elsif($option =~ m:GHOSTSIZE:i)
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E GHOSTSIZE"} = "\U$options{$option}\E";
        }
        elsif($option =~ m:DISTRIB:i)
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"} = "\U$options{$option}\E";
        }
        elsif($option =~ m:SIZE:i)
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E SIZE"} = "\U$options{$option}\E";
        }
        elsif($option =~ m:TAGS:i)
        {
          if($options{$option} =~ m/\s*^[\'\"](.*)[\'\"]$/)
          {
            $options{$option} = $1;
          }

          $options{$option} =~ s/\\/\\\\/g;
          $options{$option} =~ s/\"/\\\"/g;

          $interface_data_ref->{"\U$thorn GROUP $current_group\E TAGS"} = $options{$option};
        }
        else
        {
          &CST_error(0,"Unknown option \"$option\" in group $current_group in interface.ccl for " .
                     "of thorn $thorn\n     Perhaps you forgot a '\\' at the " .
                     "end of a continued line?\n" .
                     "The offending line is '$line'\n",'',
                     __LINE__,__FILE__);
        }
      }

      # Put in defaults
      if(! $interface_data_ref->{"\U$thorn GROUP $current_group\E GTYPE"})
      {
        $interface_data_ref->{"\U$thorn GROUP $current_group\E GTYPE"} = "SCALAR";
      }

      if (! $interface_data_ref->{"\U$thorn GROUP $current_group\E DIM"})
      {
        if ($interface_data_ref->{"\U$thorn GROUP $current_group\E GTYPE"} eq 'SCALAR')
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E DIM"} = 0;
        }
        else
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E DIM"} = 3;
        }
      }

      if(! $interface_data_ref->{"\U$thorn GROUP $current_group\E TIMELEVELS"})
      {
        $interface_data_ref->{"\U$thorn GROUP $current_group\E TIMELEVELS"} = 1;
      }

      if(! $interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"})
      {
        if ($interface_data_ref->{"\U$thorn GROUP $current_group\E GTYPE"} eq 'SCALAR')
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"} = 'CONSTANT';
        }
        else
        {
          $interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"} = 'DEFAULT';
        }
      }

      if(! $interface_data_ref->{"\U$thorn GROUP $current_group\E COMPACT"})
      {
        $interface_data_ref->{"\U$thorn GROUP $current_group\E COMPACT"} = 0;
      }

      if ($interface_data_ref->{"\U$thorn GROUP $current_group\E GTYPE"} eq "SCALAR")
      {
        my $dim = $interface_data_ref->{"\U$thorn GROUP $current_group\E DIM"};
        if ($dim && $dim ne '0')
        {
          my $message =  "Inconsistent GROUP DIM $dim for SCALAR group $current_group of thorn $thorn";
          my $hint = "The only allowed group dimension for scalar groups is '0'";
          &CST_error (0, $message, $hint, __LINE__, __FILE__);
          if ($data_ref->[$line_number+1] =~ m:\{:)
          {
            &CST_error (1, "Skipping interface block in $thorn", '',
                        __LINE__, __FILE__);
            ++ $line_number until ($data_ref->[$line_number] =~ m:\}:);
          }
          next;
        }
        $interface_data_ref->{"\U$thorn GROUP $current_group\E DIM"} = 0;

        my $distrib = $interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"};
        if ($distrib && $distrib ne 'CONSTANT')
        {
          my $message =  "Inconsistent GROUP DISTRIB $distrib for SCALAR group $current_group of thorn $thorn";
          my $hint = "The only allowed group distribution for scalar groups is 'CONSTANT'";
          &CST_error (0, $message, $hint, __LINE__, __FILE__);
          if ($data_ref->[$line_number+1] =~ m:\{:)
          {
            &CST_error (1, "Skipping interface block in $thorn", '',
                        __LINE__, __FILE__);
            ++ $line_number until ($data_ref->[$line_number] =~ m:\}:);
          }
          next;
        }
        $interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"} =
            "CONSTANT";
      }

      # Override defaults for grid functions
      if ($interface_data_ref->{"\U$thorn GROUP $current_group\E GTYPE"} eq "GF")
      {
        my $distrib = $interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"};
        if ($distrib && $distrib ne 'DEFAULT')
        {
          my $message =  "Inconsistent GROUP DISTRIB $distrib for GF group $current_group of thorn $thorn";
          my $hint = "The only allowed group distribution for grid function groups is 'DEFAULT'";
          &CST_error (0, $message, $hint, __LINE__, __FILE__);
          if ($data_ref->[$line_number+1] =~ m:\{:)
          {
            &CST_error (1, "Skipping interface block in $thorn", '',
                        __LINE__, __FILE__);
            ++ $line_number until ($data_ref->[$line_number] =~ m:\}:);
          }
          next;
        }
        $interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"} =
            "DEFAULT";
      }

      # Check that it is a known group type
      if($interface_data_ref->{"\U$thorn GROUP $current_group\E GTYPE"} !~ m:^\s*(SCALAR|GF|ARRAY)\s*$:)
      {
          $message =  "Unknown GROUP TYPE " .
          $interface_data_ref->{"\U$thorn GROUP $current_group\E GTYPE"} .
            " for group $current_group of thorn $thorn";
          $hint = "Allowed group types are SCALAR, GF or ARRAY";
          &CST_error(0,$message,$hint,__LINE__,__FILE__);
          if($data_ref->[$line_number+1] =~ m:\{:)
          {
              &CST_error(1,"Skipping interface block in $thorn","",
                         __LINE__,__FILE__);
              $line_number++ until ($data_ref->[$line_number] =~ m:\}:);
          }
        next;
      }

      # Check that it is a known distribution type
      if($interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"} !~ m:DEFAULT|CONSTANT:)
      {
          $message =  "Unknown DISTRIB TYPE " .
          $interface_data_ref->{"\U$thorn GROUP $current_group\E DISTRIB"} .
            " for group $current_group of thorn $thorn";
          $hint = "Allowed distribution types are DEFAULT or CONSTANT";
          &CST_error(0,$message,"",__LINE__,__FILE__);
          if($data_ref->[$line_number+1] =~ m:\{:)
          {
              &CST_error(1,"Skipping interface block in $thorn",'',
                         __LINE__,__FILE__);
              $line_number++ until ($data_ref->[$line_number] =~ m:\}:);
          }
        next;
      }

      # Is it a vararray?
      if($isgrouparray)
      {
        # get its size
        $interface_data_ref->{"\U$thorn GROUP $current_group\E VARARRAY_SIZE"} = $grouparray_size;
      }
        # Fill in data for the scalars/arrays/functions
        $line_number++;
        if($data_ref->[$line_number] =~ m/^\s*\{\s*$/)
        {
          $line_number++;
          while($data_ref->[$line_number] !~ m:\}:i)
          {
            @functions = split(/[^a-zA-Z_0-9]+/, $data_ref->[$line_number]);
            foreach $function (@functions)
            {
              if ($function eq $current_group)
              {
                if ($#functions == 1)
                {
                  &CST_error(1,"Group and variable '$function' in thorn " .
                             "'$thorn' should be distinct",'',
                             __LINE__,__FILE__);
                }
                else
                {
                  &CST_error(0,"The names of all variables '@functions' must " .
                             "be different from their group name " .
                             "'$current_group' in thorn '$thorn'",'',
                             __LINE__,__FILE__);
                }
              }
              $function =~ s:\s*::g;

              if($function =~ m:[^\s]+:)
              {
                if(! $known_variables{"\U$function\E"})
                {
                  $known_variables{"\U$function\E"} = 1;

                  $interface_data_ref->{"\U$thorn GROUP $current_group\E"} .= " $function";
                }
                else
                {
                  &CST_error(0,"Duplicate variable $function in thorn $thorn",'',
                             __LINE__,__FILE__);
                }
              }
            }
            $line_number++;
          }
          # Grab optional group description
          $data_ref->[$line_number] =~ m:\}\s*"([^"]*)":;
          $description = $1;
        }
        else
        {
          # If no block, create a variable with the same name as group.
          $function = $current_group;
          if(! $known_variables{"\U$function\E"})
          {
            $known_variables{"\U$function\E"} = 1;

            $interface_data_ref->{"\U$thorn GROUP $current_group\E"} .= " $function";
          }
          else
          {
            &CST_error(0,"Duplicate variable $function in thorn $thorn",'',
                       __LINE__,__FILE__);
          }

          # Decrement the line number, since the line is the first line of the next CCL statement.
          $line_number--;
        }
        $interface_data_ref->{"\U$thorn GROUP $current_group\E DESCRIPTION"} = $description;
    }
    elsif ($line =~ m/^\s*(USES\s*INCLUDE)S?\s*(SOURCE)S?\s*:\s*(.*)\s*$/i)
    {
      $interface_data_ref->{"\U$thorn USES SOURCE\E"} .= " $3";
    }
    elsif ($line =~ m/^\s*(USES\s*INCLUDE)S?\s*(HEADER)?S?\s*:\s*(.*)\s*$/i)
    {
      $interface_data_ref->{"\U$thorn USES HEADER\E"} .= " $3";
    }
    elsif ($line =~ m/^\s*(INCLUDE)S?\s*(SOURCE)S?\s*:\s*(.*)\s+IN\s+(.*)\s*$/i)
    {
      $header = $3;
      $header =~ s/ //g;
      $interface_data_ref->{"\U$thorn ADD SOURCE\E"} .= " $header";
#      print "Adding $header to $4\n";
      $interface_data_ref->{"\U$thorn ADD SOURCE $header TO\E"} = $4;
    }
    elsif ($line =~ m/^\s*(INCLUDE)S?\s*(HEADER)?S?\s*:\s*(\S*)\s+IN\s+(\S*)\s*$/i)
    {
      $header = $3;
      $header =~ s/ //g;
      $interface_data_ref->{"\U$thorn ADD HEADER\E"} .= " $header";
#      print "Adding $header to $4\n";
      $interface_data_ref->{"\U$thorn ADD HEADER $header TO\E"} = $4;
    }
    else
    {
      if($line =~ m:\{:)
      {
        &CST_error(0,'...Skipping interface block with missing keyword....','',
                   __LINE__,__FILE__);

        $line_number++ until ($data_ref->[$line_number] =~ m:\}:);
      }
      else
      {
        &CST_error(0,"Unknown line in interface.ccl for thorn $arrangement/$thorn\n\"$line\"",'',__LINE__,__FILE__) if ($line);
      }
    }
  }
}

#/*@@
#  @routine    PrintInterfaceStatistics
#  @date       Sun Sep 19 13:03:23 1999
#  @author     Tom Goodale
#  @desc
#  Prints out some statistics about a thorn's interface.ccl
#  @enddesc
#@@*/
sub PrintInterfaceStatistics
{
  my($thorn, $interface_database_ref) = @_;
  my($block);
  my($sep);

  print "           Implements: " . $interface_database_ref->{"\U$thorn IMPLEMENTS"} . "\n";

  if($interface_database_ref->{"\U$thorn INHERITS"} ne "")
  {
    print "           Inherits:  " . $interface_database_ref->{"\U$thorn INHERITS"} . "\n";
  }

  if($interface_database_ref->{"\U$thorn FRIEND"} ne "")
  {
    print "           Friend of: " . $interface_database_ref->{"\U$thorn FRIEND"} . "\n";
  }

  $sep = "           ";
  foreach $block ("Public", "Protected", "Private")
  {
    print $sep . scalar(split(" ", $interface_database_ref->{"\U$thorn $block\E GROUPS"})) . " $block";
    $sep = ", ";
  }

  print " variable groups\n";

  return;
}

1;
