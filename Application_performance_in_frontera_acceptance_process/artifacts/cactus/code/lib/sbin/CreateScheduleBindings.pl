#/*@@
#  @file      CreateScheduleBindings.pl
#  @date      Thu Sep 16 23:30:21 1999
#  @author    Tom Goodale
#  @desc
#             New schedule stuff.  Should be renamed !!!
#  @enddesc
#  @version   $Header$
#@@*/

#/*@@
#  @routine    CreateScheduleBindings
#  @date       Fri Sep 17 14:16:23 1999
#  @author     Tom Goodale
#  @desc
#              Creates the schedule bindings.
#  @enddesc
#@@*/
sub CreateScheduleBindings
{
  my($bindings_dir, $rhinterface_db, $rhschedule_db) = @_;
  my($start_dir);
  my($thorn);
  my($file_list);
  my($rsbuffer);

  if(! -d $bindings_dir)
  {
    mkdir($bindings_dir, 0755) || die "Unable to create $bindings_dir";
  }
  $start_dir = `pwd`;

  chdir $bindings_dir;

  if(! -d 'Schedule')
  {
    mkdir('Schedule', 0755) || die 'Unable to create Schedule directory';
  }

  if(! -d 'include')
  {
    mkdir('include', 0755) || die 'Unable to create include directory';
  }

  chdir 'Schedule';

  $file_list = '';

  foreach $thorn (sort split(" ", $rhinterface_db->{"THORNS"}))
  {
    $rsbuffer = &ScheduleCreateFile($thorn, $rhinterface_db, $rhschedule_db);
    &WriteFile("Schedule$thorn.c",\$rsbuffer);
    $file_list .= " Schedule$thorn.c";

    $rsbuffer = &ScheduleCreateInterfaceFile($thorn, $rhinterface_db, $rhschedule_db);
    &WriteFile("../include/${thorn}/cctk_ScheduleFunctions.h",\$rsbuffer);
    if($thorn eq "Cactus") {
        &WriteFile("../include/CactusBindings/cctk_ScheduleFunctions.h",\$rsbuffer);
    }
  }

  $rsbuffer = &ScheduleCreateBindings($rhinterface_db, $rhschedule_db);
  &WriteFile("BindingsSchedule.c", \$rsbuffer);
  $file_list .= " BindingsSchedule.c";

  $rsbuffer = &ParameterRecoveryCreateBindings($rhinterface_db, $rhschedule_db);
  &WriteFile("BindingsParameterRecovery.c", \$rsbuffer);
  $file_list .= " BindingsParameterRecovery.c";

  $line = "SRCS = $file_list\n";
  &WriteFile("make.code.defn", \$line);

  chdir "$start_dir";
}

#/*@@
#  @routine    ScheduleCreateInterfaceFile
#  @date       2006-05-11
#  @author     Erik Schnetter
#  @desc
#              Create a string containing all the data which should go
#              into a schedule interface file.
#  @enddesc
#@@*/
sub ScheduleCreateInterfaceFile
{
  my ($thorn, $rhinterface_db, $rhschedule_db) = @_;
  
  # Map groups to one of their block numbers (groups may have several
  # block numbers)
  my %group_block = ();
  
  for (my $block = 0;
       $block < $rhschedule_db->{"\U$thorn\E N_BLOCKS"};
       ++ $block)
  {
      if ($rhschedule_db->{"\U$thorn\E BLOCK_$block TYPE"} eq "GROUP")
      {
          my $group = $rhschedule_db->{"\U$thorn\E BLOCK_$block NAME"};
          $group_block{$group} = $block;
      }
  }
  
  
  
  # Process each schedule block
  
  my @data = ();
  push (@data, "#include \"$thorn/cctk_Arguments.h\"");
  
  for (my $block = 0;
       $block < $rhschedule_db->{"\U$thorn\E N_BLOCKS"};
       ++ $block)
  {
      if ($rhschedule_db->{"\U$thorn\E BLOCK_$block TYPE"} eq "FUNCTION")
      {
          my ($language, $function);
          if ($rhschedule_db->{"\U$thorn\E BLOCK_$block LANG"} =~ m:^\s*C\s*$:i)
          {
              $language = 'C';
              $function = $rhschedule_db->{"\U$thorn\E BLOCK_$block NAME"};
          }
          elsif ($rhschedule_db->{"\U$thorn\E BLOCK_$block LANG"} =~ m:^\s*(F|F77|FORTRAN|F90)\s*$:i)
          {
              $language = 'Fortran';
              $function = "CCTK_FNAME(".$rhschedule_db->{"\U$thorn\E BLOCK_$block NAME"} .")";
          }
          
          my $where = $rhschedule_db->{"\U$thorn\E BLOCK_$block WHERE"};
          # Find one outermost enclosing group iteratively (there may
          # be serveral outermost enclosing groups)
          my %been_there;       # avoid cycles
          while (! $been_there{$where} && defined $group_block{$where})
          {
              $been_there{$where} = defined;
              my $block1 = $group_block{$where};
              $where = $rhschedule_db->{"\U$thorn\E BLOCK_$block1 WHERE"};
          }
          my $is_special =
              $where eq 'CCTK_STARTUP' ||
              $where eq 'CCTK_RECOVER_PARAMETERS' ||
              $where eq 'CCTK_SHUTDOWN';
          
          push (@data, '');
          if ($language eq 'C')
          {
              push (@data, '#ifdef CCODE');
              push (@data, '#ifdef __cplusplus');
              push (@data, 'extern "C"');
              push (@data, '#endif');
              if ($is_special)
              {
                  push (@data, "int $function (void);")
              }
              else
              {
                  push (@data, "void $function (CCTK_ARGUMENTS) CCTK_ATTRIBUTE_NONNULL(1);")
              }
              push (@data, '#endif /* CCODE */');
          }
          elsif ($language eq 'Fortran')
          {
# Can't use CCTK_FNAME in header files (which is implicit in $function)
#              push (@data, '#ifdef CCODE');
#              push (@data, '#ifdef __cplusplus');
#              push (@data, 'extern "C"');
#              push (@data, '#endif');
#              if ($is_special)
#              {
#                  push (@data, "CCTK_FCALL int $function (void);")
#              }
#              else
#              {
#                  push (@data, "CCTK_FCALL void $function (CCTK_FARGUMENTS);")
#              }
#              push (@data, '#endif /* CCODE */');
# Can't declare functions in Fortran at file scope
#              push (@data, '#ifdef FCODE');
#              push (@data, '#ifdef F90CODE');
#              push (@data, '      interface');
#              if ($is_special)
#              {
#                  push (@data, "      integer function $function ()");
#                  push (@data, '         implicit none');
#                  push (@data, "      end function $function");
#              }
#              else
#              {
#                  push (@data, "      subroutine $function (CCTK_ARGUMENTS)");
#                  push (@data, '         implicit none');
#                  push (@data, '         DECLARE_CCTK_ARGUMENTS');
#                  push (@data, "      end subroutine $function");
#              }
#              push (@data, '      end interface');
#              push (@data, '#else /* ! F90CODE */');
#              push (@data, "      external $function");
#              if ($is_special)
#              {
#                  push (@data, "      integer $function");
#              }
#              push (@data, '#endif /* F90CODE */');
#              push (@data, '#endif /* FCODE */');
          }
          
      }
  }
  
  push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline
  
  return join ("\n", @data);
}



#/*@@
#  @routine    ScheduleCreateFile
#  @date       Fri Sep 17 17:34:26 1999
#  @author     Tom Goodale
#  @desc
#              Creates a string containing all the data
#              which should go into a schedule file.
#  @enddesc
#@@*/
sub ScheduleCreateFile
{
  my($thorn, $rhinterface_db, $rhschedule_db) = @_;

  my($implementation);
  my($buffer, $prototypes);
  my($recovery_buffer, $recovery_prototypes);
  my($block, $block_buffer, $block_prototype);
  my($statement, $statement_buffer, $statement_prototype);
  my(@data);

  $implementation = $rhinterface_db->{"\U$thorn\E IMPLEMENTS"};

  $buffer = $recovery_buffer = $rhschedule_db->{"\U$thorn\E FILE"};

  # Process each schedule block
  for($block = 0 ; $block < $rhschedule_db->{"\U$thorn\E N_BLOCKS"}; $block++)
  {
    if ($rhschedule_db->{"\U$thorn\E BLOCK_$block WHERE"} !~ /RECOVER_PARAMETERS/)
    {
      ($block_buffer, $block_prototype) = &ScheduleBlock($thorn, $implementation, $block,
                                                         $rhinterface_db, $rhschedule_db);
      $buffer =~ s:\@BLOCK\@$block:$block_buffer:;
      $recovery_buffer =~ s:\@BLOCK\@$block::;
      $prototypes .= "$block_prototype";
    }
    else
    {
      $block_buffer = $rhschedule_db->{"\U$thorn\E BLOCK_$block NAME"};
      if($rhschedule_db->{"\U$thorn\E BLOCK_$block LANG"} !~ m:^\s*C\s*$:i )
      {
        $block_buffer = "CCTK_FNAME ($block_buffer)";
      }

      $recovery_buffer =~ s:\@BLOCK\@$block:result =  $block_buffer();:;
      $buffer =~ s:\@BLOCK\@$block::;
      $recovery_prototypes .= "extern int $block_buffer(void);\n";
    }
  }

  # Process each schedule statement
  for($statement = 0 ; $statement < $rhschedule_db->{"\U$thorn\E N_STATEMENTS"}; $statement++)
  {
    ($statement_buffer, $statement_prototype) = &ScheduleStatement($thorn, $implementation, $statement,
                                                                   $rhinterface_db, $rhschedule_db);
    $buffer =~ s:\@STATEMENT\@$statement:$statement_buffer:;
    $recovery_buffer =~ s:\@STATEMENT\@$statement::;

    $prototypes .= $statement_prototype;
  }

  # Actually create the file contents string
  @data = ();
  push(@data, '/*@@');
  push(@data, "   \@file       Schedule$thorn.c");
  push(@data, '   @author     Automatically generated by CreateScheduleBindings.pl');
  push(@data, '   @desc');
  push(@data, '               Creates the schedule and parameter recovery bindings ');
  push(@data, "               for thorn $thorn");
  push(@data, '   @enddesc');
  push(@data, '@@*/');
  push(@data, '');
  #push(@data, "#define THORN_IS_$thorn");
  #push(@data, "#define CCTK_THORNSTRING \"$thorn\"");
  push(@data, '');
  push(@data, "#include \"$thorn/cctk.h\"");
  push(@data, "#include \"$thorn/CParameters.h\"");
  push(@data, '#include "cctki_ScheduleBindings.h"');
  push(@data, "#include \"$thorn/cctk_ScheduleFunctions.h\"");
  push(@data, '');
  #push(@data, '/* prototypes for schedule bindings functions to be registered */');
  #push(@data, '/* Note that this is a cheat, we just need a function pointer. */');
  #push(@data, $prototypes);
  #push(@data, '');
  push(@data, '/* Prototypes for Fortran schedule bindings functions to be registered */');
  push(@data, '/* Note that this is a cheat, we just need a function pointer. */');
  push(@data, join "\n", (grep /CCTK_FNAME/, (split "\n", $prototypes)));
  push(@data, '');
  push(@data, "void CCTKi_BindingsSchedule_$thorn(void);");
  push(@data, "void CCTKi_BindingsSchedule_$thorn(void)");
  push(@data, '{');
  push(@data, '  DECLARE_CCTK_PARAMETERS');

  # filter out empty lines
  foreach $line (split ("\n", $buffer))
  {
    push(@data, $line) if ($line);
  }

  push(@data, '}');

  push(@data, '');
  push(@data, '/*@@');
  push(@data, "  \@routine    CCTKi_BindingsParameterRecovery_$thorn");
  push(@data, '  @author     Automatically generated by CreateScheduleBindings.pl');
  push(@data, '  @desc');
  push(@data, "              Creates the parameter recovery bindings for thorn $thorn");
  push(@data, '  @enddesc');
  push(@data, '@@*/');
  push(@data, '');
  if (defined ($recovery_prototypes))
  {
    push(@data, '/* prototypes of parameter recovery functions to be registered. */');
    push(@data, "$recovery_prototypes");
  }
  push(@data, "int CCTKi_BindingsParameterRecovery_$thorn(void);");
  push(@data, "int CCTKi_BindingsParameterRecovery_$thorn(void)");
  push(@data, '{');
  if (defined ($recovery_prototypes))
  {
    push(@data, '  DECLARE_CCTK_PARAMETERS');
    push(@data, '  int result = 0;');
    push(@data, '');

    # filter out empty lines
    foreach $line (split ("\n", $recovery_buffer))
    {
      push(@data, $line) if ($line);
    }

    push(@data, '');
    push(@data, '  return (result);');
  }
  else
  {
    push(@data, '  /* this thorn doesn\'t define any parameter recovery routines */');
    push(@data, '  return (0);');
  }
  push(@data, '}');
  push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline

  return join ("\n", @data);
}


#/*@@
#  @routine    ScheduleCreateBindings
#  @date       Fri Sep 17 18:17:13 1999
#  @author     Tom Goodale
#  @desc
#  Creates a string containing all the data which should go into the master
#  schedule bindings file.
#  @enddesc
#@@*/
sub ScheduleCreateBindings
{
  my($rhinterface_db, $rhschedule_db) = @_;
  my(@data);

  @data = ();
  push(@data, '/*@@');
  push(@data, '   @file    BindingsSchedule.c');
  push(@data, '   @author  Automatically generated by CreateScheduleBindings.pl');
  push(@data, '   @desc');
  push(@data, '            Calls all the thorn schedule bindings file if the thorns are active.');
  push(@data, '   @enddesc');
  push(@data, '  @@*/');
  push(@data, '');

  push(@data, '#include "cctk_ActiveThorns.h"');
  push(@data, '');

  push(@data, '/* Prototypes for functions to be registered. */');
  foreach $thorn (sort split(' ', $rhinterface_db->{"THORNS"}))
  {
    push(@data, "void CCTKi_BindingsSchedule_$thorn(void);");
  }

  push(@data, '');
  push(@data, 'int CCTKi_BindingsScheduleInitialise(void);');
  push(@data, 'int CCTKi_BindingsScheduleInitialise(void)');
  push(@data, '{');
  foreach $thorn (sort split(' ', $rhinterface_db->{"THORNS"}))
  {
    push(@data, "  if(CCTK_IsThornActive(\"$thorn\"))");
    push(@data, '  {');
    push(@data, "    CCTKi_BindingsSchedule_$thorn();");
    push(@data, '  }');
  }
  push(@data, '  return 0;');
  push(@data, '}');
  push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline

  return join ("\n", @data);
}


#/*@@
#  @routine    ParameterRecoveryCreateBindings
#  @date       Tue Apr 18 13:17:13 2000
#  @author     Gabrielle Allen
#  @desc
#  Creates a string containing all the data which should go into the master
#  parameter recovery bindings file.
#  @enddesc
#@@*/
sub ParameterRecoveryCreateBindings
{
  my($rhinterface_db, $rhschedule_db) = @_;
  my(@data);

  @data = ();
  push(@data, '/*@@');
  push(@data, '   @file       BindingsParameterRecovery.c');
  push(@data, '   @author     Automatically generated by CreateScheduleBindings.pl');
  push(@data, '   @desc');
  push(@data, '               Calls all the thorn parameter recovery bindings file if the thorns are active.');
  push(@data, '   @enddesc');
  push(@data, '  @@*/');
  push(@data, '');

  push(@data, '#include "cctk_ActiveThorns.h"');
  push(@data, '');
  push(@data, '/* Prototypes for functions to be registered. */');

  foreach $thorn (sort split(' ', $rhinterface_db->{"THORNS"}))
  {
    push(@data, "int CCTKi_BindingsParameterRecovery_$thorn(void);");
  }

  push(@data, 'int CCTKi_BindingsParameterRecoveryInitialise(void);');
  push(@data, '');
  push(@data, 'int CCTKi_BindingsParameterRecoveryInitialise(void)');
  push(@data, '{');
  push(@data, '  int result;');
  push(@data, '  int retval = 0;');
  push(@data, '');
  push(@data, '  do');
  push(@data, '  {');
  foreach $thorn (sort split(' ', $rhinterface_db->{"THORNS"}))
  {
    push(@data, "    if(CCTK_IsThornActive(\"$thorn\"))");
    push(@data, '    {');
    push(@data, "      result = CCTKi_BindingsParameterRecovery_$thorn();");
    push(@data, '      if (result != 0)');
    push(@data, '        retval = result;');
    push(@data, '      if (retval > 0)');
    push(@data, '        break;');
    push(@data, '    }');
  }
  push(@data, '  } while (0);');
  push(@data, '');
  push(@data, '  return retval;');
  push(@data, '}');
  push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline

  return join ("\n", @data);
}


#/*@@
#  @routine    ScheduleBlock
#  @date       Fri Sep 17 17:37:59 1999
#  @author     Tom Goodale
#  @desc
#  Creates the code for a given schedule block
#  @enddesc
#@@*/
sub ScheduleBlock
{
  my($thorn, $implementation, $block, $rhinterface_db, $rhschedule_db) = @_;

  my($buffer, $prototype);
  my($indent, $language, $function);
  my($mem_groups);
  my($tlist);
  my($comm_groups);
  my($unused_comm_groups);
  my($trigger_groups);
  my($sync_groups);
  my(@options);
  my($tags);
  my(@before_list);
  my(@after_list);
  my(@writes_list);
  my(@reads_list);
  my(@while_list);
  my(@if_list);

  # Extract group and routine information from the databases
  ($mem_groups, $tlist) = &ScheduleSelectGroups($thorn, $implementation,
                                                $rhschedule_db->{"\U$thorn\E BLOCK_$block STOR"},
                                                $rhinterface_db);

  &ScheduleValidateTimeLevels($thorn, $implementation, $mem_groups,$tlist, $rhinterface_db);


  ($unused_comm_groups) = &ScheduleSelectGroups($thorn, $implementation,
                                                $rhschedule_db->{"\U$thorn\E BLOCK_$block COMM"},
                                                $rhinterface_db);
  if (@$unused_comm_groups)
  {
    print "No need to switch on Communication in $thorn\n";
    print "Communication automatically assigned for variables with storage\n";
  }

  $comm_groups = [];

  ($trigger_groups) = &ScheduleSelectGroups($thorn, $implementation,
                                            $rhschedule_db->{"\U$thorn\E BLOCK_$block TRIG"},
                                            $rhinterface_db);

  ($sync_groups) = &ScheduleSelectGroups($thorn, $implementation,
                                         $rhschedule_db->{"\U$thorn\E BLOCK_$block SYNC"},
                                         $rhinterface_db);

  @options = split(/,/, $rhschedule_db->{"\U$thorn\E BLOCK_$block OPTIONS"});
  $tags = $rhschedule_db->{"\U$thorn\E BLOCK_$block TAGS"};
  $tags = '' if ! defined $tags;

  @before_list = &ScheduleSelectRoutines($thorn, $implementation,
                                         $rhschedule_db->{"\U$thorn\E BLOCK_$block BEFORE"},
                                         $rhschedule_db);

  @after_list = &ScheduleSelectRoutines($thorn, $implementation,
                                        $rhschedule_db->{"\U$thorn\E BLOCK_$block AFTER"},
                                        $rhschedule_db);

  @writes_list = &ScheduleSelectRoutines($thorn, $implementation,
                                         $rhschedule_db->{"\U$thorn\E BLOCK_$block WRITES"},
                                         $rhschedule_db);

  @reads_list = &ScheduleSelectRoutines($thorn, $implementation,
                                        $rhschedule_db->{"\U$thorn\E BLOCK_$block READS"},
                                        $rhschedule_db);

  @while_list = &ScheduleSelectVars($thorn, $implementation,
                                    $rhschedule_db->{"\U$thorn\E BLOCK_$block WHILE"},
                                    $rhinterface_db);

  @if_list = &ScheduleSelectVars($thorn, $implementation,
                                 $rhschedule_db->{"\U$thorn\E BLOCK_$block IF"},
                                 $rhinterface_db);


  # Create a block so that we can have local variables.
  $buffer = "  {\n";

  # Create the timelevel array
  $buffer .= '    int cctkschedulei_tlevelarray[] = {';

  foreach $i (@$tlist)
  {
    $buffer .= "$i,";
  }

  # Add a final number to make sure we have at least one number in array
  # Can also be used to detect end of array.

  $buffer .= "0};\n";

  # add check on number of timelevels in case we were using a parameter
  for($i=0; $i < @$tlist; $i++)
  {
    $buffer .= "    if(!($$tlist[$i] >= 0 && $$tlist[$i] <= CCTK_DeclaredTimeLevels(\"$$mem_groups[$i]\")))\n";
    $buffer .= "        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n";
    $buffer .= "                   \"Tried to schedule %ld timelevels for group '%s' in schedule.ccl.\\n\"\n";
    $buffer .= "                   \"Value must be between 0 and %d (inclusive)\",\n";
    $buffer .= "                   (long)$$tlist[$i], \"$$mem_groups[$i]\", CCTK_DeclaredTimeLevels(\"$$mem_groups[$i]\"));\n";
    $buffer .= "\n";
  }

  # Start writing out the data
  if($rhschedule_db->{"\U$thorn\E BLOCK_$block TYPE"} eq "GROUP")
  {
    $prototype = '';
    $buffer .= '    CCTKi_ScheduleGroup("' . $rhschedule_db->{"\U$thorn\E BLOCK_$block NAME"} . "\",\n";
    $indent  = '                        ';
    $language = '';
  }
  elsif($rhschedule_db->{"\U$thorn\E BLOCK_$block TYPE"} eq "FUNCTION")
  {
    if($rhschedule_db->{"\U$thorn\E BLOCK_$block LANG"} =~ m:^\s*C\s*$:i )
    {
      $language = 'C';
      $function = $rhschedule_db->{"\U$thorn\E BLOCK_$block NAME"};
    }
    elsif($rhschedule_db->{"\U$thorn\E BLOCK_$block LANG"} =~ m:^\s*(F|F77|FORTRAN|F90)\s*$:i )
    {
      $language = 'Fortran';
      $function = "CCTK_FNAME(".$rhschedule_db->{"\U$thorn\E BLOCK_$block NAME"} .")";
    }
    else
    {
      if (!$rhschedule_db->{"\U$thorn\E BLOCK_$block LANG"})
      {
        $mess = "Language not specified in schedule block: " .$rhschedule_db->{"\U$thorn\E BLOCK_$block NAME"} ." in thorn: $thorn";
        &CST_error(0,$mess,'',__LINE__,__FILE__);
        return ('', '');
      }
      else
      { 
        $mess = 'Unknown language ' .$rhschedule_db->{"\U$thorn\E BLOCK_$block LANG"} ." in schedule block: ".$rhschedule_db->{"\U$thorn\E BLOCK_$block NAME"} ." in thorn: $thorn";
        &CST_error(0,$mess,'',__LINE__,__FILE__);
        return ('', '');
      }
    }
    $prototype = "extern int $function(void);\n";
    $buffer .= "    CCTKi_ScheduleFunction((void *)$function,\n";
    $indent  = '                           ';
    $buffer .= $indent;
  }
  else
  {
    $mess = 'Internal error: Unknown schedule block type ' . $rhschedule_db->{"\U$thorn\E BLOCK_$block TYPE"};
    &CST_error(0,$mess,'',__LINE__,__FILE__);
    return ('', '');
  }

  $buffer .= '"' . $rhschedule_db->{"\U$thorn\E BLOCK_$block AS"} . "\",\n";
  $buffer .= "$indent\"$thorn\",\n";
  $buffer .= "$indent\"$implementation\",\n";
  $buffer .= "$indent\"" . $rhschedule_db->{"\U$thorn\E BLOCK_$block DESCRIPTION"} . "\",\n";
  $buffer .= "$indent\"" . $rhschedule_db->{"\U$thorn\E BLOCK_$block WHERE"} . "\",\n";
  if($language ne '')
  {
    $buffer .= "$indent\"$language\",\n";
  }

  $buffer .= $indent . scalar(@$mem_groups) . ",  /* Number of STORAGE  groups   */\n";
  $buffer .= $indent . scalar(@$comm_groups) . ",  /* Number of COMM     groups   */\n";
  $buffer .= $indent . scalar(@$trigger_groups) . ",  /* Number of TRIGGERS groups   */\n";
### TR 22 Jan 2004: disabled the check for triggers of ANALYSIS routines
###                 as the absence of triggers now means to schedule
###                 such routines unconditionally
#  if (!scalar(@$trigger_groups) && $rhschedule_db->{"\U$thorn\E BLOCK_$block WHERE"} eq "CCTK_ANALYSIS")
#  {
#    $mess = "Schedule error: Scheduling at ANALYSIS in $thorn with no TRIGGERS\n";
#    $help = 'Functions or function groups scheduled in the ANALYSIS bin require TRIGGERS to be set.';
#    $help .= 'Triggers are grid variables or grid variable group names which are examined by ';
#    $help .= 'IO methods to decide whether of not execution should happen.';
#    &CST_error(0,$mess,$help,__LINE__,__FILE__);
#  }
  $buffer .= $indent . scalar(@$sync_groups)  . ", /* Number of SYNC     groups    */\n";
  $buffer .= $indent . scalar(@writes_list)   . ", /* Number of WRITES clauses     */\n";
  $buffer .= $indent . scalar(@reads_list)    . ", /* Number of READS clauses      */\n";
  $buffer .= $indent . scalar(@options)       . ", /* Number of Options            */\n";
  $buffer .= $indent . scalar(@before_list)   . ", /* Number of BEFORE   routines  */\n";
  $buffer .= $indent . scalar(@after_list)    . ", /* Number of AFTER    routines  */\n";
  $buffer .= $indent . scalar(@while_list)    . ", /* Number of WHILE    variables */\n";
  $buffer .= $indent . scalar(@if_list)       . ", /* Number of IF       variables */\n";
  $buffer .= $indent . "cctkschedulei_tlevelarray  /* Array of timelevel data for storage groups */";

  foreach $item (@$mem_groups, @$comm_groups, @$trigger_groups, @$sync_groups,
                 @writes_list, @reads_list, @options)
  {
    $buffer .= ",\n$indent\"$item\"";
  }

  $buffer .= ",\n$indent\"$tags\"";

  foreach $item (@before_list, @after_list, @while_list, @if_list)
  {
    $buffer .= ",\n$indent\"$item\"";
  }

  $buffer .= ");\n\n";

  # Close this block.
  $buffer .= "  }\n";

  return ($buffer, $prototype);
}

#/*@@
#  @routine    ScheduleStatement
#  @date       Fri Sep 17 17:38:30 1999
#  @author     Tom Goodale
#  @desc
#  Creates the code for a given schedule statement
#  @enddesc
#@@*/
sub ScheduleStatement
{
  my($thorn, $implementation, $statement, $rhinterface_db, $rhschedule_db) = @_;

  my($buffer, $prototype);
  my($groups);
  my($misc);

  # Extract the groups.
  ($groups,$misc) = &ScheduleSelectGroups($thorn, $implementation,
                                          $rhschedule_db->{"\U$thorn\E STATEMENT_$statement GROUPS"},
                                          $rhinterface_db);


  if($rhschedule_db->{"\U$thorn\E STATEMENT_$statement TYPE"} eq "STOR")
  {
    &ScheduleValidateTimeLevels($thorn, $implementation, $groups,$misc, $rhinterface_db);
    $function = "CCTKi_ScheduleGroupStorage(";
    $prototype = "";

    my $i;

    # add check on number of timelevels in case we were using a parameter
    for($i=0; $i < @$groups; $i++)
    {
      $buffer .= "  if(!($$misc[$i] >= 0 && $$misc[$i]  <= CCTK_DeclaredTimeLevels(\"$$groups[$i]\")))\n";
      $buffer .= "      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,\n";
      $buffer .= "                 \"Tried to schedule %ld timelevels for group '%s' in schedule.ccl.\\n\"\n";
      $buffer .= "                 \"Value must be between 0 and %d (inclusive)\",\n";
      $buffer .= "                 (long)$$misc[$i], \"$$groups[$i]\", CCTK_DeclaredTimeLevels(\"$$groups[$i]\"));\n";
      $buffer .= "\n";
    }

    for($i=0; $i < @$groups; $i++)
    {
      my $group     = $$groups[$i];
      my $timelevel = $$misc[$i];

      $buffer .= "  $function\"$group\",$timelevel);\n"
    }
  }
  elsif($rhschedule_db->{"\U$thorn\E STATEMENT_$statement TYPE"} eq "COMM")
  {
    print "No need to switch on Communication in $thorn\n";
    print "Communication automatically assigned for variables with storage\n";
  }
  else
  {

    $mess = "Unknown statement type '" .$rhschedule_db{"\U$thorn\E STATEMENT_$statement TYPE"};
    &CST_error(0,$mess,"",__LINE__,__FILE__);
    return ("", "");
  }

  return ($buffer, $prototype);
}


#/*@@
#  @routine    ScheduleSelectGroups
#  @date       Fri Sep 17 17:38:53 1999
#  @author     Tom Goodale
#  @desc
#  Parses a list of variable groups and selects valid ones.
#  @enddesc
#@@*/
sub ScheduleSelectGroups
{
  my($thorn, $implementation, $group_list, $rhinterface_db) = @_;
  my(@groups);
  my(@temp_list);
  my($group);
  my($other_imp, $other_thorn, $foundit, $block);
  my($misc,@misc_list);

  @temp_list = split(/[,\s\n]+/, $group_list);

  foreach $group (@temp_list)
  {
    next if($group =~ m:^\s*$:);

    $other_imp = "";

#    print "DEBUG: $thorn,$implementation - group was $group ";

    # Strip off extra miscellaneous info on group- i.e [...] bit
    $group =~ m/^([^[]+)(\[([^]]*)\])?$/;

    $group = $1;

    my $misc = $3;

    if(! defined $misc)
    {
      $misc = "";
    }

#    print ", group now is $group with misc '$misc'\n"; #DEBUG

    if($group =~ m/^(.+)::(.+)$/)
    {
      $other_imp=$1;
      $group = $2;

#      print "DEBUG: $thorn, $other_imp, $group, $misc\n";

      if(($other_imp !~ m:^\s*$thorn\s*$:i) && ($other_imp !~ m:^\s*$implementation\s*$:i))
      {
        # The name has been given completely specified but it isn't this thorn.

        if($rhinterface_db->{"IMPLEMENTATION \U$implementation\E ANCESTORS"} =~ m:\b$other_imp\b:i)
        {
          $block = "PUBLIC";
        }
        elsif($rhinterface_db->{"IMPLEMENTATION \U$implementation\E FRIENDS"} =~ m:\b$other_imp\b:i)
        {
          $block = "PROTECTED";
        }
        else
        {
          $mess = "Schedule error: Thorn $thorn - group $other_imp\:\:$group doesn't exist.";
          $help = "Check thorn $thorn inherits from implementation $other_imp";
          &CST_error(0,$mess,$help,__LINE__,__FILE__);
          next;
        }

        $rhinterface_db->{"IMPLEMENTATION \U$other_imp\E THORNS"} =~ m:(\w+):;
        $other_thorn = $1;

        if($rhinterface_db->{"\U$other_thorn\E $block GROUPS"} =~ m:\b$group\b:i)
        {
          push(@groups, "$other_imp\::$group");
          push(@misc_list, $misc);
          next;
        }
        else
        {
          $mess = "Schedule error: Thorn $thorn - group $other_imp\:\:$group doesn't exist.\n";
          &CST_error(0,$mess,"",__LINE__,__FILE__);
          next;
        }
      }
    }

    if($rhinterface_db->{"\U$thorn\E PRIVATE GROUPS"} =~ m:\b$group\b:i)
    {
      push(@groups, "$thorn\::$group");
      push(@misc_list, $misc);
    }
    elsif($rhinterface_db->{"\U$thorn\E PROTECTED GROUPS"} =~ m:\b$group\b:i)
    {
      push(@groups, "$implementation\::$group");
      push(@misc_list, $misc);
    }
    elsif($rhinterface_db->{"\U$thorn\E PUBLIC GROUPS"} =~ m:\b$group\b:i)
    {
      push(@groups, "$implementation\::$group");
      push(@misc_list, $misc);
    }
    elsif($other_imp eq "")
    {
      $foundit = 0;
      # Check ancestors and friends
      foreach $other_imp (split(" ", $rhinterface_db->{"IMPLEMENTATION \U$implementation\E ANCESTORS"}))
      {
        $rhinterface_db->{"IMPLEMENTATION \U$other_imp\E THORNS"} =~ m:(\w+):;
        $other_thorn = $1;

        if($rhinterface_db->{"\U$other_thorn\E PUBLIC GROUPS"} =~ m:\b$group\b:i)
        {
          push(@groups, "$other_imp\::$group");
          push(@misc_list, $misc);
          $foundit = 1;
          last;
        }
      }
      if(! $foundit)
      {
        foreach $other_imp (split(" ", $rhinterface_db->{"IMPLEMENTATION \U$implementation\E FRIENDS"}))
        {
          $rhinterface_db->{"IMPLEMENTATION \U$other_imp\E THORNS"} =~ m:(\w+):;
          $other_thorn = $1;

          if($rhinterface_db->{"\U$other_thorn\E PROTECTED GROUPS"} =~ m:\b$group\b:i)
          {
            push(@groups, "$other_imp\::$group");
            push(@misc_list, $misc);
            $foundit = 1;
            last;
          }
        }
      }
      if(! $foundit)
      {
        $mess = "Schedule error: Thorn $thorn - group $group doesn't exist.";
        $help = "Check $group really is in thorn $thorn. Groups from other thorns ";
        $help .= "need to be specified using <implementation>\:\:$group and ";
        $help .= "<implementation> must be inherited by your thorn.";
        &CST_error(0,$mess,$help,__LINE__,__FILE__);

      }
    }
    else
    {
      $mess = "Schedule error: Thorn $thorn - group $group doesn't exist.";
      &CST_error(0,$mess,"",__LINE__,__FILE__);

    }
  }

  return (\@groups,\@misc_list);
}

#/*@@
#  @routine    ScheduleSelectRoutines
#  @date       Fri Sep 17 17:39:29 1999
#  @author     Tom Goodale
#  @desc
#  Parses a list of schedule routines/groups.
#  FIXME - should validate
#  @enddesc
#@@*/
sub ScheduleSelectRoutines
{
  my($thorn, $implementation, $routine_list, $rhschedule_db) = @_;
  my(@routines);
  my(@temp_list);
  my($routine);

  @temp_list = split(/[,\s\n]+/, $routine_list);

  foreach $routine (@temp_list)
  {
    next if($routine =~ m:^\s*$:);

    push(@routines, $routine);

  }

  return @routines;
}

#/*@@
#  @routine    ScheduleSelectVars
#  @date       Fri Sep 17 17:39:58 1999
#  @author     Tom Goodale
#  @desc
#  Parses a list of variables
#  FIXME - should validate
#  @enddesc
#@@*/
sub ScheduleSelectVars
{
  my($thorn, $implementation, $var_list, $rhinterface_db) = @_;
  my(@vars);
  my(@temp_list);
  my($var);

  @temp_list = split(/[,\s\n]+/, $var_list);

  foreach $var (@temp_list)
  {
    next if($var =~ m:^\s*$:);

    push(@vars, $var);
  }

  return @vars;
}

#/*@@
#  @routine    ScheduleValidateTimeLevels
#  @date       Tue Apr 16 15:22:02 2002
#  @author     Tom Goodale
#  @desc
#  Validate the timelevel specifiers for a group list.
#  @enddesc
#@@*/
sub ScheduleValidateTimeLevels
{
  my($thorn, $implementation, $groups,$timelevels_list, $rhinterface_db) = @_;

  my $i;

  my $return_code;

  $return_code = 0;

  for($i = 0; $i < @$groups; $i++)
  {
    my $group = $$groups[$i];
    my $timelevels = $$timelevels_list[$i];

#    print "DEBUG: validate $thorn,$implementation,$group,$timelevels\n";

    if($timelevels !~ /^[[:alpha:]_][[:word:]]*$/ && $timelevels !~ /^\d*$/)
    {
      &CST_error(0,"Invalid timelevel specifier '$timelevels' in schedule.ccl of thorn '$thorn'","",__LINE__,__FILE__);
      $return_code++;
      next;
    }

    $group =~ m/^(.+)::(.+)$/;

    my $imp = $1;
    my $groupname = $2;

    my $allowed_timelevels;

    if($1 =~ /^$thorn$/i)
    {
      $allowed_timelevels =  $rhinterface_db->{"\U$thorn GROUP $groupname TIMELEVELS\E"};
    }
    else
    {
      @thornlist = split(" ",$rhinterface_db->{"IMPLEMENTATION \U$imp\E THORNS"});

      $allowed_timelevels =  $rhinterface_db->{"\U$thornlist[0] GROUP $groupname TIMELEVELS\E"};
    }

    # If the maximum number of timelevels is 1, the timelevel specification can be empty.
    if($timelevels eq "" && $allowed_timelevels == 1)
    {
      $$timelevels_list[$i] = 1;
      next;
    }
    elsif($timelevels eq "")
    {
      &CST_error(0,"Timelevel specifier missed for group '$group' in schedule.ccl of thorn '$thorn'\n".
                   "This group has more than one timelevel, so a timelevel specifier is mandatory.\n" .
                   "You may specify a maximum of $allowed_timelevels timelevels\n" .
                   "e.g. try $group\[$allowed_timelevels]\n" .
                   "Note that you should only activate the number of timelevels necessary\n" .
                   "for your scheme, which may be less than this maximum."
                 ,"",__LINE__,__FILE__);
      $return_code++;
    }
    elsif($timelevels =~ /^[[:alpha:]_][[:word:]]*$/ || ($timelevels > 0 && $timelevels  <= $allowed_timelevels))
    {
      next;
    }
    else
    {
      if($allowed_timelevels > 1)
      {
        &CST_error(0,"Tried to schedule $timelevels timelevels for group '$group' in schedule.ccl of thorn '$thorn'\n" .
                     "Value must be between 1 and $allowed_timelevels (inclusive)","",__LINE__,__FILE__);
      }
      else
      {
        &CST_error(0,"Tried to schedule $timelevels timelevels for group '$group' in schedule.ccl of thorn '$thorn'\n" .
                     "This variable has one timelevel only","",__LINE__,__FILE__);
      }
      $return_code++;
    }
  }

  return $return_code;
}

1;
