#! /usr/bin/perl -w
#
# /*@@
#   @file      CreateFunctionBindings.pl
#   @date      Sun Feb 16 01:18:02 2003
#   @author    Ian Hawke
#   @desc
#   Does all the work for function aliasing.
#   Something approximating documentation follows.
#   @enddesc
#   @version $Id$
# @@*/
#
# The structure of function aliasing is described by this piece of
# ASCII art.
#
#                      Calling thorn calls 'bla'
#                            /       \
#                           /         \
#                          /           \
#                         /             \
#                        /               \
#                       /                 \
#                      /                   \
#               Function                Function
#                  'bla'            'CCTK_FCALL CCTK_FNAME(bla)'
#                     |                     |
#                     |                     |
#                     |                     |
#                     |                     |
#                     |                     |
#               Function pointer       Function pointer
#               Thorn_bla_From_C       Thorn_bla_From_F
#                      \                   /
#                       \                 /
#                        \               /
#                         \             /
#                          \           /
#                           \         /
#                            \       /
#                            Thorn_bla
#
# In words. The calling thorn will link to either 'bla' or the Fortran
# equivalent, 'CCTK_FCALL CCTK_FNAME(bla)'. However, unlike the Cactus
# standard, the Fortran wrapper does not call the C function.
#
# Instead, in this case both 'bla' and the Fortran version will call a
# function pointer, which will point either to the providing function,
# or to a wrapper to the providing function.
#
# The providing function is Thorn_bla. The CST will generate wrappers
# for this function which work if called from C (Thorn_bla_From_C) or
# Fortran (Thorn_bla_From_F).
#
# Then at run time a "registration" function will set the function
# pointers Thorn_bla_From_C and Thorn_bla_From_F to either the
# providing function itself, or a wrapper for the providing function,
# as needed.
#
# The functions pointer to by the 'Thorn_bla_From_*' pointers will be
# in C. They know which language they are being called from, and also
# which language Thorn_bla has, and they know this _at CST time_. So
# they will change the arguments and calling sequence appropriately.
#
# Thus we need the following.
#
# 1) For every thorn that PROVIDES a function, generate the following
#    functions.
#
#    a) Wrappers to which Thorn_bla_From_C or Thorn_bla_From_F will
#       point, as needed.  These receive their arguments in the
#       appropriate C/Fortran style and then call Thorn_bla in the
#       (possibly different) C/Fortran style.
#
#    b) A function that "registers" all the providing functions.
#
# 2) For every function that is USEd, generate the following.
#
#    a) A function 'bla' and a Fortran wrapper for it.
#
#    b) A "registration" function that sets the appropriate function
#       pointer. This should return an error if it has already been set.
#
#    c) An "IsFunctionAliased" routine for _both_ the 'bla' and Fortran
#       functions. These can just check if the pointer is non-null.
#
# POINTS TO NOTE:
#
# - For safety, assume that any thorn that PROVIDEs a function also USEs it.
#   This will ensure that the function pointers etc. exist.
#
##########################################################################
#
# SOME EXPLANATION OF BINDING FILES GENERATED
# -------------------------------------------
#
# in <config>/bindings/Functions/AliasedFunctions.c :
# In this file the actual routines that are called by the USEing thorns are
# defined:
#
# <CCTK_TYPE> <MyAliasedFunctionName> (<ItsArguments>);
# <CCTK_TYPE> CCTK_FCALL CCTK_FNAME(<MyAliasedFunctionName>) (<ItsArguments>);
#
# Each calls a corresponding C or Fortran function pointer, respectively:
#
# static <CCTK_TYPE> (*<MyAliasedFunctionName>_C_Wrapper);
# static <CCTK_TYPE> (*<MyAliasedFunctionName>_F_Wrapper);
#
# The functions:
#
# Alias<MyAliasedFunctionName>_C(function_pointer)
# Alias<MyAliasedFunctionName>_F(function_pointer)
#
# set the pointers <MyAliasedFunctionName>_*_Wrapper to the
# function_pointer argument (for * equal to C and F respectively).
#
#
# in <config>/bindings/Functions/<ThornName>.c :
#
# For a function that is provided in C, the function
# Register_<ThornName>(void) calls
#
# Alias<MyAliasedFunctionName>_C( <MyProvidingFunctionName> )
# Alias<MyAliasedFunctionName>_F( CCTK_Wrapper_CtoF_<MyProvidingFunctionName>)
#
# The (somewhat) analogous thing is done for aliased functions that
# are provided in Fortran.
#
# in <config>/bindings/include/<ThornName>_Prototypes.h :
#
# The prototypes for USEd aliased functions are located here.
#
##########################################################################
#
# STRUCTURES
# ----------
#
# At the top level is the FunctionDatabase. This is passed in from the
# CST and is then parsed to create a slightly different structure (much
# duplication of effort, but hey). The FunctionDatabase is a hash where
# the keys are the thornnames and the values are FunctionLists.
#
# A FunctionList is another hash where the keys are the function names
# and the values are Functions (also hashes).
#
# A Function hash is probably the fundamental structure.
# The key:values are:
#
# Name            : This agrees with the key in the function list and
#                   is probably redundant, but simpler to keep here as
#                   well (plus may be some whitespace stripped). (String)
# Return Type     : The return type. This MUST be a scalar type. (String)
# Used            : Whether this function is REQUIRED/USED/not needed
#                   by this thorn (2/1/0)
# Provided        : Whether this function is PROVIDED by this thorn (1/0)
# Provider        : The name of the providing function. ONLY EXISTS if
#                   Provided = 1. (String)
# Provider Language : Language of providing function. ONLY EXISTS if
#                     Provided = 1. (String: Fortran / C)
# Arguments       : A (reference to a) list of Arguments for the Function.
# Strings         : The number of string arguments (should not be > 3) (int)
# String pointers : The number of pointer to string arguments (OBSOLETE) (int)
#
#
# An Argument is also a hash. The list of possible key/values is
#
# Name            : The (dummy) name of the argument or a reference to a
#                   Function if it's a function pointer argument (String/hash)
# Type            : The (scalar) type of the argument of return type of the
#                   function pointer argument. (String)
# Function pointer: Whether the argument is a function pointer (1/0)
# Is Array        : Whether the argument is an ARRAY (1/0)
# String          : Whether the argument is a string (1/0)

# For debugging:
my $debug = 0;
my $indent_level = 0;
sub debug_print
{
  $debug and print '#'.' 'x$indent_level."@_\n";
}


sub CreateFunctionBindings
{
  use strict;

  my($bindings_dir, $rhinterface_db) = @_;

  my($start_dir,$Function);

  my($dataout);
  my($function_db);
  my($registerfiles);

  $registerfiles = "AliasedFunctions.c IsFunctionAliased.c RegisterThornFunctions.c";

######################################################################
# Create the database
######################################################################

  $function_db = &FunctionDatabase($rhinterface_db);

  CheckRequiredFunctions($function_db);

###
# This should have returned a FunctionList (as defined above)
# for every thorn. That is, a hash with the key being the thorn name
# and the value a FunctionList
###

######################################################################
# Create all the directories
######################################################################

  if(! -d $bindings_dir)
  {
    mkdir("$bindings_dir", 0755) || die "Unable to create $bindings_dir";
  }
  $start_dir = `pwd`;

  chdir $bindings_dir;

  if(! -d "Functions")
  {
    mkdir("Functions", 0755) || die "Unable to create Functions directory";
  }

  if(! -d "include")
  {
    mkdir("include", 0755) || die "Unable to create include directory";
  }

######################################################################
# Create the appropriate files
######################################################################

###
# Create the stuff for the provided functions
###

  my $thorn;

  foreach $thorn (sort keys %{$function_db})
  {
    $debug and print "thorn is $thorn\n";

    my $localfns = keys %{$function_db->{$thorn}};
    if ($localfns)
    {
      $debug and print "  localfns = $localfns\n";
      $dataout = &ProvidedFunctions($thorn,$function_db->{"$thorn"});
      my $filename = "${thorn}_Functions.c";
      &WriteFile("Functions/$filename",\$dataout);
      $registerfiles.=" $filename";
    }
  }

###
# Create the master file to register the provided functions
###

  $dataout = &RegisterAllFunctions($function_db);
  &WriteFile("Functions/RegisterThornFunctions.c",\$dataout);

###
# Create the stuff for the used functions
###

  $dataout = &AliasedFunctions($function_db);
  &WriteFile("Functions/AliasedFunctions.c",\$dataout);

###
# Create the master IsFunctionAliased file
###

  $dataout = &IsFunctionAliased($function_db);
  &WriteFile("Functions/IsFunctionAliased.c",\$dataout);

###
# Create the master header file (i.e., all the USEd functions)
###

  $dataout = &ThornMasterIncludes($function_db);

###
# Create the prototype header file for all thorns that USE a
# function.
###

  foreach $thorn (keys %{$function_db})
  {
    $dataout = &UsesPrototypes($thorn,$function_db->{"$thorn"});
    my $filename = "${thorn}_Prototypes.h";
    &WriteFile("include/$filename",\$dataout);
  }


######################################################################
# Create the make.code.defn
######################################################################

  $dataout = "SRCS = $registerfiles";
  &WriteFile("Functions/make.code.defn",\$dataout);

  chdir $start_dir;

  return;
}

#/*@@
#  @routine    FunctionDatabase
#  @date       Sun Feb 16 01:29:48 2003
#  @author     Ian Hawke
#  @desc
#  Parses the standard bindings interface database (interface_db) into
#  the internal FunctionDatabase structure defined above.
#
#  As arguments the subroutine takes the interface_db.
#
#  It returns just the FunctionDatabase structure.
#
#  @enddesc
#@@*/

sub FunctionDatabase
{
  use strict;

  my($interface_db) = @_;

  my $thorn;
  my $FunctionDatabase={};

  foreach $thorn (split(' ',$interface_db->{'THORNS'}))
  {
    my $FunctionList={};
    my $FunctionName;
    my $ThornUses=0;
    my $ThornProvides=0;
    foreach $FunctionName (split(' ',($interface_db->{"\U$thorn FUNCTIONS\E"})))
    {
      my $Function={};
      my $nstrings;
      my $nstringptrs;
      my $warnings;
      my ($ReturnType,$Arguments,@arglist);
      $Arguments = $interface_db->{"\U${thorn} FUNCTION\E $FunctionName ARGS"};

#      &debug_print("FunctionDatabase: calling ParseArgumentsList with thorn=[$thorn] FunctionName=[$FunctionName] Arguments=[$Arguments]\n");
      ($warnings,$nstrings,$nstringptrs,@arglist)=&ParseArgumentsList($Arguments, $thorn, $FunctionName);
      $Function->{"Strings"} = $nstrings;
      $Function->{"String pointers"} = $nstringptrs;

      $ReturnType = $interface_db->{"\U${thorn} FUNCTION\E $FunctionName RET"};
      # turn 'VOID' into lower-case
      $ReturnType = lc $ReturnType;
      # turn all return types except 'void' into upper-case CCTK types
      $ReturnType = uc $ReturnType if ($ReturnType ne 'void ');

      $FunctionName =~ /([a-zA-Z][a-zA-Z0-9_]*)/;
      $Function->{"Name"}=$1;
      if ($warnings)
      {
        my $message = "The aliased function ".$Function->{"Name"}." has an error.\n".$warnings;
        &CST_error(0,$message,'',__LINE__,__FILE__);
      }
      if ( ("\L$Function->{\"Name\"}\E" eq $Function->{"Name"}) ||
           ("\U$Function->{\"Name\"}\E" eq $Function->{"Name"}) )
      {
          my $message = "An aliased function name must contain at least one upper AND one lower case letter.\n The function ".$Function->{"Name"}." does not.";
          &CST_error(0,$message,'',__LINE__,__FILE__);
      }
#
#  FIXME:: the ISO standard says that a function name should be
#          unique within the first 31 characters, which together with
#          wrapper stuff gives us only 22 characters to play with.
#          This could be gotten around with Tom's registry idea, or by
#          forcing a 22 character limit. But for the moment we just
#          comment this bit out.
#
#      if ( length($Function->{"Name"}) > 20 )
#      {
#          my $message = "An aliased function name must be less than 21 characters.\n The function ".$Function->{"Name"}." is not.";
#          &CST_error(0,$message,'',__LINE__,__FILE__);
#      }
      $Function->{"Provided"} = $interface_db->{"\U$thorn\E PROVIDES FUNCTION"} =~ / $FunctionName /;
      if ($Function->{"Provided"})
      {
        my $provider = $interface_db->{"\U$thorn PROVIDES FUNCTION\E $FunctionName WITH"};
        my $language = $interface_db->{"\U$thorn PROVIDES FUNCTION\E $FunctionName LANG"};
        $provider =~ /([a-zA-Z][a-zA-Z0-9_]*)/;
        $Function->{"Provider"}=$1;
        $language =~ /([a-zA-Z]+)/;
        $Function->{"Provider Language"}=$1;
        $ThornProvides++;
      }
      if (   $interface_db->{"\U$thorn\E USES FUNCTION"}     =~ / $FunctionName /
          || $interface_db->{"\U$thorn\E REQUIRES FUNCTION"} =~ / $FunctionName /
          || $interface_db->{"\U$thorn\E PROVIDES FUNCTION"} =~ / $FunctionName /)
      {
        $Function->{"Used"} = $interface_db->{"\U$thorn\E REQUIRES FUNCTION"} =~ / $FunctionName / ? 2 : 1;
        $ThornUses++;
      }
      else
      {
        $Function->{"Used"}=0;
      }
      $Function->{"Arguments"}=\@arglist;
      $ReturnType =~ /([a-zA-Z][a-zA-Z0-9_]*)/;
      $Function->{"Return Type"}=$1;
      if ($ReturnType =~ /:ARRAY/)
      {
          my $message = "An aliased function may not return an ARRAY type.\n The function ".$Function->{"Name"}." attempts to.";
          &CST_error(0,$message,'',__LINE__,__FILE__);
#        $Function->{"Return Type"}.="*";
      }
#  At this point we should be checking that any duplicated function
#  prototypes are consistent.
      my $KnownThorn;
      foreach $KnownThorn (sort keys %{$FunctionDatabase})
      {
        my $KnownList = $FunctionDatabase->{"$KnownThorn"};
        my $KnownFunctionName;
        foreach $KnownFunctionName (sort keys %{$KnownList})
        {
          my %KnownFunction = %{$KnownList->{$KnownFunctionName}};
          if ($KnownFunction{"Name"} eq $Function->{"Name"})
          {
            if (!($KnownFunction{"Return Type"} eq
                  $Function->{"Return Type"}))
            {
              &CST_error(0,"The prototypes for the aliased function \'".$KnownFunction{"Name"}."\'\n     given by thorns \' ".$thorn."\' and \'".$KnownThorn."\' are inconsistent.\n     The return types disagree.");
            }
            if (&CompareArguments($KnownFunction{"Arguments"},
                                 $Function->{"Arguments"}))
            {
              &debug_print("The prototypes for the aliased function \'".$KnownFunction{"Name"}."\'\n     given by thorns \'".$thorn."\' and \'".$KnownThorn."\' are inconsistent.\n     The argument lists disagree.");
              &CST_error(0,"The prototypes for the aliased function \'".$KnownFunction{"Name"}."\'\n     given by thorns \'".$thorn."\' and \'".$KnownThorn."\' are inconsistent.\n     The argument lists disagree.");
            }
          }
        }
      }
      $FunctionList->{$FunctionName}=$Function;
    }
    $FunctionDatabase->{$thorn}=$FunctionList;
  }
  return $FunctionDatabase;

}

#/*@@
#  @routine    CheckRequiredFunctions
#  @date       Fri 23 April 2004
#  @author     Thomas Radke
#  @desc
#  Checks that aliased functions which are REQUIRED by some thorn
#  are also provided by some other.
#  @enddesc
#@@*/

sub CheckRequiredFunctions
{
  use strict;
  my %FunctionDatabase = %{$_[0]};

  foreach my $thorn (keys %FunctionDatabase)
  {
    foreach my $Function (values %{$FunctionDatabase{$thorn}})
    {
      # need to check only if this thorn doesn't provide the function itself
      if ($Function->{'Used'} == 2 && ! $Function->{'Provided'})
      {
        my $is_provided = 0;

        # now go through the function database again to find the providing thorn
        foreach my $provider (keys %FunctionDatabase)
        {
          next if ($provider eq $thorn);

          foreach my $ProvidingFunction (values %{$FunctionDatabase{$provider}})
          {
            $is_provided = $Function->{'Name'} eq $ProvidingFunction->{'Name'}
                           && $ProvidingFunction->{'Provided'};
            last if ($is_provided);
          }
          last if ($is_provided);
        }
        &CST_error(0,"Aliased function \'$Function->{'Name'}\' required " .
                   "by thorn '$thorn' is not provided by any thorn in your " .
                   "thornlist\n", '',__LINE__,__FILE__)
          if (! $is_provided);
      }
    }
  }
}


#/*@@
#  @routine    ParseArgumentsList
#  @date       Sun Feb 16 01:37:55 2003
#  @author     Ian Hawke
#  @desc
#  Takes the scalar string that is the list of arguments for a
#  given function and parses it for name, return type and attributes.
#
#  The arguments to this function is just the string enclosed by braces
#  in the function declaration.
#
#  The return type is a list of Arguments as defined above, together with
#  the number of strings (no more than 4) and pointers to strings
#  (should not be any in the current version, but code is left for if/when
#  someone works out how to do this properly), and any warnings that
#  are required.
#
#  @enddesc
#@@*/

sub ParseArgumentsList
{
  use strict;

  $debug and $indent_level = 2;

  my($Arguments) = shift;
  my($Thorn) = shift;
  my($Function) = shift;

#  &debug_print("ParseArgumentsList: Arguments=[$Arguments] Thorn=[$Thorn] Function=[$Function]");

  my @ArgList=();

  my $nfptrs = 0;
  my $warnings = "";
  my @fptrargs = ();
  if ($Arguments =~ s/CCTK_FPOINTER//g)
  {
    &debug_print("$Thorn:$Function:$Arguments\n");
    while ($Arguments =~ s/(.*?)\s*(\(.*?\))(.*)/\1FPTRARGS\3/)
    {
      &debug_print("$Thorn:$Function:$Arguments\n");
      my $tempargs = $2;
      $tempargs =~ s/\((.*)\)/\1/;
      push(@fptrargs,$tempargs);
      $nfptrs++; # QUERY: This is set but never used.
    }
    &debug_print("$Thorn:$Function:$Arguments\n");
  }

  my @DummyList = split(',',$Arguments);
  my $DummyArg;

  my $nstrings = 0;
  my $nstringptrs = 0;

  $nfptrs=0;

  foreach $DummyArg (@DummyList)
  {
    &debug_print("$Thorn:$Function:$DummyArg\n");
    if ($DummyArg =~ /\S/) # ignore empty argument list
    {
      my $Arg = &ParseArgument($DummyArg, $Thorn, $Function);
      $Arg->{"Function pointer"} = 0;
      push(@ArgList,$Arg);
      if ($Arg->{"Name"} =~ /FPTRARGS/)
      {
        $Arg->{"Name"} =~ s/(.*)FPTRARGS/\1/;
        my $Name = $Arg->{"Name"};
        $Arg->{"Name"} = {"Name"=>$Name,
                          "Provided"=>0,
                          "Used"=>0,
                          "Return Type"=>$Arg->{"Type"}};
        my ($extrawarnings,$nstrings,$nstringptrs,@arglist)=&ParseArgumentsList($fptrargs[$nfptrs],$Thorn,$Function);
        $warnings .= $extrawarnings;
        $Arg->{"Name"}{"Strings"} = $nstrings;
        $Arg->{"Name"}{"String pointers"} = $nstringptrs;
        $Arg->{"Name"}{"Arguments"} = \@arglist;
        $Arg->{"Function pointer"} = 1;
        $nfptrs++;
      }
      if (!$Arg->{"Is Array"})
      {
        $nstrings += $Arg->{"String"};
      }
      else
      {
        $nstringptrs += $Arg->{"String"};
      }
      if ( ($nstrings)&&(!$Arg->{"String"}) )
      {
        $warnings .= "The argument list contains CCTK_STRINGs that are not at the end.";
      }
    }
  }

#  if ($debug)
#  {
#    print "ArgList is:\n";
#    foreach $DummyArg (@ArgList)
#    {
#      print $DummyArg->{"Type"}." ".$DummyArg->{"Name"}." ";
#    }
#    print "\n";
#  }

  if ( ($nstrings > 3) || ($nstringptrs > 3) )
  {
    $warnings .= "The argument list contains more than 3 string arguments.";
  }

  return ($warnings,$nstrings,$nstringptrs,@ArgList);
}

#/*@@
#  @routine    ParseArgument
#  @date       Sun Feb 16 01:41:08 2003
#  @author     Ian Hawke
#  @desc
#  Parses an individual Argument.
#
#  ParseArgumentsList splits the full list by commas (plus other stuff
#  for function pointers). The individual arguments are passed to here.
#
#  The input is just a string. The return is a reference to an Argument
#  as defined above.
#
#  @enddesc
#@@*/

sub ParseArgument
{
  use strict;

  $debug and $indent_level = 4;

  my($DummyArgument) = shift;
  my($Thorn) = shift;
  my($Function) = shift;

  my $Argument = {};

  my($type,$name,$fpointer,$intent,$extra);

  if ($DummyArgument =~ /FPTRARGS/)
  {
    ($type,$intent,$name,$extra) = split(' ',$DummyArgument);
    # QUERY: is $fpointer supposed to be set here?
    $fpointer = 1;
    &debug_print("$Thorn--ParseArgument: (fn pointer) type=$type name=$name");
  }
  elsif ($DummyArgument =~ s/\bARRAY\b//)
  {
    ($type,$intent,$name,$extra) = split(' ',$DummyArgument);
    $Argument->{"Is Array"} = 1;
    &debug_print("$Thorn--ParseArgument: $name is ARRAY");
  }
  else
  {
    ($type,$intent,$name,$extra) = split(' ',$DummyArgument);
    $Argument->{"Is Array"} = 0;
  }

#   if ($type =~ s/\s*(.*):ARRAY\s*/\1/)
#   {
#     $Argument->{"Is Array"} = 1;
#   }
#   else
#   {
#     $Argument->{"Is Array"} = 0;
#   }

  $Argument->{"Type"} = $type;

  if ($type =~ m/CCTK_STRING/)
  {
    $Argument->{"String"} = 1;
  }
  else
  {
    $Argument->{"String"} = 0;
  }

  $Argument->{"Name"} = $name;

  if ($name) # (meaning argument has three whitespace separated components)
  {
    if ($intent !~ /^(IN|OUT|INOUT)$/)
    {
      my $message = "Thorn $Thorn, Function $Function:\nThe intent statement must be either IN, OUT or INOUT.\nThe argument \"$DummyArgument\" has the wrong type.";
      &CST_error(0,$message,'',__LINE__,__FILE__);
    }
    $Argument->{"Intent"} = $intent;
  }
  else # (usually meaning intent is missing, so $intent holds the name and $name is empty)
  {
    $intent =~ s/FPTRARGS//; #strip off the FPTARGS if this happens for a function pointer
    my $message = "Thorn $Thorn, Function $Function:\nEvery argument must contain an intent statement of type IN, OUT or INOUT.\n The argument \"$intent\" does not.";
    &CST_error(0,$message,'',__LINE__,__FILE__);
  }

  if ($extra) # too many arguments; probably a comma missed (see PR 1886)
  {
    my $message = "Thorn $Thorn, Function $Function:\nThe argument has too many specifications: should be type - intent - name.\nThe argument \"$DummyArgument\" has too many.";
      &CST_error(0,$message,'',__LINE__,__FILE__);
  }

  if ($fpointer)
  {
    $Argument->{"Function Pointer"} = 1;
  }
  else
  {
    $Argument->{"Function Pointer"} = 0;
  }
  if ($type !~ /(\bCCTK_INT$)|(\bCCTK_REAL$)|(\bCCTK_COMPLEX$)|(\bCCTK_POINTER$)|(\bCCTK_POINTER_TO_CONST$)|(\bCCTK_STRING$)/)
  {
    my $message = "Thorn $Thorn, Function $Function:\nAn argument in an aliased function must be one of the allowed CCTK types.\nThese are CCTK_INT, CCTK_REAL, CCTK_COMPLEX, CCTK_POINTER, CCTK_POINTER_TO_CONST, or CCTK_STRING.\nThe argument ".$Argument->{"Name"}." has type \"$type\".";
    if ($type =~ /:/)
    {
      $message .= "\n(The older \"${type}ARRAY\" should be replaced with \"$type ARRAY\".)";
    }
    &CST_error(0,$message,'',__LINE__,__FILE__);
  }
  if ( ($type =~ /CCTK_STRING/) && ($Argument->{"Is Array"}) )
  {
    my $message = "An argument in an aliased function may not have the CCTK_STRING:ARRAY type.\nThe argument ".$Argument->{"Name"}." does.";
    &CST_error(0,$message,'',__LINE__,__FILE__);
  }

#  print "ParseArgument:\n";
#  print $Argument->{"Name"}." ".$Argument->{"Type"}." ".$Argument->{"Is Array"}." ".$Argument->{"Function Pointer"}."\n";

  return $Argument;
}

#/*@@
#  @routine    CompareArguments
#  @date       Sun Feb 16 01:41:08 2003
#  @author     Ian Hawke
#  @desc
#  Takes two argument lists and checks that they are the same.
#
#  Returns the number of arguments that disagree.
#
#  @enddesc
#@@*/

sub CompareArguments
{
  use strict;

  my @Arguments1 = @{$_[0]};
  my @Arguments2 = @{$_[1]};

  my $num_errors = 0;

  for (my $i=0;$i<@Arguments1;$i++) {
    my $Arg1 = $Arguments1[$i];
    my $Arg2 = $Arguments2[$i];
    &debug_print("arg1: ".$Arg1->{"Type"}." ".$Arg1->{"Intent"}." ".$Arg1->{"Name"}."\n");
    &debug_print("arg2: ".$Arg2->{"Type"}." ".$Arg2->{"Intent"}." ".$Arg2->{"Name"}."\n");
  }

  if (!($#Arguments1 == $#Arguments2))
  {
    $num_errors = 10;
  }
  else
  {
    for (my $argnum = 0; $argnum < @Arguments1; $argnum++)
    {
      my $Arg1 = $Arguments1[$argnum];
      my $Arg2 = $Arguments2[$argnum];
      if ($Arg1->{"Function pointer"})
      {
        $num_errors++ if (! $Arg2->{"Function pointer"});
      }
      elsif (
             (!($Arg1->{"Type"} eq $Arg2->{"Type"})) ||
             (!($Arg1->{"Intent"} eq $Arg2->{"Intent"}))
            )
      {

        &debug_print("Errors in arguments:\n".
          $Arg1->{"Type"}." ".$Arg2->{"Type"}."\n".
          $Arg1->{"Intent"}." ".$Arg2->{"Intent"}."\n".
          $Arg1->{"Function pointer"}." ".$Arg2->{"Function pointer"}
                    );

        $num_errors++;
      }
    }
  }

  return $num_errors;
}

#/*@@
#  @routine    FunctionDatabaseFake
#  @date       Sun Feb 16 01:45:37 2003
#  @author     Ian Hawke
#  @desc
#  A debugging routine that creates a FunctionList.
#  @enddesc
#@@*/

#sub FunctionDatabaseFake
#{
#
#
#  use strict;
#
#  my $Arg1;
#  my $Arg2;
#  my $Arg3;
#  my $Arg4;
#  my $Arg5;
#
#  my $ArgList1;
#  my $ArgList2;
#  my $ArgList3;
#
#  my $Function1;
#  my $Function2;
#  my $Function3;
#
#  my $FunctionList1;
#
#  $Arg1 = {"Function"=>"0",
#           "Type"=>"CCTK_INT",
#           "Array"=>"0",
#           "Name"=>"x"};
#  $Arg2 = {"Function"=>"0",
#           "Type"=>"CCTK_INT",
#           "Array"=>"0",
#           "Name"=>"y"};
#  $Arg3 = {"Function"=>"1",
#           "Type"=>"CCTK_REAL",
#           "Array"=>"0",
#           "Name"=>"fnptr"};
#  $Arg4 = {"Function"=>"0",
#           "Type"=>"CCTK_REAL",
#           "Array"=>"0",
#           "Name"=>"a"};
#  $Arg5 = {"Function"=>"0",
#           "Type"=>"CCTK_REAL",
#           "Array"=>"1",
#           "Name"=>"b"};
#  $ArgList1 = [$Arg1,$Arg2];
#  $ArgList2 = [$Arg3,$Arg5];
#  $ArgList3 = [$Arg4];
#  $Function1 = {"Name"=>"Sum",
#                "Aliased"=>"1",
#                "Return Type"=>"CCTK_INT",
#                "Arguments"=>$ArgList1,
#                "Thorn Uses"=>"1",
#                "Thorn Provides"=>"1",
#                "Providing Fn"=>"AddItUp",
#                "Providing Lang"=>"Fortran"};
#  $Function2 = {"Name"=>"Integrate",
#                "Aliased"=>"0",
#                "Return Type"=>"CCTK_REAL",
#                "Arguments"=>$ArgList2,
#                "Thorn Uses"=>"0",
#                "Thorn Provides"=>"0",
#                "Providing Fn"=>"General_Integrate",
#                "Providing Lang"=>"C"};
#  $Function3 = {"Name"=>"fnptr",
#                "Aliased"=>"0",
#                "Return Type"=>"CCTK_REAL",
#                "Arguments"=>$ArgList3,
#                "Thorn Uses"=>"0",
#                "Thorn Provides"=>"0",
#                "Providing Fn"=>"",
#                "Providing Lang"=>""};
#  $FunctionList1 = {"SumStuff"=>$Function1,
#                    "Integrate"=>$Function2,
#                    "fnptr"=>$Function3};
#
#  my $FunctionDatabase;
#
#  $FunctionDatabase = {"CCTK_Cactus"=>$FunctionList1,
#                       "DummyThorn"=>$FunctionList1};
#
#  return $FunctionDatabase;
#
#}

#/*@@
#  @routine    RegisterAllFunctions
#  @date       Sun Feb 16 01:49:19 2003
#  @author     Ian Hawke
#  @desc
#  The routine that prints the RegisterAllFunctions.c file
#
#  This routine checks if a thorn is active and, if so, calls
#  the routine that will register the aliased functions for that
#  thorn.
#
#  As arguments it takes a reference to the FunctionDatabase.
#  It returns a list containing the C file which is written to
#  the file in the standard way.
#
#  @enddesc
#@@*/

sub RegisterAllFunctions
{
  use strict;
  my %FunctionDatabase = %{$_[0]};
  my(@data,@prototypes,@register)=();

  # generate code to register aliased functions for all active thorns
  foreach my $thorn (sort keys %FunctionDatabase)
  {
    foreach my $Function (values %{$FunctionDatabase{$thorn}})
    {
      if ($Function && $Function->{'Provided'})
      {
        push(@prototypes,"CCTK_INT Register_${thorn}(void);");
        push(@register,"  if (CCTK_IsThornActive(\"$thorn\"))");
        push(@register,'  {');
        push(@register,"    retval += Register_${thorn}();");
        push(@register,'  }');
        last;
      }
    }
  }

  # generate code to check that all aliased functions which are REQUIRED
  # by active thorns are also provided
  my @check_required_fns = ();
  foreach my $thorn (sort keys %FunctionDatabase)
  {
    my @required = ();
    foreach my $Function (values %{$FunctionDatabase{$thorn}})
    {
      if ($Function && $Function->{'Used'} == 2 && ! $Function->{'Provided'})
      {
        # go through the function database again to find all providing thorns
        my @providing_thorns = ();
        foreach my $provider (keys %FunctionDatabase)
        {
          next if ($provider eq $thorn);

          foreach my $providing_fn (values %{$FunctionDatabase{$provider}})
          {
            if ($Function->{'Name'} eq $providing_fn->{'Name'}
                && $providing_fn->{'Provided'})
            {
              push(@providing_thorns, $provider);
              last;
            }
          }
        }

        if (! @required)
        {
          push(@required,"  if (CCTK_IsThornActive(\"$thorn\"))");
          push(@required,'  {');
        }
        push(@required,"    if (! CCTK_IsFunctionAliased(\"$Function->{Name}\"))");
        push(@required,'    {');
        push(@required,'      CCTK_Warn(1, __LINE__, __FILE__, "Bindings",');
        push(@required,'                "The aliased function ' .
                       "'$Function->{'Name'}' (required by thorn '$thorn')" .
                       ' has not been provided by any active thorn !\n"');
        push(@required,'                "Please activate one of the following '.
                       'thorns which provide this function:\n"');
        push(@required,'                "  ' . join(', ', @providing_thorns) .
                       '");');
        push(@required,'      retval++;');
        push(@required,'    }');
      }
    }
    push(@required,'  }') if (@required);
    push(@check_required_fns, @required) if (@required);
  }

  push(@data, '/*@@');
  push(@data, '   @file    RegisterThornFunctions.c');
  push(@data, '   @author  Automatically generated by CreateFunctionBindings.pl');
  push(@data, '   @desc');
  push(@data, '            Register aliased functions from active thorns');
  push(@data, '   @enddesc');
  push(@data, '  @@*/');
  push(@data, '');
  push(@data, '');

  push(@data, '#include "cctk_Flesh.h"');
  push(@data, '#include "cctk_ActiveThorns.h"');
  push(@data, '#include "cctk_WarnLevel.h"');
  push(@data, '#include "cctk_Functions.h"');
  push(@data, '');

  push(@prototypes,'CCTK_INT CCTKBindings_RegisterThornFunctions(void);');
  push(@data,@prototypes);
  push(@data,'');

  push(@data,'CCTK_INT CCTKBindings_RegisterThornFunctions(void)');
  push(@data,'{');
  push(@data,'  CCTK_INT retval;');
  push(@data,'');
  push(@data,'  retval = 0;');
  push(@data,'');

  # add the registry calls
  push(@data,@register);

  # add the provide checks if there are any function aliases required
  if (@check_required_fns)
  {
    push(@data, '');
    push(@data, '  /* verify that all aliased functions which are REQUIRED');
    push(@data, '     by active thorns are also provided by someone */');
    push(@data, @check_required_fns);
  }

  push(@data, '');
  push(@data,'  return retval;');
  push(@data,'}');
  push(@data,'');

  return join ("\n",@data);
}

#/*@@
#  @routine    AliasedFunctions
#  @date       Sun Feb 16 01:52:49 2003
#  @author     Ian Hawke
#  @desc
#  The routine that creates the AliasedFunctions.c file.
#
#  This file is the heart of function aliasing. For every aliased
#  function that is USEd this file contains:
#
#  1) Function pointers for both the C and Fortran versions
#  2) The C and Fortran wrappers of the aliased function that are
#     actually linked to by the USEing thorn at compile time. These
#     just call the function pointers in (1) (if non-NULL)
#  3) A function that checks if the aliased function has been
#     PROVIDED (i.e., that the function pointer is non-NULL)
#  4) A function that overloads the function pointer in (1)
#
#  As arguments it takes the FunctionDatabase.
#  It returns a list containing the C file which is written to
#  the file in the standard way.
#
#  @enddesc
#@@*/


sub AliasedFunctions
{
  use strict;

  my %FunctionDatabase = %{$_[0]};
  my $thornFunctionList;

  my(@data)=();

  # Header Data
  push(@data, '/*@@');
  push(@data, '   @file    AliasedFunctions.c');
  push(@data, '   @author  Automatically generated by CreateFunctionBindings.pl');
  push(@data, '   @desc');
  push(@data, '            Prototypes for the aliased functions.');
  push(@data, '   @enddesc');
  push(@data, ' @@*/');
  push(@data, '');
  push(@data, '');

  push(@data, '#include <stdlib.h>');
  push(@data, '#include <string.h>');
  push(@data, "");

  push(@data, '#include "cctk_Flesh.h"');
  push(@data, '#include "cctk_WarnLevel.h"');
  push(@data, '#include "cctk_FortranString.h"');
  push(@data, "");

  my %AliasedFunctionList = {};

  foreach my $thornFunctionKey (sort keys %FunctionDatabase)
  {
    $thornFunctionList = $FunctionDatabase{$thornFunctionKey};
    if ($thornFunctionList)
    {
      my $FunctionKey;
      foreach $FunctionKey (sort keys %{$thornFunctionList})
      {
        my $Function = $thornFunctionList->{$FunctionKey};
        if ($Function)
        {
          if (($Function->{"Provided"})||($Function->{"Used"}))
          {
            $AliasedFunctionList{$FunctionKey}=$Function;
          }
        }
      }
    }
  }

  my $FunctionKey;
  foreach $FunctionKey (sort keys %AliasedFunctionList)
  {
    my $Function = $AliasedFunctionList{$FunctionKey};
    if ($Function)
    {
      $debug and print "provided Function is ",$Function->{"Name"},"\n";
      push(@data,"/*");
      push(@data," * The function pointers to be set");
      push(@data," */");
      push(@data,"");
      push(@data,&printAliasPointers("C",$Function));
      push(@data,&printAliasPointers("Fortran",$Function));
      push(@data,"");
      push(@data,"/*");
      push(@data," * The functions that are linked to by the USEing thorn");
      push(@data," */");
      push(@data,"");
      push(@data,&printAliasPrototypes("C",$Function));
      push(@data,&printAliasToWrapper("C",$Function));
      push(@data,"");
      push(@data,&printAliasPrototypes("Fortran",$Function));
      push(@data,&printAliasToWrapper("Fortran",$Function));
      push(@data,"");
      push(@data,"/*");
      push(@data," * The functions that check if it has been PROVIDEd");
      push(@data," */");
      push(@data,"");
      push(@data,&printIsAliasedPrototypes($Function));
      push(@data,&printIsAliased($Function));
      push(@data,"");
      push(@data,"/*");
      push(@data," * The functions that overload the above function pointers.");
      push(@data," */");
      push(@data,"");
      push(@data,&printRegisterAliasedPrototypes("C",$Function));
      push(@data,&printRegisterAliased("C",$Function));
      push(@data,&printRegisterAliasedPrototypes("Fortran",$Function));
      push(@data,&printRegisterAliased("Fortran",$Function));
      push(@data,"");
    }
  }
  return join ("\n",@data);
}

#/*@@
#  @routine    printAliasPointers
#  @date       Sun Feb 16 02:03:14 2003
#  @author     Ian Hawke
#  @desc
#
#  The utility routine that prints (to the list that will be printed
#  to a file) the static function pointers in AliasedFunctions.c
#
#  As arguments it takes the type of the function (i.e., whether it
#  is the prototype of the C or Fortran wrapper) and a Function as
#  defined above. It returns the string list to be printed to a file
#  in the standard way.
#
#  @enddesc
#@@*/

sub printAliasPointers
{
  use strict;

  my ($type,%Function) = ($_[0],%{$_[1]});

  my(@data)=();

  my @args = &printArgList($type,$Function{"Arguments"});
  my $rettype = $Function{"Return Type"};
  my $name;
  if ($type eq "Fortran")
  {
    $name = $Function{"Name"}."_F_Wrapper";
  }
  else
  {
    $name = $Function{"Name"}."_C_Wrapper";
  }
  push(@data,"static $rettype (*$name) (@args) = NULL;");

  return join ("\n",@data);
}

#/*@@
#  @routine    printAliasPrototypes
#  @date       Sun Feb 16 02:05:52 2003
#  @author     Ian Hawke
#  @desc
#
#  The utility routine that prints (to the list that will be printed
#  to a file) the prototypes of the wrappers to the USEd functions
#  that appear in AliasedFunctions.c
#
#  As arguments it takes the type of the function (i.e., whether it
#  is the prototype of the C or Fortran wrapper) and a Function as
#  defined above. It returns the string list to be printed to a file
#  in the standard way.
#
#  @enddesc
#@@*/

sub printAliasPrototypes
{
  use strict;

  my ($type,%Function) = ($_[0],%{$_[1]});

  my(@data)=();

  my @args = &printArgList($type,$Function{"Arguments"});
  if ( ($type eq "Fortran")&&(($Function{"Strings"})||($Function{"String pointers"})) )
  {
    @args = @{&ConvertStringArguments($Function{"Strings"}+$Function{"String pointers"},\@args)};
#    unshift @args, "$rettype ierr,";
  }
  my $rettype = $Function{"Return Type"};
  my $name;
  if ($type eq "Fortran")
  {
    $name = "CCTK_FCALL CCTK_FNAME(".$Function{"Name"}.")";
  }
  else
  {
    $name = $Function{"Name"};
  }
  push(@data,"$rettype $name (@args);");

  return join ("\n",@data);
}

#/*@@
#  @routine    ConvertStringArguments
#  @date       Sun Feb 16 02:07:06 2003
#  @author     Ian Hawke
#  @desc
#  Normally string arguments are dealt with in a very specific
#  way (because of the inter-language issues). However, at the
#  wrapper level in AliasedFunctions.c everything is just passed
#  straight through. Rather than including more special case code
#  in the printArgList routine we just use this routine to convert
#  any string arguments to the right form.
#
#  There's probably a better way of doing this.
#
#  As arguments this takes the number of strings in the ArgList
#  and the ArgList itself. It directly alters the ArgList and returns
#  it (by reference).
#
#  @enddesc
#@@*/

sub ConvertStringArguments
{
  use strict;

  my ($nstrings,@ArgList) = ($_[0],@{$_[1]});

  my $i;
  my $stringtitle;

  if ($nstrings == 1)
  {
    $stringtitle = "ONE";
  }
  elsif ($nstrings == 2)
  {
    $stringtitle = "TWO";
  }
  elsif ($nstrings == 3)
  {
    $stringtitle = "THREE";
  }

#  print $#ArgList."\n";
#  for ($i=0; $i <= $#ArgList; $i++)
#  {
#    print $ArgList[$i]."\n";
#  }

# Pop the last argument which contains ALL the string args.

  pop(@ArgList);

  push(@ArgList,"${stringtitle}_FORTSTRING_ARG");

#  print "@{ArgList} \n";

  return \@ArgList;

}

#/*@@
#  @routine    printAliasToWrapper
#  @date       Sun Feb 16 02:11:08 2003
#  @author     Ian Hawke
#  @desc
#  This prints the function that is actually linked to by the USEing thorn.
#  It calls the function pointers (if they've actually been PROVIDEd).
#
#  The arguments to this sub are the type of the USEing routine (i.e., C
#  or Fortran) and the Function as defined above. It returns the
#  string list to be printed to a file in the standard way.
#
#  @enddesc
#@@*/

sub printAliasToWrapper
{
  use strict;

  my ($type,%Function) = ($_[0],%{$_[1]});

  my(@data)=();

  my @args = &printArgList($type,$Function{"Arguments"});
  if ( ($type eq "Fortran")&&(($Function{"Strings"})||($Function{"String pointers"})) )
  {
    @args = @{&ConvertStringArguments($Function{"Strings"}+$Function{"String pointers"},\@args)};
#    unshift @args, "$rettype ierr,";
  }
  my $rettype = $Function{"Return Type"};
  my $callname;
  my $name;

  if ($type eq "Fortran")
  {
    $name = "CCTK_FCALL CCTK_FNAME(".$Function{"Name"}.")";
    $callname = $Function{"Name"}."_F_Wrapper";
  }
  else
  {
    $name = $Function{"Name"};
    $callname = $Function{"Name"}."_C_Wrapper";
  }

  push(@data,"$rettype $name (@args)");
  push(@data,"{");
#  if ($rettype != 'void' and $type != 'Fortran')
  if ($rettype ne 'void')
  {
    push(@data,"  $rettype retval;");
    push(@data,'');
  }

  my $nstrings = "";
  my $stringargs = "";
  my $nstringptrs = "";
  my $stringptrargs = "";
  my $totalstrings = "";

  if ($type eq "Fortran")
  {
    if ($Function{"Strings"} + $Function{"String pointers"} == 1)
    {
      $nstrings = "ONE";
      $stringargs = "cctki_string1";
    }
    elsif ($Function{"Strings"} + $Function{"String pointers"} == 2)
    {
      $nstrings = "TWO";
      $stringargs = "cctki_string1,cctki_string2";
    }
    elsif ($Function{"Strings"} + $Function{"String pointers"} == 3)
    {
      $nstrings = "THREE";
      $stringargs = "cctki_string1,cctki_string2,cctki_string3";
    }
    if ($Function{"String pointers"} == 1)
    {
      $nstringptrs = "ONE";
      $stringptrargs = "cctki_pstring1";
    }
    elsif ($Function{"String pointers"} == 2)
    {
      $nstringptrs = "TWO";
      $stringptrargs = "cctki_pstring1,cctki_pstring2";
    }
    elsif ($Function{"String pointers"} == 3)
    {
      $nstringptrs = "THREE";
      $stringptrargs = "cctki_pstring1,cctki_pstring2,cctki_pstring3";
    }
  }

  if ($nstrings)
  {
    push(@data,"  ${nstrings}_FORTSTRING_CREATE(${stringargs})");
    if ($Function{"String pointers"})
    {
      push(@data,"  ${nstringptrs}_FORTSTRING_PTR($stringptrargs)");
    }
    push(@data,"");
  }

  push(@data,"  if (!${callname})");
  push(@data,"  {");
  push(@data,"    CCTK_Warn(0, __LINE__, __FILE__, \"Bindings\",\"The function ${Function{\"Name\"}} has not been provided by any active thorn.\");");
  push(@data,"  }");
  push(@data,"");

  my @seq = &printCallingSequence("",$type,$type,$Function{"Arguments"});

  if ($nstrings)
  {
    my $i;
#    print $#seq."\n";
#    for ($i=0; $i<=$#seq; $i++)
#    {
#      print $seq[$i]."\n";
#    }
#   For each, pop the argument and any necessary commas.
    pop(@seq);
#    print "nstrings ".$Function{"Strings"}."\n";
    for ($i=1; $i < $Function{"Strings"}; $i++)
    {
      pop(@seq);
      pop(@seq);
    }
    if ($Function{"Strings"})
    {
      push(@seq,$stringargs);
    }
    if ($Function{"String pointers"})
    {
      push(@seq,$stringptrargs);
    }
#    print $#seq."\n";
#    for ($i=0; $i<=$#seq; $i++)
#    {
#      print $seq[$i]."\n";
#    }
    my @freeargs=split(',',$stringargs);
    if ($rettype eq 'void')
    {
      push(@data,"  (*".$callname.")(@seq);");
      foreach $stringargs (@freeargs)
      {
        push(@data,"  free(${stringargs});");
      }
    }
    else
    {
      push(@data,"  retval = (*".$callname.")(@seq);");
      foreach $stringargs (@freeargs)
      {
        push(@data,"  free(${stringargs});");
      }
      push(@data,"  return retval;");
    }
  }
  else
  {
    if ($rettype eq 'void')
    {
      push(@data,"  (*".$callname.")(@seq);");
    }
    else
    {
      push(@data,"  retval = (*".$callname.")(@seq);");
      push(@data,"  return (retval);");
    }
  }

  push(@data,"}");

  return join ("\n",@data);
}


#/*@@
#  @routine    printIsAliasedPrototypes
#  @date       Sun Feb 16 02:14:24 2003
#  @author     Ian Hawke
#  @desc
#  Prints the prototype for the routine that checks whether an
#  aliased function has been PROVIDEd.
#
#  Takes a Function as defined above and returns the string list
#  to be printed in the standard way.
#
#  @enddesc
#@@*/

sub printIsAliasedPrototypes
{
  use strict;

  my %Function = %{$_[0]};

  my(@data)=();

  my $name;
  $name = $Function{"Name"};

  push(@data,"CCTK_INT IsAliased".$name."(void);");

  return join ("\n",@data);
}

#/*@@
#  @routine    printIsAliased
#  @date       Sun Feb 16 02:18:10 2003
#  @author     Ian Hawke
#  @desc
#
#  Prints the routine that checks if an aliased function has been
#  PROVIDEd.
#
#  Takes a Function as defined above and returns the string list
#  to be printed in the standard way.

#  @enddesc
#@@*/

sub printIsAliased
{
  use strict;

  my %Function = %{$_[0]};

  my(@data)=();

  my $name;
  $name = $Function{"Name"};

  push(@data,"CCTK_INT IsAliased${name}(void)");
  push(@data,"{");
  push(@data,"  return (${name}_C_Wrapper != NULL);");
  push(@data,"}");
  return join ("\n",@data);
}

#/*@@
#  @routine    printRegisterAliasedPrototypes
#  @date       Sun Feb 16 02:19:13 2003
#  @author     Ian Hawke
#  @desc
#
#  Prints the prototypes for the routines that will set the
#  function pointers that are defined in the AliasedFunctions.c
#  file.
#
#  The arguments to this sub are the type of the USEing routine (i.e., C
#  or Fortran) and the Function as defined above. It returns the
#  string list to be printed to a file in the standard way.
#
#  @enddesc
#@@*/

sub printRegisterAliasedPrototypes
{
  use strict;

  my ($type,%Function) = ($_[0],%{$_[1]});

  my(@data)=();

  my @args = &printArgList($type,$Function{"Arguments"});

  my $name;
  if ($type eq "Fortran")
  {
    $name = $Function{"Name"}."_F";
  }
  else
  {
    $name = $Function{"Name"}."_C";
  }
  push(@data,"CCTK_INT Alias".$name."(".$Function{"Return Type"}." (*func)(@args));");

  return join ("\n",@data);
}

#/*@@
#  @routine    printRegisterAliased
#  @date       Sun Feb 16 02:20:59 2003
#  @author     Ian Hawke
#  @desc
#
#  Prints the routines that will set the function pointers
#  that are defined in the AliasedFunctions.c file.
#
#  The arguments to this sub are the type of the USEing routine (i.e., C
#  or Fortran) and the Function as defined above. It returns the
#  string list to be printed to a file in the standard way.
#
#  The function created will return 0 for success and 1 if the
#  function was already PROVIDEd by another thorn.
#
#  @enddesc
#@@*/

sub printRegisterAliased
{
  use strict;

  my ($type,%Function) = ($_[0],%{$_[1]});

  my(@data)=();

  my @args = &printArgList($type,$Function{"Arguments"});

  my $name = $Function{"Name"};
  my $fullname;
  if ($type eq "Fortran")
  {
    $fullname = $name."_F";
  }
  else
  {
    $fullname = $name."_C";
  }
  push(@data,"CCTK_INT Alias".$fullname."(".$Function{"Return Type"}." (*func)(@args))");
  push(@data,"{");
  push(@data,"  CCTK_INT aliased = ${name}_C_Wrapper != NULL;");
  push(@data,"  if (!aliased)");
  push(@data,"  {");
  push(@data,"    ".$fullname."_Wrapper = func;");
  push(@data,"  }");
  push(@data,"  return aliased;");
  push(@data,"}");

  return join ("\n",@data);
}

#/*@@
#  @routine    IsFunctionAliased
#  @date       Sun Feb 16 02:22:32 2003
#  @author     Ian Hawke
#  @desc
#  The routine that creates the IsFunctionAliased.c file.
#
#  This file just contains one function (and the Fortran wrapper)
#  that takes a string argument. If this argument is the name of
#  an aliased function then the appropriate individual function (defined
#  in AliasedFunctions.c, see above) will be called to see if it has
#  been PROVIDEd.
#
#  As arguments it takes the FunctionDatabase.
#  It returns a list containing the C file which is written to
#  the file in the standard way.
#
#  @enddesc
#@@*/

sub IsFunctionAliased
{
  use strict;

  my %FunctionDatabase = %{$_[0]};
  my $thornFunctionList;
  my $Function;

  my(@data)=();

  # Header Data
  push(@data, '/*@@');
  push(@data, '   @file    IsFunctionAliased.c');
  push(@data, '   @author  Automatically generated by CreateFunctionBindings.pl');
  push(@data, '   @desc');
  push(@data, '            The master routine to see if the aliased functions are overloaded.');
  push(@data, '   @enddesc');
  push(@data, ' @@*/');
  push(@data, '');
  push(@data, '');

  push(@data, '#include <string.h>');
  push(@data, '#include <stdlib.h>');
  push(@data, '');
  push(@data, '#include "cctk_Flesh.h"');
  push(@data, '#include "cctk_FortranString.h"');
  push(@data, '');

  # Insert function protypes:
  my %names;
  foreach $thornFunctionList (values %FunctionDatabase)
  {
    foreach $Function (values %{$thornFunctionList})
    {
      if ($Function)
      {
        if ($Function->{"Used"})
        {
          $names{$Function->{"Name"}} = undef;
        }
      }
    }
  }
  foreach my $name (sort keys %names)
  {
    push(@data, "CCTK_INT IsAliased$name(void);");
  }

  push(@data,"CCTK_INT CCTK_IsFunctionAliased(const char *function);");
  push(@data,"");
  push(@data,"CCTK_INT CCTK_IsFunctionAliased(const char *function)");
  push(@data,"{");
  push(@data,"  CCTK_INT retval = 0;");
  push(@data,"");
  push(@data,"  (void) (function + 0);");
  push(@data,"");

  my $else = "";
  foreach my $name (sort keys %names)
  {
    push(@data, "  ${else}if (! strcmp(function, \"$name\"))");
    push(@data, "  {");
    push(@data, "    retval = IsAliased".$name."();");
    push(@data, "  }");
    $else = "else ";
  }

  push(@data,"  return retval;");

  push(@data,"}");
  push(@data, '');

  push(@data, 'void CCTK_FCALL CCTK_FNAME(CCTK_IsFunctionAliased) (int *ret, ONE_FORTSTRING_ARG);');
  push(@data, 'void CCTK_FCALL CCTK_FNAME(CCTK_IsFunctionAliased) (int *ret, ONE_FORTSTRING_ARG)');
  push(@data, '{');
  push(@data, '  ONE_FORTSTRING_CREATE(name);');
  push(@data, '  *ret = CCTK_IsFunctionAliased(name);');
  push(@data, '  free(name);');
  push(@data, '}');
  push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline

  return join ("\n",@data);
}

#/*@@
#  @routine    ThornMasterIncludes
#  @date       Sun Feb 16 02:25:43 2003
#  @author     Ian Hawke
#  @desc
#  The routine that creates the ThornMasterIncludes.h file.
#
#  Just include the appropriate prototypes so that they are available
#  to any USEing thorn.
#
#  As arguments it takes the FunctionDatabase.
#  It returns a list containing the C file which is written to
#  the file in the standard way.
#
#  @enddesc
#@@*/

sub ThornMasterIncludes
{
  use strict;

  my %function_db = %{$_[0]};
  my @thorns = keys %function_db;

  my $thorn;
  foreach $thorn (sort @thorns)
  {
    my(@data) = ();

    # Header Data
    push(@data, '/*@@');
    push(@data, '   @header  $thorn/cctk_Functions.h');
    push(@data, '   @author  Automatically generated by CreateFunctionBindings.pl');
    push(@data, '   @desc');
    push(@data, '            Prototypes for overloaded functions used by all thorns');
    push(@data, '   @enddesc');
    push(@data, '  @@*/');
    push(@data, '');
    push(@data, '');

    push(@data, '#ifndef _CCTK_FUNCTIONALIASES_H_');
    push(@data, '#define _CCTK_FUNCTIONALIASES_H_ 1');
    push(@data, '');

    push(@data, '#ifdef CCODE');
    push(@data, '#ifdef __cplusplus');
    push(@data, 'extern "C" {');
    push(@data, '#endif');
    push(@data, '  CCTK_INT CCTK_IsFunctionAliased(const char *function);');
    push(@data, '#ifdef __cplusplus');
    push(@data, '}');
    push(@data, '#endif');
    push(@data, '#endif');
    push(@data, '');

    push(@data, "#include \"${thorn}_Prototypes.h\"");
    push(@data, "#define DECLARE_CCTK_FUNCTIONS DECLARE_\U$thorn\E_FUNCTIONS");

    push(@data, '#endif  /* _CCTK_FUNCTIONALIASES_H_ */');
    push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline

    my $dataout = join ("\n",@data);
    mkdir("include/${thorn}");
    &WriteFile("include/${thorn}/cctk_Functions.h",\$dataout);
    if($thorn eq "Cactus") {
      &WriteFile("include/CactusBindings/cctk_Functions.h",\$dataout);
    }
  }
}

#/*@@
#  @routine    UsesPrototypes
#  @date       Sun Feb 16 02:27:29 2003
#  @author     Ian Hawke
#  @desc
#  For every thorn that USEs a function, create the appropriate
#  prototypes. These will be included into the ThornMasterIncludes.h file.
#
#  As arguments this takes the name of the thorn and the appropriate
#  FunctionList. It returns a list containing the C file which is written to
#  the file in the standard way.
#
#  @enddesc
#@@*/

sub UsesPrototypes
{
  use strict;

  my ($thorn,%FunctionList) = ($_[0],%{$_[1]});

  my(@data) = ();

  # Header Data
  push(@data, '/*@@');
  push(@data, "   \@header  ${thorn}_Prototypes.h");
  push(@data, '   @author  Automatically generated by CreateFunctionBindings.pl');
  push(@data, '   @desc');
  push(@data, '            Prototypes for overloaded functions used by this thorn');
  push(@data, '   @enddesc');
  push(@data, '  @@*/');
  push(@data, '');
  push(@data, '');

  push(@data, "#ifndef _\U$thorn\E_PROTOTYPES_H_");
  push(@data, "#define _\U$thorn\E_PROTOTYPES_H_  1");
  push(@data, '');

  push(@data, '#ifdef CCODE');

  push(@data, '#ifdef __cplusplus');
  push(@data, 'extern "C" {');
  push(@data, '#endif');
  push(@data, '');

  $debug and print "UsesPrototypes: thorn is $thorn\n";
  foreach my $FunctionKey (sort keys %FunctionList)
  {
    my $Function = $FunctionList{$FunctionKey};
    $debug and print "  Function is ", $Function->{"Name"},"\n";
    next if (! $Function);

    if ($Function->{"Used"})
    {
      my @cargs = &printArgList("C",$Function->{"Arguments"});
      push(@data, "$Function->{\"Return Type\"} $Function->{\"Name\"}(@cargs);");
    }
    if ($Function->{"Provided"})
    {
      my @cargs = &printArgList("C",$Function->{"Arguments"});
      push(@data, "$Function->{\"Return Type\"} $Function->{\"Provider\"}(@cargs);");
    }
    push(@data, '');
  }

  push(@data, '#ifdef __cplusplus');
  push(@data, '}');
  push(@data, '#endif');

  push(@data, '#endif /* CCODE */');
  push(@data, '');

  push(@data, '#ifdef FCODE');
  push(@data, '#ifdef F90CODE');

  push(@data, "#define DECLARE_\U$thorn\E_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \\");

  foreach my $FunctionKey (sort keys %FunctionList)
  {
    my $Function = $FunctionList{$FunctionKey};
    next if (! $Function);

    next if (! $Function->{"Used"});

    my $args = $Function->{"Arguments"};
    my $line;
    push(@data, "  interface &&\\");
    if ($Function->{"Return Type"} ne 'void')
    {
      $line = "     $Function->{\"Return Type\"} function $Function->{\"Name\"}";
    }
    else
    {
      $line = "     subroutine $Function->{\"Name\"}";
    }
    $line .= " (";
    my $sep = '';
    foreach my $arg (@$args)
    {
      $line .= $sep;
      $sep = ', ';
      my $isfunctionpointer = $arg->{"Function pointer"};
      my $name = $isfunctionpointer ? $arg->{"Name"}->{"Name"} : $arg->{"Name"};
      $line .= "$name";
    }
    $line .= ")";
    $line .= " &&\\";
    push(@data, $line);
    push(@data, "       implicit none &&\\");
    foreach my $arg (@$args)
    {
      my $isfunctionpointer = $arg->{"Function pointer"};
      my $name = $isfunctionpointer ? $arg->{"Name"}->{"Name"} : $arg->{"Name"};
      my $type = $arg->{"Type"};
      my $isarray = $arg->{"Is Array"};
      my $isstring = $arg->{"String"};
      if ($isfunctionpointer)
      {
        push(@data, "       external $name &&\\");
        push(@data, "       $type $name &&\\") if ($name ne 'void');
      }
      elsif ($isstring)
      {
        push(@data, "       character(*) $name &&\\");
      }
      else
      {
        $line = "       $type $name";
        if ($isarray)
        {
          $line .= "(*)";
        }
        $line .= " &&\\";
        push(@data, $line);
      }
    }
    if ($Function->{"Return Type"} ne 'void')
    {
      push(@data, "     end function $Function->{\"Name\"} &&\\");
    }
    else
    {
      push(@data, "     end subroutine $Function->{\"Name\"} &&\\");
    }
    push(@data, "  end interface &&\\");
  }
  push(@data, '');
  push(@data, '#else /* ! F90CODE */');
  push(@data, '');
  push(@data, "#define DECLARE_\U$thorn\E_FUNCTIONS _DECLARE_CCTK_FUNCTIONS \\");

  foreach my $FunctionKey (sort keys %FunctionList)
  {
    my $Function = $FunctionList{$FunctionKey};
    next if (! $Function);

    next if (! $Function->{"Used"});

    push(@data, "  external $Function->{\"Name\"} &&\\");
    push(@data, "  $Function->{\"Return Type\"} $Function->{\"Name\"} &&\\")
      if ($Function->{"Return Type"} ne 'void');
  }
  push(@data, '');

  push(@data, '#endif /* ! F90CODE */');
  push(@data, '#endif /* FCODE */');
  push(@data, '');

  push(@data, '#endif');
  push(@data, "\n");   # workaround for perl 5.004_04 to add a trailing newline

  return join ("\n",@data);
}

#/*@@
#  @routine    ProvidedFunctions
#  @date       Sun Feb 16 02:29:13 2003
#  @author     Ian Hawke
#  @desc
#  For every thorn that PROVIDEs a function, create the wrappers
#  to the functions. These wrappers convert the arguments to the
#  correct type depending on the language of the PROVIDEing function.
#
#  A routine is there to register these providing functions by calling
#  the function defined in AliasedFunctions.c that sets the function
#  pointers.
#
#  This function takes the thorn name and the appropriate FunctionList.
#  It returns a list containing the C file which is written to
#  the file in the standard way.
#
#  @enddesc
#@@*/

sub ProvidedFunctions
{
  use strict;

  my ($thorn,%FunctionList) = ($_[0],%{$_[1]});
  my @data;
  my $Function;

  my %WrapperFunctionList;

  # Header Data
  push(@data, '/*@@');
  push(@data, "   \@file    ${thorn}_Functions.c");
  push(@data, '   @author  Automatically generated by CreateFunctionBindings.pl');
  push(@data, '   @desc');
  push(@data, "            The wrappers for functions provided by thorn ${thorn}.");
  push(@data, '   @enddesc');
  push(@data, ' @@*/');
  push(@data, '');
  push(@data, '');

  push(@data, '#include <stdio.h>');
  push(@data, '#include <stdlib.h>');
  push(@data, '#include <string.h>');
  push(@data, '');

  push(@data,"#include \"cctk_Flesh.h\"");
  push(@data, "#include \"cctk_WarnLevel.h\"");

  push(@data,"");

  foreach my $FunctionKey (sort keys %FunctionList)
  {
    $Function = $FunctionList{$FunctionKey};
    if ($Function)
    {
      if ($Function->{"Provided"})
      {
        my $rettype = $Function->{"Return Type"};
        my $providetype = $Function->{"Provider Language"};
        my @args = &printArgList($providetype,$Function->{"Arguments"});
        my $nameC = "Alias".$Function->{"Name"}."_C";
        my $nameF = "Alias".$Function->{"Name"}."_F";
        my $provider = $Function->{"Provider"};
        if ($providetype eq "Fortran")
        {
          my @fseq = &printCallingSequence($provider,$providetype,"C",
                                           $Function->{"Arguments"});
          $WrapperFunctionList{$nameF}={"Provider"=>$provider,
                                        "Wrapper Args"=>\@args,
                                        "Calling Sequence"=>\@fseq};
          my @cargs = &printArgList("C",$Function->{"Arguments"});
          my @fargs = &printArgList("Fortran",$Function->{"Arguments"});
          $WrapperFunctionList{$nameC}={
            "Provider"=>"CCTK_Wrapper_FtoC_${provider}",
            "Wrapper Args"=>\@cargs,
            "Calling Sequence"=>\@fseq};
#                print "Calling ${provider} from C using @{fseq}\n";
          my @FnWrappers = ();
          @FnWrappers = &FunctionPointerWrappers($provider,"Fortran","C",$Function->{"Arguments"});
          if (@FnWrappers)
          {
            push(@data,@FnWrappers);
            push(@data,"");
          }
          push(@data,"extern ${rettype} CCTK_FCALL CCTK_FNAME(${provider})(@{fargs});");
          push(@data,"");
          push(@data,"static ${rettype} CCTK_Wrapper_FtoC_${provider}(@{cargs});");
          push(@data,"${rettype} CCTK_Wrapper_FtoC_${provider}(@{cargs})");
          push(@data,"{");
          my @FnPtrSets = &FunctionPointerSettings($provider,
                                                   $Function->{"Arguments"});
          if (@FnPtrSets)
          {
            push(@data,@FnPtrSets);
            push(@data,"");
          }
          if ($rettype ne 'void')
          {
            push(@data,"  return CCTK_FCALL CCTK_FNAME(${provider})(@{fseq});");
          }
          else
          {
            push(@data,"  CCTK_FCALL CCTK_FNAME(${provider})(@{fseq});");
          }
          push(@data,"}");
        }
        else # providetype is C
        {
          my @cseq = &printCallingSequence($provider,$providetype,"Fortran",
                                           $Function->{"Arguments"});
          my @cargs = &printArgList("C",$Function->{"Arguments"});
          my @fargs = &printArgList("Fortran",$Function->{"Arguments"});
          $WrapperFunctionList{$nameF}={
            "Provider"=>"CCTK_Wrapper_CtoF_${provider}",
            "Wrapper Args"=>\@fargs,
            "Calling Sequence"=>\@cseq};
          $WrapperFunctionList{$nameC}={"Provider"=>$provider,
                                        "Wrapper Args"=>\@args,
                                        "Calling Sequence"=>\@cseq};
          my @FnWrappers = ();
          @FnWrappers = &FunctionPointerWrappers($provider,"C","Fortran",$Function->{"Arguments"});
          if (@FnWrappers)
          {
            push(@data,@FnWrappers);
            push(@data,"");
          }
          push(@data,"extern ${rettype} ".$provider."(@{cargs});");
          push(@data,"static ${rettype} CCTK_Wrapper_CtoF_${provider}(@{fargs});");
          push(@data,"${rettype} CCTK_Wrapper_CtoF_${provider}(@{fargs})");
          push(@data,"{");
          my @FnPtrSets = &FunctionPointerSettings($provider,
                                                   $Function->{"Arguments"});
          if (@FnPtrSets)
          {
            push(@data,@FnPtrSets);
            push(@data,"");
          }
          if ($rettype ne 'void')
          {
            push(@data,"  return (${provider})(@{cseq});");
          }
          else
          {
            push(@data,"  ${provider}(@{cseq});");
          }
          push(@data,"}");
        }
#            print $WrapperFunctionList{$nameF}{"Provider"};
#            print "\n";
#            print @{$WrapperFunctionList{$nameF}{"Calling Sequence"}};
#            print "\n";
#            print $WrapperFunctionList{$nameC}{"Provider"};
#            print "\n";
#            print @{$WrapperFunctionList{$nameC}{"Calling Sequence"}};
#            print "\n";
      }
    }
  }

  push(@data,"");
  push(@data,"CCTK_INT Register_$thorn(void);");

  # Provide prototypes for Alias<Function Name>_[CF] functions:
  foreach my $FunctionKey (sort keys %FunctionList)
  {
    $Function = $FunctionList{$FunctionKey};
    if ($Function && $Function->{"Provided"})
    {
      push(@data,&printRegisterAliasedPrototypes("C",$Function));
      push(@data,&printRegisterAliasedPrototypes("Fortran",$Function));
    }
  }

  # Definition of Register_<Thorn>:
  push(@data,"CCTK_INT Register_$thorn(void)");
  push(@data,"{");
  push(@data,"  CCTK_INT ierr;");
  push(@data,"");
  push(@data,"  ierr = 0;");
  push(@data,"");

  foreach my $FunctionKey (sort keys %FunctionList)
  {
    $Function = $FunctionList{$FunctionKey};
    if ($Function)
    {
      if ($Function->{"Provided"})
      {
        my $type = $Function->{"Provider Language"};
        my @args = &printArgList($type,$Function->{"Arguments"});
        my $nameC = "Alias".$Function->{"Name"}."_C";
        my $nameF = "Alias".$Function->{"Name"}."_F";
        my $provider = $Function->{"Provider"};

        if ($type =~ /Fortran/)
        {
          push(@data,"  ierr += $nameF(CCTK_FNAME($provider));");
          push(@data,"  ierr += $nameC(CCTK_Wrapper_FtoC_$provider);");
        }
        else
        {
          push(@data,"  ierr += $nameF(CCTK_Wrapper_CtoF_$provider);");
          push(@data,"  ierr += $nameC($provider);");
        }
        push(@data,"  if (ierr)");
        push(@data,"  {");
        push(@data,"    CCTK_Warn(0, __LINE__, __FILE__, \"Bindings\",");
        push(@data,"              \"Function already registered!\");");
        push(@data,"  }");
      }
    }
  }

  push(@data,"  return ierr;");
  push(@data,"}");
  push(@data,"");

  return join ("\n",@data);
}

#/*@@
#  @routine    FunctionPointerWrappers
#  @date       Sun Feb 16 02:33:39 2003
#  @author     Ian Hawke
#  @desc
#  If an argument is a function pointer then it has to be treated
#  a bit differently to the standard arguments.
#
#  The assumption is made that a function passed through a function
#  pointer argument is the same language as the calling function.
#  However, this need not be the same language as the calling function.
#  So a wrapper to the function pointer argument and another local
#  function pointer are created.
#
#  When the wrapper at the thorn_ProvidedFunctions level is reached
#  it will, if necessary, set the local function pointer to point to
#  the function pointer argument and call the PROVIDEing function with
#  the wrapper. This wrapper corrects the calling sequence and calls
#  the local function pointer.
#
#  The routine creates the wrappers. As arguments it has to know
#
#  1) The name of the providing function
#  2) The type (i.e., language) of the calling function
#  3) The type (i.e., language) of the providing function
#  4) The ArgumentList of the function pointer argument.
#
#  As usual, the string list to be printed is returned.
#
#  @enddesc
#@@*/

sub FunctionPointerWrappers
{
  use strict;

  my ($provide,$calltype,$providetype,@ArgList) = ($_[0],$_[1],$_[2],@{$_[3]});

  my $nfptr = 0;

  my $Arg;

  my @data = ();

  foreach $Arg (@ArgList)
  {
    if ($Arg->{"Function pointer"})
    {
      my %Function = %{$Arg->{"Name"}};
      my $Rettype = $Arg->{"Type"};
      my $WrapperName = "CCTK_Wrap".$provide.$Function{"Name"};
      my $StaticName = "(*CCTK_Fptr${provide}${Function{\"Name\"}})";
      my @callargs = &printArgList($calltype,$Function{"Arguments"});
      my @provideargs = &printArgList($providetype,$Function{"Arguments"});
      my @provideseq = &printCallingSequence($provide,$providetype,$calltype,
                                             $Function{"Arguments"});
      push(@data,"static ${Rettype} ${StaticName}(@{provideargs});");
      push(@data,"static ${Rettype} ${WrapperName}(@{callargs});");
      push(@data,"${Rettype} ${WrapperName}(@{callargs})");
      push(@data,"{");
      if ($Rettype =~ m/void/)
      {
        push(@data,"  ${StaticName}(@{provideseq});");
      }
      else
      {
        push(@data,"  return ${StaticName}(@{provideseq});");
      }
      push(@data,"}");
      $nfptr++;
    }
  }

  return join ("\n",@data);
}

#/*@@
#  @routine    FunctionPointerSettings
#  @date       Sun Feb 16 02:41:56 2003
#  @author     Ian Hawke
#  @desc
#  The routine that sets the local function pointers (as defined above).
#  This routine takes the name of the PROVIDEing function and the
#  ArgumentList, returning the string list to be printed.
#
#  Note that this routine does not create a standalone function, but is
#  intended to be called at the right place within one that does.
#
#  @enddesc
#@@*/

sub FunctionPointerSettings
{
  use strict;

  my ($provide,@ArgList) = ($_[0],@{$_[1]});

  my $Arg;

  my @data = ();

  my $nfptrs = 0;

  foreach $Arg (@ArgList)
  {
    if ($Arg->{"Function pointer"})
    {
      my $Name = $Arg->{"Name"}{"Name"};
      push(@data,"  CCTK_Fptr${provide}${Name} = ${Name};");
      $nfptrs++;
    }
  }

  return join ("\n",@data);
}

#/*@@
#  @routine    printCallingSequence
#  @date       Sun Feb 16 02:45:17 2003
#  @author     Ian Hawke
#  @desc
#  Prints the calling sequence of a given argument list.
#
#  Needs to know:
#
#  1) Whether any function pointer arguments should be passed through
#     or whether the appropriate wrapper should be passed instead.
#  2) What the type of the PROVIDEing function is
#  3) What the type of the calling function is
#  4) The ArgumentList.
#
#  As is usual it returns the string list to be printed.
#
#  @enddesc
#@@*/

sub printCallingSequence
{
  use strict;

  my ($CallFunctionWrappers,$providetype,$calltype,@ArgList) =
      ($_[0],$_[1],$_[2],@{$_[3]});

#  print "Calling sequence. Provided by ${providetype}, called by ${calltype}\n";

  my $Arg;

  my(@data)=();

#
# I think this is going to duplicate a lot of printArg without
# quite doing the same. Oh well.
#

  my $nfptr = 0;

  for (my $i=0; $i<@ArgList; $i++)
  {
    $Arg = $ArgList[$i];
    if ($Arg->{"Function pointer"})
    {
#      my $key;
#      foreach $key (sort keys %{$Arg->{"Name"}})
#      {
#        print $key." ".$Arg->{"Name"}{$key}."\n";
#      }
#      print "Here. ${CallFunctionWrappers} ".$Arg->{"Name"}{"Name"}."\n";
      if ($CallFunctionWrappers)
      {
        push(@data,"CCTK_Wrap".$CallFunctionWrappers.$Arg->{"Name"}{"Name"});
        $nfptr++;
      }
      else
      {
        push(@data,$Arg->{"Name"}{"Name"});
      }
    }
    else
    {
      my $CallArgName=&printCallArg($providetype,$calltype,$Arg);
      push(@data,$CallArgName);
    }




    if ($i < $#ArgList)
    {
      push(@data,",");
    }

  }

  return @data;
}

#/*@@
#  @routine    printCallArg
#  @date       Sun Feb 16 02:48:12 2003
#  @author     Ian Hawke
#  @desc
#  Prints the calling sequence of an individual argument.
#
#  Needs to know that type (i.e., language) of the PROVIDEing function,
#  the calling function, and the Argument.
#
#  Returns a simple scalar string containing the calling sequence.
#
#  @enddesc
#@@*/


sub printCallArg
{
  use strict;

  my ($providetype,$calltype,%Arg) = ($_[0],$_[1],%{$_[2]});

  my $prefix = '';
  if ($providetype ne $calltype and
      not ($Arg{'Is Array'} or $Arg{'String'} or $Arg{'Intent'} =~ /OUT/))
  {
    $prefix = $calltype eq 'Fortran' ? '*' : '&';
  }
  my $varname = $Arg{"Name"};
  my $data=$prefix.$varname;

#  print "$varname $providetype $calltype $Arg{\"Is Array\"} $data\n";

  return $data;

}

#/*@@
#  @routine    printArgList
#  @date       Sun Feb 16 02:49:51 2003
#  @author     Ian Hawke
#  @desc
#  Prints the argument list of a given function.
#
#  As arguments it takes the type of the function that would be
#  calling it and returns the standard string list for printing.
#
#  @enddesc
#@@*/

sub printArgList
{
  use strict;

  my ($type,@ArgList) = ($_[0],@{$_[1]});

  my $Arg;

  my(@data)=();

  my $nstrings = 0;
  my $nstringpointers = 0;
  my @stringargnames = ();
  my @stringptrargnames = ();

  for (my $i=0; $i<@ArgList; $i++)
  {
    $Arg = $ArgList[$i];
    my ($ArgType,$ArgName)=&printArg($type,$Arg);
    if ($Arg->{"String"})
    {
      if ($Arg->{"Is Array"})
      {
        $nstringpointers++;
        push(@stringptrargnames,$ArgName);
#        if ($nstringpointers>3)
#        {
#          my $message = "An aliased function contains the argument list \n\"@{ArgList}\"\nThis must not contain more than 3 CCTK_STRING:ARRAY arguments.";
#          &CST_error(1,$message,'',__LINE__,__FILE__);
#        }
      }
      else
      {
        $nstrings++;
        push(@stringargnames,$ArgName);
#        if ($nstrings>3)
#        {
#          my @Argnames=();
#          foreach $Arg (@ArgList)
#          {
#            push(@Argnames,&printArg($type,$Arg));
#          }
#          my $message = "An aliased function contains the argument list \n\"@{Argnames}\"\nThis must not contain more than 3 CCTK_STRING arguments.";
#          &CST_error(1,$message,'',__LINE__,__FILE__);
#          $nstrings=-1;
#        }
      }
    }
    else
    {
#      if ($nstrings)
#      {
#        my @Argnames=();
#        foreach $Arg (@ArgList)
#        {
#          push(@Argnames,&printArg($type,$Arg));
#          }
#        my $message = "An aliased function contains the argument list \n\"@{Argnames}\"\nThis must have all CCTK_STRING arguments at the end of the list.";
#        &CST_error(1,$message,'',__LINE__,__FILE__);
#        $nstrings=-2;
#      }
      push(@data,($ArgType,$ArgName));
      if ($i < $#ArgList)
      {
        push(@data,",\n");
      }
    }
  }

  # Check for void arg list
  if (!@ArgList)
  {
    push(@data,"void");
  }

  push(@stringargnames,@stringptrargnames);

  if ($nstrings + $nstringpointers == 1)
  {
    push(@data,"CCTK_STRING ${stringargnames[0]}");
  }
  elsif ($nstrings + $nstringpointers == 2)
  {
    push(@data,"CCTK_STRING ${stringargnames[0]}, CCTK_STRING ${stringargnames[1]}");
  }
  elsif ($nstrings + $nstringpointers == 3)
  {
    push(@data,"CCTK_STRING ${stringargnames[0]}, CCTK_STRING ${stringargnames[1]}, CCTK_STRING ${stringargnames[2]}");
  }

  return @data;
}

#/*@@
#  @routine    printArg
#  @date       Sun Feb 16 02:51:28 2003
#  @author     Ian Hawke
#  @desc
#  Prints an individual argument.
#
#  Needs to know the type (language) of the function that calls it,
#  and the argument itself. Returns the standard string list.
#  @enddesc
#@@*/

sub printArg
{
  use strict;

  my ($type,%Arg) = ($_[0],%{$_[1]});

  my(@data)=();

  my $vartype = $Arg{"Type"};
  my $suffix = "";
  my $prefix = "";

  if ( (($type eq "Fortran")&&(!$Arg{"Function pointer"})) ||
       (($type eq "C")&&(($Arg{"Is Array"})||($Arg{"Intent"}=~/OUT/))) )
  {
    &debug_print($Arg{"Name"}." needs a *");
    $suffix = "*";
  }
  if ( ($Arg{"Intent"}=~/IN/) && (!($Arg{"Intent"}=~/OUT/)) )
  {
    $prefix = "const ";
  }

  if ($Arg{"Function pointer"})
  {
# It's a FPOINTER
#    my $key;
#    foreach $key (sort keys %{$Arg{"Name"}})
#    {
#      print $key." ".$Arg{"Name"}{$key}."\n";
#    }
    push(@data,$vartype);
    my @fptrargs = &printArgList($type,$Arg{"Name"}{"Arguments"});
#    print "\n@fptrargs\n\n";
    push(@data,"(*".$prefix.$suffix.$Arg{"Name"}{"Name"}.")(@fptrargs)");
  }
  else
  {
#    print "Argument: ".$Arg{"Name"}." ".$Arg{"Is Array"}." ".$Arg{"Intent"}."\n";

    push(@data,$prefix.$vartype.$suffix);

    push(@data, $Arg{"Name"});

  }

  return @data;

}

#/*@@
#  @routine    printFunction
#  @date       Sun Feb 16 02:44:30 2003
#  @author     Ian Hawke
#  @desc
#  A debugging function that's meant to print out the useful
#  information about a Function structure. Probably obsolete by now.
#  @enddesc
#@@*/

#sub printFunction
#{
#  use strict;
#
#  my %Function = %{$_[0]};
#
#  my @data;
#
#  if ($Function{"Aliased"})
#  {
#    push(@data, "The aliased function ");
#  }
#  else
#  {
#    push(@data, "The function pointer ");
#  }
#  push(@data, $Function{"Name"}," has return type ",$Function{"Return Type"},".\n");
#  push(@data, "It is ");
#  if (!$Function{"Thorn Uses"})
#  {
#    push(@data, "not ");
#  }
#  push(@data, "used by this thorn.\n");
#  if ($Function{"Thorn Provides"})
#  {
#    push(@data, "This thorn provides ",$Function{"Name"}," with the ",
#         $Function{"Providing Lang"}," function ",$Function{"Providing Fn"},
#         ".\n");
#  }
#  push(@data, $Function{"Name"}," has the following arguments.\n");
##    push(@data, &printArgList($Function{"Arguments"}));
#
#  return join ("\n",@data);
#}

1;
