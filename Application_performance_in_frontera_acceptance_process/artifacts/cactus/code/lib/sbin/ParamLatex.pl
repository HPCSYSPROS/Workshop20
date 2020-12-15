#!/usr/bin/perl -s

# here we say to play it safe and not let variables pop up from anywhere, except for 
# the following variables in the 'use vars' statement, that can be passed in from 
# the command.
use strict;
use vars qw($h $help $cctk_home $thornlist $outdir $verbose $debug $arrangements_dir $document_type $thorn);

#########################
# ->> ParamLatex.pl <<- #
#########################################################################
#                      standard help function                           #
#########################################################################
if ($h || $help) 
{
print <<EOC;
This program will take as input a thornlist, and outputs a latex table that contains the information in the thorns' param.ccl file(s).  This latex table can then be used as a stand-alone document, or as a section of a large "ThornGuide"         

	-cctk_home=        : root directory of Cactus 
	-arrangements_dir= : arrangements directory
	-thornlist=        : thornlist to process
	-thorn=            : specific thorn to process
	-outdir=           : where to place resulting output files
	-document_type=    : is this a self containted 'document' or 'section'

	-verbose           : verbose output
	-debug             : (=1 is debug, >1 prints more info)
	-h | -help         : this screen

Example:
	\$ perl -s ParamLatex.pl -cctk_home=/tmp/CactusCheckout/ -thornlist=/tmp/WaveToy.th -outdir=/tmp/output/
EOC
   
exit 0;
}

#/*@@
#  @file      ParamLatex.pl
#  @date      Sun Mar  3 19:05:41 CET 2002
#  @author    Ian Kelley
#  @desc 
#  This program will take as input a thornlist, and outputs a latex table    
#  that contains the information in the thorns' param.ccl file(s).  This latex
#  table can then be used as a stand-alone document, or as a section of a large
#  "ThornGuide"         
#  @enddesc 
#  @version 
#@@*/

# setup the cctk_home, if it doesn't exist, we leave it blank
$cctk_home  .= '/' if (($cctk_home !~ /\/$/) && (defined $cctk_home));

# set up the sbin dir, tacking cctk_home on the front
my $sbin_dir = "${cctk_home}lib/sbin";

##############
# REQUIRE(S) #
##############

# common procedures used to create the thornguide(s)
require "$sbin_dir/ThornUtils.pm";

# for use of create_parametere_database
require "$sbin_dir/parameter_parser.pl";
require "$sbin_dir/CSTUtils.pl";

# for reading of the thornlist routine: %thorns = &ReadThornlist($thornlist)
require "$sbin_dir/MakeUtils.pl";

#####################
# INITIAL VARIABLES #
#####################

my $start_directory = `pwd`;
chomp ($start_directory);

# table width in ''mm''
my $TABLE_WIDTH     = "160";

# maximum characters that can be in a variable (range) name until we split it out of the table
my $MAX_VAR_LENGTH  = "20";

# set some variables in ThornUtils(.pm) namespace
$ThornUtils::cctk_home       	= $cctk_home;
$ThornUtils::start_directory	= $start_directory;
$ThornUtils::verbose 		= $verbose;
$ThornUtils::debug		= $debug;

# Just declare stuff in a common place, with some brief descriptions

my %thorns;                # what we add the output of &ReadParameterDatabase to, used for
                           #  /-> used temporary when calling &ReadThornlist, which returns a hash
my %arrangements_database; # |--> simplicity of organization.
my %parameter_database;    # what we get as a return value from &create_parameter_database
my %pathsToThorns;         # hash of thorn names and absolute paths to their directories
my %arrangements;          # hash of arrangements, with links to lists of thorns
my @listOfThorns;          # a list of thorns, made by reading/creating a thornlist

# set some defaults
$document_type 	||= 'section';

#####################
# END: DECLARATIONS #
#####################

# get/setup the output directory and the arrangements directory
$outdir 		= ThornUtils::SetupOutputDirectory($outdir); 
$arrangements_dir 	= ThornUtils::GetArrangementsDir($arrangements_dir);

# determine thornlist, create one if one doesn't exist
if (defined $thornlist) {
   # provided by MakeUtils.pl, returns a hash with a list of the thorns in our thornlist
   %thorns 	 = &ReadThornlist($thornlist);
   @listOfThorns = sort keys %thorns; 
} elsif (defined $thorn) {
   @listOfThorns = $thorn;
} else {
   # we don't have a thornlist, go find all thorns in arrangements directory
   @listOfThorns = ThornUtils::CreateThornlist($arrangements_dir);
}

# this will return us a hash with keys as thorn names, and values as absolute paths to the
# thorn's directory param.ccl can be located in that path.
#   We need this information to create a parameter database using create_parameter_database
#
# We are not doing ''one'' call to parameter_database as we easily could, because there is NO WAY
# (currently) to distinguish between two identical thorns in different arrangements.  So we
# would get stuff from Alphathorns/IOHDF5 in CactusBase/IOHDF5, or/and visa-versa. 
ThornUtils::ClassifyThorns(\%arrangements, @listOfThorns);

# lets go through, one arrangement at a time
foreach my $arrangement (sort keys %arrangements) 
{
   print "\n$arrangement" if ($debug);

   # now each thorn in the given arrangement
   foreach my $thorn (@{$arrangements{$arrangement}}) 
   {
      print "\n\t$thorn" if ($debug);

      # get the path for this individual thorn
      %pathsToThorns = ThornUtils::GetThornPaths(["$arrangement/$thorn"], $arrangements_dir, "param.ccl");

      # now we create a parameter database (for one thorn), using the function 
      # &create_parameter_database which is provided by 'parameter_parser.pl'
      %parameter_database = &create_parameter_database(%pathsToThorns);

      # Go split up the information we have recieved (again, for just this one thorn)
      $arrangements_database{$arrangement}->{$thorn} = &ReadParameterDatabase(\%parameter_database);
   }
}

# just dump out the data-structure if we are in debug mode, don't create any files
if ($debug) {
   ThornUtils::Dump(\%arrangements_database);
   print "\n";
   exit 0;
} else {
   ThornUtils::ProcessAllArrangements(\%arrangements_database);
}

print "\nFinished.\n";
#############
# END: MAIN #
#############

#########################################################################
#                     BEGINNING OF SUB-ROUTINES                         #
#########################################################################

#/*@@
#  @routine   ReadParameterDatabase
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#  Hops through the parameter_database we got from &create_paramater_database, and trys 
#  to classify things in a complex data-structure (hash of hashes of (hashes/arrays)) or
#  something like that. 
#     &ReadParameterDatabase(\%parameter_database); 
#
#  It then returns this data-structure, which can then be better (and cleaner) accessed
#  when creating the latex page (as all the parsing was done here)
#  
#  Below are some examples of how to access data in the newly created structure:
#
#  [print out all the ''variables'']
#     print "\nProgram variables:";
#     foreach my $variable (sort keys %{$thorn{"variables"}}) {
#        print "\n   $variable";
#        foreach my $key (sort keys %{$thorn{"variables"}{$variable}}) {
#           print "\n      $key = $thorn{\"variables\"}->{$variable}->{$key}";
#        }
#     }
#
#  [print out which variables are shared (an array)]
#     print "\nShared group variables: ";
#     foreach (@{$thorn{"groups"}{"shares"}}) {
#        print "\n   $_"; 
#     }
#
#  NOTE: Naturally, you will have to take into account that the returned hash (%thorn) may be
#        part of a greater data-structure, so you may have to do something like:
#           foreach (sort keys %{$thorns{"CactusWave"}{"WaveToyC"}{"variables"}}) {
#  @enddesc 
#  @version 
#@@*/
sub ReadParameterDatabase
{
   my (%parameterDatabase)	= %{$_[0]};

   my %thorn;

   # only deals with one thorn at a time, as that is currently the nature of this program
   foreach (sort keys %parameterDatabase)
   {
      print "\n--> [$_] = [$parameterDatabase{$_}]" if ($debug || $verbose);

      # save he original key, as we will make it all lower case later
      # and need the original for hash keys
      my $old_key = $_;

      # just in case they have declared it, but put nothing in it.
      # if they did, don't waste our time parsing it.
      next if (! $parameterDatabase{$old_key} =~ /\w/);

      # drop the keys to lower-case, as they are all upper-case
      tr/A-Z/a-z/;       
   
      # see if we are grabbing variable lists 
      # (by scopes: global, private, restricted, shared, ...)
      if (/^(.*?)\s+(.*?)\s+(.*?)\s*?variables$/) 
      {
         my @vars = split/\s+/, ThornUtils::ToLower($parameterDatabase{$old_key});

         # remember the scopes for each variable
         foreach (@vars) { 
            $thorn{"variables"}->{$_}->{"scope"} = $2 eq "shares" ? "shared from ". ThornUtils::ToUpper($3) : $2;
         }

         @{$thorn{"groups"}->{$2}} = @vars;
         next;
      }

      # see if we can parse out some information from the keys
      # if we cannot, drop a warning and continue
      if (! /^(.*?)\s(.*?)\s(.*)$/) {
         print STDERR "\nConfused: [$_] = [$parameterDatabase{$old_key}]" if ($verbose); 
         next;
      }

      # let the thorn hash explicitly know the thorn name
      $thorn{"thorn name"} = $1;

      # add the new variable in the form of:
      #   $thorn{"variables"}->{"initialize_memory"}->{"default"} = "none";
      $thorn{"variables"}->{$2}->{$3} = $parameterDatabase{$old_key};
   } # end: foreach %parameterDatabase

   # return a reference to the newly created hash
   return \%thorn;

}

#/*@@
#  @routine   ProcessOneThorn
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     Opens/closes output file and calls the function to  'latex table' for each given thorn.
#
#     This function is called by ThornUtils::ProcessOneArrangement, the main section of this program
#     calls ThornUtils::ProcessAllArrangements which calls ProcessOneArrangement.  Things are simply 
#     split up this way so that if later things ever wanted to be changed around more, then can be.
#     This function is called by ThornUtils::ProcessOneArrangement, the main section of this program
#     calls ThornUtils::ProcessAllArrangements which calls ProcessOneArrangement.  Things are simply 
#     split up this way so that if later things ever wanted to be changed around more, then can be.
#
#  @enddesc 
#  @version 
#@@*/
sub ProcessOneThorn 
{
   # get the thorn hash
   my (%thorn)      = %{$_[0]};
   my $arrangement  = $_[1];
   my $thorn 	    = $_[2];
   my $latex_output = "";
  
   my $ofh = ThornUtils::StartDocument("param", $thorn, $outdir, $arrangement, "Parameters", $document_type);

   # go get the latex for setting the width of the columns based on the longest variable we have 
   # (contained in the vaor $paragraph_len
   if (defined $thorn{"variables"}) {
      $latex_output = ThornUtils::SetWidth($TABLE_WIDTH, ThornUtils::CleanForLatex(FindMaxVarLen($thorn{"variables"})));
   }

   # go through each group, then find all variables in the specific group type and generate the tables calling
   # the function &CreateLatexTable 
   foreach my $group_type (sort keys %{$thorn{"groups"}}) 
   {
      foreach my $variable_name (sort @{$thorn{"groups"}{$group_type}}) 
      {
         # send in the 'variable' hash, the name of the variable, and the longest variable in this thorn
         $latex_output .= &CreateLatexTable($thorn{"variables"}{$variable_name}, $variable_name);
      }
   }

   print $latex_output;
   ThornUtils::EndDocument($ofh, $document_type);
}

#/*@@
#  @routine   CreateLatexTable
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     Creates a latex table for a given variable, returns the latex code to output
#        &CreateLatexTable(<variable hash>, <variable name>);
#
#  @enddesc 
#  @version 
#@@*/
sub CreateLatexTable 
{
   # get the stuff passed in, first is a hash of the given variable, second is the name of the variable,
   # and third is the to-date longest variable in the thorn (for latex column formatting)
   my %variable      = %{$_[0]};

   # get the different properties of the variable, and clean then up so we can output good latex
   my $variable_name = ThornUtils::CleanForLatex($_[1]);
   my $scope         = ThornUtils::CleanForLatex($variable{"scope"});
   my $default       = ThornUtils::RemoveQuotes(ThornUtils::CleanForLatex($variable{"default"}));
   my $description   = ThornUtils::RemoveQuotes(ThornUtils::CleanForLatex($variable{"description"}));

   # set some vars to hold the output we create
   my $latex_output  = "";
   my $range_output  = "";
   my $extra_content = "";

   # go through each range, check to see if we can find a longer variable for latex formatting,
   # also compile the latex output for the ranges
   for (my $i = 1; $i <= $variable{"ranges"}; $i++) 
   {
      my $over_range = 1;
      my $range      = ThornUtils::RemoveQuotes(ThornUtils::CleanForLatex($variable{"range $i range"}));
      my $range_desc = ThornUtils::RemoveQuotes(ThornUtils::CleanForLatex(ThornUtils::ChopVariable($variable{"range $i description"}, $MAX_VAR_LENGTH)));

      # generate latex output for this range
      if (length($range) > $MAX_VAR_LENGTH) 
      {
         $range_output .= "\\multicolumn{1}{|p{\\maxVarWidth}|}{see [$over_range] below} & \\multicolumn{2}{p{\\paraWidth}|}{$range_desc} \\\\";
         $extra_content .= "\\noindent {\\bf [$over_range]} \\noindent \\begin{verbatim}". ThornUtils::ChopVariable($range, ($TABLE_WIDTH-10)/2)."\\end{verbatim}";

         $over_range++;
      } else {
         $range_output .= "\\multicolumn{1}{|p{\\maxVarWidth}|}{\\centering $range} & \\multicolumn{2}{p{\\paraWidth}|}{$range_desc} \\\\";
      }
   }

   # start the table (tabular) enviroment, printing out the 'name', 'scope' and 'type' of the variable
   $latex_output .= "\\noindent \\begin{tabular*}{\\tableWidth}{|c|l\@{\\extracolsep{\\fill}}r|}\n\\hline\n\\multicolumn{1}{|p{\\maxVarWidth}}{$variable_name} & {\\bf Scope:} $scope & $variable{\"type\"} \\\\";

   # print out the ranges, do nothing for shared variables
   if ($scope =~ /shared from/) {
      if ($variable{"ranges"} >= 1) {
         $latex_output .= "\\hline\n\\multicolumn{3}{|l|}{\\bf Extends ranges:}\\\\ \n\\hline";
      }
   } else {
      $latex_output .= "\\hline\n\\multicolumn{3}{|p{\\descWidth}|}{{\\bf Description:}   {\\em $description}} \\\\\n\\hline";

      if ($default !~ /\w/) {
         $default = "(none)";
      }

      # print out the range headings if we have any ranges, otherwise just print a box for the default value
      if ($variable{"ranges"} >= 1) {
         $latex_output .= "{\\bf Range} & &  {\\bf Default:} $default \\\\";
      } else {
         $latex_output .= " & & {\\bf Default:} $default \\\\";
      }
   }

   # insert the range data compiled earlier in a foreach loop
   $latex_output .= $range_output;

   # end the table
   $latex_output .= "\\hline\n\\end{tabular*}\n\n\\vspace{0.5cm}";

   # return the output, tacking on the $extra_content, which will be any variable ranges that exceeded 
   # the MAX_VAR_LENGTH restriction (for fitting into the boxes correctly)
   return $latex_output . $extra_content;
}

#/*@@
#  @routine   FindMaxVarLen
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     Finds the maximum length of all the ranges and variable names, to be used for setting the
#     length of the different table cells.
#        &FindMaxVarLen(<thorn hash>);
#
#  @enddesc 
#  @version 
#@@*/
sub FindMaxVarLen 
{
   my %thorn   = %{$_[0]};
   my $max_len = "";

   # we are going to go through each variable name and range name to see where
   # the largest text is, then we will use this for later formatting of our 
   # latex tables (so we do not get paragraph run-off
   foreach my $variable (sort keys %thorn)
   {
      # we will always take the variable length as the standard maximum length,
      # regardless of if it is longer than MAX_VAR_LENGTH
      $max_len = length($variable) > length($max_len) ? $variable : $max_len;

      for (my $i = 1; $i <= $thorn{$variable}{"ranges"}; $i++)
      {
         my $range = $thorn{$variable}{"range $i range"};

         $range =~ s/^\s*(.*?)\s*$/$1/;

         # if this new length is greater than the last one, but less than the max allowed
         # length, then we assign the maximum variable length to this new length
         if ((length($range) > length($max_len)) && (length($range) < $MAX_VAR_LENGTH)) {
            $max_len = $range;
         }
      }
   }

   return $max_len;
}
