#!/usr/bin/perl -s

#use strict;
use vars qw($cctk_home $debug $verbose $h $help $document_type $outdir $directory $thornlist $thorn);

#/*@@
#  @file      InterLatex.pl
#  @date      Sun Mar  3 19:05:41 CET 2002
#  @author    Ian Kelley
#  @desc 
#  This program will take as input a thornlist, and outputs a latex table    
#  that contains the information in the thorns' interface.ccl file(s).  This latex
#  table can then be used as a stand-alone document, or as a section of a large
#  "ThornGuide"         
#  @enddesc 
#  @version 
#@@*/

#########################
# ->> InterLatex.pl <<- #
#########################################################################
#                                                                       #
#                      standard help function                           #
#                                                                       #
#########################################################################
if ($h || $help) {
   print "--> InterLatex.pl <--\n";
   print "Options:\n";
   print "\t-thornlist=     : thornlist to process\n";
   print "\t-thorn=         : (opt) specific thorn to process\n";
   print "\t-directory=     : directory of arrangements (default arrangements/)\n";
   print "\t-outdir=        : directory to dump output files, default is .\n";
   print "\t-document_type= : create a 'document' or 'section'\n";
   print "\t-cctk_home=     : root directory of cactus installation\n";
   print "\n";
   print "\t-verbose        : gives verbose output to screen\n";
   print "\t-debug          : debug mode, does not produce output files\n";
   print "\t-h|-help        : this screen\n";
   print "Example:\n";
   print "\t\$ perl -s lib/sbin/InterLatex.pl -outdir=/tmp/ -thornlist=WaveToyC.th -document_type=document \n";
   exit 0;
}

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
require "$sbin_dir/interface_parser.pl";
require "$sbin_dir/CSTUtils.pl";

# for reading of the thornlist routine: %thorns = &ReadThornlist($thornlist)
require "$sbin_dir/MakeUtils.pl";

#####################
# INITIAL VARIABLES #
#####################

my $start_directory = `pwd`;
chomp ($start_directory);

my @valid_groups = qw(private public protected);

# fixed width of table
my $width = "150mm";

# spacing between tables
my $spacing      = "3mm";

# maximum number of cells in a table before we split the table
my $cellsintable = 6;

# set some defaults
$document_type  ||= 'section';

# set some variables in ThornUtils(.pm) namespace
$ThornUtils::cctk_home          = $cctk_home;
$ThornUtils::start_directory    = $start_directory;
$ThornUtils::verbose            = $verbose;
$ThornUtils::debug              = $debug;

##################
# INITIALIZATION #
##################
my %thorns;
my %arrangements;
my %system_database;
my %interface_database;
my %arrangements_database;
my %pathsToThorns;

my @listOfThorns;

ThornUtils::CreateSystemDatabase(\%system_database);

# get/setup the output directory and the arrangements directory
$outdir                 = ThornUtils::SetupOutputDirectory($outdir);
my $arrangements_dir    = ThornUtils::GetArrangementsDir($directory);

# determine thornlist, create one if one doesn't exist
if (defined $thornlist) {
   # provided by MakeUtils.pl, returns a hash with a list of the thorns in our thornlist
   %thorns       = &ReadThornlist($thornlist);
   @listOfThorns = sort keys %thorns;
} elsif (defined $thorn) {
   @listOfThorns = $thorn;
} else {
   # we don't have a thornlist, go find all thorns in arrangements directory
   @listOfThorns = ThornUtils::CreateThornlist($arrangements_dir);
}

# returns a hash with keys as Arrangment/Thorn and values of the location of the interface.ccl file
%pathsToThorns         =  ThornUtils::GetThornPaths(\@listOfThorns, $arrangements_dir, "interface.ccl", 1);

# run the interface parser from 'interface_parser.pl'
%interface_database    = &create_interface_database(scalar(keys %system_database), %system_database, %pathsToThorns);

# parse up the output we just got into more of a tree-like format
%arrangements_database = &ReadInterfaceDatabase(\%interface_database);

if ($debug) {
   ThornUtils::Dump(\%arrangements_database);
   print "\n";
   exit 0;
} else {
   ThornUtils::ProcessAllArrangements(\%arrangements_database);
}

print "\nFinished.\n";

#########################################################################
#                 END OF MAIN SECTION OF THE PROGRAM                    # 
#########################################################################

#########################################################################
#                     BEGINNING OF SUB-ROUTINES                         #
#########################################################################

#/*@@
#  @routine   ReadInterfaceDatabase
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     Will parse up the interfaceDatabase created by create_interface_database in 
#     interface_parser.pl.   Creates (then returns) a hash of hash of hashes, etc
#     that contains all the information split up into more of a tree structure for
#     easy iteration and printing.
#     
#     Creates something in the form of:
#        $newDatabase{"CactusWave"}->{"WaveToyC"}->{"add"}->{"header"} = "something.h somethingother.h";
#
#  @enddesc 
#  @version 
#@@*/
sub ReadInterfaceDatabase 
{
   my (%interfaceDatabase)      = %{$_[0]};

   my %newDatabase;

   foreach (sort keys %interfaceDatabase)
   {
      print STDERR "\n--> [$_] = [$interfaceDatabase{$_}]" if ($debug || $verbose);

      # save he original key, as we will make it all lower case later
      # and need the original for hash keys
      my $old_key = $_;

      # just in case they have declared it, but put nothing in it.
      # if they did, don't waste our time parsing it.
      #next if (! $interfaceDatabase{$old_key} =~ /\w/);

      # drop the keys to lower-case, as they are all upper-case
      tr/A-Z/a-z/;
      /^(.*?)\/(.*?)( |$)/;

      my ($arrangement, $thorn) = ($1, $2);

      next if /^implementation/;
      next if ($interfaceDatabase{$old_key} !~ /\w/);
   

      # try to categorize things, we have a few exceptions in here we have to deal with so that hash variables
      # are not overwriting their old values.  This whole process will give us a hash in the form of:
      #     $newDatabase{"CactusWave"}->{"WaveToyC"}->{"add"}->{"header"} = "something.h somethingother.h";
      if (/^([^\s]+) ([^\s]+) ([^\s]+) ([^\s]+) ([^\s]+) ([^\s]+)$/) {
             $newDatabase{$arrangement}->{$thorn}->{$2}->{$3}->{$4}->{$5}->{$6} = $interfaceDatabase{$old_key}; 
      } elsif (/^([^\s]+) ([^\s]+) ([^\s]+) ([^\s]+) ([^\s]+)$/) {
          if ($2 eq "add" || $2 eq "uses" || $2 eq "provides") {
             $newDatabase{$arrangement}->{$thorn}->{$2}->{$4}->{$5} = $interfaceDatabase{$old_key};
          } else {
             $newDatabase{$arrangement}->{$thorn}->{$2}->{$3}->{$4}->{$5} = $interfaceDatabase{$old_key};
          }
      } elsif (/^([^\s]+) ([^\s]+) ([^\s]+) ([^\s]+)$/) {
         if ($2 eq "group") {
               $newDatabase{$arrangement}->{$thorn}->{"group details"}->{$3}->{$4} = $interfaceDatabase{$old_key};
         } else {
               $newDatabase{$arrangement}->{$thorn}->{$2}->{$3}->{$4} = $interfaceDatabase{$old_key};
         }
      } elsif (/^([^\s]+) ([^\s]+) ([^\s]+)$/) {
          $newDatabase{$arrangement}->{$thorn}->{$2}->{$3} = $interfaceDatabase{$old_key};
      } elsif (/^([^\s]+) ([^\s]+)$/) { 
          $newDatabase{$arrangement}->{$thorn}->{$2} = $interfaceDatabase{$old_key};
      } else {
         print "\n--> Confused: [$_] = [$interfaceDatabase{$_}]";
      }
 
   }

   return %newDatabase;
}

#/*@@
#  @routine   ProcessOneThorn
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
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
   my $thorn        = $_[2];

   # open up the output file
   my $ofh = ThornUtils::StartDocument("inter", $thorn, $outdir, $arrangement, "Interfaces", $document_type);
   
   # print out the table   
   &LatexTableElement(\%thorn);

   # close the output file
   ThornUtils::EndDocument($ofh, $document_type);
}


#/*@@
#  @routine   LatexTableElement
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#    Takes whatever table element is currently reffered to by $table    
#    and prints it out into a LaTeX table.  Only nifty things it curr.  
#    does is NOT print ranges for BOOLEAN and SHARED elements.          
#
#  @enddesc 
#  @version 
#@@*/
sub LatexTableElement 
{
   # get the stuff passed in, first is a hash of the given variable, second is the name of the variable,
   # and third is the to-date longest variable in the thorn (for latex column formatting)
   my %thorn      = %{$_[0]};

   # get the different properties of the variable, and clean then up so we can output good latex
   my $variable_name = ThornUtils::CleanForLatex($_[1]);
   my @temp; 

   print "\n\\vspace\{$spacing\} \\subsection\*\{General\}";

   # print out implementations and inheriting
   PrintVar("Implements", $thorn{"implements"});
   PrintVar("Inherits",   $thorn{"inherits"}) if ($thorn{"inherits"} !~ /^$/);


   # now we are going to go through and print out stuff for each group type (public, private, protected, etc)
   my $printgridtitle = 1;
   foreach my $group_scope (@valid_groups) 
   {
      next if (! defined $thorn{$group_scope});

      my @groups = split/\s+/, $thorn{$group_scope}->{"groups"};

      next if (@groups < 1);

      print "\n\\subsection\*\{Grid Variables\}" if ($printgridtitle);

      my $gs_scope = $group_scope;
      $gs_scope    =~ tr/a-z/A-Z/;

      print "\n\\vspace\{5mm\}\\subsubsection\{$gs_scope GROUPS\}\n";
      print "\n\\vspace\{5mm\}\n\n\\begin\{tabular*\}\{$width\}\{|c|c\@\{\\extracolsep\{\\fill\}\}|rl|\} \\hline \n";
      print "~ \{\\bf Group Names\} ~ & ~ \{\\bf Variable Names\} ~  &\{\\bf Details\} ~ & ~\\\\ \n";
      print "\\hline \n";

      $printgridtitle = 0;
      my $counter     = 0;

      # where @groups is all the different groups for this thorn in this particular storage type (e.g private)
      foreach my $group (@groups) 
      {
         $group =~ tr/A-Z/a-z/;
         next if ($thorn{"group"}->{$group} =~ /^$/);
  
         # so we are splitting the table every so often, otherwise it could overrun off the page
         if (( ! ($counter % $cellsintable)) && ($counter ne 0) )
         {
            print "\\end\{tabular*\} \n\n";
            print "\n\n\\vspace\{5mm\}";
            print "\n\\vspace\{5mm\}\n\n\\begin\{tabular*\}\{$width\}\{|c|c\@\{\\extracolsep\{\\fill\}\}|rl|\} \\hline \n";
            print "~ \{\\bf Group Names\} ~ & ~ \{\\bf Variable Names\} ~  &\{\\bf Details\} ~ & ~ \\\\ \n";
            print "\\hline \n";
         }

         # print group name
         print ThornUtils::ToLower(ThornUtils::CleanForLatex($group));
         my $firstpass = 1;

         my @group_variables = split/\s+/, ThornUtils::CleanForLatex($thorn{"group"}{$group});

         $counter++;

         # now we are just dealing with one group, so we are going to go print out the details for it.
         my $var_counter = 0;
         foreach my $group_detail (sort keys %{$thorn{"group details"}->{$group}}) 
         {
            my $value = ThornUtils::CleanForLatex(ThornUtils::CleanFromC($thorn{"group details"}->{$group}->{$group_detail}));

            # print nothign as we are dealing with the same group
            if (! $firstpass) {
               print "~ &"; 
               $firstpass = 0;
            }

            my $new_counter = 1;

            # split things (group properties) onto different lines if need be
            if (@temp = split/,/, $value) 
            {
               foreach my $val (@temp) 
               {
                  if ($new_counter--) {
                     print " & $group_variables[$var_counter] & " . ExpandGroupName($group_detail) . " & $val \\\\ \n";  
                  } else {
                     print "& ~ & " .  ExpandGroupName($group_detail) ." & $val \\\\ \n";
                  }
               }
            } else {
               print "$group_variables[$var_counter] & " . ExpandGroupName($group_detail) . " & $value \\\\ \n"; 
            }

            $var_counter++;
         } #foreach %{thorn...}

         print "\\hline \n";
      } #foreach @groups

      print "\\end\{tabular*\} \n\n";
      delete $thorn{$group_scope};
   } # foreach @valid_groups

   print "\n\n\\vspace\{5mm\}";

   # now we are going to print out some extra information that is general to the thorn
   PrintHeaderOrFunction("Adds Header", $thorn{"add"}->{"header"},        $thorn{"add"});
   PrintHeaderOrFunction("Uses Header", $thorn{"uses"}->{"header"},       $thorn{"uses"});
   PrintHeaderOrFunction("Provides",    $thorn{"provides"}->{"function"}, $thorn{"provides"});
} 

#/*@@
#  @routine   PrintVar
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     Used to  print out a variable by spitting it if we need to, used for printing stuff like
#     what a group "inherits" or "implements"
#
#  @enddesc 
#  @version 
#@@*/
sub PrintVar 
{
   my $title = shift;
   my $var   = shift;

   my %alias = %{$_[0]};
   
   my @temp;

   return if ($var !~ /\w/);

   print "\n\n\\noindent \{\\bf ". ThornUtils::Translate(ThornUtils::CleanForLatex($title)) ."\}: ";

   if (@temp = split/\s+/, $var)
   {
      foreach (@temp) {
         print "\n\n" . ThornUtils::ToLower(ThornUtils::CleanForLatex($_));
      }
   } else {
      print "\n\n" . ThornUtils::ToLower(ThornUtils::CleanForLatex($var));
   }
   print "\n\\vspace\{2mm\}";
}

#/*@@
#  @routine   PrintHeaderOrFunction
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     Very similar to &PrintVar, but does not lowercase stuff and deals with optional aliasing of what
#     it is including.  This is used to print out headers/functions/etc it addes, uses or provides. 
#
#  @enddesc 
#  @version 
#@@*/
sub PrintHeaderOrFunction
{
   my $title = shift;
   my $var   = shift;

   my %alias = %{$_[0]};
   
   my @temp;

   return if ($var !~ /\w/);

   print "\n\n\\noindent \{\\bf ". ThornUtils::Translate(ThornUtils::CleanForLatex($title)) ."\}: ";

   if (@temp = split/\s+/, $var)
   {
      foreach (@temp) 
      {
         my $lc_var = $_;
         $lc_var =~ tr/A-Z/a-z/;

         print "\n\n" . ThornUtils::CleanForLatex($_);

         if ( (defined $alias{$lc_var}) && ($alias{$lc_var}->{"to"} ne $_) ) {
            print " to " . ThornUtils::CleanForLatex($alias{$lc_var}->{"to"});
         }
      }
   } else {
      my $lc_var = $var;
      $lc_var =~ tr/A-Z/a-z/;
      
      print "\n\n" . ThornUtils::CleanForLatex($var);

      if ( (defined $alias{$lc_var}) && ($alias{$lc_var}->{"to"} ne $var) ) {
         print " to " . ThornUtils::CleanForLatex($alias{$lc_var}->{"to"});
      }
   }
   print "\n\\vspace\{2mm\}";
}

#/*@@
#  @routine   ExpandGroupName
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     When we are printing out the properties of a group, sometimes we want to expand the 
#     names of stuff out to be more descriptive.  This function provides that.
#
#  @enddesc 
#  @version 
#@@*/
sub ExpandGroupName 
{
   my $name = shift;

   if ($name eq "gtype") {
       return "group type";
    } elsif ($name eq "dim") { 
       return "dimensions";
    } elsif ($name eq "distrib") { 
       return "distribution";
    } elsif ($name eq "vtype") { 
       return "variable type";
    } elsif ($name eq "timelevels") { 
       return "timelevels";
    } else {
       return ThornUtils::ToLower(ThornUtils::CleanForLatex($name));
    }
}
