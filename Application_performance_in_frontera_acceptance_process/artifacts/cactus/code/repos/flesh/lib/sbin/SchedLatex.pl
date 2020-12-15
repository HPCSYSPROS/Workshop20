#!/usr/bin/perl -s

use strict;
use vars qw($h $help $cctk_home $thornlist $outdir $verbose $debug $directory $document_type $thorn);

#/*@@
#  @file      SchedLatex.pl
#  @date      Sun Mar  3 19:05:41 CET 2002
#  @author    Ian Kelley
#  @desc 
#  This program will take as input a thornlist, and outputs a latex table    
#  that contains the information in the thorns' schedule.ccl file(s).  This latex
#  table can then be used as a stand-alone document, or as a section of a large
#  "ThornGuide"         
#  @enddesc 
#  @version 
#@@*/

#########################
# ->> SchedLatex.pl <<- #
#########################################################################
#                      Standard HELP Function                           #
#########################################################################
if ($h || $help) {
   print "--> SchedLatex.pl <--\n";
   print "   This program will take as input a thornlist, and outputs a latex table that contains the information in the thorns' schedule.ccl file(s).  This latex table can then be used as a stand-alone document, or as a section of a larger 'ThornGuide'";
   print "Options:\n";
   print "\t-thornlist=     : (opt) list specific thorns to process\n";
   print "\t-thorn=        : (opt) specific thorn to process\n";
   print "\t-directory=     : (opt) dir. of arrangements (default arrangements/)\n";
   print "\t-outdir=        : (opt) directory to dump output files, default is .\n";
   print "\t-document_type= : (opt) 'document' or 'section'\n";
   print "\n";
   print "\t-debug          : (opt) prints thorn name on each schedule \n";
   print "\t-verbose        : (opt) gives verbose output to screen\n";
   print "\t-h/-help        : (opt) this screen\n";
   print "\nExample:\n";
   print "\t\$ perl -s /lib/sbin/SchedLatex.pl -thornlist=WaveToyC.th -outdir=/tmp/ -document_type=document\n";
   exit 0;
}

#########################################################################
#                                                                       #
#  This program will take as input a thorn, and output a latex table    #
#  that contains the information in that thorns schedule.ccl file.      #
#                                                                       #
#########################################################################

##############
# REQUIRE(S) #
##############
$cctk_home .= '/' if (($cctk_home !~ /\/$/) && defined $cctk_home);

my $sbin_dir = "${cctk_home}lib/sbin";
require "$sbin_dir/ScheduleParser.pl";
require "$sbin_dir/CSTUtils.pl";

# common procedures used to create the thornguide(s)
require "$sbin_dir/ThornUtils.pm";

# for reading of the thornlist routine: %thorns = &ReadThornlist($thornlist)
require "$sbin_dir/MakeUtils.pl";

##################
# INITIALIZATION #
##################
my $TABLE_WIDTH   ||= "160";
$document_type    ||= "section";

my $start_directory = `pwd`;
chomp ($start_directory);

# set some variables in ThornUtils(.pm) namespace
$ThornUtils::cctk_home          = $cctk_home;
$ThornUtils::start_directory    = $start_directory;
$ThornUtils::verbose            = $verbose;
$ThornUtils::debug              = $debug;

my $ofh;
my %arrangements;
my %thorns;
my %schedule_database;
my %pathsToThorns;
my %arrangements_database;

my @listOfThorns;

my %var_mapping = (
	'STOR' => 'Storage',
	'LANG' => 'Language',
	'AFTER'=> 'After',
	'TRIG' => 'Triggers',
	'SYNC' => 'Sync',
);

# get/setup the output directory and the arrangements directory
$outdir                 = ThornUtils::SetupOutputDirectory($outdir);
my $arrangements_dir       = ThornUtils::GetArrangementsDir($directory);

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

# this will return us a hash with keys as thorn names, and values as absolute paths to the
# thorn's directory param.ccl can be located in that path.
#   We need this information to create a schedule database using create_schedule_database
#
# We are not doing ''one'' call to schedule_database as we easily could, because there is NO WAY
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
      %pathsToThorns = ThornUtils::GetThornPaths(["$arrangement/$thorn"], $arrangements_dir, "schedule.ccl");

      # we are selecting STDOUT so that any junk from the &create_schedule_database won't get in our output
      my $filehandle = select(STDOUT);
      %schedule_database = &create_schedule_database(%pathsToThorns);
      select($filehandle);

      $arrangements_database{$arrangement}->{$thorn} = &ReadScheduleDatabase(\%schedule_database);
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

#########################################################################
#                 END OF MAIN SECTION OF THE PROGRAM                    # 
#########################################################################

#########################################################################
#                     BEGINNING OF SUB-ROUTINES                         #
#########################################################################

#/*@@
#  @routine   ReadLatexDatabase
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#   Calls schedule_parser.pl, which will read in the schedule.ccl file  
#   from a single thorn and return all data as a %hash table, which we  
#   will then parse to put into our own %hash tables named according to 
#   the variable names.                                                 
#                                                                       
#   %(variable_name)  : any number of hashes created with their names   
#                       being the variable names, they then have $keys  
#                       (descriptions) with $values (well, values)      
#                          (e.g.) $name{"default"} = "Cactus";          
#  @enddesc 
#  @version 
#@@*/
sub ReadScheduleDatabase
{
  my (%schedule_database) = %{$_[0]};
  my ($name, $block, $var, $conditionals) = "";

  my %thorn;

  foreach (sort keys %schedule_database)
  {
     print STDERR "\n[$_] --> [$schedule_database{$_}]" if ($debug || $verbose);

     if (/^(.*?)\s(.*?)\s(.*)$/) 
     {
        ($name, $block, $var) = ($1, $2, $3); 
        chomp($schedule_database{$_});

        $thorn{$block}->{$var} = $schedule_database{$_};

        $thorn{$block}->{"CONDITIONAL"} = "0";
        $thorn{$block}->{"THORN"}       = $name;
     } else {
        print STDERR "\n\"$_ --> $schedule_database{$_}\"" if ($verbose || $verbose);
     }
   } #-- foreach %schedule_database

   # conditional blocks/statements
   $conditionals = $schedule_database{"$name FILE"};
   foreach my $key (sort keys %thorn) 
   {
      next if ($key !~ /(BLOCK|STATEMENT)\_(\d+)/);

      my $b = "\@$1\@$2";

      # try to figure out if stuff is conditional storage
      while ($conditionals =~ /\b(if|else)\b.*?\{\s*?(.*?(?:.*?\{.*?[^\{].*?\}.*?)?(?:.*?\{.*?[^\{].*?\}.*?)?(?:.*?\{.*?[^\{].*?\}.*?)?(?:.*?\{.*?[^\{].*?\}.*?)?(?:.*?\{.*?[^\{].*?\}.*?)?.*?)\}/imgs) 
      {
         my $if = $2;
         if ($if =~ /\Q$b\E/) {
            $thorn{$key}->{"CONDITIONAL"} = 1;
         }
      }
   }   

   return \%thorn;
}

#/*@@
#  @routine   ProcessOneThorn
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
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
 
   my $ofh = ThornUtils::StartDocument("schedule", $thorn, $outdir, $arrangement, "Schedule", $document_type);

   # go print out all the good stuff for any given thorn
   &CreateLatexTable(\%thorn, "$arrangement/$thorn");

   ThornUtils::EndDocument($ofh, $document_type);
}

#/*@@
#  @routine   CreateLatexTable                                                     
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#    Takes whatever table element is currently reffered to by $table    
#    and prints it out into a LaTeX table.  Only nifty things it curr.  
#    does is NOT print ranges for BOOLEAN and SHARED elements.          
#  @enddesc 
#  @version 
#@@*/
sub CreateLatexTable
{
   my %thorn      = %{$_[0]};
   my $thorn_name = $_[1];

   my $printgridtitle = 1;
   my @conditional_statements;
   my @always_statements;
   my %aliases;
   my $len;

   # categorize the storage types for STATEMENTS into conditional and always on
   foreach my $key (sort keys %thorn) 
   {
      next if ($key !~ /^STATEMENT/);

      if ($thorn{$key}->{"TYPE"} eq "STOR") 
      {
         if ($thorn{$key}->{"CONDITIONAL"} == 1) {
            push @conditional_statements, split/,/, $thorn{$key}->{"GROUPS"};
         } elsif ($thorn{$key}->{"CONDITIONAL"} == 0) {
            push @always_statements, split/,/, $thorn{$key}->{"GROUPS"};
         }
      }
   }

   # which storage type has the most elements?
   $len = @conditional_statements > @always_statements ? @conditional_statements : @always_statements;

   # print blurb about what this page is about
   print "\n\n\\noindent This section lists all the variables which are assigned storage by thorn ". ThornUtils::CleanForLatex($thorn_name) . ".  Storage can either last for the duration of the run ({\\bf Always} means that if this thorn is activated storage will be assigned, {\\bf Conditional} means that if this thorn is activated storage will be assigned for the duration of the run if some condition is met), or can be turned on for the duration of a schedule function.\n\n";

   # print out storage allocation at the top
   print "\n\\subsection\*\{Storage\}";
   if (@conditional_statements > 0 || @always_statements > 0) 
   {
      print "\n\n\\hspace\{5mm\}\n\n \\begin\{tabular*\}\{${TABLE_WIDTH}mm\}\{ll\} \n";
      print @always_statements > 0 ? "\n\{\\bf Always:\}" : "~";
      print "& ";
      print @conditional_statements > 0 ? "\{\\bf Conditional:\} \\\\ \n"  : " ~ \\\\ \n";
 
      for (my $i = 0; $i <= $len; $i++) 
      {
         print $always_statements[$i] !~ /^$/ ? ThornUtils::CleanForLatex($always_statements[$i]) : "~";
         print " & ";
         print $conditional_statements[$i] !~ /^$/ ? ThornUtils::CleanForLatex($conditional_statements[$i]) : "~";
         print "\\\\ \n";
      } 
      print "\\end\{tabular*\} \n\n";
   } else {
       print "NONE";
   }

   # print out each scheduled block
   foreach my $block (sort keys %thorn) 
   {
      next if ($block !~ /^BLOCK/);
 
      # print the title, but only once
      if ($printgridtitle) {
         print "\n\\subsection\*\{Scheduled Functions\}";
         $printgridtitle = 0;
      }
 
      print "\n\\vspace\{5mm\}\n";
      print "\n\\noindent \{\\bf " .     ThornUtils::CleanForLatex($thorn{$block}->{"WHERE"}) . "\} ";
      print $thorn{$block}->{"CONDITIONAL"} ? "  (conditional) \n" : "\n";
      print "\n\\hspace\{5mm\} " .       ThornUtils::ToLower(ThornUtils::CleanForLatex($thorn{$block}->{"NAME"})) . " \n";
      print "\n\\hspace\{5mm\}\{\\it " . ThornUtils::ToLower(ThornUtils::CleanForLatex($thorn{$block}->{"DESCRIPTION"})) . " \} \n\n";
      print "\n\\hspace\{5mm\}\n\n \\begin\{tabular*\}\{${TABLE_WIDTH}mm\}\{cll\} \n";
 
      $aliases{$thorn{$block}->{"NAME"}} = $thorn{$block}->{"AS"};
 
      # Clean up the hash before we go iterate through and print out the rest
      delete $thorn{$block}->{"CONDITIONAL"};
      delete $thorn{$block}->{"WHERE"};
      delete $thorn{$block}->{"DESCRIPTION"};
      delete $thorn{$block}->{"AS"};
      delete $thorn{$block}->{"NAME"};
      delete $thorn{$block}->{"THORN"};

      # go print out the rest of the key/value pairs
      foreach my $group_key (sort keys %{$thorn{$block}}) {
         &OutputVar($group_key, $thorn{$block}->{$group_key});
      } # foreach sort keys %{$thorn{$group}}
 
      print "\\end\{tabular*\} \n\n";
   } # foreach %blocks
 
   # delete aliases where they key equals the value
   foreach my $key (sort keys %aliases) {
      if ($key eq $aliases{$key}) {
         delete $aliases{$key};
      }
   }
 
   # print out any Aliased functions in a table
   if (scalar(keys %aliases) > 0) 
   {
      print "\n\\subsection\*\{Aliased Functions\}";
      print "\n\n\\hspace\{5mm\}\n\n \\begin\{tabular*\}\{${TABLE_WIDTH}mm\}\{ll\} \n";
      print "\n\{\\bf Alias Name:\} ~~~~~~~ & \{\\bf Function Name:\} \\\\ \n";
 
      foreach my $key (sort keys %aliases) {
         print ThornUtils::CleanForLatex($key) ." & ". ThornUtils::CleanForLatex($aliases{$key}) ." \\\\ \n";
      }
      print "\\end\{tabular*\} \n\n";
   }

   print "\n\n\\vspace\{5mm\}"; 
}

#/*@@
#  @routine   OutputVar
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     Prints out a var and description for a table element
#
#  @enddesc 
#  @version 
#@@*/
sub OutputVar 
{
   my $description = shift;
   my $value       = shift;   
 
   # print out the different storage, we split it up, because it can get long
   return if ($value =~ /^$/);

   $description = defined $var_mapping{$description} ? $var_mapping{$description} : $description;
   $description = ThornUtils::Translate($description);

   # go through and print out the values, split them onto new lines if their
   # are multiple entries separated by comas
   if (my @temp = split/,/, $value)  
   {
      print "~ & ${description}:  & "; 
 
      my $fp = 1;
      foreach my $t (@temp) 
      {
         $t = ThornUtils::ToLower(ThornUtils::CleanForLatex($t));
 
         if ($fp) {
            print "$t \\\\ \n";
            $fp = 0;
         } else {
            print "~& ~ &" . $t . "\\\\ \n";
         }
      }
   } else {
      print "~ & ${description}:  &" . ThornUtils::ToLower(ThornUtils::CleanForLatex($value)) . "\\\\ \n";
   }
}
