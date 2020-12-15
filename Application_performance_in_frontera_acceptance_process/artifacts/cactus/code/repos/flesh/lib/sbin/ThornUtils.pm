package ThornUtils;

my $parskip_set     = "0pt";
my $parskip_restore = "10pt";

#/*@@
#  @file      ThornUtils.pm
#  @date      Sun Mar  3 19:05:41 CET 2002
#  @author    Ian Kelley
#  @desc 
#  This file contains any common routines that are used by the scripts to create the
#  Cactus ThornGuide.  Currently the files that implement this package are:
#     InterLatex.pl, ParamLatex.pl, SchedLatex.pl, ThornGuide.pl
#  @enddesc 
#  @version 
#@@*/

#/*@@
#  @routine   ReadThornlist 
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#  Reads in a thornlist and returns the arrangements/thorns, strips out all the comments/etc. 
#
#  THIS FUNCTION IS PROVIDED BY "MakeUtils.pl"
#
#  @enddesc 
#  @version 
#@@*/

#/*@@
#  @routine   CreateThornlist 
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#  Takes an input directory (arrangements directory), calls &FindDirectories to get all 
#  the sub-directories (arrangements), then calls &FindDirectories once more to get all
#  the sub-sub-directories (thorns), then returns a thornlist.
#     @thorns = &CreateThornlist($arrangements_dir);
#  @enddesc 
#  @version 
#@@*/
sub CreateThornlist 
{
   my ($arrangementsDir) = shift;

   my (@foundThorns, @arrangements, @thorns);

   @arrangements = &FindDirectories($arrangementsDir);

   foreach my $arrangement (@arrangements)
   {
      @thorns = &FindDirectories("${arrangementsDir}${arrangement}");

      foreach my $thorn (@thorns) {
         push @foundThorns, "${arrangement}/${thorn}";
      }
   }

   return @foundThorns;
}

#/*@@
#  @routine   FindDirectories
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#  Grabs all directories within a given directory, minus the ones created by CVS 
#     @directories = &FindDirectories($directory);
#  @enddesc 
#  @version 
#@@*/
sub FindDirectories 
{
   my $directory = shift;

   my @directories;
   my $dirhdl;

   opendir ($dirhdl, $directory)
            or die "\nCannot open directory $directory\n";

   while (defined (my $name = readdir($dirhdl)))
   {
      if ((-d "$directory/$name")
              && ($name ne 'History') && ($name ne 'CVS')
              && !($name =~ /^\./) )      # i.e. current, parent & hidden dirs
      {
        push(@directories, $name);
      }
   }

   closedir $dirhdl;

   return @directories;
} 

#/*@@
#  @routine    GetArrangementsDir 
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  assumes $cctk_home and $start_directory are already defined
#     &GetArrangementsDir(<indir>)
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub GetArrangementsDir 
{
   my ($indir) = shift;

   if ($indir =~ /^$/) {
      $indir = "${cctk_home}arrangements";
   } elsif ($indir !~ /^\//) {
      $indir = "$start_directory/arrangements";
   } 

   $indir .= '/' if ($indir !~ /\/$/);

   return ($indir);
}

#/*@@
#  @routine    SetupOutputDirectory
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  Sets up the output directory, assigning it to the CWD if it is blank, trys to
#  resolve it if it is relative (by using $start_directory), 
#     &SetupOutputDirectory($outdir);
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub SetupOutputDirectory 
{
   my ($outdir) = shift;
 
   # set $outdir to the CWD
   if ($outdir =~ /^$/) {
      $outdir = "./";
   } 
   else 
   {
      # we have an relative path, tack $cctk_home on the front
      if ($outdir !~ /^\//) {
         $outdir = "$start_directory/$outdir";
      }

      # create the directory if it does not exist
      if (! -d $outdir) {
         mkdir($outdir, 0755)  || die "\nCannot create directory ($outdir): $!";
         print STDERR "\nCreating directory: $outdir" if ($verbose);
      }
      $outdir .= '/' if ($outdir !~ /\/$/);
   }

   print STDERR "\nOutput directory is: $outdir" if ($verbose);

   return ($outdir);
}

#/*@@
#  @routine    ClassifyThorns 
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc 
#  Goes through and classifies thorns based upon arrangement, takes a hash reference and an array
#  listing thorns (e.g. CactusWave/WaveToyC) as input, and then assigns the keys to be the 
#  arrangement names, with values as references to an array of thorn names
#
#     &ClassiftyThorns(<hash reference>, <array>)
#
#  To fill a hash by using &ClassifyThorns:
#     %arrangements;
#     &ClassifyThorns(\%arrangements, @listOfThorns);
#
#  To iterate through the hash, printing out all values:
#     foreach my $arrangement (keys %arrangements) {
#        print "\n$arrangement";
#        foreach my $thorn (@{$arrangements{$arrangement}}) {
#           print "\n\t$thorn";
#        }
#     }
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#
#@@*/
sub ClassifyThorns 
{
   my ($rhArr) = shift;

   # My most efficient code yet! (I am so proud)  

   # We add keys (arrangements) to the hash passed in ($rhArr) which are references 
   # to a list of thorns for each arrangement.
   foreach (@_) {
      /^(.*?)\/(.*?)$/;
      push @{$$rhArr{$1}}, $2;
   }
}


#/*@@
#  @routine    GetThornPaths
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  This function will cut up a list of thorns (@listOfThorns), check if a given file
#  exists in that directory ($interestingFile), if so, then it will add a key/value
#  pair to the hash (%paths) that represents the thorn name (key) and location (value)
#    %paths = &GetThornPaths(\@listOfThorns, $arrangementsDir, $interestingFile);
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub GetThornPaths 
{
   my ($rlListOfThorns)  = shift;
   my ($arrangementsDir) = shift;
   my ($interestingFile) = shift;
   my ($keepArrInKey)    = shift;

   my %paths;

   foreach (@$rlListOfThorns) {
      /^(.*?)\/(.*?)$/;
      if (($interestingFile eq "") || (-e "${arrangementsDir}${_}/${interestingFile}")) {
         my $key = $keepArrInKey ? "$1/$2" : "$1";
         $paths{$key} = "$arrangementsDir$_";
      }
   }
   return %paths;
}

#/*@@
#  @routine    StartDocument
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  Opens a file for output, and selects it as standard out. returning the old
#  filehandle.  The output filename is determined depending on how things are 
#  grouped (by thorn, arrangement, or all together) and the type of the program
#  that invoked the function (interface, schedule, or parameter)
#
#  Depending on the $docType, this will either create a 'section' to be included
#  by a  a larger .tex file, or it will create a new self-contained file.
#
#  Example, to create an 'Interfaces section' for thorn CactusBase/IOBasic
#     $ofh = &StartDocument("inter", "IOBasic", "/tmp/", "CactusBase", "Interfaces", "section");
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub StartDocument 
{
   my ($programType, $groupName, $outdir, $arrangement, $sectionName, $docType) = @_;
   my $oldfilehandle;

   my $outfile = $outdir;

   if (! defined $groupName) {
      die "\nCannot group tables for output, internal error in &StartDocument";
   }

   chdir ($start_directory) || die "\nCannot change directory to $start_directory";

   $outfile .= $arrangement !~ /^$/ ? "${arrangement}_" : "";
   $outfile .= "${groupName}_${programType}.tex";

   open(OUT, ">$outfile") || die "\nCannot open $outfile for output: $!";
   $oldfilehandle = select OUT;

   if ($docType eq 'document') {
      print "\\documentclass[12pt,a4paper]\{article\} \n";
      print "\\begin\{document\} \n\n";
      print "\n\\noindent \\section{$sectionName}\n\n";
   } elsif ($docType eq 'section') {
      print "\n\\section{$sectionName} \n\n";
   }
   print "\n\\parskip = $parskip_set\n";

   return $oldfilehandle;
}

#/*@@
#  @routine    EndDocument
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  Ends a LaTeX document (if -document) specified, and closes the OUT 
#  file, which is used to write the LaTeX table to. 
#     &EndDocument(<old file handle>, <section, document>);                  
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub EndDocument 
{
   my ($oldfilehandle, $docType) = @_;

   print "\\parskip = $parskip_restore \n";

   if ($docType eq 'document') {
      print "\\end\{document\} \n";
   }
   close OUT;

   select $oldfilehandle;
}

#/*@@
#  @routine    CleanFromC
#  @date       Aug 20 2004
#  @author     Erik Schnetter
#  @desc
#  "Interprets" a C string, i.e., takes something that would be valid
#  in C (or CST) source code inside double quotes, and return the
#  string that it represents, i.e., it mainly removes backslashes.
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub CleanFromC
{
   my $val = shift;

   # unescape special characters
   $val =~ s,\\\",\",g;

   # do not unescape backslash-letter sequences;
   # we assume that people want to see them instead of their effects

   return $val;
}

#/*@@
#  @routine    CleanForLatex
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  Cleans up our values so that latex will not give us errors.
#     $val = &CleanForLatex($val);
#  Note: Do not call ToLower or ToUpper on the result; instead,
#        transform before you clean for Latex.
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub CleanForLatex
{
   my $inval = shift;

   # at start of string, remove spaces before: "
   $inval =~ s,^\s*\",\",;

   # at end of string, remove spaces after: "
   $inval =~ s,\"\s*$,\",;

   # escape special characters
   my $outval = "";
   foreach my $char (split(//,$inval))
   {
     if ($char eq '{' or
         $char eq '}' or
         $char eq '$' or
         $char eq '_' or
         $char eq '^' or
         $char eq '&' or
         $char eq '%' or
         $char eq '#')
     {
       $outval .= '\\' . $char;
     }
     elsif ($char eq '\\')
     {
       $outval .= '{\\textbackslash}';
     }
     elsif ($char eq '~')
     {
       $outval .= '{\\textasciitilde}';
     }
     elsif ($char eq '<')
     {
       $outval .= '{\\textless}';
     }
     elsif ($char eq '>')
     {
       $outval .= '{\\textgreater}';
     }
     else
     {
       $outval .= $char;
     }
   }

   return $outval;
}

#/*@@
#  @routine    Translate
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  Takes a value, uppercases the first letter and lowercases all others
#     $val = &Translate($val); 
#  Note: Do not call this routine on a result of CleanForLatex; instead,
#        transform before you clean for Latex.
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub Translate 
{
   my $val = shift;
   my $temp = "";

   $val  =~ s/^(.)//;
   $temp = $1;
   $temp =~ tr/a-z/A-Z/;
   $val  =~ tr/A-Z/a-z/;

   return ($temp.$val);
}


#/*@@
#  @routine    ToUpper
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  Translates values passed in to upper case and returns it
#     $val = &ToUpper($val);
#  Note: Do not call this routine on a result of CleanForLatex; instead,
#        transform before you clean for Latex.
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub ToUpper
{
   my ($val) = shift;

   $val =~ tr/a-z/A-Z/;
   return $val;
}

#/*@@
#  @routine    ToLower
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  Translates values passed in to lower case and returns it
#     $val = &ToLower($val);
#  Note: Do not call this routine on a result of CleanForLatex; instead,
#        transform before you clean for Latex.
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub ToLower 
{
   my ($val) = shift;

   $val =~ tr/A-Z/a-z/;
   return $val;
}

#/*@@
#  @routine    CreateSystemDatabase
#  @date       Sun Mar  3 19:05:41 CET 2002
#  @author     Ian Kelley 
#  @desc
#  Creates a system database, used by create_interface_database.
#     &CreateSystemDatabase(<hash reference>);
#  @enddesc 
#  @calls     
#  @calledby   
#  @history 
#
#  @endhistory 
#  
#@@*/
sub CreateSystemDatabase 
{
   my ($rh) = shift;

   $$rh{"CCTK_HOME"} = $cctk_home;

   if (! defined $config_dir) {
      $$rh{"CONFIG_DIR"} = "${cctk_home}config-data";
   }
  
   if (! defined $bindings_dir) {
      $$rh{"BINDINGS_DIR"} = "${cctk_home}bindings";
   }
}

#/*@@	      
#  @routine   Dump
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#  This routine will 'dump' out some data-structure when given a reference to it. It becomes useful
#  (for debugging) to dump the contents of something like:
#     $database{"Computer"}->{"Programming Languages"}->{"Scripting"}->{"Perl"}->{"Users"} = \@users;
#  &Dump(\%database);
#
#  It should be noted that this is primary dealing with references to stuff, so if
#  this routine is used on something that is not a reference, it may not work properly.
#     (e.g. &Dump(@myarray) will only print the first element of @myarray)
#         but: &Dump(\@myarray) will print them all... etc etc
#
#  Note:  unless it isn't abundently obvious, this is all recurrsion.
#  @enddesc 
#  @calls 	Dump DumpHash DumpArray 
#  @version 
#@@*/
sub Dump 
{
   my ($ref_type) = ref($_[0]);
   
   print " [$ref_type]" if ($debug>1);
   if ($ref_type eq "HASH") {
      # we want to print out a hash, so we call the routine that 
      # knows how to properly iterate through one
      &DumpHash(@_);
   } elsif ($ref_type eq "ARRAY") {
      # routine that iterates through an array/list
      &DumpArray(@_);
   } elsif ($ref_type eq "REF") {
      # go to next level, maybe someday we get a real 'value'
      &Dump(${$_[0]}, $_[1]);
   } else {
      # scalar, or whatever, print it out.
      print "\n" . " " x $_[1] . $_[0];
   }

   return;
}

#/*@@
#  @routine   DumpHash
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#  Outputs a hash, or, if that hash contains references, calls &Dump
#  
#  Note: This is an INTERNAL FUNCTION
#  @enddesc 
#  @calledby  Dump
#  @version 
#@@*/
sub DumpHash 
{
   my ($hr, $iteration) = @_;

   foreach my $key (keys %$hr) 
   {
      print "\n" . " " x $iteration . $key;
      if (! (my $refType = ref($$hr{$key}) )) {
         print " = $$hr{$key}";
      } else { 
         &Dump($$hr{$key}, $iteration+3);
      }
   }
}

#/*@@
#  @routine   DumpArray 
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#  Outputs an array, or, if that array contains references, calls &Dump
#
#  Note: This is an INTERNAL FUNCTION
#  @enddesc 
#  @calledby  Dump
#  @version 
#@@*/
sub DumpArray 
{
   my ($ar, $iteration) = @_;

   foreach my $key (@$ar) 
   {
      if (! (my $refType = ref($key) )) {
         print "\n" . " " x $iteration . $key;
      } else { 
         &Dump($key, $iteration+3);
      }
   }
}
 
#/*@@
#  @routine   ProcessAllArrangements
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#  Will loop through all keys of the hash reference passed in, then for
#  each key, it will call &ProcessOneArrangement, this routine has been separated
#  from ProcessOneArrangement, so that ProcessOneArrangement can be called
#  independently if needed, and to add more flexiblity in the future.
#     &ProcessAllArrangements(\%arrangements);
#  @enddesc 
#  @calls     ProcessOneArrangement
#  @calledby  
#  @version 
#@@*/
sub ProcessAllArrangements
{
   # get the arrangements hash
   my ($hrArrangements) = shift;

   # go through and find  each arrangement
   foreach my $arrangement (keys %$hrArrangements) {
      &ProcessOneArrangement($$hrArrangements{$arrangement}, $arrangement);
   }
}

#/*@@
#  @routine   ProcessOneArrangement
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#  Will loop through all keys of the hash reference passed in, then for
#  each key, it will call &ProcessOneThorn in the MAIN perl file to process
#  the appropriate thorn (e.g. output latex)
#
#  It calls the &ProcessOneThorn in the MAIN perl file because depending on 
#  what type of parsing one is doing (interface, schedule, parameter), the
#  keys, sort mechanisms and datatypes in the (thorn) hash will differ.
#     &ProcessOneArrangement(\%arrangement);
#  @enddesc 
#  @calls     ProcessOneThorn (in MAIN perl file)
#  @calledby  ProcessAllArrangements
#  @version 
#@@*/
sub ProcessOneArrangement 
{
   # get the (single) arrangement hash
   my ($hrArr) = shift;
   my ($arrangement) = shift;

   foreach my $thorn (keys %$hrArr) {
      # please note, you will not get good results unless a sub-routine called
      # ProcessOneThorn exists in your main perl file.  (where this is called from)
      main::ProcessOneThorn($$hrArr{$thorn}, $arrangement, $thorn);
   }
}


#/*@@
#  @routine   AddQuotes
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     Simply takes away any enclosing quotes if they exist, and then adds enclosing quotes.
#  @enddesc 
#  @version 
#@@*/
sub AddQuotes {
   my $var = shift;
   $var =~ s/^\s*?\"(.*)\"\s*?$/$1/;
   return "\"$var\"";
}

#/*@@
#  @routine   RemoveQuotes
#  @date      Sat Apr  20 2002
#  @author    Gabrielle
#  @desc 
#     Simply takes away any enclosing quotes if they exist.
#  @enddesc 
#  @version 
#@@*/
sub RemoveQuotes {
   my $var = shift;
   $var =~ s/^\s*?\"(.*)\"\s*?$/$1/;
   return $var;
}

#/*@@
#  @routine   SetWidth
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     This function defines some latex variables for setting the width of the table we will later display.
#     It does this by defining the total table width, and then subtracting the different values for the 
#     table elements, so we can use these variables to set the widths of different cells in the table.
#
#  @enddesc 
#  @version 
#@@*/
sub SetWidth
{
   my $width       = shift;
   my $longest_var = shift;

   my $latex_output = "";

#\\newlength{\\tableWidth}

#\\newlength{\\maxVarWidth}
#\\newlength{\\paraWidth}
#\\newlength{\\descWidth}
$latex_output = <<END;

\\setlength{\\tableWidth}{${width}mm}

\\setlength{\\paraWidth}{\\tableWidth}
\\setlength{\\descWidth}{\\tableWidth}
\\settowidth{\\maxVarWidth}{$longest_var}

\\addtolength{\\paraWidth}{-\\maxVarWidth}
\\addtolength{\\paraWidth}{-\\columnsep}
\\addtolength{\\paraWidth}{-\\columnsep}
\\addtolength{\\paraWidth}{-\\columnsep}

\\addtolength{\\descWidth}{-\\columnsep}
\\addtolength{\\descWidth}{-\\columnsep}
\\addtolength{\\descWidth}{-\\columnsep}
END

return ($latex_output);
}

#/*@@
#  @routine   ChopVariable
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     If a variable has a sequence of numbers longer than the given value ($split) that does
#     not contain a space, it will add a space and newline after this sequence
#
#  @enddesc 
#  @version 
#@@*/
sub ChopVariable 
{
   my $var   = shift;
   my $split = shift;

   $var =~ s/([^\s]{$split})/$1 \n/g;
   return $var;
}


1;
