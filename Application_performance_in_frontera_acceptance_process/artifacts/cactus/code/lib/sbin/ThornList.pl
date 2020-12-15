#!/usr/bin/perl

# use vars qw($arrangements_dir $arrangements $thorns $h $help $thornlist);

#/*@@
#  @file      ThornList.pl 
#  @date      March 2001
#  @author    Ian Kelley
#  @desc 
#      Either creates a stripped down thornlist of all thorns in the arrangments
#      directory, or prints out the full paths to each thorn on a single line
#
#      Used primary by the ThornGuide Makefile
#  @enddesc 
#  @version 
#@@*/

# give help
if ($h || $help) {
   print <<EOC;
--> ThornList.pl <--

Either creates a stripped down thornlist of all thorns in the arrangments directory, or prints out the full paths to each thorn on a single line.  It can also take a combination of arrangements and thorns and will print them all out as a thornlist.

Usage:
\t\$ perl -s ThornList.pl -arrangements_dir=/tmp/Cactus/arrangements > allthorns.th
\t\$ perl -s ThornList.pl -arrangements_dir=/tmp/Cactus/arrangements -thornlist=allthorns.th
\t\$ perl -s ThornList.pl -arrangements_dir=/tmp/Cactus/arrangements -thorns="CactusBase/Boundary CactusBase/IOBasic" -arrangements="CactusTest"
EOC

exit 0;
}

my $start_directory = `pwd`;
chomp ($start_directory);

my @thorns;
my @arrangements;


if (! defined $arrangements_dir) {
   die "\nArrangements directory not defined (-arrangments_dir)";
} elsif ($arrangements_dir !~ /\/$/) {
   $arrangements_dir .= '/';
}


# if we are being passed in the arrangements or thorns manually
if ($arrangements =~ /\w/ || $thorns =~ /\w/) {
   if (defined $arrangements) {
      @arrangements = split/,/, $arrangements;

      foreach my $arrangement (@arrangements) 
      { 
         @thorns = &FindDirectories("$arrangements_dir$arrangement");

         foreach my $thorn (@thorns) {
            print "\n$arrangement/$thorn" if ($thorn ne "doc");
         }
      }
   }

   if (defined $thorns) {
      @thorns = split/,/, $thorns;

      foreach my $thorn (@thorns) {
          print "\n$thorn";
      }
   }

}
# if we are building a thornlist from thorns in $arrangements_dir
elsif (! defined $thornlist) 
{
   @arrangements = &FindDirectories($arrangements_dir);

   foreach my $arrangement (@arrangements) 
   {
      @thorns = &FindDirectories("$arrangements_dir$arrangement");

      foreach my $thorn (@thorns) {
         print "$arrangement/$thorn\n" if ($thorn ne "doc");
      }
   }
} 
# if we are printing all the thorn directories on one line
# for use by the ThornGuide makefile
else 
{
   my $thorn_paths = "";

   open (THORNLIST, "<$thornlist") || die "\nCannot open thornlist ($thornlist): $!";

   while (<THORNLIST>) 
   {
      next if /\s*?\!/;          # bypass any directives
      s/(.*?)#.*/\1/;            # read up to the first "#"
      s/\s+//g;                  # replace any spaces with nothing
      
      next if ! /\w/;            # nothing on this line?

      $thorn_paths .= "$arrangements_dir$_ ";

   }

   close (THORNLIST);

   print $thorn_paths;
}
## END: MAIN ##

#/*@@
#  @routine   FindDirectories
#  @date      Sun Mar  3 01:54:37 CET 2002
#  @author    Ian Kelley
#  @desc 
#     Grabs all directories within a given directory, minus the ones created by CVS
#  @enddesc 
#  @version 
#@@*/

sub FindDirectories 
{
   my $search_dir = shift;
   my @good_directories;
   my $dirhdl;

   opendir ($dirhdl, $search_dir)
            or die "\nCannot open directory $search_dir\n";

   while (defined (my $name = readdir($dirhdl)))
   {
      next if (! -d "$search_dir/$name");

      if (($name ne 'History') && ($name ne 'CVS')
      && !($name =~ /^\./ ) )        # i.e. current, parent & hidden dirs 
      {
         push(@good_directories, $name);
      }
   }
   closedir $dirhdl;

   return @good_directories;
} ## END :Find_Directories:
