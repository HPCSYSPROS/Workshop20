#/*@@
#  @file      CVSUpdate.pl
#  @date      Tue Nov 21 2000
#  @author    Gabrielle Allen
#  @desc
#     Updates Cactus checkout
#     (avoids problems with different versions of cvs client)
#  @enddesc
#  @version $Header$
#@@*/

my $cvs_ops="-z6 -q";
my $cvs_update_ops="-d -P";
my $svn_ops="";
# Set this to eg -r TAGNAME checkout from a TAG
my $cvs_symbolic_name="";
my $svn_symbolic_name="";

require "lib/sbin/MakeUtils.pl";

$debug = 0;
if ($debug)
{
  print "DEBUG mode: cvs/svn commands not issued\n\n";
}

print("Not updating Flesh - not in CVS anymore\n");
#$command = "cvs $cvs_ops update $cvs_update_ops $cvs_symbolic_name CONTRIBUTORS COPYRIGHT Makefile lib doc src arrangements/README";
#if ($debug)
#{
#  $this_dir = `pwd`;
#  chop($this_dir);
#  print "\nIn directory $this_dir\n";
#  print "Issuing command\n  $command\n";
#  foreach $file (`ls CVS`)
#  {
#    chop($file);
#    print "Contents of $file\n";
#    open (FILE, "<CVS/$file") || die "Could not open CVS file";
#    while (<FILE>)
#    {
#      print;
#    }
#  }
#}
#if (!$debug)
#{
#  open (CS, "$command |");
#  while (<CS>)
#  {
#    print ;
#  }
#  close (CS);
#}

$home = `pwd`;
chomp ($home);

($arrangement_dir, $thornlist) = @ARGV;

$arrangement_dir = "$home/$arrangement_dir";

if ($thornlist =~ /^$/) {
   %info = &buildthorns($arrangement_dir,"thorns-to-update");
} else {
   %info = &ReadThornlist($thornlist);
}

foreach $thorn (sort keys %info)
{
  $arrangement = $thorn;
  $arrangement =~ s/[\/]+[a-zA-Z0-9_]*//;
  
  if ( ! $visited_arrangements{"\U$arrangement\E"})
  {
    $visited_arrangements{"\U$arrangement\E"} = 1;
    if( -d "$arrangement_dir/$arrangement/doc" && -d "$arrangement_dir/$arrangement/doc/CVS")
    {
      chdir ("$arrangement_dir/$arrangement/doc") ||
        die "Cannot change to arrangement directory '$arrangement_dir/$arrangement'\n";
      print("Updating arrangement $arrangement\n");
      $command = "cvs $cvs_ops update $cvs_update_ops $cvs_symbolic_name";
      if($debug)
      {
        $this_dir = `pwd`;
        chop($this_dir);
        print "In directory $this_dir\n";
        print "Issuing command\n  $command\n";
        foreach $file (`ls CVS`)
        {
          chop($file);
          print "Contents of $file\n";
          open (FILE, "<CVS/$file") || die "Could not open CVS file";
          while (<FILE>)
          {
            print;
          }
        }
      }
      if (!$debug)
      {
        open (CS, "$command |");
        while (<CS>)
        {
          print ;
        }
      }    
    }
  }

  if( -d "$arrangement_dir/$thorn/CVS")
  {
    $command = "cvs $cvs_ops update $cvs_update_ops $cvs_symbolic_name";
  }
  if( -d "$arrangement_dir/$thorn/.svn")
  {
    $command = "svn $svn_ops update $svn_symbolic_name";
  }
  if ( ! -d "$arrangement_dir/$thorn/CVS" &&
       ! -d "$arrangement_dir/$thorn/.svn" )
  {
    print "Ignoring $thorn - no CVS or .svn directory\n";
    next;
  }

  chdir ("$arrangement_dir/$thorn") ||
    die "Cannot change to thorn directory '$arrangement_dir/$thorn'\n";
  print("Updating thorn $thorn\n");
  if($debug)
  {
    $this_dir = `pwd`;
    chop($this_dir);
    print "In directory $this_dir\n";
    print "Issuing command\n  $command\n";
    if ( -d "$arrangement_dir/$thorn/CVS" )
    {
      foreach $file (`ls CVS`)
      {
        chop($file);
        print "Contents of $file\n";
        open (FILE, "<CVS/$file") || die "Could not open CVS file";
        while (<FILE>)
        {
          print;
        }
      }
    }
  }
  if (!$debug)
  {
    open (CS, "$command |");
    while (<CS>)
    {
      print ;
    }
  }
}
  chdir $home || die "Cannot change back to Cactus home directory '$home'\n";


exit;
