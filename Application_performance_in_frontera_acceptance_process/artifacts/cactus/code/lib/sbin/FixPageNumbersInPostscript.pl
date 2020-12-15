#!/usr/bin/perl -w
#/*@@
#  @routine    FixPageNumbersInPostscript.pl
#  @date       Fri 15 Mar 2002
#  @author     Thomas Radke
#  @desc
#              Fixes the line numbers in the postscript versions of the
#              Cactus UsersGuide, ThornGuide, and MaintGuide
#  @enddesc
#  @version    $Header$
#  @history
#  @date       Sat Jul  5 17:55:07 CEST 2003
#  @author     Jonathan Thornburg <jthorn@aei.mpg.de>
#  @desc       Fix "previous line in file was end of previous page"
#              to also recognize the slightly different dvips output of
#              dvips(k) 5.92b (part of the teTeX 2.01 distribution)
#  @endhistory
#  @@*/

# $part counts the parts (chapters) in the postscript file
$part = 0;

# @part_letters lists the numbering for all the parts
# the first two parts (title page and table of contents) are skipped
@part_letters = ('', '', 'A' .. 'Z', 'a' .. 'z');

# $last_line must match EOP in order to check for a new page
$last2_line = "\n";
$last_line = "end\n";


# skip all lines in the postscript setup and prolog
while (<>)
{
  print;
  last if (/^%%EndProlog\n$/);
}


while (<>)
{
  my $previous_line_was_eop
	= (($last_line =~ /eop$/) || ($last_line =~ /eop end$/)
           || ($last2_line =~ /eop$/ && $last_line =~ /^end$/));
  if ($previous_line_was_eop && (/^%%Page: (\d+) (\d+)$/))
  {
    $part++ if ($1 == 1);
    $_ = "%%Page: ${part_letters[$part]}$1 $2\n"
      if defined(${part_letters[$part]});
  }

  print;
  $last2_line = $last_line;
  $last_line = $_;
}

1;
