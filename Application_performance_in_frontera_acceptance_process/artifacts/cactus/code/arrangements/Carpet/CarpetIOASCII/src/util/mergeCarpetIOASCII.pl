#! /usr/bin/perl -w
#/*@@
#  @file      mergeCarpetIOASCII.pl
#  @date      Tue 15 August 2006
#  @author    Thomas Radke
#  @desc
#             Perl script to merge CarpetIOASCII output files,
#             eliminating duplicate datasets.
#
#             Source code comments also by Thomas Radke :-)
#  @enddesc
#@@*/


my $help = $#ARGV < 0;
for (my $arg = 0; $arg <= $#ARGV; $arg++) {
  $help |= ($ARGV[$arg] eq '-h'    or
            $ARGV[$arg] eq '-help' or
            $ARGV[$arg] eq '--help');
}
if ($help) {
  print << "EOF";

  Usage: $0 [-h | -help | --help] <list of files>

  This script can be used to merge  CarpetIOASCII output  written before
  and after recovery. It reads one or more files in CarpetIOASCII format
  and  writes their contents  to STDOUT,  eliminating duplicate datasets
  (all but the last occurance are discarded).

    Example: $0 alp.x.asc > alp.x.asc.merged

EOF
  exit;
}

# Rauten-Feld zum Merken der Anzahl und Haeufigkeit vorhandener Datensaetze
my %datasets = ();

# Liste aller Eingabe-Dateien
my @filelist = @ARGV;

# lies zeilenweise alle Eingabe-Dateien
while (<>) {

  # falls diese Zeile einen neuen Datensatz einleitet:
  if (/^# iteration (\d+)$/) {

    # ermittle die Iterationsnummer aus der aktuellen Zeile ...
    my $iteration = $1;

    # ... und die anderen Kennwerte dieses Datensatzes aus der folgenden Zeile
    $_ = <>;
    die "Format error in file $ARGV line $.: expected '# refinement level ...'\n"
      unless (/^# refinement level (\d+)   multigrid level (\d+)   map (\d+)   component (\d+)   time level (\d+)$/);

    # vermerke den Datensatz mit seinen Parametern
    ++$datasets{"$iteration $1 $2 $3 $4 $5"};
  }
} continue {
  # setze die Zeilennummer fuer jede Eingabe-Datei wieder zurueck
  close ARGV if eof;
}

# stelle die Liste aller Eingabe-Dateien wieder her
@ARGV = @filelist;

# Flaggen-Variable, die anzeigt, ob die aktuelle Zeile ignoriert werden soll
my $discard = 0;

# lies zeilenweise alle Eingabe-Dateien
while (<>) {

  # falls diese Zeile einen neuen Datensatz einleitet:
  if (/^# iteration (\d+)$/) {

    # gib die aktuelle (Kommentar-)Zeile aus
    print;

    # ermittle die Iterationsnummer aus der aktuellen Zeile ...
    my $iteration = $1;

    # ... und die anderen Kennwerte dieses Datensatzes aus der folgenden Zeile
    $_ = <>;
    $_ =~ /^# refinement level (\d+)   multigrid level (\d+)   map (\d+)   component (\d+)   time level (\d+)$/;

    # erzeuge einen eindeutigen Schluessel fuer den aktuellen Datensatz
    my $key = "$iteration $1 $2 $3 $4 $5";

    # ueberspringe alle Datensaetze mit denselben Kennwerten bis auf den letzten
    $discard = --$datasets{"$iteration $1 $2 $3 $4 $5"};
  }

  print unless $discard;
}
