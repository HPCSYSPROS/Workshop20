#! /usr/bin/perl -w
#/*@@
#  @file      mergeCarpetIOScalar.pl
#  @date      Tue 15 August 2006
#  @author    Thomas Radke
#  @desc
#             Perl script to merge CarpetIOScalar output files,
#             eliminating duplicate timesteps.
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

  This script can be used to merge  CarpetIOScalar output  written before
  and after recovery. It reads one or more files in CarpetIOScalar format
  and  writes their contents  to STDOUT,  eliminating duplicate timesteps
  (all but the last occurance are discarded).

    Example: $0 alp.norm1.asc > alp.norm1.asc.merged

EOF
  exit;
}

# Rauten-Feld zum Merken der Anzahl und Haeufigkeit vorhandener Datensaetze
my %timesteps = ();

# Liste aller Eingabe-Dateien
my @filelist = @ARGV;

# lies zeilenweise alle Eingabe-Dateien
while (<>) {

  # falls diese Zeile eine Datenzeile ist:
  # merke den Datensatz mit seinem Zeitschritt
  $timesteps{$1} = $_ if (/^(\d+(\.\d+)?)\s+/);
}

# stelle die Liste aller Eingabe-Dateien wieder her
@ARGV = @filelist;

# lies zeilenweise alle Eingabe-Dateien
while (<>) {
  # gib alle Zeilen bis zur ersten Datenzeile aus
  last if (/^\d+(\.\d+)?\s+/);
  print;
}

# gib alle Datenzeilen aus, sortiert nach Zeitschritten
foreach my $timestep (sort {$a <=> $b} keys %timesteps) {
  print $timesteps{$timestep};
}
