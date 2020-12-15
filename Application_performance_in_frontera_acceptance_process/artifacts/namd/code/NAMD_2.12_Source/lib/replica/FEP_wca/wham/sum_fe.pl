#!/usr/bin/perl -w

open(FE_REPU_TIME, ">solv_repu_fe.dat") or die "Can't open file lyz_repu_fe.dat: $!";
open(FE_BINDING_TIME, ">solv_total_fe.dat") or die "Can't open file lyz_binding_fe.dat: $!";


my @allrepufe = <repu_*_wham_fe>;
my $sum_fe = 0;
  open(FE_CHG, "<chg_wham_fe") or die "Can't open file : $!";
NEXT:  while ($line = <FE_CHG>) {
    chomp $line;
    if (!$line) {next NEXT;}
    my @fields = split /\s+/, $line;
    my $time = $fields[1];
    my $chg_fe = $fields[2];
$sum_fe += $chg_fe;
}
close FE_CHG;

  open(FE_DISP, "<disp_wham_fe") or die "Can't open file : $!";
NEXT:  while ($line = <FE_DISP>) {
    chomp $line;
    if (!$line) {next NEXT;}
    my @fields = split /\s+/, $line;
    my $time = $fields[1];
    my $disp_fe = $fields[2];
$sum_fe += $disp_fe;
}
close FE_DISP;

  my $total_repu_fe = 0;
nextmol: foreach $i (@allrepufe) {
  open(FE_REPU, "<${i}") or die "Can't open file ${i}: $!";
NEXT:  while ($line = <FE_REPU>) {
    chomp $line;
    if (!$line) {next NEXT;}
    my @fields = split /\s+/, $line;
    my $time = $fields[1];
    my $repu_fe = $fields[2];
$total_repu_fe += $repu_fe;
$sum_fe += $repu_fe;
}
close FE_REPU;
}
printf FE_REPU_TIME "%-10.5f\n", $total_repu_fe;
printf FE_BINDING_TIME "%-10.5f\n", $sum_fe;
close FE_REPU_TIME;
