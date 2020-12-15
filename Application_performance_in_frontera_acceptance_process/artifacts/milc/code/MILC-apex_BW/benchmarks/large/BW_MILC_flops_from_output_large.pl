#!/usr/bin/perl

$fp_ops=5.17825e16;

$target_filename=$ARGV[0];
$milc_output_filename=$ARGV[1];
if(! $target_filename){
    print "\n";
    print "First argument must be name of job script output from MILC run.\n";
    print "Second argument must be name of MILC output file\n";
    print "  (of the form milc_??x??x??x???_...) (which will be in a results_* directory)\n";
    print "\n";
}
if( ! ( -f $target_filename ) ){
    print "file >$target_filename< does not exist!\n";
    print "please input name of MILC job output file as the \n";
    print "first argument.\n";
    exit;
}
$grep_start_command="grep UNIX_TIME_START $target_filename|";
$grep_end_command="grep UNIX_TIME_END $target_filename|";

open(GREP_START,$grep_start_command) or die "grepping start command returned error!\n";
open(GREP_END,$grep_end_command) or die "grepping end command returned error!\n";

$grep_start_out=<GREP_START>;
$grep_end_out=<GREP_END>;

chomp $grep_start_out;
chomp $grep_end_out;

#print "grep output: >$grep_out<<\n";

if( ! ( $grep_start_out =~ m/UNIX_TIME_START=(\d+)/ ) ){
    print "could not find START time in MILC job output file.  Exiting.\n";
    exit;
}
$unix_start_time=$1;

if( ! ( $grep_end_out =~ m/UNIX_TIME_END=(\d+)/ ) ){
    print "could not find END time in MILC job output file.  Exiting.\n";
    exit;
}
$unix_end_time=$1;

$walltime_seconds=$unix_end_time-$unix_start_time;
print "\n";
print "total time was $walltime_seconds seconds.\n";
$total_flops=$fp_ops/$walltime_seconds;
$formatted_total_flops=sprintf("%6.4G",$total_flops);
print "total application floating point operations per second: $formatted_total_flops\n\n";

if( $milc_output_filename =~ m/milc_\d+x\d+x\d+x\d+_on_(\d+)_(\d+)_\d+/ ){
    $n_nodes = $1;
    $ranks_per_node = $2;
    $flops_per_node = $total_flops / $n_nodes;
    $flops_per_rank = $flops_per_node / $ranks_per_node;
#    print " ranks=$ranks_per_node flops/node=$flops_per_node flops/rank=$flops_per_rank\n";
    $formatted_gf_per_node = sprintf("%6.4G",$flops_per_node/1.0e9);
    $formatted_gf_per_rank = sprintf("%6.4G",$flops_per_rank/1.0e6);
    print "assuming we've parsed the milc report filename right,\n";
    print "this code was run on $n_nodes computational nodes,\n";
    $total_ranks=$n_nodes * $ranks_per_node;
    print "at $ranks_per_node cores per node, on a total of $total_ranks cores.\n";
    print "The final per-node result is $formatted_gf_per_node GF/node\n";
    print "and $formatted_gf_per_rank MF/core.\n\n";
}
