#!/usr/bin/perl

# Print syntax message if not enough arguments are provided
die "Syntax: ./t_analysis.pl <timing_file>\n" if $#ARGV!=0;

# Check that the filename is the right format and can be opened
$ARGV[0]=~/\.tim$/ or die "Filename should end in '.tim'\n";
open A,$ARGV[0] or die "Can't open output file\n";

# Read in the lines of the file
$n=$steps=0;
$t_sim=$t_out=0.;
$v_mac=$t_mac=0.;
$v_fem=$t_fem=0.;
while(<A>) {

    # Skip any with comments
    next if /^#/;

    # Accumulate the timing information
    @a=split;
    $steps+=$a[1];
    $t_sim+=$a[2];$t_out+=$a[3];
    $v_mac+=$a[4];$t_mac+=$a[5];
    $v_fem+=$a[6];$t_fem+=$a[7];
    $n++;
};
close A;

# Print timing statistics
$ARGV[0]=~/(\d\d*)\.tim$/;
printf "$1 $n %d  %g %g  %g %g  %g %g  %g %g\n",$steps,$t_sim,$t_sim/$steps,$v_mac/$steps,$v_fem/$steps,1e6*$t_mac/$v_mac,1e6*$t_fem/$v_fem,
       100*$t_mac/$t_sim,100*$t_fem/$t_sim;
