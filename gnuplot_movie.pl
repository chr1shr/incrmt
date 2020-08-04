#!/usr/bin/perl
use Getopt::Std;
getopts("dhj:ln:p:tuw");
$gpm="../../shared/gpm_process";
$gpc="../../shared/make_contour";
$gpp="../../shared/gpp_info";

# Print help information if requested
if ($opt_h) {
    print "Usage: ./gnuplot_movie.pl {<options>} <filename> <suffix> {<z_min> <z_max>}\n";
    print "\nOptions:\n";
    print "-d           (Don't duplicate frames that already exist)\n";
    print "-h           (Print this information)\n";
    print "-j <n>       (Render every n frames)\n";
    print "-l           (Clean the grid up before rendering)\n";
    print "-p <n>       (Render output in parallel using n procs)\n";
    print "-n <n>       (Render up to this number of files)\n";
    print "-t           (Add tracers)\n";
    print "-u           (Add undeformed fields)\n";
    print "-w           (Render PNGs only without making a movie)\n";
    exit 0;
}

die "Need either two or four arguments" unless $#ARGV==1 || $#ARGV==3;

# Determine the number of physical cores
if(!defined $opt_p) {
    $uname=`uname`;
    if($uname=~/Linux/) {
        $nodes=`lscpu -p | egrep -v '^#' | sort -u -t, -k 2,4 | wc -l`;
        chomp $nodes;
    } elsif($uname=~/Darwin/) {
        $nodes=`sysctl -n hw.physicalcpu_max`;
        chomp $nodes;
    } else {
        $nodes=4;
    }
} else {
    $nodes=$opt_p;
}

# Make output directory
$e=@ARGV[0];
$ebase=$e;
$ebase=~s/\.odr$// or die "Filename should end in '.odr'\n";
$odir=$ebase.".frames";
mkdir $odir unless -e $odir;

# Set color range
if($#ARGV==3) {
    $cb="[@ARGV[2]:@ARGV[3]]";
} else {
    $cb="[*:*]";
}

$P=1;$a=0;$k=0;$da=$opt_j?$opt_j:1;

# Choose palette based on the field type
%pals=qw(w psy   p pag   rho byp   ds gab   dp gab   J lil   div lil);
$pal=defined $pals{$ARGV[1]}?$pals{$ARGV[1]}:"paw";
$set_pal=`$gpp gnu $pal`;

# Open header file if available
if(-e "${e}/header") {
    open A,"${e}/header" or die "Error reading header\n";
    $_=<A>;
    ($t_s,$t_e,$frn)=split;
    close A;
    $header=1;
}

# Set gnuplot header based on the filename. Strip out any directories from the filename.
my $ridx=rindex($ebase,'/')+1;
if($ridx) {
    $sim_type=substr($ebase, rindex($ebase,'/')+1);
    $sim_type=~s/\d//g;
} else {
    $sim_type=$ebase;
}
$sim_type="piston" if $sim_type=~/^piston_G/;
$sim_type="flag" if $sim_type=~/^f._/;
$sim_type="sediment" if $sim_type=~/^sed.*_/;
open A,"gp_headers/$sim_type.gnuplot" or open A,"gp_headers/default.gnuplot" or die "Can't find gnuplot header";

# Check for case of no phi field
$no_phi=1 unless -e "${e}/phi.0";

# Read the header file, skipping and altering lines as necessary
$gp="";
while(<A>) {

    # Palette and color bar substitution
    s/SET_PALETTE/$set_pal/;
    s/CBRANGE/$cb/;

    # Header substitution
    if($header) {s/^HEA://;} else {next if /^HEA:/;}

    # Tracer substitution
    if($opt_t) {s/^TRA://;} else {next if m/^TRA:/;}

    # Undeformed field substitution
    if($opt_u) {s/^UND://;} else {next if m/^UND:/;}

    # Phi field substitution
    next if $no_phi && m/LEVELSET/;
    $gp.=$_;
}
close A;

# Loop over the available frames
while(-e "${e}/$ARGV[1].$a") {

    # Terminate if the specified frame limit has been reached
    last if defined $opt_n && $a>$opt_n;

    # Prepare output and tracer filenames
    $za=sprintf "_%04d",$a;
    $of="$odir\/fr$za";
    $tfile="$e\/trace.$a";

    # Skip existing file if -d option is in use
    if($opt_d && -e $of && -M "@ARGV[0].$a" > -M $of) {
        print "$a (skipped)\n";
        $a+=$da;
        next;
    }

    # Clean up field if requested
    $lset="<$gpc $e/phi.$a - 0";
    $infile=$opt_l?"<$gpm -c $e/phi.$a $e/$ARGV[1].$a -":"$e/$ARGV[1].$a";

    # Call the utility to clean up the undeformed fields if they are in use
    if($opt_u) {
        $cmd=$no_phi?"$gpc":"$gpc -p $e/phi.$a";
        $xfield="<$cmd $e/X.$a - r -5 0.05 201";
        $yfield="<$cmd $e/Y.$a - r -5 0.05 201";
    }

    # Create temporary Gnuplot file
    print "Frame $a (thread $P)\n";

    # Make a fork to create the graph
    if(($pid[$P]=fork)==0) {
        $_=$gp;

        # File substitutions for tracers
        s/TRACERS/$tfile/ if $opt_t;

        # File substitutions for timing information
        if($header) {
            $ti=$t_s+($t_e-$t_s)*$a/$frn;
            $tif=sprintf "%.2f",$ti;
            s/TIME/$tif/;
        }

        # File substitutions for undeformed fields
        if($opt_u) {
            s/XFIELD/$xfield/;
            s/YFIELD/$yfield/;
        }

        # General file substitutions
        s/INFILE/$infile/;
        s/LEVELSET/$lset/;
        s/OUTFILE/$of/;

        # Send the plotting commands to Gnuplot
        open B,"|gnuplot -d";
        print B;
        close B;
        exit;
    }

    # Wait for one of the forked jobs to finish
    if($queue) {
        $piddone=wait;$P=1;
        $P++ while $piddone!=$pid[$P] && $P<=$nodes;
        die "PID return error!\n" if $P>$nodes;
    } else {
        $P++;$queue=1 if $P>=$nodes;
    };

    $a+=$da;
}

# Wait for all the remaining forked jobs to finish
wait foreach 1..($queue?$nodes:$P-1);

# Additional code to automatically make a movie
unless($opt_w) {
    $mf=$ebase."_".$ARGV[1];
    unlink "$mf.mov";
    system "ffmpeg -framerate 30 -y -i $odir/fr_%4d.png -preset veryslow -c:v libx265 -crf 17 -pix_fmt yuv420p -tag:v hvc1 -movflags faststart $mf.mov";

    # Delete the temporary output directory
    system "rm -rf $odir";
}
