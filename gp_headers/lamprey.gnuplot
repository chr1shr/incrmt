set pm3d
SET_PALETTE
set size ratio -1
unset surface
set view map
#xmin=-2;xmax=4;ymin=-1;ymax=1
xmin=-1.6;xmax=3.2;ymin=-0.8;ymax=0.8
set autoscale zfix
set autoscale cbfix
unset key
set xlabel 'x'
set ylabel 'y'
set cbrange CBRANGE
set term pngcairo size 1200,380 font "Helvetica, 16"
set lmargin at screen 0.08
set rmargin at screen 0.88
set tmargin at screen 0.95
set bmargin at screen 0.18
set output 'OUTFILE.png'
set multiplot
plot [xmin:xmax] [ymin:ymax] 'INFILE' matrix binary with image
unset xtics
unset ytics
unset xlabel
unset ylabel
unset clabel
unset pm3d
set surface
set style data lines
UND:set style line 1 linecolor rgb "#222222" lw 1
UND:plot [xmin:xmax] [ymin:ymax] 'XFIELD' u 1:2 ls 1
UND:plot [xmin:xmax] [ymin:ymax] 'YFIELD' u 1:2 ls 1
TRA:set style line 2 linecolor rgb "#000000"
TRA:set pointsize 0.2
TRA:plot [xmin:xmax] [ymin:ymax] 'TRACERS' binary format="%2float" u 1:2 w p ls 2 pt 7
HEA:set label 1 't = TIME' left at screen 0.0795,0.053 tc rgbcolor "#0000ff" front
set style line 3 lw 2 linecolor rgb "#000000"
plot [xmin:xmax] [ymin:ymax] 'LEVELSET' u 1:2 ls 3
unset multiplot
set output
