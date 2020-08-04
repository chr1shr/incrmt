set pm3d
SET_PALETTE
set size ratio -1
unset surface
set view map
xc=0;yc=0;xw=1.0;yw=1.0
#xc=0.5;yc=-0.5;xw=0.5;yw=0.5
xmin=xc-xw;xmax=xc+xw;ymin=yc-yw;ymax=yc+yw
#xmin=-0.5;xmax=0.5;ymin=-1;ymax=1
set autoscale zfix
set autoscale cbfix
unset key
set xtics offset 0,0.4
set xlabel 'x' offset 0,0.7
set ylabel 'y'
set cbrange CBRANGE
set term pngcairo size 780,680 font "Helvetica, 16"
set lmargin at screen 0.08
set rmargin at screen 0.86
set tmargin at screen 0.955
set bmargin at screen 0.12
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
HEA:set label 1 't = TIME' left at screen 0.112,0.03 tc rgbcolor "#0000ff" front
HEA:set label 2 '{/Symbol w}' left at screen 0.96,0.539 front
set style line 3 lw 2 linecolor rgb "#000000"
plot [xmin:xmax] [ymin:ymax] 'LEVELSET' u 1:2 ls 3
unset multiplot
set output
