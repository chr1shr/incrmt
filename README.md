# A Method for Self-Contact Within IncRMT
This branch contains the numerical implementation of a method for self-repulsion, 
allowing submersed bodies to come into contact with themselves. Previously, the 
problem with such contacts was a failure in the extrapolation of the reference 
map field. Upon nearing a portion of itself, a solid's reference map extrapolation 
procedure breaks down, being unable to evolve the reference map into a non-solid 
region. By implementing a self-contact stress which detects when a body approches 
itself, we use an additional stress tensor term to update the velocity field in 
such a way to avoid self-contact. Since we can handle collisions of a single 
body with itself, we have also developed the capability of modeling many bodies 
with a single levelset. Where previously we treated multiple bodies with their 
own levelset fields, we can now evolve one levelset containing all bodies and 
treat all contacts as self-collisions. The code here can perform a variety of 
simulations and carries out all examples and tests described in the following 
publication:



This branch modifies the existing IncRMT code to support the recent publication of 
Lara *et al.* listed above. The installation and build steps are unmodified.

## Example
The simple U example described in section 5.1 of the paper can be run with the 
following command using four threads:
```Shell
OMP_NUM_THREADS=4 ./ftest Ucurve 133
```
This is the example described in the paper using 133 grid units along the vertical 
axis, the number of horizontal grid units is scaled to create a square grid. This code 
will create a **Ucurve.odr** directy for simulation outputs. The output directory 
contains the following filetypes:
- **contdivx**.*<n>* and **contdivy**.*<n>*, the self-repulsion force vectors at 
frame *n* stored in a binary format;
- **alpha**.*<n>*, the alpha magnitude of the self-contact stress tensor at 
frame *n* stored in a binary format;
- **p**.*<n>*, the pressure field at frame *n* stored in binary format;
- **w**.*<n>*, the vorticity field at frame *n* stored in binary format;
- **X**.*<n>* and **Y**.*<n>*, the components of the reference map at frame *n*;
- **phi**.*<n>*, the levelset field at frame *n* stored in binary format;
- **header**, a small text file containing the number of simulation frames and
  the time interval simulated.

The outputs can be modified near the top of the **ftest.cc** code. To save memory,
only every 50th frame has its outputs saved. This can be modified near the 
top of the **fluid_2d.cc** code.

In Gnuplot, the alpha field can be plotted using the following command at output *n=1000*:
```Gnuplot
set view 0,90,,
set pm3d map
splot 'Ucurve.odr/alpha.1000' matrix binary
```
Using the multiplot environment, the levelset edge can be visualized:
```Gnuplot
set view 0,90,,
set pm3d map
set multiplot
splot 'Ucurve.odr/alpha.1000' matrix binary notitle
set contour base
set cntrparam levels discrete 0
splot 'Ucurve.odr/phi.1000' matrix binary notitle with line
unset multiplot
```
The self-contact forces can be overlayed with a vector field as follows:
```Gnuplot
set table $Coords
plot Ucurve.odr/contdivx.1000 using 1:2 matrix binary with table
unset table
set table $contdivx
plot Ucurve.odr/contdivx.1000 using 3 matrix binary with table
unset table
set table $contdivy
plot Ucurve.odr/contdivy.1000 using 3 matrix binary with table
unset table

set print $Vec
do for [i=1:|$Coords|] {
if (abs($contdivx[i])>0 || abs($contdivy[i])>0) {
print $Coords[i], $contdivx[i], $contdivy[i], 0, 0
}
}
unset print

set view 0,90,,
set pm3d map
set multiplot
splot 'Ucurve.odr/alpha.1000' matrix binary notitle
set contour base
set cntrparam levels discrete 0
splot 'Ucurve.odr/phi.1000' matrix binary notitle with line
unset multiplot
splot $Vec using 1:2:5:3:4:5 with vectors notitle

```
For more plotting commands, consult the [Gnuplot documention](http://www.gnuplot.info/documentation.html)

Many other types of simulation are possible with the **ftest** code, most
of which are taken from the associated publication. To see a complete list
type
```Shell
./ftest
```

## Code structure
In addition to the existing IncRMT structure, our method for self-repulsion introduces two new globally defined field using the **field** data structure. First a field storing levelset gradients, calculated at each timestep using a finite difference approximation. These are used to construct self-contact stress terms. Since this field is globally defined and object independent, all levelsets are used to calculate it, removing the necessity for multiple levelsets for many object simulations. Secondly is a field storing the magnitude of the self-contact stress tensor. While this field is not strictly necessary for simulations, it is helpful to record as a diagnostic and visualization tool. 

## Contact
For questions about the self-contact implementation, contact [Teo Lara](teolara@mit.edu).