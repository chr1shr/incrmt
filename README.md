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
treat all contacts as self-collisions. The code here can peform a variety of simulations and carries out all examples and tests described in the following publication:



This branch modifies the existing IncRMT code to support the recent publication of Lara *et al.* listed above. The installation and build steps are unmodified.

## Example
The simple three-pronged rotor example that described in Appendix C of the
paper can be run using four threads on a 128&nbsp;&times;&nbsp;128 grid
with the following command:
```Shell
OMP_NUM_THREADS=4 ./ftest simple-spin 128
```
The code will create a directory called **sspin_128.odr** for the simulation
output. (Here, the **odr** suffix stands for **o**utput **d**irecto**r**y.)
The output directory contains files of different types:

- **w**.*<n>*, the vorticity field at frame *n*;
- **phi**.*<n>*, the level set field at frame *n*;
- **X**.*<n>* and **Y**.*<n>*, the components of the reference map at frame *n*;
- **trace**.*<n>*, the fluid tracer positions at frame *n* stored in binary format;
- **track.dat**, the position of a special tracer on an arm of the rotor that tracks
  its rotation;
- **header**, a small text file containing the number of simulation frames and
  the time interval simulated.

In Gnuplot, the vorticity field at *t*&nbsp;=&nbsp;2&pi; can be plotted using
the following commands:
```Gnuplot
set pm3d map
splot 'sspin_128.odr/w.120' matrix binary
```
If FFmpeg is installed, then the following command can be used to generate a
movie:
```Shell
./gnuplot_movie.pl -t sspin_128.odr w -10 10
```
This will generate a QuickTime movie using the H.265 codec called
**sspin_128_w.mov**. Alternatively, to just make the frames without making
a movie, the command
```Shell
./gnuplot_movie.pl -t -w sspin_128.odr w -10 10
```
can be used. This will create a directory called **sspin_128.frames** that
contains the movie frames as PNG images.

Many other types of simulation are possible with the **ftest** code, most
of which are taken from the associated publication. To see a complete list
type
```Shell
./ftest
```

## Code structure
The code is structured around several C++ classes:

- The **fluid_2d** class contains the core routines for running a simulation.
  It allocates memory for the fluid and solid fields, and contains the main
  routines for updating these.

- The **field** data structure contains all of the fields required to simulate
  the fluid in one grid cell. The **fluid_2d** class allocates a
  two-dimensional array of the **field** data structure to perform the
  simulation.

- The **object** class is a pure virtual class that specifies the geometry and
  characteristics of a solid object. Many classes are derived from this, such
  as **obj_circle** describing a solid circle, and **obj_flapper** describing
  an actuated rod that can swim via a flapping motion.

- The **obj_field** class contains all of the data required to simulate one
  solid object. It is linked to a corresponding **object** type, and also
  contains information on the object's material characteristics. It contains a
  level set array for tracking the object's boundary. It allocates a
  two-dimensional array of the **s_field** data structure, which contains all
  of the simulation fields required to represent a solid.

- The **sim_type** class is a pure virtual class that specifies to global
  initial and boundary conditions of the simulation.

The simulation method requires two linear systems to be solved during each
timestep: the marker-and-cell (MAC) projection, and the approximate projection
using the finite-element method (FEM). The code contains classes that describe
these linear systems, which are used by the TGMG library. There is a hierarchy

- **mgs_base**, containing data common across all linear systems
  - **mgs_mac**, containing data that is common for all MAC systems
    - **mgs_mac_const_rho**, the MAC system for constant density simulation,
      which allows for some significant optimization
    - **mgs_mac_varying_rho**, the MAC system for varying density simulation
  - **mgs_fem**, containing data that is common for all FEM systems
    - **mgs_fem_const_rho**, the FEM system for constant density simulation,
      which allows for some significant optimization
    - **mgs_fem_varying_rho**, the FEM system for varying density simulation

During the simulation initialization, the code checks to see whether objects
with varying density are in use, and allocates the **const_rho** or
**varying_rho** class variants accordingly.

## Known issues
This code is designed to accompany the 2020 publication by Rycroft _et al._ It
is a research code and still requires additional development to make it into a
general-purpose fluid–structure simulation tool. In particular, it has the
following known issues:

- Certain parts of the code, such as the extrapolation routines, are not
  multithreaded. This results in a loss of parallel efficiency for high numbers
  of OpenMP threads.

- In cases of extreme deformation near the object boundaries, the extrapolation
  routines may cause fictitious solid to be created within the fluid. This can
  cause the simulation to terminate prematurely.

- For simulations with multiple objects the code currently creates a separate
  **obj_field** class for every one. Each **obj_field** class contains
  globally-defined fields for the reference map, stress tensor, and level set,
  even though the object may only occupy a small region of the domain. Thus the
  memory and performance do not scale well to very large numbers (*i.e.* >100)
  of objects. This could be rectified with localized allocation of the
  simulation fields for each object.

- If an object has a high component of stress tangential to its boundary, the
  current stress-blurring mechanism can cause the nearby fluid to be
  accelerated. This is only noticeable at low viscosities and high stresses,
  and will be addressed in the future.

## Contact
For questions about the code, contact [Chris Rycroft](http://seas.harvard.edu/~chr/).

## Bibliography
1. Charles S. Peskin, *The immersed boundary method*, Acta Numerica **11**,
   479–517 (2002).
   [doi:10.1017/S0962492902000077](https://doi.org/10.1017/S0962492902000077)

2. Ken Kamrin, *Stochastic and deterministic models for dense granular flow*,
   Ph.D. thesis, Massachusetts Institute of Technology (2008).
   [DSpace](http://hdl.handle.net/1721.1/43736)

3. Ken Kamrin and Jean-Christophe Nave, *An Eulerian approach to the simulation
   of deformable solids: application to finite-strain elasticity*,
   [arXiv:0901.3799](arXiv) (2009).

4. Ken Kamrin, Chris H. Rycroft, and Jean-Christophe Nave, *Reference map
   technique for finite-strain elasticity and fluid–solid interaction*, Journal
   of the Mechanics and Physics of Solids **60**, 1952–1969 (2012).
   [doi:10.1016/j.jmps.2012.06.003](https://doi.org/10.1016/j.jmps.2012.06.003)

5. Boris Valkov, Chris H. Rycroft, and Ken Kamrin, *Eulerian method for
   multiphase interactions of soft solid bodies in fluids*, Journal of Applied
   Mechanics **82**, 041011 (2015).
   [doi:10.1115/1.4029765](https://doi.org/10.1115/1.4029765)

6. James W. Demmel, *Applied Numerical Linear Algebra*, SIAM (1997).
   [doi:10.1137/1.9781611971446](https://doi.org/10.1137/1.9781611971446)

7. William L. Briggs, Van Emden Henson, and Steve F. McCormick, *A Multigrid
   Tutorial, Second Edition*, SIAM (2000).
   [doi:10.1137/1.9780898719505](https://doi.org/10.1137/1.9780898719505)

8. James A. Sethian, *Level Set Methods and Fast Marching Methods*, Cambridge
   University Press (1999).
   [ISBN:9780521645577](https://www.cambridge.org/us/academic/subjects/mathematics/computational-science/level-set-methods-and-fast-marching-methods-evolving-interfaces-computational-geometry-fluid-mechanics-computer-vision-and-materials-science-2nd-edition?format=PB&isbn=9780521645577)
(1999).
