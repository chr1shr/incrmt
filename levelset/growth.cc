#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#include "level++.hh"

/** Class containing a function "velocity" for moving the level set. */
struct vel_class {
    /** Computes the normal interface velocity at a position.
     * \param[in] (x,y) the position (in terms of real coordinates) to
     *                  consider.
     * \return The normal interface velocity. */
    double velocity(double x,double y) {
        return cos(6*atan2(y,x));
    }
};

int main() {

    // Set up grid dimensions and size
    const int m=101,n=101;
    const double ax=-1,ay=-1,bx=1,by=1;
    const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);

    // Initialize the levelset class. The last two arguments control the
    // size of the band to compute.
    levelset ls(m,n,ax,ay,dx,dy,-4*(dx+dy),4*(dx+dy));

    // Initialize the phi function to be equal to some arbitrary boundary. The
    // only thing that matters is the position of the contour phi(x)=0, which
    // will be used as the starting point of the fast marching method.
    double x,y,*phi=ls.phi;
    int i,j,ij;
    for(ij=j=0;j<n;j++) {
        y=ay+dy*j;
        for(i=0;i<m;i++,ij++) {
            x=ax+dx*i;
            phi[ij]=sqrt(x*x+y*y)-0.5;
        }
    }

    // Output this initial function in the Gnuplot "matrix binary" format.
    // Within Gnuplot, type "help binary matrix" for a description of this
    // format. To plot in Gnuplot, type "splot 'phi' mat bin". To enable the
    // contour plot, do "set contour", "unset surface", "set view map", and then
    // repeat the plotting command. In addition type "help cntrparam" to tweak
    // the contour layout.
    output("phi.orig",ls);

    // Do the fast marching method to solve |grad(phi)|^2=1
    ls.build_band();

    // Move the levelset, and output the updated phi field
    vel_class vc;
    output("phi.0",ls);
    char buf[128];
    for(int k=1;k<=10;k++) {
        for(int i=0;i<4;i++) ls.move(&vc,0.01);
        sprintf(buf,"phi.%d",k);
        output(buf,ls);
    }
}
