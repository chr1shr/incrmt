#define _USE_MATH_DEFINES
#include <cmath>

#include "common.hh"
#include "sim_type.hh"
#include "fluid_2d.hh"

/** Remaps tracers that flow out of the simulation on the right hand side,
 * reintroducing them on the left hand side.
 * \param[in,out] (x,y) the position of the tracer, which will be remapped
 *                      during the function if necessary. */
void sim_type::tracer_remap(double &x,double &y) {
    if(x<ax) x+=(bx-ax),y=rnd(ay,by);
    if(x>bx) x-=(bx-ax),y=rnd(ay,by);
    if(y<ay) y+=(by-ay),x=rnd(ax,bx);
    if(y>by) y-=(by-ay),x=rnd(ax,bx);
}

/** Specifies the initial velocity.
 * \param[in] (x,y) the position.
 * \param[out] (u,v) the velocity. */
void sim_four_vortices::velocity(double x, double y, double &u, double &v) {
    const double tx=x-.7,ty=y-.7,sx=x+.7,sy=y+.7;

    // Prescribe four vortices
    v=+50*tx*exp(-10*(tx*tx+ty*ty))
      -50*sx*exp(-10*(sx*sx+ty*ty))
      +50*sx*exp(-10*(sx*sx+sy*sy))
      -50*tx*exp(-10*(tx*tx+sy*sy));
}

/** Specifies the initial velocity.
 * \param[in] (x,y) the position.
 * \param[out] (u,v) the velocity. */
void sim_four_rolls::velocity(double x, double y, double &u, double &v) {
    u=sin(M_PI*x)*cos(M_PI*y);
    v=-cos(M_PI*x)*sin(M_PI*y);
}

/** Specifies the initial velocity.
 * \param[in] (x,y) the position.
 * \param[out] (u,v) the velocity. */
void sim_velocity_pulses::velocity(double x, double y, double &u, double &v) {
    u=v=0;
    for(double k=2,s=1,z=-5/6.;z<1;z+=1/3.,k+=2,s=-s)
        pulse(s,k,x-z,y-z,u,v);
}

/** Specifies a boundary condition of constant horizontal flow from the left.
 */
void sim_horiz_flow::boundary() {
    for (field *fp=fm-2*ml,*fe=fm+(n+2)*ml;fp<fe;fp+=ml)
        fp[-1].u=fp[-2].u=U,
        fp[-1].v=fp[-2].v=0;
}

/** Sets additional constants for the horizontal flow simulation.
 * \param[in] f2d a reference to the parent fluid_2d class. */
void sim_horiz_flow::start(fluid_2d &f2d) {
    m=f2d.m;
    n=f2d.n;
    ml=f2d.ml;
    fm=f2d.fm;
    f2d.zpr=true;
}
