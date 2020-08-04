#ifndef SIM_TYPE_HH
#define SIM_TYPE_HH

#include <cstring>
#include <fstream>
#include <cmath>

#include "mat.hh"
#include "fields.hh"
#include "common.hh"

class fluid_2d;

/** The default scale setting the size of the boundary repulsion region. */
const double w_scale_default=3.;

/** The default size of the boundary repulsion force. */
const double w_force_default=20.;

/** \brief A structure describing the simulation geometry, initial conditions,
 * and boundary conditions. */
struct sim_type {
    /** The lower bound in the x direction. */
    const double ax;
    /** The upper bound in the x direction. */
    const double bx;
    /** The lower bound in the y direction. */
    const double ay;
    /** The upper bound in the y direction. */
    const double by;
    /** The scale setting the size of the boundary repulsion region. */
    const double w_scale;
    /** The size of the boundary repulsion force. */
    const double w_force;
    /** The periodicity in the x direction. */
    const bool x_prd;
    /** The periodicity in the y direction. */
    const bool y_prd;
    virtual ~sim_type() {}
    /** Sets up the simulation geometry to be default [-1,1]^2 box with
     * periodic boundary conditions. */
    sim_type() : ax(-1), bx(1), ay(-1), by(1), w_scale(w_scale_default),
      w_force(w_force_default), x_prd(true), y_prd(true) {}
    /** Sets up the simulation geometry.
     * \param[in] (ax_,bx_) the coordinate bounds in the x direction.
     * \param[in] (ay_,by_) the coordinate bounds in the y direction.
     * \param[in] (x_prd_,y_prd_) the periodicity in the two directions. */
    sim_type(double ax_,double bx_,double ay_,double by_,bool x_prd_,
        bool y_prd_) : ax(ax_), bx(bx_), ay(ay_), by(by_), w_scale(w_scale_default),
      w_force(w_force_default), x_prd(x_prd_), y_prd(y_prd_) {}
    /** Sets up the simulation geometry, and uses custom values for the
     * boundary repulsion parameters.
     * \param[in] (ax_,bx_) the coordinate bounds in the x direction.
     * \param[in] (ay_,by_) the coordinate bounds in the y direction.
     * \param[in] (x_prd_,y_prd_) the periodicity in the two directions.
     * \param[in] w_scale_ the scale setting the size of the boundary repulsion
     *                     region.
     * \param[in] w_force_ the boundary repulsion force. */
    sim_type(double ax_,double bx_,double ay_,double by_,bool x_prd_,
        bool y_prd_,double w_scale_,double w_force_) : ax(ax_), bx(bx_),
      ay(ay_), by(by_), w_scale(w_scale_), w_force(w_force_), x_prd(x_prd_),
      y_prd(y_prd_) {}
    /** Specifies the initial velocity of the fluid. This default function
     * initializes the velocity as zero everywhere, but it can be overridden by
     * the derived classes.
     * \param[in] (x,y) the position.
     * \param[out] (u,v) the velocity. */
    virtual void velocity(double x, double y,double &u, double &v) {
      u=v=0;
    }
    /** Applies boundary conditions to the simulation fields in the parent
     * fluid_2d class. This default function does nothing but it can be
     * overridden. */
    virtual void boundary() {}
    virtual void tracer_remap(double &x,double &y);
    /** Copies required data from the parent fluid_2d class at the start of the
     * simulation. This default function does nothing but it can be overridden.
     * \param[in] f2d a reference to the parent fluid_2d class. */
    virtual void start(fluid_2d& f2d) {}
};

/** \brief A structure describing a simulation where the initial velocity has
 * four vortices at random positions. */
struct sim_four_vortices : sim_type {
    sim_four_vortices() : sim_type(-1,1,-1,1,false,false) {}
    virtual void velocity(double x, double y,double &u, double &v);
};

/** \brief A structure describing a simulation with the four-roll velocity
 * field, a frequently-used incompressible test velocity field. */
struct sim_four_rolls : sim_type {
    sim_four_rolls() : sim_type(-1,1,-1,1,true,true) {}
    virtual void velocity(double x, double y,double &u, double &v);
};

/** \brief A class describing a simulation with initial velocity field made up
 * of pulses with different length scales, useful for performing convergence
 * tests. */
class sim_velocity_pulses : public sim_type {
    public:
        sim_velocity_pulses() : sim_type(-1,1,-1,1,true,true) {}
        virtual void velocity(double x, double y,double &u, double &v);
    private:
        /** Calculates a incompressible Gaussian pulse velocity field.
         * \param[in] s the size of the pulse.
         * \param[in] k the inverse square length scale of the pulse.
         * \param[in] (x,y) the position.
         * \param[out] (u,v) the velocity. */
        inline void pulse(double s,double k,double x,double y,double &u,double &v) {
            double f=s*exp(-k*(2-cos(x*M_PI)-cos(y*M_PI)));
            u-=sin(y*M_PI)*f;
            v+=sin(x*M_PI)*f;
        }
};

/** \brief A class describing a simulation with constant horizontal flow. */
class sim_horiz_flow : public sim_type {
  public:
        /** The horizontal flow velocity. */
        const double U;
        /** Initializes the class constants.
         * \param[in] U_ the horizontal flow velocity.
         * \param[in] (ax_,bx_) the coordinate bounds in the x direction.
         * \param[in] (ay_,by_) the coordinate bounds in the y direction. */
        sim_horiz_flow(double U_,double ax_,double bx_,double ay_,double by_)
          : sim_type(ax_,bx_,ay_,by_,false,true), U(U_) {}
        virtual void start(fluid_2d& f2d);
        virtual void boundary();
        /** Initializes the velocity to be constant horizontal flow.
         * \param[in] (x,y) the position.
         * \param[out] (u,v) the velocity. */
        virtual void velocity(double x,double y,double &u,double &v) {
          u=U;
          v=0;
        }
  private:
        /** The number of grid cells in the horizontal direction. (Copied from
         * the parent fluid_2d classs.) */
        int m;
        /** The number of grid cells in the vertical direction. (Copied from
         * the parent fluid_2d class.) */
        int n;
        /** The memory step length, taking into account ghost point allocation.
         * (Copied from the parent fluid_2d class.) */
        int ml;
        /** A pointer to the main field data structure in the parent fluid_2d
         * class. */
        field *fm;
};

#endif
