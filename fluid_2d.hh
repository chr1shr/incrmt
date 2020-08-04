#ifndef FLUID_2D_HH
#define FLUID_2D_HH

#include <cstdio>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>

#include "int_box.hh"
#include "fields.hh"
#include "object.hh"
#include "sim_type.hh"
#include "mgs_mac.hh"
#include "mgs_fem.hh"
#include "obj_field.hh"
#include "tgmg.hh"
#include "bi_interp.hh"

/** The initial memory allocation for the number of solid objects. */
const int f2d_init_objects=16;

/** The maximum memory allocation for the number of solid objects. */
const int f2d_max_objects=65536;

#ifdef _OPENMP
#include "omp.h"
#endif

/** \brief A class to carry out a 2D elasticity simulation. */
class fluid_2d {
    public:
        /** The number of grid cells in the horizontal direction. */
        const int m;
        /** The number of grid cells in the vertical direction. */
        const int n;
        /** The total number of grid cells. */
        const int mn;
        /** The number of cell-cornered points in the pressure
         * finite-element problem in the horizontal direction. */
        const int m_fem;
        /** The number of cell-cornered points in the pressure
         * finite-element problem in the vertical direction. */
        const int n_fem;
        /** The memory step length, taking into account ghost point allocation. */
        const int ml;
        /** The vertical grid size, taking into account ghost point allocation. */
        const int nl;
        /** The number of tracers. */
        const int ntrace;
        /** The file output flags. */
        const unsigned int fflags;
        /** The periodicity in the x direction. */
        const bool x_prd;
        /** The periodicity in the y direction. */
        const bool y_prd;
        /** The lower bound in the x direction. */
        const double ax;
        /** The upper bound in the x direction. */
        const double bx;
        /** The lower bound in the y direction. */
        const double ay;
        /** The upper bound in the y direction. */
        const double by;
        /** The grid spacing in the x direction. */
        const double dx;
        /** The grid spacing in the y direction. */
        const double dy;
        /** The inverse grid spacing in the x direction. */
        const double xsp;
        /** The inverse grid spacing in the y direction. */
        const double ysp;
        /** The square inverse grid spacing in the x direction. */
        const double xxsp;
        /** The square inverse grid spacing in the y direction. */
        const double yysp;
        /** The fluid viscosity. */
        const double visc;
        /** The fluid density. */
        const double rhof;
        /** The transition width multiplier. */
        const double tw;
        /** The pinning width multiplier. */
        const double pinw;
        /** The transition width. */
        const double eps;
        /** The reciprocal of the transition width. */
        const double epsinv;
        /** The scale (in multiples of eps) over which to apply the boundary
         * repulsion. */
        const double w_scale;
        /** The size of the boundary repulsion. */
        const double w_force;
        /** The filename of the output directory. */
        const char *filename;
        /** An array containing the simulation fields. */
        field* const fbase;
        /** A pointer to the (0,0) grid cell in the field array. */
        field* const fm;
        /** An array containing the tracer positions. */
        double* const tm;
        /** An array for the source terms used during the algebraic
         * multigrid solve.*/
        double* const src;
        /** The regular timestep to be used. */
        double dt_reg;
        /** The current simulation time. */
        double time;
        /** The current frame number. */
        int f_num;
        /** The total number of solid objects in the simulation. */
        int n_obj;
        /** The total memory allocation for solid objects in the simulation. */
        int c_obj;
        /** Whether the density is constant or not. */
        bool const_rho;
        /** A constant to enable a zero pressure boundary condition on the right hand
         * side. */
        bool zpr;
        /** Whether to use the simplified fluid viscosity or not. */
        bool simple_fluid_visc;
        /** Whether to apply periodic boundary conditions on the solids, mainly
         * intended for a full solid simulation. */
        bool solid_prd_bc;
        /** An pointer to the configuration of the first linear system to be solved
         * using the multigrid method. */
        mgs_mac_base *ms_mac;
        /** An pointer to the configuration of the second linear system to be solved
         * using the multigrid method. */
        mgs_fem_base *ms_fem;
        fluid_2d(const int m_,const int n_,sim_type &sim_t_,const double visc_,
                 const double rhof_,const double tw_,const double pinw_,const int ntrace_,
                 const unsigned int fflags_,const char *filename_);
        ~fluid_2d();
        void initialize(double ev_mult,double dt_pad,double dt_ev_pad);
        void setup_multigrid();
        void add_object(object *o,mat_const &mc);
        void choose_ev_and_dt(double ev_mult,double dt_pad,double dt_ev_pad,bool verbose);
        void solve(double duration,int frames);
        void save_header(double duration,int frames);
        void step_forward(const double& dt);
        void compute_stress();
        template<bool left>
        inline void solid_stress(obj_field *op,int ij,double &s1s,double &s2s,double &sfrac,obj_field **collt,int &k);
        template<bool left>
        inline void fluid_stress(field* fp,double &s1f,double &s2f);
        template<bool left>
        inline void collision_stress(int ij,double &s1,double &s2,obj_field **collt,int k);
        void init_fields();
        void init_fields_bend();
        void read_fields(const int& type);
        void init_fields_from_files();
        void write_files(int k);
        void tangential_diffusion();
        void pin();
        void init_tracers();
        void update_tracers(int i,const double& dt);
        void output(const char *prefix,const int mode,const int sn);
        void output(const char *prefix,const int mode,const int sn,int_box &ib);
        void output_tracers(const int sn);
        void tracer_update_1(double *tp);
        void tracer_update_2(double *tp,double dt);
        void momentum(double &momx,double &momy);
        void set_rho();
        void set_edge_irho(double *irho);
        void set_corner_irho(double *irho);
        void l_comparison(fluid_2d &f2d,double *err);
        /** Chooses a timestep size that is the largest value smaller than dt_reg,
         * such that a given interval length is a perfect multiple of this timestep.
         * \param[in] interval the interval length to consider.
         * \param[out] adt the timestep size.
         * \return The number of timesteps the fit into the interval. */
        inline int timestep_select(double interval,double &adt) {
            int l=static_cast<int>(interval/dt_reg)+1;
            adt=interval/l;
            return l;
        }
 private:
        /** An array of pointers to obj_field classes, which contain the fields for
         * each solid in the simulation. */
        obj_field **olist;
        /** A pointer to the end of olist array, frequently used by for-loops over
         * all of the obj_field classes. */
        obj_field **oe;
        /** An array of a pointers to obj_field classes involved in collisions at
         * grid point. */
        obj_field ***coll;
        /** A reference the to the simulation type class, containing information about the
         * geometry and initial conditions. */
        sim_type &sim_t;
        /** Temporary storage for used during the output routine. */
        float *buf;
        /** A bicubic interpolation class, used for computing bicubic interpolations
         * of the velocity for tracer updates. */
        bicubic_interp bic;
        /** An array for the solution to the MAC projection. */
        double *smac;
        /** An array for the finite-element solution. */
        double *sfem;
        void normal_derivatives_velocity(const double& dt);
        void normal_derivatives_refmap(const double& dt);
        void compute_tangential_stability_derivatives();
        void wall_force(double x,double y,double phiv,double &fx,double &fy);
        void compute_half_time_edge_velocities(const double& dt);
        void compute_half_time_edge_reference_map(const double& dt);
        template<bool tang>
        void godunov_upwinding_set();
        void mac_projection();
        void finite_element_projection(const double& dt);
        void copy_pressure();
        double average_pressure();
        void update_reference_map(const double& dt);
        void compute_ustar(const double& dt);
        void update_velocity(const double& dt);
        void set_boundaries();
        void set_boundaries_solid();
        void edge_boundary_conditions();
        void edge_boundary_conditions_solid_pre();
        void edge_boundary_conditions_solid_post();
        void fem_source_term_conditions();
        void acceleration(field *fp,int ij,double &accx,double &accy,double x,double y,bool pressure);
        double min_phi(int ij);
        obj_field* min_phi_obj(int ij);
        obj_field* min_corner_phi_obj(int ij);
        double vol_change(int ij);
        double mono_diff(const double& f0,const double& f1,const double& f2,
                         const double& f3,const double& f4);
        inline double delta_lim(const double& f0,const double& f1,const double& f2);
        inline double delta_f(const double& f0,const double& f1,const double& f2);
        template<bool tang>
        inline void godunov_set_lr(int ij);
        template<bool tang>
        inline void godunov_set_du(int ij);
        inline double delta_func(double x) {
            const double c1=0.5*epsinv,c3=M_PI*epsinv;
            return c1*(1.+cos(c3*x));
        }
        inline double delta_func_wall(const double& x) {
            return (1.-x)/x;
        }
        /** Checks whether all of the densities in the problem are the same, in
         * which case the simplified multigrid systems can be used.
         * \return Whether all the densities are the same. */
        inline bool same_rho() {
            for(obj_field **op=olist; op<oe; op++)
                if(fabs((*op)->rhos-rhof)>small_number) return false;
            return true;
        }
        /** Computes the vorticity at the lower left corner of a grid cell, based
         * on centered finite differences of the velocity.
         * \param[in] fp a pointer to the grid cell. */
        inline double vorticity(field *fp) {
            return 0.5*(xsp*(fp->v-fp[-1].v+fp[-ml].v-fp[-ml-1].v)
                       -ysp*(fp->u+fp[-1].u-fp[-ml].u-fp[-ml-1].u));
        }
        /** Computes the divergence at the lower left corner of a grid cell, based
         * on centered finite differences of the velocity.
         * \param[in] fp a pointer to the grid cell. */
        inline double divergence(field *fp) {
            return 0.5*(xsp*(fp->u-fp[-1].u+fp[-ml].u-fp[-ml-1].u)
                       +ysp*(fp->v+fp[-1].v-fp[-ml].v-fp[-ml-1].v));
        }
        /** Computes an interpolation of the phi field at the lower left corner
         * of a grid cell.
         * \param[in] phip a pointer to the phi value in the grid cell. */
        inline double corner_phi(double *phip) {
            return phip[-ml-1]+phip[-ml]+phip[-1]+*phip;
        }
        bool inside_object(double tx,double ty);
        inline void bicubic_vel(double tx,double ty,double &vx,double &vy);
        void tracer_cell(double tx,double ty,int &i,int &j,double &x,double &y);
};

#endif
