#ifndef OBJ_FIELD_HH
#define OBJ_FIELD_HH

#define _USE_MATH_DEFINES
#include <cmath>

#include "int_box.hh"
#include "fields.hh"
#include "object.hh"
#include "level++.hh"

class fluid_2d;

/** \brief A structure containing material constants for a solid. */
struct mat_const {
    /** The shear modulus. */
    double G;
    /** The solid density. */
    double rhos;
    /** The multiplier to apply to the extra viscosity in the transition
     * region. */
    double ev_trans_mult;
    /** Whether the solid is initialized to its own velocity, or inherits it
     * from the fluid. */
    bool set_velocity;
    mat_const(double G_,double rhos_,double ev_trans_mult_,bool set_velocity_=true)
        : G(G_), rhos(rhos_), ev_trans_mult(ev_trans_mult_),
            set_velocity(set_velocity_) {}
};

/** \brief A class describing all data needed to simulate a single solid. */
class obj_field {
    public:
        /** The number of fields to extrapolate, set to two for the two
         * components of the reference map. */
        static const int c=2;
        /** The memory step length, taking into account ghost point allocation. */
        const int ml;
        /** The vertical grid size, taking into account ghost point allocation. */
        const int nl;
        /** The total number of grid points, including ghost regions. */
        const int mem;
        /** The shear modulus. */
        const double G;
        /** The solid density. */
        const double rhos;
        /** The multiplier to apply to the extra viscosity in the transition
         * region. */
        const double ev_trans_mult;
        /** The transition width. */
        const double eps;
        /** The level set pinning width. */
        const double pinw;
        /** Whether the solid is initialized to its own velocity, or inherits it
         * from the fluid. */
        const bool set_velocity;
        /** An array containing the solid fields. */
        s_field* const sbase;
        /** A pointer to the (0,0) grid cell in the solid field array. */
        s_field* const sm;
        /** A class for representing the solid interface using a level set. */
        levelset ls;
        /** A pointer to the phi array held within the levelset class. */
        double* const phi_base;
        /** A pointer to the (0,0) grid cell in the phi array. */
        double* const phi;
        /** A pointer to the indicator field array within the levelset class. */
        int* const cbase;
        /** A pointer to the (0,0) grid cell in the levelset indicator field. */
        int* const cc;
        /** A pointer to the object class containing details about the object's
         * geometry and other attributes. */
        object* obj;
        /** The extra viscosity to apply in the solid. */
        double ex_visc;
        obj_field(fluid_2d &f,object* obj_,mat_const mc);
        void set_ex_visc(double tfac,double &sws_max,double &ex_visc_max,double &rhos_div_ev_min);
        /** The class destructor frees the dynamically allocated memory. */
        ~obj_field() {
            delete [] sbase;
        }
        inline void acc(int ij,double *p,double w,double x,double y) {
            s_field &s=sbase[ij];
            *p+=w*s.nX;p[1]+=w*x*s.nX;p[2]+=w*y*s.nX;
            p[3]+=w*s.nY;p[4]+=w*x*s.nY;p[5]+=w*y*s.nY;
        }
        inline void set(int ij,double *p,double al,double be,double ga) {
            sbase[ij].nX=*p*al+p[1]*be+p[2]*ga;
            sbase[ij].nY=p[3]*al+p[4]*be+p[5]*ga;
        }
        /** Extrapolates the reference map from inside the solid to the
         * neighboring grid points. */
        inline void extrapolate() {
            ls.extrapolate_fields(*this);
        }
        inline double coll_func(double x) {
            return x>eps?0:0.5-x*tf1;
        }
        /** Calculates the smoothed Heaviside transition function.
         * \param[in] x the function argument.
         * \return The function value. */
        inline double trans_func(double x) {
            return x>eps?0:trans_func_in(x);
        }
        /** Calculates the smoothed Heaviside transition function, assuming
         * that the given argument corresponds to being inside the solid or the
         * blur zone.
         * \param[in] x the function argument.
         * \return The function value. */
        inline double trans_func_in(double x) {
            return x<-eps?1:0.5-x*tf1-tf2*sin(tf3*x);
        }
        /** Calculates the derivative of the smoothed Heaviside transition
         * function multiplied by epsilon, assuming that the given argument
         * corresponds to being inside the solid or the blur zone.
         * \param[in] x the function argument.
         * \return The function value. */
        inline double tderiv_func_in(double x) {
            return x<-eps?0:0.5*(1+cos(M_PI*x/eps));
        }
        void pin();
        /** Initializes the reference map and level set values at a grid point
         * \param[in] ij the grid point index.
         * \param[in] (x,y) the position of this grid point. */
        void init(int ij,double x,double y) {
            obj->transform(x,y,sm[ij].X,sm[ij].Y);
            phi[ij]=obj->phi(sm[ij].X,sm[ij].Y);
        }
        void bound(int_box &ib);
        void mom_contrib(int ij,double &mu,double &mv,double &rho,double &sfrac);
        void rho_contrib(int ij,double &rho,double &sfrac);
    private:
        /** Constants that appear in the transition function calculations. */
        const double tf1,tf2,tf3;
};

#endif
