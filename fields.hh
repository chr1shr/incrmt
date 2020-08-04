#ifndef FLUID_2D_FIELDS_HH
#define FLUID_2D_FIELDS_HH

#include <cmath>

/** Data structure for storing the fields at grid points. */
struct field {
    /** The horizontal velocity. */
    double u;
    /** The vertical velocity. */
    double v;
    /** The pressure. */
    double p;
    /** The current density. */
    double rho;
    /** Gradient of velocity, stored at the cell-center. */
    double ux,uy,vx,vy;
    /** Velocities extrapolated to edges. */
    double ud,vd,ul,vl,ur,vr,uu,vu;
    /** Temporary storage for intermediate quantities. */
    double c0,c1,c2,c3;
    /** Solid stress. */
    double s11,s12,s21,s22;
    /** Computes the speed, calculated using the Euclidean norm of the velocity
     * vector.
     * \return The speed. */
    inline double spd() {
        return sqrt(u*u+v*v);
    }
    /** Sets the field values to apply a fixed velocity boundary condition,
     * based on another field object.
     * \param[in] f the field object to use. */
    inline void fixed_bc(field &f) {
        u=-f.u;v=-f.v;
    }
    /** Sets the field values to apply a periodic velocity boundary condition,
     * based on another field object.
     * \param[in] f the field object to use. */
    inline void prd_bc(field &f) {
        u=f.u;v=f.v;
    }
    /** Calculates the velocity based on a linear extrapolation from two
     * neighboring field objects.
     * \param[in] (f1,f2) the two field values to use. */
    inline void extrap(field &f1, field &f2) {
        u = 2*f1.u-f2.u;
        v = 2*f1.v-f2.v;
    }
};

struct s_field {
    /** The reference map. */
    double X,Y;
    /** The new reference map. */
    double nX,nY;
    /** Gradient of reference map, stored at the cell-center. */
    double Xx,Xy,Yx,Yy;
    /** Reference maps extrapolated to edges. */
    double Xd,Yd,Xl,Yl,Xr,Yr,Xu,Yu;
    /** Temporary storage for intermediate quantities. */
    double r0,r1,r2,r3;
    /** Sets the reference map to apply a periodic velocity boundary condition,
     * based on another field object.
     * \param[in] f the field object to use.
     * \param[in] (delx,dely) the displacement to apply. */
    inline void prd_bc(s_field &f, const double& delx, const double& dely) {
        X=f.X+delx;Y=f.Y+dely;
    }
    inline void pre_update_ref_map(double cX,double cY) {
        nX=X+cX;nY=Y+cY;
    }
    inline void set_half_timestep() {
        X=0.5*(X+nX);
        Y=0.5*(Y+nY);
    }
    inline void set_full_timestep() {
        X=nX;Y=nY;
    }
};

#endif
