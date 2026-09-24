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
    /** Velocities extrapolated to edges. */
    double ux,uy,vx,vy;
    /** Gradient of velocity, stored at the cell-center. */
    double ud,vd,ul,vl,ur,vr,uu,vu;
    /** Solid stress. */
    double s11,s12,s21,s22;
    /** Temporary storage for intermediate quantities. */
    double c0,c1,c2,c3;
    /** Total phi field */
    double combinedphi;
    /** Total phi field with larger narrow band for selfcontact stress*/
    double detection_phi_contact;
    /** Gradient total phi field. */
    double gradphidownx,gradphidowny,gradphileftx,gradphilefty;   
    /** Normal direction at cell centres */ 
    double n_1,n_2;
    /** Self-Contact stress tensor components. */
    double s11c,s12c,s21c,s22c;
    /** Friction Force Field */
    double Ffx,Ffy;
    /** Friction elastic Force Field */
    double Ff_e_x,Ff_e_y;    
    /** Friction Trial Force Field */
    double F_trial_x,F_trial_y;

    // Diagnostic fields (not needed, just exist for diagnostics)
    /** ustar and vstar */
    double ustar,vstar;
    /** divergence contact x direction, this is contact repulsion force*/
    double contdivx,contdivy;
    /** alpha self-repulsion field.*/
    double alpha_down,alpha_left;
    /** lambda eigenvector self-contact direction vectors.*/
    double lambda_x_down,lambda_y_down,lambda_x_left,lambda_y_left;
    /** lambda eigenvector self-contact direction vectors, centred.*/
    double lambda_1,lambda_2;
    /** Estimated Distance Field */
    double EstDist;
    /** Estimated Distance Field Weight Function */
    double EstDistWeight;
    /** Curvature Field */
    double curvature;
    /** Difference in velocity field computed via integral */
    double relative_velocity_tan_field, relative_velocity_norm_field;
    /**Componenets of Normal Force */
    double F_norm_x, F_norm_y;
    /**Magnitude Normal Force and Friction Force */
    double F_norm_mag, F_fric_mag;

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
