#ifndef BI_INTERP_HH
#define BI_INTERP_HH

#include "vec.hh"
#include "fields.hh"

/** \brief A class to perform bicubic interpolation of a velocity field or
 * reference map field. */
class bicubic_interp {
    public:
        /** The number of grid points in the x direction. */
        const int m;
        /** The number of grid points in the y direction. */
        const int n;
        /** The lower x coordinate. */
        const double ax;
        /** The lower y coordinate. */
        const double ay;
        /** The reciprocal of the x grid spacing. */
        const double xsp;
        /** The reciprocal of the y grid spacing. */
        const double ysp;
        /** A pointer the main fluid_2d fields. */
        field* const fbase;
        /** A pointer to the fields defining a solid. */
        s_field* const sbase;
        /** Temporary workspace for computing the interpolation. */
        vec a[16];
        bicubic_interp(int m_,int n_,double ax_,double bx_,double ay_,double by_,field *fbase,s_field *sbase);
        vec f(double x,double y);
        vec f_grad_f(double x,double y,vec &fx,vec &fy);
    private:
        void table_setup(int i,int j,int ij);
        void compute_x(int i,int ij,vec &c0,vec &c1,vec &c2,vec &c3);
        /** Returns the reference map field at a grid point.
         * \param[in] ij the grid point index.
         * \return The reference map. */
        inline vec sf(int ij) {
            return vec(sbase[ij].X,sbase[ij].Y);
        }
        /** Returns the velocity at a grid point.
         * \param[in] ij the grid point index.
         * \return The velocity. */
        inline vec ff(int ij) {
            return vec(fbase[ij].u,fbase[ij].v);
        }
        /** Calculates a cubic interpolant value.
         * \param[in] ap a pointer to the cubic coefficients.
         * \param[in] y the function argument. */
        inline vec yl(vec *ap,double y) {
            return *ap+y*(ap[1]+y*(ap[2]+y*ap[3]));
        }
        /** Calculates the first derivative of the cubic interpolant.
         * \param[in] ap a pointer to the cubic coefficients.
         * \param[in] y the function argument.*/
        inline vec dyl(vec *ap,double y) {
           return ap[1]+y*(2.0*ap[2]+3.0*y*ap[3]);
        }
        void grid_index(double &x,double &y);
        void fill_ad(vec *ap,vec c1,vec c2,vec c3);
        void fill_au(vec *ap,vec c0,vec c1,vec c2);
        void fill_a(vec *ap,vec c0,vec c1,vec c2,vec c3);
};

#endif
