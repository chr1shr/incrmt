#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <cstring>

#include "mgs_common.hh"

/** The constructor for the base multigrid solution class sets up constants
 * and memory that are used in all multigrid problems.
 * \param[in] (m_,n_) the dimensions of the solution grid.
 * \param[in] (x_prd_,y_prd_) the periodicity in the x and y directions.
 * \param[in] zpr_ whether to apply a pressure boundary condition on the left
 *                 edge of the domain. */
mgs_base::mgs_base(int m_,int n_,bool x_prd_,bool y_prd_,bool zpr_) : m(m_),
    n(n_), mn(m*n), x_prd(x_prd_), y_prd(y_prd_), zpr(zpr_), wc_time(0.),
    z_prev(mg_lin_t_interp?new double[mn]:NULL), z(new double[mn]) {

    // For reproducibility, set the initial solution guess to zero
    for(double *zp=z,*ze=z+mn;zp<ze;zp++) *zp=0.;
}

/** The class destructor frees the dynamically allocated memory. */
mgs_base::~mgs_base() {
    delete [] z;
    if(mg_lin_t_interp) delete [] z_prev;
}

/** If the option is enabled, then this routine calculates an initial guess for
 * the solution to the linear system using linear interpolation from the
 * previous two solutions.
 * \param[in] t the current simulation time. */
void mgs_base::initial_guess(double t) {
    if(!mg_lin_t_interp) return;

    // Check the number of previous solutions
    if(scount<2) {

        // If there aren't enough solutions to do linear interpolation, then
        // just record the time and/or previous solution in preparation for
        // subsequent calls
        if(scount==1) {
            memcpy(z_prev,z,mn*sizeof(double));
            t0=t;
        } else t_prev=t;
        scount++;
    } else {

        // Linear interpolation will only be useful if the simulation time has
        // advanced. If it hasn't, then just return and do nothing.
        if(fabs(t0-t_prev)<small_number) {
            fputs("Multigrid routine called for duplicate times\n",stderr);
            return;
        }

        // Perform the linear interpolation, and update the times for the
        // previous solves
        double a=1./(t0-t_prev),b=(t-t_prev)*a,c=(t-t0)*a;
#pragma omp parallel for
        for(int ij=0;ij<mn;ij++) {
            double zz=z[ij];
            z[ij]=b*zz-c*z_prev[ij];
            z_prev[ij]=zz;
        }
        t_prev=t0;
        t0=t;
    }
}

#include "tgmg.hh"

// Explicit instantiation
#include "tgmg.cc"
template class tgmg_base<tgmg_level<double,double>,double,double>;
