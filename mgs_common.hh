#ifndef FLUID_2D_MGS_COMMON_HH
#define FLUID_2D_MGS_COMMON_HH

class fluid_2d;

#include "common.hh"
#include "tgmg.hh"

/** Whether to use linear interpolation in time to create the initial guess for
 * the multigrid solvers. This takes some wall clock time, but can reduce the
 * number of V-cycles required by the multigrid cycles. */
static const bool mg_lin_t_interp=false;

/** The scaling factor to apply to the multigrid convergence limit, which is
 * computed for each problem based on machine epsilon. */
const double mgs_accuracy=1e4;

/** \brief A base class describing all linear systems to be solved by the
 * multigrid library. */
struct mgs_base {
    /** The number of gridpoints in the x direction. */
    const int m;
    /** The number of gridpoints in the y direction. */
    const int n;
    /** Total number of gridpoints. */
    const int mn;
    /** The periodicity in the x direction. */
    const bool x_prd;
    /** The periodicity in the y direction. */
    const bool y_prd;
    /** A constant to enable a zero pressure boundary condition on the right
     * hand side. */
    bool zpr;
    /** The Gauss-Seidel mode to use. */
    static const char gs_mode=0;
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    double acc;
    /** The simulation time corresponding to the solution stored in the z_prev
     * array. */
    double t_prev;
    /** The simulation time corresponding to the solution in the z array. */
    double t0;
    /** The total wall clock time spent on solving this multigrid problem. */
    double wc_time;
    /** An array for storing a previous solution, which is used to initialize a
     * better starting guess for the multigrid solve, using linear
     * interpolation. */
    double* const z_prev;
    /** An array for computing the solution. */
    double* const z;
    /** The total number of multigrid solves that have been performed. This is
     * used by the initial guess function, since the code only starts using
     * linear interpolation after two solves have been performed so that there
     * are enough solutions to interpolate from. */
    int scount;
    /** A structure for predicting the number of V-cycles to perform before
     * checking the error threshold has been reached. */
    tgmg_predict tp;
    /** Returns the average number of V-cycles based on internal counters in
     * the tgmg_predict class. Those counters are reset in this call.
     * \return The average number of V-cycles. */
    inline float avg_v_cycles() {
        wc_time=0;
        return tp.avg_iters(true);
    }
    mgs_base(int m_,int n_,bool x_prd_,bool y_prd_,bool zpr_);
    virtual ~mgs_base();
    void initial_guess(double t);
    /** Solves the linear system using the multigrid library, recording the
     * time taken.
     * \param[in] f a reference to the parent fluid_2d class. */
    inline void solve(fluid_2d &f) {
        double t0=wtime();
        solve_internal(f);
        wc_time+=wtime()-t0;
    }
    virtual void solve_internal(fluid_2d &f)=0;
};

#endif
