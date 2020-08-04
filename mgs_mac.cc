#include "mgs_mac.hh"
#include "fluid_2d.hh"
#include "tgmg.hh"

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
mgs_mac_const_rho::mgs_mac_const_rho(const fluid_2d &f)
    : mgs_mac_base(f.m,f.n,f.x_prd,f.y_prd,f.zpr,f.xxsp,f.yysp,f.dt_reg),
        mg(*this,f.src,z) {
    mg.setup();
}

/** An alternative constructor that independently sets up the multigrid
 * class for testing purposes.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (dx,dy) the grid spacings in the x and y directions,
 *                    respectively.
 * \param[in] z_ a pointer to the solution array. */
mgs_mac_const_rho::mgs_mac_const_rho(int m_, int n_, bool x_prd_, bool y_prd_,
        double dx, double dy, double *src)
    : mgs_mac_base(m_,n_,x_prd_,y_prd_,false,1./(dx*dx),1./(dy*dy),1.),
        mg(*this,src,z) {
    mg.setup();
}

/** Solves the linear system using the multigrid library.
 * \param[in] f a reference to the parent fluid_2d class. */
void mgs_mac_const_rho::solve_internal(fluid_2d &f) {
    initial_guess(f.time);
    if(!mg.solve_v_cycle(tp)) {
        fputs("V-cycle failed to converge in MAC problem\n",stderr);
        exit(1);
    }
}

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
mgs_mac_varying_rho::mgs_mac_varying_rho(const fluid_2d &f)
    : mgs_mac_base(f.m,f.n,f.x_prd,f.y_prd,f.zpr,f.xxsp,f.yysp,f.dt_reg),
        irho(new double[2*mn+n+m]), irhor(irho+2*mn), irhou(irhor+n),
        mg(*this,f.src,z) {}

/** Calculates the result of multiplying the matrix A by the solution vector at
 * one grid point. The diagonal term in the matrix A is omitted.
 * \param[in] i the horizontal index of the grid point.
 * \param[in] ij the overall index of the grid point. */
double mgs_mac_varying_rho::mul_a(int i,int ij) {
    double *w=z+ij;
    const double el=ep(i,ij,l),er=ep(i,ij,r),ed=ep(i,ij,d),eu=ep(i,ij,u);
    return xxsp*(i>0?(i<m-1?el*w[-1]+er*w[1]:el*w[-1]+(x_prd?er*w[1-m]:0)):(x_prd?el*w[m-1]:0)+er*w[1])
          +yysp*(ij>=m?(ij<mn-m?ed*w[-m]+eu*w[m]:ed*w[-m]+(y_prd?eu*w[m-mn]:0)):(y_prd?ed*w[mn-m]:0)+eu*w[m]);
}

/** Solves the linear system using the multigrid library.
 * \param[in] f a reference to the parent fluid_2d class. */
void mgs_mac_varying_rho::solve_internal(fluid_2d &f) {

    // Compute the inverse edge densities. Since these enter in the linear system,
    // the multigrid setup routine must be called again.
    f.set_edge_irho(irho);
    mg.setup();

    // Create an initial guess and call the multigrid solver
    initial_guess(f.time);
    if(!mg.solve_v_cycle(tp)) {
        fputs("V-cycle failed to converge in MAC problem\n",stderr);
        exit(1);
    }
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<mgs_mac_const_rho,double,double>;
template void tgmg_base<mgs_mac_const_rho,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<mgs_mac_const_rho,double,double>::clear_z();
template void tgmg_base<mgs_mac_const_rho,double,double>::output_res(char const*,double,double,double,double);
template double tgmg_base<mgs_mac_const_rho,double,double>::mds();
template class tgmg<mgs_mac_varying_rho,double,double>;
template void tgmg_base<mgs_mac_varying_rho,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<mgs_mac_varying_rho,double,double>::clear_z();
template void tgmg_base<mgs_mac_varying_rho,double,double>::output_res(char const*,double,double,double,double);
template double tgmg_base<mgs_mac_varying_rho,double,double>::mds();
