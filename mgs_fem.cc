#include "mgs_fem.hh"
#include "fluid_2d.hh"
#include "tgmg.hh"

/** Initializes the class with the constants for the finite-element multigrid
 * problem.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (x_prd_,y_prd_) the periodicity in the x and y directions.
 * \param[in] zpr_ whether a pressure boundary condition is used on the right
 *                 edge of the domain.
 * \param[in] (dx,dy) the grid spacings in the x and y directions. */
mgs_fem_base::mgs_fem_base(int m_,int n_,bool x_prd_,bool y_prd_,
        bool zpr_,double dx,double dy) :
    mgs_base(m_,n_,x_prd_,y_prd_,zpr_), dydx(dy/dx), dxdy(dx/dy),
    fm(1./3.*(dxdy+dydx)), fm_inv(1.0/fm), fey(1./3.*(-2*dxdy+dydx)),
    hey(0.5*fey), fex(1./3.*(-2*dydx+dxdy)), hex(0.5*fex),
    fc(-1./6.*(dxdy+dydx)), acc(tgmg_accuracy(4*fm,mgs_accuracy)),
    pfix(new bool[mn]) {

    for(bool *pp=pfix,*pe=pfix+mn;pp<pe;pp++) *pp=false;
    if(zpr) for(int ij=m-1;ij<mn;ij+=m) pfix[ij]=true;
}

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
mgs_fem_const_rho::mgs_fem_const_rho(const fluid_2d &f) :
    mgs_fem_base(f.m_fem,f.n_fem,f.x_prd,f.y_prd,f.zpr,f.dx,f.dy),
    mg(*this,f.src,z) {
    mg.setup();
}

/** An alternative constructor that independently sets up the multigrid
 * class for testing purposes.
 * \param[in] (m_,n_) the dimensions of the grid.
 * \param[in] (dx,dy) the grid spacings in the x and y directions,
 *                    respectively.
 * \param[in] z_ a pointer to the solution array. */
mgs_fem_const_rho::mgs_fem_const_rho(int m_, int n_, bool x_prd_,
        bool y_prd_, double dx, double dy, double *src)
    : mgs_fem_base(m_,n_,x_prd_,y_prd_,false,dx,dy), mg(*this,src,z) {
    mg.setup();
}

/** Calculates the result of multiplying the matrix A by the solution vector at
 * one grid point. The diagonal term in the matrix A is omitted.
 * \param[in] i the horizontal index of the grid point.
 * \param[in] ij the overall index of the grid point. */
double mgs_fem_const_rho::mul_a(int i,int ij) {
    if(pfix[ij]) return 0;
    double *w=z+ij;

    // Since most gridpoints are interior, deal with this case
    // first
    if(i>0&&i<m-1&&ij>=m&&ij<mn-m)
        return fc*(w[-m-1]+w[-m+1]+w[m-1]+w[m+1])
               +fey*(w[-m]+w[m])+fex*(w[-1]+w[1]);

    // Compute constants and memory shifts for x-periodicity
    int sl=-1,sr=1;
    double lmu=1,rmu=1,afex=fex,afey=fey,ans;
    if(i==0) {
        sl+=m;
        if(!x_prd) lmu=0,afey=hey;
    } else if(i==m-1) {
        sr-=m;
        if(!x_prd) rmu=0,afey=hey;
    }

    // Assemble terms, taking into account y-periodicity
    if(ij>=m) ans=fc*(lmu*w[-m+sl]+rmu*w[-m+sr])+afey*w[-m];
    else if(y_prd) ans=fc*(lmu*w[mn-m+sl]+rmu*w[mn-m+sr])+afey*w[mn-m];
    else ans=0,afex=hex;
    if(ij<mn-m) ans+=fc*(lmu*w[m+sl]+rmu*w[m+sr])+afey*w[m];
    else if(y_prd) ans+=fc*(lmu*w[m-mn+sl]+rmu*w[m-mn+sr])+afey*w[m-mn];
    else afex=hex;
    return ans+afex*(lmu*w[sl]+rmu*w[sr]);
}

/** Solves the linear system using the multigrid library.
 * \param[in] f a reference to the parent fluid_2d class. */
void mgs_fem_const_rho::solve_internal(fluid_2d &f) {

    // Set the pressure boundary condition on the right edge if it is in use
    if(zpr) for(int ij=m-1;ij<mn;ij+=m) f.src[ij]=0;

    // Create an initial guess and call the multigrid solver
    initial_guess(f.time);
    if(!mg.solve_v_cycle(tp)) {
        fputs("V-cycle failed to converge in FEM problem\n",stderr);
        exit(1);
    }
}

/** The constructor sets many internal constants from the parent fluid_2d
 * class.
 * \param[in] f a reference to the parent fluid_2d class. */
mgs_fem_varying_rho::mgs_fem_varying_rho(const fluid_2d &f)
  : mgs_fem_base(f.m_fem,f.n_fem,f.x_prd,f.y_prd,f.zpr,f.dx,f.dy),
    irho(new double[mn]), mg(*this,f.src,z) {}

/** Solves the linear system using the multigrid library.
 * \param[in] f a reference to the parent fluid_2d class. */
void mgs_fem_varying_rho::solve_internal(fluid_2d &f) {

    // Create the reciprocal of the density field. Since this enters into the
    // linear system, the multigrid setup routine must be called again.
    f.set_corner_irho(irho);
    mg.setup();

    // Set the pressure boundary condition on the right edge if it is in use
    if(zpr) for(int ij=m-1; ij<mn; ij+=m) f.src[ij]=0;

    // Create an initial guess and call the multigrid solver
    initial_guess(f.time);
    if(!mg.solve_v_cycle(tp)) {
        fputs("V-cycle failed to converge in FEM problem\n",stderr);
        exit(1);
    }
}

/** Calculates the result of multiplying the matrix A by the solution vector at
 * one grid point. The diagonal term in the matrix A is omitted.
 * \param[in] i the horizontal index of the grid point.
 * \param[in] ij the overall index of the grid point. */
double mgs_fem_varying_rho::mul_a(int i,int ij) {

    // If this point is
    if(pfix[ij]) return 0;

    // Since most grid points are in the interior, deal with this case first
    double *w=z+ij;
    const double cul=cp(i,ij,ul),cdl=cp(i,ij,dl),
                 cur=cp(i,ij,ur),cdr=cp(i,ij,dr);
    if(i>0&&i<m-1&&ij>=m&&ij<mn-m)
        return cul*(w[-1]*hex+w[+m-1]*fc+w[+m]*hey)
              +cdl*(w[-1]*hex+w[-m-1]*fc+w[-m]*hey)
              +cur*(w[+1]*hex+w[+m+1]*fc+w[+m]*hey)
              +cdr*(w[+1]*hex+w[-m+1]*fc+w[-m]*hey);

    // Compute constants and memory shifts for x-periodicity
    int sl=-1,sr=1;
    double lmu=1,rmu=1,ans,
           hex_l=hex*(cdl+cul),hex_r=hex*(cdr+cur),
           hey_d=hey*(cdl+cdr),hey_u=hey*(cul+cur);
    if(i==0) {
        sl+=m;
        if(!x_prd) lmu=0,hey_d=hey*cdr,hey_u=hey*cur;
    } else if(i==m-1) {
        sr-=m;
        if(!x_prd) rmu=0,hey_d=hey*cdl,hey_u=hey*cul;
    }

    // Assemble terms, taking into account y-periodicity
    if(ij>=m) ans=fc*(lmu*w[-m+sl]*cdl+rmu*w[-m+sr]*cdr)+hey_d*w[-m];
    else if(y_prd) ans=fc*(lmu*w[mn-m+sl]*cdl+rmu*w[mn-m+sr]*cdr)+hey_d*w[mn-m];
    else ans=0,hex_l=hex*cul,hex_r=hex*cur;
    if(ij<mn-m) ans+=fc*(lmu*w[m+sl]*cul+rmu*w[m+sr]*cur)+hey_u*w[m];
    else if(y_prd) ans+=fc*(lmu*w[m-mn+sl]*cul+rmu*w[m-mn+sr]*cur)+hey_u*w[m-mn];
    else hex_l=hex*cdl,hex_r=hex*cdr;
    return ans+hex_l*lmu*w[sl]+hex_r*rmu*w[sr];
}

// Explicit instantiation
#include "tgmg.cc"
template class tgmg<mgs_fem_const_rho,double,double>;
template void tgmg_base<mgs_fem_const_rho,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<mgs_fem_const_rho,double,double>::clear_z();
template void tgmg_base<mgs_fem_const_rho,double,double>::output_res(char const*,double,double,double,double);
template double tgmg_base<mgs_fem_const_rho,double,double>::mds();
template class tgmg<mgs_fem_varying_rho,double,double>;
template void tgmg_base<mgs_fem_varying_rho,double,double>::output(char const*,double*,double,double,double,double);
template void tgmg_base<mgs_fem_varying_rho,double,double>::clear_z();
template void tgmg_base<mgs_fem_varying_rho,double,double>::output_res(char const*,double,double,double,double);
template double tgmg_base<mgs_fem_varying_rho,double,double>::mds();
