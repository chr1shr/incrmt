#ifndef FLUID_2D_MGS_FEM_HH
#define FLUID_2D_MGS_FEM_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "mgs_common.hh"
#include "tgmg.hh"

class fluid_2d;

/** \brief A base class for describing common features of all finite-element
 * linear systems. */
struct mgs_fem_base : public mgs_base {
    /** The ratio of vertical and horizontal grid spacings, which appears in
     * the finite-element stencils. */
    const double dydx;
    /** The ratio of horizontal and vertical grid spacings, which appears in
     * the finite-element stencils. */
    const double dxdy;
    /** The central term in the finite-element stencil. */
    const double fm;
    /** The reciprocal of the central term in the finite-element stencil, which
     * is stored separately because it is used frequently. */
    const double fm_inv;
    /** The vertical term in the finite-element stencil. */
    const double fey;
    /** Half of the vertical term in the finite-element stencil. */
    const double hey;
    /** The horizontal term in the finite-element stencil. */
    const double fex;
    /** Half of the horizontal term in the finite-element stencil. */
    const double hex;
    /** The corner term in the finite-element stencil. */
    const double fc;
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    const double acc;
    /** The array for setting essential (Dirichlet) boundary conditions. */
    bool* const pfix;
    mgs_fem_base(int m_,int n_,bool x_prd_,bool y_prd_,bool zpr_,double dx_,double dy_);
    /** The class destructor frees the dynamically allocated memory for
     * applying the Dirichlet boundary conditions. */
    ~mgs_fem_base() {
        delete [] pfix;
    }
    inline bool not_l(int i) {return x_prd||i>0;}
    inline bool not_r(int i) {return x_prd||i<m-1;}
    inline bool not_lr(int i) {return x_prd||(i>0&&i<m-1);}
    inline bool not_d(int ij) {return y_prd||ij>=m;}
    inline bool not_u(int ij) {return y_prd||ij<mn-m;}
    inline bool not_du(int ij) {return y_prd||(ij>=m&&ij<mn-m);}
    inline bool not_dl(int i,int ij) {return not_d(ij)&&not_l(i);}
    inline bool not_dr(int i,int ij) {return not_d(ij)&&not_r(i);}
    inline bool not_ul(int i,int ij) {return not_u(ij)&&not_l(i);}
    inline bool not_ur(int i,int ij) {return not_u(ij)&&not_r(i);}
};

/** \brief A class describing the finite-element linear system for a simulation
 * with constant density, which allows for some substantial simplifications. */
struct mgs_fem_const_rho : public mgs_fem_base {
    mgs_fem_const_rho(const fluid_2d &f);
    mgs_fem_const_rho(int m_,int n_,bool x_prd_,bool y_prd_,double dx,double dy,double *src);
    ~mgs_fem_const_rho() {}
    inline double a_dl(int i,int ij) {return pfix[ij]?0:(not_d(ij)&&not_l(i)?fc:0);}
    inline double a_dr(int i,int ij) {return pfix[ij]?0:(not_d(ij)&&not_r(i)?fc:0);}
    inline double a_ul(int i,int ij) {return pfix[ij]?0:(not_u(ij)&&not_l(i)?fc:0);}
    inline double a_ur(int i,int ij) {return pfix[ij]?0:(not_u(ij)&&not_r(i)?fc:0);}
    inline double a_dc(int i,int ij) {return pfix[ij]?0:(not_d(ij)?(not_lr(i)?fey:hey):0);}
    inline double a_uc(int i,int ij) {return pfix[ij]?0:(not_u(ij)?(not_lr(i)?fey:hey):0);}
    inline double a_cl(int i,int ij) {return pfix[ij]?0:(not_l(i)?(not_du(ij)?fex:hex):0);}
    inline double a_cr(int i,int ij) {return pfix[ij]?0:(not_r(i)?(not_du(ij)?fex:hex):0);}
    inline double a_cc(int i,int ij) {
        return (not_lr(i)?2:1)*(not_du(ij)?2:1)*fm;
    }
    inline double inv_cc(int i,int ij,double v) {
        return (not_lr(i)?0.5:1)*(not_du(ij)?0.5:1)*fm_inv*v;
    }
    double mul_a(int i,int ij);
    virtual void solve_internal(fluid_2d &f);
    tgmg<mgs_fem_const_rho,double,double> mg;
};

/** \brief A class describing the finite-element linear system for a simulation
 * with varying density. */
struct mgs_fem_varying_rho : public mgs_fem_base {
    /** The reciprocal of the density field. */
    double* const irho;
    mgs_fem_varying_rho(const fluid_2d &f);
    /** The class destructor frees the dynamically allocated memory for the
     * reciprocal density field. */
    ~mgs_fem_varying_rho() {
        delete [] irho;
    }
    /** The inverse rho at left,right,down,up edges. */
    enum cell_center {dl,dr,ul,ur};
    /** cp stands for the "corner pointer that points to the ij
     * position in (dl,dr,ul,ur)-direction" */
    inline double cp(int i,int ij,int dir) {
        switch(dir) {
            case dl: return irho[ij];
            case dr: return x_prd&&i==m-1?irho[ij+1-m]:irho[ij+1];
            case ul: return y_prd&&ij>=mn-m?irho[ij+m-mn]:irho[ij+m];
        }
        return irho[ij+(x_prd&&i==m-1?1-m:1)+(y_prd&&ij>=mn-m?m-mn:m)];
    }
    inline double a_dl(int i,int ij) {return pfix[ij]?0:(not_dl(i,ij)?fc*cp(i,ij,dl):0);}
    inline double a_dr(int i,int ij) {return pfix[ij]?0:(not_dr(i,ij)?fc*cp(i,ij,dr):0);}
    inline double a_ul(int i,int ij) {return pfix[ij]?0:(not_ul(i,ij)?fc*cp(i,ij,ul):0);}
    inline double a_ur(int i,int ij) {return pfix[ij]?0:(not_ur(i,ij)?fc*cp(i,ij,ur):0);}
    inline double a_dc(int i,int ij) {
        const double rho_dc=(not_l(i)?cp(i,ij,dl):0)+(not_r(i)?cp(i,ij,dr):0);
        return pfix[ij]?0:(not_d(ij)?rho_dc*hey:0);
    }
    inline double a_uc(int i,int ij) {
        const double rho_uc=(not_l(i)?cp(i,ij,ul):0)+(not_r(i)?cp(i,ij,ur):0);
        return pfix[ij]?0:(not_u(ij)?rho_uc*hey:0);
    }
    inline double a_cl(int i,int ij) {
        const double rho_cl=(not_d(ij)?cp(i,ij,dl):0)+(not_u(ij)?cp(i,ij,ul):0);
        return pfix[ij]?0:(not_l(i)?rho_cl*hex:0);
    }
    inline double a_cr(int i,int ij) {
        const double rho_cr=(not_d(ij)?cp(i,ij,dr):0)+(not_u(ij)?cp(i,ij,ur):0);
        return pfix[ij]?0:(not_r(i)?rho_cr*hex:0);
    }
    inline double a_cc(int i,int ij) {
        return (!not_d(ij)?((not_l(i)?cp(i,ij,ul):0)+(not_r(i)?cp(i,ij,ur):0))
             :(!not_u(ij))?((not_l(i)?cp(i,ij,dl):0)+(not_r(i)?cp(i,ij,dr):0))
                          :((not_l(i)?(cp(i,ij,ul)+cp(i,ij,dl)):0)+(not_r(i)?(cp(i,ij,ur)+cp(i,ij,dr)):0)))*fm;
    }
    inline double inv_cc(int i,int ij,double v) {return v/a_cc(i,ij);}
    double mul_a(int i,int ij);
    virtual void solve_internal(fluid_2d &f);
    tgmg<mgs_fem_varying_rho,double,double> mg;
};

#endif
