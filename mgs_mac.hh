#ifndef FLUID_2D_MGS_MAC_HH
#define FLUID_2D_MGS_MAC_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "mgs_common.hh"
#include "tgmg.hh"

class fluid_2d;

/** \brief A base class describing common features of all MAC linear systems. */
struct mgs_mac_base : public mgs_base {
    /** The inverse horizontal grid spacing squared. */
    const double xxsp;
    /** The inverse vertical grid spacing squared. */
    const double yysp;
    /** Threshold on L_2 norm of residual to terminate the multigrid solve. */
    const double acc;
    /** Sets up the class constants for the MAC linear system.
     * \param[in] (m_,n_) the number of grid points in the horizontal and
     *                    vertical directions.
     * \param[in] (x_prd_,y_prd_) the periodicity in horizontal and vertical
     *                            directions.
     * \param[in] zpr_ whether to apply a pressure boundary condition on the
     *                 right edge of the domain.
     * \param[in] (xxsp_,yysp_) the square reciprocal grid spacings.
     * \param[in] dt_reg the simulation timestep. */
    mgs_mac_base(int m_,int n_,bool x_prd_,bool y_prd_,bool zpr_,double xxsp_,double yysp_,double dt_reg)
        : mgs_base(m_, n_, x_prd_, y_prd_, zpr_), xxsp(xxsp_), yysp(yysp_),
        acc(tgmg_accuracy(2*(xxsp+yysp)*dt_reg,mgs_accuracy)) {}
    inline double a_dl(int i,int ij) {return 0;}
    inline double a_dr(int i,int ij) {return 0;}
    inline double a_ul(int i,int ij) {return 0;}
    inline double a_ur(int i,int ij) {return 0;}
};

/** \brief A class describing the MAC linear system for a simulation with
 * constant density, which allows for some substantial simplifications. */
struct mgs_mac_const_rho : public mgs_mac_base {
    mgs_mac_const_rho(const fluid_2d &f);
    mgs_mac_const_rho(int m_, int n_, bool x_prd_, bool y_prd_,double dx, double dy, double *src);
    ~mgs_mac_const_rho() {}
    inline double a_dc(int i,int ij) {return y_prd||ij>=m?yysp:0;}
    inline double a_uc(int i,int ij) {return y_prd||ij<mn-m?yysp:0;}
    inline double a_cl(int i,int ij) {return x_prd||i>0?xxsp:0;}
    inline double a_cr(int i,int ij) {return x_prd||i<m-1?xxsp:0;}
    inline double a_cc(int i,int ij) {
        return (x_prd?-2*xxsp:(i>0?(i<m-1?-2*xxsp:(zpr?-3*xxsp:-xxsp)):-xxsp))
              +(y_prd?-2*yysp:(ij>=m?(ij<mn-m?-2*yysp:-yysp):-yysp));}
    inline double inv_cc(int i,int ij,double v) {return v/a_cc(i,ij);}
    inline double mul_a(int i,int ij) {
        double *w=z+ij;
        return xxsp*(i>0?(i<m-1?w[1]+w[-1]:w[-1]+(x_prd?w[1-m]:0)):(x_prd?w[m-1]:0)+w[1])
              +yysp*(ij>=m?(ij<mn-m?w[m]+w[-m]:w[-m]+(y_prd?w[m-mn]:0)):(y_prd?w[mn-m]:0)+w[m]);
    }
    virtual void solve_internal(fluid_2d &f);
    tgmg<mgs_mac_const_rho,double,double> mg;
};

/** \brief A class describing the MAC linear system for a simulation with
 * varying density. */
struct mgs_mac_varying_rho : public mgs_mac_base {
    /** An array for storing the reciprocal densities on edges. */
    double* const irho;
    /** An array for storing the reciprocal densities on the right edge of the
     * grid. */
    double* const irhor;
    /** An array for storing the reciprocal densities on the top edge of the
     * grid. */
    double* const irhou;
    mgs_mac_varying_rho(const fluid_2d &f);
    /** The class destructor frees the dynamically allocated memory for the
     * reciprocal densities. */
    ~mgs_mac_varying_rho() {
        delete [] irho;
    }
    /** The inverse rho at left, right, down, up edges. */
    enum edge {l,r,d,u};
    inline double ep(int i,int ij,int dir) {
        switch(dir) {
            case l: return irho[2*ij];
            case r: return i<m-1?irho[2*ij+2]:irhor[ij/m];
            case d: return irho[2*ij+1];
        }
        return ij<mn-m?irho[2*(ij+m)+1]:irhou[i];
    }
    inline double a_cl(int i,int ij) {return (x_prd||i>0)?xxsp*ep(i,ij,l):0;}
    inline double a_cr(int i,int ij) {return (x_prd||i<m-1)?xxsp*ep(i,ij,r):0;}
    inline double a_dc(int i,int ij) {return (y_prd||ij>=m)?yysp*ep(i,ij,d):0;}
    inline double a_uc(int i,int ij) {return (y_prd||ij<mn-m)?yysp*ep(i,ij,u):0;}
    inline double a_cc(int i,int ij) {
        const double el=ep(i,ij,l),er=ep(i,ij,r),ed=ep(i,ij,d),eu=ep(i,ij,u);
        return -xxsp*(((x_prd||(i>0))?el:0)+(x_prd?er:(i<m-1?er:(zpr?2*er:0))))
                     -yysp*(((y_prd||(ij>=m))?ed:0)+((y_prd||(ij<mn-m))?eu:0));
    }
    inline double inv_cc(int i,int ij, double v) {return v/a_cc(i, ij);}
    double mul_a(int i,int ij);
    virtual void solve_internal(fluid_2d &f);
    tgmg<mgs_mac_varying_rho,double,double> mg;
};

#endif
