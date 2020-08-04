#define _USE_MATH_DEFINES
#include <cmath>

#include "obj_field.hh"
#include "fluid_2d.hh"
#include "level++.hh"

/** Constructs the object field class, used to store all of the information to
 * represent a single solid object.
 * \param[in] f a reference to the parent fluid_2d class.
 * \param[in] obj_ a pointer to the object class containing details about the
 *                 object's geometry and other attributes.
 * \param[in] mc the material constants. */
obj_field::obj_field(fluid_2d &f,object* obj_,mat_const mc) :
    ml(f.ml), nl(f.nl), mem(ml*nl), G(mc.G), rhos(mc.rhos),
    ev_trans_mult(mc.ev_trans_mult), eps(f.eps), pinw(f.pinw*f.dx),
    set_velocity(mc.set_velocity), sbase(new s_field[mem]), sm(sbase+2*ml+2),
    ls(ml,nl,f.ax-1.5*f.dx,f.ay-1.5*f.dy,f.dx,f.dy,-eps-4.5*f.dx,eps+4.5*f.dx),
    phi_base(ls.phi), phi(phi_base+2*ml+2), cbase(ls.s), cc(cbase+2*ml+2),
    obj(obj_), tf1(0.5/eps), tf2(0.5/M_PI), tf3(M_PI/eps) {}

/** Calculates the extra stabilising viscosity in this solid, and updates
 * several metrics used by the parent fluid_2d class in timestep selection.
 * \param[in] tfac the prefactor to use on the extra viscosity, incorporating
 *                 grid spacing and the constant multiplier.
 * \param[in,out] sws_max the maximum shear wave speed in the fluid_2d
 *                        simulation, updated if this object has a higher
 *                        value.
 * \param[in,out] ex_visc_max the maximum extra viscosity in the fluid_2d
 *                            simulation, updated if this object has a higher
 *                            value.
 * \param[in,out] rho_div_ev_min the minimum value of the solid density divided
 *                               by the extra viscosity, updated if this object
 *                               has a lower value. */
void obj_field::set_ex_visc(double tfac,double &sws_max,double &ex_visc_max,double &rhos_div_ev_min) {

    // Calculate the shear wave speed
    double sws=sqrt(G/rhos);
    if(sws>sws_max) sws_max=sws;

    // Set extra viscosity
    ex_visc=tfac*rhos*sws;
    if(ex_visc>ex_visc_max) ex_visc_max=ex_visc;

    // Set timestep factor, taking into account that in the transition region,
    // the extra viscosity contribution is multiplied
    double r_d_ev=rhos/(ex_visc*(ev_trans_mult+1));
    if(r_d_ev<rhos_div_ev_min) rhos_div_ev_min=r_d_ev;
}

/** Computes the contribution to the momentum at a grid point from this solid,
 * plus the density and solid fraction, weighting by the transition function.
 * If the set_velocity option is false, the contributions from this solid are
 * skipped.
 * \param[in] ij the grid point to consider.
 * \param[in,out] (mu,mv) the momentum at this grid point to add the
 *                        contribution to.
 * \param[in,out] rho the density at this grid point to add the
 *                    contribution to.
 * \param[in,out] sfrac the solid fraction at this grid point to add the
 *                      contribution to. */
void obj_field::mom_contrib(int ij,double &mu,double &mv,double &rho,double &sfrac) {
    if(set_velocity) {

        // Compute the velocity and the transition function at this grid point
        double tf=trans_func(phi[ij]),uu,vv;
        obj->velocity(sm[ij].X,sm[ij].Y,uu,vv);

        // Add the contributions from this solid
        sfrac+=tf;
        rho+=tf*rhos;
        mu+=tf*rhos*uu;
        mv+=tf*rhos*vv;
    }
}

/** Computes the contribution to the density and solid fraction at a grid point
 * from this solid.
 * \param[in] ij the grid point to consider (relative to the full grid with
 *               ghost regions).
 * \param[in,out] rho the density at this grid point to add the
 *                    contribution to.
 * \param[in,out] sfrac the solid fraction at this grid point to add the
 *                      contribution to. */
void obj_field::rho_contrib(int ij,double &rho,double &sfrac) {
    double tf=trans_func(phi_base[ij]);
    sfrac+=tf;
    rho+=tf*rhos;
}

/** Pins the level set to match the reference map field. */
void obj_field::pin() {
    s_field *sp=sbase;
    for(double *p=phi_base,*pe=phi_base+mem;p<pe;p++,sp++)
        if(*p<pinw&&*p>-pinw) *p=obj->phi(sp->X,sp->Y);
    ls.build_band();
}

/** Computes the grid point bounds for this solid, and extends and integer box
 * structure to include these bounds.
 * \param[in,out] ib the integer box to extend. */
void obj_field::bound(int_box &ib) {
    for(int *cp=cbase,j=0;j<nl;j++)
        for(int i=0;i<ml;i++,cp++)
          if(*cp<4) ib.bound(i,j);
}

// Explicit instantiation
template class traverse<tra_positive,obj_field>;
