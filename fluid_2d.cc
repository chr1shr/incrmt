#include <cstring>
#include <cmath>
#include <vector>

#include "mat.hh"
#include "common.hh"
#include "sim_type.hh"
#include "fluid_2d.hh"

/** The class constructor sets up constants the control the geometry and the
 * simulation, dynamically allocates memory for the fields, and calls the
 * routine to initialize the fields.
 * \param[in] (m_,n_) the number of grid points to use in the horizontal and
 *                    vertical directions.
 * \param[in] sim_t_ a reference to the sim_type class that controls the initial
 *                   conditions and boundary conditions.
 * \param[in] visc_ the fluid viscosity.
 * \param[in] rhof_ the fluid density.
 * \param[in] tw_ the transition width, measured in multiples of the minimum
 *                grid spacing.
 * \param[in] pinw_ the level set pinning width, measured in multiples of the
 *                  minimum grid spacing.
 * \param[in] ntrace_ the number of tracers to use.
 * \param[in] fflags_ the types of file output to save.
 * \param[in] filename_ the filename of the output directory. */
fluid_2d::fluid_2d(const int m_,const int n_,sim_type &sim_t_,
                   const double visc_,const double rhof_,const double tw_,
                   const double pinw_,const int ntrace_,const unsigned int fflags_,
                   const char *filename_)
    : m(m_), n(n_), mn(m_*n_), m_fem(sim_t_.x_prd?m:m+1),
    n_fem(sim_t_.y_prd?n:n+1), ml(m+4), nl(n+4), ntrace(ntrace_),
    fflags(fflags_), x_prd(sim_t_.x_prd), y_prd(sim_t_.y_prd), ax(sim_t_.ax),
    bx(sim_t_.bx), ay(sim_t_.ay), by(sim_t_.by), dx((bx-ax)/m), dy((by-ay)/n),
    xsp(1/dx), ysp(1/dy), xxsp(xsp*xsp), yysp(ysp*ysp), visc(visc_),
    rhof(rhof_), tw(tw_), pinw(pinw_), eps(tw*(dx<dy?dx:dy)), epsinv(1/eps),
    w_scale(sim_t_.w_scale), w_force(sim_t_.w_force), filename(filename_),
    fbase(new field[ml*nl]), fm(fbase+2*ml+2), tm(new double[ntrace<<2]),
    src(new double[m_fem*n_fem]), time(0.0), f_num(0), n_obj(0),
    c_obj(f2d_init_objects), zpr(false), simple_fluid_visc(false),
    solid_prd_bc(false), olist(new obj_field*[c_obj]), oe(olist),
    coll(new obj_field**[omp_get_max_threads()]), sim_t(sim_t_),
    buf(new float[m>123?m+5:128]),
    bic(ml,nl,ax-1.5*dx,bx+1.5*dx,ay-1.5*dy,by+1.5*dy,fbase,NULL) {}

/** The class destructor frees the dynamically allocated memory. */
fluid_2d::~fluid_2d() {

    // Free the thread-based collision buffer
#pragma omp parallel
    {
        delete [] coll[omp_get_thread_num()];
    }

    // Free the object-related data
    for(int i=n_obj-1;i>=0;i--) delete olist[i];
    delete [] coll;
    delete [] olist;

    // Free the multigrid solvers and temporary workspace
    delete ms_fem;
    delete ms_mac;
    delete [] buf;
    delete [] src;

    // Free the tracer memory and main grid data structure
    delete [] tm;
    delete [] fbase;
}

/** Adds an object to the simulation, creating a new obj_field class containing
 * all of the required solid fields.
 * \param[in] o a pointer to the object to add.
 * \param[in] mc the material constants for the object. */
void fluid_2d::add_object(object *o,mat_const &mc) {

    // Allocate more memory in the obj_field list and coll list if required
    if(n_obj==c_obj) {
        c_obj<<=1;
        if(c_obj>f2d_max_objects)
            fatal_error("Maximum limit of objects reached",1);
        obj_field **olist_=new obj_field*[c_obj];
        memcpy(olist_,olist,n_obj*sizeof(obj_field*));
        delete [] olist;
        olist=olist_;oe=olist+n_obj;
    }

    // Create the new obj_field class on the list
    olist[n_obj++]=new obj_field(*this,o,mc);
    oe++;
}

/** Initializes the simulation. It chooses the timestep, initializes the
 * simulation fields, sets up collision detection workspace, sets up the
 * multigrid solvers, and initializes the tracers.
 * \param[in] ev_mult the padding factor for the extra viscosity, which should
 *                    be larger than 1.
 * \param[in] dt_pad the padding factor for the timestep for the physical
 *                   terms, which should be smaller than 1.
 * \param[in] dt_ev_pad the padding factor for the timestep for the extra
 *                      viscosity terms, which should be smaller than 1. */
void fluid_2d::initialize(double ev_mult,double dt_pad,double dt_ev_pad) {

    // Set up multithreaded collision workspace
#pragma omp parallel
    {
        coll[omp_get_thread_num()]=new obj_field*[c_obj];
    }

    // Choose extra viscosity and timestep
    choose_ev_and_dt(ev_mult,dt_pad,dt_ev_pad,true);

    // Do any setup for the BCs or the objects
    sim_t.start(*this);
    for(obj_field **op=olist;op<oe;op++) (*op)->obj->start(*this,*op);

    // Set up the multigrid classes
    setup_multigrid();

    // Initialize the fields and the tracers
    init_fields();
    if(ntrace>0) init_tracers();
}

/* Sets up the multigrid classes, using the simplified versions if the density
 * is constant */
void fluid_2d::setup_multigrid() {
    if(same_rho()) {
        ms_mac=new mgs_mac_const_rho(*this);
        ms_fem=new mgs_fem_const_rho(*this);
        const_rho=true;
        puts("# Constant density");
    } else {
        ms_mac=new mgs_mac_varying_rho(*this);
        ms_fem=new mgs_fem_varying_rho(*this);
        const_rho=false;
        puts("# Non-constant density");
    }

    // Set pointers to solution arrays in multigrid classes
    smac=ms_mac->z;
    sfem=ms_fem->z;
}

/** Chooses the extra viscosity and timestep based on material parameters.
 * \param[in] ev_mult the padding factor for the extra viscosity, which should
 *                    be larger than 1.
 * \param[in] dt_pad the padding factor for the timestep for the physical
 *                   terms, which should be smaller than 1.
 * \param[in] dt_ev_pad the padding factor for the timestep for the extra
 *                      viscosity terms, which should be smaller than 1.
 * \param[in] verbose whether to print out messages to the screen. */
void fluid_2d::choose_ev_and_dt(double ev_mult,double dt_pad,double dt_ev_pad,bool verbose) {
    double hmax=dx>dy?dx:dy,hmin=dx>dy?dy:dx,
                sws_max=0,ex_visc_max=0,rhos_div_ev_min=big_number;

    // Scan the objects to determine the minimal timestep information for each,
    // based on their densities and shear moduli
    for(obj_field **op=olist;op<oe;op++)
        (*op)->set_ex_visc(ev_mult*hmax,sws_max,ex_visc_max,rhos_div_ev_min);

    // Calculate the three timestep restrictions
    int ca;
    double dt0=hmin/sws_max*dt_pad,
           dt1=0.5*rhos_div_ev_min/(xxsp+yysp)*dt_ev_pad,
           dt2=0.5*rhof/(visc*(xxsp+yysp))*dt_pad;

    // Choose the minimum of the three timestep restrictions
    if(dt0<dt1) {
        if(dt0<dt2) {ca=0;dt_reg=dt0;}
        else {ca =2;dt_reg=dt2;}
    } else {
        if(dt1<dt2) {ca=1;dt_reg=dt1;}
        else {ca=2;dt_reg=dt2;}
    }

    // Print information if requested
    if(verbose) {
        const char mno[]="",myes[]=" <-- use this";
        printf("# Ex. visc. multiplier    : %g\n"
               "# dt padding (phys)       : %g\n"
               "# dt padding (EV)         : %g\n"
               "# Max shear wave speed    : %g (L/T)\n"
               "# Max extra viscosity     : %g (M/LT)\n"
               "# Fluid viscosity         : %g (M/LT)\n"
               "# Shear wave dtmax        : %g (T)%s\n"
               "# Ex. visc. dtmax         : %g (T)%s\n"
               "# Fluid visc. dtmax       : %g (T)%s\n"
               "# Chosen dtmax            : %g (T)\n",
               ev_mult,dt_pad,dt_ev_pad,sws_max,ex_visc_max,visc,dt0,ca==0?myes:mno,
               dt1,ca==1?myes:mno,dt2,ca==2?myes:mno,dt_reg);
    }
}

/** Initializes the simulation fields. */
void fluid_2d::init_fields() {

    // Set up the reference maps for any objects
#pragma omp parallel for collapse(2)
    for(obj_field **op=olist;op<oe;op++) {
        for(int j=-2;j<n+2;j++) {
            double y=ay+dy*(0.5+j);
            for(int i=-2;i<m+2;i++)
                (*op)->init(i+ml*j,ax+dx*(0.5+i),y);
        }
    }

    // Set up the level sets for any objects
#pragma omp parallel for
    for(obj_field **op=olist;op<oe;op++)
        (*op)->ls.build_band();

    // Set up the global velocity field
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        field *fp=fm+ml*j;
        double y=ay+dy*(0.5+j);
        for(int i=0;i<m;i++,fp++) {
            double mu=0,mv=0,rho=0,sfrac=0;

            // Set the reference map and the level set for each object,
            // and gather their contributions to the velocity field
            for(obj_field **op=olist;op<oe;op++)
                (*op)->mom_contrib(i+ml*j,mu,mv,rho,sfrac);

            // Compute the velocity using a weighted average of fluid and solid
            // momenta
            if(sfrac>=1) {
                fp->u=mu/rho;
                fp->v=mv/rho;
            } else {
                double uu,vv,irho=1/((1-sfrac)*rhof+rho);
                sim_t.velocity(ax+dx*(0.5+i),y,uu,vv);
                fp->u=((1-sfrac)*rhof*uu+mu)*irho;
                fp->v=((1-sfrac)*rhof*vv+mv)*irho;
            }
            fp->p=0;
        }
        fp->p=0;
    }

    // Set the final line of the cell-cornered pressure field
    field *fp=fm+ml*n,*fe=fp+(m+1);
    while (fp<fe) (fp++)->p=0;

    // Now that the primary grid points are set up, initialize the ghost
    // points according to the boundary conditions
    set_boundaries();

    // Do an initial stress computation
    compute_stress();

    // Set the fluid/solid density
    set_rho();
}

/** Carries out the simulation for a specified time interval, periodically
 * saving the output.
 * \param[in] duration the simulation duration.
 * \param[in] frames the number of frames to save. */
void fluid_2d::solve(double duration,int frames) {
    double t0,t1,t2,adt;
    int l=timestep_select(duration/frames,adt);

    // Save header file, output the initial fields and record the initial wall
    // clock time
    save_header(duration,frames);
    if(f_num==0) write_files(0),puts("# Output frame 0");
    t0=wtime();

    // Loop over the output frames
    for(int k=1;k<=frames;k++) {

        // Perform the simulation steps
        for(int j=0;j<l;j++) step_forward(adt);

        // Output the fields
        t1=wtime();
        write_files(k+f_num);

        // Print diagnostic information. If the 8192 output flag is set, then
        // output in a machine-readable format. Otherwise, output in a
        // human-readable format.
        t2=wtime();
        if(fflags&8192) {
            printf("%d %d %.8g %.8g  %d %.8g %d %.8g\n",k+f_num,
                   l,t1-t0,t2-t1,ms_mac->tp.vcount_extra(),
                   ms_mac->wc_time,ms_fem->tp.vcount_extra(),ms_fem->wc_time);
            ms_mac->avg_v_cycles();ms_fem->avg_v_cycles();
        } else printf("# Output frame %d [%d, %.8g s, %.8g s] {MAC %.2f, FEM %.2f}\n",k+f_num,
                      l,t1-t0,t2-t1,ms_mac->avg_v_cycles(),ms_fem->avg_v_cycles());
        t0=t2;
    }
    f_num+=frames;
}

/** Steps the simulation fields forward.
 * \param[in] dt the time step to use. */
void fluid_2d::step_forward(const double& dt) {
    double ntime=time+dt;

    // Update the simulation time,and compute the tracer intermediate step
    update_tracers(1,dt);

    // Compute the edge velocities and reference maps that are needed for the
    // stable normal derivative calculation
    normal_derivatives_velocity(dt);
    normal_derivatives_refmap(dt);

    // Select which velocities to use at each boundary, using the Godunov
    // upwinding procedure
    godunov_upwinding_set<true>();
    edge_boundary_conditions();

    // Store the tangential stability information in the temporary q arrays
    compute_tangential_stability_derivatives();

    // Use the tangential-stability derivatives to compute edge velocities
    // and the reference map at a half-timestep
    compute_half_time_edge_velocities(dt);
    compute_half_time_edge_reference_map(dt);

    // Select which velocities to use at each boundary,using the Godunov
    // upwinding procedure
    godunov_upwinding_set<false>();
    edge_boundary_conditions();

    // Apply the MAC projection to correct the edge velocities to satisfy a
    // discrete divergence criterion
    time+=0.5*dt;
    mac_projection();

    // Update the reference map. This computes the full update, but leaves the
    // half timestep update in the X and Y fields for now.
    update_reference_map(dt);

    // Pin the level set to match the reference map (half timestep)
    pin();

    // Reset the ghost points according to the boundary conditions
    set_boundaries();

    // Compute stress (half timestep)
    compute_stress();

    // Compute u_star
    compute_ustar(dt);

    // Perform the finite element projection step
    time=ntime;
    finite_element_projection(dt);

    // Copy the pressure back into the main data structure
    copy_pressure();

    // Update the velocity. This also completes the reference map update.
    update_velocity(dt);

    // Pin the level set to match the reference map (end of timestep)
    pin();

    // Reset the ghost points according to the boundary conditions
    set_boundaries();

    // Compute the stress tensor using the reference map
    compute_stress();

    // Complete the tracer update
    update_tracers(2,dt);
}

/** Pins the level sets so that they are consistent with the reference map
 * fields. */
void fluid_2d::pin() {
#pragma omp parallel for
    for(obj_field **op=olist;op<oe;op++) (*op)->pin();
    if(!const_rho) set_rho();
}

/** Calculates the stress of the fluid and the reference map using the
 * nonlinear elasticity model. */
void fluid_2d::compute_stress() {

    // Perform any object-specific calculations prior to the stress computation
    for(obj_field **op=olist;op<oe;op++) (*op)->obj->pre_stress_setup(time);

#pragma omp parallel for
    for(int j=0;j<n+1;j++) {
        obj_field **collt=coll[omp_get_thread_num()];
        int ij=j*ml,k;
        field *fp=fm+ij,*fe=fp+(m+1);
        while (fp<fe) {

            // Set left edge stress, by first computing any solid stress
            // components
            k=0;
            double s1s=0,s2s=0,sfrac=0;
            for(obj_field **op=olist;op<oe;op++)
                solid_stress<true>(*op,ij,s1s,s2s,sfrac,collt,k);

            // If the solid fraction is greater than 1, the scale the stress
            // values down by this fraction. Otherwise, compute a fluid stress
            // to add in for the remaining stress.
            if(sfrac>=1) {
                fp->s11=s1s/sfrac;
                fp->s21=s2s/sfrac;
            } else {
                double s1f,s2f;
                fluid_stress<true>(fm+ij,s1f,s2f);
                fp->s11=(1-sfrac)*s1f+s1s;
                fp->s21=(1-sfrac)*s2f+s2s;
            }

            // If more than solid is present, then compute the collision stress
            // contribution. The indices of the solids involved will be in the
            // collision table.
            if(k>1) collision_stress<true>(ij,fp->s11,fp->s21,collt,k);

            // Set the down edge stress, by first computing and solid stress
            // components
            k=0;
            s1s=s2s=sfrac=0;
            for(obj_field **op=olist;op<oe;op++)
                solid_stress<false>(*op,ij,s1s,s2s,sfrac,collt,k);

            // If the solid fraction is greater than 1, the scale the stress
            // values down by this fraction. Otherwise, compute a fluid stress
            // to add in for the remaining stress.
            if(sfrac>=1) {
                fp->s12=s1s/sfrac;
                fp->s22=s2s/sfrac;
            } else {
                double s1f,s2f;
                fluid_stress<false>(fm+ij,s1f,s2f);
                fp->s12=(1-sfrac)*s1f+s1s;
                fp->s22=(1-sfrac)*s2f+s2s;
            }

            // If more than solid is present, then compute the collision stress
            // contribution. The indices of the solids involved will be in the
            // collision table.
            if(k>1) collision_stress<false>(ij,fp->s12,fp->s22,collt,k);
            fp++;ij++;
        }
    }
}

/** Computes the stress contribution at an edge, due to collisions between
 * objects.
 * \param[in] ij the grid point to consider
 * \param[in] (s1,s2) the stress components to contribute to.
 * \param[in] collt a pointer to the collision table to use, containing pointers
 *                  to the objects that are present at this location.
 * \param[in] k the number of objects on the collision table. */
template<bool left>
inline void fluid_2d::collision_stress(int ij,double &s1,double &s2,obj_field** collt,int k) {
    for(int i=0;i<k-1;i++) {
        double *phi1=collt[i]->phi+ij,p1,p1x,p1y;

        // Compute the edge-based level set and its gradient
        if(left) {
            p1=0.5*(*phi1+phi1[-1]);
            p1x=xsp*(*phi1-phi1[-1]);
            p1y=0.25*ysp*(phi1[ml]+phi1[ml-1]-phi1[-ml]-phi1[-ml-1]);
        } else {
            p1=0.5*(*phi1+phi1[-ml]);
            p1x=0.25*xsp*(phi1[1]+phi1[-ml+1]-phi1[-1]-phi1[-ml-1]);
            p1y=ysp*(*phi1-phi1[-ml]);
        }

        // Loop over all other objects, considering pairs (i,j) where i<j
        for(int j=i+1;j<k;j++) {
            double *phi2=collt[j]->phi+ij,p2,p2x,p2y;

            // Compute the edge-based level set and its gradient
            if(left) {
                p2=0.5*(*phi2+phi2[-1]);
                p2x=xsp*(*phi2-phi2[-1]);
                p2y=0.25*ysp*(phi2[ml]+phi2[ml-1]-phi2[-ml]-phi2[-ml-1]);
            } else {
                p2=0.5*(*phi2+phi2[-ml]);
                p2x=0.25*xsp*(phi2[1]+phi2[-ml+1]-phi2[-1]-phi2[-ml-1]);
                p2y=ysp*(*phi2-phi2[-ml]);
            }

            // Assemble the collision stress using the difference between the
            // two level sets
            double norx=(p2x-p1x),nory=(p2y-p1y),nn=norx*norx+nory*nory,
                   fac=4.*std::min(collt[i]->coll_func(p1),collt[j]->coll_func(p2))
                         *(collt[i]->G+collt[j]->G);
            fac=nn<small_number?0.:fac/nn;
            if(left) {
                s1-=fac*0.5*(norx*norx-nory*nory);
                s2-=fac*norx*nory;
            } else {
                s1-=fac*norx*nory;
                s2-=fac*0.5*(nory*nory-norx*norx);
            }
        }
    }
}

/** Calculates the solid stress for a given object at an edge grid point.
 * \param[in] op a pointer to the object.
 * \param[in,out] (s1s,s2s) the stress components to add to.
 * \param[in,out] sfrac the solid fraction to add to.
 * \param[in] collt a pointer to the collision table. If the solid is present,
 *                  then the object pointer is added to the table.
 * \param[in,out] the number of solids on the collision table. */
template<bool left>
inline void fluid_2d::solid_stress(obj_field *op,int ij,double &s1s,double &s2s,double &sfrac,obj_field** collt,int &k) {
    double phiv=0.5*(op->phi[ij]+op->phi[ij-(left?1:ml)]);
    if(phiv>op->eps) return;
    collt[k++]=op;
    double Xx,Yx,Xy,Yy,X,Y;
    s_field *sp=op->sm+ij;

    // Compute the edge-based reference map and deformation gradient tensor
    if(left) {
        X=0.5*(sp->X+sp[-1].X);
        Y=0.5*(sp->Y+sp[-1].Y);
        Xx=xsp*(sp->X-sp[-1].X);
        Yx=xsp*(sp->Y-sp[-1].Y);
        Xy=0.25*ysp*(sp[ml].X+sp[ml-1].X-sp[-ml].X-sp[-ml-1].X);
        Yy=0.25*ysp*(sp[ml].Y+sp[ml-1].Y-sp[-ml].Y-sp[-ml-1].Y);
    } else {
        X=0.5*(sp->X+sp[-ml].X);
        Y=0.5*(sp->Y+sp[-ml].Y);
        Xx=0.25*xsp*(sp[1].X+sp[-ml+1].X-sp[-1].X-sp[-ml-1].X);
        Yx=0.25*xsp*(sp[1].Y+sp[-ml+1].Y-sp[-1].Y-sp[-ml-1].Y);
        Xy=ysp*(sp->X-sp[-ml].X);
        Yy=ysp*(sp->Y-sp[-ml].Y);
    }

    // Compute solid stress, by first evaluating the deformation gradient
    // tensor, applying an active component F_a from the object's actuate()
    // function. The stress is then computed using F F^T.
    mat F;
    sym_mat sigma;
    F=mat(Yy,-Xy,-Yx,Xx);
    op->obj->actuate(X,Y,F);
    sigma=F.AAT();

    // Store required components, and add ex_visc contribution
    double sfac=(1/3.0)*(sigma.trace()+1),
           tf=op->trans_func_in(phiv),G=op->G*tf,
           ex_visc=op->ex_visc*tf*(1+op->ev_trans_mult*op->tderiv_func_in(phiv));
    field *fp=fm+ij;
    sfrac+=tf;
    if(left) {
        s1s+=G*(sigma.a-sfac);
        s2s+=G*sigma.b;
        s1s+=ex_visc*xsp*(fp->u-fp[-1].u);
        s2s+=ex_visc*xsp*(fp->v-fp[-1].v);
    } else {
        s1s+=G*sigma.b;
        s2s+=G*(sigma.d-sfac);
        s1s+=ex_visc*ysp*(fp->u-fp[-ml].u);
        s2s+=ex_visc*ysp*(fp->v-fp[-ml].v);
    }
}

/** Calculates the fluid stress at a cell-centered grid point.
 * \param[in] fp a pointer to the grid point.
 * \param[out] (s1f,s2f) the stress components to store.
 * \param[in] left whether to compute the left edge stress. If not, the down
 *                 edge stress is computed. */
template<bool left>
void fluid_2d::fluid_stress(field* fp,double& s1f,double &s2f) {
    if(simple_fluid_visc) {

        // If the simple fluid viscosity model is being used, then the term
        // for (nabla v)^T is omitted, since it should be zero in the
        // incompressible limit, in the bulk of the fluid.
        field *fp2=fp-(left?1:ml);
        double fac=visc*(left?xsp:ysp);
        s1f=fac*(fp->u-fp2->u);
        s2f=fac*(fp->v-fp2->v);
    } else {

        // Otherwise, the full viscosity model is used where (nabla v)^T term
        // is included. This additional term is required to accurately model
        // the stress in the blur zone.
        double ux,vx,uy,vy;
        if(left) {
            ux=xsp*(fp->u-fp[-1].u);
            vx=xsp*(fp->v-fp[-1].v);
            uy=0.25*ysp*(fp[ml].u+fp[ml-1].u-fp[-ml].u-fp[-ml-1].u);
            s1f=2*visc*ux;
            s2f=visc*(vx+uy);
        } else {
            vx=0.25*xsp*(fp[1].v+fp[-ml+1].v-fp[-1].v-fp[-ml-1].v);
            uy=ysp*(fp->u-fp[-ml].u);
            vy=ysp*(fp->v-fp[-ml].v);
            s1f=visc*(vx+uy);
            s2f=2*visc*vy;
        }
    }
}

/** Calculates the fourth-order monotonicity-limited derivative using five
 * consecutive field values.
 * \param[in] (f0,f1,f2,f3,f4) the field pointers along the input axis.
 * \return The derivative. */
double fluid_2d::mono_diff(const double& f0,const double& f1,const double& f2,const double& f3,const double& f4) {
    double tdrc=f3-f1;
    double s=(1/6.)*(4*tdrc-delta_f(f0,f1,f2)-delta_f(f2,f3,f4));
    double t=delta_lim(f1,f2,f3);
    return tdrc>0?std::min(s,t):std::max(s,-t);
}

/** Calculates the delta_f field at the field fp, which is an intermediate
 * quantity in the fourth-order monotonicity-limited derivative.
 * \param[in] f0 the value at left/down.
 * \param[in] f1 the value at center.
 * \param[in] f2 the value at right/upper.
 * \return The delta_f value. */
inline double fluid_2d::delta_f(const double& f0,const double& f1,const double& f2) {
    double drc=0.5*(f2-f0);
    return drc>0?std::min(drc,delta_lim(f0,f1,f2))
                :std::max(drc,-delta_lim(f0,f1,f2));
}

/** Calculates the delta_lim field at the field fp, which is an intermediate
 * quantity in the fourth-order monotonicity-limited derivative.
 * \param[in] f0 the value at left/down.
 * \param[in] f1 the value at center.
 * \param[in] f2 the value at right/upper.
 * \return The delta_lim value. */
inline double fluid_2d::delta_lim(const double &f0,const double &f1,
        const double &f2) {
    double drm=f1-f0,drp=f2-f1;
    return drm>0?(drp>0?2*std::min(drm,drp):0)
        :(drp>0?0:-2*std::max(drm,drp));
}

/** Applies the Godunov upwinding procedure to select field values at
 * left/right edge, based on the horizontal velocity. This function has two
 * variants: one for the tangential stability calculation, and one for the main
 * Godunov update.
 * \param[in] ij the grid point index. */
template<bool tang>
inline void fluid_2d::godunov_set_lr(int ij) {
    field *fp=fm+ij,&f=*fp,&fr=fp[1];
    if(f.ur+fr.ul>0&&f.ur>0) {
        fr.ul=f.ur;
        fr.vl=f.vr;
        for(obj_field **op=olist;op<oe;op++)
            (*op)->sm[ij+1].Xl=(*op)->sm[ij].Xr,(*op)->sm[ij+1].Yl=(*op)->sm[ij].Yr;
    } else if(f.ur+fr.ul<0&&fr.ul<0) {
        if(tang) f.ur=fr.ul;
    } else {
        if(tang) f.ur=0.5*(f.ur+fr.ul);
        fr.ul=0;
        fr.vl=0.5*(f.vr+fr.vl);
        for(obj_field **op=olist;op<oe;op++) {
            (*op)->sm[ij+1].Xl=0.5*((*op)->sm[ij].Xr+(*op)->sm[ij+1].Xl);
            (*op)->sm[ij+1].Yl=0.5*((*op)->sm[ij].Yr+(*op)->sm[ij+1].Yl);
        }
    }
}

/** Applies the Godunov upwinding procedure to select field values at down/up
 * edge, based on the vertical velocity. This function has two variants: one
 * for the tangential stability calculation, and one for the main Godunov
 * update.
 * \param[in] ij the grid point index. */
template<bool tang>
inline void fluid_2d::godunov_set_du(int ij) {
    field *fp=fm+ij,&f=*fp,&fu=fp[ml];
    if(f.vu+fu.vd>0&&f.vu>0) {
        fu.vd=f.vu;
        fu.ud=f.uu;
        for(obj_field **op=olist;op<oe;op++)
            (*op)->sm[ij+ml].Xd=(*op)->sm[ij].Xu,(*op)->sm[ij+ml].Yd=(*op)->sm[ij].Yu;
    } else if(f.vu+fu.vd<0&&fu.vd<0) {
        if(tang) f.vu=fu.vd;
    } else {
        if(tang) f.vu=0.5*(f.vu+fu.vd);
        fu.vd=0;
        fu.ud=0.5*(f.uu+fu.ud);
        for(obj_field **op=olist;op<oe;op++) {
            (*op)->sm[ij+ml].Xd=0.5*((*op)->sm[ij].Xu+(*op)->sm[ij+ml].Xd);
            (*op)->sm[ij+ml].Yd=0.5*((*op)->sm[ij].Yu+(*op)->sm[ij+ml].Yd);
        }
    }
}

/** Sets the fields in the boundaries according to the boundary conditions. */
void fluid_2d::set_boundaries() {
    if(solid_prd_bc) set_boundaries_solid();

    // Set left and right ghost values
    if(x_prd) {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-2].prd_bc(fp[m-2]);fp[-1].prd_bc(fp[m-1]);
            fp[m].prd_bc(*fp);fp[m+1].prd_bc(fp[1]);
        }
    } else {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-2].fixed_bc(fp[1]);fp[-1].fixed_bc(*fp);
            if(zpr) {
                fp[m].extrap(fp[m-1],fp[m-2]);
                fp[m+1].extrap(fp[m],fp[m-1]);
            } else {
                fp[m].fixed_bc(fp[m-1]);fp[m+1].fixed_bc(fp[m-2]);
            }
        }
    }

    // Set top and bottom ghost values
    const int tl=2*ml,g=n*ml;
    if(y_prd) {
        for(field *fp=fm-2,*fe=fp+m+4;fp<fe;fp++) {
            fp[-tl].prd_bc(fp[g-tl]);fp[-ml].prd_bc(fp[g-ml]);
            fp[g].prd_bc(*fp);fp[g+ml].prd_bc(fp[ml]);
        }
    } else {
        for(field *fp=fm-2,*fe=fp+m+4;fp<fe;fp++) {
            fp[-tl].fixed_bc(fp[ml]);fp[-ml].fixed_bc(*fp);
            fp[g].fixed_bc(fp[g-ml]);fp[g+ml].fixed_bc(fp[g-tl]);
        }
    }

    // Apply any simulation-specific boundary conditions
    sim_t.boundary();
}

/** Sets the values of the edge velocities to match the boundary conditions. */
void fluid_2d::edge_boundary_conditions() {
    if(solid_prd_bc) edge_boundary_conditions_solid_pre();

    // Set the ghost edge velocities on the left and right edges
    if(x_prd) {
        int ij=0;
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml,ij+=ml) {
            fp[-1].ur=fp[m-1].ur;
            fp[-1].vr=fp[m-1].vr;
            godunov_set_lr<true>(ij-1);
            fp[m].ul=fp->ul;
            fp[m].vl=fp->vl;
            fp[m-1].ur=fp[-1].ur;
        }
    } else {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            if(zpr) {
                fp[-1].ur=fp->ul;
                fp[m].ul=fp[m-1].ur;
                fp[m].vl=fp[m-1].vr;
            } else {
                fp->ul=fp->vl=0;
                fp[-1].ur=0;
                fp[m].ul=fp[m].vl=0;
                fp[m-1].ur=0;
            }
        }
    }

    // Set the ghost edge velocities on the bottom and top edges
    const int g=n*ml;
    if(y_prd) {
        int ij=0;
        for(field *fp=fm,*fe=fm+m;fp<fe;fp++,ij++) {
            fp[-ml].vu=fp[g-ml].vu;
            fp[-ml].uu=fp[g-ml].uu;
            godunov_set_du<true>(ij-ml);
            fp[g].vd=fp->vd;
            fp[g].ud=fp->ud;
            fp[g-ml].vu=fp[-ml].vu;
        }
    } else {
        for(field *fp=fm,*fe=fm+m;fp<fe;fp++) {
            fp->ud=fp->vd=0;
            fp[-ml].vu=0;
            fp[g].ud=fp[g].vd=0;
            fp[g-ml].vu=0;
        }
    }
    if(solid_prd_bc) edge_boundary_conditions_solid_post();
}

/** If the full solid simulation is being used, then this routine is applied
 * before the upwinding procedure, to set the extrapolated edge-based reference
 * map values based on the boundary conditions. */
void fluid_2d::edge_boundary_conditions_solid_pre() {
    const double lx=bx-ax,ly=by-ay;
    for(obj_field **op=olist;op<oe;op++) {
        for(s_field *sp=(*op)->sm,*se=sp+n*ml;sp<se;sp+=ml) {
            sp[-1].Xr=sp[m-1].Xr-lx;
            sp[-1].Yr=sp[m-1].Yr;
        }
        const int g=n*ml;
        for(s_field *sp=(*op)->sm,*se=sp+m;sp<se;sp++) {
            sp[-ml].Xu=sp[g-ml].Xu;
            sp[-ml].Yu=sp[g-ml].Yu-ly;
        }
    }
}

/** If the full solid simulation is being used, then this routine is applied
 * after the upwinding procedure, to set the upwinded edge-based reference map
 * values based on the boundary conditions. */
void fluid_2d::edge_boundary_conditions_solid_post() {
    const double lx=bx-ax,ly=by-ay;
    for(obj_field **op=olist;op<oe;op++) {
        for(s_field *sp=(*op)->sm,*se=sp+n*ml;sp<se;sp+=ml) {
            sp[m].Xl=sp->Xl+lx;
            sp[m].Yl=sp->Yl;
        }
        const int g=n*ml;
        for(s_field *sp=(*op)->sm,*se=sp+m;sp<se;sp++) {
            sp[g].Xd=sp->Xd;
            sp[g].Yd=sp->Yd+ly;
        }
    }
}

/** Sets the fields in the boundaries according to the boundary conditions. */
void fluid_2d::set_boundaries_solid() {
    const double lx=bx-ax,ly=by-ay;
    for(obj_field **op=olist;op<oe;op++) {

        // Set left and right ghost values
        for(s_field *sp=(*op)->sm,*se=sp+n*ml;sp<se;sp+=ml) {
            sp[-2].prd_bc(sp[m-2],-lx,0);sp[-1].prd_bc(sp[m-1],-lx,0);
            sp[m].prd_bc(*sp,lx,0);sp[m+1].prd_bc(sp[1],lx,0);
        }

        // Set top and bottom ghost values
        const int tl=2*ml,g=n*ml;
        for(s_field *sp=(*op)->sm-2,*se=sp+m+4;sp<se;sp++) {
            sp[-tl].prd_bc(sp[g-tl],0,-ly);sp[-ml].prd_bc(sp[g-ml],0,-ly);
            sp[g].prd_bc(*sp,0,ly);sp[g+ml].prd_bc(sp[ml],0,ly);
        }
    }
}

/** Sets boundary conditions for the FEM source term computation,taking into
 * account periodicity */
void fluid_2d::fem_source_term_conditions() {

    // Set left and right ghost values
    if(x_prd) {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-1].c0=fp[m-1].c0;
            fp[-1].c1=fp[m-1].c1;
        }
    } else if(zpr) {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-1].c0=fp[-1].u;fp[m].c0=0;
            fp[-1].c1=fp[-1].v;fp[m].c1=0;
        }
    } else {
        for(field *fp=fm,*fe=fm+n*ml;fp<fe;fp+=ml) {
            fp[-1].c0=fp[m].c0=0;
            fp[-1].c1=fp[m].c1=0;
        }
    }

    // Set top and bottom ghost values
    const int xl=x_prd?m:m+1,g=n*ml;
    if(y_prd) {
        for(field *fp=fm-1,*fe=fm+xl;fp<fe;fp++) {
            fp[-ml].c0=fp[g-ml].c0;
            fp[-ml].c1=fp[g-ml].c1;
        }
    } else {
        for(field *fp=fm-1,*fe=fm+xl;fp<fe;fp++) {
            fp[-ml].c0=fp[g].c0=0;
            fp[-ml].c1=fp[g].c1=0;
        }
    }
}

/** Do the pre-computation step to evaluate the normal derivatives of the
 * velocity.
 * \param[in] dt the size of the timestep. */
void fluid_2d::normal_derivatives_velocity(const double& dt) {
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {

            // Create references to the fields in the neighboring gridpoints
            field *fp=fm+(ml*j+i);
            field &f=*fp,&fl=fp[-1],&fr=fp[1];
            field &fd=fp[-ml],&fu=fp[ml];

            // Calculate monotonicity-limited derivatives of velocities
            f.ux=xsp*mono_diff(fp[-2].u,fl.u,f.u,fr.u,fp[2].u);
            f.vx=xsp*mono_diff(fp[-2].v,fl.v,f.v,fr.v,fp[2].v);
            f.uy=ysp*mono_diff(fp[-2*ml].u,fd.u,f.u,fu.u,fp[2*ml].u);
            f.vy=ysp*mono_diff(fp[-2*ml].v,fd.v,f.v,fu.v,fp[2*ml].v);

            // Extrapolate to left edge velocity
            f.ul=f.u-0.5*(dx+dt*f.u)*f.ux;
            f.vl=f.v-0.5*(dx+dt*f.u)*f.vx;

            // Extrapolate to right edge velocity
            f.ur=f.u+0.5*(dx-dt*f.u)*f.ux;
            f.vr=f.v+0.5*(dx-dt*f.u)*f.vx;

            // Extrapolate to bottom edge velocity
            f.ud=f.u-0.5*(dy+dt*f.v)*f.uy;
            f.vd=f.v-0.5*(dy+dt*f.v)*f.vy;

            // Extrapolate to top edge velocity
            f.uu=f.u+0.5*(dy-dt*f.v)*f.uy;
            f.vu=f.v+0.5*(dy-dt*f.v)*f.vy;
        }
    }
}

/** Do the pre-computation step to evaluate the normal derivatives of the
 * reference map.
 * \param[in] dt the timestep. */
void fluid_2d::normal_derivatives_refmap(const double& dt) {
#pragma omp parallel for collapse(2)
    for(obj_field **op=olist;op<oe;op++) {
        for(int j=0;j<n;j++) {
            for(int i=0;i<m;i++) {

                // Create references to the fields in the neighboring gridpoints
                field &f=fm[ml*j+i];
                s_field *sp=(*op)->sm+(ml*j+i);
                s_field &s=*sp,&sl=sp[-1],&sr=sp[1],&sd=sp[-ml],&su=sp[ml];

                // Calculate monotonicity-limited derivatives of (X,Y) fields
                s.Xx=xsp*mono_diff(sp[-2].X,sl.X,s.X,sr.X,sp[2].X);
                s.Yx=xsp*mono_diff(sp[-2].Y,sl.Y,s.Y,sr.Y,sp[2].Y);
                s.Xy=ysp*mono_diff(sp[-2*ml].X,sd.X,s.X,su.X,sp[2*ml].X);
                s.Yy=ysp*mono_diff(sp[-2*ml].Y,sd.Y,s.Y,su.Y,sp[2*ml].Y);

                // Extrapolate to left edge reference map
                s.Xl=s.X-0.5*(dx+dt*f.u)*s.Xx;
                s.Yl=s.Y-0.5*(dx+dt*f.u)*s.Yx;

                // Extrapolate to right edge reference map
                s.Xr=s.X+0.5*(dx-dt*f.u)*s.Xx;
                s.Yr=s.Y+0.5*(dx-dt*f.u)*s.Yx;

                // Extrapolate to bottom edge reference map
                s.Xd=s.X-0.5*(dy+dt*f.v)*s.Xy;
                s.Yd=s.Y-0.5*(dy+dt*f.v)*s.Yy;

                // Extrapolate to top edge reference map
                s.Xu=s.X+0.5*(dy-dt*f.v)*s.Xy;
                s.Yu=s.Y+0.5*(dy-dt*f.v)*s.Yy;
            }
        }
    }
}

/** Store the tangential stability information in the temporary q arrays. */
void fluid_2d::compute_tangential_stability_derivatives() {
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            field *fp=fm+(ml*j+i);
            field &f=*fp,&fl=fp[-1],&fr=fp[1],&fu=fp[ml],&fd=fp[-ml];
            double u_adv=0.5*xsp*(fr.ul+f.ul),v_adv=0.5*ysp*(fu.vd+f.vd);
            f.c0=u_adv*(f.ur-fl.ur);
            f.c1=u_adv*(fr.vl-f.vl);
            f.c2=v_adv*(fu.ud-f.ud);
            f.c3=v_adv*(f.vu-fd.vu);

            // temporary quantities for the reference map
            for(obj_field **op=olist;op<oe;op++) {
                s_field *sp=(*op)->sm+(ml*j+i),&s=*sp,&sr=sp[1],&su=sp[ml];
                s.r0=u_adv*(sr.Xl-s.Xl);
                s.r1=u_adv*(sr.Yl-s.Yl);
                s.r2=v_adv*(su.Xd-s.Xd);
                s.r3=v_adv*(su.Yd-s.Yd);
            }
        }
    }
}

/** Use the tangential-stability derivatives to compute the edge velocities
 *  at a half-timestep */
void fluid_2d::compute_half_time_edge_velocities(const double& dt) {
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        double y=ay+dy*(0.5+j);
        for(int i=0;i<m;i++) {
            const int ij=ml*j+i;
            double x=ax+dx*(0.5+i),accx,accy;

            // Create references to the fields in the neighboring gridpoints
            field *fp=fm+ij,&f=*fp;
            acceleration(fp,ij,accx,accy,x,y,true);

            // Extrapolate to left
            f.ul=f.u+0.5*(-dx*f.ux+dt*(-f.u*f.ux-f.c2+accx));//
            f.vl=f.v+0.5*(-dx*f.vx+dt*(-f.u*f.vx-f.c3+accy));//     |           |
                                                             //   --+- (uu,vu) -+--
            // Extrapolate to right                          //     |           |
            f.ur=f.u+0.5*( dx*f.ux+dt*(-f.u*f.ux-f.c2+accx));// (ul,vl) fp  (ur,vr)
            f.vr=f.v+0.5*( dx*f.vx+dt*(-f.u*f.vx-f.c3+accy));//     |           |
                                                             //   --+- (ud,vd) -+--
            // Extrapolate to bottom                         //     |           |
            f.ud=f.u+0.5*(-dy*f.uy+dt*(-f.v*f.uy-f.c0+accx));
            f.vd=f.v+0.5*(-dy*f.vy+dt*(-f.v*f.vy-f.c1+accy));

            // Extrapolate to top
            f.uu=f.u+0.5*( dy*f.uy+dt*(-f.v*f.uy-f.c0+accx));
            f.vu=f.v+0.5*( dy*f.vy+dt*(-f.v*f.vy-f.c1+accy));
        }
    }
}

/** Computes the acceleration at a grid point.
 * \param[in] fp a pointer to the grid point to consider.
 * \param[in] ij the index of the grid point to consider.
 * \param[out] (accx,accy) the components of acceleration.
 * \param[in] (x,y) the physical position of the grid point.
 * \param[in] pressure whether to include the contribution from the pressure
 *                     gradient. This is used in the predictor step, but not in
 *                     the computation of u_star since the pressure term is
 *                     dealt with in the projection. */
void fluid_2d::acceleration(field *fp,int ij,double &accx,double &accy,double x,double y,bool pressure) {
    field &f=*fp,&fr=fp[1],&fu=fp[ml];
    double accx_=0,accy_=0;

    // Calculate net force due to stress imbalance
    accx=xsp*(fr.s11-f.s11)+ysp*(fu.s12-f.s12);
    accy=xsp*(fr.s21-f.s21)+ysp*(fu.s22-f.s22);

    // Calculate pressure gradient
    if(pressure) {
        accx-=0.5*xsp*(fp[ml+1].p+fr.p-fu.p-f.p);
        accy-=0.5*ysp*(fp[ml+1].p-fr.p+fu.p-f.p);
    }

    // Apply boundary forces and additional accerelations
    for(obj_field **op=olist;op<oe;op++) {
        const double &phiv=(*op)->phi[ij];
        s_field *sp=((*op)->sm)+ij;
        if(phiv<(*op)->eps&&!solid_prd_bc) {
            wall_force(x,y,phiv,accx,accy);
            (*op)->obj->accel(x,y,sp->X,sp->Y,phiv,accx_,accy_);
        }
    }

    // All previous contributions are forces. Now divide by rho to get the
    // acceleration.
    accx/=f.rho;accy/=f.rho;

    // Add the extra accelerations, which are not normalized by 1/rho
    accx+=accx_;
    accy+=accy_;
}

/** Calculates the wall forces when object close to four edges.
 * \param[in] (x,y) the position of the current field pointing to.
 * \param[in] phiv the level set value.
 * \param[out] (fx,fy) the computed wall forces. */
void fluid_2d::wall_force(double x,double y,double phiv,double &fx,double &fy) {

    // Strength of the wall force function,w:support range;
    const double w=w_scale*eps,winv=1/w,k_rep=w_force*winv;

    // Left and right walls
    if(x<=ax+w) {
        fx+=k_rep*delta_func_wall((x-ax)*winv)
                 *(phiv>-dx?(phiv>dx?0:0.5*(1-phiv*xsp)):1);
    } else if(x>=bx-w) {
        fx-=k_rep*delta_func_wall((bx-x)*winv)
                 *(phiv>-dx?(phiv>dx?0:0.5*(1-phiv*xsp)):1);
    }

    // Bottom and top walls
    if(y<=ay+w) {
        fy+=k_rep*delta_func_wall((y-ay)*winv)
                 *(phiv>-dy?(phiv>dy?0:0.5*(1-phiv*ysp)):1);
    } else if(y>=by-w) {
        fy-=k_rep*delta_func_wall((by-y)*winv)
                 *(phiv>-dy?(phiv>dy?0:0.5*(1-phiv*ysp)):1);
    }
}

/** Use the tangential-stability derivatives to compute X,Y at the edge
 *  in a half-timestep */
void fluid_2d::compute_half_time_edge_reference_map(const double& dt) {
#pragma omp parallel for collapse(2)
    for(obj_field **op=olist;op<oe;op++)
        for(int j=0;j<n;j++) for(int i=0;i<m;i++) {

            // Create references to the fields
            s_field &s=(*op)->sm[ml*j+i];
            field &f=fm[ml*j+i];

            // Extrapolate to left
            double Xdt=dt*(-f.u*s.Xx-f.v*s.Xy),Ydt=dt*(-f.u*s.Yx-f.v*s.Yy);
            s.Xl=s.X+0.5*(-dx*s.Xx+Xdt);
            s.Yl=s.Y+0.5*(-dx*s.Yx+Ydt);

            // Extrapolate to right
            s.Xr=s.X+0.5*( dx*s.Xx+Xdt);
            s.Yr=s.Y+0.5*( dx*s.Yx+Ydt);

            // Extrapolate to bottom
            s.Xd=s.X+0.5*(-dy*s.Xy+Xdt);
            s.Yd=s.Y+0.5*(-dy*s.Yy+Ydt);

            // Extrapolate to top
            s.Xu=s.X+0.5*( dy*s.Xy+Xdt);
            s.Yu=s.Y+0.5*( dy*s.Yy+Ydt);
    }
}

/** Godunov upwinding procedure */
template<bool tang>
void fluid_2d::godunov_upwinding_set() {
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            int ij=i+ml*j;
            if(i<m-1) godunov_set_lr<tang>(ij);
            if(j<n-1) godunov_set_du<tang>(ij);
        }
    }
}

/** Computes the density field. */
void fluid_2d::set_rho() {
#pragma omp parallel for
    for(int ij=0;ij<ml*nl;ij++) {
        double rho=0,sfrac=0;

        // Evaluate the density contributions from each solid
        for(obj_field **op=olist;op<oe;op++) (*op)->rho_contrib(ij,rho,sfrac);

        // If the solid fraction is greater than 1, then renormalize the
        // density. If it's less than one, then pad the remaining fraction with
        // the fluid density.
        fbase[ij].rho=sfrac>=1?rho/sfrac:(1-sfrac)*rhof+rho;
    }
}

/** Computes the edge density field. */
void fluid_2d::set_edge_irho(double *irho) {

    // Deal with the bulk of the grid
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        double *ip=irho+2*j*m;
        for(field *fp=fm+j*ml,*fe=fp+m;fp<fe;fp++) {
            *(ip++)=0.5/(fp->rho)+0.5/(fp[-1].rho);
            *(ip++)=0.5/(fp->rho)+0.5/(fp[-ml].rho);
        }
    }

    // Deal with the right edge
    double *ip=irho+2*mn;
    for(field *fp=fm+m-1,*fe=fm+n*ml;fp<fe;fp+=ml)
        *(ip++)=0.5/(fp->rho)+0.5/(fp[1].rho);

    // Deal with the top edge
    for(field *fp=fm+(n-1)*ml,*fe=fp+m;fp<fe;fp++)
        *(ip++)=0.5/(fp->rho)+0.5/(fp[ml].rho);
}

/** Computes the inverse density field at grid corners. */
void fluid_2d::set_corner_irho(double *irho) {
#pragma omp parallel for
    for(int j=0;j<n_fem;j++) {
        field *fp=fm-ml-1+j*ml;
        for(double *ip=irho+j*m_fem,*ie=ip+m_fem;ip<ie;ip++,fp++)
            *ip=1./fp->rho;
    }
}

/** Applies the marker-and-cell (MAC) projection, correcting the edge-based
 * velocities in order to satisfy a discrete incompressibility condition on
 * each grid cell. */
void fluid_2d::mac_projection() {

    // Compute the source term for the projection, based on any local
    // divergence contribution in each cell
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        field *fp=fm+j*ml,*fe=fp+m;
        double *srp=src+j*m;
        while (fp<fe) {
            *(srp++)=xsp*(fp[1].ul-fp->ul)+ysp*(fp[ml].vd-fp->vd);
            fp++;
        }
    }

    // Solve the linear system to find the q field from which correction to the
    // velocities can be computed
    ms_mac->solve(*this);

    if(const_rho) {

        // Correct the edge velocities using the q field. This simpler version
        // assumes that the density is constant, and just uses gradients of the
        // q field.
#pragma omp parallel for
        for(int j=0;j<n;j++) {
            double *sop=smac+j*m;
            for(field *fp=fm+j*ml,*fe=fp+m;fp<fe;fp++,sop++) {
                if(j>0) fp->vd-=(*sop-sop[-m])*ysp;
                else if(y_prd) {
                    fp->vd-=(*sop-sop[mn-m])*ysp;
                    fp[n*ml].vd=fp->vd;
                }
                if(fp>fm+j*ml) fp->ul-=(*sop-sop[-1])*xsp;
                else if(x_prd) {
                    fp->ul-=(*sop-sop[m-1])*xsp;
                    fp[m].ul=fp->ul;
                }
            }
        }
    } else {

        // For the case of non-constant density, this version uses gradients of
        // the q field and the density information to correct the velocity
        double *irho=static_cast<mgs_mac_varying_rho*>(ms_mac)->irho;
#pragma omp parallel for
        for(int j=0;j<n;j++) {
            double *sop=smac+j*m,*ip=irho+2*m*j;
            for(field *fp=fm+j*ml,*fe=fp+m;fp<fe;fp++,sop++,ip+=2) {
                if(j>0) fp->vd-=ip[1]*(*sop-sop[-m])*ysp;
                else if(y_prd) {
                    fp->vd-=ip[1]*(*sop-sop[mn-m])*ysp;
                    fp[n*ml].vd=fp->vd;
                }
                if(fp>fm+j*ml) fp->ul-=*ip*(*sop-sop[-1])*xsp;
                else if(x_prd) {
                    fp->ul-=*ip*(*sop-sop[m-1])*xsp;
                    fp[m].ul=fp->ul;
                }
            }
        }
    }
}

/** Updates the reference map. It computes and stores the update at the next
 * timestep, but temporarily sets the field to have the half-timestep values in
 * order to get a more accurate acceleration computation.
 * \param[in] dt the timestep to use. */
void fluid_2d::update_reference_map(const double& dt) {

    // Compute the reference map at the end of the timestep, for all grid
    // points inside objects
    double hx=0.5*dt*xsp,hy=0.5*dt*ysp;
#pragma omp parallel for collapse(2)
    for(obj_field **op=olist;op<oe;op++)
        for(int j=0;j<n;j++) {
            int *cp=(*op)->cc+j*ml;
            s_field *sp=(*op)->sm+j*ml;
            for(field *fp=fm+j*ml,*fe=fp+m;fp<fe;fp++,cp++,sp++) {
                field &f=*fp,&fr=fp[1],&fu=fp[ml];
                s_field &s=*sp,&sr=sp[1],&su=sp[ml];
                if(*cp<4) s.pre_update_ref_map(-hx*(fr.ul+f.ul)*(sr.Xl-s.Xl)-hy*(fu.vd+f.vd)*(su.Xd-s.Xd),
                                               -hx*(fr.ul+f.ul)*(sr.Yl-s.Yl)-hy*(fu.vd+f.vd)*(su.Yd-s.Yd));
        }
    }

    // Extrapolate the reference map values from the object interiors to a
    // neighborhood outside them
#pragma omp parallel for
    for(obj_field **op=olist;op<oe;op++) (*op)->extrapolate();

    // Set the main field to contain the half-timestep reference map values
#pragma omp parallel for collapse(2)
    for(obj_field **op=olist;op<oe;op++)
        for(int j=0;j<n;j++) {
            int *cp=(*op)->cc+j*ml;
            for(s_field *sp=(*op)->sm+j*ml,*se=sp+m;sp<se;sp++,cp++)
                if(*cp<6) sp->set_half_timestep();
    }
}

/** Computes the intermediate velocity u_star, prior to performing the
 * finite-element projection.
 * \param[in] dt the timestep to use. */
void fluid_2d::compute_ustar(const double& dt) {
    double hx=0.5*dt*xsp,hy=0.5*dt*ysp;
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        double y=ay+dy*(.5+j),x=ax+0.5*dx,accx,accy;
        int ij=j*ml;
        for(field *fp=fm+ij,*fe=fp+m;fp<fe;fp++,x+=dx,ij++) {
            field &f=*fp,&fr=fp[1],&fu=fp[ml];
            acceleration(fp,ij,accx,accy,x,y,false);
            f.c0=f.u-hx*(fr.ul+f.ul)*(fr.ul-f.ul)-hy*(fu.vd+f.vd)*(fu.ud-f.ud)+dt*accx;
            f.c1=f.v-hx*(fr.ul+f.ul)*(fr.vl-f.vl)-hy*(fu.vd+f.vd)*(fu.vd-f.vd)+dt*accy;
        }
    }
}

/** Calculates the source term, then performs the finite-element projection
 * Note that pressure is stored at the corner of each cell. The index of
 * cell-center matches to the left-down corner.
 * \param[in] dt the timestep to use. */
void fluid_2d::finite_element_projection(const double& dt) {
    fem_source_term_conditions();

    // Set the source term
    double hx=0.5*dy/dt,hy=0.5*dx/dt;
#pragma omp parallel for
    for(int j=0;j<n_fem;j++) {
        double *srp=src+j*m_fem;
        for(field *fp=fm+j*ml,*fe=fp+m_fem;fp<fe;fp++) {
            *(srp++)=hx*(fp[-ml-1].c0+fp[-1].c0-fp[-ml].c0-fp->c0)
                    +hy*(fp[-ml-1].c1-fp[-1].c1+fp[-ml].c1-fp->c1);
        }
    }

    // Solve the finite-element problem
    ms_fem->solve(*this);
}

/** Copies the pressure back into the main data structure, subtracting off the
 * mean and taking into account the boundary conditions. */
void fluid_2d::copy_pressure() {
    double pavg=zpr?0.:average_pressure();
#pragma omp parallel for
    for(int j=0;j<n+1;j++) {
        // Set a pointer to the row to copy. If the domain is
        // y-periodic and this is the last line,then
        double *sop=y_prd&&j==n?sfem:sfem+j*m_fem,*soe=sop+m;
        field *fp=fm+j*ml;
        while (sop<soe) (fp++)->p=*(sop++)-pavg;
        fp->p=x_prd?fp[-m].p:*sop-pavg;
    }
}

/** Computes the average pressure that has been computed using the
 * finite-element method.
 * \return The pressure. */
double fluid_2d::average_pressure() {
    double pavg=0;

    // Loop over the finite-element grid. If a dimension is non-periodic, then
    // multiply the pressure contribution of the end points by 1/2.
#pragma omp parallel for reduction(+:pavg)
    for(int j=0;j<n_fem;j++) {
        double *sop=sfem+j*m_fem,*soe=sop+(m_fem-1);
        double prow=*(sop++)*(x_prd?1:0.5);
        while (sop<soe) prow+=*(sop++);
        prow+=*sop*(x_prd?1:0.5);
        if(!y_prd&&(j==0||j==n)) prow*=0.5;
        pavg+=prow;
    }
    return pavg*(1./mn);
}

/** Computes the velocity at the end of the timestep, using the intermediate
 * velocity and the pressure gradient. It also sets the reference map at the
 * end of the timestep, using the previously computed value.
 * \param[in] dt the timestep to use. */
void fluid_2d::update_velocity(const double& dt) {

    // Compute the velocity at the end of the timestep
#pragma omp parallel for
    for(int j=0;j<n;j++) {
        for(field *fp=fm+j*ml,*fe=fp+m;fp<fe;fp++) {
            const double rhoinv=1./fp->rho,
                         u_star=fp->c0,v_star=fp->c1;
            fp->u=u_star-0.5*dt*rhoinv*xsp*(fp[ml+1].p+fp[1].p-fp[ml].p-fp->p);
            fp->v=v_star-0.5*dt*rhoinv*ysp*(fp[ml+1].p-fp[1].p+fp[ml].p-fp->p);
        }
    }

    // Set the reference map field at the end of the timestep, using the
    // previously stored values
#pragma omp parallel for collapse(2)
    for(obj_field **op=olist;op<oe;op++)
        for(int j=0;j<n;j++)
            for(s_field *sp=(*op)->sm+j*ml,*se=sp+m;sp<se;sp++)
                sp->set_full_timestep();
}

/** Sets up the fluid tracers by initializing them at random positions. */
void fluid_2d::init_tracers() {
    for(double *tp=tm;tp<tm+(ntrace<<2);) {

        // Create a random position vector within the simulation region
        double tx=rnd(ax,bx),ty=rnd(ay,by);

        // If this position is within the solid, skip it
        if(!inside_object(tx,ty)) {*(tp++)=tx;*tp=ty;tp+=3;}
    }
}

/** Determines whether a given position vector is inside any object,based on
 * whether the bilinear interpolation of the any phi field is smaller than the
 * corresponding transition width.
 * \param[in] (tx,ty) the position vector.
 * \return True if the vector is in inside an object,false otherwise. */
bool fluid_2d::inside_object(double tx,double ty) {
    int i,j;
    double x,y,fx,fy;

    // Determine which grid cell the position vector
    tracer_cell(tx,ty,i,j,x,y);
    fx=1-x;fy=1-y;

    // Loop over all of the obj_fields and check the bilinear interpolation of
    // each phi value
    for(obj_field **op=olist;op<oe;op++) {
        double *p=(*op)->phi+(i+ml*j);
        if(fy*(*p*fx+p[1]*x)+y*(p[ml]*fx+p[ml+1]*x)<(*op)->eps) return true;
    }
    return false;
}

/** Calculates which grid cell a tracer is within.
 * \param[in] (tx,y) the tracer position.
 * \param[out] (i,j) the grid cell indices.
 * \param[out] (x,y) the tracer's fractional position within the grid cell. */
void fluid_2d::tracer_cell(double tx,double ty,int &i,int &j,double &x,double &y) {

    // Find which grid cell the tracer is in
    x=(tx-ax)*xsp-0.5;
    y=(ty-ay)*ysp-0.5;
    i=static_cast<int>(x);
    if(i<-1) i=-1;
    else if(i>m-1) i=m-1;
    j=static_cast<int>(y);
    if(j<-1) j=-1;
    else if(j>n-1) j=n-1;

    // Compute the tracer's fractional position with the grid cell
    x-=i;
    y-=j;
}

/** Computes the velocity at a point using a bilinear interpolation of the
 * velocity field.
 * \param[in] (tx,ty) the position at which to interpolate.
 * \param[out] (vx,vy) the interpolated velocity. */
inline void fluid_2d::bicubic_vel(double tx,double ty,double &vx,double &vy) {
    vec v=bic.f(tx,ty);
    vx=v.x;vy=v.y;
}

/** Computes the intermediate step of the (second-order) improved Euler
 * method.
 * \param[in] tp a pointer to the tracer information to work with. The first
 *               two doubles correspond to the tracer's current position. The
 *               intermediate step will be stored in the next two positions. */
void fluid_2d::tracer_update_1(double *tp) {
    bicubic_vel(*tp,tp[1],tp[2],tp[3]);
}

/** Completes the improved Euler method, using the previously computed
 * intermediate step.
 * \param[in] tp a pointer to the tracer information to work with. Upon
 *               completion, the first two entries will contain the updated
 *               tracer position. */
void fluid_2d::tracer_update_2(double *tp,double dt) {
    double tx=*tp+dt*tp[2],ty=tp[1]+dt*tp[3],vx,vy;
    bicubic_vel(tx,ty,vx,vy);
    *tp+=0.5*dt*(tp[2]+vx);
    tp[1]+=0.5*dt*(tp[3]+vy);
}

/** Moves the tracers according to the bilinear interpolation of the fluid
 * velocity.
 * \param[in] dt the timestep to use. */
void fluid_2d::update_tracers(int i,const double& dt) {

    if(i==1) {

        // On the first pass, the routine computes the intermediate step
        // that will be needed in the improved Euler method
        for(double *tp=tm,*te=tm+(ntrace<<2);tp<te;tp+=4)
            tracer_update_1(tp);
    } else {

        // On the second pass, the routine updates the tracer position using
        // the improved Euler method, and remaps the position if necessary
        for(double *tp=tm,*te=tm+(ntrace<<2);tp<te;tp+=4) {
            tracer_update_2(tp,dt);
            sim_t.tracer_remap(*tp,tp[1]);
        }
    }

    // Update special tracers associated with objects
    for(obj_field **op=olist;op<oe;op++) (*op)->obj->update_tracers(i,*this,dt);
}

/** Finds the minimum value of the level set fields at a given gridpoint.
 * \param[in] ij the gridpoint index to consider, relative to the base memory
 *               structure including ghost regions.
 * \return The minimum phi value. */
double fluid_2d::min_phi(int ij) {
    double mphi=big_number;
    for(obj_field **op=olist;op<oe;op++)
        if((*op)->phi_base[ij]<mphi) mphi=(*op)->phi_base[ij];
    return mphi;
}

/** Finds the obj_field corresponding to the minimum of the level set fields at
 * a given gridpoint.
 * \param[in] ij the gridpoint index to consider, relative to the base memory
 *               structure including ghost regions.
 * \return A pointer to the obj_field. */
obj_field* fluid_2d::min_phi_obj(int ij) {
    double mphi=(*olist)->phi_base[ij];
    obj_field *mop=*olist;
    for(obj_field **op=olist+1;op<oe;op++)
        if((*op)->phi_base[ij]<mphi) {mphi=(*op)->phi_base[ij];mop=*op;}
    return mop;
}

/** Finds the obj_field corresponding to the minimum of the level set fields at
 * a given corner point.
 * \param[in] ij the gridpoint index to consider, relative to the base memory
 *               structure including ghost regions. The corner is in the bottom
 *               left of this grid cell.
 * \return A pointer to the obj_field. */
obj_field* fluid_2d::min_corner_phi_obj(int ij) {
    double mphi=corner_phi((*olist)->phi_base+ij);
    obj_field *mop=*olist;
    for(obj_field **op=olist+1;op<oe;op++)
        if((*op)->phi_base[ij]<mphi) {mphi=corner_phi((*op)->phi_base+ij);mop=*op;}
    return mop;
}

// Explicit instantiation
template void fluid_2d::godunov_upwinding_set<true>();
template void fluid_2d::godunov_upwinding_set<false>();
