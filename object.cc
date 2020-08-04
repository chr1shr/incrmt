#define _USE_MATH_DEFINES
#include <cmath>

#include "common.hh"
#include "sim_type.hh"
#include "obj_field.hh"
#include "fluid_2d.hh"

/** Calculates the initial velocity at a position, in terms of a constant
 * velocity plus an angular velocity around the object center.
 * \param[in] (X,Y) the reference map at the position.
 * \param[out] (u,v) the velocity. */
void obj_standard::velocity(double X,double Y,double &u,double &v) {
    u=uu-omega*Y;
    v=vv+omega*X;
}

/** Computes the level set function of the circle as a function of the
 * reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return The level set function. */
double obj_circle::phi(double X,double Y) {
    return sqrt(X*X+Y*Y)-cr;
}

/** Applies a gravitational acceleration to the circle.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phiv the level set value.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_circle::accel(double x,double y,double X,double Y,
                       double phiv,double &accx,double &accy) {
    accy-=grav*op->trans_func_in(phiv);
}

/** Applies a gravitational acceleration to the square.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phiv the level set value.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_square::accel(double x,double y,double X,double Y,
                       double phiv,double &accx,double &accy) {
    accy-=grav*op->trans_func_in(phiv);
}

/** Initializes the triangle shape, by computing normal vectors to its three
 * sides.
 * \param[in] (cx_,cy_) the center of the shape.
 * \param[in] cr_ the distance from the center to the triangle vertices.
 * \param[in] theta the rotation to apply to the triangle.
 * \param[in] grav_ the gravitational acceleration to apply. */
obj_triangle::obj_triangle(double cx_,double cy_,double cr_,double theta,double grav_) :
    obj_standard(cx_,cy_,cr_*0.5), grav(grav_) {
    for(int i=0;i<3;i++) {
        nor[2*i]=cos(theta+2*M_PI/3*i);
        nor[2*i+1]=sin(theta+2*M_PI/3*i);
    }
}

/** Applies a gravitational acceleration to the triangle.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phiv the level set value.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_triangle::accel(double x,double y,double X,double Y,
                         double phiv,double &accx,double &accy) {
    accy-=grav*op->trans_func_in(phiv);
}

/** Applies a gravitational acceleration to the rounded rod.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phiv the level set value.
 * \param[in,out] (accx,accy) the acceleration vector to add to. */
void obj_rounded_rod::accel(double x,double y,double X,double Y,
                            double phiv,double &accx,double &accy) {
    accy-=grav*op->trans_func_in(phiv);
}

/** Computes the level set function of the squaree as a function of position.
 * \param[in] (x,y) the position to consider.
 * \return The level set function. */
double obj_square::phi(double X,double Y) {
    X=fabs(X)-cr;
    Y=fabs(Y)-cr;
    return X<0||Y<0?std::max(X,Y):sqrt(X*X+Y*Y);
}

/** Computes the level set function of the triangle as a function of the
 * reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return The level set function. */
double obj_triangle::phi(double X,double Y) {
    return sp(Y,-X,sp(X,Y,0)>0?(sp(X,Y,1)>0?0:2):(sp(X,Y,2)>0?1:0))-cr;
}

/** Sets up the seven-pointed star shape, by computing constants for applying
 * the pivot forces.
 * \param[in] f2d a reference to the parent fluid_2d class. */
void obj_seven_point::start_extra(fluid_2d &f2d) {

    // Compute the threshold for quickly determining if a point is inside the
    // pivot region
    piv_thresh=piv_r+op->eps;
    piv_thresh*=piv_thresh;

    // Set the anchoring stiffness so that the natural frequency is ten times
    // the numerical stability limit
    K_stiff=0.01/(f2d.dt_reg*f2d.dt_reg);
    printf("# Anchoring K_stiff       : %g (1/T^2)\n",K_stiff);
}

/** Calculates the level set function of the seven-pointed star as a function
 * of the reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return The level set function. */
double obj_seven_point::phi(double X,double Y) {
    const double fac=2*M_PI/7.;
    double theta=fabs(atan2(Y,X)),R=sqrt(X*X+Y*Y);
    theta=fabs(theta-fac*int(theta/fac+0.5));
    X=R*cos(theta)-cr;
    Y=R*sin(theta);
    return X*cos(0.25*fac)-Y*sin(0.25*fac)>0?sqrt(X*X+Y*Y):Y*cos(0.25*fac)+X*sin(0.25*fac);
}

/** Calculates the rotation of the pivot prior to a stress computation.
 * \param[in] time the current simulation time. */
void obj_seven_point::pre_stress_setup(double time) {
    double theta=th_max*(1-cos(time*omega));
    cth=cos(theta);sth=sin(theta);
}

/** Applies an acceleration due to a rotating anchor.
 * \param[in] (x,y) the position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phiv the level set value.
 * \param[in,out] (accx,acy) the acceleration vector to add to. */
void obj_seven_point::accel(double x,double y,double X,double Y,
                            double phiv,double &accx,double &accy) {
    double delx=x-cx,dely=y-cy,rsq=delx*delx+dely*dely;
    if(rsq<piv_thresh) {
        double rx=X*cth-Y*sth,ry=X*sth+Y*cth,
               K=K_stiff*op->trans_func_in(sqrt(rsq)-piv_r);
        accx-=K*(x-rx);
        accy-=K*(y-ry);
    }
}

/** Sets up the simple spinner shape, by computing constants for applying the
 * pivot forces.
 * \param[in] f2d a reference to the parent fluid_2d class. */
void obj_simple_spin::start_extra(fluid_2d &f2d) {
    open_track_file();
    setup_bicubic(f2d);

    // Compute the threshold for quickly determining if a point is inside the
    // pivot region
    piv_thresh=piv_r+op->eps;
    piv_thresh*=piv_thresh;

    // Set the anchoring stiffness so that the natural frequency is ten times
    // the numerical stability limit
    K_stiff=0.1/(f2d.dt_reg*f2d.dt_reg);
    printf("# Anchoring K_stiff       : %g (1/T^2)\n",K_stiff);
}

/** Calculates the level set function of the simple spinner as a function of
 * the reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return The level set function. */
double obj_simple_spin::phi(double X,double Y) {
    double th=atan2(Y,X),c=cos(3*th);
    return sqrt(X*X+Y*Y)-c0-c*(c1+c2*c);
}

/** Calculates the rotation of the pivot prior to a stress computation.
 * \param[in] time the current simulation time. */
void obj_simple_spin::pre_stress_setup(double time) {
    double ang=time*omega,theta=ang>2*M_PI?0:th_max*(1-cos(ang));
    cth=cos(theta);sth=sin(theta);
}

/** Applies an acceleration due to a rotating anchor.
 * \param[in] (x,y) the position.
 * \param[in] (X,Y) the reference map at this position.
 * \param[in] phiv the level set value.
 * \param[in,out] (accx,acy) the acceleration vector to add to. */
void obj_simple_spin::accel(double x,double y,double X,double Y,
                            double phiv,double &accx,double &accy) {
    double rsq=x*x+y*y;
    if(rsq<piv_thresh) {
        double rx=X*cth-Y*sth,ry=X*sth+Y*cth,
               K=K_stiff*op->trans_func_in(sqrt(rsq)-piv_r);
        accx-=K*(x-rx);
        accy-=K*(y-ry);
    }
}

/** Writes the position of the tracer to a file, along with its rotation angle
 * (taking into account winding).
 * \param[in] time the current simulation time. */
void obj_simple_spin::file_output(double time) {
    double ang=atan2(tm[1],*tm);
    correct_tracers();

    // Check for a change in winding number
    int q=fabs(ang)>0.5*M_PI?(ang>0?1:-1):0;
    if(q==1&&prev_q==-1) wind--;
    else if(q==-1&&prev_q==1) wind++;
    prev_q=q;

    // Assemble the output buffer and save it to file
    *buf=time;
    buf[1]=*tm;
    buf[2]=tm[1];
    buf[3]=ang+2*M_PI*wind;
    fwrite(buf,sizeof(double),4,fp);
    count=1;
}

/** Sets up the piston, by computing constants for applying the pivot forces,
 * and opening the tracer tracking file.
 * \param[in] f2d a reference to the parent fluid_2d class. */
void obj_piston::start_extra(fluid_2d &f2d) {
    open_track_file();
    setup_bicubic(f2d);

    // Compute the threshold for quickly determining if a point is inside the
    // anchored region
    anch_thresh=anch_r+op->eps;
    anch_thresh*=anch_thresh;

    // Set the anchoring stiffness so that the natural frequency is ten times
    // the numerical stability limit
    K_stiff=0.01/(f2d.dt_reg*f2d.dt_reg);
    printf("# Anchoring K_stiff       : %g (1/T^2)\n",K_stiff);
}

/** Calculates the level set function of the piston as a function of the
 * reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return The level set function. */
double obj_piston::phi(double X,double Y) {
    double gx=fabs(X)-lx,gy=fabs(Y)-ly;
    return gx>gy?gx:gy;
}

/** Sets the piston anchor y position as a function of time.
 * \param[in] time the current simulation time. */
void obj_piston::pre_stress_setup(double time) {
    time*=idur;
    cy=y_start+(time>=1?len:len*time*time*(3-2*time));
}

/** Applies a rotating anchor force.
 * \param[in] (x,y) the current position.
 * \param[in] (X,Y) the reference map at the position.
 * \param[in] phiv the level set value at the position.
 * \param[out] (accx,accy) the acceleration due to the anchor. */
void obj_piston::accel(double x,double y,double X,double Y,
                       double phiv,double &accx,double &accy) {
    double delx=x-cx,dely=y-cy;
    double rsq=delx*delx+dely*dely;
    if(rsq<anch_thresh) {
        double K=K_stiff*op->trans_func_in(sqrt(rsq)-anch_r);
        accx-=K*(x-X);
        accy-=K*(dely-Y);
    }
}

/** Computes the level set function of the rounded rod as a function of the
 * reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return The level set function. */
double obj_rounded_rod::phi(double X,double Y) {
    double fY=fabs(Y)-cl;
    return (fY<0?fabs(X):sqrt(fY*fY+X*X))-cr;
}

/** Computes the level set function of the flapper as a function of the
 * reference map.
 * \param[in] (X,Y) the reference map position to consider.
 * \return The level set function. */
double obj_flapper::phi(double X,double Y) {
    double fX=fabs(X)-cl;
    return (fX<0?fabs(Y):sqrt(fX*fX+Y*Y))-cr;
}

/** Computes the current contraction amount of the flapper.
 * \param[in] time the current simulation time. */
void obj_flapper::pre_stress_setup(double time) {
    double a=sin(omega*time);
    a*=a*a;a*=a;a*=a;
    con=-amp*a;
}

/** Modifies the deformation gradient tensor to apply the flapper swimming
 * motion.
 * \param[in] (X,Y) the reference map position to consider.
 * \param[in,out] F the deformation gradient tensor to modify. */
void obj_flapper::actuate(double X,double Y,mat &F) {
    double fX=fabs(X)-al,ph=(fX<0?fabs(Y):sqrt(fX*fX+Y*Y))-ar;
    if(ph<op->eps) {
        double fac=exp(con*Y*op->trans_func_in(ph));
        F.a*=fac;F.b*=1./fac;
        F.c*=fac;F.d*=1./fac;
    }
}

/** Opens the file for storing tracers, and sets up the bicubic interpolation
 * classes for correcting the tracers.
 * \param[in] f2d a reference to the parent fluid_2d class. */
void obj_track::start_extra(fluid_2d &f2d) {
    open_track_file();
    setup_bicubic(f2d);
}

/** Creates the bicubic interpolation class of the object reference map, for
 * periodically correcting the tracers to match.
 * \param[in] f2d a reference to the parent fluid_2d class. */
void obj_track::setup_bicubic(fluid_2d &f2d) {
    double ddx=1.5*f2d.dx,ddy=1.5*f2d.dy;
    bic=new bicubic_interp(f2d.ml,f2d.nl,f2d.ax-ddx,f2d.bx+ddx,f2d.ay-ddy,f2d.by+ddy,NULL,op->sbase);
}

/** Opens the file for storing tracers. */
void obj_track::open_track_file() {
    char *fn=new char[strlen(filename)+16];
    sprintf(fn,"%s.odr/track.dat",filename);
    fp=safe_fopen(fn,"wb");
}

/** Writes the tracer positions to the output file. Output to the file is
 * flushed, so that it can analyzed on-the-fly.
 * \param[in] k the frame output number.
 * \param[in] time the current simulation time. */
void obj_track::write(int k,double time) {
    file_output(time);
    fflush(fp);
}

/** Writes the tracer positions to the output file.
 * \param[in] time the current simulation time. */
void obj_track::file_output(double time) {
    correct_tracers();
    *buf=time;
    double *bp=buf+1;
    for(double *tp=tm;tp<te;tp+=6) {*(bp++)=*tp; *(bp++)=tp[1];}
    fwrite(buf,sizeof(double),rec_size,fp);
    count=1;
}

/** Updates the tracer positions, integrating them with the improved Euler
 * method.
 * \param[in] i the step in the Euler method.
 * \param[in] f2d a reference to the parent fluid_2d class.
 * \param[in] dt the timestep to use. */
void obj_track::update_tracers(int i,fluid_2d &f2d,double dt) {
    if(i==1) {
        for(double *tp=tm;tp<te;tp+=6) f2d.tracer_update_1(tp);
    } else {
        for(double *tp=tm;tp<te;tp+=6) f2d.tracer_update_2(tp,dt);
        if(count==every) file_output(f2d.time);
        else count++;
    }
}

/** Corrects the tracer positions to exactly agree with the bicubic
 * interpolation of the reference map field. If the tracer correction fails,
 * then it is switched on from this point onward. */
void obj_track::correct_tracers() {
    if(!t_correction) return;
    for(double *tp=tm; tp<te; tp+=6) {
        if(!correct_tracer(tp[4],tp[5],*tp,tp[1])) {
          printf("# Tracer correction failed for (X,Y)=(%g,%g), (x,y)=(%g,%g)\n",
                 tp[4],tp[5],*tp,tp[1]);
          t_correction=false;
        }
    }
}

/** Corrects tracer position to be anchored to the reference map field.
 * \param[in] (xref,yref) the reference map position to anchor to.
 * \param[in,out] (xx,yy) the tracer position to be corrected.
 * \return True if the anchoring was successful, false otherwise. */
bool obj_track::correct_tracer(double xref,double yref,double &xx,double &yy) {
    const double tolnewt=small_number*small_number*100;
    const double bigstep=1;
    vec f,fx,fy;
    double det,delx,dely;

    int count=0,iter=0;
    while(count<4) {

        // Check for too many iterations
        if(iter>100) return false;

        // Compute function values and Jacobian
        f=bic->f_grad_f(xx,yy,fx,fy)-vec(xref,yref);

        // Check for convergence
        if(f.x*f.x+f.y*f.y<tolnewt) count++;

        // Bail out if the determinant is within the machine epsilon
        det=fx.x*fy.y-fy.x*fx.y;
        if(det<small_number&&det>-small_number) return false;

        // Compute update and bail out if there's a huge step
        delx=( fy.y*f.x-fy.x*f.y)/det;
        dely=(-fx.y*f.x+fx.x*f.y)/det;
        if(delx*delx+dely*dely>bigstep) return false;

        // Apply updates
        xx-=delx;yy-=dely;iter++;
    }
    return true;
}

/** Sets up the object representing a flag.
 * \param[in] (cx_,cy_) the translation to apply to the flag position.
 * \param[in] flag_l_ the length of the flag.
 * \param[in] flag_w_ the width of the flag.
 * \param[in] perturb_v_ the initial velocity perturbation to apply to the flag.
 * \param[in] ntrace the number of tracers to initialize along the flag's length.
 * \param[in] filename_ the name of the f*/
obj_flag::obj_flag(double cx_,double cy_,double flag_l_,double flag_w_,
                   double flag_a_,double perturb_v_,int ntrace,const char* filename_) :
    obj_track(ntrace, 100, filename_), cx(cx_), cy(cy_), flag_l(flag_l_),
    flag_w(flag_w_), flag_a(flag_a_), perturb_v(perturb_v_) {

        // Initialize tracers along length of the flag
        double step=ntrace==1?0.:flag_l/(ntrace-1),ex=flag_l,*tp=tm;
        for(;tp<te;tp+=6,ex-=step) {
          *tp=ex;tp[1]=0;
          transform(*tp,tp[1],tp[4],tp[5]);
        }

        // Check that the anchoring region is within the flag's blurred region.
        // Otherwise its influence will be cut off.
        if(flag_a>=0.5*flag_w) {
          fputs("flag_a should be smaller than 0.5*flag_w to avoid cutting off its influence\n",stderr);
          exit(1);
        }
    }

/** Computes the level set function for the flag, as a rectangle with
 * semicircular end caps.
 * \param[in] (X,Y) the reference map.
 * \return The level set function. */
double obj_flag::phi(double X,double Y) {
	return (X<0?sqrt(X*X+Y*Y)
               :(X>flag_l?sqrt((X-flag_l)*(X-flag_l)+Y*Y):fabs(Y)))-0.5*flag_w;
}

/** Sets up the flag object, initializing the tracer output file, and
 * calculating constants.
 * \param[in] f2d a reference to the parent fluid_2d class. */
void obj_flag::start_extra(fluid_2d &f2d) {
    open_track_file();
    setup_bicubic(f2d);

    // Compute the threshold for quickly determining if a point is inside the
    // anchor region
    flag_thresh=flag_a+op->eps;
    flag_thresh*=flag_thresh;

    // Set the anchoring stiffness so that the natural frequency is ten times
    // the numerical stability limit
    K_stiff=0.01/(f2d.dt_reg*f2d.dt_reg);
    printf("# Anchoring K_stiff       : %g (1/T^2)\n",K_stiff);
}

/** Compute the anchor force at position (x, y) with radius r
 * \param[in] (x,y) the center of the anchor.
 * \param[in] (X,Y) the reference map of the anchor.
 * \param[in] phiv the level set value at this position.
 * \param[out] (accx,accy) the anchor force. */
void obj_flag::accel(double x,double y,double X,double Y,
                     double phiv,double &accx,double &accy) {
    double delx=x-cx,dely=y-cy,rsq=delx*delx+dely*dely;
    if(rsq<flag_thresh) {
        double K=K_stiff*op->trans_func_in(sqrt(rsq)-flag_a);
        accx-=K*(x-X);
        accy-=K*(y-Y);
    }
}

/** Sets the flag velocity to have small perturbation in the tail.
 * \param[in] (X,Y) the reference map at which to compute the velocity.
 * \param[out] (u,v) the computed velocity. */
void obj_flag::velocity(double X,double Y,double &u,double &v) {
    X-=flag_l;
    double rsq=X*X+Y*Y;
    u=0.;
    v=rsq<flag_thresh?perturb_v*op->trans_func_in(sqrt(rsq)-flag_a):0.;
}

/** Initializes the circle object that can track its position and rotation.
 * \param[in] (cx_,cy_) the center of the circle.
 * \param[in] cr_ the radius of the circle.
 * \param[in] (uu_,vv_) the initial velocity of the circle.
 * \param[in] omega_ the initial angular velocity of the circle.
 * \param[in] filename_ the name of the file to write tracking information to. */
obj_circle_track::obj_circle_track(double cx_,double cy_,double cr_,double uu_,double vv_,
                                   double omega_,const char* filename_)
    : obj_track(2, 32, 4, filename_), cx(cx_), cy(cy_), cr(cr_), uu(uu_),
    vv(vv_), omega(omega_), prev_q(0), wind(0) {

    // Initialize the two tracers
    *tm=cx-0.5*cr;tm[1]=cy;transform(*tm,tm[1],tm[4],tm[5]);
    tm[6]=cx+0.5*cr;tm[7]=cy;transform(tm[6],tm[7],tm[10],tm[11]);
}

/** Computes the level set function as a function of position.
 * \param[in] (x,y) the position to consider.
 * \return The level set function. */
double obj_circle_track::phi(double X,double Y) {
    return sqrt(X*X+Y*Y)-cr;
}

/** Calculates the initial velocity, as a combination of a translational motion
 * and an angular spin.
 * \param[in] (X,Y) the reference map at the position.
 * \param[out] (u,v) the velocity. */
void obj_circle_track::velocity(double X,double Y,double &u,double &v) {
    u=uu-omega*Y;
    v=vv+omega*X;
}

/** Saves the tracer-related output to file.
 * \param[in] time the current simulation time. */
void obj_circle_track::file_output(double time) {
    double ang=atan2(tm[5]-tm[1],tm[4]-*tm);
    correct_tracers();

    // Check for a change in winding number
    int q=fabs(ang)>0.5*M_PI?(ang>0?1:-1):0;
    if(q==1&&prev_q==-1) wind--;
    else if(q==-1&&prev_q==1) wind++;
    prev_q=q;

    // Assemble the output buffer and save it to file
    *buf=time;
    buf[1]=*tm;
    buf[2]=tm[1];
    buf[3]=ang+2*M_PI*wind;
    fwrite(buf,sizeof(double),rec_size,fp);
    count=1;
}

/** Computes the level set function as a function of position.
 * \param[in] (x,y) the position to consider.
 * \return The level set function. */
double obj_lamprey::phi(double X,double Y) {
    if(X>0) return dis*Y*Y/(w1*w1)+X-dis;
    double wl=w1-(w2-w1)*X/la;
    return std::max(-la-X,Y*Y/(wl*wl)-1);
}

/** Computes the overall amplitude of swimming, and the shift of the
 * contraction wave down the lamprey body.
 * \param[in] the current simulation time. */
void obj_lamprey::pre_stress_setup(double time) {
    sh=time/T;
    double o=2*sh;
    con=o<1?A*o*o*(3-2*o):A;
    sh*=la;
}

/** Apply the actuation to create the lamprey's swimming motion.
 * \param[in] (X,Y) the reference map at this point.
 * \param[in,out] F the deformation gradient tensor to modify. */
void obj_lamprey::actuate(double X,double Y,mat &F) {
    double dx=X+(Y>0?1.5*la:la)+sh,
           wl=w1-(w2-w1)*X/la,fac=fabs(Y)/wl,
           st=exp(con*(fac<1?fac*(2-fac):1)
              *op->trans_func(-std::min(-X,la*0.15-fabs(dx-la*(int(dx/la)+0.5)))));
    F.a*=st;F.b*=1./st;
    F.c*=st;F.d*=1./st;
}

/** Computes the level set function as a function of position.
 * \param[in] (X,Y) the position to consider.
 * \return The level set function. */
double obj_hoop::phi(double X,double Y) {
    return fabs(sqrt(X*X+Y*Y)-cc)-cr;
}
