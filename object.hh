#ifndef OBJECT_HH
#define OBJECT_HH

#include <cstdio>
#define _USE_MATH_DEFINES
#include <cmath>

#include "common.hh"
#include "mat.hh"
#include "bi_interp.hh"

class obj_field;
class fluid_2d;

/** \brief Base class to define an solid object within the simulation. */
class object {
    public:
        virtual ~object() {}
        virtual double phi(double X,double Y) = 0;
        /** Calculates the initial reference map at a position. This default
         * function uses the identity mapping, but it can be overridden to
         * apply different transforms.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            X=x;Y=y;
        }
        /** Calculates the initial velocity at a position. This default
         * function sets the velocity to be zero, but it can be overridden.
         * \param[in] (X,Y) the reference map at the position.
         * \param[out] (u,v) the velocity. */
        virtual void velocity(double X,double Y,double &u,double &v) {
            u=v=0.;
        }
        /** Modifies the deformation gradient tensor to apply actuation. This
         * default function applies no modification, but it can be overridden.
         * \param[in] (X,Y) the reference map at the position.
         * \param[in,out] F the deformation gradient tensor to modify. */
        virtual void actuate(double X,double Y,mat &F) {}
        /** Links the parent obj_field class to this object, and performs
         * any additional initialization.
         * \param[in] f2d a reference to the parent fluid_2d class.
         * \param[in] op_ a pointer to the obj_field class representing this
         *                object. */
        inline void start(fluid_2d &f2d,obj_field *op_) {
            op=op_;start_extra(f2d);
        }
        /** Performs additional initialization. This default function does
         * nothing, but it can be overridden.
         * \param[in] f2d a reference to the parent fluid_2d class. */
        virtual void start_extra(fluid_2d &f2d) {}
        /** Some objects contain special tracers. This function is called to
         * update their positions. This default function does nothing, but it
         * can be overridden.
         * \param[in] i the RK step number in the improved Euler method to
         *              perform.
         * \param[in] f2d a reference to the parent fluid_2d class.
         * \param[in] dt the timestep. */
        virtual void update_tracers(int i,fluid_2d &f2d,double dt) {}
        /** Performs any object-related output when a simulation frame is being
         * stored. This default function does nothing, but it can be
         * overridden.
         * \param[in] k the frame number being stored.
         * \param[in] time the current simulation time. */
        virtual void write(int k,double time) {}
        /** Applies extra accelerations to the object, such as due to an
         * anchor. This default function does nothing, but it can be
         * overridden.
         * \param[in] (x,y) the current position.
         * \param[in] (X,Y) the reference map at this position.
         * \param[in] phiv the level set value at this position.
         * \param[in,out] (accx,accy) the acceleration vector to add to. */
        virtual void accel(double x,double y,double X,double Y,
                           double phiv,double &accx,double &accy) {}
        virtual void pre_stress_setup(double time) {}
    protected:
        /** A pointer to the corresponding obj_field class. */
        obj_field *op;
};

/** \brief A special object type representing a solid filling the entire
 * domain. */
struct obj_full : public object {
    virtual ~obj_full() {}
    virtual double phi(double X,double Y) {return -1;}
};

/** \brief A base class for many objects, allowing for them to be translated
 * to different positions and initialized with velocity and spin. */
struct obj_standard : public object {
    /** The x position of the object center. */
    const double cx;
    /** The y position of the object center. */
    const double cy;
    /** A length scale for the object size. */
    const double cr;
    /** The initial horizontal velocity of the object. */
    const double uu;
    /** The initial vertical velocity of the object. */
    const double vv;
    /** The initial angular velocity of the object. */
    const double omega;
    /** Initializes the class constants, setting the velocity parameters to
     * zero.
     * \param[in] (cx_,cy_) the position of the object center.
     * \param[in] cr_ the object length scale. */
    obj_standard(double cx_,double cy_,double cr_) :
        cx(cx_), cy(cy_), cr(cr_), uu(0.), vv(0.), omega(0.) {}
    /** Initalizes the class constants.
     * \param[in] (cx_,cy_) the position of the object center.
     * \param[in] cr_ the object length scale.
     * \param[in] (uu_,vv_) the initial object velocity.
     * \param[in] (omega_) the initial object angular velocity. */
    obj_standard(double cx_,double cy_,double cr_,double uu_,double vv_,double omega_) :
        cx(cx_), cy(cy_), cr(cr_), uu(uu_), vv(vv_), omega(omega_) {}
    /** Calculates the initial reference map at a position, applying a
     * translation by the (cx,cy) vector.
     * \param[in] (x,y) the position to consider.
     * \param[out] (X,Y) the coordinates of the reference map. */
    virtual void transform(double x,double y,double &X,double &Y) {
        X=x-cx;Y=y-cy;
    }
    virtual void velocity(double X,double Y,double &u,double &v);
};

/** \brief A class to define a circle that experiences an accerelation in the
 * negative y direction (e.g. due to gravity). */
struct obj_circle : public obj_standard {
    /** The acceleration to apply in the negative y direction. */
    const double grav;
    obj_circle(double cx_,double cy_,double cr_,double grav_=0) :
        obj_standard(cx_,cy_,cr_), grav(grav_) {}
    obj_circle(double cx_,double cy_,double cr_,double uu_,double vv_,double omega_,double grav_=0) :
        obj_standard(cx_,cy_,cr_,uu_,vv_,omega_), grav(grav_) {}
    virtual double phi(double X,double Y);
    virtual void accel(double x,double y,double X,double Y,
                       double phiv,double &accx,double &accy);
};

/** \brief A class to define a square that experiences an accerelation in the
 * negative y direction (e.g. due to gravity). */
struct obj_square : public obj_standard {
    /** The acceleration to apply in the negative y direction. */
    const double grav;
    obj_square(double cx_,double cy_,double cr_,double grav_=0) :
        obj_standard(cx_,cy_,cr_), grav(grav_) {}
    obj_square(double cx_,double cy_,double cr_,double uu_,double vv_,double omega_,double grav_=0) :
        obj_standard(cx_,cy_,cr_,uu_,vv_,omega_), grav(grav_) {}
    virtual double phi(double X,double Y);
    virtual void accel(double x,double y,double X,double Y,
                       double phiv,double &accx,double &accy);
};

/** \brief A class to define an equilateral triangle. */
class obj_triangle : public obj_standard {
    public:
        const double grav;
        obj_triangle(double cx_,double cy_,double cr_,double theta,double grav_);
        virtual double phi(double X,double Y);
        virtual void accel(double x,double y,double X,double Y,
                           double phiv,double &accx,double &accy);
    private:
        /** The normal vectors for the sides of triangle. */
        double nor[6];
        /** Computes the scalar product between the reference map and one of
         * the normal vectors.
         * \param[in] (X,Y) the reference map.
         * \param[in] k the normal vector to use.
         * \return The scalar product. */
        inline double sp(double X,double Y,int k) {
            return X*nor[2*k]+Y*nor[2*k+1];
        }
};

/** \brief A class to define a seven-pointed star. */
class obj_seven_point : public object {
    public:
        /** The x position of the object center. */
        const double cx;
        /** The y position of the object center. */
        const double cy;
        /** The distance from the center to the star tips. */
        const double cr;
        /** The radius over which to apply the pivoting force. */
        const double piv_r;
        /** Half of the maximum theta value to spin the star to. */
        const double th_max;
        /** The angular velocity of the applied spinning term. */
        const double omega;
        /** The stiffness constant for the pivot. */
        double K_stiff;
        obj_seven_point(double cx_,double cy_,double cr_,double piv_r_,double th_max_,
                double omega_) :
            cx(cx_), cy(cy_), cr(cr_), piv_r(piv_r_), th_max(th_max_), omega(omega_) {}
        /** Calculates the initial reference map at a position, applying a
         * translation by the (cx,cy) vector.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            X=x-cx;Y=y-cy;
        }
        virtual void start_extra(fluid_2d &f2d);
        virtual void accel(double x,double y,double X,double Y,
                           double phiv,double &accx,double &accy);
        virtual double phi(double X,double Y);
        virtual void pre_stress_setup(double time);
    private:
        /** A threshold used for rapid computations of the pivot region. */
        double piv_thresh;
        /** The cosine of the current angle of the spinner. */
        double cth;
        /** The sine of the current angle of the spinner. */
        double sth;
};

/** \brief A class to define a rounded rod that experiences a gravitational
* pull. */
struct obj_rounded_rod : public object {
    /** The x position of the rod center. */
    const double cx;
    /** The y position of the rod center. */
    const double cy;
    /** The radius of the rod end caps. */
    const double cr;
    /** The half-length of the rod (not counting the end caps). */
    const double cl;
    /** The gravitational acceleration to apply to the rod. */
    const double grav;
    obj_rounded_rod(double cx_,double cy_,double cr_,double cl_,double grav_) :
        cx(cx_), cy(cy_), cr(cr_), cl(cl_), grav(grav_) {}
    /** Calculates the initial reference map at a position, applying a
     * translation by the (cx,cy) vector.
     * \param[in] (x,y) the position to consider.
     * \param[out] (X,Y) the coordinates of the reference map. */
    virtual void transform(double x,double y,double &X,double &Y) {
        X=x-cx;Y=y-cy;
    }
    virtual void accel(double x,double y,double X,double Y,
                       double phiv,double &accx,double &accy);
    virtual double phi(double X,double Y);
};

/** \brief A class to define a filament that can swim with a jellyfish-like
 * motion. */
struct obj_flapper : public object {
    public:
        /** The x position of the flapper center. */
        const double cx;
        /** The y position of the flapper center. */
        const double cy;
        /** The radius of the flapper end caps. */
        const double cr;
        /** The half-length of the flapper (not counting the end caps). */
        const double cl;
        /** The radius of the end caps of the actuated region. */
        const double ar;
        /** The half-length of the actuated region. */
        const double al;
        /** The maximum flapping amplitude. */
        const double amp;
        /** The angular frequency of the flapping motion. */
        const double omega;
        obj_flapper(double cx_,double cy_,double cr_,double cl_,double ar_,
                    double al_,double theta,double m_stretch,double T) :
            cx(cx_), cy(cy_), cr(cr_), cl(cl_), ar(ar_), al(al_),
            amp(log(m_stretch)/ar), omega(2*M_PI/T), cth(cos(theta)),
            sth(sin(theta)) {}
        /** Calculates the initial reference map at a position, applying a
         * translation by the (cx,cy) vector, and applying a rotation.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            x-=cx;y-=cy;
            X=x*cth+y*sth;
            Y=-x*sth+y*cth;
        }
        virtual double phi(double X,double Y);
        virtual void pre_stress_setup(double time);
        virtual void actuate(double X,double Y,mat &F);
    private:
        /** Cosine of the flapper rotation angle. */
        const double cth;
        /** Sine of the flapper rotation angle. */
        const double sth;
        /** The current contraction level of the flapper. */
        double con;
        /** A pointer to the corresponding obj_field class. */
        obj_field *op;
};

/** \brief A class that can track one or more special tracers that move with
 * the object. */
struct obj_track : public object {
    public:
        /** The filename to write the tracer information to. */
        const char* filename;
        /** The number of tracers. */
        const int ntrace;
        /** The frequency at which to store the tracer information. */
        const int every;
        /** The timestep counter, used to determine when to save output snapshots.
         */
        int count;
        /** Whether to periodically correct the tracers so that they exactly align
         * with the reference map. */
        bool t_correction;
        /** Initializes the obj_track class, setting up constants and
         * allocating memory for the tracers.
         * \param[in] ntrace_ the number of tracers to use.
         * \param[in] every_ the number of timesteps after which to output
         *                   data.
         * \param[in] filename_ the name of the output file to write to.
         * \param[in] t_correction_ whether to periodically correct the tracer
         *                          positions to exactly align with the
         *                          reference map. */
        obj_track(int ntrace_,int every_,const char* filename_,bool t_correction_=true) :
            filename(filename_), ntrace(ntrace_), every(every_), count(1),
            t_correction(t_correction_), rec_size(1+(ntrace<<1)),
            tm(new double[6*ntrace]), te(tm+6*ntrace), fp(NULL),
            buf(new double[rec_size]), bic(NULL) {}
        /** Initializes the obj_track class, setting up constants and
         * allocating memory for the tracers. This special constructor adjusts
         * the size of the output data, in cases where special tracers analysis
         * is performed.
         * \param[in] ntrace_ the number of tracers to use.
         * \param[in] every_ the number of timesteps after which to output
         *                   data.
         * \param[in] rec_size_ the size of each output record (in
         *                      double-precision floating point numbers).
         * \param[in] filename_ the name of the output file to write to.
         * \param[in] t_correction_ whether to periodically correct the tracer
         *                          positions to exactly align with the
         *                          reference map. */
        obj_track(int ntrace_,int every_,int rec_size_,const char* filename_,bool t_correction_=true) :
            filename(filename_), ntrace(ntrace_), every(every_), count(1),
            t_correction(t_correction_), rec_size(rec_size_),
            tm(new double[6*ntrace]), te(tm+6*ntrace), fp(NULL),
            buf(new double[rec_size]), bic(NULL) {}
        /** The class constructor closes the output file and frees the
         * dynamically allocated memory. */
        virtual ~obj_track() {
            if(bic!=NULL) delete bic;
            if(fp!=NULL) fclose(fp);
            delete [] buf;
            delete [] tm;
        }
        virtual void start_extra(fluid_2d &f2d);
        virtual void update_tracers(int i,fluid_2d& f2d,double dt);
        virtual void write(int k,double time);
        virtual void file_output(double time);
        void open_track_file();
        void correct_tracers();
        bool correct_tracer(double xref,double yref,double &xx,double &yy);
        void setup_bicubic(fluid_2d &f2d);
    protected:
        /** The amount of data store per snapshot in the output file. */
        const size_t rec_size;
        /** Memory for storing the tracer positions. */
        double* const tm;
        /** A pointer to the end of tracer memory. */
        double* const te;
        /** The file handle for the output. */
        FILE *fp;
        /** A buffer for constructing the output snapshots. */
        double *buf;
        /** A bicubic interpolation interpolations of the velocity field. */
        bicubic_interp* bic;
};

/** \brief Class to represent a flexible piston that is anchored at one end. */
struct obj_piston : public obj_track {
    public:
        /** The half-width of the piston. */
        const double lx;
        /** The half-height of the piston. */
        const double ly;
        /** The x position of the piston anchor. */
        const double cx;
        /** The starting y position of the piston. */
        const double y_start;
        /** The total vertical displacement of the piston anchor. */
        const double len;
        /** The reciprocal time for the anchor to move the complete
         * displacement. */
        const double idur;
        /** The radius of the anchored region. */
        const double anch_r;
        /** The current y position of the piston anchor. */
        double cy;
        /** The size of the anchoring acceleration. */
        double K_stiff;
        obj_piston(double lx_,double ly_,double cx_,double y_start_,double y_end_,
                   double duration_,double anch_r_,const char* filename_) :
            obj_track(1,100,filename_), lx(lx_), ly(ly_), cx(cx_),
            y_start(y_start_), len(y_end_-y_start), idur(1./duration_),
            anch_r(anch_r_) {

            // Set up the tracer position at the mirror image of the anchor
            // center
            *tm=-cx_;
            tm[1]=y_start;
            transform(*tm,tm[1],tm[4],tm[5]);
        }
        /** Calculates the initial reference map at a position, applying a
         * translation by the (cx,cy) vector, and applying a rotation.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            X=x;Y=y-y_start;
        }
        virtual void start_extra(fluid_2d &f2d);
        virtual void accel(double x,double y,double X,double Y,
                           double phiv,double &accx,double &accy);
        virtual void pre_stress_setup(double time);
        virtual double phi(double X,double Y);
    private:
        /** A threshold for quickly determining if a grid point is within
         * the anchor. */
        double anch_thresh;
};

/** \brief Class to represent a simple rotating spinner, primarily used for
 * resolution and performance tests. */
struct obj_simple_spin : public obj_track {
    public:
        /** The scale factor that sets the overall size of the spinner. */
        const double scf,c0,c1,c2;
        /** The coefficient of the third Fourier mode in the R(theta) function
         * describing the spinner's boundary. */
        const double a1;
        /** The coefficient of the sixth Fourier mode in the R(theta) function
         * describing the spinner's boundary. */
        const double a2;
        /** The radius of the pivot region. */
        const double piv_r;
        /** Half of the maximum theta value that the spinner attains. */
        const double th_max;
        /** The angular velocity prefactor in the spinner's rotation function. */
        const double omega;
        /** The stiffness constant for the pivot. */
        double K_stiff;
        /** The previous quadrant that the tracer was in. */
        int prev_q;
        /** The winding number of the tracer. */
        int wind;
        obj_simple_spin(double r_scale_,double a1_,double a2_,double piv_r_,
                        double th_max_,double omega_,const char* filename_) :
            obj_track(1,10,filename_), scf(r_scale_/(1+a1_+a2_)), c0(scf*(1-a2_)),
            c1(scf*a1_), c2(2*scf*a2_), a1(a1_),
            a2(a2_), piv_r(piv_r_), th_max(th_max_), omega(omega_), prev_q(0),
            wind(0) {
            *tm=0.9*r_scale_;
            tm[1]=0;
            transform(*tm,tm[1],tm[4],tm[5]);
        }
        virtual void start_extra(fluid_2d &f2d);
        virtual void accel(double x,double y,double X,double Y,
                           double phiv,double &accx,double &accy);
        virtual double phi(double X,double Y);
        virtual void file_output(double time);
        virtual void pre_stress_setup(double time);
    private:
        /** A threshold used for rapid computations of the pivot region. */
        double piv_thresh;
        /** The cosine of the current angle of the spinner. */
        double cth;
        /** The sine of the current angle of the spinner. */
        double sth;
};

/** Class to represent an filament anchored at one end, which can be used as a
 * simplified model of a flag. */
struct obj_flag : public obj_track {
    public:
        /** The x translation to apply to the flag. */
        const double cx;
        /** The y translation to apply to the flag. */
        const double cy;
        /** The length of the flag. */
        const double flag_l;
        /** The width of the flag. */
        const double flag_w;
        /** The radius of the flag anchoring region. */
        const double flag_a;
        /** The initial velocity perturbation to apply to the flag. */
        const double perturb_v;
        /** The anchoring acceleration constant. */
        double K_stiff;
        obj_flag(double cx_,double cy_,double flag_l_,double flag_w_,double flag_a_,
                 double perturb_v_,int ntrace,const char* filename_);
        virtual void start_extra(fluid_2d &f2d);
        /** Calculates the initial reference map at a position, applying a
         * translation by the (cx,cy) vector.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            X=x-cx;Y=y-cy;
        }
        virtual void velocity(double X,double Y,double &u,double &v);
        virtual double phi(double X,double Y);
        virtual void accel(double x,double y,double X,double Y,
                           double phiv,double &accx,double &accy);
    private:
        /** The threshold to use to quickly determine if a grid point is
         * anchored or not. */
        double flag_thresh;
};

/** \brief Class to represent a circle that can also track its position and
 * rotation. */
struct obj_circle_track : public obj_track {
    /** The x position of the circle center. */
    const double cx;
    /** The y position of the circle center. */
    const double cy;
    /** The radius of the circle. */
    const double cr;
    /** The initial horizontal velocity of the circle. */
    const double uu;
    /** The initial vertical velocity of the circle. */
    const double vv;
    /** The initial angular velocity of the circle. */
    const double omega;
    /** The previous quadrant that the tracer was in. */
    int prev_q;
    /** The winding number of the tracer. */
    int wind;
    obj_circle_track(double cx_,double cy_,double cr_,double uu_,
            double vv_,double omega_, const char* filename_);
    virtual double phi(double X,double Y);
    /** Calculates the initial reference map at a position, applying a
     * translation by the (cx,cy) vector.
     * \param[in] (x,y) the position to consider.
     * \param[out] (X,Y) the coordinates of the reference map. */
    virtual void transform(double x,double y,double &X,double &Y) {
        X=x-cx;Y=y-cy;
    }
    virtual void velocity(double X,double Y,double &u,double &v);
    virtual void file_output(double time);
};

/** \brief Class to represent a swimming lamprey, using the geometry and
 * actuation of Tytell et al. (PNAS, 2010). */
struct obj_lamprey : public object {
        /** The x position of the lamprey center. */
        const double cx;
        /** The y position of the lamprey center. */
        const double cy;
        /** The length of the head of the lamprey. */
        const double dis;
        /** The length of the body of the lamprey. */
        const double la;
        /** The half-width of the lamprey at the body--head connection. */
        const double w1;
        /** The half-width of the lamprey at the tail. */
        const double w2;
        /** The period of the lamprey swimming motion. */
        const double T;
        /** The swimming amplitude. */
        const double A;
        /** Initializes the swimming lamprey.
         * \param[in] (cx_,cy_) the center of the lamprey.
         * \param[in] theta the lamprey rotation angle.
         * \param[in] le_ the total length of the lamprey.
         * \param[in] dis_ the length of the head of the lamprey.
         * \param[in] w1_ the half-width of the lamprey at the body--head
         *                connection.
         * \param[in] w2_ the half-width of the lamprey at the tail.
         * \param[in] T_ the period of the lamprey swimming motion.
         * \param[in] mlam_ the maximum contraction to apply during swimming.
         */
        obj_lamprey(double cx_,double cy_,double theta,double le_,
                    double dis_,double w1_,double w2_,double T_,double mlam) :
            cx(cx_), cy(cy_), dis(dis_), la(le_-dis), w1(w1_), w2(w2_),
            T(T_), A(log(mlam)), cth(cos(theta)), sth(sin(theta)) {}
        /** Calculates the initial reference map at a position, applying a
         * translation by the (cx,cy) vector, and a rotation.
         * \param[in] (x,y) the position to consider.
         * \param[out] (X,Y) the coordinates of the reference map. */
        virtual void transform(double x,double y,double &X,double &Y) {
            x-=cx;y-=cy;
            X=x*cth+y*sth;
            Y=-x*sth+y*cth;
        }
        virtual double phi(double X,double Y);
        virtual void pre_stress_setup(double time);
        virtual void actuate(double X,double Y,mat &F);
    private:
        /** Cosine of the lamprey rotation angle. */
        const double cth;
        /** Sine of the lamprey rotation angle. */
        const double sth;
        /** The current contraction level of the lamprey. */
        double con;
        /** The shift of the contraction wave in the lamprey. */
        double sh;
};

/** \brief Class to represent a circular hoop. */
struct obj_hoop : public object {
    /** The average of the inner and outer radii of the hoop. */
    const double cc;
    /** The thickness of the hoop. */
    const double cr;
    /** Gravity to apply to the hoop in the negative y direction. */
    const double grav;
    /** Initializes the circular hoop object.
     * \param[in] (r_in,r_out) the inner and outer radii of the hoop.
     * \param[in] grav the gravitational acceleration to apply to the hoop. */
    obj_hoop(double r_in,double r_out,double grav_=0) :
        cc(0.5*(r_in+r_out)), cr((r_out-r_in)*0.5), grav(grav_) {}
    virtual double phi(double X,double Y);
};

#endif
