#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

#include "common.hh"
#include "object.hh"
#include "obj_field.hh"
#include "fluid_2d.hh"

/** Checks to see if two strings are equal.
 * \param[in] (p1,p2) the two strings.
 * \return True if the are equal, false if they are not. */
inline bool se(const char* p1,const char *p2) {
    return strcmp(p1,p2)==0;
}

/** Checks to ensure that the number of command-line arguments matches the
 * number that is expected.
 * \param[in] argc the number of command-line arguments.
 * \param[in] req the number of extra command-line arguments expected. */
void check_cmd(int argc,int req) {
    if(argc!=req+2) {
        fprintf(stderr,"Selected option requires %d argument%s (%d supplied)\n",req,req==1?"":"s",argc-2);
        exit(1);
    }
}

/** Checks that the grid resolution provided is in a reasonable range.
 * \param[in] m the grid resolution. */
void check_grid(int m) {
    if(m<0) {
        fputs("The grid resolution must be positive\n",stderr);
        exit(1);
    }
    if(m<20||m>65536) {
        fprintf(stderr,"Supplied grid resolution of %d is too %s\n",m,m<20?"small":"large");
        exit(1);
    }
}

/** Prints an error message with the command-line arguments for a simulation.
 * \param[in] na the name of the simulation.
 * \param[in] ar the command-line arguments. */
void cmd_args(const char* na,const char *ar) {
    fprintf(stderr,"./ftest %s %s\n",na,ar);
    exit(1);
}

int main(int argc,char **argv) {

    // Check for at least one command-line argument. If not, provide a list
    if(argc<2) {
        fputs("Syntax: ./ftest <sim_type> {extra options}\n\n"
              "Simulation types:\n"
              "simple-spin     simple three-pronged rotor\n"
              "circ-move       circle moving through fluid\n"
              "square          rotating square\n"
              "seven-star      rotating seven-pointed star\n"
              "piston          flexible paddle moving through cavity\n"
              "fluid           fluid only\n"
              "full           solid only\n"
              "flapper         single flapping swimmer\n"
              "flappers        multiple flapping swimmers\n"
              "conc-hoops      concentric hoops\n"
              "multi-drop      multiple objects sedimenting\n"
              "flag            flapping flag\n",stderr);
        return 1;
    }

    // Some default settings
    int m=200,n=200;             // Simulation resolution
    double T=4;                  // Simulation duration
    int num_frames=200;          // Number of output frames
    double tw=2.5;               // Transition width (in grid spacings)
    double pinw=2.5;             // Pinning width (in grid spacings)
    int ntrace=512;              // Number of tracers
    bool set_velocity=true;      // Whether solids initialize their own velocity, or
                                 // inherit the fluid's velocity
    bool solid_prd_bc=false;     // Whether to apply periodic boundary conditions to
                                 // the reference map fields

    // File output flags
    // 1 - u                128 - solid Y
    // 2 - v                256 - solid phi
    // 4 - speed            512 - only subregions for solids
    // 8 - pressure         1024 - add ghost regions
    // 16 - vorticity       2048 - difference fields
    // 32 - density         4096 - divergence
    // 64 - solid X         8192 - parsable timing information
    unsigned int fflags=8|16|64|128|256|512;

    // Default fluid/solid parameters
    double rhof=1;            // Fluid density
    double rhos=1;            // Solid density
    double visc=1e-3;         // Fluid viscosity
    double G=24;              // Shear modulus

    // Specify padding factors for extra viscosity and timestep
    double ev_mult=0.8;       // The multiplier to apply to the extra viscosity
    double ev_trans_mult=1;   // The multiplier to apply to the additional
                              // extra viscosity term in the transition region
    double dt_pad=0.4;        // The padding to apply to the timestep restriction
                              // due to physical terms
    double dt_ev_pad=0.8;     // The padding to apply to the timestep restriction
                              // due to the extra viscosity

    // Set up different simulation examples
    sim_type *sim_t=NULL;
    char *fn=argv[1],*fa=NULL;
    std::vector<object*> obj;

    // Parse first argument to determine simulation type
    if(se(argv[1],"simple-spin")) {

        // Check for the right number of command-line arguments
        if(argc<3||argc>4)
            cmd_args("simple-spin","<grid resolution> [<transition width>]");

        // Open output filename
        fn=fa=new char[128];
        if(argc==3) {
            sprintf(fa,"sspin_%s",argv[2]);
            num_frames=240;
        } else {
            sprintf(fa,"sspin_%s_tw%s",argv[2],argv[3]);
            tw=atof(argv[3]);

            // If the transition width is provided, then switch to parseable
            // timing information, and include divergence fields
            fflags|=8192|4096;
        }

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        G=48;visc=1e-2;rhos=3;
        T=4*M_PI;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_simple_spin(0.8,0.5,0.125,0.2,M_PI,1,fa));
        sim_t=new sim_type(-1,1,-1,1,false,false);
    } else if(se(argv[1],"circ-move")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circ-move","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=4;visc=1e-3;
        T=3.;num_frames=300;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_track(-0.5,-0.5,0.25,2,4,24,argv[1]));
        sim_t=new sim_type(-1,1,-1,1,false,false);
    } else if(se(argv[1],"square")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("square","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=3;

        // Create object
        obj.push_back(new obj_square(-0.02,0,0.5,0.01,0,2));
    } else if(se(argv[1],"seven-star")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("seven-star","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=3;
        T=4*M_PI;num_frames=600;

        // Create object
        obj.push_back(new obj_seven_point(0,0,0.62,0.16,M_PI,1));
    } else if(se(argv[1],"piston")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("piston","<shear modulus>");
        check_cmd(argc,1);

        // Assemble output filename
        fn=fa=new char[128];
        sprintf(fa,"piston_G%s",argv[2]);

        // Read the shear modulus and check it is positive
        G=atof(argv[2]);
        if(G<=0) {
            fputs("The shear modulus must be positive\n",stderr);
            return 1;
        }

        // Set grid resolution and various simulation constants
        m=320;n=800;
        T=20;num_frames=1000;
        rhos=2.;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_piston(0.8,0.2,0.6,0.4,4.6,10,0.15,fa));
        sim_t=new sim_type(-1,1,0,5,false,false);
    } else if(se(argv[1],"fluid")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("fluid","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        visc=1e-3;
        T=20;num_frames=200;
        ntrace=256;

        // Set initial velocity to be the "pulses" field
        sim_t=new sim_velocity_pulses;

    } else if(se(argv[1],"full")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("full","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=1.;G=3.;visc=1e-3;
        T=1;num_frames=200;

        // Switch off tracers, since there is no fluid to initialize them in.
        // Set several special flags appropriate for solid-only simulation.
        ntrace=0;
        set_velocity=false;
        solid_prd_bc=true;

        // Create full solid, at set initial velocity to be the "pulses" field
        obj.push_back(new obj_full());
        sim_t=new sim_velocity_pulses;
    } else if(se(argv[1],"flapper")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("flapper","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=4;G=10;visc=5e-4;
        T=30;num_frames=600;

        // Creat the object, and set the large domain with periodic boundary
        // conditions
        obj.push_back(new obj_flapper(0,-0.8,0.026,0.25,0.021,0.14,0,2.2,8));
        sim_t=new sim_type(-1.5,1.5,-1.5,1.5,false,false);
    } else if(se(argv[1],"flappers")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("flappers","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=4;G=10;visc=5e-4;
        T=30;num_frames=1200;

        // Insert non-overlapping flappers
        const int cnum=28;
        double x,y,wx=0.25+0.026+0.01,wy=0.026+0.01,dx,dy;
        int i=0,j;
        bool good;
        while(i<cnum) {

            // Choose random position for flapper and check it is in range
            x=rnd(-1,1);
            y=rnd(-1,1);
            if(x+wx>0.98 || x-wx<-0.98 || y+wy>0.98 || y-wy<-0.98) continue;
            good=true;

            // Loop over previously added flappers and check for overlaps
            for(j=0;j<i;j++) {
                obj_flapper *oc=reinterpret_cast<obj_flapper*>(obj[j]);
                dx=x-oc->cx;dy=y-oc->cy;
                if(fabs(dx)<2*wx&&fabs(dy)<2*wy) {good=false;break;}
            }

            // If this flapper doesn't overlap with a previous one, then add it
            if(good) {
                i++;
                obj.push_back(new obj_flapper(x,y,0.026,0.25,0.021,0.14,0,2.2,rnd(7,9)));
            }
        }

        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-1,1,-1,1,false,false);
    } else if(se(argv[1],"multi-drop")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("multi-drop","<grid resolution> <type>\n\n"
                             "Type is 0 for circles and 1 for squares");
        bool square=atoi(argv[3])==1;

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=3;G=2;visc=1e-3;
        T=25;num_frames=1200;

        // Insert non-overlapping shapes
        const int cnum=42;
        double x,y,r,dx,dy,dr;
        int i=0,j;
        bool good;
        while(i<cnum) {

            // Choose a random position and radius for the shape and check it
            // is in range
            x=rnd(-1,1);
            y=rnd(-1,1);
            r=rnd(0.1,0.25);
            if(x+r>1||x-r<-1||y+r>1||y-r<-1) continue;
            good=true;

            // Loop over previously added shapes and check for overlaps
            for(j=0;j<i;j++) {
                obj_circle *oc=reinterpret_cast<obj_circle*>(obj[j]);
                dx=x-oc->cx;dy=y-oc->cy;dr=oc->cr+r;
                if(fabs(dx)<dr&&fabs(dy)<dr) {good=false;break;}
            }

            // If the shape doesn't overlap with a previous one, then add it
            if(good) {
                i++;
                square?obj.push_back(new obj_square(x,y,r-0.05,0,0,rnd(-5,5)))
                      :obj.push_back(new obj_circle(x,y,r-0.05,0,0,rnd(-5,5)));
            }
        }

        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-1,1,-1,1,false,false);
    } else if(se(argv[1],"flag")) {

        // Three arguments (matching Connell & Yue, JFM 2007), plus a fourth
        // for our additional degree of freedom
        // - Mass ratio (mu)
        // - Dimensionless bending rigidity (k_B)
        // - Reynolds number (Re)
        // - Aspect ratio (ar)
        // - [Perturbation velocity to apply to tip [optional]]
        if(argc<6||argc>7) {
            fprintf(stderr,"Selected option requires 4 or 5 arguments (%d supplied)\n",argc-2);
            exit(1);
        }
        double mu=atof(argv[2]),k_B=atof(argv[3]),Re=atof(argv[4]),
               ar=atof(argv[5]),perturb_v=argc==7?atof(argv[6]):0;

        // Grid and flag dimensions
        n=608;m=3*n;
        double flag_L=1,flag_w=1./ar,U=1,
               flag_a=flag_w*0.25,flag_ntrace=17;

        // Make output filename
        fn=fa=new char[strlen(argv[1])+128];
        argc==6?sprintf(fa,"fl_%s_%s_%s_%s",argv[2],argv[3],argv[4],argv[5])
               :sprintf(fa,"fl_%s_%s_%s_%s_%s",argv[2],argv[3],argv[4],argv[5],argv[6]);

        // Create flag
        sim_t=new sim_horiz_flow(U,-1,5,-1,1);
        obj.push_back(new obj_flag(0,0,flag_L,flag_w,flag_a,perturb_v,flag_ntrace,fa));

        // Set constants
        rhof=1;rhos=mu*rhof*flag_L/flag_w;
        visc=U*flag_L*rhof/Re;
        if(argc==7) {
            T=160;num_frames=80;
        } else {
            T=30;num_frames=900;
        }
        G=4*k_B/(flag_w*flag_w*flag_w);

        // Print diagnostic information
        printf("# Shear modulus                : %g (M/LT^2)\n"
               "# Flag width                   : %g dy\n"
               "# Solid density                : %g\n#\n",G,0.5*flag_w*n,rhos);
    } else {
        fprintf(stderr,"Simulation type '%s' not recognized\n",argv[1]);
        return true;
    }

    // Create a directory for output if one hasn't been created yet
    char *fout=new char[strlen(fn)+5];
    sprintf(fout,"%s.odr",fn);
    mkdir(fout,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Create default simulation type if none has been set up
    if(sim_t==NULL) sim_t=new sim_type();

    // Construct the simulation class and add objects
    fluid_2d f2d(m,n,*sim_t,visc,rhof,tw,pinw,ntrace,fflags,fout);
    f2d.solid_prd_bc=solid_prd_bc;
    mat_const mc(G,rhos,ev_trans_mult,set_velocity);
    for(unsigned int k=0;k<obj.size();k++) f2d.add_object(obj[k],mc);

    // Initialize the simulation and set the extra viscosity and timestep based
    // on the padding factors
    f2d.initialize(ev_mult,dt_pad,dt_ev_pad);

    // Run the simulation
    f2d.solve(T,num_frames);

    // Free dynamically allocated objects
    delete sim_t;
    for(unsigned int k=0;k<obj.size();k++) delete obj[k];
    delete [] fout;
    if(fa!=NULL) delete [] fa;
}
