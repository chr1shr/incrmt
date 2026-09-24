#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>

// Activate for Windows
#if (defined(_WIN32) || defined(__WIN32__))
    #include <direct.h>
#endif
// test change friction branch

#include "common.hh"
#include "object.hh"
#include "obj_field.hh"
#include "fluid_2d.hh"

#if (defined(_WIN32) || defined(__WIN32__))
    #define mkdir(A, B) _mkdir(A)
#endif

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
              "circle_hoop     circle moving inside hoop\n"
              "jellyfish_2D    circle in square deforming\n"
              "square          rotating square\n"
              "rod             rounded rod experiencing constant fluid velocity\n"
              "rvemrve         U shape being acclerated to the left \n"
              "UcurveMany      Aseveral U shapes being acclerated to the left \n"              
              "accordiansingle one accordian being acclerated to the left \n"
              "accordian_with_circles   multiple accordian and circles acclerated to the left \n"       
              "twocircles      two circles part of the same levelset\n" 
              "seven-star      rotating seven-pointed star\n"
              "piston          flexible paddle moving through cavity\n"
              "fluid           fluid only\n"
              "full            solid only\n"
              "flapper         single flapping swimmer\n"
              "flappers        multiple flapping swimmers\n"
              "lamprey         a swimming lampry\n"
              "hoop            hoop acted on by gravity\n"              
              "conc-hoops      concentric hoops, inside one is spinning\n"
              "roll-hoop       hoop squased between two rectangles that are shearing\n"   
              "multi-square    two squares moving past each other\n"
              "2-rod_push     two vertically sliding rectangles\n"
              "2-rod_push_one_levelset     two vertically sliding rectangles made from one levelset\n"
              "Square_ramp_slide    Square on a plane being pushed with gravity at an angle\n"
              "multi-drop      multiple objects sedimenting\n"
              "flag            flapping flag\n",stderr);
        return 1;
    }

    // Some default settings
    int m=200,n=200;             // Simulation resolution
    double T=4;                  // Simulation duration
    int num_frames=200;          // Number of output frames
    double tw=2.5;               // Transition width (in grid spacings)
    double pinw=tw;             // Pinning width (in grid spacings)
    int ntrace=512;              // Number of tracers
    bool set_velocity=true;      // Whether solids initialize their own velocity, or
                                 // inherit the fluid's velocity
    bool solid_prd_bc=false;     // Whether to apply periodic boundary conditions to
                                 // the reference map fields

    // File output flags
    // 1 - u                128 - solid Y
    // 2 - v                256 - solid phi
    // 4 - speed            512 - only subregions for solids - make sure to turn off if making grid plots of reference map
    // 8 - pressure         1024 - add ghost regions
    // 16 - vorticity       2048 - difference fields
    // 32 - density         4096 - divergence
    // 64 - solid X         8192 - parsable timing information
    // 16384 - Diagnostic fields, gradphi(xy), contdiv(xy), alpha, sigma terms
    // 32768 - Friction Diagnosics sigma_nt, friction forces, combined phi
    // 65536 - phi fields of 2 objects and estimated phi field and weight function

    unsigned int fflags= 1|2|8|16|64|128|256|16384|32768|65536; //1|2|4|8|16|64|128|256|16384|32768|65536; //1|2|4|256|16384|32768|65536 //1|2|4|8|16|64|128|256; //8|16|64|128|256|512;

    // Default fluid/solid parameters
    double rhof=1;            // Fluid density
    double rhos=1;            // Solid density
    double visc=1e-3;         // Fluid viscosity
    double G=24; //24;        // Shear modulus
    double mufr = 0; //24;     // Friction coefficient
    double solid_visc = 0; // solid viscosity, note that this will only apply if bigger than the extra viscosity, rho_s * shear wave speed

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
        G=24;visc=1e-2;rhos=3;
        T=4*M_PI;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_simple_spin(0.8,0,0.5,0.15,M_PI,1,fa));
        //obj.push_back(new obj_simple_spin(0.8,0.5,0.125,0.2,M_PI,1,fa));
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

    } else if(se(argv[1],"circle_hoop_100")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_100","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5* 0.5;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"circle_hoop_150")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_150","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5* 1.5/2.0;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"circle_hoop_200")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_200","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"circle_hoop_300")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_300","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5*1.5;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"circle_hoop_400")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_400","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5*2;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"circle_hoop_500")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_500","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5*5.0/2.0;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"circle_hoop_600")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_600","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5*3;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"circle_hoop_800")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_800","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5*4;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"circle_hoop_1000")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_1000","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5*5;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"circle_hoop_1200")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("circle_hoop_1200","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 5.0;
        double ywidth = xwidth;
        m=n=atoi(argv[2]);check_grid(m);
        G=60;rhos=100;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m;
        solid_visc = 0; //solid viscosity
        T=50;num_frames=1000;
        tw = 2.5*6;

        // circle and outer radii
        double c_in = 2.64;
        double r_out = 4.0;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_circle_hoop(c_in,3.0,r_out,-3*M_PI*0.2/100));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);

    } else if(se(argv[1],"square")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("square","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=3;
        num_frames=1;

        // Create object
        obj.push_back(new obj_square(-0.02,0,0.5,0.01,0,2));

    } else if(se(argv[1],"rod")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("rod","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=atoi(argv[2]);check_grid(m);
        n=ceil(8/0.4*m);check_grid(n);
        rhos=1;visc=1e-3;
        T=3.;num_frames=5;

        // Create object
        obj.push_back(new obj_rounded_rod(0,0,0.05,2,2));
        // Set initial velocity to be constant velocity
        sim_t=new sim_horiz_half_flow(0,-0.2,0.2,-4,4);

    } else if(se(argv[1],"multi-square")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("multi-square","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 0.6;
        double ywidth = 1.2;
        m=atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        rhos=1;visc=1e-3; //solid density and viscosity
        T=1;num_frames=500;
        // Create object
        obj.push_back(new obj_square(0,0.51,0.5,0.2));
        obj.push_back(new obj_square(0,-0.51,0.5,0.2));

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);

  

    } else if(se(argv[1],"roll-hoop"))  {
        // Meant for 400 grid cells
        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("roll-hoop","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 1.5*6.0;
        double ywidth = 3.15;
        m=1.5*atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=600;rhos=200;
        mufr=0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=200;num_frames=10000; //T=26;num_frames=5400;
        // tw=4; pinw=4; //make transition zone larger

        // Create object
        obj.push_back(new obj_tworectangles(0,2.2875,5.475,0.6,0)); //two rectangles
        obj.push_back(new obj_hoop(0.7425,1.4925,0)); //hoop

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);  

    } else if(se(argv[1],"roll-hoop_mufr1"))  {
        // Meant for 400 grid cells
        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("roll-hoop_mufr1","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 1.5*6.0;
        double ywidth = 3.15;
        m=1.5*atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=600;rhos=200;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=200;num_frames=10000; //T=26;num_frames=5400;
        // tw=4; pinw=4; //make transition zone larger

        // Create object
        obj.push_back(new obj_tworectangles(0,2.2875,5.475,0.6,0)); //two rectangles
        obj.push_back(new obj_hoop(0.7425,1.4925,0)); //hoop

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);  

    } else if(se(argv[1],"roll-hoop_mufr10"))  {
        // Meant for 400 grid cells
        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("roll-hoop_mufr10","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 1.5*6.0;
        double ywidth = 3.15;
        m=1.5*atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=600;rhos=200;
        mufr=10.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=200;num_frames=6000; //T=26;num_frames=5400;
        // tw=4; pinw=4; //make transition zone larger

        // Create object
        obj.push_back(new obj_tworectangles(0,2.2875,5.475,0.6,0)); //two rectangles
        obj.push_back(new obj_hoop(0.7425,1.4925,0)); //hoop

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);  

    } else if(se(argv[1],"roll-ucurve_right"))  {
        // Meant for 400 grid cells
        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("roll-roll_ucurve_right","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 1.5*6.0;
        double ywidth = 3.15;
        m=1.5*atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=600;rhos=200; //200 produced strange crash at 1233
        mufr=0.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=400;num_frames=20000; //T=26;num_frames=5400;
        // tw=4; pinw=4; //make transition zone larger

        // Create object
        obj.push_back(new obj_tworectangles(0,2.2875,7.0,0.6,0)); //two rectangles
        //obj.push_back(new obj_hoop(0.7425,1.4925,0)); //hoop
        obj.push_back(new obj_Ucurve(0,0,1.1175,0.375,3.4,0));

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);  

    } else if(se(argv[1],"roll-ucurve_mufr1_right"))  {
        // Meant for 400 grid cells
        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("roll-roll_ucurve_mufr1_right","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 1.5*6.0;
        double ywidth = 3.15;
        m=1.5*atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=600;rhos=200;
        mufr=1.0;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=400;num_frames=40000; //T=26;num_frames=5400;
        // tw=4; pinw=4; //make transition zone larger

        // Create object
        obj.push_back(new obj_tworectangles(0,2.2875,7.0,0.6,0)); //two rectangles
        //obj.push_back(new obj_hoop(0.7425,1.4925,0)); //hoop
        obj.push_back(new obj_Ucurve(0,0,1.1175,0.375,3.4,0));

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);  

    } else if(se(argv[1],"Square_ramp_slide"))  {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("Square_ramp_slide","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 6.0;
        double ywidth = 1.575;
        m=atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=600;rhos=50;
        mufr=0.5;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=450;num_frames=12000; //T=26;num_frames=5400;
        // tw=4; pinw=4; //make transition zone larger

        // Create object
        obj.push_back(new obj_rectangle(0,-0.75,5.475,0.3,0,0,true)); //bottom plate
        obj.push_back(new obj_rectangle(0,0.75,3*0.6,0.3,1.5,0.02*M_PI,false)); //push down top square

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);  

    } else if(se(argv[1],"Square_ramp_slide_N120"))  { // designed for N=120

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("Square_ramp_slide_N120","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 4.0;
        double ywidth = 1.75;
        m=atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=500;rhos=50;
        mufr=0.5;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=450;num_frames=12000; //T=26;num_frames=5400;
        // tw=4; pinw=4; //make transition zone larger

        // Create object
        obj.push_back(new obj_rectangle(0,-0.833,3.417,0.333,0,0,true)); //bottom plate
        obj.push_back(new obj_rectangle(-1.5,0.833,1.333,0.333,1.0,0.02*M_PI,false)); //push down top square

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);  

    } else if(se(argv[1],"Square_ramp_slide_N300"))  { // designed for N=120

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("Square_ramp_slide_N300","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 4.0;
        double ywidth = 1.75;
        m=atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=500;rhos=50;
        mufr=0.5;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=450;num_frames=12000; //T=26;num_frames=5400;
        // tw=4; pinw=4; //make transition zone larger

        // Create object
        obj.push_back(new obj_rectangle(0,-0.833,3.417,0.333,0,0,true)); //bottom plate
        obj.push_back(new obj_rectangle(-1.5,0.833,1.333,0.333,1.0,0.02*M_PI,false)); //push down top square

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);  

    } else if(se(argv[1],"Square_ramp_slide_N440"))  {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("Square_ramp_slide_N440","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 4.0;
        double ywidth = 1.75;
        m=atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=500;rhos=50;
        mufr=0.5;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=450;num_frames=12000; //T=26;num_frames=5400;
        // tw=4; pinw=4; //make transition zone larger

        // Create object
        obj.push_back(new obj_rectangle(0,-0.833,3.417,0.333,0,0,true)); //bottom plate
        obj.push_back(new obj_rectangle(-1.5,0.833,1.333,0.333,1.0,0.02*M_PI,false)); //push down top square

        sim_t=new sim_horiz_flow(0,-xwidth,xwidth,-ywidth,ywidth);  

    } else if(se(argv[1],"Ucurve_Friction_testing")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("Ucurve_Friction_testing","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 0.18; //0.17667;
        m=atoi(argv[2]);check_grid(m);
        n=ceil(1/(2*xwidth)*m);check_grid(n);
        G=6;rhos=1;visc=1e-3; //solid density and viscosity
        mufr=0.5;
        T=4;num_frames=500;
        // Create object
        obj.push_back(new obj_Ucurve(0,0,0.1,0.14,0.7,2));
        sim_t=new sim_type(-xwidth,xwidth,-0.12,0.88,false,false);

    } else if(se(argv[1],"Ucurve")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("Ucurve","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth = 7.0;
        double ywidth = 2.5;
        m=atoi(argv[2]);check_grid(m);
        n=ceil(ywidth/xwidth*m);check_grid(n);
        G=600;rhos=200;
        mufr=0.5;
        visc=0.7*sqrt(G/rhos)*2*xwidth/m; //fluid viscosity
        solid_visc = 0; //solid viscosity
        T=500;num_frames=5000; //T=26;num_frames=5400;
        // Create object
        obj.push_back(new obj_Ucurve(0,0,1.25,1,8.25,1.0));
        sim_t=new sim_type(-xwidth,xwidth,-ywidth,ywidth,false,false);
    } else if(se(argv[1],"accordian_with_circles")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("accordian_with_circles","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth=1;
        double yheight=1.3;
        n=atoi(argv[2]);check_grid(n);
        m=ceil(xwidth/yheight*n);check_grid(m);
        G=20;rhos=1.0;rhof=1.0;visc=1e-2; //solid density and viscosity
        mufr = 0;
        T=80.;num_frames=40000;
        // Create objects
        obj.push_back(new obj_accordianmultiple(-0.5,0,1.2,0.14,0.4,-1.0));
        // obj.push_back(new obj_accordianmultiple(0, 0, 1.2, 0.14, 0.4, -0.9));
        obj.push_back(new obj_accordianmultiple(0.5,0,1.2,0.14,0.4,-1.0));
        obj.push_back(new obj_circle(0.7, 0.4, 0.15, -1.0));
        obj.push_back(new obj_circle(-0.7, 0, 0.15, -1.0));
        obj.push_back(new obj_circle(0.7, -0.4, 0.15, -1.0));
        obj.push_back(new obj_circle(-0.3, 0.4, 0.15, -1.0));
        obj.push_back(new obj_circle(0.3, 0, 0.15, -1.0));
        obj.push_back(new obj_circle(-0.3, -0.4, 0.15, -1.0));        
        ///obj.push_back(new obj_circle(0.4, -0.79, 0.15, -0.9));
        ///obj.push_back(new obj_circle(0, -0.65, 0.15, -0.9));
        ///obj.push_back(new obj_circle(-0.4, -0.79, 0.15, -0.9));    

        // Set initial velocity to be constant velocity
        sim_t=new sim_horiz_half_flow(0,-xwidth,xwidth,-yheight,yheight);

    } else if(se(argv[1],"accordian_with_circles_600")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("accordian_with_circles_600","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth=1;
        double yheight=1.3;
        n=atoi(argv[2]);check_grid(n);
        m=ceil(xwidth/yheight*n);check_grid(m);
        G=10;rhos=1.0;rhof=1.0;visc=5e-3; //solid density and viscosity
        mufr = 0;
        T=80.;num_frames=40000;
        // Create objects
        obj.push_back(new obj_accordianmultiple(-0.5,0,1.2,0.14,0.4,-1.0));
        // obj.push_back(new obj_accordianmultiple(0, 0, 1.2, 0.14, 0.4, -0.9));
        obj.push_back(new obj_accordianmultiple(0.5,0,1.2,0.14,0.4,-1.0));
        obj.push_back(new obj_circle(0.7, 0.4, 0.15, -1.0));
        obj.push_back(new obj_circle(-0.7, 0, 0.15, -1.0));
        obj.push_back(new obj_circle(0.7, -0.4, 0.15, -1.0));
        obj.push_back(new obj_circle(-0.3, 0.4, 0.15, -1.0));
        obj.push_back(new obj_circle(0.3, 0, 0.15, -1.0));
        obj.push_back(new obj_circle(-0.3, -0.4, 0.15, -1.0));        
        ///obj.push_back(new obj_circle(0.4, -0.79, 0.15, -0.9));
        ///obj.push_back(new obj_circle(0, -0.65, 0.15, -0.9));
        ///obj.push_back(new obj_circle(-0.4, -0.79, 0.15, -0.9));    

        // Set initial velocity to be constant velocity
        sim_t=new sim_horiz_half_flow(0,-xwidth,xwidth,-yheight,yheight);

   

    } else if(se(argv[1],"accordian_with_circles_mufr10")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("accordian_with_circles_mufr10","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth=1;
        double yheight=1.3;
        n=atoi(argv[2]);check_grid(n);
        m=ceil(xwidth/yheight*n);check_grid(m);
        G=5;rhos=1.0;rhof=1.0;visc=1e-3; //solid density and viscosity
        mufr = 1.0;
        T=80.;num_frames=40000;
        // Create objects
        obj.push_back(new obj_accordianmultiple(-0.5,0,1.2,0.14,0.4,-1.0));
        // obj.push_back(new obj_accordianmultiple(0, 0, 1.2, 0.14, 0.4, -0.9));
        obj.push_back(new obj_accordianmultiple(0.5,0,1.2,0.14,0.4,-1.0));
        obj.push_back(new obj_circle(0.7, 0.4, 0.15, -1.0));
        obj.push_back(new obj_circle(-0.7, 0, 0.15, -1.0));
        obj.push_back(new obj_circle(0.7, -0.4, 0.15, -1.0));
        obj.push_back(new obj_circle(-0.3, 0.4, 0.15, -1.0));
        obj.push_back(new obj_circle(0.3, 0, 0.15, -1.0));
        obj.push_back(new obj_circle(-0.3, -0.4, 0.15, -1.0));     
        ///obj.push_back(new obj_circle(0.4, -0.79, 0.15, -0.9));
        ///obj.push_back(new obj_circle(0, -0.65, 0.15, -0.9));
        ///obj.push_back(new obj_circle(-0.4, -0.79, 0.15, -0.9));    

        // Set initial velocity to be constant velocity
        sim_t=new sim_horiz_half_flow(0,-xwidth,xwidth,-yheight,yheight);
    } else if(se(argv[1],"twocircles")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("twocircles","<grid resolution>");

        double Rc = 50;
        double yheight = 102.0;
        double xwidth = 0.75*yheight;

        // Set grid resolution and various simulation constants
        n=atoi(argv[2]);check_grid(n);
        m=ceil(xwidth/yheight*n);check_grid(m);
        // radii and seperation of disks
        double seperation = Rc + 1.3*2.5*2*yheight/n;
        solid_visc = 500;
        rhos=15;
        G=solid_visc*solid_visc/rhos;
        mufr=0.5;
        visc=0.7*solid_visc/rhos*2*yheight/n; //fluid viscosity
        rhof = 1;
        double grav = 0; //0.01*G/rhos*1/Rc;
        T=16;num_frames=4000;
        // Create objects
        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_twocircles(Rc,seperation,grav,0,2*yheight/n));
        
        // Set initial velocity to be constant velocity
        sim_t=new sim_horiz_half_flow(0,-xwidth,xwidth,-yheight,yheight);



    } else if(se(argv[1],"twocircles_shear_1021_slide_018")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("twocircles_shear_1021_slide_018","<grid resolution>");

        double Rc = 50;
        double yheight = 102.0;
        double xwidth = 0.75*yheight;

        // Set grid resolution and various simulation constants
        n=atoi(argv[2]);check_grid(n);
        m=ceil(xwidth/yheight*n);check_grid(m);
        // radii and seperation of disks
        double seperation = Rc + 1.3*2.5*2*yheight/n;
        solid_visc = 5000; //500;
        rhos = 500;
        G=solid_visc*solid_visc/rhos;
        mufr=0.5;
        visc=0.1; //0.7*solid_visc/rhos*2*yheight/n; //viscosity
        rhof = 1;
        double grav = 0.01*G/rhos*1/Rc;
        T=60;num_frames=6000;
        // Create objects
        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_twocircles(Rc,seperation,grav,0.18,2*yheight/n)); //
        
        // Set initial velocity to be constant velocity
        sim_t=new sim_horiz_half_flow(0,-xwidth,xwidth,-yheight,yheight);


    } else if(se(argv[1],"hoop")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("hoop","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=4;visc=1e-3;
        T=3.;num_frames=200;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_hoop(0.625,0.05,1));
        sim_t=new sim_type(-1,1,-1,1,false,false,1,0.5);

    } else if(se(argv[1],"conc-hoops")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("conc-hoops","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=4;visc=1e-3;
        T=3.;num_frames=200;

        // Create object and set non-periodic boundary conditions
        obj.push_back(new obj_hoop(0.625,0.05,1));
        obj.push_back(new obj_hoop(0.625,0.05,1));
        sim_t=new sim_type(-1,1,-1,1,false,false,1,0.5);

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
        visc=1e-5;
        T=1;num_frames=500;
        ntrace=256; 

        // Set initial velocity to be the "pulses" field

        sim_t=new sim_velocity_pulses;

    } else if(se(argv[1],"full")) {

        // Check for the right number of command-line arguments
        if(argc!=3) cmd_args("full","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=1.;G=3.;visc=1e-5;
        T=1;num_frames=200;

        // Switch off tracers, since there is no fluid to initialize them in.
        // Set several special flags appropriate for solid-only simulation.
        ntrace=0;
        set_velocity=false;
        solid_prd_bc=true;

        // Create full solid, at set initial velocity to be the "pulses" field
        obj.push_back(new obj_full());
        sim_t=new sim_velocity_pulses;

    } else if(se(argv[1],"flappers_n12")) {

        if(argc!=3) cmd_args("flappers_n12","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=3;G=5;visc=1e-4;
        mufr=0.0; // friction coefficient

        T=48;num_frames=12000;

        // Insert non-overlapping flappers
        const int cnum=12;
        double x,y,wx=0.25+0.026+0.01,wy=0.026+0.01,dx,dy;

        int i=0,j;
        bool good;

        srand(42); // same random seed for all trials

        while(i<cnum) {

            // Choose random position for flapper and check it is in range
            x=rnd(-1,1);
            y=rnd(-1,1);
            if(x+wx>0.99 || x-wx<-0.99 || y+wy>0.99 || y-wy<-0.99) continue;
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
                double theta = (rand() % 2 == 0) ? 0 : M_PI;
                obj.push_back(new obj_flapper(x,y,0.026,0.25,0.021,0.14,theta,2.2,rnd(7,9)));
            }
        }

        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-1.0,1.0,-1.0,1.0,false,false);

    } else if(se(argv[1],"flappers_n16")) {

        if(argc!=3) cmd_args("flappers","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=3;G=5;visc=1e-4;
        mufr=0.0; // friction coefficient

        T=48;num_frames=12000;

        // Insert non-overlapping flappers
        const int cnum=16;
        double x,y,wx=0.25+0.026+0.01,wy=0.026+0.01,dx,dy;

        int i=0,j;
        bool good;

        srand(42); // same random seed for all trials

        while(i<cnum) {

            // Choose random position for flapper and check it is in range
            x=rnd(-1,1);
            y=rnd(-1,1);
            if(x+wx>0.99 || x-wx<-0.99 || y+wy>0.99 || y-wy<-0.99) continue;
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
                double theta = (rand() % 2 == 0) ? 0 : M_PI;
                obj.push_back(new obj_flapper(x,y,0.026,0.25,0.021,0.14,theta,2.2,rnd(7,9)));
            }
        }

        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-1.0,1.0,-1.0,1.0,false,false);

    } else if(se(argv[1],"flappers_n20")) {

        if(argc!=3) cmd_args("flappers","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=3;G=5;visc=1e-4;
        mufr=0.0; // friction coefficient

        T=48;num_frames=12000;

        // Insert non-overlapping flappers
        const int cnum=20;
        double x,y,wx=0.25+0.026+0.01,wy=0.026+0.01,dx,dy;

        int i=0,j;
        bool good;

        srand(42); // same random seed for all trials

        while(i<cnum) {

            // Choose random position for flapper and check it is in range
            x=rnd(-1,1);
            y=rnd(-1,1);
            if(x+wx>0.99 || x-wx<-0.99 || y+wy>0.99 || y-wy<-0.99) continue;
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
                double theta = (rand() % 2 == 0) ? 0 : M_PI;
                obj.push_back(new obj_flapper(x,y,0.026,0.25,0.021,0.14,theta,2.2,rnd(7,9)));
            }
        }

        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-1.0,1.0,-1.0,1.0,false,false);

    } else if(se(argv[1],"flappers_no_friction")) {

        if(argc!=3) cmd_args("flappers_no_friction","<grid resolution>");

        // Set grid resolution and various simulation constants
        m=n=atoi(argv[2]);check_grid(m);
        rhos=4;G=10;visc=5e-4;
        mufr=0.0; // friction coefficient

        T=30;num_frames=8000;

        // Insert non-overlapping flappers
        const int cnum=20;
        double x,y,wx=0.25+0.026+0.01,wy=0.026+0.01,dx,dy;

        int i=0,j;
        bool good;

        srand(42); // same random seed for all trials

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
                double theta = (rand() % 2 == 0) ? 0 : M_PI;
                obj.push_back(new obj_flapper(x,y,0.026,0.25,0.021,0.14,theta,2.2,rnd(7,9)));
            }
        }

        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-1.0,1.0,-1.0,1.0,false,false);

    } else if(se(argv[1],"flagella_n12_friction")) {

        if(argc!=3) cmd_args("flagella_n12_friction","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth=1.5;
        double yheight=1.5;
        m=atoi(argv[2]);check_grid(n); //input number of grid cells in x
        n=ceil(yheight/(xwidth)*m);check_grid(m);
        rhos=3;G=5;visc=1e-4;
        mufr=1.0; // friction coefficient

        T=48;num_frames=12000;

        // Insert non-overlapping flagellas
        const int cnum=12;
        double x,y,wx=0.45+0.03+0.005,wy=0.1+0.005,dx,dy;

        int i=0,j;
        bool good;

        srand(24); // same random seed for all trials

        while(i<cnum) {

            // Choose random position for flagella and check it is in range
            x=rnd(-1.5,1.5);
            y=rnd(-1.5,1.5);
            if(x+wx>0.99*1.5 || x-wx<-0.99*1.5 || y+wy>0.99*1.5 || y-wy<-0.99*1.5) continue;
            good=true;

            // Loop over previously added flagella and check for overlaps
            for(j=0;j<i;j++) {
                obj_flapper *oc=reinterpret_cast<obj_flapper*>(obj[j]);
                dx=x-oc->cx;dy=y-oc->cy;
                if(fabs(dx)<2*wx&&fabs(dy)<2*wy) {good=false;break;}
            }

            // If this flagella doesn't overlap with a previous one, then add it
            if(good) {
                i++;
                double theta = (rand() % 2 == 0) ? 0 : M_PI;
                 obj.push_back(new obj_flagella(x,y,0.03,0.35,0.1,0.027,0.32,theta,rnd(2.1,2.3),rnd(8,12)));
            }
        }
        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-xwidth,xwidth,-yheight,yheight,false,false);

    } else if(se(argv[1],"flagella_n16_friction")) {

        if(argc!=3) cmd_args("flagella_n16_friction","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth=1.5;
        double yheight=1.5;
        m=atoi(argv[2]);check_grid(n); //input number of grid cells in x
        n=ceil(yheight/(xwidth)*m);check_grid(m);
        rhos=3;G=5;visc=1e-4;
        mufr=1.0; // friction coefficient

        T=48;num_frames=12000;

        // Insert non-overlapping flagellas
        const int cnum=16;
        double x,y,wx=0.45+0.03+0.005,wy=0.1+0.005,dx,dy;

        int i=0,j;
        bool good;

        srand(24); // same random seed for all trials

        while(i<cnum) {

            // Choose random position for flagella and check it is in range
            x=rnd(-1.5,1.5);
            y=rnd(-1.5,1.5);
            if(x+wx>0.99*1.5 || x-wx<-0.99*1.5 || y+wy>0.99*1.5 || y-wy<-0.99*1.5) continue;
            good=true;

            // Loop over previously added flagella and check for overlaps
            for(j=0;j<i;j++) {
                obj_flapper *oc=reinterpret_cast<obj_flapper*>(obj[j]);
                dx=x-oc->cx;dy=y-oc->cy;
                if(fabs(dx)<2*wx&&fabs(dy)<2*wy) {good=false;break;}
            }

            // If this flagella doesn't overlap with a previous one, then add it
            if(good) {
                i++;
                double theta = (rand() % 2 == 0) ? 0 : M_PI;
                 obj.push_back(new obj_flagella(x,y,0.03,0.35,0.1,0.025,0.2,theta,rnd(2.0,2.2),rnd(8,12)));
            }
        }
        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-xwidth,xwidth,-yheight,yheight,false,false);

    } else if(se(argv[1],"flagella_n20_friction")) {

        if(argc!=3) cmd_args("flagella_n20_friction","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth=1.5;
        double yheight=1.5;
        m=atoi(argv[2]);check_grid(n); //input number of grid cells in x
        n=ceil(yheight/(xwidth)*m);check_grid(m);
        rhos=3;G=5;visc=1e-4;
        mufr=1.0; // friction coefficient

        T=48;num_frames=12000;

        // Insert non-overlapping flagellas
        const int cnum=20;
        double x,y,wx=0.45+0.03+0.005,wy=0.1+0.005,dx,dy;

        int i=0,j;
        bool good;

        srand(24); // same random seed for all trials

        while(i<cnum) {

            // Choose random position for flagella and check it is in range
            x=rnd(-1.5,1.5);
            y=rnd(-1.5,1.5);
            if(x+wx>0.99*1.5 || x-wx<-0.99*1.5 || y+wy>0.99*1.5 || y-wy<-0.99*1.5) continue;
            good=true;

            // Loop over previously added flagella and check for overlaps
            for(j=0;j<i;j++) {
                obj_flapper *oc=reinterpret_cast<obj_flapper*>(obj[j]);
                dx=x-oc->cx;dy=y-oc->cy;
                if(fabs(dx)<2*wx&&fabs(dy)<2*wy) {good=false;break;}
            }

            // If this flagella doesn't overlap with a previous one, then add it
            if(good) {
                i++;
                double theta = (rand() % 2 == 0) ? 0 : M_PI;
                 obj.push_back(new obj_flagella(x,y,0.03,0.35,0.1,0.025,0.2,theta,rnd(2.0,2.2),rnd(8,12)));
            }
        }
        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-xwidth,xwidth,-yheight,yheight,false,false);
    
    } else if(se(argv[1],"lamprey_n24")) {

        if(argc!=3) cmd_args("lamprey_n24","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth=1.0;
        double yheight=1.0;
        m=atoi(argv[2]);check_grid(n); //input number of grid cells in x
        n=ceil(yheight/(xwidth)*m);check_grid(m);
        rhos=3;G=5;visc=1e-4;
        mufr=0.0; // friction coefficient

        T=50;num_frames=10000;

        double length = 0.5;
        double head_length = 0.0;
        double width = 0.05; 

        const int cnum=24;
        double x,y,wx=1.02*(length+head_length),wy=1.02*width,dx,dy;

        int i=0,j;
        bool good;

        srand(60); // same random seed for all trials

        while(i<cnum) {

            // Choose random position for flapper and check it is in range
            x=rnd(0.99*(-xwidth+length),0.99*xwidth);
            y=rnd(0.99*(-yheight+width),0.99*(yheight-width));
            good=true;

            for(j=0;j<i;j++) {
                obj_lamprey *oc=reinterpret_cast<obj_lamprey*>(obj[j]);
                dx=x-oc->cx;dy=y-oc->cy;
                if(fabs(dx)<2*wx&&fabs(dy)<2*wy) {good=false;break;}
            }

            if(good) {
                i++;
                obj.push_back(new obj_lamprey(x,y,0,length,head_length,width,width,rnd(3,6),rnd(1.0,1.5)));
            }

        }

        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-xwidth,xwidth,-yheight,yheight,false,false);
    
    } else if(se(argv[1],"lamprey_n16")) {

        if(argc!=3) cmd_args("lamprey_n16","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth=1.0;
        double yheight=1.0;
        m=atoi(argv[2]);check_grid(n); //input number of grid cells in x
        n=ceil(yheight/(xwidth)*m);check_grid(m);
        rhos=3;G=5;visc=1e-4;
        mufr=0.0; // friction coefficient

        T=50;num_frames=10000;

        double length = 0.5;
        double head_length = 0.0;
        double width = 0.05; 

        const int cnum=16;
        double x,y,wx=1.02*(length+head_length),wy=1.02*width,dx,dy;

        int i=0,j;
        bool good;

        srand(60); // same random seed for all trials

        while(i<cnum) {

            // Choose random position for flapper and check it is in range
            x=rnd(0.99*(-xwidth+length),0.99*xwidth);
            y=rnd(0.99*(-yheight+width),0.99*(yheight-width));
            good=true;

            for(j=0;j<i;j++) {
                obj_lamprey *oc=reinterpret_cast<obj_lamprey*>(obj[j]);
                dx=x-oc->cx;dy=y-oc->cy;
                if(fabs(dx)<2*wx&&fabs(dy)<2*wy) {good=false;break;}
            }

            if(good) {
                i++;
                obj.push_back(new obj_lamprey(x,y,0,length,head_length,width,width,rnd(3,6),rnd(1.0,1.5)));
            }

        }

        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-xwidth,xwidth,-yheight,yheight,false,false);

    } else if(se(argv[1],"lamprey_n8")) {

        if(argc!=3) cmd_args("lamprey_n8","<grid resolution>");

        // Set grid resolution and various simulation constants
        double xwidth=1.0;
        double yheight=1.0;
        m=atoi(argv[2]);check_grid(n); //input number of grid cells in x
        n=ceil(yheight/(xwidth)*m);check_grid(m);
        rhos=3;G=5;visc=1e-4;
        mufr=0.0; // friction coefficient

        T=50;num_frames=10000;

        double length = 0.5;
        double head_length = 0.0;
        double width = 0.05; 

        const int cnum=8;
        double x,y,wx=1.02*(length+head_length),wy=1.02*width,dx,dy;

        int i=0,j;
        bool good;

        srand(60); // same random seed for all trials

        while(i<cnum) {

            // Choose random position for flapper and check it is in range
            x=rnd(0.99*(-xwidth+length),0.99*xwidth);
            y=rnd(0.99*(-yheight+width),0.99*(yheight-width));
            good=true;

            for(j=0;j<i;j++) {
                obj_lamprey *oc=reinterpret_cast<obj_lamprey*>(obj[j]);
                dx=x-oc->cx;dy=y-oc->cy;
                if(fabs(dx)<2*wx&&fabs(dy)<2*wy) {good=false;break;}
            }

            if(good) {
                i++;
                obj.push_back(new obj_lamprey(x,y,0,length,head_length,width,width,rnd(3,6),rnd(1.0,1.5)));
            }

        }

        // Initialize non-periodic boundary conditions
        sim_t=new sim_type(-xwidth,xwidth,-yheight,yheight,false,false);

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

    mat_const mc(G,rhos,mufr,ev_trans_mult,solid_visc,set_velocity);
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
