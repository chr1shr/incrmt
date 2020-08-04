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

int main(int argc,char **argv) {

    // Check for at least one command-line argument
    if(argc!=2) {
        fputs("Syntax: ./sediment <grid_size>\n",stderr);
        return 1;
    }

    // Parse the grid size
    int m=atoi(argv[1]),n=2*m,mn=m*n;
    if(m<=0||m>65536) {
        fputs("Grid size out of range\n",stderr);
        return 1;
    }

    // Some default settings
    const int num_frames=1600,      // Number of output frames
              ntrace=512;           // Number of tracers
    const double T=160,             // Simulation duration
                 tw=2.5,            // Transition width (in grid spacings)
                 pinw=2.5,          // Pinning width (in grid spacings)
                 dx=1.5/m,dy=3./n;  // Grid spacings
    const bool set_velocity=true;   // Whether solids initialize their own
                                    // velocity, or inherit the fluid's
                                    // velocity

    // File output flags
    // 1 - u                128 - solid Y
    // 2 - v                256 - solid phi
    // 4 - speed            512 - only subregions for solids
    // 8 - pressure         1024 - add ghost regions
    // 16 - vorticity       2048 - difference fields
    // 32 - density         4096 - divergence
    // 64 - solid X         8192 - parsable timing information
    unsigned int fflags=8|16|32|64|128|256|512;

    // Default fluid/solid parameters
    const double rhof=1,            // Fluid density
                 visc=2e-3,         // Fluid viscosity
                 G=1.3,             // Shear modulus

    // Specify padding factors for extra viscosity and timestep
                 ev_mult=0.8,       // The multiplier to apply to the extra viscosity
                 ev_trans_mult=1,   // The multiplier to apply to the
                                    // additional extra viscosity term in the
                                    // transition region
                 dt_pad=0.4,        // The padding to apply to the timestep
                                    // restriction due to physical terms
                 dt_ev_pad=0.8;     // The padding to apply to the timestep
                                    // restriction due to the extra viscosity

    // Set up different simulation examples
    char fout[64];
    sprintf(fout,"sed_%d.odr",m);
    mkdir(fout,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

    // Create vectors for the inserted objects and their corresponding material
    // properties
    std::vector<object*> obj;
    std::vector<mat_const> mc;

    // Create status array for finding object overlaps
    bool *f=new bool[mn],*fp;
    for(fp=f;fp<f+mn;fp++) *fp=false;

    // Constants about objects to add
    const int cnum=20;              // Total number of objects
    const double r=0.044,           // Half-width of object
                 le=0.32,           // Half-length of object
                 pad=2.5*sqrt(dx*dx+dy*dy); // Overlap padding

    // Randomly insert objects into the domain, checking for overlaps
    double x,y,rhos;
    int i,j,k=0;
    object* nobj;
    while(k<cnum) {

        // Pick a random position and density for the object
retry:  x=rnd(-0.75,0.75);
        y=rnd(-1.5,1.5);
        if(x+r>0.73||x-r<-0.73||y+le+r>1.48||y-le-r<-1.48) goto retry;
        rhos=rnd(0.4,1.6);

        // Create the candidate object and see if it overlaps with any
        // previously inserted object
        nobj=new obj_rounded_rod(x,y,r,le,(rhos-1)/rhos);
        for(fp=f,j=0;j<n;j++) {
            double py=-1.5+(j+0.5)*dy-y;
            for(i=0;i<m;i++,fp++) {
                double px=-0.75+(i+0.5)*dx-x;
                if(*fp&&nobj->phi(px,py)<pad) {delete nobj;goto retry;}
            }
        }

        // If the object doesn't overlap with a previous one, then update the
        // overlap status array
        for(fp=f,j=0;j<n;j++) {
            double py=-1.5+(j+0.5)*dy-y;
            for(i=0;i<m;i++,fp++) {
                double px=-0.75+(i+0.5)*dx-x;
                if(nobj->phi(px,py)<pad) *fp=true;
            }
        }

        // Store the object for later insertion into the simulation
        printf("Object %d: (x,y)=(%g,%g)\n",k++,x,y);
        obj.push_back(nobj);
        mc.push_back(mat_const(G,rhos,ev_trans_mult,set_velocity));
    }

    // Construct the simulation class and add objects
    sim_type sim_t=sim_type(-0.75,0.75,-1.5,1.5,false,false);
    fluid_2d f2d(m,n,sim_t,visc,rhof,tw,pinw,ntrace,fflags,fout);
    for(unsigned int k=0;k<obj.size();k++) f2d.add_object(obj[k],mc[k]);

    // Initialize the simulation and set the extra viscosity and timestep based
    // on the padding factors
    f2d.initialize(ev_mult,dt_pad,dt_ev_pad);

    // Run the simulation
    f2d.solve(T,num_frames);

    // Free dynamically allocated objects and memory
    for(unsigned int k=0;k<obj.size();k++) delete obj[k];
    delete [] f;
}
