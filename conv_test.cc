#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector>
#include <string>

#include "common.hh"
#include "object.hh"
#include "obj_field.hh"
#include "fluid_2d.hh"

int main(int argc,char **argv) {

    // Check for one command-line argument; if not provided, then print out a
    // help message
    if (argc!=2) {
        fputs("Syntax: ./conv_test <test_number>\n\n"
              "Test type:\n"
              "0: fluid only (mu=1e-3)\n"
              "1: fluid only (mu=4e-3)\n"
              "2: solid only (rho_s=1)\n"
              "3: solid square (rho_s=1, mu=1e-3)\n"
              "4: solid circle (rho_s=3, mu=1e-3)\n"
              "5: solid circle (rho_s=1, mu=1e-3)\n"
              "6: solid circle (rho_s=1, mu=4e-3)\n"
              "7: solid circle (rho_s=1, mu=1.6e-2)\n\n"
              "Additional options, assembled bitwise:\n"
              "8: constant extra viscosity\n"
              "16: extra viscosity timestep multiplier halved\n"
              "32: extra viscosity x1.5\n"
              "64: increase max resolution from 2520 to 5040\n",stderr);
        return 1;
    }

    // Read test number
    int type=atoi(argv[1]),ctype=type&7;
    if(type<0||type>127) {
        fputs("Test number out of range\n",stderr);
        return 1;
    }

    // Some default settings
    const double tw=2.5,        // Transition width (in grid spacings)
                 pinw=2.5;      // Pinning width (in grid spacings)
    const int ntrace=512;  // Number of tracers

    // File output flags
    // 1 - u                128 - solid Y
    // 2 - v                256 - solid phi
    // 4 - speed            512 - only subregions for solids
    // 8 - pressure         1024 - add ghost regions
    // 16 - vorticity       2048 - difference fields
    // 32 - density         4096 - divergence
    // 64 - solid X         8192 - parsable timing information
    unsigned int fflags=8|2048|(ctype==0||ctype==1?0:64|128|256);

    // Default fluid/solid parameters
    double rhof=1,                                              // Fluid density
           rhos[8]={1,1,1,1,3,1,1,1},                           // Solid density
           visc[8]={1e-3,4e-3,1e-20,1e-3,1e-3,1e-3,4e-3,16e-3}, // Fluid viscosity
           G=1,                                                 // Shear modulus
           ev_trans_mult=1;

    // Specify padding factors for extra viscosity and timestep
    double ev_mult=0.4,dt_pad=0.4,dt_ev_pad=0.8;

    // Try some alternative padding factors and timesteps
    if(type&16) dt_ev_pad=0.4;
    if(type&32) ev_mult=0.6;

    // Initialize objects. Note that currently this only works for objects that
    // don't require direct access to the simulation itself.
    std::vector<object*> obj;
    switch(ctype) {
        case 0: case 1:
          break;
        case 2:
          obj.push_back(new obj_full());break;
        case 3:
          obj.push_back(new obj_square(-0.1,0,0.6,0,0,0));break;
        case 4: case 5: case 6: case 7:
          obj.push_back(new obj_circle(-0.1,0,0.6,0,0,0));
    }

    // Set up list of resolutions to examine. The coarser grids should all be
    // factors of the finest grid.
    const int nsteps=50,nsim_max=12;
    const double interval=0.02;
    int nsim=9;
    int mtab1[9]={2520,1260,840,630,504,420,360,315,280},
        mtab2[12]={5040,2520,1680,1260,1008,840,720,630,560,504,420,360},
        *mtab=mtab1;
    if(type&64) {nsim=12;mtab=mtab2;}

    // Allocate per-simulation arrays
    char buf[nsim_max][32],cbuf[32];
    fluid_2d* f2d[nsim_max];
    int i=0,j,j0,k,l[nsim_max];
    double t0,t1,adt[nsim_max],err[32*(nsteps+1)*(nsim_max-1)],*ep=err;

    // Use initial velocity field that has six vortices with different length
    // scales
    sim_velocity_pulses sim_t;
    mat_const mc(G,rhos[ctype],ev_trans_mult,false);

    // Construct the simulation classes
    for(k=0;k<nsim;k++) {

        // Create output directory
        sprintf(buf[k],"cnv_%02d_%d.odr",type,k);
        mkdir(buf[k],S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);

        // Create simulation and add objects. For the full solid simulation,
        // disable the tracers and enable the solid periodic BCs.
        f2d[k]=new fluid_2d(mtab[k],mtab[k],sim_t,visc[ctype],
                            rhof,tw,pinw,ctype==2?0:ntrace,fflags,buf[k]);
        for(unsigned int o=0;o<obj.size();o++) f2d[k]->add_object(obj[o],mc);
        if(ctype==2) f2d[k]->solid_prd_bc=true;

        // Initialize the simulation. For if the 8 bit is set in the type,
        // scale the extra viscosity to match the coarsest grid.
        f2d[k]->initialize(type&8?(ev_mult*mtab[k]/mtab[nsim-1]):ev_mult,
                           dt_pad*mtab[nsim-1]/mtab[k],dt_ev_pad*mtab[nsim-1]/mtab[k]);
        l[k]=f2d[k]->timestep_select(interval,adt[k]);
        f2d[k]->save_header(interval*nsteps,nsteps);
        printf("# m=%d, steps=%d, dt=%g\n#\n",mtab[k],l[k],adt[k]);
    }

    // Open the output file
    sprintf(cbuf,"cnv_%02d.dat",type);
    FILE *fp=safe_fopen(cbuf,"w");

    // Step all simulations forward simultaneously and compute the L-norms
    // of velocity and pressure to the finest grid
    while(true) {

        // Perform the simulation output, and compute the L-norm comparisons
        fprintf(fp,"%d %g",i,i*interval);
        for(k=1;k<nsim;k++) {
          f2d[k]->l_comparison(*(*f2d),ep);
          for(double *ee=ep+32;ep<ee;ep+=8)
            fprintf(fp," %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g",
                    *ep,ep[1],ep[2],ep[3],ep[4],ep[5],ep[6],ep[7]);
          f2d[k]->write_files(i);
        }
        fputc('\n',fp);
        fflush(fp);

        // Perform the integration timesteps for each simulation
        printf("Timestep %d\n",i);
        if(i==nsteps) break;
        t0=wtime();
        for(k=0;k<nsim;k++) {

          // Perform the integration timesteps for the kth simulation, and output
          // diagnostic information every eight steps.
          for(j0=j=0;j<l[k];j++) {
            printf("Grid %d, step %d of %d\r",mtab[k],j,l[k]);fflush(stdout);
            f2d[k]->step_forward(adt[k]);
            if(j>0&&(j%8==0||j==l[k]-1)) {
              t1=wtime();
              printf("Grid %d, step %d of %d [%.6g s/step] {MAC %.2f,FEM %.2f}\n",
                     mtab[k],j,l[k],(t1-t0)/(j-j0),f2d[k]->ms_mac->avg_v_cycles(),
                     f2d[k]->ms_fem->avg_v_cycles());
              t0=t1;j0=j;
            }
          }
        }
        i++;
    }
    fclose(fp);

    // Output the second file, with the results ordered according to the grid
    // sizes
    sprintf(cbuf,"cnv_%02d.dig",type);
    fp=safe_fopen(cbuf,"w");
    for(k=1;k<nsim;k++) {
        fprintf(fp,"%d %g",mtab[k],2./mtab[k]);
        for(i=0;i<=nsteps;i++) {
            ep=err+32*((nsim-1)*i+(k-1));
            for(double *ee=ep+32;ep<ee;ep+=8)
                fprintf(fp," %.14g %.14g %.14g %.14g %.14g %.14g %.14g %.14g",
                        *ep,ep[1],ep[2],ep[3],ep[4],ep[5],ep[6],ep[7]);
        }
        fputc('\n',fp);
    }

    // Close the output file and free dynamically allocated objects
    fclose(fp);
    for(k=nsim-1;k>=0;k--) delete f2d[k];
    for(unsigned int k=0;k<obj.size();k++) delete obj[k];
}
