#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "common.hh"
#include "tgmg.hh"
#include "mgs_mac.hh"
#include "mgs_fem.hh"

// Grid dimensions
const double ax=-1,bx=1,ay=-1,by=1;

// Safe grid limit
const int grid_limit=32768;

// Minimum time for testing V-cycles
const double min_wtime=1.;

// Template for setting up the multigrid hierarchy and solving the problem
template<class T>
void setup_and_do(T& mg,bool solve,double dx,double dy) {
    mg.setup();

    if (solve) {
        mg.solve_v_cycle();

        // Output the fields using the Gnuplot binary matrix format
        mg.output_b("b.0",ax,dx,ay,dy);
        mg.output_z("z.0",ax,dx,ay,dy);
        mg.output_res("r.0",ax,dx,ay,dy);
    } else {
        double t0=wtime(),t1;
        int k=0;

        // Perform V-cycles until the testing time has been reached
        do {
            mg.v_cycle();
            k++;t1=wtime();
        } while(t1<t0+min_wtime);
        printf("%d V-cycles performed in %g s\n"
               "%g ms per V-cycle\n",k,t1-t0,1e3*(t1-t0)/k);
    }
}

int main(int argc,char **argv) {

    // Periodicity in the x and y directions
    bool x_prd=false,y_prd=false;

    // Check for periodicity flags
    int c=1;
    while (c+2<argc) {
        if(strcmp(argv[c],"-px")==0) x_prd=true;
        else if(strcmp(argv[c],"-py")==0) y_prd=true;
        else break;
        c++;
    }

    // Check for the correct number of command-line arguments
    if(c+2>argc||c+3<argc) {
        fputs("Syntax: ./mg_test [-px] [-py] <case> <x gridpoints> [<y gridpoints>]\n\n"
              "Case=0 : MAC correction (solve)\n"
              "Case=1 : MAC correction (time V-cycle)\n"
              "Case=2 : Finite-element projection (solve)\n"
              "Case=3 : Finite-element projection (time V-cycle)\n\n"
              "The -px and -py flags enable periodicity in the x and y directions.\n"
              "If y gridpoints aren't specified, then the code uses a square grid.\n",stderr);
        return 1;
    }

    // Check case number
    int ca=atoi(argv[c]);
    if (ca < 0 or ca > 3) {
        fputs("Case number should be between 0 and 3\n",stderr);
        return 1;
    }

    // Check horizontal grid size
    int m=atoi(argv[c+1]);
    if(m<=8||m>grid_limit) {
        fputs("x gridpoints out of range\n",stderr);
        return 1;
    }

    // Check vertical grid size
    int n;
    if(argc==c+2) n=m;
    else {
        n=atoi(argv[c+2]);
        if(n<=8||n>grid_limit) {
            fputs("y gridpoints out of range\n",stderr);
            return 1;
        }
    }

    // Set constants and allocate memory for the grids
    int ij,mn=m*n,i1=m/5,j1=n/3,i2=(6*m)/7,j2=(2*n)/3;
    double dx=(bx-ax)/(x_prd?m:m-1),
           dy=(by-ay)/(y_prd?n:n-1),
           *b=new double[mn],*z=new double[mn];

    // Set up the RHS array, putting two delta function sources in the grid
    for(ij=0;ij<mn;ij++) b[ij]=0;
    b[i1+j1*m]=1;
    b[i2+j2*m]=-1;

    // Solve the multigrid problem
    if(ca<2) {
        mgs_mac_const_rho ms(m,n,x_prd,y_prd,dx,dy,b);
        setup_and_do(ms.mg,ca==0,dx,dy);
    } else {
        mgs_fem_const_rho ms(m,n,x_prd,y_prd,dx,dy,b);
        setup_and_do(ms.mg,ca==2,dx,dy);
    }

    // Delete dynamically allocated memory
    delete [] z;
    delete [] b;
}
