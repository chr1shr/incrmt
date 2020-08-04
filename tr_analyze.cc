#include <cstdio>
#include <cstdlib>
#include <cstring>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

#include "common.hh"

int main(int argc, char **argv) {

    // Check for either two or three command line arguments
    if(argc!=2) {
        fputs("Syntax: ./tr_analyze <binary tracer file>\n",stderr);
        return 1;
    }

    // Open the file and read sets of 35 entries
    FILE *fp=fopen(argv[argc-1],"rb");
    std::vector<double*> p;
    int n=0,i;
    p.push_back(new double[35]);
    while(fread(p[n],sizeof(double),35,fp)==35) {
        p.push_back(new double[35]);
        n++;
    }

    // Perform Fourier transform of the flag tail marker's y position
    double T,s,c,z,ome,nor=1./static_cast<float>(n-3*n/4);
    for(T=0.02;T<10;T+=0.02) {
        ome=2*M_PI/T;
        s=c=0;
        for(int i=3*n/4;i<n;i++) {
            z=ome*p[i][0];
            s+=sin(z)*p[i][2];
            c+=cos(z)*p[i][2];
        }
        printf("%.14g %.14g\n",T,sqrt(s*s+c*c)*nor);
    }

    // Free dynamically allocated memory
    for(i=n;i>=0;i--) delete [] p[i];
}
