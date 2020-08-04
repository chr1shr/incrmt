#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "common.hh"

int main(int argc, char **argv) {

    // Check for either two or three command line arguments
    if(argc<2||argc>3) {
        fputs("Syntax: ./tr_unpack [-t] <binary tracer file>\n",stderr);
        return 1;
    }

    // If three are provided, then the second one needs to be "-t"
    if(argc==3) {
        if(strcmp(argv[1],"-t")!=0) {
            fputs("Error reading command-line arguments\n",stderr);
        }
    }

    // Open the file and read sets of 35 entries
    FILE *fp=fopen(argv[argc-1],"rb");
    double buf[35],*bp,ltime=0.;
    while(fread(buf,sizeof(double),35,fp)==35) {
        if(argc==2) {
            for(bp=buf; bp<buf+33; bp+=3)
                printf("%g %g %g ",*bp,bp[1],bp[2]);
            printf("%g %g\n",buf[33],buf[34]);
        }
        ltime=*buf;
    }

    // If the "-t" option is supplied, then print the final time in the file
    if(argc==3) printf("%g\n",ltime);
}
