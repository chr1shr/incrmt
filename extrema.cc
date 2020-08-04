#include <cstdio>
#include <cstring>
#include <sys/stat.h>

#include "common.hh"
#include "gp_matrix.hh"

// Maximum safe limit on the number of simulation frames to search for
const int frames_max=65536;

int main(int argc,char **argv) {

    // Check for the correct number of command-line arguments
    if(argc!=3) {
        fputs("Usage: ./extrema <output_dir> <suffix>\n",stderr);
        return 1;
    }

    // Loop over the available frames
    char buf[256];
    struct stat fe;
    for(int i=0;i<=frames_max;i++) {

        // Check that the file exists
        sprintf(buf,"%s/%s.%d",argv[1],argv[2],i);
        if(stat(buf,&fe)!=0) return 0;

        // Read in the file and calculate a variety of norms of the field values
        gp_matrix gp(buf);
        float l1,l2,l4=gp.lp_norm(4.),l8=gp.lp_norm(8.),lmin,lmax;
        gp.field_range(lmin,lmax);
        gp.l1_l2_norm(l1,l2);
        printf("%d %g %g %g %g %g\n",i,l1,l2,l4,l8,lmax>-lmin?lmax:-lmin);
    }
}
