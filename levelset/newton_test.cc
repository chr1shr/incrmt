#include <cstdio>
#include <cmath>
#include <sys/types.h>
#include <sys/stat.h>

#include "level++.hh"

int main() {
    const int m=101,n=101;
    const double ax=-1,ay=-1,bx=1,by=1;
    const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);

    levelset ls(m,n,ax,ay,dx,dy,-1e10,1e10);
    double x0,y0,x,y,*phi=ls.phi,*vp=ls.vp;
    int i,j,ij;

    for(ij=j=0;j<n;j++) {
        y=ay+dy*j;
        for(i=0;i<m;i++,ij++) {
            x=ax+dx*i;
            phi[ij]=sqrt((x-0.08)*(x-0.08)+y*y)-0.4+0.18*sin(11*x)+0.13*sin(13.5*y);
    //        phi[ij]=sqrt(x*x+(abs(y-dx*0.5)-0.4)*(abs(y-dx*0.5)-0.4))-0.3;
        //    phi[ij]=-sqrt(x*x+y*y)+0.6;
            vp[ij]=0;
        }
    }

    ls.build_band();

    for(ij=j=0;j<n;j++) {
        y=ay+dy*j;
        for(i=0;i<m;i++,ij++) {
            x=ax+dx*i;
            phi[ij]=sqrt((x-0.08)*(x-0.08)+y*y)-0.4+0.18*sin(11*x)+0.13*sin(13.5*y);
    //        phi[ij]=sqrt(x*x+(abs(y-dx*0.5)-0.4)*(abs(y-dx*0.5)-0.4))-0.3;
        //    phi[ij]=-sqrt(x*x+y*y)+0.6;
            vp[ij]=0;
        }
    }

    output("phi",ls);

    for(ij=j=0;j<n;j++) {
        for(i=0;i<m;i++,ij++) {
            if(fabs(phi[ij])<2*(dx+dy)) {
                if(ls.newton(i,ij,x0,y0,x,y)) printf("%d %d %d %g %g %g %g\n",i,j,ij,x0,y0,x,y);
                else fprintf(stderr,"Problem at %d %d\n",i,j);
            }
        }
    }
}
