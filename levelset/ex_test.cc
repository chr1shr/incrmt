#include <cstdio>
#include <cmath>
#include <ctime>
#include <sys/types.h>
#include <sys/stat.h>

#include "level++.hh"

clock_t ci;
const double pi=3.1415926535897932384626433832795;

inline double fphi(double x,double y) {return sqrt(x*x+y*y)-2;}
inline double fu(double x,double y) {return x+y;}

int main(int argc,char **argv) {
    if(argc!=2&&argc!=3) {
        fputs("Syntax: ./ex_test <grid_size> <trials>\n",stderr);
        return 1;
    }
    clock_t c1,c2,c3;
    double sctime=0,bbtime=0,extime=0;
    int n=atof(argv[1]),nn=n*n;
    int trials(argc==3?atoi(argv[2]):1);
    double a=-(1-1.0/n)*pi,b=-a;
    double d=(b-a)/(n-1);

    levelset ls(n,n,a,a,d,d,-5,5);
//    levelset ls(n,n,a,a,d,d,-5*d*sqrt(2),2.5*d*sqrt(2));trials*=4;

    double x,y,temp,*phi=ls.phi;
    double *u=new double[nn];
    double l1=0,l2=0,linf=0,lsinf=0;
    double lp1=0,lp2=0,lpinf=0,lpsinf=0;
    int i,j,ij;
    ex_single e(u);

    for(int tr=0;tr<trials;tr++) {
        for(ij=j=0;j<n;j++) {
            y=a+d*j;
            for(i=0;i<n;i++,ij++) {
                x=a+d*i;
                phi[ij]=sqrt(x*x+y*y)-2;
                //u[ij]=x+0.5*y;
                u[ij]=(x+0.5*y)*(x+0.5*y);
                if(phi[ij]>0) u[ij]=0;
            }
        }

        output("ua",0,u,n,n,a,a,d,d);
        output("phia",0,ls);
        c1=clock();
        ls.build_band();
        c2=clock();
        sctime+=double(ci-c1)/CLOCKS_PER_SEC;
        bbtime+=double(c2-ci)/CLOCKS_PER_SEC;

    //    ls.check_sorted();

    /*    for(ij=j=0;j<n;j++) {
            y=a+d*j;
            for(i=0;i<n;i++,ij++) {
                x=a+d*i;
                phi[ij]=sqrt(x*x+y*y)-2;
            }
        }*/

        ls.extrapolate_fields(e);
        c3=clock();
        extime+=double(c3-c2)/CLOCKS_PER_SEC;
    }

    output("ue",0,u,n,n,a,a,d,d);
    output("phi",0,ls);
    return 0;

    for(ij=j=0;j<n;j++) {
        y=a+d*j;
        for(i=0;i<n;i++,ij++) {
            x=a+d*i;
            temp=d*d*(i==0||i==n-1?0.5:1)*(j==0||j==n-1?0.5:1);
            u[ij]-=(x+0.5*y)*(x+0.5*y);
            phi[ij]-=sqrt(x*x+y*y)-2;
            l1+=temp*fabs(u[ij]);
            l2+=temp*u[ij]*u[ij];
            if(fabs(u[ij])>linf) linf=fabs(u[ij]);
            if(ls.s[ij]==4&&fabs(u[ij])>lsinf) lsinf=fabs(u[ij]);
            lp1+=temp*fabs(phi[ij]);
            lp2+=temp*phi[ij]*phi[ij];
            if(fabs(phi[ij])>lpinf) lpinf=fabs(phi[ij]);
            if(ls.s[ij]==4&&fabs(phi[ij])>lpsinf) lpsinf=fabs(phi[ij]);
        }
    }

    output("phid",0,ls);
    output("ud",0,u,n,n,a,a,d,d);
    printf("%d %g %g %g %g %g %g %g %g %g %g %g\n",n,l1,sqrt(l2),linf,lsinf,lp1,sqrt(lp2),lpinf,sctime/trials,bbtime/trials, extime/trials,lpsinf);

    delete [] u;
}
