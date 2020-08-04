#ifndef LEVELPP_TRAVERSE_HH
#define LEVELPP_TRAVERSE_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "config.hh"
#include "helper.hh"
#include "levelset.hh"

template<class T,class E>
class traverse {
    public:
        traverse(levelset &ls_,E &e_);
        ~traverse() {};
        void extrapolate(int i,int ij);
        void extrapolate_staggered(int i,int ij);
    private:
        levelset &ls;
        E &e;
        T z;
        const int m,n,mn,me,ne;
        const double dx,dy,xsp,ysp;
        unsigned int bps,&bpc,sgs,&sgc;
        int &r1,&mem_p1;
        int *s;
        int *ss;
        unsigned int *bp,*sg;
        double *phi,*vp,*&p1;
        double staggered_phi(int ie,int je,int ij);
        void set_staggered(int ie,int je);
};

struct tra_positive {
    inline bool cm(double a,double b) {return a<b;}
    inline int r(int a) {return a;}
};

struct tra_negative {
    inline bool cm(double a,double b) {return a>b;}
    inline int r(int a) {return 7-a;}
};

template<class T,class E>
traverse<T,E>::traverse(levelset &ls_,E &e_) : ls(ls_), e(e_), m(ls_.m),
    n(ls_.n), mn(ls_.mn), me(ls_.me), ne(ls_.ne), dx(ls_.dx), dy(ls_.dy),
    xsp(ls_.xsp), ysp(ls_.ysp), bps(ls_.bpc), bpc(ls_.bpc), sgs(ls_.sgc),
    sgc(ls_.sgc), r1(ls_.r1),
    mem_p1(ls_.mem_p1), s(ls_.s), ss(ls_.ss), bp(ls_.bp), sg(ls_.sg),
    phi(ls_.phi), vp(ls_.vp), p1(ls_.p1) {
    r1=0;
}

template<class T,class E>
void traverse<T,E>::extrapolate_staggered(int i,int ij) {
    int sij[8],sn,j=ij/m,swi,k,l=1;
    double sphi[4],swd;
    unsigned int vls=0;
    bp[ij]=bpc;

    // Assemble a mask describing the 3 by 3 arrangement of grid points
    if(ij>=m) {
        if(i>0) {if(bp[ij-m-1]==bpc) vls|=1;} else vls|=1;
        if(bp[ij-m]==bpc) vls|=2;
        if(i<m-1) {if(bp[ij-m+1]==bpc) vls|=4;} else vls|=4;
    } else vls|=7;
    if(ij<mn-m) {
        if(i>0) {if(bp[ij+m-1]==bpc) vls|=32;} else vls|=32;
        if(bp[ij+m]==bpc) vls|=64;
        if(i<m-1) {if(bp[ij+m+1]==bpc) vls|=128;} else vls|=128;
    } else vls|=224;
    if(i>0) {if(bp[ij-1]==bpc) vls|=8;} else vls|=8;
    if(i<m-1) {if(bp[ij+1]==bpc) vls|=16;} else vls|=16;

    // Test for blocks that have been completed by filling in this grid point
    if((vls&11)==11) {*sij=i;sij[1]=j;*sphi=staggered_phi(i,j,ij);sn=1;} else sn=0;
    if((vls&22)==22) {sij[2*sn]=i+1;sij[2*sn+1]=j;sphi[sn++]=staggered_phi(i+1,j,ij+1);}
    if((vls&104)==104) {sij[2*sn]=i;sij[2*sn+1]=j+1;sphi[sn++]=staggered_phi(i,j+1,ij+m);}
    if((vls&208)==208) {sij[2*sn]=i+1;sij[2*sn+1]=j+1;sphi[sn++]=staggered_phi(i+1,j+1,ij+m+1);}

    // Do a bubble sort on the list of phi values that have been filled in.
    // This list will usually have one or two elements so a bubble sort seems
    // most efficient.
    bool srt;l=0;
    do {
        srt=false;
        for(k=0;k+1<sn;k++) if(z.cm(sphi[k+1],sphi[k])) {
            swi=sij[2*k];sij[2*k]=sij[2*k+2];sij[2*k+2]=swi;
            swi=sij[2*k+1];sij[2*k+1]=sij[2*k+3];sij[2*k+3]=swi;
            swd=sphi[k];sphi[k]=sphi[k+1];sphi[k+1]=swd;
            srt=true;
        }
        l++;
    } while(l<sn&&srt==true);

    // Scan points in the order of the staggered phi values
    for(k=0;k<sn;k++) set_staggered(sij[2*k],sij[2*k+1]);
}

template<class T,class E>
double traverse<T,E>::staggered_phi(int ie,int je,int ij) {
    int di=ie==0?1:(ie==m?-1:0),dj=je==0?1:(je==n?-1:0);
    if(di==0&&dj==0) return 0.25*(phi[ij-m-1]+phi[ij-m]+phi[ij-1]+phi[ij]);

    ij+=di+m*dj;
    return (0.5+dj)*((0.5+di)*phi[ij-m-1]+(0.5-di)*phi[ij-m])
          +(0.5-dj)*((0.5+di)*phi[ij-1]+(0.5-di)*phi[ij]);
}

template<class T,class E>
void traverse<T,E>::set_staggered(int ie,int je) {
    int ije=ie+me*je,di,dj,k,ij=ie+m*je;
    if(sg[ije]==sgc+1) return;

    // Set bounds on which to fit
    const int r=2;
    int ilo=ie<=r?-ie:-r,
        ihi=ie+r>=me?me-1-ie:r,
        jlo=je<=r?-je:-r,
        jhi=je+r>=ne?ne-1-je:r;

    // Initialize accumulators
    int co=0;
    double s=0,sx=0,sy=0,sxx=0,sxy=0,syy=0,w,sphi=staggered_phi(ie,je,ij);
    while(3*e.c>mem_p1) ls.add_memory_p1();
    for(k=0;k<3*e.c;k++) p1[k]=0;

    for(dj=jlo;dj<=jhi;dj++) for(di=ilo;di<=ihi;di++) {
        if(sg[ije+dj*me+di]==sgc+1) {
            w=staggered_phi(ie+di,je+dj,ij+dj*m+di);
            if(z.cm(w,sphi)) {
                w=sphi-w;
                s+=w;
                sx+=w*di;
                sy+=w*dj;
                sxx+=w*di*di;
                sxy+=w*di*dj;
                syy+=w*dj*dj;
                e.acc(ije+dj*me+di,p1,w,di,dj);
            }
        }
    }

    double det=s*(sxy*sxy-sxx*syy)+sx*sx*syy-2*sx*sxy*sy+sy*sy*sxx,al,be,ga;
    if(fabs(det)<1e-10) {
        for(dj=-7;dj<=7;dj++) {
            for(di=-7;di<=7;di++) {
                printf("%s",sg[ije+dj*me+di]==sgc+1?"*":".");
            }
            puts("");
        }
        printf("%d %d co=%d, s=%g, det=%g\n\n",ie,je,co,s,det);
        al=1./s;
        be=0;
        ga=0;
        exit(1);
    } else {
        det=1/det;
        al=(sxy*sxy-sxx*syy)*det;
        be=(sx*syy-sxy*sy)*det;
        ga=(sxx*sy-sxy*sx)*det;
    }
    e.set(ije,p1,al,be,ga);
    sg[ije]=sgc+1;
}

template<class T,class E>
void traverse<T,E>::extrapolate(int i,int ij) {
    int di,dj,j=ij/m,k;

    if(z.r(s[ij])<4) {
        bp[ij]=bpc+1;
        return;
    }

    // Set bounds on which to fit
    double s,sx,sy,sxx,sxy,syy,w,det,al,be,ga;
    int r=3;
    do {
        int ilo=i<=r?-i:-r,
            ihi=i+r>=m?m-1-i:r,
            jlo=j<=r?-j:-r,
            jhi=j+r>=n?n-1-j:r;

        // Initialize accumulators
        s=0;sx=0;sy=0;sxx=0;sxy=0;syy=0;
        while(3*e.c>mem_p1) ls.add_memory_p1();
        for(k=0;k<3*e.c;k++) p1[k]=0;

        for(dj=jlo;dj<=jhi;dj++) for(di=ilo;di<=ihi;di++) {
            if(bp[ij+dj*m+di]==bpc+1) {
                w=phi[ij+dj*m+di];
                if(z.cm(w,phi[ij])) {
                    w=phi[ij]-w;
                    s+=w;
                    sx+=di*w;
                    sy+=dj*w;
                    sxx+=di*di*w;
                    sxy+=di*dj*w;
                    syy+=dj*dj*w;
                    e.acc(ij+dj*m+di,p1,w,di,dj);
                }
            }
        }

        det=s*(sxy*sxy-sxx*syy)+sx*sx*syy-2*sx*sxy*sy+sy*sy*sxx;
        r++;
    } while(fabs(det)<1e-10&&r<=7);

    if(fabs(det)<1e-10) {
        for(dj=-7;dj<=7;dj++) {
            for(di=-7;di<=7;di++) {
                printf("%s",bp[ij+dj*m+di]==bpc+1?"*":".");
            }
            puts("");
        }
        printf("%d %d s=%g, det=%g [reg]\n\n",i,j,s,det);
        al=1./s;
        be=0;
        ga=0;
        exit(1);
    } else {
        det=1/det;
        al=(sxy*sxy-sxx*syy)*det;
        be=(sx*syy-sxy*sy)*det;
        ga=(sxx*sy-sxy*sx)*det;
    }
    e.set(ij,p1,al,be,ga);
    bp[ij]=bpc+1;
}

/** Extrapolates all of the currently registered fields. */
template<class E>
void levelset::extrapolate_fields(E& e) {
    const double s_point=-1.2*sqrt(dx*dx+dy*dy);
    int i,ij,l=u_start;

    // Check that the mask counter has enough space for mn separate entries;
    // otherwise, reset it
    if(bpc+1<1) reset_back_pointers();

    // Mark all the points on the interior of the band from which values can be
    // extrapolated from
    if (l==u_end) return;
    while(phi[u[l]]<s_point) {
        bp[u[l]]=bpc+1;
        l=ul_up(l);
        if (l==u_end) return;
    }

    // Traverse upwards through the band and extrapolate values
    traverse<tra_positive,E> tr1(*this,e);
    while(l!=u_end) {
        ij=u[l];i=ij%m;
        tr1.extrapolate(i,ij);
        l=ul_up(l);
    }
    bpc++;
}

/** Extrapolates all of the currently registered fields. */
template<class E>
void levelset::extrapolate_fields_reverse(E& e) {
    const double s_point=1.2*sqrt(dx*dx+dy*dy);
    int i,ij,l=u_end;

    // Check that the mask counter has enough space for mn separate entries;
    // otherwise, reset it
    if(bpc+1<1) reset_back_pointers();

    // Mark all the points on the interior of the band from which values can be
    // extrapolated from
    if (l==u_start) return;
    while(phi[u[l]]>s_point) {
        l=ul_down(l);
        bp[u[l]]=bpc+1;
        if (l==u_start) return;
    }

    // Traverse upwards through the band and extrapolate values
    traverse<tra_negative,E> tr1(*this,e);
    while(l!=u_start) {
        l=ul_down(l);
        ij=u[l];i=ij%m;
        tr1.extrapolate(i,ij);
    }
    bpc++;
}

/** Extrapolates all of the currently registered fields. */
template<class E>
void levelset::extrapolate_staggered_fields(E& e) {
    const double s_point=-1.2*sqrt(dx*dx+dy*dy);
    int i,ij,l=u_start;

    // Check that the mask counters
    if(bpc+1<1) reset_back_pointers();
    if(sgc+1<1) reset_staggered_back_pointers();

    // Slow setup of staggered mask
    for(ij=0;ij<mne;ij++) if(ss[ij]==0) sg[ij]=sgc+1;

    // Mark all the points on the interior of the band from which values can be
    // extrapolated from
    if (l==u_end) return;
    while(phi[u[l]]<s_point) {
        bp[u[l]]=bpc+1;
        l=ul_up(l);
        if (l==u_end) return;
    }

    // Traverse upwards through the band and extrapolate values
    traverse<tra_positive,E> tr1(*this,e);
    while(l!=u_end) {
        ij=u[l];i=ij%m;
        tr1.extrapolate_staggered(i,ij);
        l=ul_up(l);
    }
    bpc++;
    sgc++;
}

/** Extrapolates all of the currently registered fields. */
template<class E>
void levelset::extrapolate_staggered_fields_reverse(E& e) {
    const double s_point=1.2*sqrt(dx*dx+dy*dy);
    int i,ij,l=u_end;

    // Check that the mask counters
    if(bpc+1<1) reset_back_pointers();
    if(sgc+1<1) reset_staggered_back_pointers();

    // Slow setup of staggered mask
    for(ij=0;ij<mne;ij++) if(ss[ij]==7) sg[ij]=sgc+1;

    // Mark all the points on the interior of the band from which values can be
    // extrapolated from
    if (l==u_start) return;
    while(phi[u[l]]>s_point) {
        l=ul_down(l);
        bp[u[l]]=bpc+1;
        if (l==u_start) return;
    }

    // Traverse upwards through the band and extrapolate values
    traverse<tra_negative,E> tr1(*this,e);
    while(l!=u_start) {
        l=ul_down(l);
        ij=u[l];i=ij%m;
        tr1.extrapolate_staggered(i,ij);
    }
    bpc++;
    sgc++;
}

#endif
