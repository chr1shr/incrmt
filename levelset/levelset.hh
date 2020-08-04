#ifndef LEVELPP_LEVELSET_HH
#define LEVELPP_LEVELSET_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "config.hh"
#include "helper.hh"
#include "heap.hh"
#include "sorting.hh"

class levelset {
    public:
        levelset(int m_,int n_,double ax_,double ay_,double dx_,double dy_,double phi_down_=-1e30,double phi_up_=1e30);
        ~levelset();
        const int m,n,mn,me,ne,mne;
        const double ax,ay,bx,by,dx,dy;
        const double xsp,ysp;
        const double dx2,dy2;
        const double phi_down,phi_up;
        int u_mem,u_start,u_mid,u_end;
        unsigned int bpc;
        unsigned int sgc;
        int num_n,num_p;
        int q1,r1;
        int mem_t1,mem_p1;
        int mem_up,mem_down;
        int *s;
        int *ss;
        unsigned int *bp;
        unsigned int *sg;
        double *phi;
        double *vp;
        int *u;
        int *h_up;
        int *h_down;
        int *t1;
        double *p1;
        ls_quicksort sort;
        bool failsafe;
        double bicubic(double x,double y,double& phix,double& phiy);
        template<class T>
        T bicubic(T *fe,double x,double y);
        double bilinear(double *fe,double x,double y);
        bool newton(int i,int ij,double &x0,double &y0,double &px,double &py);
        bool newton_failsafe(int i,int ij,double &x0,double &y0,double &px,double &py);
        void add_positive(int ij);
        void add_negative(int ij);
        void build_band();
        template<class m_class>
        void move(m_class* mc,double dt);
        template<class m_class>
        void move_and_smooth(m_class* mc,double dt,double smth);
        void clean();
        template<class E>
        void extrapolate_fields(E &e);
        template<class E>
        void extrapolate_fields_reverse(E &e);
        template<class E>
        void extrapolate_staggered_fields(E &e);
        template<class E>
        void extrapolate_staggered_fields_reverse(E &e);
        double centered_deriv(int i,int ij);
        void check_sorted();
        void check_counts();
        void check_divided();
        void count_s_array();
        void output_band(const char* filename);
        void add_memory_p1();
    private:
        inline int ul_up(int i) {return i==u_mem-1?0:i+1;}
        inline int ul_down(int i) {return i==0?u_mem-1:i-1;}
        void sort_boundary(int uv,int u_min,int u_max);
        void sort_boundary_wrap();
        int search_down_simple(int uv,int u_min,double tphi);
        int search_down(int iu,double iphi);
        int search_up(int iu,double iphi);
        bool check_up(int i,int ij);
        bool check_down(int i,int ij);
        bool check_neg(int i,int ij);
        bool check_pos(int i,int ij);
        void change_to_neg(int i,int ij);
        void change_to_pos(int i,int ij);
        int bisect(int lu,int uu,double iphi);
        void reset_back_pointers();
        void reset_staggered_back_pointers();
        void check_two_points(int i,int j);
        void add_u_memory_init();
        void add_memory_t1();
        void post_move_update(double dt);
        double laplacian(int i,int ij);
};

inline void output(const char* filename,levelset &ls) {
    output(filename,ls.phi,ls.m,ls.n,ls.ax,ls.ay,ls.dx,ls.dy);
}

inline void output(const char* filename,int index,levelset &ls) {
    char buffer[256];
    sprintf(buffer,"output/%s.%d",filename,index);
    output(buffer,ls);
}

template<class m_class>
void levelset::move(m_class* mc,double dt) {
    int l,i,ij;
    double x0,y0,x,y;

    q1=0;

    // Construct velocity extensions for positive points on the interface,
    // and store the position of the first non-interface point in l_p
    l=u_mid;
    while(l!=u_end) {
        ij=u[l];i=ij%m;
        if(s[ij]==4) {
            if(q1==mem_t1) add_memory_t1();
            vp[ij]=newton(i,ij,x0,y0,x,y)?mc->velocity(x,y):mc->velocity(x0,y0);
            t1[q1++]=ij;
        } else {
            s[ij]=7;
            if(ij>=m&&s[ij-m]==6) s[ij-m]=7;
            if(i>0&&s[ij-1]==6) s[ij-1]=7;
            if(i<m-1&&s[ij+1]==6) s[ij+1]=7;
            if(ij<mn-m&&s[ij+m]==6) s[ij+m]=7;
        }
        l=ul_up(l);
    }

    // Construct velocity extensions for negative points on the interface,
    // and store the position of the first non-interface point in l_n
    l=u_mid;
    while(l!=u_start) {
        l=ul_down(l);
        ij=u[l];i=ij%m;
        if(s[ij]==3) {
            if(q1==mem_t1) add_memory_t1();
            vp[ij]=newton(i,ij,x0,y0,x,y)?mc->velocity(x,y):mc->velocity(x0,y0);
            t1[q1++]=ij;
        } else {
            s[ij]=0;
            if(ij>=m&&s[ij-m]==1) s[ij-m]=0;
            if(i>0&&s[ij-1]==1) s[ij-1]=0;
            if(i<m-1&&s[ij+1]==1) s[ij+1]=0;
            if(ij<mn-m&&s[ij+m]==1) s[ij+m]=0;
        }
    }
    post_move_update(dt);
}

template<class m_class>
void levelset::move_and_smooth(m_class* mc,double dt,double smth) {
    int l,i,ij;
    double x0,y0,x,y;

    q1=0;

    // Construct velocity extensions for positive points on the interface,
    // and store the position of the first non-interface point in l_p
    l=u_mid;
    while(l!=u_end) {
        ij=u[l];i=ij%m;
        if(s[ij]==4) {
            if(q1==mem_t1) add_memory_t1();
            vp[ij]=newton(i,ij,x0,y0,x,y)?mc->velocity(x,y):mc->velocity(x0,y0);
            vp[ij]+=smth*laplacian(i,ij);
            t1[q1++]=ij;
        } else {
            s[ij]=7;
            if(ij>=m&&s[ij-m]==6) s[ij-m]=7;
            if(i>0&&s[ij-1]==6) s[ij-1]=7;
            if(i<m-1&&s[ij+1]==6) s[ij+1]=7;
            if(ij<mn-m&&s[ij+m]==6) s[ij+m]=7;
        }
        l=ul_up(l);
    }

    // Construct velocity extensions for negative points on the interface,
    // and store the position of the first non-interface point in l_n
    l=u_mid;
    while(l!=u_start) {
        l=ul_down(l);
        ij=u[l];i=ij%m;
        if(s[ij]==3) {
            if(q1==mem_t1) add_memory_t1();
            vp[ij]=newton(i,ij,x0,y0,x,y)?mc->velocity(x,y):mc->velocity(x0,y0);
            vp[ij]+=smth*laplacian(i,ij);
            t1[q1++]=ij;
        } else {
            s[ij]=0;
            if(ij>=m&&s[ij-m]==1) s[ij-m]=0;
            if(i>0&&s[ij-1]==1) s[ij-1]=0;
            if(i<m-1&&s[ij+1]==1) s[ij+1]=0;
            if(ij<mn-m&&s[ij+m]==1) s[ij+m]=0;
        }
    }
    post_move_update(dt);
}

/** Carries out a bicubic interpolation of any function at a given point.
 * \param[in] fe a pointer to the function to interpolate.
 * \param[in] (x,y) the position to interpolate at.
 * \return the interpolated value. */
template<class T>
T levelset::bicubic(T *fe,double x,double y) {
    int i,j,ij;
    const double dxf=0.5,dyf=0.5,d2f=0.25;
    T plbx,prbx,pltx,prtx;
    T plby,prby,plty,prty;
    T plbxy,prbxy,pltxy,prtxy;
    double me,mf;
    double e=(x-ax)*xsp;
    double f=(y-ay)*ysp;
    i=int(e);j=int(f);
#include "built/bicubic2.cc"
}

#endif
