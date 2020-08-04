#ifndef LEVELPP_HEAP_HH
#define LEVELPP_HEAP_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "config.hh"

class levelset;

template<class d_option>
class heap {
    public:
        levelset &ls;
        const int m,n,mn;
        const double dx,dy,xsp,ysp,dx2,dy2,dxbydy,dx2bydy2;
        unsigned int bpc;
        int &mem;
        int *&hp;
        int *s;
        unsigned int *bp;
        double *phi;
        int w;
        d_option z;
        heap(levelset &ils,int &imem,int *&ihp);
        void march(double limit);
        void full_march();
        void stack_test(int i,int ij);
        void stack_test2(int i,int ij);
        void update(int i,int ij);
        void scan(int i,int ij);
        void trickle(int c,int ij,double tphi);
        void introduce(int ij,double tphi);
        void add(int ij,double tphi);
        void calc_tentative();
        void build_heap();
        void check_heap();
        double compute_eikonal(int i,int ij,int &stat);
    private:
        void add_memory();
        void check_gridpoint(int i,int ij);
        void check_gridpoint2(int i,int ij);
        bool look_l(int i,int ij,double &e);
        bool look_r(int i,int ij,double &e);
        bool look_d(int ij,double &e);
        bool look_u(int ij,double &e);
        double threeway1(bool h1,bool h2,bool v,double &a,double &c,double &b);
        double threeway2(bool h,bool v1,bool v2,double &b,double &a,double &c);
        double fourway(double &a,double &b,double &c,double &d);
        double corner(bool h,bool v,double &a,double &b);
        double corner_bail(bool h,bool v,double &a,double &b);
        double calc(int i,int ij);
        void list(int i,int ij);
};

class dir_positive {
    public:
        levelset &ls;
        const int o,p,q;
        dir_positive(levelset &ils,int *&s);
        inline bool cm(double a,double b) {return a<b;};
        inline bool in(int ij) {return ss[ij]<6;};
        inline double ct(double a,double b) {return a+b;};
        inline void add(int ij);
    private:
        int *&ss;
};

class dir_negative {
    public:
        levelset &ls;
        const int o,p,q;
        dir_negative(levelset &ils,int *&s);
        inline bool cm(double a,double b) {return a>b;};
        inline bool in(int ij) {return ss[ij]>1;};
        inline double ct(double a,double b) {return a-b;};
        inline void add(int ij);
    private:
        int *&ss;
};

#endif
