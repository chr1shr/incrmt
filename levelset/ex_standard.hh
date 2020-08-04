#ifndef LEVELPP_EX_STANDARD_HH
#define LEVELPP_EX_STANDARD_HH

#include "config.hh"

class ex_single {
    public:
        const int c;
        ex_single(double *f_) : c(1), f(f_) {}
        inline void acc(int ij,double *p,double w,double x,double y) {
            *p+=w*f[ij];
            p[1]+=w*x*f[ij];
            p[2]+=w*y*f[ij];
        }
        inline void set(int ij,double *p,double al,double be,double ga) {
            f[ij]=*p*al+p[1]*be+p[2]*ga;
        }
    private:
        double* f;
};

class ex_array {
    public:
        int c;
        ex_array() {}
        ex_array(int c_,...);
        void setup_fields(int c_,...);
        inline void acc(int ij,double *p,double w,double x,double y) {
            for(int i=0;i<c;i++) {
                *(p++)+=w*f[i][ij];
                *(p++)+=w*x*f[i][ij];
                *(p++)+=w*y*f[i][ij];
            }
        }
        inline void set(int ij,double *p,double al,double be,double ga) {
            for(int i=0;i<c;i++,p+=3)
                f[i][ij]=*p*al+p[1]*be+p[2]*ga;
        }
    private:
        double *f[max_extrapolation_fields];
};

#endif
