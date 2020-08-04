#ifndef MAT_HH
#define MAT_HH

#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "vec.hh"

struct sym_mat;

struct mat {
    public:
        double a,b,c,d;
        mat() {}
        mat(double a_) : a(a_),b(0),c(0),d(a_) {}
        mat(double a_,double b_,double c_,double d_) : a(a_),b(b_),c(c_), d(d_) {}
        inline mat operator+ (mat p) {return mat(a+p.a,b+p.b,c+p.c,d+p.d);}
        inline mat operator- (mat p) {return mat(a-p.a,b-p.b,c-p.c,d-p.d);}
        inline mat operator* (double e) {return mat(a*e,b*e,c*e,d*e);}
        inline mat operator*(mat e) {return mat(a*e.a+b*e.c,a*e.b+b*e.d,c*e.a+d*e.c,c*e.b+d*e.d);}
        inline vec operator* (vec e) {return vec(a*e.x+b*e.y,c*e.x+d*e.y);}
        inline mat operator/ (double e) {
            double ei=1/e;
            return mat(a*ei,b*ei,c*ei,d*ei);
        }
        inline void operator+= (mat p) {
            a+=p.a;b+=p.b;
            c+=p.c;d+=p.d;
        }
        inline void operator-= (mat p) {
            a-=p.a;b-=p.b;
            c-=p.c;d-=p.d;
        }
        inline void operator*= (double e) {
            a*=e;b*=e;c*=e;d*=e;
        }
        inline void operator*=(mat e) {
            double pa=a,pc=c;
            a=a*e.a+b*e.c;
            b=pa*e.b+b*e.d;
            c=c*e.a+d*e.c;
            d=pc*e.b+d*e.d;
        }
        inline void operator/= (double e) {
            a/=e;b/=e;c/=e;d/=e;
        }
        inline void set(double a_,double b_,double c_,double d_) {
            a=a_;b=b_;c=c_;d=d_;
        }
        inline double devmod() {
            return sqrt(0.5*(a-d)*(a-d)+b*b+c*c);
        }
        inline double trace() {return a+d;}
        inline double det() {return a*d-b*c;}
        inline double mod_sq() {return a*a+b*b+c*c+d*d;}
        inline mat transpose() {return mat(a,c,b,d);}
        inline mat inverse() {
            double idet=1.0/det();
            return mat(idet*d,-idet*b,-idet*c,idet*a);
        }
        inline mat inv_transpose() {
            double idet=1.0/det();
            return mat(idet*d,-idet*c,-idet*b,idet*a);
        }
        inline void remove_trace() {
            double rem=0.5*(a+d);
            a-=rem;d-=rem;
        }
        inline void remove_trace3() {
            double rem=(a+d+1)/3.;
            a-=rem;d-=rem;
        }
        inline sym_mat AAT();
        inline mat operator*(sym_mat e);
        inline sym_mat ATA();
        inline sym_mat ATDA(double l1,double l2);
        inline sym_mat ATSA(sym_mat m);
        inline sym_mat sym();
        inline double frob_norm() {return a*a+b*b+c*c+d*d;};
        inline void print_mat() {printf(" [%g %g %g %g]",a,b,c,d);}
        void sym_eigenvectors(double &l1,double &l2,mat &Lam);
};

inline mat operator*(const double e,mat f) {
    return mat(e*f.a,e*f.b,e*f.c,e*f.d);
}

inline mat operator/(double e,mat f) {
    double idet=e/(f.a*f.d-f.b*f.c);
    return mat(f.d*idet,-f.b*idet,-f.c*idet,f.a*idet);
}

inline mat operator-(mat f) {
    return mat(-f.a,-f.b,-f.c,-f.d);
}

struct sym_mat {
    public:
        double a,b,d;
        sym_mat() {};
        sym_mat(double a_) : a(a_),b(0),d(a_) {};
        sym_mat(double a_,double b_,double d_) : a(a_),b(b_),d(d_) {};
        inline sym_mat operator+ (sym_mat p) {return sym_mat(a+p.a,b+p.b,d+p.d);}
        inline sym_mat operator- (sym_mat p) {return sym_mat(a-p.a,b-p.b,d-p.d);}
        inline sym_mat operator* (double e) {return sym_mat(a*e,b*e,d*e);}
        inline sym_mat operator*(sym_mat e) {return sym_mat(a*e.a+b*e.b,a*e.b+b*e.d,b*e.b+d*e.d);}
        inline mat operator*(mat e) {return mat(a*e.a+b*e.c,a*e.b+b*e.d,b*e.a+d*e.c,b*e.b+d*e.d);}
        inline vec operator*(vec e) {return vec(a*e.x+b*e.y,b*e.x+d*e.y);}
        inline sym_mat operator/ (double e) {
            double ei=1/e;
            return sym_mat(a*ei,b*ei,d*ei);
        }
        inline void operator+= (sym_mat p) {a+=p.a;b+=p.b;d+=p.d;}
        inline void operator-= (sym_mat p) {a-=p.a;b-=p.b;d-=p.d;}
        inline void operator*= (double e) {a*=e;b*=e;d*=e;}
        inline void operator/= (double e) {a/=e;b/=e;d/=e;}
        inline void set(double a_,double b_,double d_) {a=a_;b=b_;d=d_;}
        inline double devmod() {
            return sqrt(0.5*(a-d)*(a-d)+2*b*b);
        }
        inline double trace() {return a+d;}
        inline double invariant2() {return a*d-b*b;}
        inline double det() {return a*d-b*b;}
        inline double mod_sq() {return a*a+2*b*b+d*d;}
        inline sym_mat transpose() {return sym_mat(a,b,d);}
        inline sym_mat inverse() {
            double idet=1.0/det();
            return sym_mat(idet*d,-idet*b,idet*a);
        }
        inline sym_mat inv_transpose() {
            return inverse();
        }
        inline void remove_trace() {
            double rem=0.5*(a+d);
            a-=rem;d-=rem;
        }
        inline void remove_trace3() {
            double rem=(1./3.)*(a+d+1);
            a-=rem;d-=rem;
        }
        inline double frob_norm() {return a*a+2*b*b+d*d;};
        inline sym_mat MTA_plus_AM(mat m);
        inline void print_sym_mat() {printf(" [%g %g %g %g]",a,b,b,d);}
        void eigenvectors(double &l1,double &l2,mat &Lam);

        double& operator()(int i) {return i==0?a:(i==1?b:d);};
        const double& operator()(int i) const {return i==0?a:(i==1?b:d);};
};

inline sym_mat mat::sym() {
    return sym_mat(a,0.5*(b+c),d);
}

inline sym_mat sym_mat::MTA_plus_AM(mat m) {
    return sym_mat(2*(m.a*a+m.c*b),m.trace()*b+m.b*a+m.c*d,2*(m.b*b+m.d*d));
}

inline sym_mat mat::AAT() {
    return sym_mat(a*a+b*b,a*c+b*d,c*c+d*d);
}

inline sym_mat mat::ATA() {
    return sym_mat(a*a+c*c,a*b+c*d,b*b+d*d);
}

inline sym_mat mat::ATDA(double l1,double l2) {
    return sym_mat(a*a*l1+c*c*l2,a*b*l1+c*d*l2,b*b*l1+d*d*l2);
}

inline sym_mat mat::ATSA(sym_mat s) {
    double e=a*s.a+c*s.b,f=a*s.b+c*s.d,
           g=b*s.a+d*s.b,h=b*s.b+d*s.d;
    return sym_mat(a*e+c*f,b*e+d*f,b*g+d*h);
}

inline mat mat::operator*(sym_mat e) {
        return mat(a*e.a+b*e.b,a*e.b+b*e.d,c*e.a+d*e.b,c*e.b+d*e.d);
}

inline sym_mat operator*(const double e,sym_mat f) {
    return sym_mat(e*f.a,e*f.b,e*f.d);
}

inline sym_mat operator/(double e,sym_mat f) {
    double idet=e/(f.a*f.d-f.b*f.b);
    return sym_mat(f.d*idet,-f.b*idet,f.a*idet);
}

inline sym_mat operator-(sym_mat f) {
    return sym_mat(-f.a,-f.b,-f.d);
}

#endif
