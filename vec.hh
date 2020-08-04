#ifndef VEC_HH
#define VEC_HH

#include <cmath>

struct vec {
    public:
        double x,y;
        vec() {}
        vec(double x_) : x(x_),y(x_) {}
        vec(double x_,double y_) : x(x_),y(y_) {}
        inline vec operator+ (vec p) {
            return vec(x+p.x,y+p.y);
        }
        inline vec operator- (vec p) {
            return vec(x-p.x,y-p.y);
        }
        inline vec operator* (double a) {
            return vec(x*a,y*a);
        }
        inline vec operator/ (double a) {
            return vec(x/a,y/a);
        }
        inline void operator+= (vec p) {
            x+=p.x;y+=p.y;
        }
        inline void operator-= (vec p) {
            x-=p.x;y-=p.y;
        }
        inline void operator*= (double a) {
            x*=a;y*=a;
        }
        inline void operator/= (double a) {
            x/=a;y/=a;
        }
        inline operator double() {
            return double(x*x+y*y);
        }
        inline double magnitude() {
            return sqrt(x*x+y*y);
        }
        inline void normalize() {
            double mag=x*x+y*y;
            if(mag>0) {
                mag=1./sqrt(mag);
                x*=mag;y*=mag;
            }
        }
};

inline vec operator*(const double e,vec f) {
    return vec(e*f.x,e*f.y);
}

inline vec operator-(vec f) {
    return vec(-f.x,-f.y);
}

inline double mod_sq(vec a) {return a.x*a.x+a.y*a.y;}

#endif
