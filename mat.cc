#include "mat.hh"

/** Calculates the eigenvectors of the matrix, assuming it is symmetric so that
 * b=c.
 * \param[out] l1 the first eigenvalue (l1>=l2).
 * \param[out] l2 the second eigenvalue.
 * \param[out] Lam the corresponding eigenvectors. */
void mat::sym_eigenvectors(double &l1,double &l2,mat &Lam) {
    if(fabs(b)<1e-12*(fabs(a)+fabs(d))) {
        if(a>d) {
            l1=a;l2=d;
            Lam.set(1,0,0,1);
            return;
        } else {
            l1=d;l2=a;
            Lam.set(0,1,1,0);
            return;
        }
    }
    double det=sqrt((a-d)*(a-d)+4*b*b);
    l1=0.5*(a+d+det);
    l2=0.5*(a+d-det);
    double r=b,s=a-l1;
    double nor1=1/sqrt(r*r+s*s);
    Lam.set(r*nor1,-s*nor1,s*nor1,r*nor1);
}

/** Calculates the eigenvectors of the matrix.
 * \param[out] l1 the first eigenvalue (l1>=l2).
 * \param[out] l2 the second eigenvalue.
 * \param[out] Lam the corresponding eigenvectors. */
void sym_mat::eigenvectors(double &l1,double &l2,mat &Lam) {
    if(fabs(b)<1e-12*(fabs(a)+fabs(d))) {
        if(a>d) {
            l1=a;l2=d;
            Lam.set(1,0,0,1);
            return;
        } else {
            l1=d;l2=a;
            Lam.set(0,1,1,0);
            return;
        }
    }
    double det=sqrt((a-d)*(a-d)+4*b*b);
    l1=0.5*(a+d+det);
    l2=0.5*(a+d-det);
    double r=b,s=a-l1;
    double nor1=1/sqrt(r*r+s*s);
    Lam.set(r*nor1,-s*nor1,s*nor1,r*nor1);
}
