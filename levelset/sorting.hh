#ifndef LEVELPP_SORTING_HH
#define LEVELPP_SORTING_HH

#include "config.hh"

class ls_quicksort {
    public:
        ls_quicksort(double *&iphi) : phi(iphi) {};
        void sort(int *a,int k);
    private:
        double *&phi;
        void quicksort_internal(int *a,int k);
};

class ls_smoothsort {
    public:
        ls_smoothsort(double *&iphi) : phi(iphi) {};
        void sort(int *a,int n);
    private:
#include "built/leonardo.hh"
        double *&phi;
        void sift(int *a,int l,int r);
        void trinkle(int *a,int p,int l,int r);
        void semitrinkle(int *a,int &p,int &l,int &r);
};

#endif
