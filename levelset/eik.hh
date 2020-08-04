#ifndef LEVELPP_EIK_HH
#define LEVELPP_EIK_HH

class eik_l2norm {
    private:
        const int m,n,mn;
        const double nfac;
        double b;
    public:
        eik_l2norm(levelset &ls) : m(ls.m), n(ls.n), mn(ls.mn),
                    nfac((ls.bx-ls.ax)*(ls.by-ls.ay)/(4*(m-1)*(n-1))),b(0) {};
        inline void contrib(int i,int ij,double c) {
            int d=i==0||i==m-1?1:2;
            if(ij>=m&&ij<mn-m) d*=2;
            b+=(c-1)*(c-1)*d;
        }
        inline double output() {return sqrt(b*nfac);}
};

class eik_l1norm {
    private:
        const int m,n,mn;
        const double nfac;
        double b;
    public:
        eik_l1norm(levelset &ls) : m(ls.m), n(ls.n), mn(ls.mn),
                    nfac((ls.bx-ls.ax)*(ls.by-ls.ay)/(4*(m-1)*(n-1))),b(0) {};
        inline void contrib(int i,int ij,double c) {
            int d=i==0||i==m-1?1:2;
            if(ij>=m&&ij<mn-m) d*=2;
            b+=(c<1?1-c:c-1)*d;
        }
        inline double output() {return b*nfac;}
};

class eik_max {
    private:
        double b;
    public:
        eik_max(levelset &ls) : b(0) {};
        inline void contrib(int i,int ij,double c) {c=c<1?1-c:c-1;if(c>b) b=c;}
        inline double output() {return b;}
};

#endif
