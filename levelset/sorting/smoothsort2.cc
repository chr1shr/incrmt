#include <cstdio>
#include <cstdlib>
#include <cmath>

const int nn=2000;
int a[nn];

inline void sm_sift(int *v,int b,int c,int r) {
    int rm,r2,d,e;
    while(b>=3) {
        rm=r-1;e=b-c-1;r2=rm-e;
        if(v[r2]<v[rm]) {
            if(v[r]<v[rm]) {
                d=v[rm];v[rm]=v[r];v[r]=d;
                b=e;c=c-b-1;
                r=rm;
            } else return;
        } else {
            if(v[r]<v[r2]) {
                d=v[r2];v[r2]=v[r];v[r]=d;
                b=c;c=e;
                r=r2;
            } else return;
        }
    }
}

inline void sm_trinkle(int *v,int p,int b,int c,int r) {
    int r2,r3,d;
    while(p>0) {
        while((p&1)==0) {
            if((p&2)==0) {p>>=2;c=b+c+1;b=b+c+1;}
            else {p>>=1;d=b;b=b+c+1;c=d;break;}
        }
        r3=r-b;
        if(p>1&&v[r3]>v[r]) {
            p--;
            if(b==1) {
                d=v[r3];v[r3]=v[r];v[r]=d;r=r3;
            } else {
                r2=r-b+c;
                if(v[r2]<v[r-1]) {
                    r2=r-1;d=c;c=b-c-1;b=d;p<<=1;
                }
                if(v[r3]<v[r2]) {
                    d=v[r];v[r]=v[r2];v[r2]=d;r=r2;
                    d=c;c=b-c-1;b=d;
                    sm_sift(v,b,c,r);return;
                } else {
                    d=v[r];v[r]=v[r3];v[r3]=d;r=r3;
                }
            }
        } else break;
    }
    sm_sift(v,b,c,r);
}

inline void sm_semitrinkle(int *v,int &p,int &b,int &c,int &r) {
    int r1=r-c;
    if (v[r1]>v[r]) {
        int d=v[r];v[r]=v[r1];v[r1]=d;sm_trinkle(v,p,b,c,r1);
    }
}

void smoothsort(int *v,int n) {
    int q,r=0,p=1,b=1,c=1,d;
    for(q=1;q<n;q++) {
        if((p&7)==3) {
            sm_sift(v,b,c,r);
            p=(p+1)>>2;c=b+c+1;b=b+c+1;
        } else {
            if(q+c<n) sm_sift(v,b,c,r);
            else sm_trinkle(v,p,b,c,r);
            do{
                if(b>3) {
                    p<<=2;b=b-c-1;c=c-b-1;
                } else {
                    p<<=1;d=c;c=b-c-1;b=d;break;
                }
            } while(b>1);
            p++;
        }
        r++;
    }
    sm_trinkle(v,p,b,c,r);
    while(q>1) {
        q--;p--;
        if(b==1) {
            r--;
            while((p&1)==0) {
                if((p&2)==0) {
                    p>>=2;c=b+c+1;b=b+c+1;
                } else {
                    p>>=1;d=b;b=b+c+1;c=d;break;
                }
            }
        } else {
            r=r-b+c;

            if(p>0) sm_semitrinkle(v,p,b,c,r);

            b=b-c-1;
            p=(p<<1)+1;

            r=r+b;sm_semitrinkle(v,p,c,b,r);
            c=c-b-1;
            p=(p<<1)+1;
        }
    }
}

int main() {
    int o,p;
    for (p=0;p<1;p++) {
        for(o=0;o<nn;o++) a[o]=rand()%1000;
        smoothsort(a,nn);
    }
    for(o=0;o<nn;o++) printf("%d %d\n",o,a[o]);
}
