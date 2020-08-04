#include <cstdio>
#include <cstdlib>
#include <cmath>

const int n=2000;
int a[n];

inline void printit() {
    for(int i=0;i<n-1;i++) printf("%d,",a[i]);
    printf("%d\n",a[n-1]);
}

void sm_sift(int b,int c,int r) {
    int rm,r2,d,e;
    while(b>=3) {
        rm=r-1;e=b-c-1;r2=rm-e;
        if(a[r2]<a[rm]) {
            if(a[r]<a[rm]) {
                d=a[rm];a[rm]=a[r];a[r]=d;
                b=e;c=c-b-1;
                r=rm;
            } else return;
        } else {
            if(a[r]<a[r2]) {
                d=a[r2];a[r2]=a[r];a[r]=d;
                b=c;c=e;
                r=r2;
            } else return;
        }
    }
}

inline void sm_trinkle(int p,int b,int c,int r) {
    int r2,r3,d;
    while(p>0) {
        while((p&1)==0) {
            if((p&2)==0) {p>>=2;c=b+c+1;b=b+c+1;}
            else {p>>=1;d=b;b=b+c+1;c=d;break;}
        }
        r3=r-b;
        if(p>1&&a[r3]>a[r]) {
            p--;
            if(b==1) {
                d=a[r3];a[r3]=a[r];a[r]=d;r=r3;
            } else {
                r2=r-b+c;
                if(a[r2]<a[r-1]) {
                    r2=r-1;d=c;c=b-c-1;b=d;p<<=1;
                }
                if(a[r3]<a[r2]) {
                    d=a[r];a[r]=a[r2];a[r2]=d;r=r2;
                    d=c;c=b-c-1;b=d;
                    sm_sift(b,c,r);return;
                } else {
                    d=a[r];a[r]=a[r3];a[r3]=d;r=r3;
                }
            }
        } else break;
    }
    sm_sift(b,c,r);
}

inline void sm_semitrinkle(int &p,int &b,int &c,int &r) {
    int r1=r-c;
    if (a[r1]>a[r]) {
        int d=a[r];a[r]=a[r1];a[r1]=d;sm_trinkle(p,b,c,r1);
    }
}

void smoothsort() {
    int q,r=0,p=1,b=1,c=1,d;
    for(q=1;q<n;q++) {
        if((p&7)==3) {
            sm_sift(b,c,r);
            p=(p+1)>>2;c=b+c+1;b=b+c+1;
        } else {
            if(q+c<n) sm_sift(b,c,r);
            else sm_trinkle(p,b,c,r);
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
    sm_trinkle(p,b,c,r);
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

            if(p>0) sm_semitrinkle(p,b,c,r);

            b=b-c-1;
            p=(p<<1)+1;

            r=r+b;sm_semitrinkle(p,c,b,r);
            c=c-b-1;
            p=(p<<1)+1;
        }
    }
}

int main() {
    int o,p;
    for (p=0;p<1;p++) {
        for(o=0;o<n;o++) a[o]=rand()%1000;
        smoothsort();
    }
    for(o=0;o<n;o++) printf("%d %d\n",o,a[o]);
}
