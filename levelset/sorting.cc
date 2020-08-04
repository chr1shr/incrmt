#include "sorting.hh"

void ls_quicksort::sort(int *a,int k) {
    int a1;
    if (k<3) {
        if (k==2&&phi[a[0]]>phi[a[1]]) {a1=a[0];a[0]=a[1];a[1]=a1;}
    } else {
        if (k>3) quicksort_internal(a,k);
        else {
            a1=a[0];int a2=a[1],a3=a[2];
            double p1=phi[a1],p2=phi[a2],p3=phi[a3];
            if(p1<p2) {
                if(p3<p2) {
                    if(p1<p3) a[1]=a3;
                    else {a[0]=a3;a[1]=a1;}
                    a[2]=a2;
                }
            } else {
                if(p3>p1) {
                    a[0]=a2;a[1]=a1;a[2]=a3;
                } else {
                    if(p2>p3) a[0]=a3;
                    else {a[0]=a2;a[1]=a3;}
                    a[2]=a1;
                }
            }
        }
    }
}

void ls_quicksort::quicksort_internal(int *a,int k) {
    k--;int pivot=k>>1;
    int h,x=a[pivot],i=0,j=k-1;
    double px=phi[x];a[pivot]=a[k];

    while(phi[a[i]]<px&&i<j-1) i++;
    while(phi[a[j]]>px&&i+1<j) j--;

    while(i+1<j) {
        h=a[i];a[i++]=a[j];a[j--]=h;
        while(phi[a[i]]<px&&i<j-1) i++;
        while(phi[a[j]]>px&&i+1<j) j--;
    }

    if(i==j) {
        if(phi[a[j]]<px) j++;
        a[k]=a[j];a[j]=x;
    } else {
        if(phi[a[i]]<px) {
            if(phi[a[j]]<px) j++;
            a[k]=a[j];a[j]=x;
        } else {
            if(phi[a[j]]<px) {
                a[k]=a[i];a[i]=a[j];a[j]=x;
            } else {
                a[k]=a[i];a[i]=x;j=i;
            }
        }
    }

    sort(a,j);
    sort(a+(j+1),k-j);
}

inline void ls_smoothsort::sift(int *a,int l,int r) {
    int rm,r2,d;
    while(l>=2) {
        rm=r-1;r2=rm-smc[l-1];
        if(phi[a[r2]]<phi[a[rm]]) {
            if(phi[a[r]]<phi[a[rm]]) {
                d=a[rm];a[rm]=a[r];a[r]=d;
                l-=2;r=rm;
            } else return;
        } else {
            if(phi[a[r]]<phi[a[r2]]) {
                d=a[r2];a[r2]=a[r];a[r]=d;
                l--;r=r2;
            } else return;
        }
    }
}

void ls_smoothsort::trinkle(int *a,int p,int l,int r) {
    int r2,r3,d;
    while(p>0) {
        while((p&1)==0) {
            if((p&2)==0) {p>>=2;l+=2;}
            else {p>>=1;l++;break;}
        }
        r3=r-smb[l];
        if(p>1&&phi[a[r3]]>phi[a[r]]) {
            p--;
            if(l<=1) {
                d=a[r3];a[r3]=a[r];a[r]=d;r=r3;
            } else {
                r2=r-smd[l];
                if(phi[a[r2]]<phi[a[r-1]]) {
                    r2=r-1;l--;p<<=1;
                }
                if(phi[a[r3]]<phi[a[r2]]) {
                    d=a[r];a[r]=a[r2];a[r2]=d;r=r2;
                    l--;sift(a,l,r);return;
                } else {
                    d=a[r];a[r]=a[r3];a[r3]=d;r=r3;
                }
            }
        } else break;
    }
    sift(a,l,r);
}

void ls_smoothsort::semitrinkle(int *a,int &p,int &l,int &r) {
    int r1=r-smc[l];
    if (phi[a[r1]]>phi[a[r]]) {
        int d=a[r];a[r]=a[r1];a[r1]=d;trinkle(a,p,l,r1);
    }
}

void ls_smoothsort::sort(int *a,int n) {
    int q,r=0,p=1,l=1;
    for(q=1;q<n;q++) {
        if((p&7)==3) {
            sift(a,l,r);
            p=(p+1)>>2;l+=2;
        } else {
            if(q+smc[l]<n) sift(a,l,r);else trinkle(a,p,l,r);
            if (l==1) {p<<=1;l=0;} else {p<<=l-1;l=1;}
            p++;
        }
        r++;
    }
    trinkle(a,p,l,r);
    while(q>1) {
        q--;p--;
        if(l<=1) {
            r--;
            while((p&1)==0) {
                if((p&2)==0) {p>>=2;l+=2;}
                else {p>>=1;l+=1;break;}
            }
        } else {
            r=r-smd[l];

            if(p>0) semitrinkle(a,p,l,r);
            l--;p=(p<<1)+1;

            r=r+smc[l];semitrinkle(a,p,l,r);
            l--;p=(p<<1)+1;
        }
    }
}

#include "built/leonardo.cc"
