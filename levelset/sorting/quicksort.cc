#include <cstdio>
#include <cstdlib>
#include <cmath>

const int n=100000;
int a[n];
int b[n];

void quicksort_internal(int lo,int hi);

inline void quicksort(int i,int j) {
    int k=j-i;
    if (k<2) {
        if (i<j&&a[i]>a[j]) {k=a[i];a[i]=a[j];a[j]=k;}
    } else {
        if (k>2) quicksort_internal(i,j);
        else {
            int a1=a[i],a2=a[i+1],a3=a[i+2];
            if(a1<a2) {
                if(a3<a2) {
                    if(a1<a3) a[i+1]=a3;
                    else {a[i]=a3;a[i+1]=a1;}
                    a[i+2]=a2;
                }
            } else {
                if(a3>a1) {
                    a[i]=a2;a[i+1]=a1;a[i+2]=a3;
                } else {
                    if(a2>a3) a[i]=a3;
                    else {a[i]=a2;a[i+1]=a3;}
                    a[i+2]=a1;
                }
            }
        }
    }
}

void quicksort_internal(int lo,int hi) {
    int h,i=lo,j=hi-1;
    int pivot=(lo+hi)>>1;

    int x=a[pivot];a[pivot]=a[hi];

    while(a[i]<x&&i<j-1) i++;
    while(a[j]>x&&i+1<j) j--;

    while(i+1<j) {
        h=a[i];a[i++]=a[j];a[j--]=h;
        while(a[i]<x&&i<j-1) i++;
        while(a[j]>x&&i+1<j) j--;
    }

    if(i==j) {
        if(a[j]<x) j++;
        a[hi]=a[j];a[j]=x;
    } else {
        if(a[i]<x) {
            if(a[j]<x) j++;
            a[hi]=a[j];a[j]=x;
        } else {
            if(a[j]<x) {
                a[hi]=a[i];a[i]=a[j];a[j]=x;
            } else {
                a[hi]=a[i];a[i]=x;j=i;
            }
        }
    }
    quicksort(lo,j-1);
    quicksort(j+1,hi);
}

int main() {
    int o,p;
    for (p=0;p<2000;p++) {
        for(o=0;o<n;o++) a[o]=rand()%5+o;
        quicksort(0,n-1);
    }
    //for(o=0;o<n;o++) printf("%d %d\n",o,a[o]);
}
