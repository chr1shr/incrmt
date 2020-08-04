#include <cstdio>
#include <cstdlib>
#include <cmath>

const int n=100000;
int c[n];

void quicksort_internal(int *a,int hi);

inline void quicksort(int *a,int k) {
    int a1;
    if (k<3) {
        if (k==2&&a[0]>a[1]) {a1=a[0];a[0]=a[1];a[1]=a1;}
    } else {
        if (k>3) quicksort_internal(a,k);
        else {
            a1=a[0];int a2=a[1],a3=a[2];
            if(a1<a2) {
                if(a3<a2) {
                    if(a1<a3) a[1]=a3;
                    else {a[0]=a3;a[1]=a1;}
                    a[2]=a2;
                }
            } else {
                if(a3>a1) {
                    a[0]=a2;a[1]=a1;a[2]=a3;
                } else {
                    if(a2>a3) a[0]=a3;
                    else {a[0]=a2;a[1]=a3;}
                    a[2]=a1;
                }
            }
        }
    }
}

void quicksort_internal(int *a,int k) {
    k--;int pivot=k>>1;
    int h,x=a[pivot],i=0,j=k-1;a[pivot]=a[k];

    while(a[i]<x&&i<j-1) i++;
    while(a[j]>x&&i+1<j) j--;

    while(i+1<j) {
        h=a[i];a[i++]=a[j];a[j--]=h;
        while(a[i]<x&&i<j-1) i++;
        while(a[j]>x&&i+1<j) j--;
    }

    if(i==j) {
        if(a[j]<x) j++;
        a[k]=a[j];a[j]=x;
    } else {
        if(a[i]<x) {
            if(a[j]<x) j++;
            a[k]=a[j];a[j]=x;
        } else {
            if(a[j]<x) {
                a[k]=a[i];a[i]=a[j];a[j]=x;
            } else {
                a[k]=a[i];a[i]=x;j=i;
            }
        }
    }

    quicksort(a,j);
    quicksort(a+(j+1),k-j);
}

int main() {
    int o,p;
    for (p=0;p<2000;p++) {
        for(o=0;o<n;o++) c[o]=rand()%5+o;
        quicksort(c,n);
    }
//    for(o=0;o<n;o++) printf("%d %d\n",o,c[o]);
}
