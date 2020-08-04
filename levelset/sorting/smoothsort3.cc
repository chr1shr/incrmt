#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "leonardo.cc"

const int nn=100000;
int a[nn];

inline void sm_sift(int *v,int l,int r) {
	int rm,r2,d;
	while(l>=2) {
		rm=r-1;r2=rm-smc[l-1];
		if(v[r2]<v[rm]) {
			if(v[r]<v[rm]) {
				d=v[rm];v[rm]=v[r];v[r]=d;
				l-=2;r=rm;
			} else return;
		} else {
			if(v[r]<v[r2]) {
				d=v[r2];v[r2]=v[r];v[r]=d;
				l--;r=r2;
			} else return;
		}
	}
}

inline void sm_trinkle(int *v,int p,int l,int r) {
	int r2,r3,d;
	while(p>0) {
		while((p&1)==0) {
			if((p&2)==0) {p>>=2;l+=2;}
			else {p>>=1;l++;break;}
		}
		r3=r-smb[l];
		if(p>1&&v[r3]>v[r]) {
			p--;
			if(l<=1) {
				d=v[r3];v[r3]=v[r];v[r]=d;r=r3;
			} else {
				r2=r-smd[l];
				if(v[r2]<v[r-1]) {
					r2=r-1;l--;p<<=1;
				}
				if(v[r3]<v[r2]) {
					d=v[r];v[r]=v[r2];v[r2]=d;r=r2;
					l--;sm_sift(v,l,r);return;
				} else {
					d=v[r];v[r]=v[r3];v[r3]=d;r=r3;
				}
			}
		} else break;
	}
	sm_sift(v,l,r);
}

inline void sm_semitrinkle(int *v,int &p,int &l,int &r) {
	int r1=r-smc[l];
	if (v[r1]>v[r]) {
		int d=v[r];v[r]=v[r1];v[r1]=d;sm_trinkle(v,p,l,r1);
	}
}

void smoothsort(int *v,int n) {
	int q,r=0,p=1,l=1;
	for(q=1;q<n;q++) {
		if((p&7)==3) {
			sm_sift(v,l,r);
			p=(p+1)>>2;l+=2;
		} else {
			if(q+smc[l]<n) sm_sift(v,l,r);else sm_trinkle(v,p,l,r);
			if (l==1) {p<<=1;l=0;} else {p<<=l-1;l=1;}
			p++;
		}
		r++;
	}
	sm_trinkle(v,p,l,r);
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

			if(p>0) sm_semitrinkle(v,p,l,r);
			l--;p=(p<<1)+1;

			r=r+smc[l];sm_semitrinkle(v,p,l,r);
			l--;p=(p<<1)+1;
		}
	}
}

int main() {
	int o,p;
	for (p=0;p<2000;p++) {
		for(o=0;o<nn;o++) a[o]=rand()%5+o;
		smoothsort(a,nn);
	}
//	for(o=0;o<nn;o++) printf("%d %d\n",o,a[o]);
}
