// Level++, a levelset library in C++
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : January 20th 2009

#ifndef LEVELPP_HELPER_HH
#define LEVELPP_HELPER_HH

#include <cstdio>
#include <cstdlib>

#include "config.hh"

void levelpp_fatal_error(const char *p,const int status);

template<class f_type>
void output(const char* filename,f_type *array,int m,int n,double ax,double ay,double dx,double dy);

template <class T>
inline void output(const char* filename,int index,T *q,double m,double n,double ax,double ay,double dx,double dy) {
    char buffer[256];
    sprintf(buffer,"output/%s.%d",filename,index);
    output(buffer,q,m,n,ax,ay,dx,dy);
}

#endif
