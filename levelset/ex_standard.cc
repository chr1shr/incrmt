#include "ex_standard.hh"
#include "traverse.hh"

#include <cstdarg>

ex_array::ex_array(int c_,...) : c(c_) {
    va_list vl;
    va_start(vl,c_);
    for(double **fp=f,**fe=f+c;fp<fe;) *(fp++)=va_arg(vl,double*);
    va_end(vl);
}

void ex_array::setup_fields(int c_,...) {
    c=c_;
    va_list vl;
    va_start(vl,c_);
    for(double **fp=f,**fe=f+c;fp<fe;) *(fp++)=va_arg(vl,double*);
    va_end(vl);
}

// Explicit instantiation
template class traverse<tra_positive,ex_single>;
template class traverse<tra_negative,ex_single>;
template class traverse<tra_positive,ex_array>;
template class traverse<tra_negative,ex_array>;
