#ifndef COMMON_HH
#define COMMON_HH

#include <cstdio>
#include <cstdlib>
#include <limits>

#ifdef _OPENMP
#include "omp.h"
#else
#include <ctime>
#endif

/** A small number that is a few orders of magnitude larger than machine
 * epsilon, and is used as a threshold for various iterations in the code. */
const double small_number=100.*std::numeric_limits<double>::epsilon();

/** A large number that is used as an initial value in searches for minimal
 * values. */
const double big_number=std::numeric_limits<double>::max();

FILE* safe_fopen(const char* filename,const char* mode);
void fatal_error(const char *p,int code);

/** Calculates a random number uniformly distributed in a given range.
 * \param[in] (a,b) the range.
 * \return The computed number. */
inline double rnd(double a,double b) {
    return a+((b-a)/RAND_MAX)*static_cast<double>(rand());
}

// Set up timing routine. If code was compiled with OpenMP, then use the
// accurate wtime function. Otherwise use the clock function in the ctime
// library.
#ifdef _OPENMP
inline double wtime() {return omp_get_wtime();}
#else
inline double wtime() {return (1./CLOCKS_PER_SEC)*static_cast<double>(clock());}
#endif

#endif
