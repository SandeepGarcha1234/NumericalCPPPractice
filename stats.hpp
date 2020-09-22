#ifndef STATS.H
#define STATS.H

#include "matrix.hpp"

double normcdf(double x);
double norminv(double x);
Matrix randuniform(int n);
double boxMuller();
Matrix boxMuller(int N);


#endif