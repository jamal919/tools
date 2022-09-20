

#ifndef MPI_REDUCERS_H

#define MPI_REDUCERS_H

#include <complex>

using namespace std;

extern int P_sum(int value);
extern size_t P_sum(size_t value);
extern double P_sum(double value);

extern double P_average(double value);

extern double P_norm(double value);

extern double P_max(double value);

extern double P_min(double value);

extern complex<double> P_sum(complex<double> value);

extern int P_maxloc(double & value, int & n);

extern int P_minloc(double & value, int & n);



#endif










