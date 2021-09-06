#ifndef WLS_HH
#define WLS_HH

namespace ivanp {

void wls(
  const double* A, // matrix of functions values: p (row) then x (col)
  const double* y, // measured values
  const double* u, // variances
  unsigned nx, // number of measured values
  unsigned np, // number of parameters
  double* p, // fitted functions coefficients (parameters)
  double* cov = nullptr // covariance matrix
);

}

#endif
