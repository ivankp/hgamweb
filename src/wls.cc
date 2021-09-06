#include "wls.hh"
#include "linalg.hh"
#include <cstring>

// Cowan p. 97
// https://en.wikipedia.org/wiki/Generalized_least_squares

// p = (At V^-1 A)^-1 At V^-1 y

namespace ivanp {
using namespace ivanp::linalg;

void wls(
  const double* A, // matrix of functions values: p (row) then x (col)
  const double* y, // measured values
  const double* u, // variances
  unsigned nx, // number of measured values
  unsigned np, // number of parameters
  double* p, // fitted functions coefficients (parameters)
  double* cov // covariance matrix
) {
  const unsigned N = utn(np);
  double* const L = new double[N+nx];
  double* const V = L + N; // V^-1

  // ================================================================
  // V^-1 = (u^2)^-1
  // LL = At V^-1 A
  // p = At V^-1 y
  // solve p = L^-1 p
  // solve p = LT^-1 p
  // ================================================================

  // V^-1 = (u^2)^-1
  if (u) {
    for (unsigned i=nx; i; ) { --i; V[i] = 1./u[i]; }
  } else {
    for (unsigned i=nx; i; ) { --i; V[i] = 1.; }
  }

  // LL = At V^-1 A
  for (unsigned i=0, r1=0, k=0; i<np; ++i, r1+=nx) {
    for (unsigned j=0, r2=0; j<=i; ++j, r2+=nx, ++k) {
      double& l = L[k] = 0;
      for (unsigned x=nx; x; ) {
        --x;
        l += A[r1+x] * A[r2+x] * V[x];
      }
    }
  }

  cholesky(L,N);

  // p = At V^-1 y
  for (unsigned i=nx; i; ) { --i; V[i] *= y[i]; }

  for (unsigned i=0, row=0; i<np; ++i, row+=nx) {
    double& pi = p[i] = 0;
    for (unsigned x=nx; x; ) {
      --x;
      pi += A[row+x] * V[x];
    }
  }

  solve_triang  (L,p,np); // solve p = L^-1 p
  solve_triang_T(L,p,np); // solve p = LT^-1 p

  if (cov) {
    inv_triang(L,np); // invert
    LT_L(L,np); // multiply LT by L
    memcpy(cov,L,N*sizeof(double));
  }

  delete[] L;
}

}
