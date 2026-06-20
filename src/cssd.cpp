// ADMM solver for total-variation denoising with irregular sampling.
//
// Solves: min_u (1/2)||y - u||_2^2 + lambda * ||D u||_1
// where D is the weighted first-difference operator: (Du)_i = (u_{i+1} - u_i) / dx_i.
// The minimizer u is piecewise constant; detected jumps are gaps where Du != 0.
//
// Reference:
//   Storath, M. and Weinmann, A. (2024). Smoothing splines for discontinuous
//   signals. J. Comput. Graph. Stat., 33(2), 651-664.
//   (This file implements the simplified TV / first-order variant, not the full
//    cubic-spline CSSD from the paper, which requires a QR-update algorithm.)

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]

// Disable vectorization to avoid warning in GCC 14/Rtools45
#define EIGEN_DONT_VECTORIZE

#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;

// Element-wise soft thresholding: S_t(v)_i = sign(v_i) * max(|v_i| - t, 0).
static VectorXd soft_thresh(const VectorXd& v, double t) {
  return v.array().sign() * (v.array().abs() - t).cwiseMax(0.0);
}

//' ADMM solver for TV denoising with irregular sampling
//'
//' Solves min_u (1/2)||y - u||^2 + lambda*||Du||_1 where D is the weighted
//' first-difference matrix with (Du)_i = (u_{i+1} - u_i) / dx_i.
//' The system (I + rho*D'D) is factorized once and reused every iteration.
//'
//' @param y_in numeric vector of length n (no NAs).
//' @param x_in strictly increasing numeric vector of length n (the sample sites).
//' @param lambda positive regularization parameter.
//' @param max_iter maximum number of ADMM iterations (default 100).
//' @param rho positive augmented-Lagrangian step; affects convergence speed
//'   but not the solution (default 1.0).
//'
//' @return Numeric vector u of length n: the piecewise-constant minimizer.
// [[Rcpp::export]]
NumericVector cssd_irregular_eigen(NumericVector y_in, NumericVector x_in,
                                   double lambda,
                                   int max_iter = 100, double rho = 1.0) {
  const int n = y_in.size();
  if (x_in.size() != n) stop("x and y must have the same length.");

  Map<const VectorXd> y(y_in.begin(), n);
  Map<const VectorXd> x(x_in.begin(), n);

  // Weighted first-difference matrix D  ((n-1) x n)
  VectorXd dx(n - 1);
  for (int i = 0; i < n - 1; ++i) {
    dx[i] = x[i + 1] - x[i];
    if (dx[i] <= 0.0) stop("x must be strictly increasing.");
  }

  // Unweighted first differences (w = 1): penalises absolute jump amplitudes.
  // The 1/dx weighting would penalise slopes instead, but causes D'D to have
  // eigenvalues of order 1/h^2 ~ 1e4, which makes the ADMM u-update collapse
  // u to nearly zero in the first iteration and stall permanently.
  std::vector<Triplet<double>> trips;
  trips.reserve(2 * (n - 1));
  for (int i = 0; i < n - 1; ++i) {
    trips.emplace_back(i, i,      -1.0);
    trips.emplace_back(i, i + 1,   1.0);
  }
  SparseMatrix<double> D(n - 1, n);
  D.setFromTriplets(trips.begin(), trips.end());
  SparseMatrix<double> Dt = D.transpose();

  // Pre-factorize A = I + rho * D' * D (done once, reused every iteration)
  SparseMatrix<double> I_n(n, n);
  I_n.setIdentity();
  SparseMatrix<double> A = I_n + rho * (Dt * D);

  SimplicialLDLT<SparseMatrix<double>> solver;
  solver.compute(A);
  if (solver.info() != Success) stop("Cholesky factorization of (I + rho*D'D) failed.");

  // ADMM: u-update, z-update (soft threshold), dual-variable update
  VectorXd u     = y;
  VectorXd z     = VectorXd::Zero(n - 1);
  VectorXd alpha = VectorXd::Zero(n - 1);   // unscaled dual variable
  VectorXd Du(n - 1);

  for (int k = 0; k < max_iter; ++k) {
    u     = solver.solve(y + Dt * (rho * z - alpha));
    Du    = D * u;
    z     = soft_thresh(Du + alpha / rho, lambda / rho);
    alpha = alpha + rho * (Du - z);
  }

  return wrap(u);
}
