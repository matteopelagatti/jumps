// It contains the functions for the Alternating time-warping algorithm
// for cubic splines with discontinuities.

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]

// Disable vectorization to avoid warning in GCC 14/Rtools45
#define EIGEN_DONT_VECTORIZE

#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;

// Helper to invert a 2x2 matrix
inline Matrix2d inv2x2(const Matrix2d& M) {
  double det = M(0,0)*M(1,1) - M(0,1)*M(1,0);
  double invDet = 1.0 / det;
  Matrix2d res;
  res <<  M(1,1) * invDet, -M(0,1) * invDet,
          -M(1,0) * invDet,  M(0,0) * invDet;
  return res;
}

// Projection onto the simplex.
// Writes the result into `out`; uses `work` as a scratch sort buffer.
// Fuses the cumsum and threshold loops to avoid a separate sv[] allocation.
void project_simplex_inplace(const VectorXd& v, double budget,
                              VectorXd& out, VectorXd& work) {
  // budget <= 0: the only feasible point of {x >= 0, sum(x) <= budget} is 0.
  // The general algorithm below relies on a *strict* inequality to find the
  // threshold theta and never fires when budget == 0 exactly, leaving theta
  // at its initial value of 0 and letting v.cwiseMax(0) through unconstrained.
  if (budget <= 0.0) {
    out.setZero();
    return;
  }
  if (v.sum() <= budget && (v.array() >= 0).all()) {
    out = v.cwiseMax(0.0);
    return;
  }
  work = v;
  std::sort(work.data(), work.data() + work.size(), std::greater<double>());
  double cumsum = 0.0, theta = 0.0;
  for (int i = 0; i < (int)work.size(); ++i) {
    cumsum += work[i];
    double val = (cumsum - budget) / (i + 1.0);
    if (work[i] > val) theta = val;
  }
  out = (v.array() - theta).cwiseMax(0.0);
}

//' Alternating time-warping algorithm for smoothing splines with jumps
//'
//' @param x_in vector of observations for the x variable
//' @param y_in vector of observations for the y variable
//' @param lambda positive smoothing constant
//' @param M positive value for the maximum sum of the time extensions
//' @param max_iter positive integer with the maximum number of iterations
//' @param tol small positive number with tolerance for numerical convergence
//' @param learning_rate positive number: initial and fallback step size for the
//'   Barzilai-Borwein adaptive scheme (default 0.1).
//' @param gamma_init optional numeric vector of length n-1 with a warm-start
//'   value for gamma; ignored (zero initialisation used) when its length differs
//'   from n-1.
//'
//' @returns A list with the following slots:
//' \describe{
//'   \item{f}{The smoothing spline with jumps computed at x.}
//'   \item{gamma}{The vector of added times (time jumps).}
//'   \item{df}{Derivative of f.}
//'   \item{edf}{Number of effective degrees of freedom.}
//'   \item{rss}{Residual sum of squares with residuals = y - f.}
//'   \item{gcv}{Generalized cross-validation.}
//'   \item{loglik}{Gaussian log-likelihood evaluated at the unbiased residual variance sigma^2 = rss/(n - edf).}
//'   \item{aic}{Akaike information criterion.}
//'   \item{aicc}{Akaike information criterion with bias correction.}
//'   \item{bic}{Bayesian information criterion.}
//'   \item{hq}{Hannan-Quinn information criterion.}
//' }
// [[Rcpp::export]]
List solve_jump_spline_fast(NumericVector x_in, NumericVector y_in,
                            double lambda, double M,
                            int max_iter = 100, double tol = 1e-4,
                            double learning_rate = 0.1,
                            NumericVector gamma_init = NumericVector(0)) {

  Map<VectorXd> x(as<Map<VectorXd> >(x_in));
  Map<VectorXd> y(as<Map<VectorXd> >(y_in));
  int n = x.size();

  VectorXd delta_x(n - 1);
  for (int i = 0; i < n - 1; ++i) delta_x[i] = x[i+1] - x[i];

  VectorXd f_curr = y;
  VectorXd df_curr = VectorXd::Zero(n);
  VectorXd rhs = VectorXd::Zero(2*n);
  for (int i = 0; i < n; ++i) rhs(2*i) = y(i);

  std::vector<Matrix2d> D_blocks(n);
  std::vector<Matrix2d> U_blocks(n - 1);
  std::vector<Matrix2d> D_tilde_inv(n);

  // Pre-allocate all loop temporaries once to avoid per-iteration heap churn.
  MatrixXd z_mat(2, n);
  VectorXd x_sol(2 * n);
  VectorXd grad_gamma(n - 1);
  VectorXd gamma_unconstrained(n - 1);
  VectorXd gamma_new(n - 1);
  VectorXd sort_buf(n - 1);
  VectorXd gamma_prev(n - 1);     // Barzilai–Borwein: previous iterate
  VectorXd grad_prev(n - 1);      // Barzilai–Borwein: previous gradient
  double alpha = learning_rate;   // BB step size (adapts each iteration)

  // Initialise gamma: from warm start (project onto feasible set) or zero.
  VectorXd gamma(n - 1);
  if (gamma_init.size() == n - 1) {
    gamma = as<VectorXd>(gamma_init);
    project_simplex_inplace(gamma, M, gamma_new, sort_buf);
    gamma = gamma_new;
  } else {
    gamma.setZero();
  }
  gamma_prev.setZero();
  grad_prev.setZero();

  double edf_final = 0.0;

  // --- Optimization loop ---
  for (int iter = 0; iter < max_iter; ++iter) {

    // 1. Build matrices
    for (int i = 0; i < n; ++i) D_blocks[i] << 1.0, 0.0, 0.0, 0.0;

    for (int i = 0; i < n - 1; ++i) {
      double tau = delta_x[i] + gamma[i];
      double t2 = tau*tau, t3 = t2*tau;

      double lk11 = lambda*12./t3, lk12 = lambda*6./t2;
      double lk22  = lambda*4./tau, lk22_off = lambda*2./tau;

      D_blocks[i](0,0) += lk11; D_blocks[i](0,1) += lk12;
      D_blocks[i](1,0) += lk12; D_blocks[i](1,1) += lk22;

      D_blocks[i+1](0,0) += lk11; D_blocks[i+1](0,1) -= lk12;
      D_blocks[i+1](1,0) -= lk12; D_blocks[i+1](1,1) += lk22;

      Matrix2d U;
      U << -lk11, lk12, -lk12, lk22_off;
      U_blocks[i] = U;
    }

    // 2. Forward Elimination
    D_tilde_inv[0] = inv2x2(D_blocks[0]);
    z_mat.col(0) = D_tilde_inv[0] * rhs.segment<2>(0);

    for (int i = 1; i < n; ++i) {
      Matrix2d L_fact = D_tilde_inv[i-1] * U_blocks[i-1];
      D_blocks[i] -= U_blocks[i-1].transpose() * L_fact;
      D_tilde_inv[i] = inv2x2(D_blocks[i]);
      z_mat.col(i) = D_tilde_inv[i] *
                     (rhs.segment<2>(2*i) - U_blocks[i-1].transpose() * z_mat.col(i-1));
    }

    // 3. Backward Substitution
    x_sol.segment<2>(2*(n-1)) = z_mat.col(n-1);
    for (int i = n - 2; i >= 0; --i) {
      x_sol.segment<2>(2*i) = z_mat.col(i) -
                               (D_tilde_inv[i] * U_blocks[i]) * x_sol.segment<2>(2*(i+1));
    }

    for (int i = 0; i < n; ++i) {
      f_curr(i)  = x_sol(2*i);
      df_curr(i) = x_sol(2*i + 1);
    }

    // 4. Gradient + Barzilai–Borwein projected step
    for (int i = 0; i < n - 1; ++i) {
      double tau = delta_x[i] + gamma[i];
      double D = f_curr(i+1) - f_curr(i);
      double S = df_curr(i) + df_curr(i+1);
      double P = df_curr(i)*df_curr(i) + df_curr(i)*df_curr(i+1) + df_curr(i+1)*df_curr(i+1);
      double t2 = tau*tau, t3 = t2*tau, t4 = t3*tau;
      grad_gamma(i) = lambda * (-36.*D*D/t4 + 24.*D*S/t3 - 4.*P/t2);
    }

    // BB1 step: reuse gamma_unconstrained / gamma_new as scratch to avoid
    // allocating ds and dg vectors (both are overwritten immediately after).
    if (iter > 0) {
      gamma_unconstrained.noalias() = gamma - gamma_prev;   // ds
      gamma_new.noalias()           = grad_gamma - grad_prev; // dg
      double ss = gamma_unconstrained.squaredNorm();
      double sg = gamma_unconstrained.dot(gamma_new);
      if (ss > 0.0 && sg > 1e-10 * ss) alpha = ss / sg;
    }
    gamma_prev = gamma;
    grad_prev  = grad_gamma;

    gamma_unconstrained.noalias() = gamma - alpha * grad_gamma;
    project_simplex_inplace(gamma_unconstrained, M, gamma_new, sort_buf);

    if ((gamma_new - gamma).norm() < tol) break;
    gamma = gamma_new;
  }

  // --- 5. Takahashi selected inversion for EDF ---
  Matrix2d Sigma_next = D_tilde_inv[n-1];
  edf_final = Sigma_next(0,0);

  for (int i = n - 2; i >= 0; --i) {
    Matrix2d P = D_tilde_inv[i] * U_blocks[i];
    Matrix2d Sigma_curr = D_tilde_inv[i] + P * Sigma_next * P.transpose();
    edf_final += Sigma_curr(0,0);
    Sigma_next = Sigma_curr;
  }

  // --- 6. Statistical criteria ---
  double rss = (y - f_curr).squaredNorm();
  double mse = rss / n;

  double gcv = mse / std::pow(1.0 - edf_final / n, 2);

  // Log-likelihood at sigma^2 = rss/(n - edf): does not collapse to +Inf
  // for near-perfect fits the way the concentrated (sigma^2 = rss/n) form does.
  double df_resid   = n - edf_final;
  double sigma2_hat = (df_resid > 0) ? rss / df_resid : mse;
  double log_lik    = -0.5 * n * std::log(2 * M_PI)
                    - 0.5 * n * std::log(sigma2_hat)
                    - 0.5 * rss / sigma2_hat;

  double aic  = n * std::log(mse) + 2 * edf_final;
  double aicc = aic + (2 * edf_final * (edf_final + 1)) / (n - edf_final - 1);
  double bic  = n * std::log(mse) + std::log((double)n) * edf_final;
  double hq   = n * std::log(mse) + 2 * std::log(std::log((double)n)) * edf_final;

  return List::create(
    Named("f")      = f_curr,
    Named("gamma")  = gamma,
    Named("df")     = df_curr,
    Named("edf")    = edf_final,
    Named("rss")    = rss,
    Named("gcv")    = gcv,
    Named("loglik") = log_lik,
    Named("aic")    = aic,
    Named("aicc")   = aicc,
    Named("bic")    = bic,
    Named("hq")     = hq
  );
}
