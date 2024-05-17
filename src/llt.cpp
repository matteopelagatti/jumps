#include <Rcpp.h>
using namespace Rcpp;

//' Kalman filtering and smoothing for local linear trend plus noise
//' 
//' It uses the power of C++, scalar computation and pointers
//' to run the Kalman filter, the smoother and compute the log-likelihood.
//' The R user has to supply many vectors that in most cases will be
//' overwritten by the llt() function since they are passed by reference.
//' All passed parameters must be numerical (floating point) vectors:
//' any other kind of variable may cause serious problems to the stability
//' of your system. Passing vectors of integers will make the computations fail.
//'  
//' @param y vector of n observations
//' @param var_eps vector of n variances for the observation noises
//' @param var_eta vector of n variances for the level disturbances
//' @param var_zeta vector of n variances for the slope disturbances
//' @param cov_eta_zeta vector of n covariances between level and slope disturbances
//' @param a1 vector of n+1 one-step-ahead prediction of the level; the first
//' element is the initial condition for the level at time t=1, the other elements
//' are arbitrary and will be overwritten
//' @param a2 vector of n+1 one-step-ahead prediction of the slope; the first
//' element is the initial condition for the slope at time t=1, the other elements
//' are arbitrary and will be overwritten
//' @param p11 vector of n+1 one-step-ahead prediction error variance of the level;
//' the first element is the initial condition for the level at time t=1, the other elements
//' are arbitrary and will be overwritten
//' @param p12 vector of n+1 one-step-ahead prediction covariances for level and slope;
//' the first element is the initial condition for the slope at time t=1, the other elements
//' are arbitrary and will be overwritten
//' @param p22 vector of n+1 one-step-ahead prediction error variance of the slope;
//' the first element is the initial condition for the level at time t=1, the other elements
//' are arbitrary and will be overwritten
//' @param k1 vector of the n Kalman gains for the level equation; values are
//' arbitrary and will be overwritten;
//' @param k2 vector of the n Kalman gains for the slope equation; values are
//' arbitrary and will be overwritten;
//' @param i vector of the n innovations; values are
//' arbitrary and will be overwritten;
//' @param f vector of the n innovatoin variances;values are
//' arbitrary and will be overwritten;
//' @param r1 vector of the n+1 smoothers (Th.5.4 in Pelagatti, 2015)
//' for the level equation; values are arbitrary and will be overwritten;
//' @param r2 vector of the n+1 smoothers (Th.5.4 in Pelagatti, 2015) for
//' the slope equation; values are arbitrary and will be overwritten;
//' @param n11 vector of the n+1 variance smoothers (Th.5.4 in Pelagatti, 2015) for the level
//' equation; values are arbitrary and will be overwritten;
//' @param n12 vector of the n+1 covariance smoothers (Th.5.4 in Pelagatti, 2015) for the level
//' and slope; values are arbitrary and will be overwritten;
//' @param n22 vector of the n+1 variance smoothers (Th.5.4 in Pelagatti, 2015) for the slope
//' equation; values are arbitrary and will be overwritten;
//' @param e vector of the n+1 observation error smoothers (Th.5.4 in Pelagatti, 2015);
//' values are arbitrary and will be overwritten;
//' @param d vector of the n+1 observation error variance smoothers (Th.5.4 in Pelagatti, 2015);
//' values are arbitrary and will be overwritten;
//' @param w_ NULL (default) or vector of n weights for the effect of observation y_t
//' on the estimation of the hp filter (with jumps) at time t;
//' values are arbitrary and will be overwritten;
//' @return The value of the Gaussian log-likelihood net of the -log(2*pi)*n/2 part
//' that can be added if needed.
// [[Rcpp::export]]
 double llt(NumericVector y,
            NumericVector var_eps,
            NumericVector var_eta,
            NumericVector var_zeta,
            NumericVector cov_eta_zeta,
            NumericVector a1,
            NumericVector a2,
            NumericVector p11,
            NumericVector p12,
            NumericVector p22,
            NumericVector k1,
            NumericVector k2,
            NumericVector i,
            NumericVector f,
            NumericVector r1,
            NumericVector r2,
            NumericVector n11,
            NumericVector n12,
            NumericVector n22,
            NumericVector e,
            NumericVector d,
            Nullable<NumericVector> w_ = R_NilValue) {
   unsigned int n = y.length();
   NumericVector w;
   if (w_.isNotNull()) w = w_;
   // Kalman filter (only from t|t-1 to t+1|t)
   i[0] = y[0] - a1[0];
   f[0] = p11[0] + var_eps[0];
   double loglik = log(f[0]) + i[0]*i[0]/f[0];
   k1[0] = (p11[0] + p12[0]) / f[0];
   k2[0] = p12[0] / f[0];
   for (unsigned int t = 0; t < n; ++t) {
     // one-step-prediction
     a1[t+1] = a1[t] + a2[t] + k1[t]*i[t];
     a2[t+1] =         a2[t] + k2[t]*i[t];
     // their variances/covariance
     if (NumericVector::is_na(y[t])) { // in case of missing y[t]
       p11[t+1] = p11[t] + 2*p12[t] + p22[t] + var_eta[t];
       p12[t+1] = p12[t] + p22[t] + cov_eta_zeta[t];
       p22[t+1] = p22[t] + var_zeta[t];
     } else {
       p11[t+1] = p11[t] + 2*p12[t] + p22[t] + var_eta[t] - k1[t]*k1[t]*f[t];
       p12[t+1] = p12[t] + p22[t] + cov_eta_zeta[t] - k1[t]*k2[t]*f[t];
       p22[t+1] = p22[t] + var_zeta[t] - k2[t]*k2[t]*f[t];
     }
     if (t < n-1) {
       if (NumericVector::is_na(y[t+1])) { // in case of missing y[t+1]
         // innovation
         i[t+1] = 0;
         f[t+1] = R_PosInf;
         // Kalman gain
         k1[t+1] = 0;
         k2[t+1] = 0;
       } else {
         // innovation
         i[t+1] = y[t+1] - a1[t+1];
         f[t+1] = p11[t+1] + var_eps[t+1];
         // Kalman gain
         k1[t+1] = (p11[t+1] + p12[t+1]) / f[t+1];
         k2[t+1] = p12[t+1] / f[t+1];
         // Gaussian negative log-lik kernel
         loglik += log(f[t+1]) + i[t+1]*i[t+1]/f[t+1];
       }
     }
   }
   // smoother
   r1[n] = 0;
   r2[n] = 0;
   n11[n] = 0;
   n12[n] = 0;
   n22[n] = 0;
   e[n-1] = i[n-1]/f[n-1];
   d[n-1] = 1/f[n-1];
   double k;
   double ksq;
   double kk;
   for (unsigned int t = n-1; t >= 0; --t) {
     k = (1 - k1[t]);
     ksq = k*k;
     kk = k2[t]*k;
     r1[t] = i[t]/f[t] + k*r1[t+1] - k2[t]*r2[t+1];
     r2[t] = r1[t+1] + r2[t+1];
     n11[t] = ksq*n11[t+1] - 2*kk*n12[t+1] + k2[t]*k2[t]*n22[t+1] + 1/f[t];
     n12[t] = k*(n11[t+1] + n12[t+1]) - k2[t]*(n12[t+1] + n22[t+1]);
     n22[t] = n11[t+1] + 2*n12[t+1] + n22[t+1];
     if (w_.isNotNull()) {
       w[t] = p11[t]*(1/f[t] + k1[t]*k1[t]*n11[t] + 2*k1[t]*k2[t]*n12[t] +
         k2[t]*k2[t]*n22[t] - k1[t]*n11[t] - k2[t]*n12[t]) -
         p12[t]*(k1[t]*(n11[t] + n12[t]) + k2[t]*(n12[t] + n22[t]));
     }
     if (t == 0) break;
     e[t-1] = i[t-1]/f[t-1] - k1[t-1]*r1[t] - k2[t-1]*r2[t];
     d[t-1] = 1/f[t-1] + k1[t-1]*k1[t-1]*n11[t] +
       2*k1[t-1]*k2[t-1]*n12[t] + k2[t-1]*k2[t-1]*n22[t];
   }
   return -loglik/2.0;
 }
 