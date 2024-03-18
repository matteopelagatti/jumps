#include <Rcpp.h>
using namespace Rcpp;

//' Internal function for computing scores w/r to regression coefficients
//' 
//' This function, not intended for end-users, implements the following
//' recursions needed in computing scores with respect to regression
//' coefficients:
//' \deqn{D a^{(1)}_{t+1} = D a^{(1)}_{t} + D a^{(2)}_{t} - k^{(1)}_t x_t -
//' k^{(1)}_t D a^{(1)}_{t}}
//' \deqn{D a^{(2)}_{t+1} = a^{(2)}_{t} - k^{(2)}_t x_t - k^{(2)}_t Da^{(1)}_{t}}
//' where \eqn{a^{(1)}_{t}}, \eqn{a^{(2)}_{t}} are the one-step-ahead Kalman filtered
//' state variables, and \eqn{k^{(1)}_{t}}, \eqn{k^{(2)}_{t}} the respective
//' Kalman gain elements. The symbol $D$ represent the partial derivative with
//' respect to the regression coefficients and $x_t$ is the vector of regressors.
//' All variables are passed by reference and, so, no output is needed.
//' 
//' @param k1 numeric vector of n elements with the Kalman gain sequence for
//' the first state variable;
//' @param k2 numeric vector of n elements with the Kalman gain sequence for
//' the second state variable;
//' @param X numeric matrix of dimension \eqn{n\times k} with the regressors;
//' @param A1 numeric matrix of dimension \eqn{n\times k} that, after calling
//' the function will contain the sequence of gradients \eqn{D a^{(1)}_t}; the
//' first row must be of zero values;
//' @param A2 numeric matrix of dimension \eqn{n\times k} that, after calling
//' the function will contain the sequence of gradients \eqn{D a^{(2)}_t}; the
//' first row must be of zero values;
//' @returns It does not return anything as it writes on the A1 and A2 matrices
//' passed as reference.
// [[Rcpp::export]]
void da(NumericVector k1,
       NumericVector k2,
       NumericMatrix X,
       NumericMatrix A1,
       NumericMatrix A2) {
 unsigned int n = A1.nrow();
 for (unsigned int t = 0; t < n-1; ++t) {
   A1(t+1, _) = (1-k1[t])*A1(t, _) + A2(t, _) - k1[t]*X(t, _);
   A2(t+1, _) = A2(t, _) - k2[t]*X(t, _) - k2[t]*A1(t, _);
 }
}