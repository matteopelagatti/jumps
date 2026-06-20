// CSSD: Cubic Smoothing Spline with Discontinuities
//
// Implements the exact global minimiser of
//
//   min_{f,J}  p * sum_i ((y_i - f(x_i))/delta_i)^2
//            + (1-p) * integral (f''(t))^2 dt
//            + gamma * |J|
//
// from Storath & Weinmann (2024), JCGS 33(2):651-664.
//
// Algorithm (Sections 2.2-2.5):
//   1. Restrict discontinuity search to midpoints of data sites (Lemma 1).
//   2. Compute all interval energies E[l,r] in O(N^2) via incremental
//      Givens-rotation QR updates on a fixed 5x5 subsystem (Theorem 3).
//   3. Recover the globally optimal discontinuity set by dynamic programming
//      (Bellman equation 7).
//   4. Fit the piecewise cubic spline on each optimal interval using the
//      same QR machinery to get Hermite control points.

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp17)]]
#define EIGEN_DONT_VECTORIZE
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;

// ---------------------------------------------------------------------------
// V and W: the two 2x2 halves of U_i from eq.(9).
//   U_i = beta * [ 2*sqrt(3)/d^{3/2},  sqrt(3)/d^{1/2}, -2*sqrt(3)/d^{3/2},  sqrt(3)/d^{1/2} ]
//                [           0,          1/sqrt(d),                  0,         -1/sqrt(d)       ]
//   V_i = U_i[:,0:1]  acts on [f_i, f'_i]
//   W_i = U_i[:,2:3]  acts on [f_{i+1}, f'_{i+1}]
// ---------------------------------------------------------------------------
static void make_VW(double d, double beta,
                    double V[2][2], double W[2][2]) {
    const double sq3 = std::sqrt(3.0);
    const double sqd  = std::sqrt(d);
    const double d32  = sqd * d;      // d^{3/2}

    V[0][0] =  beta * 2.0 * sq3 / d32;
    V[0][1] =  beta *       sq3 / sqd;
    V[1][0] =  0.0;
    V[1][1] =  beta / sqd;

    W[0][0] = -beta * 2.0 * sq3 / d32;
    W[0][1] =  beta *       sq3 / sqd;
    W[1][0] =  0.0;
    W[1][1] = -beta / sqd;
}

// ---------------------------------------------------------------------------
// Givens rotation: use row ia as pivot to zero entry (ib, col).
// Operates on a (nrows x ncols) matrix M and matching RHS vector b,
// both stored in a single (nrows x (ncols+1)) augmented block passed
// as raw double arrays in row-major order with stride (ncols+1).
// ---------------------------------------------------------------------------
static void givens_zero(double* aug, int nrows, int ncols_aug,
                         int ia, int ib, int col) {
    double a = aug[ia * ncols_aug + col];
    double bv = aug[ib * ncols_aug + col];
    if (std::abs(bv) < 1e-14 * (std::abs(a) + 1e-300)) return;
    double r = std::hypot(a, bv);
    double c = a / r, s = bv / r;
    for (int j = col; j < ncols_aug; ++j) {
        double tmp         = c * aug[ia * ncols_aug + j] + s * aug[ib * ncols_aug + j];
        aug[ib * ncols_aug + j] = -s * aug[ia * ncols_aug + j] + c * aug[ib * ncols_aug + j];
        aug[ia * ncols_aug + j] = tmp;
    }
}

// ---------------------------------------------------------------------------
// QR state: the 2x2 upper-triangular R_bot and 2-vector z_bot that encode
// the "active" bottom-right block needed for the incremental update.
// ---------------------------------------------------------------------------
struct QRState {
    double R00, R01, R11;   // upper-triangular 2x2
    double z0,  z1;         // corresponding RHS
    double energy;          // accumulated sum of squared residuals
};

// ---------------------------------------------------------------------------
// Initialise QR state for the first two points {x_l, x_{l+1}}.
// Builds the 4x4 augmented system A^(2) and QR-decomposes it via Givens.
// The energy is 0 (square full-rank system).
// ---------------------------------------------------------------------------
static QRState init_qr(double yl, double yl1,
                        double d, double al, double al1, double beta) {
    // 4x5 augmented [M | b]
    double aug[4 * 5] = {};
    double V[2][2], W[2][2];
    make_VW(d, beta, V, W);

    // row 0: data at x_l
    aug[0 * 5 + 0] = al;
    aug[0 * 5 + 4] = al * yl;

    // rows 1-2: smoothness on [x_l, x_{l+1}]
    for (int k = 0; k < 2; ++k) {
        aug[(1+k)*5 + 0] = V[k][0];
        aug[(1+k)*5 + 1] = V[k][1];
        aug[(1+k)*5 + 2] = W[k][0];
        aug[(1+k)*5 + 3] = W[k][1];
    }

    // row 3: data at x_{l+1}
    aug[3 * 5 + 2] = al1;
    aug[3 * 5 + 4] = al1 * yl1;

    // QR via Givens: eliminate sub-diagonal entries column by column
    for (int col = 0; col < 4; ++col)
        for (int row = col + 1; row < 4; ++row)
            givens_zero(aug, 4, 5, col, row, col);

    QRState st;
    st.R00    = aug[2*5 + 2];
    st.R01    = aug[2*5 + 3];
    st.R11    = aug[3*5 + 3];
    st.z0     = aug[2*5 + 4];
    st.z1     = aug[3*5 + 4];
    st.energy = 0.0;   // square system -> zero residual
    return st;
}

// ---------------------------------------------------------------------------
// Extend the interval by one point x_{r+1}.
//   d     = x_{r+1} - x_r
//   alpha = sqrt(p) / delta_{r+1}
//   yval  = y_{r+1}
// Updates st in-place (new R_bot, z_bot, energy).
// ---------------------------------------------------------------------------
static void qr_extend(QRState& st,
                       double d, double alpha, double yval, double beta) {
    // 5x5 augmented [M | b]
    //   cols 0-1 : [f_r, f'_r] from R_bot
    //   cols 2-3 : [f_{r+1}, f'_{r+1}] (new)
    //   col  4   : RHS
    double aug[5 * 5] = {};
    double V[2][2], W[2][2];
    make_VW(d, beta, V, W);

    // rows 0-1: existing R_bot
    aug[0*5+0] = st.R00;  aug[0*5+1] = st.R01;  aug[0*5+4] = st.z0;
    aug[1*5+1] = st.R11;                          aug[1*5+4] = st.z1;

    // rows 2-3: smoothness rows for [x_r, x_{r+1}]
    for (int k = 0; k < 2; ++k) {
        aug[(2+k)*5+0] = V[k][0];
        aug[(2+k)*5+1] = V[k][1];
        aug[(2+k)*5+2] = W[k][0];
        aug[(2+k)*5+3] = W[k][1];
    }

    // row 4: data at x_{r+1}
    aug[4*5+2] = alpha;
    aug[4*5+4] = alpha * yval;

    // Givens: eliminate sub-diagonal column by column (4 cols, 5 rows)
    for (int col = 0; col < 4; ++col)
        for (int row = col + 1; row < 5; ++row)
            givens_zero(aug, 5, 5, col, row, col);

    // New R_bot = 2x2 block in rows 2-3, cols 2-3
    st.R00 = aug[2*5+2];
    st.R01 = aug[2*5+3];
    st.R11 = aug[3*5+3];
    st.z0  = aug[2*5+4];
    st.z1  = aug[3*5+4];
    double resid = aug[4*5+4];
    st.energy += resid * resid;
}

// ---------------------------------------------------------------------------
// Solve for Hermite control points [f_i, f'_i] on an interval.
// Returns a (2m)-vector: [f_{start}, f'_{start}, f_{start+1}, f'_{start+1}, ..., f_{end}, f'_{end}].
// ---------------------------------------------------------------------------
static VectorXd solve_hermite(const VectorXd& y, const VectorXd& x,
                                const VectorXd& delta,
                                int start, int end, double p) {
    int m    = end - start + 1;
    int s    = 3 * m - 2;
    double beta = std::sqrt(1.0 - p);

    // Build (s x 2m) matrix A and (s)-vector rhs
    MatrixXd A  = MatrixXd::Zero(s, 2 * m);
    VectorXd rh = VectorXd::Zero(s);

    int row = 0;
    double alpha0 = std::sqrt(p) / delta[start];
    A(row, 0)  = alpha0;
    rh[row]    = alpha0 * y[start];
    ++row;

    for (int i = start; i < end; ++i) {
        double d = x[i+1] - x[i];
        double V[2][2], W[2][2];
        make_VW(d, beta, V, W);
        int col = 2 * (i - start);
        for (int k = 0; k < 2; ++k) {
            A(row+k, col)   = V[k][0];
            A(row+k, col+1) = V[k][1];
            A(row+k, col+2) = W[k][0];
            A(row+k, col+3) = W[k][1];
        }
        row += 2;
        double ai1 = std::sqrt(p) / delta[i+1];
        A(row, col+2) = ai1;
        rh[row]       = ai1 * y[i+1];
        ++row;
    }

    // Solve via Eigen's ColPivHouseholderQR (handles rank-deficient edge cases)
    return A.colPivHouseholderQr().solve(rh);
}

// ---------------------------------------------------------------------------
// Evaluate a piecewise cubic spline at arbitrary query points xq.
// hermite[i] = [f_i, f'_i] for each data point in the interval.
// ---------------------------------------------------------------------------
static void eval_hermite_interval(const VectorXd& xq,
                                   const VectorXd& xi, const VectorXd& hermite,
                                   VectorXd& out) {
    int nq = xq.size();
    int n  = xi.size();
    for (int k = 0; k < nq; ++k) {
        double xv = xq[k];
        // Find containing sub-interval via binary search
        int seg = std::lower_bound(xi.data(), xi.data() + n, xv) - xi.data();
        seg = std::max(0, std::min(n - 2, seg - 1));
        double d  = xi[seg+1] - xi[seg];
        double t  = xv - xi[seg];
        double fi  = hermite[2*seg],      fpi  = hermite[2*seg+1];
        double fi1 = hermite[2*seg+2],    fpi1 = hermite[2*seg+3];
        double c2 = (3.0*(fi1-fi)/d - 2.0*fpi - fpi1) / d;
        double c3 = (2.0*(fi-fi1)/d +     fpi + fpi1) / (d*d);
        out[k] = fi + fpi*t + c2*t*t + c3*t*t*t;
    }
}

// ---------------------------------------------------------------------------
// Main exported function.
//
// @param y_in  numeric vector length n
// @param x_in  strictly increasing numeric vector length n
// @param p     stiffness in (0,1)
// @param gamma jump-count penalty > 0
// @param delta_in  noise std dev vector length n (or length 1 = common sigma)
//
// @return list with:
//   f_hat     : fitted values at each x_i
//   disc_locs : midpoint x-values where discontinuities are placed
//   energy    : optimal target function value (excluding gamma*|J| already added)
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
List cssd_potts(NumericVector y_in, NumericVector x_in, double p, double gamma,
                NumericVector delta_in) {
    int n = y_in.size();
    if (n < 2) stop("Need at least 2 observations.");
    if (x_in.size() != n) stop("x and y must have the same length.");

    Map<const VectorXd> y(y_in.begin(), n);
    Map<const VectorXd> x(x_in.begin(), n);

    // Expand scalar delta to vector
    VectorXd delta(n);
    if (delta_in.size() == 1) {
        delta.setConstant(delta_in[0]);
    } else if (delta_in.size() == n) {
        for (int i = 0; i < n; ++i) delta[i] = delta_in[i];
    } else {
        stop("delta must have length 1 or length(y).");
    }

    double beta = std::sqrt(1.0 - p);

    // -----------------------------------------------------------------------
    // Step 1: compute all interval energies E[l][r] for 0 <= l <= r < n
    // (0-indexed). E[l][r] = energy of fitting a cubic spline to {x_l,...,x_r}
    // using parameters (p, delta).
    // Stored in a flat lower-triangular array: idx = r*n + l.
    // -----------------------------------------------------------------------
    const double INF = std::numeric_limits<double>::infinity();
    std::vector<double> E(n * n, 0.0);

    for (int l = 0; l < n - 1; ++l) {
        double al  = std::sqrt(p) / delta[l];
        double al1 = std::sqrt(p) / delta[l+1];
        double d   = x[l+1] - x[l];

        QRState st = init_qr(y[l], y[l+1], d, al, al1, beta);
        E[( l+1)*n + l] = st.energy;   // E[l, l+1]

        for (int r = l + 2; r < n; ++r) {
            double dr   = x[r] - x[r-1];
            double alphr = std::sqrt(p) / delta[r];
            qr_extend(st, dr, alphr, y[r], beta);
            E[r*n + l] = st.energy;    // E[l, r]
        }
    }
    // E[l][l] = 0 for all l (single point, trivially zero energy)

    // -----------------------------------------------------------------------
    // Step 2: Dynamic programming (Bellman eq. 7, 1-indexed → here 0-indexed).
    // F[r] = min_{l=0,...,r} { E[l,r] + gamma + F[l-1] }
    // with F[-1] = -gamma.
    // Z[r] = argmin l (for backtracking).
    // -----------------------------------------------------------------------
    std::vector<double> F(n, 0.0);
    std::vector<int>    Z(n, 0);

    // r = 0: only choice is l = 0 (single point)
    // F[0] = E[0,0] + gamma + F[-1] = 0 + gamma + (-gamma) = 0
    F[0] = 0.0;
    Z[0] = 0;

    for (int r = 1; r < n; ++r) {
        double best    = INF;
        int    best_l  = 0;
        double F_prev  = -gamma;   // F[l-1] with l=0 → F[-1] = -gamma
        for (int l = 0; l <= r; ++l) {
            double val = E[r*n + l] + gamma + F_prev;
            if (val < best) { best = val; best_l = l; }
            F_prev = (l < n - 1) ? F[l] : INF;
        }
        F[r]   = best;
        Z[r]   = best_l;
    }

    // -----------------------------------------------------------------------
    // Step 3: Backtrack to find optimal intervals.
    // -----------------------------------------------------------------------
    std::vector<std::pair<int,int>> intervals;   // (start, end) indices
    {
        int r = n - 1;
        while (r >= 0) {
            int l = Z[r];
            intervals.emplace_back(l, r);
            r = l - 1;
        }
        std::reverse(intervals.begin(), intervals.end());
    }

    // -----------------------------------------------------------------------
    // Step 4: Fit cubic spline on each optimal interval.
    // -----------------------------------------------------------------------
    VectorXd f_hat(n);
    std::vector<double> disc_locs;

    for (int k = 0; k < (int)intervals.size(); ++k) {
        int st_idx = intervals[k].first;
        int en_idx = intervals[k].second;

        if (st_idx == en_idx) {
            // Single point: fitted value = observed value
            f_hat[st_idx] = y[st_idx];
        } else {
            VectorXd herm = solve_hermite(y, x, delta, st_idx, en_idx, p);
            // Evaluate at all interior data points
            for (int i = st_idx; i <= en_idx; ++i) {
                if (i == st_idx) {
                    f_hat[i] = herm[0];
                } else if (i == en_idx) {
                    f_hat[i] = herm[2*(en_idx - st_idx)];
                } else {
                    int seg = i - st_idx;
                    f_hat[i] = herm[2*seg];
                }
            }
        }

        // Record discontinuity at the midpoint before this interval (k > 0)
        if (k > 0) {
            int prev_end = intervals[k-1].second;
            disc_locs.push_back((x[prev_end] + x[st_idx]) / 2.0);
        }
    }

    // Convert disc_locs to NumericVector
    NumericVector disc_r(disc_locs.begin(), disc_locs.end());

    return List::create(
        Named("f_hat")     = wrap(f_hat),
        Named("disc_locs") = disc_r,
        Named("energy")    = F[n-1] + gamma,   // add back -gamma offset
        Named("n_disc")    = (int)disc_locs.size()
    );
}

// ---------------------------------------------------------------------------
// Evaluate the piecewise cubic spline at arbitrary query points.
// Called after cssd_potts() to get smooth curves on a fine grid.
//
// @param y_in, x_in    training data (sorted)
// @param f_hat_in      fitted values at training points (from cssd_potts)
// @param disc_locs_in  discontinuity locations (from cssd_potts)
// @param xq_in         query points to evaluate at
// @param p             stiffness parameter
// @param delta_in      noise std vector
// ---------------------------------------------------------------------------
// [[Rcpp::export]]
NumericVector cssd_potts_predict(NumericVector y_in, NumericVector x_in,
                                  NumericVector f_hat_in,
                                  NumericVector disc_locs_in,
                                  NumericVector xq_in,
                                  double p, NumericVector delta_in) {
    int n  = y_in.size();
    int nq = xq_in.size();

    Map<const VectorXd> y(y_in.begin(), n);
    Map<const VectorXd> x(x_in.begin(), n);
    Map<const VectorXd> xq(xq_in.begin(), nq);

    VectorXd delta(n);
    if (delta_in.size() == 1)
        delta.setConstant(delta_in[0]);
    else
        for (int i = 0; i < n; ++i) delta[i] = delta_in[i];

    // Reconstruct interval boundaries from disc_locs
    // disc_locs contains midpoints between segments; convert to data-index boundaries
    std::vector<int> seg_start = {0};
    for (double dl : disc_locs_in) {
        // find the first data point index > dl
        int idx = (int)(std::upper_bound(x.data(), x.data() + n, dl) - x.data());
        if (idx > seg_start.back() && idx < n)
            seg_start.push_back(idx);
    }
    seg_start.push_back(n);   // sentinel

    int n_segs = (int)seg_start.size() - 1;
    VectorXd out(nq);
    out.setZero();

    for (int s = 0; s < n_segs; ++s) {
        int l = seg_start[s], r = seg_start[s+1] - 1;
        double xl = x[l], xr = x[r];

        // Query points belonging to this segment
        std::vector<int> qidx;
        for (int k = 0; k < nq; ++k) {
            double xv = xq[k];
            bool in_seg = (s == 0)            ? xv <= xr :
                          (s == n_segs - 1)   ? xv >= xl :
                                                (xv >= xl && xv <= xr);
            // tie-break at disc locations: assign to right segment
            if (s > 0 && std::abs(xv - xl) < 1e-12) in_seg = true;
            if (in_seg) qidx.push_back(k);
        }
        if (qidx.empty()) continue;

        if (l == r) {
            // Single-point segment: constant extrapolation
            for (int k : qidx) out[k] = y[l];
            continue;
        }

        VectorXd herm = solve_hermite(y, x, delta, l, r, p);
        // sub-interval x values
        VectorXd xi_seg = x.segment(l, r - l + 1);
        VectorXd xq_seg(qidx.size());
        for (int i = 0; i < (int)qidx.size(); ++i) xq_seg[i] = xq[qidx[i]];

        VectorXd vals(qidx.size());
        eval_hermite_interval(xq_seg, xi_seg, herm, vals);
        for (int i = 0; i < (int)qidx.size(); ++i) out[qidx[i]] = vals[i];
    }

    return wrap(out);
}
