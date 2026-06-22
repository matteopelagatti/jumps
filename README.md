# Hodrick-Prescott filter with automatically selected jumps

This R package implements our novel method to supplement the classical HP filter
with jumps and, possibly, regressors. The method is based on the following state-space
representation

$$y_t = x_t^\top \beta + \mu_t + \varepsilon_t$$

$$\mu_{t+1} = \mu_t + \nu_t$$

$$\nu_{t+1} = \nu_t + \zeta_t,$$

where $y_t$ is the observable time series, $\mu_t$ is the level component,
$\nu_t$ is the slope component, $\varepsilon_t$ and $\zeta_t$ are white noise sequences
with variances $\sigma^2_\varepsilon$ and $\sigma^2_\zeta$, respectively.
The smoother, that is, the linear projection of $\mu_t$ on the span of the observations
$\{y_1,\ldots,y_n\}$, coincides with the HP filter, where the smoothing constant $\lambda$
is given by $\sigma^2_\varepsilon / \sigma^2_\zeta$. Finally, $x_t$ is a vector of regressors,
and $\beta$ is a vector of regression coefficients. These regressors are mainly used to
model seasonal patterns in the data and should have a zero mean to not alter the interpretation
of the HP filter as a trend extractor.

# Smoothing splines with discontinuities

This part of the package implements two novel methods for fitting cubic smoothing splines
that allow for a finite number of discontinuities at unknown locations.
The key idea is to replace the physical spacing $\delta_i = x_{i+1} - x_i$ between consecutive
observations with a latent spacing $\tau_i = \delta_i + \gamma_i$, where $\gamma_i \ge 0$ is a
non-negative *domain-warping* variable. A large $\gamma_i$ stretches the latent domain between
$x_i$ and $x_{i+1}$, compressing a smooth latent transition into a tiny physical interval and
thereby producing an apparent jump. The warping variables are constrained by a budget:
$\gamma_i \ge 0$ and $\sum_i \gamma_i \le M$.

The objective function to be minimised over $f$ and $\boldsymbol{\gamma}$ is

$$J(f, \boldsymbol{\gamma}) = \sum_{i=1}^n (y_i - f(x_i))^2 + \lambda \sum_{i=1}^{n-1} \int_0^{\tau_i} \left(\frac{d^2 f_i(s)}{ds^2}\right)^2 ds$$

subject to $\gamma_i \ge 0$ and $\sum_i \gamma_i \le M$, where the roughness penalty is
evaluated in the latent domain. For $\boldsymbol{\gamma} = \mathbf{0}$ the problem reduces
to the standard cubic smoothing spline.

Two estimation methods are provided.

**MLE via Kalman filter and smoother.** Following Wecker and Ansley (1983), the cubic
smoothing spline coincides with the signal extracted from a state-space model whose transition
matrix depends on the spacing $\delta_i$. Replacing $\delta_i$ with $\tau_i = \delta_i + \gamma_i$
yields a state-space model for the discontinuous case. The likelihood is maximised jointly over
$(\sigma, \lambda, \boldsymbol{\gamma})$ using constrained optimisation (CCSA/NLopt) with
analytical scores derived from the Kalman filter and smoother. The smoothed state provides
the estimated function. This approach is implemented in `ssj_mle()` (for fixed $\lambda$
and $M$) and `auto_ssj_mle()` (which selects $M$ automatically by EBIC).

**Alternating Time-Warping (ATW) algorithm.** As an alternative, a two-step block coordinate
descent algorithm with $O(n \log n)$ complexity per iteration is provided. The *f-step* solves
the banded linear system $(\mathbf{S} + \lambda \mathbf{R})\mathbf{z} = \mathbf{y}^*$ via a
block-tridiagonal Cholesky decomposition in $O(n)$ operations. The *$\gamma$-step* performs a
gradient descent on the roughness energy followed by a projection onto the simplex
$\{\boldsymbol{\gamma} \ge 0,\, \sum_i \gamma_i \le M\}$. The smoothing parameter $\lambda$
can be supplied or piloted from a preliminary `ssj_mle()` fit; the budget $M$ is selected
by EBIC over a grid. This approach is implemented in `ssj_atw()` (fixed $\lambda$ and $M$)
and `auto_ssj_atw()` (automatic selection of $M$).