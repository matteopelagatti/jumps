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