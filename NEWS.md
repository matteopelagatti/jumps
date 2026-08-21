# jumps 2.0.0

This is a major release: the package now covers smoothing splines with
discontinuities alongside the original Hodrick-Prescott filter with jumps, which is
why the version jumps from 1.0.0 to 2.0.0.

## New features

* New family of functions for cubic smoothing splines with discontinuities
  (CSSD), implementing the dynamic-programming algorithm of Storath and
  Weinmann (2024): `ssj_cssd_dp()`, `auto_ssj_cssd()` (K-fold
  cross-validation of `p` and `gamma`) and `predict_cssd_dp()`.
  `ssj_cssd()` provides a fast piecewise-constant (total-variation)
  approximation for larger problems.

* New family of smoothing splines with jumps based on time warping:
  `ssj_mle()`/`auto_ssj_mle()` (smoothing parameter estimated by maximum
  likelihood), `ssj_kfs()`/`auto_ssj_kfs()` (fixed smoothing parameter) and
  `ssj_atw()`/`auto_ssj_atw()` (fast alternating time-warping algorithm).
  All of them support the extended BIC (`ebic`) for selecting the jump budget.

* The `ssj_*` functions now return an object of class `ssj`, with `print` and
  `plot` methods, mirroring what `hpj()` already offered. The object is still a
  plain list underneath, so existing code that reaches into its slots keeps
  working; it gains a `method` code, the response `y` as the estimator saw it,
  and the `call`. `plot()` reproduces Figure 4 of the paper: the fitted curve
  over the data, vertical lines at the declared discontinuities and a lower
  panel with the estimated domain warps; pass `prob` for a confidence band.
  No title is drawn unless one is passed through `main`.

* New `predict()` method for `ssj` objects. For the domain-warping estimators
  it evaluates the cubic the model actually implies -- the state carries both
  the level and the slope at every sampling site, so the fit is a cubic Hermite
  in the latent coordinate -- instead of interpolating the fitted values. It was
  validated against the exact model-implied interpolation, obtained by inserting
  the query points as missing observations and re-running the Kalman smoother:
  the two agree to 1e-12. Beyond the range of the data the fit is extended
  linearly. CSSD objects are passed to `predict_cssd_dp()`, and the
  piecewise-constant TV fit gets a step interpolation.

* The `ssj_*` fitters gained a `smoothed_slope` slot, the derivative of the fit
  with respect to the latent (warped) coordinate, which is what `predict()`
  uses.

* `ssj_mle()`, `ssj_kfs()`, `ssj_fix()` and `ssj_atw()`, and their `auto_*`
  wrappers, gained a `weights` argument. The observation variance becomes
  `lambda*sigma^2/w_i` (for the alternating algorithm the data term becomes
  `sum w_i (y_i-f_i)^2`), so a point with weight `w_i` counts as the mean of
  `w_i` independent observations. This makes it correct to collapse tied values
  of `x` to their mean and pass the replicate counts as weights, which is the
  right way to handle a covariate with many repeated values: verified to
  reproduce the fit on the raw replicates exactly for `ssj_atw()` and to 1e-9
  and 1e-5 for `ssj_fix()` and `ssj_kfs()`. Unit weights reproduce the previous
  results exactly, so nothing existing changes. See the caveat in `?ssj_mle`:
  because that function estimates the noise level, collapsing loses the
  within-group scatter that identifies it.

* `ssj_mle()`, `ssj_kfs()` and `ssj_fix()` now collapse tied values of `x` to
  their weighted mean automatically, as `ssj_atw()` already did, turning the
  replicate counts into observation weights and -- for `ssj_mle()`, which
  estimates the noise level -- feeding the within-group scatter back into the
  likelihood. Previously they fitted the raw replicated data, where the
  intervals between replicates have zero length and a warp placed on one of
  them is a discontinuity between two observations at the same `x`, which is
  not identified. On the electricity-demand application that made the fit
  depend on the order of the rows, with effective degrees of freedom anywhere
  between 8 and 79. The `x`, `smoothed_level` and `gamma` slots are now
  returned on the distinct sites, as they already were for `ssj_atw()`.

* `llt()` and `llt_delta()` returned `NA` for the whole likelihood when the
  FIRST observation was missing: the `t = 0` initialisation lacked the
  missing-value guard that the loop applies to every later observation. Any
  series starting with `NA` was affected.

* `ssj_mle()` and `auto_ssj_mle()` gained `ss_within`, the pooled within-group
  sum of squares discarded when replicates are collapsed to their means.
  Supplying it alongside the replicate counts in `weights` restores the term of
  the likelihood that collapsing removes, so the noise level -- and hence
  `lambda` -- is estimated as if the raw replicates had been supplied: on a test
  case the estimated observation variance goes from 0.086 back to 0.021 against
  a raw-data value of 0.021, at three different noise levels and to four or five
  significant digits. It changes nothing at a fixed smoothing parameter, and
  nothing when omitted.

* `ssj_atw()` already collapsed tied values of `x`, but gave every collapsed
  point unit weight, silently treating a mean of six observations as no more
  informative than a single one. It now carries the replicate counts as
  weights, which changes the effective degrees of freedom on data with ties
  (3.52 rather than 2.84 on a test case).

* New function `ssj_breaks()`, which applies the detection rule of the paper --
  a warp counts as a discontinuity when it clears both 5% of the total warping
  and 30% of the mean sampling interval -- and returns the declared breaks with
  their location, warp and jump size. Breaks are dated at the first observation
  of the new regime. An absolute threshold can be given through `thresh`, as in
  Section 4.2 of the paper.

* New dataset `load_temp` with the hourly electricity load of Northern Italy
  and the temperature of Milan, and `demand_temp`, the filtered analysis
  sample used in the demand-versus-temperature application. Both are shipped
  as `.rda`, so neither the 11 MB raw market file nor the Open-Meteo API call
  is needed to reproduce the analysis. The script that derives `demand_temp`
  is in `data-raw/demand_temp.R`.

* The introduction vignette gained two sections: one on the smoothing splines
  with jumps, and one working through the electricity demand versus
  temperature application with base graphics only.

* **The `ebic_xi` defaults are now calibrated rather than taken from theory**,
  and they differ by estimator: 0.3 for `ssj_mle()`, 0.15 for `ssj_kfs()` and
  0.25 for `ssj_atw()`, against 1, 1 and 0.3 before. Each is the value that
  minimises worst-case regret in root mean squared error over a factorial
  design of 48 cells -- four sample sizes, three noise levels and four counts
  of discontinuities, 17,400 replications -- shipped as
  `inst/scripts/xi_calibration.R`.

  Fits selected by `auto_ssj_mle()`, `auto_ssj_kfs()` and `auto_ssj_atw()` will
  therefore change: the smaller penalties admit more discontinuities. Pass
  `ebic_xi` explicitly to reproduce earlier results.

  The three differ for real reasons rather than by accident. `ssj_kfs()` sits
  below `ssj_mle()` because `auto_ssj_kfs()` holds `lambda` fixed across the
  budget grid while `auto_ssj_mle()` re-estimates it, and being able to re-tune
  the smoothness gives more scope to overfit. For `ssj_atw()` the optimum falls
  steadily with the sample size -- roughly 0.32, 0.22, 0.15, 0.12 at
  `n` = 50, 100, 200, 400 -- while the total charge per declared discontinuity
  stays flat at about 3.9, which yields the rule
  `xi(n) = (3.93 - 0.419*log(n)) / (2*log(n-1))`. That rule is documented and
  is worth using when `n` is far from a few hundred; the constant is the
  default because it behaves acceptably everywhere and a default that changes
  with the data is easy to overlook.

  All three sit below the Chen and Chen (2008) consistency floor of about 0.5,
  which is expected rather than a defect: `df` is the trace of the smoothing
  matrix and grows with each activated warp, so the plain BIC term already pays
  part of the search cost -- about 0.08, 0.17 and 0.21 on the `xi` scale for
  the three estimators.

* `plot()` for `ssj` objects gains `use_ggplot` and **defaults to it**, as
  `plot.hpj()` already did. The figure is built with ggplot2 and *returned*
  rather than drawn, so it can be extended, themed or saved with `ggsave()`;
  the fitted curve and the warp profile become facets of a single object rather
  than two separate plots, so the x axes stay aligned and one value comes back.
  Pass `use_ggplot = FALSE` for base graphics, which draws directly and returns
  `NULL`.

  One consequence worth knowing: the plot is now a value, so inside a loop or a
  script it has to be `print`ed to appear, exactly as any ggplot2 figure does.

* New `auto_ssj_atw_gcv()` and `auto_ssj_atw_cv()`, which choose the smoothing
  constant for the alternating algorithm without recourse to `ssj_mle()`, then
  hand it to `auto_ssj_atw()` for the warping budget as before. The first
  minimises generalized cross-validation over `lambda` at zero budget, at one
  fit per candidate; the second cross-validates over `lambda` and the budget
  jointly and keeps only `lambda`. Both accept `weights`, warn when the
  selected value is an endpoint of the grid, and return the search path.

  The two are not equivalent. Cross-validation is the reliable one: the
  `lambda` it selects tracks the maximum-likelihood value closely and its
  boundary bound in 0 to 23% of samples on the calibration design. Generalized
  cross-validation is around two orders of magnitude cheaper but degenerates
  exactly where the method is aimed -- with three clear jumps and
  `sigma = 0.05` its floor bound in 83% of samples at `n = 100` and in every
  sample at `n = 200`, because at zero budget the spline can only represent a
  discontinuity by bending, so the criterion improves without limit as `lambda`
  falls. Accuracy is nevertheless within a few percent of the alternatives,
  since the floor is a serviceable value; but it is a floor, not a minimum.

  The default `lambda` grid is built in the normalised coordinates
  `ssj_atw()` works in and mapped back, so it follows the scale of `x` -- a
  grid fixed in absolute terms would suit only one scaling, `lambda` having the
  units of `x^3`.

* `ssj_atw()` gains a `pars` slot, holding the unbiased residual standard
  deviation and the `lambda` it was given, so the smoothing constant can be
  recovered from a fit as it already could from the other estimators.

## Bug fixes

* **`ssj_kfs()` and `ssj_fix()` frequently failed to converge**, returning a
  scale parameter many orders of magnitude too large and a log-likelihood far
  below the attainable one. Both started the optimiser at a fixed `sdy/10`,
  which bears no relation to the optimum: when `lambda` is small the optimal
  `sigma` is large -- 37 against a start of 0.1 on the test case -- and the
  optimiser diverged. Against `ssj_mle()` run on the same data at the same
  `lambda`, where the two fit exactly the same model and must agree,
  `ssj_kfs()` reached a worse optimum in **37 of 40** test replications, with
  gaps of 800 log-likelihood units; `ssj_fix()` was worse still, at around
  7,000,000 units. `auto_ssj_kfs()` propagated the failure along the whole
  grid through its warm start.

  Both now start from the profile solution for `sigma`, which needs a single
  filter pass: with `lambda` fixed, `sigma` is a pure scale, so the profile
  log-likelihood is maximised at `mean(v^2/F)` evaluated at `sigma = 1`. That
  lands essentially on the optimum, and it is followed by a short multi-start
  in the same spirit as the grid of `lambda` starts already used by
  `ssj_mle()`. `ssj_kfs()` now matches `ssj_mle()` to 1e-10 at the same
  `lambda`, and no longer falls short in any of the 40 replications.

  Users of `ssj_kfs()`, `auto_ssj_kfs()`, `ssj_fix()` or `auto_ssj_fix()`
  should re-run anything fitted with an earlier version.

* **`ssj_mle()`, `ssj_kfs()` and `ssj_atw()` reported the log-likelihood of an
  internally rescaled copy of `y`, not of the data as supplied**, and with it
  the absolute values of `aic`, `aicc`, `bic`, `hq` and `ebic`. All three
  divide `y` by its standard deviation before fitting, and the Jacobian of that
  map was never added back, so the reported value sat below the true one by
  `nobs * log(sd(y))` -- 12.5 units on a 60-point test case. The symptom was
  that multiplying `y` by 1000 left the reported log-likelihood unchanged,
  where it should fall by `nobs * log(1000)`.

  All three now report the likelihood on the scale of the data supplied.
  `ssj_mle()` uses `nobs + df_within` factors when the within-group scatter has
  been restored, that being the number of raw observations. `ssj_fix()` and the
  `hpfj*()` family never rescaled `y` and are unchanged.

  **This does not change any model selection.** The offset is identical at
  every `maxsum` for a given data set, so it cancels in the comparisons the
  `auto_*` functions make; the budget they choose, the fitted curve and `df`
  are all unaffected. What changes is the reported `loglik` and the absolute
  level of the criteria, which are now comparable across estimators: on the
  identical fitted model `ssj_kfs()` and `ssj_fix()` now report a `bic` of
  41.088 against 41.060, where before it was 66.157 against 41.060.

  A caveat that remains, and is not a defect of the above: these fitters
  approximate the diffuse initial conditions with a large proper prior
  (`1e5` times the data variance), and the reported likelihood depends on that
  arbitrary multiplier -- by about `0.5*log(10)` per decade per diffuse state.
  Absolute criterion values therefore carry an additive constant fixed by that
  convention. It is common to all of `ssj_mle()`, `ssj_kfs()` and `ssj_fix()`,
  but `ssj_fix()` does not rescale `x`, so its slope prior differs slightly
  from the others'; that is the residual 0.03 above. Exact diffuse
  initialisation would remove the constant altogether.

* **The effective degrees of freedom were wrong for every Kalman-based
  estimator** (`ssj_mle()`, `ssj_kfs()`, `ssj_fix()`, `hpfj()`, `hpfj_fix()`,
  `hpfjx()`), and with them `aic`, `aicc`, `bic`, `hq` and `ebic`, and hence the
  jump budget selected by the `auto_*` functions. `df` is the trace of the
  smoothing matrix, obtained by summing the observation weights computed inside
  `llt()` and `llt_delta()` with the recursion of Koopman and Harvey (2003).
  Those weights carried the wrong time index on the backward variance smoother:
  `N` was taken one step too early. Since `N` grows going backwards, the
  subtracted terms overshot and the weights could fall well outside `(0, 1]` and
  sum to a negative number -- reported `df` values such as -579 on a smooth fit
  of 60 observations. The weights are now exactly the diagonal of the hat
  matrix, verified against a central-difference computation of it for equally
  and unequally spaced designs, with and without observation weights, and with
  missing values (which correctly receive weight zero). `hpj` objects return
  the corrected weights in their `weights` slot. `ssj_atw()` was never
  affected: its `df` comes from the trace of `A^-1 S` in the alternating
  algorithm.

  Note for users reproducing published results: the formula in Appendix B.5 of
  Maranzano and Pelagatti (2025) carries the same index error, and information
  criteria computed with earlier versions of this package will differ.

* `hpfj_fix()` reported the wrong value for `pars["gamma"]`: it squared the
  first additional standard deviation instead of the estimated `sqrt(gamma)`.

* `plot.hpj()` drew an inverted confidence band when `use_ggplot = TRUE` and
  `prob` was given (`ymin` and `ymax` were swapped, because the negative
  quantile from `qnorm()` was not taken in absolute value).

* `plot.hpj()` now falls back to base graphics, with a warning, when
  **ggplot2** is not installed, and only uses **xts** when that package is
  available.

* `BIC.hpj()` failed with "argument to 'which' is not logical" when one of the
  objects passed through `...` had more than one class. It now uses
  `inherits()`.

* `predict_cssd_dp()` returned its values in the order of `sort(xnew)`
  rather than in the order of `xnew`, silently mis-aligning predictions when
  the query points were not sorted.

* `predict_cssd_dp()` returned exactly 0 for query points falling between
  the last data site of one segment and the first data site of the next, that
  is, precisely around each detected discontinuity. Segment membership is now
  determined by cutting at the midpoint of each gap.

* `auto_ssj_cssd()` could not be called without an explicit `x` argument,
  although its body handled `x = NULL`; `x` now defaults to `NULL`. It also
  sorts and de-duplicates `x` before cross-validation, as the other CSSD
  functions do.

* `ssj_atw()` produced an all-`NaN` fit when `y` contained missing values,
  which its documentation claimed to allow. It now raises an informative error
  pointing to `ssj_kfs()`/`ssj_mle()`.

* `hpfj_fix()` used `var(y)` without `na.rm = TRUE` and initialised the filter
  at `y[1]`, both of which failed on series with missing values.

* The examples of `ssj_atw()` and `auto_ssj_atw()` failed with "'x' and 'y'
  lengths differ", because those functions collapse tied `x` values. All the
  spline examples now plot against the `x` slot of the returned object.

## Documentation and packaging

* The `NAMESPACE` file is now generated by **roxygen2**, so that exports and
  imports can no longer drift apart from the source (`hpfj_fix()` was exported
  by hand but had no `@export` tag).

* Argument and return-value documentation has been made uniform across all
  functions, and several mismatches with the code have been fixed: the
  `maxit`/`maxit_nm` argument of `auto_ssj_cssd()`, the `std_devs` slot of
  `hpj()` objects (actually named `sigmas`), the documented default of the
  `ic` argument of `auto_ssj_atw()`, and the dangling link to a
  non-existent `cssd_cv_lambda()`.

* `predict_cssd_dp()` extrapolated beyond the end knots of a segment with
  that segment's end cubic. The reference implementation of Storath and Weinmann
  extends each segment linearly from its end knots (`linext_pp`), which is the
  convention CSSD is defined with; a cubic extrapolation diverges quickly across
  the gaps around a discontinuity. It now extends linearly.

* `predict()` on a domain-warping fit extrapolated beyond the data using the
  physical slope just inside the end interval, that is the latent slope times
  that interval's dilation `tau/delta`. The warps live strictly between
  observations, so outside the sample the dilation is 1. On a boundary interval
  carrying a warp the factor can be several tens -- one simulated replication
  had `tau/delta = 40`, turning a slope of 6.6 into 267 and predicting -11.7
  where the signal was 0. It now extrapolates with the bare latent slope.

* The documentation of `ssj_cssd()` no longer claims that the difference
  operator is weighted by `1/dx`; the implementation deliberately uses the
  unweighted operator.

* `hpj()` now raises an explicit error when `xreg` is supplied, instead of
  documenting an argument that was not implemented.

* The default `maxsum` of `ssj_mle()` was `sd(y)/mean(diff(x))`, a formula
  inherited from `ssj_fix()`, where the budget constrains disturbance standard
  deviations. For the domain-warping estimators the budget is expressed in the
  units of `x`, so the old default was on the wrong scale: on `faithful` it
  evaluated to about 1053 against an x-range of 3.5, leaving the warping
  effectively unconstrained. It is now `max(x) - min(x)`, matching both
  `ssj_kfs()` and the upper end of the grid searched by `auto_ssj_mle()`.
  Calls that pass `maxsum` explicitly are unaffected.

* `ssj_mle()`, `ssj_kfs()` and `ssj_fix()` gave different results depending on
  the order in which the observations were supplied, even when the covariate
  had no ties. They computed `var(diff(y))` on the *unsorted* response and used
  it to scale the diffuse prior on the slope, so a permutation of the rows
  changed the starting point of the Kalman filter (by a factor of 6.5 on the
  `demand_temp` data). It is now computed after sorting by `x`. Results for
  input that was already sorted by `x` are unchanged.

* `ssj_kfs()` silently ignored its own default `maxsum`. The default,
  `max(x) - min(x)`, refers to `x`, and R forces such a default lazily, at its
  first use; that happened only after `x` had been rescaled to `[0, 1]`, so the
  budget evaluated to 1 instead of the range of `x`. The defaults of
  `ssj_kfs()`, `ssj_mle()` and `ssj_atw()` are now forced while `x` is still on
  its original scale. Calls that pass `maxsum` explicitly were never affected.

* The examples of `ssj_atw()`, `ssj_kfs()` and `ssj_mle()` used `maxsum = 150`
  on the `faithful` data. For the domain-warping estimators `maxsum` is a
  budget in the units of `x`, whose range there is about 3.5, so the budget was
  some forty times the whole domain and the resulting fits were degenerate (the
  Kalman variant even reported negative effective degrees of freedom). The
  examples now use a budget of 1.

* The reference for the HP filter with jumps has been updated to the published
  article: Maranzano and Pelagatti (2025), *Economic Modelling* 150, 107132,
  <doi:10.1016/j.econmod.2025.107132>, replacing the FEEM working paper in
  `DESCRIPTION`, `inst/CITATION`, the package help page, the README and the
  introduction vignette.

* Package-level documentation is now generated from `R/jumps-package.R`, and
  the reference papers, replication data and rendered vignettes are no longer
  shipped inside the source tarball.
