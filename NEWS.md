# jumps 2.0.0

This is a major release: the package now covers smoothing splines with
discontinuities alongside the original Hodrick-Prescott filter with jumps, which is
why the version jumps from 1.0.0 to 2.0.0.

## New features

* New family of functions for cubic smoothing splines with discontinuities
  (CSSD), implementing the dynamic-programming algorithm of Storath and
  Weinmann (2024): `ssj_cssd_potts()`, `auto_ssj_cssd()` (K-fold
  cross-validation of `p` and `gamma`) and `predict_cssd_potts()`.
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

* The default `ebic_xi` of `ssj_atw()` and `auto_ssj_atw()` changed from 1 to
  0.3. The alternating algorithm evaluates the log-likelihood at the unbiased
  residual variance rather than maximising the exact Gaussian likelihood, so it
  is not on the same scale as `ssj_mle()`/`ssj_kfs()` and needs a weaker
  extended-BIC search penalty to reach comparable sparsity. On the simulated
  signal of the paper the new default recovers exactly the three true
  discontinuities, where the old one selected none.

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

## Bug fixes

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

* `predict_cssd_potts()` returned its values in the order of `sort(xnew)`
  rather than in the order of `xnew`, silently mis-aligning predictions when
  the query points were not sorted.

* `predict_cssd_potts()` returned exactly 0 for query points falling between
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
