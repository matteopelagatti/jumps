# cran-comments

## Submission summary

This is a major release of an existing CRAN package (previous version 1.0,
published 2025-03-24). It adds two new families of functions for smoothing
splines with discontinuities, a new data set, and fixes a number of bugs in the
existing code, several of which change numerical results; `NEWS.md` lists them
all and flags the ones that do.

The methods are described in two papers, both cited in `DESCRIPTION`:

* the Hodrick-Prescott filter with jumps in Maranzano and Pelagatti (2025),
  *Economic Modelling* 150, 107132, <doi:10.1016/j.econmod.2025.107132>;
* the smoothing splines with discontinuities in Pelagatti, Maranzano and
  Gatto Ronchieri (2026), <doi:10.48550/arXiv.2610.33137>.

## Test environments

* local Windows 11, R 4.5.2 (x86_64-w64-mingw32)
* win-builder, R-devel (x86_64-w64-mingw32)        -- TODO paste result
* win-builder, R-release (x86_64-w64-mingw32)      -- TODO paste result
* macOS builder, R-release (aarch64-apple-darwin20) -- TODO paste result

## R CMD check results

0 errors | 0 warnings | 3 notes

`CRAN incoming feasibility`, `R code for possible problems`, `tests`,
`re-building of vignette outputs` and the PDF manual all pass cleanly.

Notes seen locally and why they are not actionable:

* "Examples with CPU (user + system) or elapsed time > 5s: plot.ssj" -- the
  example draws the figure six times to document the plotting options. The
  model fit itself takes 0.02s; the time is ggplot2's namespace load (about
  3.7s on the check machine) plus five ggplot renders. We can reduce the
  number of renders if you would prefer the example under the threshold.
* "unable to verify current time" -- the check machine had no access to the
  time service used by that test.
* "Skipping checking math rendering: package 'V8' unavailable" -- V8 is not
  installed on the check machine.

## Reverse dependencies

None: `jumps` has no reverse depends, imports or suggests on CRAN, so no
other package is affected by the changes in numerical results.
