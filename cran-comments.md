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

* local Windows 11, R 4.5.2 (x86_64-w64-mingw32):
  0 errors | 0 warnings | 2 notes
* win-builder, R 4.6.1 release (x86_64-w64-mingw32):
  0 errors | 0 warnings | 1 note
* win-builder, R-devel (x86_64-w64-mingw32)         -- TODO paste result
* macOS builder, R 4.6.1 Patched (aarch64-apple-darwin23):
  the check could not run; see 'macOS builder' below

## R CMD check results

No errors and no warnings on either platform. Examples, tests, re-building
of the vignettes and the PDF manual all pass; so do `CRAN incoming
feasibility` and `R code for possible problems`.

The notes:

* "Possibly misspelled words in DESCRIPTION: Gatto, Ronchieri, Storath,
  Weinmann" (win-builder). All four are surnames: Giulio Cesare Gatto
  Ronchieri, a contributor and co-author of the preprint cited above, and
  Martin Storath and Andreas Weinmann, whose algorithm the CSSD functions
  implement.
* "unable to verify current time" (local only) -- that machine had no access
  to the time service used by the test.
* "Skipping checking math rendering: package 'V8' unavailable" (local only)
  -- V8 is not installed on that machine.

## macOS builder

The macOS builder could not check the package. It stopped after six seconds
at "checking package dependencies", reporting that none of the dependencies
-- 'Rcpp', 'nloptr', 'RcppEigen', 'knitr', 'rmarkdown', 'testthat',
'ggplot2', 'xts' -- was installed on the build machine. All eight are plain
CRAN packages, all are present locally and on win-builder, and the same
check step passes in both of those environments, so this appears to be an
unpopulated library on that service rather than a problem with the package.
We will retry it.

## Reverse dependencies

None: `jumps` has no reverse depends, imports or suggests on CRAN, so no
other package is affected by the changes in numerical results.
