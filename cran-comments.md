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

The macOS builder could not check the package, on two attempts a day apart
(build IDs 1793189011-e99585e12f34b669 and 1793206562-7e713ac5eaa851f9,
profile r-release-macosx-arm64, R 4.6.1). Both stopped within fifteen
seconds at "checking package dependencies" reporting

    Packages required but not available: 'Rcpp', 'nloptr', 'RcppEigen'
    VignetteBuilder package required for checking but not installed: 'knitr'

with the remaining suggested packages -- 'rmarkdown', 'testthat', 'ggplot2',
'xts' -- also reported unavailable. All eight are plain CRAN packages, all
are present locally and on win-builder, and that same check step passes in
both of those environments, so the build machine's library appears to be
unpopulated rather than the package having a dependency problem. No macOS
result is therefore claimed here.

## Reverse dependencies

None: `jumps` has no reverse depends, imports or suggests on CRAN, so no
other package is affected by the changes in numerical results.
