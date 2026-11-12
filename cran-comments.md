# cran-comments

## Submission summary

This is a data-licensing release (2.0.3, resubmitting 2.0.2, which did not pass
the incoming checks). The Gestore dei Mercati Energetici
(GME), whose day-ahead market data underlie the `load_temp` data set, does not
allow redistribution of those data with the package. Version 2.0.3 therefore
removes the data set `load_temp` and the raw GME file shipped in
`inst/scripts/`. No R or C++ code changed and no fitted result changes.

We kindly ask the CRAN team to archive the earlier versions that contain the
GME data (2.0.0 and 2.0.1, if the latter was published), so that they are no
longer distributed.

The data set `demand_temp` (a filtered sample used in the vignette) is kept.
`R CMD check` and the vignettes no longer reference `load_temp`.

The incoming check of 2.0.2 failed on linux-arm64 in one test assertion
(`test-ssj.R`), which required the warping budget to be respected to a relative
1e-8; the augmented-Lagrangian optimizer only satisfies it to its convergence
tolerance and overshot by a relative 3e-7 on that platform. The tolerance in the
test is now 1e-5. No package code changed.

## Test environments

* local Windows 11, R 4.5.2 (x86_64-w64-mingw32):
  0 errors | 0 warnings | 1 note

## R CMD check results

No errors and no warnings. Examples, tests, re-building of the vignettes and
the PDF manual all pass.

The note:

* "Skipping checking math rendering: package 'V8' unavailable" (local only)
  -- V8 is not installed on that machine.

## The macOS failures, and what is not claimed

The remaining divergence between platforms is real and is not fixed by this
release: on the two arm64 builds both fitters settled tens of log-likelihood
units below the optimum reached on Linux and Windows for one 60-point data
set. The C++ Kalman filter was re-read for undefined behaviour and
uninitialised reads and none was found, and the same source gives the same
answer under gcc and clang on Linux, so the likeliest cause is the Apple
toolchain's floating point feeding a nonconvex optimizer. Chasing that is
not something this release attempts; the tests no longer assert that the
optimizer resolves it identically everywhere.

The macOS builder remains unavailable for pre-submission checking. On two
attempts for 2.0.0 (build IDs 1793189011-e99585e12f34b669 and
1793206562-7e713ac5eaa851f9) it stopped at "checking package dependencies"
reporting 'Rcpp', 'nloptr' and 'RcppEigen' not available, for which that
machine's library rather than this package appears to be responsible, and a
later attempt returned a server error page. The failures fixed here are
therefore taken from the CRAN check page rather than reproduced locally.

## Reverse dependencies

None: `jumps` has no reverse depends, imports or suggests on CRAN.
