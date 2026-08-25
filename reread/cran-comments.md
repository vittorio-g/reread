## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new release, so there is a "New submission" NOTE (checking CRAN
  incoming feasibility). This is expected for a first submission.
* The "checking HTML version of manual" NOTE ("no command 'tidy' found") is
  local to the build machine (HTML Tidy is not installed); it does not arise on
  systems where Tidy is available.

## Test environments

* Windows 11, R 4.6.0 (local), R CMD check --as-cran

## Downstream dependencies

There are currently no downstream dependencies.

## Notes

* The package is a faithful port of a companion browser tool; its scores are
  validated to reproduce that tool (deterministic components exactly, the
  permutation baseline to Monte Carlo error).
