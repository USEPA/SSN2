This is an initial submission. It is an update to the `SSN` package,
which was archived alongside the retirement of `rgdal`, `rgeos`, and
`maptools`. Thank you.

-------

## Submission

This is an initial submission.

## R CMD check results

Here is the output from `devtools::check(manual = TRUE)` on
the Windows 10 x64 operating system

0 errors | 0 warnings | 2 notes

Note 1: `lsndata` subdirectory is 18.1 Mb. This subdirectory contains the
`.ssn` object required to run examples when written to the a user's temporary
directory. This file was included with `SSN`.

Note 2: Found a detritus in the temp directory. So examples that use the `lsn`
subdirectory write to a user's temporary directory. `SSN` also had examples
that were written to a user's temporary directory.

## Downstream dependencies

There are no downstream dependencies.
