Thank you very much for all the time and effort the CRAN team puts into maintaining
packages and assuring their high quality.

## Revised CRAN Submission

* Size of tarball: 7201753 bytes A CRAN package should not be larger than 5 MB.
  Please reduce the size.
    * The CRAN package is now smaller than 5 MB (as indicated by `devtools::check()`)
    
* Please add `\value` to .Rd files regarding exported methods and explain the functions results in the documentation. Please write about the structure of the output (class) and also what the output means. (If a function does not return a value, please document that too, e.g.
`\value{No return value, called for side effects}` or similar) Missing Rd-tags: ssn_initial.Rd: `\value` ssn_update_path.Rd: `\value`
    * We have added `value` tags to `ssn_initial.Rd` and `ssn_update_path.Rd`
    
* `\dontrun{}` should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in `\dontrun{}` adds the comment ("# Not run:") as a warning for the user. Does not seem necessary.
Please replace `\dontrun{}` with `\dontrun{}`. Please unwrap the examples if they are executable in < 5 sec, or replace `dontrun{}` with `\donttest{}`.
    * We have removed the cases of `\dontrun{}`
    
* You write information messages to the console that cannot be easily suppressed.
It is more R like to generate objects that can be used to extract the information a user is interested in, and then `print()` that object.
Instead of `print()/cat()` rather use `message()/warning()` or
`if(verbose)cat(..) (or maybe stop())` if you really have to write text to the console. (except for print, summary, interactive functions) -> `R/createBinaryID.R`; `R/ssn_create_distmat.R`, ...
    * We have removed the use of `print()` and `cat()` (outside of print, summary, interactive functions) and replaced with `message()` or `warning()`

    
* Please make sure that you do not change the user's options, par or working directory. If you really have to do so within functions, please ensure with an *immediate* call of on.exit() that the settings are reset when the function is exited.
e.g.: -> `R/ssn_subset.R`; `R/ssn_write.R`; `R/createBinaryID.R` ...
`oldwd <- getwd() # code line i`
`on.exit(setwd(oldwd)) # code line i+1`
`...`
`setwd(...) # somewhere after`
`...`
e.g.:
If you're not familiar with the function, please check ?on.exit. This function makes it possible to restore options before exiting a function even if the function breaks. Therefore it needs to be called immediately after the option change within a function.
    * We have ensured that `on.exit(setwd(oldwd))` occurs immediately after the first call to `getwd()` within a function


-------

## Submission Notes

This is an initial submission.

## R CMD check results

Here is the output from a submission to win-builder checking the 
development, current, and old versions of R (accessed by running
`devtools::check_win_devel()`, `devtools::check_win_release()`,
and `devtools::check_win_oldrelease()`, respectively). 

0 errors | 0 warnings | 1 note

Note 1: `checking CRAN incoming feasibility ... NOTE`
  `Maintainer: 'Michael Dumelle <Dumelle.Michael@epa.gov>'`

## Downstream dependencies

There are no downstream dependencies. We have informed maintainers of downstream
dependendies of `SSN` that `SSN2` will be released shortly and have provided
them with the package structure.

-------

## Initial CRAN Submission Notes

This is an initial submission. We decide to archive `SSN`
alongside the retirement of `rgdal`, `rgeos`, and
`maptools` and create `SSN2`. We created a new package (rather than update 
`SSN`) for several reasons, the first and third of which are most important:

1. `SSN2` looks and behaves very differently from `SSN`. While making updates
to prepare for the retirement of `rgdal`, `rgeos`, and `maptools`, we decided
to go a step further, tearing down and totally rebuilding `SSN`. We changed the
underlying stream network objects from S4 to S3, added new model-fitting routines,
completely changed function and argument names, and added several other new 
features. `SSN2` also makes use of more R generics and incorporates tidyverse
infrastructure, which gives it a totally different, more modern feel than
`SSN` had.

2. The maintainers have changed. Michael Dumelle is now the maintainer, and
changing the package names makes this more explicit to users.

3. USEPA (Michael Dumelle's employer) requires technical review of research
projects before they are made public. They have reviewed and cleared `SSN2`
as a "new product", but have not cleared and reviewed any `SSN`-related material.
Keeping the name `SSN` should initiative another review request, which could 
take months, as the review will encompass all of `SSN`'s history.

4. We plan to publish new research articles using `SSN2` and believe the name
change is important so that users know which publications reference which 
`SSN` behaviors.

For these reasons, we would very much prefer to use the name `SSN2` and keep
`SSN` archived. If this is at all possible, we would really appreciate it.
Thank you again so much for your time.

