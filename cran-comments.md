Thank you very much for all the time and effort the CRAN team puts into maintaining
packages and assuring their high quality.

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

-------

## Submission

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
