#' @export
# use generics to export tidy
generics::tidy

#' @export
# use generics to export glance
generics::glance

#' @export
# use generics to export augment
generics::augment

#' @export
spmodel::AICc

#' @export
spmodel::covmatrix

#' @export
spmodel::dispersion_initial

#' @export
spmodel::dispersion_params

#' @export
spmodel::glances

#' @export
spmodel::loocv

#' @export
spmodel::pseudoR2

#' @export
spmodel::randcov_initial

#' @export
spmodel::randcov_params

#' @export
spmodel::varcomp

# logit function
logit <- function(x) {
  if (x < 0 | x > 1) {
    stop("logit argument must be between zero and one", call. = FALSE)
  }
  log(x / (1 - x))
}

# expit function
expit <- function(x) {
  1 / (1 + exp(-x))
}

# remove type class
remove_covtype <- function(class_string) {
  sub("^[^_]*_", "", class_string)
}

# CRAN release questions
release_questions <- function() {
  c(
    "Have you changed version numbers in DESCRIPTION, CITATION, and README?",
    "Have you run pkgdown::build_site() and committed?"
  )
}
