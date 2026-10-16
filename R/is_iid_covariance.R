is_iid_covariance <- function(tailup, taildown, euclid, random = NULL) {
  is.null(random) && inherits(tailup, "tailup_none") &&
    inherits(taildown, "taildown_none") && inherits(euclid, "euclid_none")
}
