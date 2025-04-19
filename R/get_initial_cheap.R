get_initial_cheap <- function(params_object) {
  list(
    tailup_initial = params_object$tailup,
    taildown_initial = params_object$taildown,
    euclid_initial = params_object$euclid
  )
}
