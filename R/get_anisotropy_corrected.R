#' Correct anisotropy argument if initial values specified
#'
#' @param anisotropy Anisotropy argument
#' @param initial_object Initial object
#'
#' @return A corrected anisotropy argument. Values specified in the initial object
#'   take precedence over the anisotropy argument.
#' @noRd
get_anisotropy_corrected <- function(anisotropy, initial_object) {
  if (inherits(initial_object$euclid_initial, "euclid_none")) {
    anisotropy <- FALSE
  } else {
    if (anisotropy) {
      if (all(c("rotate", "scale") %in% names(initial_object$euclid_initial$initial))) {
        is_rotate_zero <- (!is.na(initial_object$euclid_initial$initial[["rotate"]])) && initial_object$euclid_initial$initial[["rotate"]] == 0
        is_rotate_known <- initial_object$euclid_initial$is_known[["rotate"]]
        is_scale_one <- (!is.na(initial_object$euclid_initial$initial[["scale"]])) && initial_object$euclid_initial$initial[["scale"]] == 1
        is_scale_known <- initial_object$euclid_initial$is_known[["scale"]]
        if (is_rotate_zero && is_rotate_known && is_scale_one && is_scale_known) {
          anisotropy <- FALSE
        }
      }
    } else {
      if ("rotate" %in% names(initial_object$euclid_initial$initial)) {
        is_rotate_zero <- (!is.na(initial_object$euclid_initial$initial[["rotate"]])) && initial_object$euclid_initial$initial[["rotate"]] == 0
        is_rotate_known <- initial_object$euclid_initial$is_known[["rotate"]]
        if (!is_rotate_zero || !is_rotate_known) {
          anisotropy <- TRUE
        }
      }

      if ("scale" %in% names(initial_object$euclid_initial$initial)) {
        is_scale_one <- (!is.na(initial_object$euclid_initial$initial[["scale"]])) && initial_object$euclid_initial$initial[["scale"]] == 1
        is_scale_known <- initial_object$euclid_initial$is_known[["scale"]]
        if (!is_scale_one || !is_scale_known) {
          anisotropy <- TRUE
        }
      }
    }
  }
  anisotropy
}
