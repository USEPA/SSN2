#' Print values
#'
#' @description Print fitted model objects and summaries.
#'
#' @param x A fitted model object from [ssn_lm()], a fitted model object from [ssn_glm()],
#'   or output from \code{summary(x)} or or \code{anova(x)}.
#' @param digits The number of significant digits to use when printing.
#' @param signif.stars Logical. If \code{TRUE}, significance stars are printed for each coefficient
#' @param ... Other arguments passed to or from other methods.
#'
#' @return Printed fitted model objects and summaries with formatting.
#'
#' @name print.SSN2
#' @method print ssn_lm
#' @export
print.ssn_lm <- function(x, digits = max(3L, getOption("digits") - 3L),
                         ...) {
  cat("\nCall:\n", paste(deparse(x$call),
    sep = "\n",
    collapse = "\n"
  ), "\n\n", sep = "")

  cat("\n")

  cat("Coefficients (fixed):\n")
  print.default(format(coef(x, type = "fixed"), digits = digits),
    print.gap = 2L,
    quote = FALSE
  )

  cat("\n")

  cat("Coefficients (covariance):\n")
  x_tailup <- coef(x, type = "tailup")
  x_tailup <- data.frame(
    Effect = class(x_tailup),
    Parameter = names(x_tailup),
    Estimate = as.vector(x_tailup),
    is_known = x$is_known$tailup
  )
  x_taildown <- coef(x, type = "taildown")
  x_taildown <- data.frame(
    Effect = class(x_taildown),
    Parameter = names(x_taildown),
    Estimate = as.vector(x_taildown),
    is_known = x$is_known$taildown
  )
  x_euclid <- coef(x, type = "euclid")
  x_euclid <- data.frame(
    Effect = class(x_euclid),
    Parameter = names(x_euclid),
    Estimate = as.vector(x_euclid),
    is_known = x$is_known$euclid
  )
  x_nugget <- coef(x, type = "nugget")
  x_nugget <- data.frame(
    Effect = "nugget",
    Parameter = names(x_nugget),
    Estimate = as.vector(x_nugget),
    is_known = x$is_known$nugget
  )
  x_cov <- rbind(x_tailup, x_taildown, x_euclid, x_nugget)
  x_cov$Effect <- gsub("_", " ", x_cov$Effect)
  x_cov$Parameter <- gsub("de", "de (parsill)", x_cov$Parameter)
  logi1 <- x_cov$Estimate == 0 | x_cov$Estimate == Inf | (x_cov$Estimate == 1 & x_cov$Parameter == "scale")
  logi2 <- x_cov$is_known
  x_cov <- x_cov[!(logi1 & logi2), , drop = FALSE]
  x_cov <- x_cov[, -which(names(x_cov) == "is_known"), drop = FALSE]


  if (!is.null(x$random)) {
    x_rand <- coef(x, type = "randcov")
    x_rand <- data.frame(
      Effect = "random",
      Parameter = names(x_rand),
      Estimate = as.vector(x_rand),
      is_known = x$is_known$randcov
    )
    logi1 <- x_rand$Estimate == 0
    logi2 <- x_rand$is_known
    x_rand <- x_rand[!(logi1 & logi2), , drop = FALSE]
    x_rand <- x_rand[, -which(names(x_rand) == "is_known"), drop = FALSE]
    x_cov <- rbind(x_cov, x_rand)
  }
  print.data.frame(format(x_cov, digits = digits),
    print.gap = 2L,
    quote = FALSE,
    row.names = FALSE
  )

  cat("\n")

  #
  # tailup_coef <- coef(x, type = "tailup")
  # cat(paste("\nCoefficients (", gsub("tailup_", "", class(tailup_coef)), " tailup covariance):\n", sep = ""))
  # print.default(format(tailup_coef, digits = digits),
  #               print.gap = 2L,
  #               quote = FALSE
  # )
  #
  # taildown_coef <- coef(x, type = "taildown")
  # cat(paste("\nCoefficients (", gsub("taildown_", "", class(taildown_coef)), " taildown covariance):\n", sep = ""))
  # print.default(format(taildown_coef, digits = digits),
  #               print.gap = 2L,
  #               quote = FALSE
  # )
  #
  # euclid_coef <- coef(x, type = "euclid")
  # if (!x$anisotropy) {
  #   euclid_coef <- euclid_coef[-which(names(euclid_coef) %in% c("rotate", "scale"))]
  # } # class gets dropped here
  # cat(paste("\nCoefficients (", gsub("euclid_", "", class(coef(x, type = "euclid"))), " Euclidean covariance):\n", sep = ""))
  # print.default(format(euclid_coef, digits = digits),
  #               print.gap = 2L,
  #               quote = FALSE
  # )
  #
  # nugget_coef <- coef(x, type = "nugget")
  # cat(paste("\nCoefficients (", "nugget covariance):\n", sep = ""))
  # print.default(format(nugget_coef, digits = digits),
  #               print.gap = 2L,
  #               quote = FALSE
  # )
  #
  # # cat("\n")
  #
  # if (length(coef(x, type = "randcov"))) {
  #   cat("Coefficients (random effects):\n")
  #   print.default(format(coef(x, type = "randcov"), digits = digits),
  #                 print.gap = 2L,
  #                 quote = FALSE
  #   )
  #
  #   cat("\n")
  # }

  invisible(x)
}

#' @rdname print.SSN2
#' @method print ssn_glm
#' @export
print.ssn_glm <- function(x, digits = max(3L, getOption("digits") - 3L),
                          ...) {
  cat("\nCall:\n", paste(deparse(x$call),
    sep = "\n",
    collapse = "\n"
  ), "\n\n", sep = "")

  cat("\n")

  cat("Coefficients (fixed):\n")
  print.default(format(coef(x, type = "fixed"), digits = digits),
    print.gap = 2L,
    quote = FALSE
  )

  cat("\n")

  cat("Coefficients (covariance):\n")
  x_tailup <- coef(x, type = "tailup")
  x_tailup <- data.frame(
    Effect = class(x_tailup),
    Parameter = names(x_tailup),
    Estimate = as.vector(x_tailup),
    is_known = x$is_known$tailup
  )
  x_taildown <- coef(x, type = "taildown")
  x_taildown <- data.frame(
    Effect = class(x_taildown),
    Parameter = names(x_taildown),
    Estimate = as.vector(x_taildown),
    is_known = x$is_known$taildown
  )
  x_euclid <- coef(x, type = "euclid")
  x_euclid <- data.frame(
    Effect = class(x_euclid),
    Parameter = names(x_euclid),
    Estimate = as.vector(x_euclid),
    is_known = x$is_known$euclid
  )
  x_nugget <- coef(x, type = "nugget")
  x_nugget <- data.frame(
    Effect = "nugget",
    Parameter = names(x_nugget),
    Estimate = as.vector(x_nugget),
    is_known = x$is_known$nugget
  )
  x_dispersion <- coef(x, type = "dispersion")
  x_dispersion <- data.frame(
    Effect = "dispersion",
    Parameter = names(x_dispersion),
    Estimate = as.vector(x_dispersion),
    is_known = x$is_known$dispersion
  )
  x_cov <- rbind(x_tailup, x_taildown, x_euclid, x_nugget, x_dispersion)
  x_cov$Effect <- gsub("_", " ", x_cov$Effect)
  x_cov$Parameter <- gsub("de", "de (parsill)", x_cov$Parameter)
  logi1 <- x_cov$Estimate == 0 | x_cov$Estimate == Inf | (x_cov$Estimate == 1 & x_cov$Parameter == "scale")
  logi2 <- x_cov$is_known
  x_cov <- x_cov[!(logi1 & logi2), , drop = FALSE]
  x_cov <- x_cov[, -which(names(x_cov) == "is_known"), drop = FALSE]


  if (!is.null(x$random)) {
    x_rand <- coef(x, type = "randcov")
    x_rand <- data.frame(
      Effect = "random",
      Parameter = names(x_rand),
      Estimate = as.vector(x_rand),
      is_known = x$is_known$randcov
    )
    logi1 <- x_rand$Estimate == 0
    logi2 <- x_rand$is_known
    x_rand <- x_rand[!(logi1 & logi2), , drop = FALSE]
    x_rand <- x_rand[, -which(names(x_rand) == "is_known"), drop = FALSE]
    x_cov <- rbind(x_cov, x_rand)
  }
  print.data.frame(format(x_cov, digits = digits),
    print.gap = 2L,
    quote = FALSE,
    row.names = FALSE
  )

  cat("\n")

  invisible(x)
}

#' @rdname print.SSN2
#' @method print summary.ssn_lm
#' @export
print.summary.ssn_lm <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 signif.stars = getOption("show.signif.stars"),
                                 ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  cat("\nResiduals:\n")
  resQ <- c(
    min(x$residuals$response), quantile(x$residuals$response, p = c(0.25, 0.5, 0.75), na.rm = TRUE),
    max(x$residuals$response)
  )
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)

  # pasting the fixed coefficient summary
  cat("\nCoefficients (fixed):\n")
  coefs_fixed <- x$coefficients$fixed
  # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  colnames(coefs_fixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  printCoefmat(coefs_fixed, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  # pasting the generalized r squared
  if (x$pseudoR2 != 0) {
    cat("\nPseudo R-squared: ")
    cat(formatC(x$pseudoR2, digits = digits))
    cat("\n")
  }

  cat("\nCoefficients (covariance):\n")
  x_tailup <- x$coefficients$params_object$tailup
  x_tailup <- data.frame(
    Effect = class(x_tailup),
    Parameter = names(x_tailup),
    Estimate = as.vector(x_tailup),
    is_known = x$is_known$tailup
  )
  x_taildown <- x$coefficients$params_object$taildown
  x_taildown <- data.frame(
    Effect = class(x_taildown),
    Parameter = names(x_taildown),
    Estimate = as.vector(x_taildown),
    is_known = x$is_known$taildown
  )
  x_euclid <- x$coefficients$params_object$euclid
  x_euclid <- data.frame(
    Effect = class(x_euclid),
    Parameter = names(x_euclid),
    Estimate = as.vector(x_euclid),
    is_known = x$is_known$euclid
  )
  x_nugget <- x$coefficients$params_object$nugget
  x_nugget <- data.frame(
    Effect = "nugget",
    Parameter = names(x_nugget),
    Estimate = as.vector(x_nugget),
    is_known = x$is_known$nugget
  )
  x_cov <- rbind(x_tailup, x_taildown, x_euclid, x_nugget)
  x_cov$Effect <- gsub("_", " ", x_cov$Effect)
  x_cov$Parameter <- gsub("de", "de (parsill)", x_cov$Parameter)
  logi1 <- x_cov$Estimate == 0 | x_cov$Estimate == Inf | (x_cov$Estimate == 1 & x_cov$Parameter == "scale")
  logi2 <- x_cov$is_known
  x_cov <- x_cov[!(logi1 & logi2), , drop = FALSE]
  x_cov <- x_cov[, -which(names(x_cov) == "is_known"), drop = FALSE]


  if (!is.null(x$coefficients$params_object$randcov)) {
    x_rand <- x$coefficients$params_object$randcov
    x_rand <- data.frame(
      Effect = "random",
      Parameter = names(x_rand),
      Estimate = as.vector(x_rand),
      is_known = x$is_known$randcov
    )
    logi1 <- x_rand$Estimate == 0
    logi2 <- x_rand$is_known
    x_rand <- x_rand[!(logi1 & logi2), , drop = FALSE]
    x_rand <- x_rand[, -which(names(x_rand) == "is_known"), drop = FALSE]
    x_cov <- rbind(x_cov, x_rand)
  }
  print.data.frame(format(x_cov, digits = digits),
    print.gap = 2L,
    quote = FALSE,
    row.names = FALSE
  )

  cat("\n")

  invisible(x)
}

#' @rdname print.SSN2
#' @method print summary.ssn_glm
#' @export
print.summary.ssn_glm <- function(x,
                                  digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
  # pasting the formula call
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  # pasting the residual summary
  cat("\nResiduals:\n")
  resQ <- c(
    min(x$residuals$response), quantile(x$residuals$response, p = c(0.25, 0.5, 0.75), na.rm = TRUE),
    max(x$residuals$response)
  )
  names(resQ) <- c("Min", "1Q", "Median", "3Q", "Max")
  print(resQ, digits = digits)

  # pasting the fixed coefficient summary
  cat("\nCoefficients (fixed):\n")
  coefs_fixed <- x$coefficients$fixed
  # colnames(coefs_fixed) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  colnames(coefs_fixed) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  printCoefmat(coefs_fixed, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)

  # pasting the generalized r squared
  if (x$pseudoR2 != 0) {
    cat("\nPseudo R-squared: ")
    cat(formatC(x$pseudoR2, digits = digits))
    cat("\n")
  }

  cat("\nCoefficients (covariance):\n")
  x_tailup <- x$coefficients$params_object$tailup
  x_tailup <- data.frame(
    Effect = class(x_tailup),
    Parameter = names(x_tailup),
    Estimate = as.vector(x_tailup),
    is_known = x$is_known$tailup
  )
  x_taildown <- x$coefficients$params_object$taildown
  x_taildown <- data.frame(
    Effect = class(x_taildown),
    Parameter = names(x_taildown),
    Estimate = as.vector(x_taildown),
    is_known = x$is_known$taildown
  )
  x_euclid <- x$coefficients$params_object$euclid
  x_euclid <- data.frame(
    Effect = class(x_euclid),
    Parameter = names(x_euclid),
    Estimate = as.vector(x_euclid),
    is_known = x$is_known$euclid
  )
  x_nugget <- x$coefficients$params_object$nugget
  x_nugget <- data.frame(
    Effect = "nugget",
    Parameter = names(x_nugget),
    Estimate = as.vector(x_nugget),
    is_known = x$is_known$nugget
  )
  x_dispersion <- x$coefficients$params_object$dispersion
  x_dispersion <- data.frame(
    Effect = "dispersion",
    Parameter = names(x_dispersion),
    Estimate = as.vector(x_dispersion),
    is_known = x$is_known$dispersion
  )
  x_cov <- rbind(x_tailup, x_taildown, x_euclid, x_nugget, x_dispersion)
  x_cov$Effect <- gsub("_", " ", x_cov$Effect)
  x_cov$Parameter <- gsub("de", "de (parsill)", x_cov$Parameter)
  logi1 <- x_cov$Estimate == 0 | x_cov$Estimate == Inf | (x_cov$Estimate == 1 & x_cov$Parameter == "scale")
  logi2 <- x_cov$is_known
  x_cov <- x_cov[!(logi1 & logi2), , drop = FALSE]
  x_cov <- x_cov[, -which(names(x_cov) == "is_known"), drop = FALSE]


  if (!is.null(x$coefficients$params_object$randcov)) {
    x_rand <- x$coefficients$params_object$randcov
    x_rand <- data.frame(
      Effect = "random",
      Parameter = names(x_rand),
      Estimate = as.vector(x_rand),
      is_known = x$is_known$randcov
    )
    logi1 <- x_rand$Estimate == 0
    logi2 <- x_rand$is_known
    x_rand <- x_rand[!(logi1 & logi2), , drop = FALSE]
    x_rand <- x_rand[, -which(names(x_rand) == "is_known"), drop = FALSE]
    x_cov <- rbind(x_cov, x_rand)
  }
  print.data.frame(format(x_cov, digits = digits),
    print.gap = 2L,
    quote = FALSE,
    row.names = FALSE
  )

  cat("\n")
  invisible(x)
}

#' @rdname print.SSN2
#' @method print anova.ssn_lm
#' @export
print.anova.ssn_lm <- function(x, digits = max(getOption("digits") - 2L, 3L),
                               signif.stars = getOption("show.signif.stars"), ...) {
  cat(attr(x, "heading")[1])
  cat("\n")
  cat(attr(x, "heading")[2])
  cat("\n")
  if ("Pr(>Chi2)" %in% colnames(x)) {
    P.values <- TRUE
    has.Pvalue <- TRUE
  } else {
    P.values <- FALSE
    has.Pvalue <- FALSE
  }
  printCoefmat(x, digits = digits, signif.stars = signif.stars, P.values = P.values, has.Pvalue = has.Pvalue, ...)
}

#' @rdname print.SSN2
#' @method print anova.ssn_glm
#' @export
print.anova.ssn_glm <- print.anova.ssn_lm
