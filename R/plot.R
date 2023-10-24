#' Plot fitted model diagnostics
#'
#' @description Plot fitted model diagnostics such as residuals vs fitted values,
#'   quantile-quantile, scale-location, Cook's distance, residuals vs leverage,
#'   and Cook's distance vs leverage.
#'
#' @param x A fitted model object from [ssn_lm()] or [ssn_glm()].
#' @param which An integer vector taking on values between 1 and 6, which indicates
#'   the plots to return. Available plots are described in Details. If \code{which}
#'   has length greater than one, additional plots are stepped through in order
#'   using \code{<Return>}. The default is \code{which = c(1, 2)}
#' @param ... Other arguments passed to other methods.
#'
#' @details For all fitted model objects,, the values of \code{which} make the
#'   corresponding plot:
#'   \itemize{
#'     \item 1: Standardized residuals vs fitted values (of the response)
#'     \item 2: Normal quantile-quantile plot of standardized residuals
#'     \item 3: Scale-location plot of standardized residuals
#'     \item 4: Cook's distance
#'     \item 5: Standardized residuals vs leverage
#'     \item 6: Cook's distance vs leverage
#'   }
#'
#' @return No return value. Function called for plotting side effects.
#'
#' @name plot.SSN2
#' @method plot ssn_lm
#' @order 1
#' @export
#'
#' @seealso [plot.Torgegram()]
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#'
#' ssn_mod <- ssn_lm(
#'   formula = Summer_mn ~ ELEV_DEM,
#'   ssn.object = mf04p,
#'   tailup_type = "exponential",
#'   additive = "afvArea"
#' )
#' plot(ssn_mod, which = 1)
plot.ssn_lm <- function(x, which, ...) {
  if (missing(which)) {
    which <- c(1, 2)
  }

  if (any(!(which %in% 1:10))) {
    stop("Values of which can only take on 1, 2, 3, 4, 5, or 6.", call. = FALSE)
  }

  # setting old graphical parameter value
  oldpar <- par(no.readonly = TRUE)
  # setting exit handler
  on.exit(par(ask = oldpar$ask), add = TRUE)
  # set ask
  if (length(which) > 1) {
    par(ask = TRUE)
  }


  cal <- x$call
  if (!is.na(m.f <- match("formula", names(cal)))) {
    cal <- cal[c(1, m.f)]
    names(cal)[2L] <- ""
  }
  cc <- deparse(cal, 80)
  nc <- nchar(cc[1L], "c")
  abbr <- length(cc) > 1 || nc > 75
  sub.caption <- if (abbr) {
    paste(substr(cc[1L], 1L, min(75L, nc)), "...")
  } else {
    cc[1L]
  }

  # plot 1
  if (1 %in% which) {
    plot(
      x = fitted(x),
      y = rstandard(x),
      main = "Standardized Residuals vs Fitted",
      xlab = "Fitted values",
      ylab = "Standardized residuals",
      ...
    )
    title(sub = sub.caption)
    abline(h = 0, lty = 3, col = "gray")
  }

  # plot 2
  if (2 %in% which) {
    qqnorm(rstandard(x), main = "Normal Q-Q", ylab = "Standardized residuals", ...)
    qqline(rstandard(x), ...)
    title(sub = sub.caption)
  }


  # plot 3
  if (3 %in% which) {
    plot(
      x = fitted(x),
      y = sqrt(abs(rstandard(x))),
      main = "Scale-Location",
      xlab = "Fitted values",
      ylab = as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Standardized residuals")))),
      ...
    )
    title(sub = sub.caption)
  }




  # plot 4
  if (4 %in% which) {
    plot(
      x = seq_len(x$n),
      y = cooks.distance(x),
      xlab = "Obs. Number",
      ylab = "Cook's distance",
      main = "Cook's Distance",
      type = "h",
      ...
    )
    title(sub = sub.caption)
  }


  # plot 5
  if (5 %in% which) {
    plot(
      x = hatvalues(x),
      y = rstandard(x),
      xlab = "Leverage",
      ylab = "Standardized residuals",
      main = "Standardized residuals vs Leverage",
      ...
    )
    title(sub = sub.caption)
    abline(h = 0, v = 0, lty = 3, col = "gray")
  }

  # plot 6
  if (6 %in% which) {
    plot(
      x = hatvalues(x),
      y = cooks.distance(x),
      xlab = "Leverage",
      ylab = "Cook's distance",
      main = "Cook's dist vs Leverage",
      ...
    )
    title(sub = sub.caption)
  }
}

#' @rdname plot.SSN2
#' @method plot ssn_glm
#' @export
plot.ssn_glm <- plot.ssn_lm



#' Plot Torgegram
#'
#' @description Plot Torgegram
#'
#' @param x A Torgegram object from [Torgegram()].
#' @param type The type of semivariogram. Can take character values that are a subset
#'   of objects in \code{x}. The default is \code{names(x)}.
#' @param separate When \code{type} is length greater than one, whether each
#'   \code{type} be placed in a separate plot. The default is \code{FALSE}.
#' @param ... Other arguments passed to other methods.
#'
#' @return No return value. Function called for plotting side effects.
#'
#' @name plot.Torgegram
#' @method plot Torgegram
#' @export
#'
#' @seealso [plot.SSN2]
#'
#' @examples
#' # Copy the mf04p .ssn data to a local directory and read it into R
#' # When modeling with your .ssn object, you will load it using the relevant
#' # path to the .ssn data on your machine
#' copy_lsn_to_temp()
#' temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
#' mf04p <- ssn_import(temp_path, overwrite = TRUE)
#'
#' tg <- Torgegram(Summer_mn ~ 1, mf04p)
#' plot(tg)
plot.Torgegram <- function(x, type, separate = FALSE, ...) {
  if (missing(type)) {
    type <- names(x)
  }

  if (any(!type %in% names(x))) {
    stop("\"type\" values must be contained in the names of x, the Torgegram object.", call. = FALSE)
  }

  x <- lapply(type, function(y) {
    x_sub <- x[[y]]
    x_sub[["type"]] <- y
    x_sub
  })

  x <- do.call("rbind", x)
  x$type <- droplevels(factor(x$type, levels = c("flowcon", "flowuncon", "euclid")))
  # scale to [1, 3]
  x$cex <- (x$np - min(x$np)) / (max(x$np) - min(x$np)) * 2 + 1
  col_key <- c("flowcon" = "#000000", "flowuncon" = "#E69F00", euclid = "#56B4E9")
  x$col <- col_key[as.character(x$type)]

  # cb_b8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
  #                                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  dotlist <- list(...)
  names_dotlist <- names(dotlist)

  # set defaults
  if (!"main" %in% names_dotlist) {
    dotlist$main <- "Torgegram"
  }

  if (!"xlab" %in% names_dotlist) {
    dotlist$xlab <- "Distance"
  }

  if (!"ylab" %in% names_dotlist) {
    dotlist$ylab <- "Semivariance"
  }

  if (!"pch" %in% names_dotlist) {
    dotlist$pch <- 19
  }


  # hard code to FALSE if length type only one
  if (length(type) == 1) {
    separate <- FALSE
  }

  if (separate) {
    # setting old graphical parameter value
    oldpar <- par(no.readonly = TRUE)
    # setting exit handler
    on.exit(par(ask = oldpar$ask), add = TRUE)
    par(ask = TRUE)

    tg_split <- split(x, x$type)

    invisible(lapply(tg_split, function(p) {
      do.call("plot", args = c(list(
        x = p$dist, y = p$gamma, cex = p$cex,
        xlim = c(0, max(x$dist)), ylim = c(0, max(x$gamma) * 1.6),
        col = p$col
      ), dotlist))
      legend(0, max(x$gamma) * 1.6, legend = levels(droplevels(p$type)), col = col_key[levels(droplevels(p$type))], pch = dotlist$pch)
    }))
  } else {
    do.call("plot", args = c(list(
      x = x$dist, y = x$gamma, cex = x$cex, ylim = c(0, max(x$gamma) * 1.6),
      col = x$col
    ), dotlist))
    legend(0, max(x$gamma) * 1.6, legend = levels(x$type), col = col_key[levels(x$type)], pch = dotlist$pch)
  }
}
