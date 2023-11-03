#' @title Split a prediction dataset in an \code{SSN} object
#'
#' @description The \command{splitPrediction} function is used to
#'   split prediction sets in an \code{SSN} object into smaller
#'   prediction sets. It returns a \code{SSN} object with additional
#'   prediction sets based on equal interval splits, a factor,
#'   integer, character or logical column stored within the prediction
#'   set, or a logical expression.
#'
#' @param ssn An \code{SSN} object.
#' @param predpts A character string representing the name of the
#'   prediction dataset
#' @param size_predpts numeric value representing the size of the new
#'   prediction sets. The existing prediction set is split equally to
#'   produce multiple prediction sets of this size
#' @param by character string representing the column name of type
#'   factor, integer, character or logical that the split will be
#'   based on
#' @param subset logical expression indicating which elements or rows
#'   to keep; missing values are taken as \code{FALSE}
#' @param id_predpts character string representing the new prediction
#'   dataset name. This value is only specified when the subset method is
#'   used
#' @param keep logical value indicating whether the original
#'   prediction dataset should be retained in the \code{SSN}
#'   object. Default is \code{TRUE}
#' @param drop_levels logical value indicating whether empty factor
#'   levels should be dropped in the \code{by} column when the new
#'   prediction dataset(s) are created. Default is \code{FALSE}
#' @param overwrite logical indicating whether the new prediction
#'   dataset shapefile should be deleted in the .ssn directory if it
#'   already exists. Default = \code{FALSE}
#'
#' @details Three methods have been provided to split prediction sets:
#'   \code{size_predpts}, \code{by}, and \code{subset}. The
#'   \code{size_predpts} method is used to split the existing prediction
#'   set into multiple equally-sized prediction sets. Note that the
#'   final prediction set may be smaller in size than the others if
#'   the total number of predictions is not evenly divisible by
#'   \code{size_predpts}. The \code{by} method is used if the prediction
#'   set is to be split into multiple new prediction sets based on an
#'   existing column of type factor, integer, character, or
#'   logical. The \code{subset} method is used to create one new
#'   prediction set based on a logical expression.
#'
#'   When more than one prediction dataset is created the prediction
#'   dataset names will be appended with a hyphen and prediction
#'   dataset number if more than one prediction dataset is
#'   created. For example, when "preds" is split using
#'   \code{size_predpts}, the new names will be "preds-1", "preds-2", and
#'   so forth.
#'
#'   When \code{keep=FALSE}, the prediction dataset is removed from
#'   the \code{SSN} object stored in memory, but is not deleted from
#'   the .ssn directory specified in \code{ssn$path}.
#'
#'   Note that, only one
#'   method may be specified when the \command{ssn_split_predpts}
#'   function is called. The distance matrices for the new prediction
#'   datasets must be created using the \code{ssn_create_distmat} before
#'   predictions can be made.
#'
#' @return returns the \code{SSN} specified in \code{ssn}, with one or more new prediction
#'   sets. Shapefiles of the new prediction sets are written to the
#'   .ssn directory designated in ssn$path.
#'
#' @name ssn_split_predpts
#' @export
#' @examples
#' ## Import SSN object
#' copy_lsn_to_temp() ## Only needed for this example
#' ssn <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
#'   predpts = c("pred1km.shp", "Knapp", "CapeHorn"),
#'   overwrite = TRUE
#' )
#'
#' ## Split predictions into size_predpts 200
#' ssn1 <- ssn_split_predpts(ssn, "CapeHorn",
#'   size_predpts = 200,
#'   keep = FALSE, overwrite = TRUE
#' )
#' names(ssn1$preds)
#' nrow(ssn1$preds[["CapeHorn-1"]])
#'
#' ## Split predictions using by method
#' ssn$preds$pred1km$net.fac <- as.factor(ssn$preds$pred1km$netID)
#' ssn2 <- ssn_split_predpts(ssn, "pred1km",
#'   by = "net.fac",
#'   overwrite = TRUE
#' )
#' names(ssn2$preds)
#'
#' ## Split predictions using subset method
#' ssn3 <- ssn_split_predpts(ssn, "pred1km",
#'   subset = ratio > 0.5,
#'   id_predpts = "RATIO_05", overwrite = TRUE
#' )
#' names(ssn3$preds)
ssn_split_predpts <- function(ssn, predpts, size_predpts, by,
                              subset, id_predpts, keep = TRUE,
                              drop_levels = FALSE,
                              overwrite = FALSE) {
  suppressWarnings({
    ## Check arguments
    if (missing(ssn) || missing(predpts)) {
      stop("Arguments ssn and predpts must be specified")
    }
    if (sum(c(
      !missing(size_predpts), !missing(by),
      !missing(subset)
    )) != 1) {
      stop("Only one of 'size_predpts', 'by' and 'subset' can be input to ssn_split_predpts")
    }
    if (!missing(id_predpts) && missing(subset)) {
      stop("Input 'id_predpts' should only be specified if input 'subset' is used")
    }
    if (!missing(subset) && missing(id_predpts)) {
      stop("Input 'id_predpts' must be specified if input 'subset' is used")
    }
    if (!predpts %in% names(ssn$preds)) {
      stop("predpts no found in ssn")
    }

    ## Get vector of pids in obs and all sets of preds
    pid.vec <- ssn_get_netgeometry(ssn$obs,
      netvars = "all",
      reformat = TRUE
    )
    pid.vec <- pid.vec$pid

    for (i in seq_len(length(ssn$preds))) {
      tmp <- ssn_get_netgeometry(ssn$preds[[i]],
        netvars = "all",
        reformat = TRUE
      )
      pid.vec <- append(pid.vec, tmp$pid)
    }

    ## Check for duplicated pids
    if (sum(duplicated(pid.vec)) > 0) {
      stop("duplicated pid values are not allowed in sites and/or prediction sites")
    }

    ## Get max pid
    max.pid <- max(pid.vec)

    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(ssn$path)

    ## REPLACE PID VALUES IF NECESSARY
    if (keep == TRUE) {
      ## Replace network geometry columns in predpts with values from
      ## netgeometry
      ng.cols <- c("netID", "rid", "upDist", "ratio", "pid", "locID")
      ind <- colnames(ssn$preds[[predpts]]) %in% ng.cols
      ssn$preds[[predpts]] <- ssn$preds[[predpts]][, !ind]
      ng.df <- ssn_get_netgeometry(ssn$preds[[predpts]],
        netvars = "all",
        reformat = TRUE
      )
      colnames(ng.df) <- ng.cols

      ## Avoid merging issues with sf & data.frame
      tmp <- ssn$preds[[predpts]]
      for (m in seq_len(length(ng.cols))) {
        tmp[, ng.cols[m]] <- ng.df[, ng.cols[m]]
      }

      ## Modify pid values in new prediction subsets
      tmp.pid <- tmp$pid + max.pid
      if (sum(tmp.pid %in% pid.vec) > 0) stop("Duplicate pids exist")
      tmp$pid <- tmp.pid

      pid.vec <- append(pid.vec, tmp.pid)
      rm(tmp.pid)

      ## Recreate netgeometry
      ind <- colnames(tmp) == "netgeometry"
      tmp <- tmp[, !ind]
      tmp <- create_netgeometry(tmp, "point")

      ## Replace predpts with updated data.frame
      ssn$preds[[predpts]] <- tmp
    }
    ## -------------------------------------------------------
    ## Split based on 'size_predpts'
    ## -------------------------------------------------------
    if (!missing(size_predpts)) {
      ## Get number of preds and number of chunks
      n.points <- nrow(ssn$preds[[predpts]])
      nchunks <- round((n.points / size_predpts) + 0.5)
      if (nchunks > 1000) {
        stop(paste(
          "Specified value of size_predpts would result in ",
          predpts,
          " being split into more than 1000 prediction sets",
          sep = ""
        ))
      }
      ## Get start row for each chunk
      chunks <- seq.int(1, n.points, size_predpts)

      ## Add a start row for the final phantom chunk
      if (tail(chunks, 1) != n.points) {
        chunks <- c(chunks, n.points + 1)
      }

      for (j in seq_len(length(chunks) - 1)) {
        id_predpts2 <- paste(predpts, "-", j, sep = "")

        ## Get chunk
        sub.data <-
          ssn$preds[[predpts]][chunks[j]:(chunks[j + 1] - 1), , drop = FALSE]


        ## write subset of prediction points to file
        fe.old <- file.exists(paste0(ssn$path, "/", id_predpts2, ".shp"))
        if (overwrite == TRUE & fe.old == TRUE) {
          st_delete(paste0(ssn$path, "/", id_predpts2, ".shp"),
            quiet = TRUE
          )
        }
        if (overwrite == FALSE & fe.old == TRUE) {
          stop(paste0("overwrite = FALSE and ", id_predpts2, ".shp already exists"))
        }

        ## Write shapefile
        st_write(sub.data, paste0(id_predpts2, ".shp"), quiet = TRUE)

        ## Add to existing ssn
        new.index <- length(ssn$preds) + 1
        ssn$preds[[new.index]] <- sub.data
        names(ssn$preds)[[new.index]] <- id_predpts2
      }
      ## Create vector of subset prediction names
      new.predids <- paste(predpts, "-", 1:(length(chunks) - 1), sep = "")

      ## -------------------------------------------------------
      ## Split based on subset
      ## -------------------------------------------------------
    } else if (!missing(subset)) {
      e <- substitute(subset)
      ## e <- substitute(CUMDRAINAG > 64)

      tmp.df <- ssn_get_data(ssn, predpts)

      values <- eval(e, ssn$preds[[predpts]], parent.frame())
      if (!is.logical(values)) {
        stop("Input 'subset' must evaluate to a logical value")
      }

      sub.data <- ssn$preds[[predpts]][values, ]

      ## write subset of prediction points to file
      fe.old <- file.exists(paste0(ssn$path, "/", id_predpts, ".shp"))
      if (overwrite == TRUE & fe.old == TRUE) {
        st_delete(paste0(ssn$path, "/", id_predpts, ".shp"),
          quiet = TRUE
        )
      }
      if (overwrite == FALSE & fe.old == TRUE) {
        stop(paste0("overwrite = FALSE and ", id_predpts, ".shp already exists"))
      }

      ## write subset of prediction points to file
      st_write(sub.data, paste0(id_predpts, ".shp"), quiet = TRUE)

      ## Add to existing ssn
      new.index <- length(ssn$preds) + 1
      ssn$preds[[new.index]] <- sub.data
      names(ssn$preds)[[new.index]] <- id_predpts

      new.predids <- id_predpts
      ## -------------------------------------------------------
      ## Split based on 'by'
      ## -------------------------------------------------------
    } else {
      if (!(by %in% colnames(ssn$preds[[predpts]]))) {
        stop(paste("Could not find column named ", by,
          " in point.coords entry of specified prediction points set",
          sep = ""
        ))
      }
      ## Extract predpts data.frame
      tmp.df <- ssn_get_data(ssn, predpts)
      current.class <- class(as.data.frame(tmp.df)[, by])

      ## Check that by is factor, integer, character or logical value
      if (!current.class %in% c("integer", "factor", "character", "logical")) {
        stop("'by' column must be an integer, factor, character, or logical")
      }
      ## Get levels to split by
      levels <- unique(as.data.frame(tmp.df)[, by])

      new.predids <- c()
      for (k in seq_len(length(levels))) {
        id_predpts <- paste(predpts, "-", by, "-", levels[k], sep = "")
        new.predids <- c(new.predids, id_predpts)

        relevent <- as.data.frame(tmp.df)[, by] == levels[k]
        sub.data <- ssn$preds[[predpts]][relevent, ]

        if (current.class %in% "factor" & drop_levels == TRUE) {
          sub.data[, by] <- droplevels(as.data.frame(sub.data)[, by])
        }

        ## write subset of prediction points to file
        fe.old <- file.exists(paste0(ssn$path, "/", id_predpts, ".shp"))
        if (overwrite == TRUE & fe.old == TRUE) {
          st_delete(paste0(ssn$path, "/", id_predpts, ".shp"),
            quiet = TRUE
          )
        }
        if (overwrite == FALSE & fe.old == TRUE) {
          stop(paste0("overwrite = FALSE and ", id_predpts, ".shp already exists"))
        }

        ## write subset of prediction points to file
        st_write(sub.data, paste0(id_predpts, ".shp"), quiet = TRUE)

        ## Add to existing ssn
        new.index <- length(ssn$preds) + 1
        ssn$preds[[new.index]] <- sub.data
        names(ssn$preds)[[new.index]] <- id_predpts
      }
    }
    if (keep == FALSE) {
      ind <- names(ssn$preds) == predpts
      ssn$preds <- ssn$preds[!ind]
    }

    return(ssn)
  })
}
