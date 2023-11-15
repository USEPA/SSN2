#' @title Subset an \code{SSN} object
#'
#' @description Returns an \code{SSN} object that has been subset based on a
#'   logical expression.
#' @param ssn An \code{SSN} object.
#' @param path The filepath to the .ssn folder, in string format,
#'   where the subset \code{SSN} will be saved.
#' @param subset A logical expression indicating which features to keep.
#' @param clip If \code{TRUE}, create a subset of the edges and
#'   prediction sites, based on the same logical expression used to
#'   subset the observed sites.  Default = \code{FALSE}.
#' @param overwrite If \code{TRUE}, overwrite the folder specified in
#'   \code{path} if it exists. Default = FALSE.
#'
#' @details This function creates a subset of the original \code{SSN}
#'   object based on a logical expression defined in the \code{subset}
#'   argument. The \code{subset} argument is treated as an expression
#'   within \code{ssn_subset()} and so the full argument is not a
#'   string; although values in factor or character format will still
#'   require quotes (see examples). If \code{clip = TRUE}, then the
#'   columns referred to in \code{subset} must be present in the edges
#'   and all of the prediction datasets (if present in the \code{SSN}
#'   object). Note that features with missing values in the \code{subset}
#'   expression are treated as false and are not included in the
#'   subset \code{SSN} object.
#'
#'   Once the subset \code{SSN} object has been written to the local
#'   directory, it is re-imported using
#'   \code{\link[SSN2]{ssn_import}}. During this process, the
#'   binaryID.db is recreated. If distance matrices exist in the
#'   original \code{SSN} object, they are not copied or recalculated
#'   for the new \code{SSN} object. Users will need to run the
#'   \code{\link[SSN2]{ssn_create_distmat}} to create the distance
#'   matrices before fitting models to the data in the subset
#'   \code{SSN}.
#'
#' @return an object of class \code{SSN}, which is stored locally in the .ssn
#'   directory specified in \code{path}. It also creates and
#'   stores an SQLite database, binaryID.db, within the .ssn
#'   directory.
#'
#' @name ssn_subset
#' @export
#' @examples
#' ## Import SSN object
#' copy_lsn_to_temp() ## Only needed for this example
#' mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
#'   predpts = c("pred1km.shp", "Knapp"),
#'   overwrite = TRUE
#' )
#'
#' ## Subset SSN observations, edges, and prediction sites on network 1
#' ssn.sub1 <- ssn_subset(mf04p,
#'   path = paste0(tempdir(), "/subset1.ssn"),
#'   subset = netID == 1, clip = TRUE,
#'   overwrite = TRUE
#' )
#'
#' ## Subset SSN observations, removing two sites
#' ssn.sub2 <- ssn_subset(mf04p,
#'   path = paste0(tempdir(), "/subset2.ssn"),
#'   subset = !COMID %in% c("23519461", "23519365"),
#'   overwrite = TRUE
#' )
ssn_subset <- function(ssn, path, subset, clip = FALSE, overwrite = FALSE) {
  file <- path

  suppressWarnings({
    if (!file.exists(file)) {
      dir.create(file)
    } else {
      if (overwrite == FALSE) stop("file exists and overwrite = FALSE")
      if (overwrite == TRUE) {
        unlink(file, recursive = TRUE)
        dir.create(file)
      }
    }

    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(file)

    ssn.tmp <- ssn
    ssn.tmp$path <- getwd()
    pred.len <- length(ssn.tmp$preds)
    if (pred.len > 0) {
      pred.names.vec <- attributes(ssn.tmp$preds)$names
    }

    ## ##Check to see if attribute exists
    s <- deparse(substitute(subset))
    ## ##s<- deparse(substitute(netID == 2))

    ## If the subset expression contains special netgeometry columns
    if (grepl(
      pattern = "netID|rid|upDist|ratio|pid|locID",
      s
    ) == TRUE) {
      ## Identify which netgeometry columns are used
      data.cols <- c("netID", "rid", "upDist", "ratio", "pid", "locID")
      netg.obs <- ssn_get_netgeometry(ssn.tmp$obs, reformat = TRUE)
      colnames(netg.obs) <- data.cols

      ind <- colnames(ssn.tmp$obs) %in% data.cols
      colnames(ssn.tmp$obs)[ind] <- paste0("._", colnames(ssn.tmp$obs)[ind], "_")

      ssn.tmp$obs <- cbind(ssn.tmp$obs, netg.obs)
      rm(ind)
    }

    ## Select observations based on subset expression
    ind <- eval(substitute(subset), ssn.tmp$obs)
    ## ind<- eval(substitute(netID == 2), ssn.tmp$obs)
    ind.na <- is.na(ind)
    ind[ind.na] <- FALSE
    rm(ind.na)

    if (sum(ind) == 0) {
      stop("No records were selected based on subset expression")
    }

    ## Subset observations
    ssn.tmp$obs <- ssn.tmp$obs[ind, ]
    rm(ind)

    ## Fix netgeometry column names if necessary
    if (exists("netg.obs")) {
      ind <- colnames(ssn.tmp$obs) %in% data.cols
      ssn.tmp$obs <- ssn.tmp$obs[, !ind]

      ind2 <- colnames(ssn.tmp$obs) %in% paste0("._", data.cols, "_")
      colnames(ssn.tmp$obs)[ind2] <- substr(
        colnames(ssn.tmp$obs)[ind2],
        3, nchar(colnames(ssn.tmp$obs[ind2]))
      )
      colnames(ssn.tmp$obs)[ind2] <- substr(
        colnames(ssn.tmp$obs)[ind2], 1,
        nchar(colnames(ssn.tmp$obs[ind2])) - 1
      )
      rm(ind, ind2)
    }

    ## Write sites shapefile
    st_write(ssn.tmp$obs, paste0(file, "/sites.shp"), quiet = TRUE)


    if (clip == FALSE) {
      ssn.files <- list.files(ssn$path)
      for (i in seq_len(length(ssn.files))) {
        fn.old <- file.path(ssn$path, ssn.files[i])
        if (basename(fn.old) != "distance") {
          if (substr(basename(fn.old), 1, 5) != "sites") {
            fn.new <- file.path(ssn.tmp$path, ssn.files[i])
            file.copy(fn.old, fn.new, overwrite = TRUE)
          }
        }
      }
      rm(fn.old, fn.new)

      ## Clip everything based on subset expression
    } else {
      ## If subset expression depends on netgeometry columns then extract and
      ## rename them
      if (exists("netg.obs")) {
        netg.edges <- ssn_get_netgeometry(ssn.tmp$edges, reformat = TRUE)
        colnames(netg.edges) <- data.cols[1:3]
        ind <- colnames(ssn.tmp$edges) %in% data.cols[1:3]
        colnames(ssn.tmp$edges)[ind] <- paste0("._", colnames(ssn.tmp$edges)[ind], "_")
        ssn.tmp$edges <- cbind(ssn.tmp$edges, netg.edges)
        rm(ind)
      }

      ## Subset edges
      ind.edges <- eval(substitute(subset), ssn.tmp$edges)
      ## ind.edges<- eval(substitute(netID == 2), ssn.tmp$edges)
      ind.na <- is.na(ind.edges)
      ind.edges[ind.na] <- FALSE
      rm(ind.na)

      if (sum(ind.edges) == 0) {
        stop("No edges have have been selected based on subset expression")
      }

      if (exists("netg.edges")) {
        ind <- colnames(ssn.tmp$edges) %in% data.cols[1:3]
        ssn.tmp$edges <- ssn.tmp$edges[, !ind]

        ind2 <- colnames(ssn.tmp$edges) %in% paste0("._", data.cols[1:3], "_")
        colnames(ssn.tmp$edges)[ind2] <- substr(
          colnames(ssn.tmp$edges)[ind2],
          3, nchar(colnames(ssn.tmp$edges[ind2]))
        )
        colnames(ssn.tmp$edges)[ind2] <- substr(
          colnames(ssn.tmp$edges)[ind2], 1,
          nchar(colnames(ssn.tmp$edges[ind2])) - 1
        )
        rm(ind, ind2)
      }

      ## Subset edges
      edges.sub <- ssn.tmp$edges[ind.edges, ]

      ## Save subset of edges
      st_write(edges.sub, paste0(ssn.tmp$path, "/edges.shp"), quiet = TRUE)


      ## Subset prediction points
      if (pred.len > 0) {
        for (i in 1:pred.len) {
          pred.name <- attributes(ssn.tmp$preds)$names[i]

          if (exists("netg.obs")) {
            netg.pred <- ssn_get_netgeometry(ssn.tmp$preds[[pred.name]], reformat = TRUE)
            colnames(netg.pred) <- data.cols
            ind <- colnames(ssn.tmp$preds[[pred.name]]) %in% data.cols
            colnames(ssn.tmp$preds[[pred.name]])[ind] <-
              paste0("._", colnames(ssn.tmp$pred[[pred.name]])[ind], "_")
            ssn.tmp$preds[[pred.name]] <- cbind(ssn.tmp$preds[[pred.name]], netg.pred)
            rm(ind)
          }

          ind.preds <- eval(substitute(subset), ssn.tmp$preds[[pred.name]])
          ## ind.preds <- eval(substitute(netID == 2), ssn.tmp$preds[[pred.name]])
          ind.na <- is.na(ind.preds)
          ind.preds[ind.na] <- FALSE
          rm(ind.na)

          ## Write subset of predictions
          if (sum(ind.preds) > 0) {
            if (exists("netg.pred")) {
              ind <- colnames(ssn.tmp$preds[[pred.name]]) %in% data.cols
              ssn.tmp$preds[[pred.name]] <- ssn.tmp$preds[[pred.name]][, !ind]

              ind2 <- colnames(ssn.tmp$preds[[pred.name]]) %in% paste0("._", data.cols, "_")
              colnames(ssn.tmp$preds[[pred.name]])[ind2] <- substr(
                colnames(ssn.tmp$preds[[pred.name]])[ind2],
                3, nchar(colnames(ssn.tmp$preds[[pred.name]][ind2]))
              )
              colnames(ssn.tmp$preds[[pred.name]])[ind2] <- substr(
                colnames(ssn.tmp$preds[[pred.name]])[ind2], 1,
                nchar(colnames(ssn.tmp$preds[[pred.name]][ind2])) - 1
              )
              rm(ind, ind2)
            }
            ## Subset predictions
            preds.sub <- ssn.tmp$preds[[pred.name]][ind.preds, ]
            st_write(preds.sub, paste0(ssn.tmp$path, "/", pred.name, ".shp"),
              quiet = TRUE
            )
            rm(preds.sub)
          } else {
            ind.rm <- pred.names.vec == pred.name
            pred.names.vec <- pred.names.vec[!ind.rm]
          }
          rm(ind.preds)
        }
      }

      ## Get list of unique netIDs
      ind.dup <- !duplicated(edges.sub$netID)
      netID.list <- edges.sub$netID[ind.dup]

      # copy netID files
      for (i in seq_len(length(netID.list))) {
        fn.old <- file.path(ssn$path, paste("netID", netID.list[i], ".dat", sep = ""))
        fn.new <- file.path(ssn.tmp$path, paste("netID", netID.list[i], ".dat", sep = ""))
        file.copy(fn.old, fn.new, overwrite = TRUE)
      }
      rm(fn.new, fn.old)
    }

    ## Import subset SSN
    if (pred.len == 0) {
      ssn.tmp <- ssn_import(ssn.tmp$path, overwrite = TRUE)
    }

    if (pred.len > 0) {
      ## pred.names.vec <- attributes(ssn.tmp$preds)$names
      ssn.tmp <- ssn_import(ssn.tmp$path, overwrite = TRUE)

      for (j in seq_len(length(pred.names.vec))) {
        ssn.tmp <- ssn_import_predpts(ssn.tmp, pred.names.vec[j])
      }
    }

    return(ssn.tmp)
  })
}
