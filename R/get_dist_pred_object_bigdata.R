# get_dist_pred_object_bigdata <- function(newdata_list_row, network_index_obs, pid_obs, ssn.object, newdata_name, params_object, additive, anisotropy) {
#
#
#   newdata_netgeom <- ssn_get_netgeom(newdata_list_row)
#   network_index_pred <- newdata_netgeom$NetworkID
#   pid_pred <- newdata_netgeom$pid
#
#
#
#   dist_pred_matlist <- get_dist_pred_matlist_bigdata(network_index_obs, pid_obs, network_index_pred, pid_pred, ssn.object, newdata_name, params_object, additive, newdata_list_row)
#
#   # see whether euclid is none to avoid unnecessary computations
#   euclid_none <- inherits(params_object$euclid, "euclid_none")
#   if (euclid_none) {
#     # not needed if no euclid covariance
#     dist_pred_matlist <-c(dist_pred_matlist, list(euclid_pred = NULL))
#     # dist_pred_matlist$euclid_matix <- NULL does not return anything
#   } else {
#     # store euclid coordinates (in original data order)
#     obs_coords <- sf::st_coordinates(ssn.object$obs)
#     pred_coords <- sf::st_coordinates(newdata_list_row)
#     if (anisotropy) {
#       dist_pred_matlist$.xcoord_obs <- obs_coords[, 1, drop = TRUE]
#       dist_pred_matlist$.ycoord_obs <- obs_coords[, 2, drop = TRUE]
#       dist_pred_matlist$.xcoord_pred <- pred_coords[, 1, drop = TRUE]
#       dist_pred_matlist$.ycoord_pred <- pred_coords[, 2, drop = TRUE]
#     } else {
#       x_diffs <- t(outer(obs_coords[, 1, drop = TRUE], pred_coords[, 1, drop = TRUE], "-"))
#       y_diffs <- t(outer(obs_coords[, 2, drop = TRUE], pred_coords[, 2, drop = TRUE], "-"))
#       euclid_pred_mat <- sqrt(x_diffs^2 + y_diffs^2)
#       rownames(euclid_pred_mat) <- pid_pred
#       colnames(euclid_pred_mat) <- pid_obs
#       dist_pred_matlist$euclid_pred_mat <- Matrix::Matrix(euclid_pred_mat, sparse = TRUE)
#     }
#   }
#
#   dist_pred_matlist
# }
#
# get_dist_pred_matlist_bigdata <- function(network_index_obs, pid_obs, network_index_pred, pid_pred, ssn.object, newdata_name, params_object, additive, newdata_list_row) {
#
#   # see whether tailup and taildown are none to avoid unnecssary computations
#   tailup_none <- inherits(params_object$tailup, "tailup_none")
#   taildown_none <- inherits(params_object$taildown, "taildown_none")
#
#   # return all NULL if they are both none (no stream distance needed)
#   if (tailup_none && taildown_none) {
#     dist_pred_matlist <- list(
#       distjunc_pred = NULL,
#       mask_pred = NULL,
#       a_pred = NULL,
#       b_pred = NULL,
#       hydro_pred = NULL,
#       w_pred = NULL
#     )
#   } else {
#     # get dist junction matrices as a list (for efficiency, do things
#     # network by network and then combine so zeroes populate accordingly)
#     distjunc_pred_matlist <- get_distjunc_pred_matlist_bigdata(network_index_obs, pid_obs, network_index_pred, pid_pred, ssn.object, newdata_name)
#     # get other matrices as a list
#     dist_pred_matlist <- list(
#       distjunca_pred = distjunc_pred_matlist$distjunca,
#       distjuncb_pred = distjunc_pred_matlist$distjuncb,
#       mask_pred = get_mask_pred_matlist_bigdata(network_index_pred, network_index_obs, distjunc_pred_matlist),
#       a_pred = get_a_pred_matlist_bigdata(distjunc_pred_matlist),
#       b_pred = get_b_pred_matlist_bigdata(distjunc_pred_matlist),
#       hydro_pred = get_hydro_pred_matlist_bigdata(distjunc_pred_matlist)
#     )
#
#     # if only taildown covariance, do not need additive matrix
#     if (tailup_none) {
#       # store additive matrix as NULL
#       dist_pred_matlist <- c(dist_pred_matlist, list(w_pred = NULL))
#     } else {
#       obs_additive <- ssn.object$obs[[additive]]
#       pred_additive <- newdata_list_row[[additive]]
#       dist_pred_matlist$w_pred <- get_w_pred_matlist_bigdata(
#         obs_additive,
#         pred_additive,
#         dist_pred_matlist$b_pred,
#         dist_pred_matlist$mask_pred
#       )
#     }
#   }
#   dist_pred_matlist
# }
#
# get_distjunc_pred_matlist_bigdata <- function(network_index_obs, pid_obs, network_index_pred, pid_pred, ssn.object, newdata_name) {
#
#
#   network_index_obs <- as.numeric(as.character(network_index_obs))
#   network_index_pred <- as.numeric(as.character(network_index_pred))
#   network_index_vals <- sort(unique(c(network_index_obs, network_index_pred)))
#
#   # get network pid
#   # might require ordering
#   ord_pid_obs <- order(as.numeric(pid_obs))
#   pid_obs_orig <- pid_obs
#   pid_obs <- pid_obs[ord_pid_obs]
#   distmata <- matrix(0, nrow = 1, ncol = length(pid_obs))
#   rownames(distmata) <- pid_pred
#   colnames(distmata) <- pid_obs
#   distmatb <- distmata
#
#   # find distance junction prediction matrices (as a list) separately for each
#   # network index
#
#   # operate differently if observations are raw prediction data or induced
#   # by NA values in the response
#   if (newdata_name == ".missing") {
#     # on the disk, distance matrices are stored by network
#     workspace_name <- paste("dist.net", x, ".bmat", sep = "")
#     # path to the distance matrices on disk
#     path <- file.path(ssn.object$path, "distance", "obs", workspace_name)
#     # check to see if the file exists on the disk
#     if (!file.exists(path)) {
#       stop("Unable to locate required distance matrix", call. = FALSE)
#     }
#     distmat <- fm.open(path)
#     rownames_val <- rownames(distmat)
#     # find observations that are used to build model
#     which_obs <- rownames_val %in% pid_obs[network_index_obs == x]
#     # find prediction observations
#     which_pred <- colnames(distmat) %in% pid_pred
#     distmata <- distmat[which_obs, which_pred, drop = FALSE]
#     distmatb <- distmat[which_pred, which_obs, drop = FALSE]
#     close(distmat)
#   } else {
#     for (x in network_index_vals) { # order by net ID
#       # on the disk, distance matrices are stored by network
#       workspace.name.a <- paste("dist.net", x,
#                                 ".a.bmat",
#                                 sep = ""
#       )
#       workspace.name.b <- paste("dist.net", x,
#                                 ".b.bmat",
#                                 sep = ""
#       )
#       # path to the distance matrices on disk
#       path.a <- file.path(
#         ssn.object$path,
#         "distance", newdata_name, workspace.name.a
#       )
#       # check to see if the file exists on the disk
#       if (!file.exists(path.a)) {
#         stop("Unable to locate required distance matrix", call. = FALSE)
#       }
#       path.b <- file.path(
#         ssn.object$path,
#         "distance", newdata_name, workspace.name.b
#       )
#       # check to see if the file exists on the disk
#       if (!file.exists(path.b)) {
#         stop("Unable to locate required distance matrix", call. = FALSE)
#       }
#       # distance matrix a
#       distjunc_fm <- fm.open(path.a)
#       rows_keep <- pid_pred[rownames(distmata) %in% rownames(distjunc_fm)]
#       cols_keep <- pid_obs[colnames(distmata) %in% colnames(distjunc_fm)]
#       rows_keep_distjunc_fm <- which(rownames(distjunc_fm) %in% rows_keep)
#       cols_keep_distjunc_fm <- which(colnames(distjunc_fm) %in% cols_keep)
#       if (length(rows_keep_distjunc_fm) > 0 && length(cols_keep_distjunc_fm)) {
#         rows_keep_distjunc <- which(rownames(distmata) %in% rows_keep)
#         cols_keep_distjunc <- which(colnames(distmata) %in% cols_keep)
#         distjunc_fm_sub <- distjunc_fm[rows_keep_distjunc_fm, cols_keep_distjunc_fm]
#         # don't need ordering here like we do with local indexing
#         distmata[rows_keep_distjunc, cols_keep_distjunc] <- distjunc_fm_sub
#       }
#       close(distjunc_fm)
#       # distance matrix b
#       distjunc_fm <- fm.open(path.b)
#       rows_keep <- pid_pred[rownames(distmatb) %in% rownames(distjunc_fm)]
#       cols_keep <- pid_obs[colnames(distmatb) %in% colnames(distjunc_fm)]
#       rows_keep_distjunc_fm <- which(rownames(distjunc_fm) %in% rows_keep)
#       cols_keep_distjunc_fm <- which(colnames(distjunc_fm) %in% cols_keep)
#       if (length(rows_keep_distjunc_fm) > 0 && length(cols_keep_distjunc_fm)) {
#         rows_keep_distjunc <- which(rownames(distmatb) %in% rows_keep)
#         cols_keep_distjunc <- which(colnames(distmatb) %in% cols_keep)
#         distjunc_fm_sub <- distjunc_fm[rows_keep_distjunc_fm, cols_keep_distjunc_fm]
#         # don't need ordering here like we do with local indexing
#         distmatb[rows_keep_distjunc, cols_keep_distjunc] <- distjunc_fm_sub
#       }
#       # do we need to order by netid/pid here?
#       close(distjunc_fm)
#     }
#   }
#   distjunc_pred_matlist <- list(distjunca = distmata, distjuncb = t(distmatb))
#
# }
#
# get_mask_pred_matlist_bigdata <- function(network_index_pred, network_index_obs, distjunc_pred_matlist) {
#   mask_mat <- outer(network_index_pred, network_index_obs, FUN = "==")
#   rownames(mask_mat) <- rownames(distjunc_pred_matlist$distjunca)
#   colnames(mask_mat) <- colnames(distjunc_pred_matlist$distjunca)
#   1 * mask_mat
# }
#
# get_a_pred_matlist_bigdata <- function(distjunc_pred_matlist) {
#   pmax(distjunc_pred_matlist$distjunca, t(distjunc_pred_matlist$distjuncb))
# }
#
# get_b_pred_matlist_bigdata <- function(distjunc_pred_matlist) {
#   pmin(distjunc_pred_matlist$distjunca, t(distjunc_pred_matlist$distjuncb))
# }
#
# get_hydro_pred_matlist_bigdata <- function(distjunc_pred_matlist) {
#   distjunc_pred_matlist$distjunca + t(distjunc_pred_matlist$distjuncb)
# }
#
# get_w_pred_matlist_bigdata <- function(obs_additive, pred_additive, b_pred_matlist, mask_pred_matlist) {
#   obs_additive_matrix_val <- matrix(replicate(length(pred_additive), obs_additive, simplify = "array"), nrow = 1)
#   pred_additive_matrix_val <- matrix(replicate(length(obs_additive), pred_additive, simplify = "array"), ncol = 1)
#   w_matrix_val <- pmin(obs_additive_matrix_val, t(pred_additive_matrix_val)) / pmax(obs_additive_matrix_val, t(pred_additive_matrix_val))
#   # sqrt(w_matrix_val) # bug obviously
#   sqrt(w_matrix_val) * (b_pred_matlist == 0) * mask_pred_matlist
# }
