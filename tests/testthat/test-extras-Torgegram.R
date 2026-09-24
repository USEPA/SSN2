skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("EACF changes only the statistic for each Torgegram pair type", {
  types <- c("flowcon", "flowuncon", "euclid")
  for (reversed in c(FALSE, TRUE)) {
    ssn <- mf04p
    if (reversed) ssn$obs <- ssn$obs[rev(seq_len(NROW(ssn$obs))), , drop = FALSE]
    residuals <- residuals(lm(Summer_mn ~ ELEV_DEM, data = ssn$obs))
    products <- outer(residuals, residuals)
    products <- products[upper.tri(products)]
    distances <- get_dist_object(ssn, get_Torgegram_initial_object(types), NULL, FALSE)
    hydro <- as.matrix(distances$hydro_mat * distances$mask_mat)
    connected <- as.matrix(distances$b_mat) == 0
    euclid <- as.matrix(distances$euclid_mat)
    pair_distances <- list(
      flowcon = (hydro * connected)[upper.tri(hydro)],
      flowuncon = (hydro * !connected)[upper.tri(hydro)],
      euclid = euclid[upper.tri(euclid)]
    )
    for (partition in list(NULL, ~as.factor(netID))) {
      within_partition <- if (is.null(partition)) {
        rep(TRUE, length(products))
      } else {
        mat <- outer(ssn$obs$netID, ssn$obs$netID, "==")
        mat[upper.tri(mat)]
      }
      for (cutoff in list(NULL, 10000)) {
        for (cloud in c(FALSE, TRUE)) {
          args <- list(formula = Summer_mn ~ ELEV_DEM, ssn.object = ssn,
                       type = types, bins = 5, cutoff = cutoff,
                       cloud = cloud, partition_factor = partition)
          semivariogram <- do.call(Torgegram, args)
          acov <- do.call(Torgegram, c(args, list(eacf = TRUE)))
          expect_identical(names(acov), names(semivariogram))
          expect_identical(attr(acov, "cloud"), cloud)
          expect_true(attr(acov, "eacf"))
          expect_null(attr(semivariogram, "eacf"))
          for (type in types) {
            shared <- if (cloud) "dist" else c("bins", "dist", "np")
            expect_identical(acov[[type]][shared], semivariogram[[type]][shared])
            expect_named(acov[[type]], if (cloud) c("dist", "acov") else c("bins", "dist", "acov", "np"))
            dist <- pair_distances[[type]][within_partition]
            values <- products[within_partition]
            limit <- if (is.null(cutoff)) max(dist) / 2 else cutoff
            keep <- dist > 0 & dist <= limit
            if (cloud) {
              expected <- values[keep]
            } else {
              bin <- ceiling(dist[keep] / limit * 5)
              expected <- vapply(seq_len(5), function(i) {
                if (any(bin == i)) mean(values[keep][bin == i]) else NA_real_
              }, numeric(1))
            }
            expect_equal(as.numeric(acov[[type]]$acov), expected, tolerance = 1e-12)
          }
        }
      }
    }
  }
})

test_that("Euclidean EACF includes the same cross-network pairs as the semivariogram", {
  coords <- sf::st_coordinates(mf04p$obs)
  distances <- as.matrix(dist(coords))
  network <- ssn_get_netgeom(mf04p$obs)$NetworkID
  cross_network <- outer(network, network, "!=")
  cutoff <- max(distances)
  keep <- upper.tri(distances) & distances > 0 & distances <= cutoff
  expect_true(any(cross_network[keep]))
  acov <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p, type = "euclid",
                   cloud = TRUE, cutoff = cutoff, eacf = TRUE)
  expect_equal(NROW(acov$euclid), sum(keep))
  for (eacf in c(FALSE, TRUE)) {
    expect_error(Torgegram(Summer_mn ~ 1, mf04p, type = "crossnet", eacf = eacf),
                 "All elements of type must be")
  }
})

test_that("EACF plots support negative values and both output forms", {
  withr::local_pdf(NULL)
  for (cloud in c(FALSE, TRUE)) {
    acov <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p,
                     type = c("flowcon", "flowuncon", "euclid"), cloud = cloud, eacf = TRUE)
    expect_true(any(unlist(lapply(acov, function(x) x$acov)) < 0, na.rm = TRUE))
    expect_invisible(plot(acov))
  }
  expect_error(Torgegram(Summer_mn ~ 1, mf04p, robust = TRUE, eacf = TRUE),
               "robust = TRUE is not available")
  expect_error(Torgegram(Summer_mn ~ 1, mf04p, eacf = NA), "eacf must be")
})
