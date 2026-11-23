test_that("torgegram works", {
  withr::local_pdf(NULL)
  # mf04p's only zero-hydro-distance pairs are cross-network (masked to zero
  # by mask_mat, not genuine coincident sites), so this must stay silent --
  # regression guard for the mask-vs-genuine-zero distinction in
  # check_torgegram_zero_distance()
  expect_no_warning(tg <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p))
  expect_s3_class(tg, "Torgegram")
  tg_names <- names(tg)
  expect_identical(c("flowcon", "flowuncon"), tg_names)
  tg_names_flowcon <- names(tg$flowcon)
  expect_identical(tg_names_flowcon, c("bins", "dist", "gamma", "np"))
  expect_s3_class(tg$flowcon, "data.frame")
  expect_equal(unlist(tg$flowcon[1, c("dist", "gamma", "np")]), c("dist" = 409.648, "gamma" = 0.167, "np" = 11), tolerance = 0.01)
  tg_names_flowuncon <- names(tg$flowuncon)
  expect_identical(tg_names_flowuncon, c("bins", "dist", "gamma", "np"))
  expect_s3_class(tg$flowuncon, "data.frame")
  expect_equal(unlist(tg$flowuncon[1, c("dist", "gamma", "np")]), c("dist" = 487.595, "gamma" = 2.739, "np" = 6), tolerance = 0.01)
  expect_invisible(plot(tg))
})

test_that("torgegram works euclid", {
  withr::local_pdf(NULL)
  expect_no_warning(tg <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p, type = "euclid"))
  expect_s3_class(tg, "Torgegram")
  tg_names <- names(tg)
  expect_identical(c("euclid"), tg_names)
  tg_names_euclid <- names(tg$euclid)
  expect_identical(tg_names_euclid, c("bins", "dist", "gamma", "np"))
  expect_equal(unlist(tg$euclid[1, c("dist", "gamma", "np")]), c("dist" = 545.053, "gamma" = 0.742, "np" = 40), tolerance = 0.01)
  expect_invisible(plot(tg))
})

test_that("torgegram works (partition factor)", {
  withr::local_pdf(NULL)
  expect_no_warning(
    tg <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p, partition_factor = ~ as.factor(netID), type = c("flowcon", "flowuncon", "euclid"))
  )
  expect_s3_class(tg, "Torgegram")
  tg_names <- names(tg)
  expect_identical(c("flowcon", "flowuncon", "euclid"), tg_names)
  tg_names_flowcon <- names(tg$flowcon)
  expect_identical(tg_names_flowcon, c("bins", "dist", "gamma", "np"))
  expect_s3_class(tg$flowcon, "data.frame")
  expect_equal(unlist(tg$flowcon[1, c("dist", "gamma", "np")]), c("dist" = 409.648, "gamma" = 0.167, "np" = 11), tolerance = 0.01)
  tg_names_flowuncon <- names(tg$flowuncon)
  expect_identical(tg_names_flowuncon, c("bins", "dist", "gamma", "np"))
  expect_s3_class(tg$flowuncon, "data.frame")
  expect_equal(unlist(tg$flowuncon[1, c("dist", "gamma", "np")]), c("dist" = 487.595, "gamma" = 2.739, "np" = 6), tolerance = 0.01)
  tg_names_euclid <- names(tg$euclid)
  expect_identical(tg_names_euclid, c("bins", "dist", "gamma", "np"))
  expect_equal(unlist(tg$euclid[1, c("dist", "gamma", "np")]), c("dist" = 310.580, "gamma" = 0.526, "np" = 17), tolerance = 0.01)
  expect_invisible(plot(tg))
})

test_that("torgegram robust works", {
  withr::local_pdf(NULL)
  expect_no_warning(tg <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p, robust = TRUE))
  expect_s3_class(tg, "Torgegram")
  tg_names <- names(tg)
  expect_identical(c("flowcon", "flowuncon"), tg_names)
  tg_names_flowcon <- names(tg$flowcon)
  expect_identical(tg_names_flowcon, c("bins", "dist", "gamma", "np"))
  expect_s3_class(tg$flowcon, "data.frame")
  expect_equal(unlist(tg$flowcon[1, c("dist", "gamma", "np")]), c("dist" = 409.648, "gamma" = 0.0365, "np" = 11), tolerance = 0.01)
  tg_names_flowuncon <- names(tg$flowuncon)
  expect_identical(tg_names_flowuncon, c("bins", "dist", "gamma", "np"))
  expect_s3_class(tg$flowuncon, "data.frame")
  expect_equal(unlist(tg$flowuncon[1, c("dist", "gamma", "np")]), c("dist" = 487.595, "gamma" = 2.2045, "np" = 6), tolerance = 0.01)
  expect_invisible(plot(tg))
})

test_that("torgegram cloud works", {
  withr::local_pdf(NULL)
  expect_no_warning(tg <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p, cloud = TRUE))
  expect_s3_class(tg, "Torgegram")
  tg_names <- names(tg)
  expect_identical(c("flowcon", "flowuncon"), tg_names)
  tg_names_flowcon <- names(tg$flowcon)
  expect_identical(tg_names_flowcon, c("dist", "gamma"))
  expect_s3_class(tg$flowcon, "data.frame")
  expect_equal(unlist(tg$flowcon[1, c("dist", "gamma")]), c("dist" = 1.9630e3, "gamma" = 7.8004e-4), tolerance = 0.01)
  tg_names_flowuncon <- names(tg$flowuncon)
  expect_identical(tg_names_flowuncon, c("dist", "gamma"))
  expect_s3_class(tg$flowuncon, "data.frame")
  expect_equal(unlist(tg$flowuncon[1, c("dist", "gamma")]), c("dist" = 1.2033e2, "gamma" = 5.2123e-3), tolerance = 0.01)
  expect_invisible(plot(tg))
})

test_that("warn_torgegram_zero_distance() warns only when observation coordinates are genuinely duplicated", {
  # checking hydro/euclid distance vectors directly for zeros is ambiguous
  # as hydro_mat is forced to zero for masked-out cross-network pairs; checking
  # physical coincidence directly in the coordinates sidesteps this entirely.
  no_dups <- sf::st_as_sf(data.frame(x = c(0, 1, 2), y = c(0, 1, 2)), coords = c("x", "y"))
  expect_no_warning(warn_torgegram_zero_distance(no_dups))

  has_dups <- sf::st_as_sf(data.frame(x = c(0, 1, 0), y = c(0, 1, 0)), coords = c("x", "y"))
  expect_warning(
    warn_torgegram_zero_distance(has_dups),
    "Zero distances observed between at least one pair"
  )
})

