# Copy the mf04p .ssn data to a local directory and read it into R
# When modeling with your .ssn object, you will load it using the relevant
# path to the .ssn data on your machine
copy_lsn_to_temp()
temp_path <- paste0(tempdir(), "/MiddleFork04.ssn")
mf04p <- ssn_import(
  temp_path,
  predpts = c("pred1km", "CapeHorn", "Knapp"),
  overwrite = TRUE
)

test_that("torgegram works", {
  tg <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p)
  tg_names <- names(tg$flowcon)
  expect_s3_class(tg$flowcon, "data.frame")
  expect_equal(names(tg$flowcon), tg_names)
  expect_s3_class(tg$flowuncon, "data.frame")
  expect_equal(names(tg$flowuncon), tg_names)
  expect_false("euclid" %in% names(tg))
  expect_invisible(plot(tg))
})

test_that("torgegram works euclid", {
  tg <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p, type = "euclid")
  tg_names <- names(tg$euclid)
  expect_s3_class(tg$euclid, "data.frame")
  expect_equal(names(tg$euclid), tg_names)
  expect_false(any(c("flowcon", "flowuncon") %in% names(tg)))
  expect_invisible(plot(tg))
})

test_that("torgegram works (partition factor)", {
  tg <- Torgegram(Summer_mn ~ ELEV_DEM, mf04p, partition_factor = ~ as.factor(netID))
  tg_names <- names(tg$flowcon)
  expect_s3_class(tg$flowcon, "data.frame")
  expect_equal(names(tg$flowcon), tg_names)
  expect_s3_class(tg$flowuncon, "data.frame")
  expect_equal(names(tg$flowuncon), tg_names)
  expect_false("euclid" %in% names(tg))
  expect_invisible(plot(tg))
})
