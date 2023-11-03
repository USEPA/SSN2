test_that("ssn manipulation functions work", {

  expect_no_error(copy_lsn_to_temp())

  mf04p <- ssn_import(paste0(tempdir(), "/MiddleFork04.ssn"),
                      predpts = c("pred1km", "Knapp"),
                      overwrite = TRUE
  )
  expect_s3_class(mf04p, "SSN")


  expect_no_error(ssn_import_predpts(mf04p, predpts = "CapeHorn"))

  expect_no_error(ssn_update_path(mf04p, mf04p$path))

  expect_message(ssn_create_distmat(mf04p, predpts = c("pred1km")))
  expect_no_error(ssn_create_distmat(mf04p, overwrite = TRUE))
  expect_no_error(ssn_create_distmat(mf04p, predpts = c("pred1km", "Knapp"),
                     overwrite = TRUE, among_predpts = TRUE, only_predpts = TRUE))



  expect_no_error(ssn_subset(mf04p, path = paste0(tempdir(), "/subset1.ssn"),
            subset = netID == 1, clip = TRUE, overwrite = TRUE))


  expect_no_error(ssn_split_predpts(mf04p, "pred1km", size_predpts = 200,
                            keep = FALSE, overwrite = TRUE))
  mf04p$preds$pred1km$net.fac <- as.factor(mf04p$preds$pred1km$netID)
  expect_no_error(ssn_split_predpts(mf04p, "pred1km", by = "net.fac", overwrite = TRUE))
  expect_no_error(ssn_split_predpts(mf04p, "pred1km", subset = ratio > 0.5,
                                    id_predpts = "RATIO_05", overwrite = TRUE))

  expect_no_error(ssn_get_stream_distmat(mf04p))

  expect_no_error(ssn_write(mf04p, path = paste0(tempdir(), "/tempSSN.ssn"),
                            overwrite = TRUE))

  expect_no_error(ssn_put_data(ssn_get_data(mf04p), mf04p))

})
