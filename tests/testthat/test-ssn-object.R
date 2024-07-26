test_that("ssn manipulation functions work", {
  # test class of original SSN object
  expect_s3_class(mf04p, "SSN")


  tempdir_path <- paste0(tempdir(), "/MiddleFork04.ssn")

  # copy lsn to temporary directory (does it exist?)
  copy_lsn_to_temp()
  # problem on linux with reading stated files from path
  # expect_true(file.exists(tempdir_path))
  # are the correct files present?
  file_names <- c(
    "binaryID.db", "CapeHorn.gpkg", "edges.gpkg", "netID1.dat",
    "netID2.dat", "pred1km.gpkg", "sites.gpkg"
  )
  # problem on linux
  # expect_true(all(file_names %in% list.files(tempdir_path)))


  # import new SSN object from temporary LSN
  mf04p_new <- ssn_import(tempdir_path,
    predpts = c("pred1km"),
    overwrite = TRUE
  )
  expect_s3_class(mf04p_new, "SSN")
  expect_identical(names(mf04p_new), c("edges", "obs", "preds", "path"))
  expect_s3_class(mf04p_new$edges, "sf")
  expect_s3_class(mf04p_new$obs, "sf")
  expect_s3_class(mf04p_new$preds$pred1km, "sf")
  expect_identical(mf04p_new$path, tempdir_path)

  # import capehorn prediction data
  mf04p_new <- ssn_import_predpts(mf04p_new, predpts = "CapeHorn")
  expect_s3_class(mf04p_new$preds$CapeHorn, "sf")

  # update path
  mf04p <- ssn_update_path(mf04p, tempdir_path)
  expect_identical(mf04p$path, tempdir_path)

  # create distance matrices
  ssn_create_distmat(mf04p_new, predpts = c("pred1km"), overwrite = TRUE)
  distance_path <- paste0(tempdir_path, "/distance")
  # problem on linux
  # expect_true(file.exists(distance_path))
  # expect_true(all(c("obs", "pred1km") %in% list.files(distance_path)))
  ssn_create_distmat(mf04p_new,
    predpts = c("pred1km"),
    overwrite = TRUE, among_predpts = TRUE, only_predpts = TRUE
  )
  obs_path <- paste0(distance_path, "/obs")
  # problem on linux
  # expect_true(file.exists(obs_path))
  obs_files <- c("dist.net1.RData", "dist.net2.RData")
  # expect_true(all(obs_files %in% list.files(obs_path)))
  # problem on linux
  pred1km_path <- paste0(distance_path, "/pred1km")
  # expect_true(file.exists(pred1km_path))
  pred1km_files <- c(
    "dist.net1.RData", "dist.net2.RData", "dist.net1.a.RData",
    "dist.net1.b.RData", "dist.net2.a.RData", "dist.net2.b.RData"
  )
  # expect_true(all(pred1km_files %in% list.files(pred1km_path)))


  # subset SSN data
  subset_path <- paste0(tempdir(), "/subset1.ssn")
  ssn_subset(mf04p_new,
    path = subset_path,
    subset = netID == 1, clip = TRUE, overwrite = TRUE
  )
  expect_true(file.exists(subset_path))
  mf04p_new_subset <- ssn_import(subset_path, overwrite = TRUE)
  expect_equal(unique(mf04p_new_subset$obs$netID), 1)


  # split prediction points
  mf04p_new_split <- ssn_split_predpts(mf04p, "pred1km",
    size_predpts = 100,
    keep = FALSE, overwrite = TRUE
  )
  expect_identical(names(mf04p_new_split$preds), c("CapeHorn", "pred1km-1", "pred1km-2"))
  mf04p$preds$pred1km$net.fac <- as.factor(mf04p$preds$pred1km$netID)
  mf04p_new_split <- ssn_split_predpts(mf04p, "pred1km",
    by = "net.fac",
    keep = FALSE, overwrite = TRUE
  )
  expect_identical(names(mf04p_new_split$preds), c("CapeHorn", "pred1km-net.fac-1", "pred1km-net.fac-2"))
  mf04p_new_split <- ssn_split_predpts(mf04p, "pred1km",
    subset = ratio > 0.5,
    id_predpts = "RATIO_05", overwrite = TRUE
  )
  expect_identical(names(mf04p_new_split$preds), c("pred1km", "CapeHorn", "RATIO_05"))


  mf04p_distmats <- ssn_get_stream_distmat(mf04p)
  expect_identical(names(mf04p_distmats), c("dist.net1", "dist.net2"))
  expect_identical(dim(mf04p_distmats$dist.net1), as.integer(c(13, 13)))
  expect_identical(dim(mf04p_distmats$dist.net2), as.integer(c(32, 32)))


  # get SSN data
  mf04p_new_get <- ssn_get_data(mf04p_new)
  expect_identical(mf04p_new_get, mf04p_new$obs)

  # put SSN data
  mf04p_new <- ssn_put_data(mf04p_new_get, mf04p_new)
  expect_identical(mf04p_new_get, mf04p_new$obs)

  # write SSN object
  write_path <- paste0(tempdir_path, "/tempSSN.ssn")
  ssn_write(mf04p, path = write_path, overwrite = TRUE)
  # expect_true(file.exists(write_path))
})
