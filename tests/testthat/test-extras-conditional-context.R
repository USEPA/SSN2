skip_on_cran()
skip_if_not(identical(Sys.getenv("SSN2_RUN_EXTRAS"), "true"),
            "set SSN2_RUN_EXTRAS = 'true' to run the extras suite")

test_that("random designs preserve tuple identities and base conditioning membership", {
  observed <- data.frame(g = factor(c("a:b", "a", "a:b", "a")),
    h = factor(c("c", "b:c", "c", "b:c")), p = factor(c("x", "x", "y", "y")),
    x = c(0, 2, -1, 3))
  newdata <- observed[c(2, 3), ]
  newdata$g <- c("new", "a:b")
  params <- c("1 | g:h" = .4, "x | g:h" = .2, "1 | g" = 0)
  full <- get_randcov_design(observed, params, ~p, newdata)
  data <- rbind(observed, newdata)
  group <- outer(as.character(data$g), as.character(data$g), "==") &
    outer(as.character(data$h), as.character(data$h), "==")
  partition <- outer(as.character(data$p), as.character(data$p), "==")
  reference <- (.4 + .2 * outer(data$x, data$x)) * group * partition
  Z <- rbind(full$observed, full$newdata)
  expect_equal(as.matrix(Z %*% Matrix::Diagonal(x = full$variance) %*% t(Z)),
               reference, ignore_attr = TRUE)
  expect_identical(anyDuplicated(rownames(full$info)), 0L)
  base <- get_randcov_design(observed, params, ~p, newdata, index = 1:2)
  expect_true(any(base$info$observed & !base$info$conditioned))
  expect_true(any(!base$info$observed & !base$info$conditioned))
  expect_equal(NROW(base$observed), 2)
  expect_true(all(c("level_tuple", "partition_tuple") %in% names(base$info)))
  expect_equal(NCOL(get_randcov_design(observed, NULL)$observed), 0)
  observed$g[1] <- NA
  expect_error(get_randcov_design(observed, params), "must not be missing")
})
