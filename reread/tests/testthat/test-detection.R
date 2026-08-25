test_that("reread() runs and returns a well-formed object", {
  set.seed(1)
  clean <- simulate_clean(12, 6, 250, seed = 1)
  dat   <- inject_careless(clean$data, prevalence = 0.2, seed = 2)
  fit   <- reread(dat$data)
  expect_s3_class(fit, "reread")
  expect_length(fit$flagged, nrow(dat$data))
  expect_type(fit$flagged, "logical")
  expect_true(all(c("rc", "longstring", "person_total", "eta", "prob", "flagged")
                  %in% names(fit$scores)))
  expect_false(anyNA(fit$scores$eta))
  expect_true(all(fit$prob >= 0 & fit$prob <= 1))
})

test_that("rc points the right way and the ensemble separates careless", {
  clean <- simulate_clean(12, 6, 300, seed = 3)
  dat   <- inject_careless(clean$data, prevalence = 0.25, seed = 4)
  fit   <- reread(dat$data)
  # low rc = careless -> negative correlation with the careless label
  expect_lt(cor(fit$scores$rc, dat$labels), -0.3)
  # high eta = careless -> positive
  expect_gt(cor(fit$scores$eta, dat$labels), 0.3)
})

test_that("clean data flags almost nobody", {
  clean <- simulate_clean(12, 6, 300, seed = 5)
  fit   <- reread(clean$data)
  expect_lt(fit$flagged_fraction, 0.03)
})

test_that("manual prevalence override flags the requested share", {
  clean <- simulate_clean(10, 6, 200, seed = 6)
  dat   <- inject_careless(clean$data, prevalence = 0.2, seed = 7)
  fit   <- reread(dat$data, prevalence = 0.15)
  expect_equal(fit$n_flagged, round(0.15 * 200), tolerance = 1)
})

test_that("person-total can be cut by the same calibration", {
  clean <- simulate_clean(10, 6, 250, seed = 11)
  dat   <- inject_careless(clean$data, prevalence = 0.2, seed = 12)
  ens   <- reread(dat$data)
  pt    <- reread(dat$data, score = "person_total")
  expect_identical(ens$score, "ensemble")
  expect_identical(pt$score, "person_total")
  # the alternative must actually decide, not silently fall back to the ensemble
  expect_gt(pt$n_flagged, 0)
  expect_false(identical(ens$flagged, pt$flagged))
  # the score column is the one that was cut
  expect_true(all(is.finite(pt$scores$eta)))
  # quantities defined only for the fitted ensemble are withheld, not faked
  expect_true(all(is.na(pt$scores$prob)))
  expect_true(all(is.na(pt$scores$eta_se)))
  expect_false(any(pt$scores$borderline))
  expect_true(all(is.finite(ens$scores$prob)))
})

test_that("the cut is equivariant to the scale of the score it is given", {
  clean <- simulate_clean(10, 6, 250, seed = 13)
  dat   <- inject_careless(clean$data, prevalence = 0.25, seed = 14)
  eta   <- ensemble_score(dat$data)$eta
  a <- auto_flag(eta)
  b <- auto_flag(3.7 * eta + 12)
  expect_identical(a$flagged, b$flagged)
})
