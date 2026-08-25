test_that("proportion rescaling is a no-op on a homogeneous scale", {
  clean <- simulate_clean(10, 6, 200, seed = 8)
  a <- rc_index(clean$data, rescale = "proportion", seed = 1)
  b <- rc_index(clean$data, rescale = "none", seed = 1)
  # correlation matrix is unchanged, so rc is identical up to permutation seed
  expect_gt(cor(as.numeric(a), as.numeric(b)), 0.999)
})

test_that("reverse-keyed items do not collapse rc", {
  clean <- simulate_clean(12, 6, 250, seed = 9)
  # rc must still separate careless despite ~40% reverse-keyed items
  dat <- inject_careless(clean$data, prevalence = 0.25, seed = 10)
  rc  <- rc_index(dat$data)
  expect_lt(cor(as.numeric(rc), dat$labels), -0.3)
})

test_that("straight-liners (degenerate) get the most-careless rc", {
  clean <- simulate_clean(10, 6, 200, seed = 11)
  d <- clean$data
  d[1, ] <- 3L  # a pure straight-liner
  rc <- rc_index(d)
  expect_true(is.finite(rc[1]))
  expect_equal(as.numeric(rc[1]), min(as.numeric(rc)))
})

test_that("gauss_mixture recovers a known mixing proportion", {
  set.seed(12)
  x <- c(stats::rnorm(170, 0, 1), stats::rnorm(30, 5, 1))
  g <- gauss_mixture(x)
  expect_gt(g$bic1 - g$bic2, 2)          # two components clearly better
  expect_equal(g$pi, 0.15, tolerance = 0.05)
})

test_that("auto_flag flags only the far tail of unimodal (clean) scores", {
  set.seed(13)
  eta <- stats::rnorm(300, -1, 0.8)      # one mode, no careless component
  fl <- auto_flag(eta)
  # The cut is unconditional since the structure gate was removed, so a clean
  # sample is not guaranteed to yield nothing: it yields its own 2.5-sigma tail,
  # which for a normal score is about 0.6% of respondents. What must hold is that
  # the flagged set stays a tail rather than becoming a subpopulation, and that
  # the two-component diagnostic still reports the absence of one.
  expect_false(fl$two_component)
  expect_lt(fl$flagged_fraction, 0.03)
})

test_that("benchmark ranks the ensemble first", {
  skip_if_not_installed("pROC")
  clean <- simulate_clean(12, 6, 300, seed = 14)
  dat   <- inject_careless(clean$data, prevalence = 0.25, seed = 15)
  b <- benchmark_indices(dat$data, dat$labels)
  expect_equal(b$method[1], "reread")
  expect_gt(b$auc[b$method == "reread"], b$auc[b$method == "mahalanobis_d2"])
})
