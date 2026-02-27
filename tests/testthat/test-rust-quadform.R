# Tests for the P2 quadratic form SE Rust functions
#
# These tests compare the Rust quad_form_se() and cumulative_quad_form_se()
# output against the R implementation for correctness.

# Helper: build a crosspred scenario and extract Xpred, coef, vcov
# Returns list(xpred_lagspec, xpred_cumul, vcov, predvar, predlag_spec, predlag_cumul)
build_crosspred_scenario <- function(dt, config, use_group = TRUE) {
  group_arg <- if (use_group) dt$city else NULL

  if (config == "C1") {
    cb <- crossbasis(dt$temp, lag = c(0, 15),
                     argvar = list(fun = "lin"),
                     arglag = list(fun = "poly", degree = 4),
                     group = group_arg)
  } else if (config == "C2") {
    cb <- crossbasis(dt$temp, lag = c(0, 21),
                     argvar = list(fun = "ns", df = 5),
                     arglag = list(fun = "ns", df = 4),
                     group = group_arg)
  }

  # Fit model
  mod <- glm(death ~ cb + ns(time, 7 * 14) + dow, data = dt, family = quasipoisson())

  # Extract coef and vcov for the cross-basis terms
  ind <- grep("cb", names(coef(mod)))
  beta <- coef(mod)[ind]
  V <- vcov(mod)[ind, ind]

  # Build prediction matrices
  lag <- attr(cb, "lag")
  at <- pretty(dt$temp, 20)
  predvar <- at
  predlag_spec <- seqlag(lag, 1)  # may be non-integer
  predlag_cumul <- seqlag(lag)     # integer lags for cumulative

  # Lag-specific Xpred
  Xpred_spec <- mkXpred("cb", cb, at, predvar, predlag_spec, cen = NULL)

  # Cumulative Xpred (integer lags)
  Xpred_cumul <- mkXpred("cb", cb, at, predvar, predlag_cumul, cen = NULL)

  list(
    xpred_spec = Xpred_spec,
    xpred_cumul = Xpred_cumul,
    coef = beta,
    vcov = V,
    predvar = predvar,
    predlag_spec = predlag_spec,
    predlag_cumul = predlag_cumul,
    lag = lag
  )
}

# Load benchmark data
pkg_root <- normalizePath(file.path(test_path(), "..", ".."))
dt_10mb <- readRDS(file.path(pkg_root, "benchmarks", "data", "scale_10mb.rds"))

# Use only first 3 cities for speed in tests
dt_small <- dt_10mb[dt_10mb$city %in% levels(dt_10mb$city)[1:3], ]
dt_small$city <- droplevels(dt_small$city)

# ---- Lag-Specific SE Tests ----

test_that("P2 quad_form_se matches R for C1 lag-specific SE", {
  scenario <- build_crosspred_scenario(dt_small, "C1")

  # R reference: rowSums((Xpred %*% V) * Xpred)
  r_se <- sqrt(pmax(0, rowSums((scenario$xpred_spec %*% scenario$vcov) * scenario$xpred_spec)))

  # Rust
  rust_se <- as.numeric(quad_form_se(scenario$xpred_spec, scenario$vcov))

  expect_equal(rust_se, r_se, tolerance = 1e-10)
})

test_that("P2 quad_form_se matches R for C2 lag-specific SE", {
  scenario <- build_crosspred_scenario(dt_small, "C2")

  # R reference
  r_se <- sqrt(pmax(0, rowSums((scenario$xpred_spec %*% scenario$vcov) * scenario$xpred_spec)))

  # Rust
  rust_se <- as.numeric(quad_form_se(scenario$xpred_spec, scenario$vcov))

  expect_equal(rust_se, r_se, tolerance = 1e-10)
})

# ---- Cumulative SE Tests ----

test_that("P2 cumulative_quad_form_se matches R for C1 cumulative SE", {
  scenario <- build_crosspred_scenario(dt_small, "C1")
  n_at <- length(scenario$predvar)
  n_lag <- length(scenario$predlag_cumul)
  Xpred <- scenario$xpred_cumul
  V <- scenario$vcov

  # R reference: cumulative loop
  Xpredall <- 0
  cumse_r <- matrix(0, n_at, n_lag)
  for (i in seq(n_lag)) {
    ind <- seq(n_at) + n_at * (i - 1)
    Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
    cumse_r[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% V) * Xpredall)))
  }

  # Rust
  rust_result <- cumulative_quad_form_se(Xpred, V, as.integer(n_at), as.integer(n_lag))
  rust_cumse <- rust_result[, 1:n_lag, drop = FALSE]

  expect_equal(as.numeric(rust_cumse), as.numeric(cumse_r), tolerance = 1e-10)
})

test_that("P2 cumulative_quad_form_se matches R for C2 cumulative SE", {
  scenario <- build_crosspred_scenario(dt_small, "C2")
  n_at <- length(scenario$predvar)
  n_lag <- length(scenario$predlag_cumul)
  Xpred <- scenario$xpred_cumul
  V <- scenario$vcov

  # R reference
  Xpredall <- 0
  cumse_r <- matrix(0, n_at, n_lag)
  for (i in seq(n_lag)) {
    ind <- seq(n_at) + n_at * (i - 1)
    Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
    cumse_r[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% V) * Xpredall)))
  }

  # Rust
  rust_result <- cumulative_quad_form_se(Xpred, V, as.integer(n_at), as.integer(n_lag))
  rust_cumse <- rust_result[, 1:n_lag, drop = FALSE]

  expect_equal(as.numeric(rust_cumse), as.numeric(cumse_r), tolerance = 1e-10)
})

# ---- Overall SE Tests ----

test_that("P2 overall SE matches R for C1 (final accumulation)", {
  scenario <- build_crosspred_scenario(dt_small, "C1")
  n_at <- length(scenario$predvar)
  n_lag <- length(scenario$predlag_cumul)
  Xpred <- scenario$xpred_cumul
  V <- scenario$vcov

  # R reference: overall SE = SE of fully accumulated Xpredall
  Xpredall <- 0
  for (i in seq(n_lag)) {
    ind <- seq(n_at) + n_at * (i - 1)
    Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
  }
  allse_r <- sqrt(pmax(0, rowSums((Xpredall %*% V) * Xpredall)))

  # Rust: overall SE is in the last column (n_lag + 1)
  rust_result <- cumulative_quad_form_se(Xpred, V, as.integer(n_at), as.integer(n_lag))
  rust_allse <- as.numeric(rust_result[, n_lag + 1])

  expect_equal(rust_allse, allse_r, tolerance = 1e-10)
})

test_that("P2 overall SE matches R for C2 (final accumulation)", {
  scenario <- build_crosspred_scenario(dt_small, "C2")
  n_at <- length(scenario$predvar)
  n_lag <- length(scenario$predlag_cumul)
  Xpred <- scenario$xpred_cumul
  V <- scenario$vcov

  # R reference
  Xpredall <- 0
  for (i in seq(n_lag)) {
    ind <- seq(n_at) + n_at * (i - 1)
    Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
  }
  allse_r <- sqrt(pmax(0, rowSums((Xpredall %*% V) * Xpredall)))

  # Rust
  rust_result <- cumulative_quad_form_se(Xpred, V, as.integer(n_at), as.integer(n_lag))
  rust_allse <- as.numeric(rust_result[, n_lag + 1])

  expect_equal(rust_allse, allse_r, tolerance = 1e-10)
})

# ---- pmax(0,...) Clamping Tests ----

test_that("P2 quad_form_se handles near-singular vcov (pmax clamping)", {
  # Create a near-singular vcov that could produce negative quadratic forms
  # due to floating-point arithmetic
  set.seed(123)
  p <- 5
  n <- 10

  # Create a rank-deficient vcov by constructing from fewer vectors than p
  A <- matrix(rnorm(p * 2), p, 2)  # rank 2 matrix
  V_singular <- A %*% t(A)  # rank-2 vcov (singular)

  # Add tiny noise to make it nearly singular
  V_near_singular <- V_singular + diag(p) * 1e-15

  X <- matrix(rnorm(n * p), n, p)

  # R reference (with pmax clamping)
  r_se <- sqrt(pmax(0, rowSums((X %*% V_near_singular) * X)))

  # Rust
  rust_se <- as.numeric(quad_form_se(X, V_near_singular))

  # Both should produce non-negative results (no NaN from sqrt of negative)
  expect_true(all(is.finite(rust_se)))
  expect_true(all(rust_se >= 0))
  expect_equal(rust_se, r_se, tolerance = 1e-10)
})

test_that("P2 cumulative_quad_form_se handles near-singular vcov", {
  set.seed(456)
  p <- 4
  n_at <- 3
  n_lag <- 5

  # Near-singular vcov
  A <- matrix(rnorm(p * 2), p, 2)
  V_near_singular <- A %*% t(A) + diag(p) * 1e-15

  Xpred <- matrix(rnorm(n_at * n_lag * p), n_at * n_lag, p)

  # R reference
  Xpredall <- 0
  cumse_r <- matrix(0, n_at, n_lag)
  for (i in seq(n_lag)) {
    ind <- seq(n_at) + n_at * (i - 1)
    Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
    cumse_r[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% V_near_singular) * Xpredall)))
  }

  # Rust
  rust_result <- cumulative_quad_form_se(Xpred, V_near_singular, as.integer(n_at), as.integer(n_lag))
  rust_cumse <- rust_result[, 1:n_lag, drop = FALSE]

  expect_true(all(is.finite(as.numeric(rust_cumse))))
  expect_true(all(as.numeric(rust_cumse) >= 0))
  expect_equal(as.numeric(rust_cumse), as.numeric(cumse_r), tolerance = 1e-10)
})

# ---- Small Manual Test ----

test_that("P2 quad_form_se correct on tiny manual example", {
  # Simple 2x2 case for verification
  X <- matrix(c(1, 2, 3, 4), 2, 2)  # column-major: row1=[1,3], row2=[2,4]
  V <- matrix(c(2, 1, 1, 3), 2, 2)  # column-major: V = [[2,1],[1,3]]

  # Manual: row1 = [1,3], x'Vx = 1*2*1 + 1*1*3 + 3*1*1 + 3*3*3 = 2+3+3+27=35
  # Actually: x'Vx = [1,3] %*% [[2,1],[1,3]] %*% [1,3]' = [1,3] %*% [5,10] = 5+30=35
  # sqrt(35) = 5.916...
  # row2 = [2,4], x'Vx = [2,4] %*% [[2,1],[1,3]] %*% [2,4]' = [2,4] %*% [8,14] = 16+56=72
  # sqrt(72) = 8.485...

  r_se <- sqrt(pmax(0, rowSums((X %*% V) * X)))
  rust_se <- as.numeric(quad_form_se(X, V))

  expect_equal(rust_se, r_se, tolerance = 1e-10)
  expect_equal(rust_se[1], sqrt(35), tolerance = 1e-10)
  expect_equal(rust_se[2], sqrt(72), tolerance = 1e-10)
})

test_that("P2 cumulative_quad_form_se correct on tiny manual example", {
  set.seed(99)
  n_at <- 2
  n_lag <- 3
  p <- 2

  X <- matrix(rnorm(n_at * n_lag * p), n_at * n_lag, p)
  V <- crossprod(matrix(rnorm(p * p), p, p))

  # R reference
  Xpredall <- 0
  cumse_r <- matrix(0, n_at, n_lag)
  for (i in seq(n_lag)) {
    ind <- seq(n_at) + n_at * (i - 1)
    Xpredall <- Xpredall + X[ind, , drop = FALSE]
    cumse_r[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% V) * Xpredall)))
  }
  allse_r <- sqrt(pmax(0, rowSums((Xpredall %*% V) * Xpredall)))

  # Rust
  rust_result <- cumulative_quad_form_se(X, V, as.integer(n_at), as.integer(n_lag))

  expect_equal(as.numeric(rust_result[, 1:n_lag]), as.numeric(cumse_r), tolerance = 1e-10)
  expect_equal(as.numeric(rust_result[, n_lag + 1]), allse_r, tolerance = 1e-10)
})

# ---- Full crosspred Equivalence Test ----

test_that("P2 functions produce results matching full crosspred for C1", {
  # Run actual crosspred and compare SE matrices
  cb <- crossbasis(dt_small$temp, lag = c(0, 15),
                   argvar = list(fun = "lin"),
                   arglag = list(fun = "poly", degree = 4),
                   group = dt_small$city)
  mod <- glm(death ~ cb + ns(time, 7 * 14) + dow, data = dt_small, family = quasipoisson())
  pred <- crosspred(cb, mod, cumul = TRUE)

  # Get the SE matrices from crosspred
  matse_r <- pred$matse      # lag-specific SE
  cumse_r <- pred$cumse      # cumulative SE
  allse_r <- pred$allse      # overall SE

  # Now compute using Rust functions with the same inputs
  ind <- grep("cb", names(coef(mod)))
  beta <- coef(mod)[ind]
  V <- vcov(mod)[ind, ind]

  lag <- attr(cb, "lag")
  at <- pred$predvar
  predlag_spec <- seqlag(lag, pred$bylag)
  predlag_cumul <- seqlag(lag)

  # Lag-specific
  Xpred_spec <- mkXpred("cb", cb, at, at, predlag_spec, cen = pred$cen)
  rust_matse_vec <- as.numeric(quad_form_se(Xpred_spec, V))
  rust_matse <- matrix(rust_matse_vec, length(at), length(predlag_spec))

  expect_equal(as.numeric(rust_matse), as.numeric(matse_r), tolerance = 1e-10)

  # Cumulative
  Xpred_cumul <- mkXpred("cb", cb, at, at, predlag_cumul, cen = pred$cen)
  n_at <- length(at)
  n_lag <- length(predlag_cumul)
  rust_cum_result <- cumulative_quad_form_se(Xpred_cumul, V, as.integer(n_at), as.integer(n_lag))

  expect_equal(as.numeric(rust_cum_result[, 1:n_lag]), as.numeric(cumse_r), tolerance = 1e-10)

  # Overall
  expect_equal(as.numeric(rust_cum_result[, n_lag + 1]), as.numeric(allse_r), tolerance = 1e-10)
})

test_that("P2 functions produce results matching full crosspred for C2", {
  cb <- crossbasis(dt_small$temp, lag = c(0, 21),
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4),
                   group = dt_small$city)
  mod <- glm(death ~ cb + ns(time, 7 * 14) + dow, data = dt_small, family = quasipoisson())
  pred <- crosspred(cb, mod, cumul = TRUE)

  matse_r <- pred$matse
  cumse_r <- pred$cumse
  allse_r <- pred$allse

  ind <- grep("cb", names(coef(mod)))
  beta <- coef(mod)[ind]
  V <- vcov(mod)[ind, ind]

  lag <- attr(cb, "lag")
  at <- pred$predvar
  predlag_spec <- seqlag(lag, pred$bylag)
  predlag_cumul <- seqlag(lag)

  # Lag-specific
  Xpred_spec <- mkXpred("cb", cb, at, at, predlag_spec, cen = pred$cen)
  rust_matse_vec <- as.numeric(quad_form_se(Xpred_spec, V))
  rust_matse <- matrix(rust_matse_vec, length(at), length(predlag_spec))

  expect_equal(as.numeric(rust_matse), as.numeric(matse_r), tolerance = 1e-10)

  # Cumulative
  Xpred_cumul <- mkXpred("cb", cb, at, at, predlag_cumul, cen = pred$cen)
  n_at <- length(at)
  n_lag <- length(predlag_cumul)
  rust_cum_result <- cumulative_quad_form_se(Xpred_cumul, V, as.integer(n_at), as.integer(n_lag))

  expect_equal(as.numeric(rust_cum_result[, 1:n_lag]), as.numeric(cumse_r), tolerance = 1e-10)

  # Overall
  expect_equal(as.numeric(rust_cum_result[, n_lag + 1]), as.numeric(allse_r), tolerance = 1e-10)
})
