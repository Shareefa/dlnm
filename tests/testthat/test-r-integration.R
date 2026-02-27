# Tests for R integration: Rust dispatch + fallback in crossbasis() and crosspred()
#
# These tests verify:
# 1. crossbasis() dispatches to Rust and produces correct results
# 2. crossbasis() falls back to R when Rust is unavailable
# 3. crosspred() dispatches to Rust for SE computation
# 4. Partial availability: Rust P1 + R P2 fallback works
# 5. Full pipeline equivalence: crossbasis -> glm -> crosspred -> crossreduce

# Load benchmark data
pkg_root <- normalizePath(file.path(test_path(), "..", ".."))
dt_10mb <- readRDS(file.path(pkg_root, "benchmarks", "data", "scale_10mb.rds"))

# Use a reasonable subset for quick tests
dt_test <- dt_10mb[dt_10mb$city %in% levels(dt_10mb$city)[1:3], ]
dt_test$city <- droplevels(dt_test$city)

# Helper: temporarily replace a function in the dlnm namespace, run code, then restore
with_disabled_rust <- function(fn_names, code) {
  ns <- asNamespace("dlnm")
  originals <- list()
  for (fn_name in fn_names) {
    originals[[fn_name]] <- get(fn_name, envir = ns)
    unlockBinding(fn_name, ns)
    assign(fn_name, function(...) stop("Rust unavailable"), envir = ns)
  }
  on.exit({
    for (fn_name in fn_names) {
      unlockBinding(fn_name, ns)
      assign(fn_name, originals[[fn_name]], envir = ns)
      lockBinding(fn_name, ns)
    }
  })
  force(code)
}

# ---- crossbasis dispatch tests ----

test_that("crossbasis() dispatches to Rust and produces correct results for C2", {
  # Build crossbasis via the main function (which dispatches to Rust)
  cb_auto <- crossbasis(dt_test$temp, lag = c(0, 21),
                        argvar = list(fun = "ns", df = 5),
                        arglag = list(fun = "ns", df = 4),
                        group = dt_test$city)

  # Build crossbasis via pure R fallback
  lag <- c(0, 21)
  basisvar <- do.call("onebasis", list(x = dt_test$temp, fun = "ns", df = 5))
  basislag <- do.call("onebasis", list(x = seqlag(lag), fun = "ns", df = 4, intercept = TRUE))
  cb_r <- matrix(0, nrow = nrow(dt_test), ncol = ncol(basisvar) * ncol(basislag))
  for (v in seq_len(ncol(basisvar))) {
    mat <- as.matrix(tsModel::Lag(basisvar[, v], seqlag(lag), group = dt_test$city))
    for (l in seq_len(ncol(basislag))) {
      cb_r[, ncol(basislag) * (v - 1) + l] <- mat %*% basislag[, l]
    }
  }

  # The matrix values should match within tolerance
  expect_equal(as.numeric(unclass(cb_auto)), as.numeric(cb_r), tolerance = 1e-10)

  # Attributes should be preserved
  expect_true(inherits(cb_auto, "crossbasis"))
  expect_true(inherits(cb_auto, "matrix"))
  expect_equal(attr(cb_auto, "df"), c(ncol(basisvar), ncol(basislag)))
  expect_equal(attr(cb_auto, "lag"), lag)
  expect_true(!is.null(attr(cb_auto, "argvar")))
  expect_true(!is.null(attr(cb_auto, "arglag")))
  expect_true(!is.null(dimnames(cb_auto)))
})

test_that("crossbasis() R fallback works when Rust is unavailable", {
  cb_fallback <- with_disabled_rust("fused_crossbasis", {
    crossbasis(dt_test$temp, lag = c(0, 21),
               argvar = list(fun = "ns", df = 5),
               arglag = list(fun = "ns", df = 4),
               group = dt_test$city)
  })

  # Build reference (with Rust)
  cb_ref <- crossbasis(dt_test$temp, lag = c(0, 21),
                       argvar = list(fun = "ns", df = 5),
                       arglag = list(fun = "ns", df = 4),
                       group = dt_test$city)

  # Fallback should produce same results as Rust
  expect_equal(as.numeric(unclass(cb_fallback)), as.numeric(unclass(cb_ref)), tolerance = 1e-10)
  expect_true(inherits(cb_fallback, "crossbasis"))
})

test_that("crossbasis() works with group=NULL (single time series)", {
  dt_single <- dt_test[dt_test$city == levels(dt_test$city)[1], ]

  cb <- crossbasis(dt_single$temp, lag = c(0, 21),
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4))

  expect_true(inherits(cb, "crossbasis"))
  expect_equal(nrow(cb), nrow(dt_single))
  expect_equal(ncol(cb), 5 * 4)  # nv=5, nl=4
})

test_that("crossbasis() handles non-time-series (matrix) input without Rust", {
  # When x is a matrix (not time series), Rust dispatch is skipped
  # and R fallback is used
  dt_single <- dt_test[dt_test$city == levels(dt_test$city)[1], ]
  lag <- c(0, 21)

  # Create a lagged matrix manually
  x_mat <- as.matrix(tsModel::Lag(dt_single$temp, seqlag(lag)))

  cb <- crossbasis(x_mat, lag = lag,
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4))

  expect_true(inherits(cb, "crossbasis"))
  expect_equal(nrow(cb), nrow(dt_single))
})

# ---- crosspred dispatch tests ----

test_that("crosspred() dispatches to Rust for lag-specific SE", {
  cb <- crossbasis(dt_test$temp, lag = c(0, 21),
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4),
                   group = dt_test$city)
  mod <- glm(death ~ cb + ns(time, 7 * 14) + dow, data = dt_test, family = quasipoisson())

  # crosspred with Rust SE
  pred <- crosspred(cb, mod)

  # Compute reference SE via R
  ind <- grep("cb", names(coef(mod)))
  beta <- coef(mod)[ind]
  V <- vcov(mod)[ind, ind]
  lag <- attr(cb, "lag")
  at <- pred$predvar
  predlag <- seqlag(lag, pred$bylag)
  Xpred <- mkXpred("cb", cb, at, at, predlag, cen = pred$cen)

  r_matse <- matrix(sqrt(pmax(0, rowSums((Xpred %*% V) * Xpred))),
                    length(at), length(predlag))

  expect_equal(as.numeric(pred$matse), as.numeric(r_matse), tolerance = 1e-10)
})

test_that("crosspred() dispatches to Rust for cumulative SE", {
  cb <- crossbasis(dt_test$temp, lag = c(0, 21),
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4),
                   group = dt_test$city)
  mod <- glm(death ~ cb + ns(time, 7 * 14) + dow, data = dt_test, family = quasipoisson())

  pred <- crosspred(cb, mod, cumul = TRUE)

  # Compute reference cumulative SE via R
  ind <- grep("cb", names(coef(mod)))
  V <- vcov(mod)[ind, ind]
  lag <- attr(cb, "lag")
  at <- pred$predvar
  predlag_cumul <- seqlag(lag)
  Xpred <- mkXpred("cb", cb, at, at, predlag_cumul, cen = pred$cen)
  n_at <- length(at)
  n_lag <- length(predlag_cumul)

  Xpredall <- 0
  cumse_r <- matrix(0, n_at, n_lag)
  for (i in seq(n_lag)) {
    ind_r <- seq(n_at) + n_at * (i - 1)
    Xpredall <- Xpredall + Xpred[ind_r, , drop = FALSE]
    cumse_r[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% V) * Xpredall)))
  }
  allse_r <- sqrt(pmax(0, rowSums((Xpredall %*% V) * Xpredall)))

  expect_equal(as.numeric(pred$cumse), as.numeric(cumse_r), tolerance = 1e-10)
  expect_equal(as.numeric(pred$allse), as.numeric(allse_r), tolerance = 1e-10)
})

test_that("crosspred() R fallback works when Rust SE functions are unavailable", {
  cb <- crossbasis(dt_test$temp, lag = c(0, 21),
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4),
                   group = dt_test$city)
  mod <- glm(death ~ cb + ns(time, 7 * 14) + dow, data = dt_test, family = quasipoisson())

  # Get reference with Rust
  pred_ref <- crosspred(cb, mod, cumul = TRUE)

  # Disable Rust SE functions and run crosspred
  pred_fb <- with_disabled_rust(c("quad_form_se", "cumulative_quad_form_se"), {
    crosspred(cb, mod, cumul = TRUE)
  })

  expect_equal(as.numeric(pred_fb$matse), as.numeric(pred_ref$matse), tolerance = 1e-10)
  expect_equal(as.numeric(pred_fb$cumse), as.numeric(pred_ref$cumse), tolerance = 1e-10)
  expect_equal(as.numeric(pred_fb$allse), as.numeric(pred_ref$allse), tolerance = 1e-10)
})

# ---- Partial Availability Tests ----

test_that("Partial availability: Rust P1 crossbasis + R P2 SE fallback works", {
  # Rust crossbasis available, Rust SE functions disabled
  cb <- crossbasis(dt_test$temp, lag = c(0, 21),
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4),
                   group = dt_test$city)

  # Save reference with full Rust
  mod <- glm(death ~ cb + ns(time, 7 * 14) + dow, data = dt_test, family = quasipoisson())
  pred_full <- crosspred(cb, mod, cumul = TRUE)

  # Disable P2 (SE) but keep P1 (crossbasis) available
  # Build a new crossbasis (uses Rust P1), then crosspred falls back to R for SE
  pred_partial <- with_disabled_rust(c("quad_form_se", "cumulative_quad_form_se"), {
    cb2 <- crossbasis(dt_test$temp, lag = c(0, 21),
                      argvar = list(fun = "ns", df = 5),
                      arglag = list(fun = "ns", df = 4),
                      group = dt_test$city)
    mod2 <- glm(death ~ cb2 + ns(time, 7 * 14) + dow, data = dt_test, family = quasipoisson())
    crosspred(cb2, mod2, cumul = TRUE)
  })

  # Should produce same results
  expect_equal(as.numeric(pred_partial$matse), as.numeric(pred_full$matse), tolerance = 1e-10)
  expect_equal(as.numeric(pred_partial$cumse), as.numeric(pred_full$cumse), tolerance = 1e-10)
  expect_equal(as.numeric(pred_partial$allse), as.numeric(pred_full$allse), tolerance = 1e-10)
  expect_equal(as.numeric(pred_partial$allfit), as.numeric(pred_full$allfit), tolerance = 1e-10)
})

# ---- Full Pipeline Equivalence Tests ----

test_that("Full pipeline equivalence: crossbasis -> glm -> crosspred -> crossreduce for C2", {
  # ---- Pure R pipeline ----
  r_results <- with_disabled_rust(
    c("fused_crossbasis", "quad_form_se", "cumulative_quad_form_se"),
    {
      cb_r <- crossbasis(dt_test$temp, lag = c(0, 21),
                         argvar = list(fun = "ns", df = 5),
                         arglag = list(fun = "ns", df = 4),
                         group = dt_test$city)
      mod_r <- glm(death ~ cb_r + ns(time, 7 * 14) + dow, data = dt_test, family = quasipoisson())
      pred_r <- crosspred(cb_r, mod_r, cumul = TRUE)
      red_r <- crossreduce(cb_r, mod_r, type = "overall")
      list(cb = cb_r, mod = mod_r, pred = pred_r, red = red_r)
    }
  )

  # ---- Rust pipeline ----
  cb_rust <- crossbasis(dt_test$temp, lag = c(0, 21),
                        argvar = list(fun = "ns", df = 5),
                        arglag = list(fun = "ns", df = 4),
                        group = dt_test$city)
  mod_rust <- glm(death ~ cb_rust + ns(time, 7 * 14) + dow, data = dt_test, family = quasipoisson())
  pred_rust <- crosspred(cb_rust, mod_rust, cumul = TRUE)
  red_rust <- crossreduce(cb_rust, mod_rust, type = "overall")

  # ---- Compare 4 outputs ----
  # 1. Crossbasis matrix
  expect_equal(as.numeric(unclass(cb_rust)), as.numeric(unclass(r_results$cb)), tolerance = 1e-10,
               info = "Crossbasis matrix values differ")

  # 2. Model coefficients
  expect_equal(as.numeric(coef(mod_rust)), as.numeric(coef(r_results$mod)), tolerance = 1e-10,
               info = "Model coefficients differ")

  # 3. Predictions (crosspred)
  expect_equal(as.numeric(pred_rust$matfit), as.numeric(r_results$pred$matfit), tolerance = 1e-10,
               info = "Prediction matfit differs")
  expect_equal(as.numeric(pred_rust$matse), as.numeric(r_results$pred$matse), tolerance = 1e-10,
               info = "Prediction matse differs")
  expect_equal(as.numeric(pred_rust$allfit), as.numeric(r_results$pred$allfit), tolerance = 1e-10,
               info = "Prediction allfit differs")
  expect_equal(as.numeric(pred_rust$allse), as.numeric(r_results$pred$allse), tolerance = 1e-10,
               info = "Prediction allse differs")
  expect_equal(as.numeric(pred_rust$cumfit), as.numeric(r_results$pred$cumfit), tolerance = 1e-10,
               info = "Prediction cumfit differs")
  expect_equal(as.numeric(pred_rust$cumse), as.numeric(r_results$pred$cumse), tolerance = 1e-10,
               info = "Prediction cumse differs")

  # 4. Reduced estimates (crossreduce)
  expect_equal(as.numeric(red_rust$fit), as.numeric(r_results$red$fit), tolerance = 1e-10,
               info = "Crossreduce fit differs")
  expect_equal(as.numeric(red_rust$se), as.numeric(r_results$red$se), tolerance = 1e-10,
               info = "Crossreduce se differs")
  expect_equal(as.numeric(red_rust$coefficients), as.numeric(r_results$red$coefficients), tolerance = 1e-10,
               info = "Crossreduce coefficients differ")
})

# ---- Function signature preservation tests ----

test_that("crossbasis() return type and class are preserved", {
  cb <- crossbasis(dt_test$temp, lag = c(0, 21),
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4),
                   group = dt_test$city)

  expect_s3_class(cb, "crossbasis")
  expect_s3_class(cb, "matrix")
  expect_true(is.matrix(cb))

  # Check all expected attributes are present
  expect_true(!is.null(attr(cb, "df")))
  expect_true(!is.null(attr(cb, "range")))
  expect_true(!is.null(attr(cb, "lag")))
  expect_true(!is.null(attr(cb, "argvar")))
  expect_true(!is.null(attr(cb, "arglag")))
  expect_true(!is.null(attr(cb, "group")))
})

test_that("crosspred() return type is preserved", {
  cb <- crossbasis(dt_test$temp, lag = c(0, 21),
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4),
                   group = dt_test$city)
  mod <- glm(death ~ cb + ns(time, 7 * 14) + dow, data = dt_test, family = quasipoisson())
  pred <- crosspred(cb, mod, cumul = TRUE)

  expect_s3_class(pred, "crosspred")

  # Check required components
  expect_true(!is.null(pred$matfit))
  expect_true(!is.null(pred$matse))
  expect_true(!is.null(pred$allfit))
  expect_true(!is.null(pred$allse))
  expect_true(!is.null(pred$cumfit))
  expect_true(!is.null(pred$cumse))
  expect_true(!is.null(pred$matRRfit))
  expect_true(!is.null(pred$allRRfit))
  expect_true(!is.null(pred$cumRRfit))
})

# ---- crosspred without cumul ----

test_that("crosspred() works correctly without cumul=TRUE", {
  cb <- crossbasis(dt_test$temp, lag = c(0, 21),
                   argvar = list(fun = "ns", df = 5),
                   arglag = list(fun = "ns", df = 4),
                   group = dt_test$city)
  mod <- glm(death ~ cb + ns(time, 7 * 14) + dow, data = dt_test, family = quasipoisson())

  pred <- crosspred(cb, mod, cumul = FALSE)

  # Should have allfit/allse but not cumfit/cumse
  expect_true(!is.null(pred$allfit))
  expect_true(!is.null(pred$allse))
  expect_null(pred$cumfit)
  expect_null(pred$cumse)

  # Compare with R reference
  pred_fb <- with_disabled_rust(c("quad_form_se", "cumulative_quad_form_se"), {
    crosspred(cb, mod, cumul = FALSE)
  })

  expect_equal(as.numeric(pred$matse), as.numeric(pred_fb$matse), tolerance = 1e-10)
  expect_equal(as.numeric(pred$allse), as.numeric(pred_fb$allse), tolerance = 1e-10)
})
