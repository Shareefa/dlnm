# Tests for the P1 fused cross-basis Rust kernel
#
# These tests compare the Rust fused_crossbasis() output against the R
# implementation (Lag() + matrix multiply) for correctness.

# Helper: compute cross-basis the R way (reference implementation)
crossbasis_r_reference <- function(basisvar, basislag, lag, group = NULL) {
  nv <- ncol(basisvar)
  nl <- ncol(basislag)
  n <- nrow(basisvar)
  cb <- matrix(0, nrow = n, ncol = nv * nl)
  for (v in seq_len(nv)) {
    mat <- as.matrix(tsModel::Lag(basisvar[, v], seqlag(lag), group = group))
    for (l in seq_len(nl)) {
      cb[, nl * (v - 1) + l] <- mat %*% basislag[, l]
    }
  }
  cb
}

# Helper: build cross-basis for a config using the full crossbasis() function
build_crossbasis_config <- function(dt, config, use_group = TRUE) {
  group_arg <- if (use_group) dt$city else NULL

  if (config == "C1") {
    crossbasis(dt$temp, lag = c(0, 15),
               argvar = list(fun = "lin"),
               arglag = list(fun = "poly", degree = 4),
               group = group_arg)
  } else if (config == "C2") {
    crossbasis(dt$temp, lag = c(0, 21),
               argvar = list(fun = "ns", df = 5),
               arglag = list(fun = "ns", df = 4),
               group = group_arg)
  } else if (config == "C3") {
    crossbasis(dt$temp, lag = c(0, 40),
               argvar = list(fun = "bs", df = 6),
               arglag = list(fun = "ns", df = 4),
               group = group_arg)
  } else if (config == "C4") {
    crossbasis(dt$temp, lag = c(0, 30),
               argvar = list(fun = "ps", df = 10),
               arglag = list(fun = "ps", df = 5),
               group = group_arg)
  } else if (config == "C5") {
    crossbasis(dt$temp, lag = c(0, 60),
               argvar = list(fun = "ps", df = 15),
               arglag = list(fun = "ps", df = 8),
               group = group_arg)
  }
}

# Load benchmark data (find project root relative to test dir)
pkg_root <- normalizePath(file.path(test_path(), "..", ".."))
dt_10mb <- readRDS(file.path(pkg_root, "benchmarks", "data", "scale_10mb.rds"))

# ---- Config Tests: Rust fused_crossbasis matches R for C1-C5 ----

test_that("P1 Rust fused_crossbasis matches R for C1 (lin x poly4, lag 0-15)", {
  lag <- c(0, 15)
  basisvar <- do.call("onebasis", list(x = dt_10mb$temp, fun = "lin"))
  basislag <- do.call("onebasis", list(x = seqlag(lag), fun = "poly", degree = 4, intercept = TRUE))

  r_result <- crossbasis_r_reference(basisvar, basislag, lag, group = dt_10mb$city)

  # Build group start/end indices (1-based for R, will be adjusted in Rust wrapper)
  group_factor <- dt_10mb$city
  group_rle <- rle(as.integer(group_factor))
  group_ends <- cumsum(group_rle$lengths)
  group_starts <- c(1L, group_ends[-length(group_ends)] + 1L)

  rust_result <- fused_crossbasis(
    basisvar, basislag,
    as.integer(lag[1]), as.integer(lag[2]),
    as.integer(group_starts), as.integer(group_ends)
  )

  expect_equal(as.numeric(rust_result), as.numeric(r_result), tolerance = 1e-10)
})

test_that("P1 Rust fused_crossbasis matches R for C2 (ns5 x ns4, lag 0-21)", {
  lag <- c(0, 21)
  basisvar <- do.call("onebasis", list(x = dt_10mb$temp, fun = "ns", df = 5))
  basislag <- do.call("onebasis", list(x = seqlag(lag), fun = "ns", df = 4, intercept = TRUE))

  r_result <- crossbasis_r_reference(basisvar, basislag, lag, group = dt_10mb$city)

  group_factor <- dt_10mb$city
  group_rle <- rle(as.integer(group_factor))
  group_ends <- cumsum(group_rle$lengths)
  group_starts <- c(1L, group_ends[-length(group_ends)] + 1L)

  rust_result <- fused_crossbasis(
    basisvar, basislag,
    as.integer(lag[1]), as.integer(lag[2]),
    as.integer(group_starts), as.integer(group_ends)
  )

  expect_equal(as.numeric(rust_result), as.numeric(r_result), tolerance = 1e-10)
})

test_that("P1 Rust fused_crossbasis matches R for C3 (bs6 x ns4, lag 0-40)", {
  lag <- c(0, 40)
  basisvar <- do.call("onebasis", list(x = dt_10mb$temp, fun = "bs", df = 6))
  basislag <- do.call("onebasis", list(x = seqlag(lag), fun = "ns", df = 4, intercept = TRUE))

  r_result <- crossbasis_r_reference(basisvar, basislag, lag, group = dt_10mb$city)

  group_factor <- dt_10mb$city
  group_rle <- rle(as.integer(group_factor))
  group_ends <- cumsum(group_rle$lengths)
  group_starts <- c(1L, group_ends[-length(group_ends)] + 1L)

  rust_result <- fused_crossbasis(
    basisvar, basislag,
    as.integer(lag[1]), as.integer(lag[2]),
    as.integer(group_starts), as.integer(group_ends)
  )

  expect_equal(as.numeric(rust_result), as.numeric(r_result), tolerance = 1e-10)
})

test_that("P1 Rust fused_crossbasis matches R for C4 (ps10 x ps5, lag 0-30)", {
  lag <- c(0, 30)
  basisvar <- do.call("onebasis", list(x = dt_10mb$temp, fun = "ps", df = 10))
  basislag <- do.call("onebasis", list(x = seqlag(lag), fun = "ps", df = 5, intercept = TRUE))

  r_result <- crossbasis_r_reference(basisvar, basislag, lag, group = dt_10mb$city)

  group_factor <- dt_10mb$city
  group_rle <- rle(as.integer(group_factor))
  group_ends <- cumsum(group_rle$lengths)
  group_starts <- c(1L, group_ends[-length(group_ends)] + 1L)

  rust_result <- fused_crossbasis(
    basisvar, basislag,
    as.integer(lag[1]), as.integer(lag[2]),
    as.integer(group_starts), as.integer(group_ends)
  )

  expect_equal(as.numeric(rust_result), as.numeric(r_result), tolerance = 1e-10)
})

test_that("P1 Rust fused_crossbasis matches R for C5 (ps15 x ps8, lag 0-60)", {
  lag <- c(0, 60)
  basisvar <- do.call("onebasis", list(x = dt_10mb$temp, fun = "ps", df = 15))
  basislag <- do.call("onebasis", list(x = seqlag(lag), fun = "ps", df = 8, intercept = TRUE))

  r_result <- crossbasis_r_reference(basisvar, basislag, lag, group = dt_10mb$city)

  group_factor <- dt_10mb$city
  group_rle <- rle(as.integer(group_factor))
  group_ends <- cumsum(group_rle$lengths)
  group_starts <- c(1L, group_ends[-length(group_ends)] + 1L)

  rust_result <- fused_crossbasis(
    basisvar, basislag,
    as.integer(lag[1]), as.integer(lag[2]),
    as.integer(group_starts), as.integer(group_ends)
  )

  expect_equal(as.numeric(rust_result), as.numeric(r_result), tolerance = 1e-10)
})

# ---- NA Handling ----

test_that("P1 Rust fused_crossbasis handles NA values correctly", {
  # Inject NAs into temp
  dt_na <- data.frame(dt_10mb)
  set.seed(123)
  na_idx <- sample(nrow(dt_na), 500)
  dt_na$temp[na_idx] <- NA

  lag <- c(0, 21)
  basisvar <- do.call("onebasis", list(x = dt_na$temp, fun = "ns", df = 5))
  basislag <- do.call("onebasis", list(x = seqlag(lag), fun = "ns", df = 4, intercept = TRUE))

  r_result <- crossbasis_r_reference(basisvar, basislag, lag, group = dt_na$city)

  group_factor <- dt_na$city
  group_rle <- rle(as.integer(group_factor))
  group_ends <- cumsum(group_rle$lengths)
  group_starts <- c(1L, group_ends[-length(group_ends)] + 1L)

  rust_result <- fused_crossbasis(
    basisvar, basislag,
    as.integer(lag[1]), as.integer(lag[2]),
    as.integer(group_starts), as.integer(group_ends)
  )

  # NA positions should match
  expect_equal(is.na(as.numeric(rust_result)), is.na(as.numeric(r_result)))
  # Non-NA values should match
  non_na <- !is.na(as.numeric(r_result))
  expect_equal(as.numeric(rust_result)[non_na], as.numeric(r_result)[non_na], tolerance = 1e-10)
})

# ---- Group Boundary Tests ----

test_that("P1 Rust fused_crossbasis respects group boundaries", {
  # Use a small subset with known group boundaries
  # Take first 3 cities (each 5114 rows)
  dt_small <- dt_10mb[dt_10mb$city %in% levels(dt_10mb$city)[1:3], ]
  dt_small$city <- droplevels(dt_small$city)

  lag <- c(0, 21)
  basisvar <- do.call("onebasis", list(x = dt_small$temp, fun = "ns", df = 5))
  basislag <- do.call("onebasis", list(x = seqlag(lag), fun = "ns", df = 4, intercept = TRUE))

  r_result <- crossbasis_r_reference(basisvar, basislag, lag, group = dt_small$city)

  group_factor <- dt_small$city
  group_rle <- rle(as.integer(group_factor))
  group_ends <- cumsum(group_rle$lengths)
  group_starts <- c(1L, group_ends[-length(group_ends)] + 1L)

  rust_result <- fused_crossbasis(
    basisvar, basislag,
    as.integer(lag[1]), as.integer(lag[2]),
    as.integer(group_starts), as.integer(group_ends)
  )

  # Check boundary rows specifically — first 21 rows of each city should have some NAs
  # due to lag window extending beyond city start
  city_boundaries <- group_starts
  for (start in city_boundaries) {
    # First row of each city with lag 0-21: rows 0-20 from city start should have NAs
    # Check that R and Rust agree on which rows are NA
    check_rows <- start:(min(start + 20, nrow(r_result)))
    expect_equal(
      is.na(rust_result[check_rows, ]),
      is.na(r_result[check_rows, ]),
      info = paste("Group boundary at row", start)
    )
  }

  # Overall equality
  expect_equal(as.numeric(rust_result), as.numeric(r_result), tolerance = 1e-10)
})

# ---- Single Group (No Grouping) ----

test_that("P1 Rust fused_crossbasis works with single group (no grouping)", {
  # Use data from a single city
  dt_single <- dt_10mb[dt_10mb$city == levels(dt_10mb$city)[1], ]

  lag <- c(0, 21)
  basisvar <- do.call("onebasis", list(x = dt_single$temp, fun = "ns", df = 5))
  basislag <- do.call("onebasis", list(x = seqlag(lag), fun = "ns", df = 4, intercept = TRUE))

  # R reference without group
  r_result <- crossbasis_r_reference(basisvar, basislag, lag, group = NULL)

  # Rust: single group = entire dataset
  group_starts <- 1L
  group_ends <- nrow(basisvar)

  rust_result <- fused_crossbasis(
    basisvar, basislag,
    as.integer(lag[1]), as.integer(lag[2]),
    as.integer(group_starts), as.integer(group_ends)
  )

  expect_equal(as.numeric(rust_result), as.numeric(r_result), tolerance = 1e-10)
})

# ---- Small Manual Test ----

test_that("P1 Rust fused_crossbasis produces correct output on tiny example", {
  # Tiny manual test for sanity
  set.seed(42)
  x <- rnorm(20)
  g <- factor(c(rep("a", 10), rep("b", 10)))

  basisvar <- do.call("onebasis", list(x = x, fun = "ns", df = 3))
  basislag <- do.call("onebasis", list(x = 0:5, fun = "ns", df = 2, intercept = TRUE))

  r_result <- crossbasis_r_reference(basisvar, basislag, c(0, 5), group = g)

  group_rle <- rle(as.integer(g))
  group_ends <- cumsum(group_rle$lengths)
  group_starts <- c(1L, group_ends[-length(group_ends)] + 1L)

  rust_result <- fused_crossbasis(
    basisvar, basislag,
    0L, 5L,
    as.integer(group_starts), as.integer(group_ends)
  )

  expect_equal(as.numeric(rust_result), as.numeric(r_result), tolerance = 1e-10)
})
