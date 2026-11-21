set.seed(123)
S <- 100
n <- 10
rate <- 1.5
y <- rpois(n, rate)
rate_draws <- matrix(rgamma(S * n, shape = 50, rate = 50 / rate), nrow = S)
yrep <- matrix(rpois(S * n, rate_draws), nrow = S)
log_lik <- matrix(dpois(rep(y, each = S), rate_draws, log = TRUE), nrow = S)
group <- rep(c("s1", "s2"), each = n / 2)


test_that("weighted_pit_intervals computes weighted CDF steps", {
  yrep_ex <- matrix(c(0, 1, 1, 3,
                      2, 2, 0, 1), nrow = 4)
  w <- matrix(rep(c(0.1, 0.2, 0.3, 0.4), 2), nrow = 4)
  out <- weighted_pit_intervals(c(1, 2), yrep_ex, w)
  expect_equal(out$lower, c(0.1, 0.7))
  expect_equal(out$upper, c(0.6, 1.0))
})

test_that("weighted_pit_intervals gives point interval when no draw equals y", {
  yrep_ex <- matrix(c(1, 2, 3, 4), nrow = 4)
  w <- matrix(1 / 4, nrow = 4, ncol = 1)
  out <- weighted_pit_intervals(5, yrep_ex, w)
  expect_equal(out$lower, 1)
  expect_equal(out$upper, 1)
  out <- weighted_pit_intervals(0, yrep_ex, w)
  expect_equal(out$lower, 0)
  expect_equal(out$upper, 0)
})

test_that("bin_pit_intervals distributes interval mass by overlap proportion", {
  br <- c(0, 0.25, 0.5, 0.75, 1)
  out <- bin_pit_intervals(lower = c(0, 0.25), upper = c(0.5, 0.5), breaks = br)
  expect_equal(out$expected, c(0.5, 1.5, 0, 0))
  expect_equal(out$var_randomization, c(0.25, 0.25, 0, 0))
})

test_that("bin_pit_intervals handles unequal bin widths", {
  br <- c(0, 0.3, 1)
  out <- bin_pit_intervals(lower = c(0, 0.25), upper = c(0.5, 0.5), breaks = br)
  expect_equal(out$expected, c(0.6 + 0.2, 0.4 + 0.8))
  expect_equal(out$var_randomization, c(0.6 * 0.4 + 0.2 * 0.8,
                                        0.4 * 0.6 + 0.8 * 0.2))
})

test_that("bin_pit_intervals assigns point intervals left-closed/right-open, 1 in last bin", {
  br <- c(0, 0.25, 0.5, 1)
  out <- bin_pit_intervals(lower = c(0, 0.25, 1), upper = c(0, 0.25, 1), breaks = br)
  expect_equal(out$expected, c(1, 1, 1))
  expect_equal(out$var_randomization, c(0, 0, 0))
})

test_that("bin_pit_intervals gives zero-width overlap at internal breaks", {
  out <- bin_pit_intervals(lower = 0, upper = 0.3, breaks = c(0, 0.3, 1))
  expect_equal(out$expected, c(1, 0))
  expect_equal(out$var_randomization, c(0, 0))
})

test_that("bin_pit_intervals matches a brute-force all-bin scan", {
  set.seed(7)
  n_i <- 200
  lo <- runif(n_i)
  hi <- pmin(lo + runif(n_i, 0, 0.3), 1)
  point_idx <- sample(n_i, 20)
  hi[point_idx] <- lo[point_idx]
  br <- c(0, 0.13, 0.31, 0.55, 0.72, 0.9, 1)
  fast <- bin_pit_intervals(lo, hi, br)

  B <- length(br) - 1
  E <- V <- numeric(B)
  for (i in seq_len(n_i)) {
    w <- hi[i] - lo[i]
    if (w == 0) {
      j <- findInterval(lo[i], br, all.inside = TRUE)
      E[j] <- E[j] + 1
    } else {
      for (j in seq_len(B)) {
        p <- max(0, min(hi[i], br[j + 1]) - max(lo[i], br[j])) / w
        E[j] <- E[j] + p
        V[j] <- V[j] + p * (1 - p)
      }
    }
  }
  expect_equal(fast$expected, E, tolerance = 1e-12)
  expect_equal(fast$var_randomization, V, tolerance = 1e-12)
})

test_that("loo_pit_calibration returns the expected structure", {
  cal <- suppressWarnings(loo_pit_calibration(y, yrep, log_lik, group = group))
  expect_s3_class(cal, "loo_pit_calibration")
  expect_named(cal, c("bins", "intervals", "diagnostics", "breaks"))

  expect_equal(nrow(cal$bins), 2 * 10)
  expect_named(cal$bins,
               c("group", "bin", "lower", "upper", "expected", "var_randomization"))
  expect_equal(cal$bins$lower, rep(seq(0, 0.9, by = 0.1), 2))
  expect_equal(cal$bins$upper, rep(seq(0.1, 1, by = 0.1), 2))

  # expected counts sum to the number of observations in each group
  expect_equal(as.numeric(tapply(cal$bins$expected, cal$bins$group, sum)),
               c(5, 5), tolerance = 1e-10)
  # randomization variance is bounded by the expected count
  expect_true(all(cal$bins$var_randomization >= 0))
  expect_true(all(cal$bins$var_randomization <= cal$bins$expected + 1e-12))

  expect_equal(nrow(cal$intervals), n)
  expect_named(cal$intervals,
               c("group", "y", "pit_lower", "pit_upper", "pareto_k", "n_eff"))
  expect_equal(cal$intervals$y, y)
  expect_true(all(cal$intervals$pit_lower >= 0 & cal$intervals$pit_lower <= 1))
  expect_true(all(cal$intervals$pit_upper >= cal$intervals$pit_lower))
  expect_true(all(cal$intervals$pit_upper <= 1))

  expect_equal(nrow(cal$diagnostics), 2)
  expect_equal(as.character(cal$diagnostics$group), c("s1", "s2"))
  expect_equal(cal$diagnostics$n_obs, c(5, 5))
})

test_that("loo_pit_calibration intervals match manual PSIS-weighted computation", {
  cal <- suppressWarnings(loo_pit_calibration(y, yrep, log_lik, group = group))
  for (g in 1:2) {
    idx <- ((g - 1) * 5 + 1):(g * 5)
    psis_g <- suppressWarnings(psis(-log_lik[, idx]))
    w_g <- weights(psis_g, log = FALSE)
    lo <- colSums(w_g * sweep(yrep[, idx], 2, y[idx], "<"))
    hi <- lo + colSums(w_g * sweep(yrep[, idx], 2, y[idx], "=="))
    expect_equal(cal$intervals$pit_lower[idx], lo)
    expect_equal(cal$intervals$pit_upper[idx], hi)
    expect_equal(cal$intervals$pareto_k[idx], psis_g$diagnostics$pareto_k)
  }
})

test_that("loo_pit_calibration is deterministic and draws no random numbers", {
  set.seed(42)
  cal_a <- suppressWarnings(loo_pit_calibration(y, yrep, log_lik, group = group))
  seed_after_first_call <- .Random.seed
  cal_b <- suppressWarnings(loo_pit_calibration(y, yrep, log_lik, group = group))
  expect_identical(cal_a, cal_b)
  expect_identical(.Random.seed, seed_after_first_call)
})

test_that("loo_pit_calibration works without groups and with unequal breaks", {
  br <- c(0, 0.2, 0.5, 0.8, 0.95, 1)
  cal <- suppressWarnings(loo_pit_calibration(y, yrep, log_lik, breaks = br))
  expect_equal(nrow(cal$diagnostics), 1)
  expect_equal(as.character(cal$diagnostics$group), "all")
  expect_equal(nrow(cal$bins), length(br) - 1)
  expect_equal(sum(cal$bins$expected), n, tolerance = 1e-10)
  expect_equal(cal$breaks, br)
})

test_that("loo_pit_calibration validates its inputs", {
  expect_error(loo_pit_calibration(c(-1, y[-1]), yrep, log_lik), "y")
  expect_error(loo_pit_calibration(y + 0.5, yrep, log_lik), "y")
  expect_error(loo_pit_calibration(y, yrep[, -1], log_lik), "yrep")
  expect_error(loo_pit_calibration(y, yrep, log_lik[, -1]), "log_lik")
  ll_na <- log_lik
  ll_na[1, 1] <- NA
  expect_error(loo_pit_calibration(y, yrep, ll_na), "missing")
  expect_error(loo_pit_calibration(y, yrep, log_lik, group = rep("a", 3)), "group")
  group_na <- group
  group_na[1] <- NA
  expect_error(loo_pit_calibration(y, yrep, log_lik, group = group_na), "group")
  expect_error(loo_pit_calibration(y, yrep, log_lik, breaks = c(0.1, 0.5, 1)),
               "start at 0 and end at 1")
  expect_error(loo_pit_calibration(y, yrep, log_lik, breaks = c(0, 0.5, 0.9)),
               "start at 0 and end at 1")
  expect_error(loo_pit_calibration(y, yrep, log_lik, breaks = c(0, 0.5, 0.4, 1)),
               "sorted")
  expect_error(loo_pit_calibration(y, yrep, log_lik, breaks = c(0, 0.5, 0.5, 1)),
               "duplicated")
  expect_error(loo_pit_calibration(y, yrep, log_lik, breaks = 0.5),
               "breaks")
})

test_that("print.loo_pit_calibration works", {
  cal <- suppressWarnings(loo_pit_calibration(y, yrep, log_lik, group = group))
  expect_output(print(cal), "LOO-PIT calibration data for discrete predictive distributions")
  expect_output(print(cal), "PSIS diagnostics by group")
  expect_invisible(print(cal))
})
