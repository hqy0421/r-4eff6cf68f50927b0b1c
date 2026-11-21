#' LOO-PIT calibration data for discrete predictive distributions
#'
#' The `loo_pit_calibration()` function computes the data needed for
#' leave-one-out (LOO) probability integral transform (PIT) calibration plots
#' of discrete (e.g., count) predictive distributions. For discrete data the
#' PIT value of an observation is not unique, and calibration plots are
#' usually based on *randomized* PIT values that change every time they are
#' generated. `loo_pit_calibration()` is a deterministic alternative: the
#' probability mass of each PIT interval is distributed over the histogram
#' bins (Czado, Gneiting, and Held, 2009), yielding *expected* bin counts,
#' and the conditional variance of each bin count under randomization is
#' reported alongside. No random numbers are drawn, so repeated calls give
#' identical results.
#'
#' @export
#' @param y A vector of observed counts (non-negative integers), one per
#'   observation.
#' @param yrep An `S` by `N` matrix of draws from the predictive
#'   distribution (counts), where `S` is the number of posterior draws and
#'   `N = length(y)` is the number of observations.
#' @param log_lik An `S` by `N` matrix of pointwise log-likelihood values
#'   with the same dimensions as `yrep`, used to compute the LOO importance
#'   weights with [psis()].
#' @param group An optional grouping vector of length `N` (e.g., store ids).
#'   The PSIS smoothing and normalization are carried out separately within
#'   each group and all summaries are reported per group. If `NULL`
#'   (default), all observations form a single group.
#' @param breaks A numeric vector of bin edges for the PIT histogram. Must be
#'   strictly increasing, start at 0, and end at 1. Bins may have unequal
#'   widths. The default is 10 equal-width bins.
#' @param r_eff Vector of relative effective sample size estimates containing
#'   one element per observation, passed to [psis()]. See [psis()] for
#'   details. The default is 1.
#' @template cores
#'
#' @return An object of class `"loo_pit_calibration"`, which is a named list
#'   with the following components:
#' \describe{
#'   \item{`bins`}{
#'   A data frame with one row per group and bin and columns `group`, `bin`
#'   (bin number), `lower` and `upper` (bin edges), `expected` (the expected
#'   count: the sum over the group's observations of the proportions
#'   \eqn{p_{ij}} of each PIT interval overlapping the bin), and
#'   `var_randomization` (the conditional variance
#'   \eqn{\sum_i p_{ij}(1-p_{ij})} of the bin count if the PIT values were
#'   instead randomized uniformly within their intervals).
#'   }
#'   \item{`intervals`}{
#'   A data frame with one row per observation and columns `group`, `y`,
#'   `pit_lower` (\eqn{P(Y < y)}), `pit_upper` (\eqn{P(Y \leq y)}), and the
#'   PSIS diagnostics `pareto_k` and `n_eff`.
#'   }
#'   \item{`diagnostics`}{
#'   A data frame with one row per group: the number of observations, the
#'   largest Pareto \eqn{k}, the number of Pareto \eqn{k} estimates above the
#'   sample-size-specific [threshold][pareto-k-diagnostic], and the smallest
#'   PSIS effective sample size.
#'   }
#'   \item{`breaks`}{
#'   The bin edges used.
#'   }
#' }
#'
#' @details For observation \eqn{i} with normalized LOO weights
#'   \eqn{w_i^{(s)}} (obtained with [psis()] from `log_lik` within its
#'   group), the PIT interval is
#'   \deqn{\left[\hat{P}(Y < y_i),\; \hat{P}(Y \leq y_i)\right], \quad
#'   \hat{P}(Y \leq y) = \sum_s w_i^{(s)}\, I\!\left(y^{(s)} \leq y\right).}
#'   If the interval has positive width, its probability mass is split across
#'   the bins in proportion to the overlap of the interval with each bin. If
#'   the interval is a single point (no predictive draw equals \eqn{y_i}),
#'   the point is assigned to the single bin containing it, with bins treated
#'   as left-closed/right-open, except that a point at 1 falls in the last
#'   bin. Only the bins overlapping an interval are visited when
#'   accumulating, so the cost per observation does not grow with the number
#'   of bins, and unequal bin widths are handled exactly.
#'
#' @references
#' Czado, C., Gneiting, T., and Held, L. (2009). Predictive model assessment
#' for count data. *Biometrics*, 65(4), 1254--1261.
#' \doi{10.1111/j.1541-0420.2009.01191.x}
#'
#' @template loo-and-psis-references
#'
#' @seealso [psis()], [pareto-k-diagnostic], [E_loo()]
#'
#' @examples
#' # Two "stores" with different Poisson rates and many zero counts
#' set.seed(123)
#' S <- 500 # posterior draws
#' n <- 30  # observations per store
#' group <- rep(c("store_1", "store_2"), each = n)
#' rate <- rep(c(0.5, 2), each = n)
#' y <- rpois(2 * n, rate)
#'
#' # stand-in posterior draws of each store's rate, predictive draws, log-lik
#' rate_draws <- matrix(rgamma(S * 2 * n, shape = 100, rate = 100 / rate), nrow = S)
#' yrep <- matrix(rpois(S * 2 * n, rate_draws), nrow = S)
#' log_lik <- matrix(dpois(rep(y, each = S), rate_draws, log = TRUE), nrow = S)
#'
#' cal <- loo_pit_calibration(y, yrep, log_lik, group = group)
#' cal
#' head(cal$bins)
#' head(cal$intervals)
#'
#' # bins may have unequal widths
#' loo_pit_calibration(y, yrep, log_lik, group = group,
#'                     breaks = c(0, 0.25, 0.5, 0.75, 0.9, 1))
#'
loo_pit_calibration <- function(y,
                                yrep,
                                log_lik,
                                group = NULL,
                                breaks = seq(0, 1, by = 0.1),
                                r_eff = 1,
                                cores = getOption("mc.cores", 1)) {
  checkmate::assert_integerish(y, lower = 0, any.missing = FALSE, min.len = 1)
  N <- length(y)
  checkmate::assert_matrix(yrep, mode = "numeric", any.missing = FALSE, ncols = N)
  checkmate::assert_matrix(log_lik, mode = "numeric", any.missing = FALSE,
                           nrows = nrow(yrep), ncols = N)
  validate_log_lik(log_lik)
  breaks <- validate_pit_breaks(breaks)
  group <- prepare_pit_group(group, N)
  r_eff <- prepare_psis_r_eff(r_eff, len = N)
  cores <- loo_cores(cores)

  n_bins <- length(breaks) - 1L
  k_threshold <- ps_khat_threshold(nrow(yrep))
  idx_list <- split(seq_len(N), group)
  group_levels <- names(idx_list)

  bins_list <- vector("list", length(idx_list))
  intervals_list <- vector("list", length(idx_list))
  diagnostics_list <- vector("list", length(idx_list))

  for (g in seq_along(idx_list)) {
    idx <- idx_list[[g]]
    n_g <- length(idx)
    psis_g <- psis(-log_lik[, idx, drop = FALSE], r_eff = r_eff[idx], cores = cores)
    w_g <- weights(psis_g, log = FALSE)
    pit <- weighted_pit_intervals(y[idx], yrep[, idx, drop = FALSE], w_g)
    binned <- bin_pit_intervals(pit$lower, pit$upper, breaks)

    bins_list[[g]] <- data.frame(
      group = factor(group_levels[g], levels = group_levels),
      bin = seq_len(n_bins),
      lower = breaks[-length(breaks)],
      upper = breaks[-1L],
      expected = binned$expected,
      var_randomization = binned$var_randomization
    )
    intervals_list[[g]] <- data.frame(
      group = factor(rep.int(group_levels[g], n_g), levels = group_levels),
      y = y[idx],
      pit_lower = pit$lower,
      pit_upper = pit$upper,
      pareto_k = psis_g$diagnostics$pareto_k,
      n_eff = psis_g$diagnostics$n_eff
    )
    k_g <- psis_g$diagnostics$pareto_k
    diagnostics_list[[g]] <- data.frame(
      group = factor(group_levels[g], levels = group_levels),
      n_obs = n_g,
      max_pareto_k = max(k_g),
      n_pareto_k_high = sum(k_g > k_threshold),
      min_n_eff = min(psis_g$diagnostics$n_eff)
    )
  }

  structure(
    list(
      bins = do.call(rbind, bins_list),
      intervals = do.call(rbind, intervals_list),
      diagnostics = do.call(rbind, diagnostics_list),
      breaks = breaks
    ),
    class = c("loo_pit_calibration", "list")
  )
}

#' @rdname loo_pit_calibration
#' @export
#' @param x An object returned by [loo_pit_calibration()].
#' @param digits An integer passed to [base::round()] when printing.
#' @param ... Ignored.
#' @return `x`, invisibly.
#'
print.loo_pit_calibration <- function(x, digits = 2, ...) {
  cat("LOO-PIT calibration data for discrete predictive distributions\n")
  cat(
    "Observations:", sum(x$diagnostics$n_obs),
    "| Groups:", nrow(x$diagnostics),
    "| Bins:", length(x$breaks) - 1L,
    "\n\nPSIS diagnostics by group:\n"
  )
  print(x$diagnostics, digits = digits, row.names = FALSE)
  invisible(x)
}


# internal ----------------------------------------------------------------

#' PSIS-weighted PIT intervals for discrete observations
#'
#' @noRd
#' @param y Vector of `n` observed counts.
#' @param yrep `S` by `n` matrix of predictive draws.
#' @param w `S` by `n` matrix of normalized (columns sum to 1) weights.
#' @return A list with vectors `lower` (\eqn{P(Y < y)}) and `upper`
#'   (\eqn{P(Y \leq y)}), each of length `n`.
#'
weighted_pit_intervals <- function(y, yrep, w) {
  less <- sweep(yrep, MARGIN = 2, STATS = y, FUN = "<")
  equal <- sweep(yrep, MARGIN = 2, STATS = y, FUN = "==")
  lower <- colSums(w * less)
  upper <- lower + colSums(w * equal)
  # the values are weighted sums of indicators with normalized weights, so
  # they lie in [0, 1] up to floating-point error; clamp rare excursions
  lower <- pmin(pmax(lower, 0), 1)
  upper <- pmin(pmax(upper, 0), 1)
  list(lower = lower, upper = upper)
}

#' Distribute PIT intervals over bins and accumulate moments
#'
#' @noRd
#' @param lower,upper Vectors of PIT interval endpoints in `[0, 1]` with
#'   `lower <= upper`.
#' @param breaks Strictly increasing vector of `B + 1` bin edges from 0 to 1.
#' @return A list with vectors of length `B`: `expected` (sum over
#'   observations of the proportion of each interval overlapping each bin)
#'   and `var_randomization` (sum of `p * (1 - p)` over observations, the
#'   conditional variance of the bin counts if one PIT value per observation
#'   were drawn uniformly from its interval).
#'
#' @details Only the bins overlapping an interval are visited: the first and
#'   last such bin are located with [findInterval()] and the contribution of
#'   the bins in between is computed directly, so the cost per observation
#'   does not grow with the number of bins.
#'
bin_pit_intervals <- function(lower, upper, breaks) {
  n_bins <- length(breaks) - 1L
  # first and last bin touched by each interval; bins are left-closed and
  # right-open, with the rightmost bin also containing 1, which is exactly
  # the convention of findInterval(all.inside = TRUE)
  bin_lo <- findInterval(lower, breaks, all.inside = TRUE)
  bin_hi <- findInterval(upper, breaks, all.inside = TRUE)
  n_touch <- bin_hi - bin_lo + 1L

  # expand to one row per (observation, touched bin) pair
  obs_id <- rep.int(seq_along(lower), n_touch)
  bin_id <- rep.int(bin_lo, n_touch) + sequence(n_touch) - 1L

  lo <- lower[obs_id]
  hi <- upper[obs_id]
  width <- hi - lo
  overlap <- pmin(hi, breaks[bin_id + 1L]) - pmax(lo, breaks[bin_id])

  p <- numeric(length(obs_id))
  is_interval <- width > 0
  p[is_interval] <- overlap[is_interval] / width[is_interval]
  # a point interval places all its mass in the single located bin
  p[!is_interval] <- 1
  p <- pmin(pmax(p, 0), 1)

  agg <- rowsum(cbind(p, p * (1 - p)), group = bin_id)
  present <- as.integer(rownames(agg))
  expected <- numeric(n_bins)
  var_randomization <- numeric(n_bins)
  expected[present] <- agg[, 1]
  var_randomization[present] <- agg[, 2]
  list(expected = expected, var_randomization = pmax(var_randomization, 0))
}

#' Validate the `breaks` argument
#'
#' @noRd
#' @param breaks User's `breaks` argument.
#' @return `breaks`, invisibly, if no error is thrown.
#'
validate_pit_breaks <- function(breaks) {
  checkmate::assert_numeric(
    breaks,
    lower = 0,
    upper = 1,
    any.missing = FALSE,
    min.len = 2,
    sorted = TRUE,
    unique = TRUE
  )
  if (breaks[1] != 0 || breaks[length(breaks)] != 1) {
    stop("'breaks' must start at 0 and end at 1.", call. = FALSE)
  }
  invisible(breaks)
}

#' Validate and factorize the `group` argument
#'
#' @noRd
#' @param group User's `group` argument.
#' @param N Number of observations.
#' @return A factor of length `N` without unused levels.
#'
prepare_pit_group <- function(group, N) {
  if (is.null(group)) {
    group <- rep.int("all", N)
  }
  checkmate::assert_vector(group, len = N, any.missing = FALSE)
  group <- as.factor(group)
  if (anyNA(group)) {
    stop("'group' must not contain missing values.", call. = FALSE)
  }
  droplevels(group)
}
