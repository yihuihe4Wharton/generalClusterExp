##############################################################
# Estimand assembly and the variance engine.
#
# Given the per-regime identification weights (weights1, weights0) produced by
# one of the weight constructors, this module
#   1. enumerates the well-defined population average potential outcomes
#      (marginal and non-marginal), Section 3 of the companion manuscript;
#   2. forms each Hajek LW point estimate and its per-unit contribution column
#      \hat V_ij = u_i * b_ij * (Y_ij - \hat mu);
#   3. assembles the d x d HAC covariance \hat Sigma over all estimands using
#      the kernel rule appropriate to the (method, design) pair;
#   4. reads off variances of every potential outcome and of the direct /
#      indirect / total / overall causal effects (contrasts of the estimands).
##############################################################

#' Hajek linear-weighted estimate and contribution column
#'
#' Computes the Hajek LW point estimate \eqn{\hat\mu = \sum_i u_i b_i Y_i /
#' \sum_i u_i b_i} for one estimand and its per-unit centred contribution
#' column \eqn{V_i = u_i b_i (Y_i - \hat\mu)}.
#'
#' NOTE: following the definition of \eqn{\hat\Sigma^K} in Section 5 of the
#' companion manuscript, the contribution is NOT divided by the Hajek
#' denominator \eqn{D = \sum_i u_i b_i} (which converges to 1; see the
#' estimator definition in Section 4). The point estimate itself is the
#' properly normalised Hajek ratio. If the finite-sample denominator D departs
#' noticeably from 1, dividing V by D yields a sharper variance.
#'
#' @param b per-unit identification weight column.
#' @param Y observed outcomes.
#' @param u per-unit weights.
#' @return A list with the point estimate \code{mu}, the contribution column
#'   \code{V}, and a flag \code{ok} (FALSE when the denominator is 0 or
#'   non-finite).
#' @keywords internal
#' @examples
#' set.seed(1)
#' Y <- rnorm(10); b <- runif(10); u <- rep(1 / 10, 10)
#' est <- generalClusterExp:::.lw_estimand(b, Y, u)
#' est$mu               # Hajek ratio sum(u * b * Y) / sum(u * b)
#' sum(est$V)           # centred contributions sum to 0
.lw_estimand <- function(b, Y, u) {
  den <- sum(u * b)
  if (!is.finite(den) || den == 0) {
    return(list(mu = NA_real_, V = rep(0, length(Y)), ok = FALSE))
  }
  mu <- sum(u * b * Y) / den
  list(mu = mu, V = u * b * (Y - mu), ok = TRUE)
}

#' Canonical potential-outcome specifications
#'
#' Canonical list of the six population average potential outcomes
#' (Section 3). Each entry has \code{regime} (1 for phi_1 / weights1, 0 for
#' phi_0 / weights0), \code{w} (\code{NA} for the marginal
#' \eqn{\bar Y(\phi)}; 0/1 for the non-marginal \eqn{\bar Y(w, \phi)}) and a
#' display \code{label}.
#'
#' @return A named list of specification lists.
#' @keywords internal
#' @examples
#' specs <- generalClusterExp:::.po_specs()
#' names(specs)
#' vapply(specs, function(s) s$label, character(1))
.po_specs <- function() {
  list(
    Ybar_phi1     = list(regime = 1, w = NA_integer_, label = "Ybar(phi1)"),
    Ybar_phi0     = list(regime = 0, w = NA_integer_, label = "Ybar(phi0)"),
    Ybar_w1_phi1  = list(regime = 1, w = 1L,          label = "Ybar(1,phi1)"),
    Ybar_w0_phi1  = list(regime = 1, w = 0L,          label = "Ybar(0,phi1)"),
    Ybar_w1_phi0  = list(regime = 0, w = 1L,          label = "Ybar(1,phi0)"),
    Ybar_w0_phi0  = list(regime = 0, w = 0L,          label = "Ybar(0,phi0)")
  )
}

#' Causal-effect contrast specifications
#'
#' Causal-effect contrasts over the potential-outcome keys (Section 3):
#' \itemize{
#'   \item overall \eqn{CE^O = \bar Y(\phi_1) - \bar Y(\phi_0)}
#'   \item direct \eqn{CE^D(\phi_c) = \bar Y(1, \phi_c) - \bar Y(0, \phi_c)}
#'   \item indirect \eqn{CE^I(w) = \bar Y(w, \phi_1) - \bar Y(w, \phi_0)}
#'   \item total \eqn{CE^T = \bar Y(1, \phi_1) - \bar Y(0, \phi_0)}
#' }
#'
#' @return A named list; each entry has the potential-outcome keys \code{plus}
#'   and \code{minus} and a display \code{label}.
#' @keywords internal
#' @examples
#' effs <- generalClusterExp:::.effect_specs()
#' names(effs)
#' effs$CE_total
.effect_specs <- function() {
  list(
    CE_overall      = list(plus = "Ybar_phi1",    minus = "Ybar_phi0",    label = "Overall CE^O(phi1,phi0)"),
    CE_direct_phi1  = list(plus = "Ybar_w1_phi1", minus = "Ybar_w0_phi1", label = "Direct CE^D(phi1)"),
    CE_direct_phi0  = list(plus = "Ybar_w1_phi0", minus = "Ybar_w0_phi0", label = "Direct CE^D(phi0)"),
    CE_indirect_w1  = list(plus = "Ybar_w1_phi1", minus = "Ybar_w1_phi0", label = "Indirect CE^I(1;phi1,phi0)"),
    CE_indirect_w0  = list(plus = "Ybar_w0_phi1", minus = "Ybar_w0_phi0", label = "Indirect CE^I(0;phi1,phi0)"),
    CE_total        = list(plus = "Ybar_w1_phi1", minus = "Ybar_w0_phi0", label = "Total CE^T(phi1,phi0)")
  )
}

#' Well-definedness of a potential outcome
#'
#' Decides whether a potential outcome is well-defined from the regime
#' own-treatment probabilities. Marginal outcomes are always well-defined;
#' \eqn{\bar Y(w, \phi)} needs \eqn{P_\phi(W = w) > 0} for every contributing
#' unit.
#'
#' @param spec one specification from \code{\link{.po_specs}}.
#' @param p1_own,p0_own per-unit own-treatment probabilities under phi_1 and
#'   phi_0 (see \code{\link{.own_treatment_prob1}}).
#' @return \code{TRUE} or \code{FALSE}.
#' @keywords internal
#' @examples
#' spec <- generalClusterExp:::.po_specs()$Ybar_w1_phi0   # Ybar(1, phi0)
#' p1 <- rep(0.5, 4)
#' generalClusterExp:::.is_well_defined(spec, p1, p0_own = rep(0.2, 4))
#' ## not defined when phi_0 never treats a unit:
#' generalClusterExp:::.is_well_defined(spec, p1, p0_own = rep(0, 4))
.is_well_defined <- function(spec, p1_own, p0_own) {
  if (is.na(spec$w)) return(TRUE)
  p_own <- if (spec$regime == 1) p1_own else p0_own
  pw <- if (spec$w == 1L) p_own else (1 - p_own)
  all(pw > 0)
}

#' Per-unit identification weight column for one estimand
#'
#' Builds the per-unit identification weight column b for one estimand:
#' \itemize{
#'   \item marginal: \eqn{b_i = wts_i}
#'   \item non-marginal: \eqn{b_i = wts_i \cdot 1(W_i = w) / P_\phi(W_i = w)}
#'     (Theorem 3.6(b))
#' }
#'
#' @param spec one specification from \code{\link{.po_specs}}.
#' @param weights1,weights0 per-regime identification weights.
#' @param W_Y 0/1 vector of individual-level treatments.
#' @param p1_own,p0_own per-unit own-treatment probabilities under phi_1 and
#'   phi_0.
#' @return A numeric weight vector.
#' @keywords internal
#' @examples
#' spec <- generalClusterExp:::.po_specs()$Ybar_w1_phi1   # Ybar(1, phi1)
#' W_Y <- c(1, 0, 1, 0)
#' ## = weights1 * 1(W = 1) / P_phi1(W = 1)
#' generalClusterExp:::.estimand_weight(spec,
#'   weights1 = rep(1, 4), weights0 = rep(1, 4), W_Y = W_Y,
#'   p1_own = rep(0.5, 4), p0_own = rep(0.2, 4))
.estimand_weight <- function(spec, weights1, weights0, W_Y, p1_own, p0_own) {
  wts <- if (spec$regime == 1) weights1 else weights0
  if (is.na(spec$w)) return(wts)
  p_own <- if (spec$regime == 1) p1_own else p0_own
  pw <- if (spec$w == 1L) p_own else (1 - p_own)
  ind <- as.numeric(W_Y == spec$w)
  b <- numeric(length(wts))
  pos <- pw > 0
  b[pos] <- wts[pos] * ind[pos] / pw[pos]
  b
}

#' Assemble the HAC covariance for a (method, design) pair
#'
#' Assembles the d x d HAC covariance for the supplied contribution matrix
#' \code{V}, following the variance rule for the (method, design) pair:
#' \itemize{
#'   \item \code{crn}: \eqn{\hat\Sigma_2 = V^\top K_3^+ V} (Theorem 5.9);
#'   \item \code{mrn} / \code{iptw}, Bernoulli:
#'     \eqn{\hat\Sigma_1 = \hat\Sigma^{K_1} \vee \hat\Sigma^{K_2}}, the Loewner
#'     max of the within- and overlap-kernel covariances (Theorem 5.4);
#'   \item \code{mrn} / \code{iptw}, complete:
#'     \eqn{\hat\Sigma_3 = \hat\Sigma_1 - \hat M} (Theorem 5.11(ii)).
#' }
#'
#' @param V n_unit x d contribution matrix.
#' @param A n_unit x n_unit adjacency matrix (with self-loops).
#' @param C integer vector of cluster ids.
#' @param W_C 0/1 vector of cluster-level treatments.
#' @param method \code{"mrn"}, \code{"iptw"} or \code{"crn"}.
#' @param design \code{"bernoulli"} or \code{"complete"}.
#' @return A list with the d x d matrix \code{Sigma} and a describing string
#'   \code{variance_type}.
#' @keywords internal
#' @examples
#' set.seed(1)
#' A <- Matrix::Matrix(outer(1:12, 1:12, function(i, j) abs(i - j) <= 1) * 1,
#'                     sparse = TRUE)
#' C <- rep(1:3, each = 4); W_C <- c(1, 0, 1)
#' V <- matrix(rnorm(24), 12, 2) / 12   # contributions of two estimands
#' generalClusterExp:::.assemble_Sigma(V, A, C, W_C,
#'                                     method = "mrn", design = "bernoulli")
#' generalClusterExp:::.assemble_Sigma(V, A, C, W_C,
#'                                     method = "crn", design = "bernoulli")
.assemble_Sigma <- function(V, A, C, W_C, method, design) {
  if (method == "crn") {
    K3 <- .kernel_K3_plus(A)
    return(list(Sigma = .sigma_from_kernel(V, K3),
                variance_type = "Sigma2: K3+ kernel (Theorem 5.9)"))
  }

  ker <- .cluster_kernels(A, C)
  S_within  <- .sigma_from_kernel(V, ker$K1)
  S_overlap <- .sigma_from_kernel(V, ker$K2)
  Sigma1 <- .lowner_max(S_within, S_overlap)

  if (design == "complete") {
    M <- .bias_correction_matrix(A, C, W_C, V)
    return(list(Sigma = Sigma1 - M,
                variance_type = "Sigma3: bias-corrected, complete randomization (Theorem 5.11(ii))"))
  }
  list(Sigma = Sigma1,
       variance_type = "Sigma1: Loewner-max, Bernoulli randomization (Theorem 5.4)")
}

#' Estimation engine: weights to tidy output
#'
#' Main engine shared by \code{\link{mrn}}, \code{\link{iptw}} and
#' \code{\link{crn}}: takes the per-regime identification weights, enumerates
#' the well-defined potential outcomes, forms the Hajek LW point estimates and
#' contribution columns, assembles the HAC covariance, and returns the tidy
#' output lists.
#'
#' @param A n_unit x n_unit adjacency matrix (with self-loops).
#' @param C integer vector of cluster ids.
#' @param Y numeric vector of observed outcomes.
#' @param W_C 0/1 vector of cluster-level treatments.
#' @param W_Y 0/1 vector of individual-level treatments.
#' @param weights1,weights0 per-regime identification weights (from one of the
#'   weight constructors).
#' @param p1_own,p0_own per-unit own-treatment probabilities under phi_1 and
#'   phi_0.
#' @param unit_weights optional per-unit weights (default: uniform).
#' @param method \code{"mrn"}, \code{"iptw"} or \code{"crn"}.
#' @param design \code{"bernoulli"} or \code{"complete"}.
#' @return See \code{\link{mrn}}.
#' @keywords internal
#' @examples
#' set.seed(1)
#' A <- Matrix::Matrix(outer(1:12, 1:12, function(i, j) abs(i - j) <= 1) * 1,
#'                     sparse = TRUE)
#' C <- rep(1:3, each = 4); W_C <- c(1, 0, 1)
#' W_Y <- rbinom(12, 1, ifelse(W_C[C] == 1, 0.5, 0.2))
#' Y <- rnorm(12, mean = 2 * as.vector(A %*% W_Y))
#' w <- generalClusterExp:::.iptw_weights(A, C, W_C,
#'                                        distC = "binom",
#'                                        paramsC = list(prob = 0.7))
#' out <- generalClusterExp:::.run_engine(A, C, Y, W_C, W_Y,
#'   weights1 = w$weights1, weights0 = w$weights0,
#'   p1_own = rep(0.5, 12), p0_own = rep(0.2, 12),
#'   unit_weights = NULL, method = "iptw", design = "bernoulli")
#' out$causal_effects
.run_engine <- function(A, C, Y, W_C, W_Y,
                        weights1, weights0,
                        p1_own, p0_own,
                        unit_weights, method, design) {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  N <- length(C)
  if (is.null(unit_weights)) unit_weights <- rep(1 / N, N)

  po_specs <- .po_specs()
  keys <- names(po_specs)

  # 1. well-definedness + point estimates + contribution columns
  well_defined <- vapply(po_specs, .is_well_defined, logical(1),
                         p1_own = p1_own, p0_own = p0_own)
  mu      <- stats::setNames(rep(NA_real_, length(keys)), keys)
  ok      <- stats::setNames(rep(FALSE, length(keys)), keys)
  Vcols   <- vector("list", length(keys)); names(Vcols) <- keys

  for (key in keys) {
    if (!well_defined[[key]]) next
    b <- .estimand_weight(po_specs[[key]], weights1, weights0, W_Y, p1_own, p0_own)
    est <- .lw_estimand(b, Y, unit_weights)
    mu[key]    <- est$mu
    ok[key]    <- est$ok
    Vcols[[key]] <- est$V
  }

  # 2. covariance over the estimable potential outcomes
  est_keys <- keys[ok]
  d <- length(est_keys)
  Sigma_full <- matrix(NA_real_, length(keys), length(keys),
                       dimnames = list(keys, keys))
  variance_type <- NA_character_
  if (d > 0) {
    V <- do.call(cbind, Vcols[est_keys])
    asm <- .assemble_Sigma(V, A, C, W_C, method, design)
    variance_type <- asm$variance_type
    Sigma_full[est_keys, est_keys] <- asm$Sigma
  }

  # 3. potential-outcome table
  po_var <- stats::setNames(rep(NA_real_, length(keys)), keys)
  for (key in est_keys) po_var[key] <- max(Sigma_full[key, key], 0)
  potential_outcomes <- data.frame(
    estimand     = vapply(po_specs, function(s) s$label, character(1)),
    key          = keys,
    estimate     = mu,
    variance     = po_var,
    se           = sqrt(po_var),
    well_defined = ok,
    row.names    = NULL,
    stringsAsFactors = FALSE
  )

  # 4. causal-effect table (contrasts of the estimands)
  eff_specs <- .effect_specs()
  ce_rows <- lapply(names(eff_specs), function(name) {
    sp <- eff_specs[[name]]
    avail <- ok[[sp$plus]] && ok[[sp$minus]]
    estv <- varv <- NA_real_
    if (avail) {
      estv <- mu[[sp$plus]] - mu[[sp$minus]]
      cc <- stats::setNames(numeric(length(keys)), keys)
      cc[sp$plus] <- 1; cc[sp$minus] <- -1
      cc <- cc[est_keys]
      varv <- max(as.numeric(t(cc) %*% Sigma_full[est_keys, est_keys] %*% cc), 0)
    }
    data.frame(effect = sp$label, key = name, estimate = estv,
               variance = varv, se = sqrt(varv), well_defined = avail,
               row.names = NULL, stringsAsFactors = FALSE)
  })
  causal_effects <- do.call(rbind, ce_rows)

  list(
    method             = method,
    design             = design,
    variance_type      = variance_type,
    causal_effects     = causal_effects,
    potential_outcomes = potential_outcomes,
    Sigma              = Sigma_full,
    weights1           = weights1,
    weights0           = weights0
  )
}
