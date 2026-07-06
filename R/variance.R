##############################################################
# Estimand assembly and the variance engine.
#
# Given the per-regime identification weights (weights1, weights0) produced by
# one of the weight constructors, this module
#   1. enumerates the well-defined population average potential outcomes
#      (marginal and non-marginal), Section 4-5 of main.tex;
#   2. forms each Hajek LW point estimate and its per-unit contribution column
#      \hat V_ij = u_i * b_ij * (Y_ij - \hat mu);
#   3. assembles the d x d HAC covariance \hat Sigma over all estimands using
#      the kernel rule appropriate to the (method, design) pair;
#   4. reads off variances of every potential outcome and of the direct /
#      indirect / total / overall causal effects (contrasts of the estimands).
##############################################################

# Hajek LW estimate mu and per-unit centred contribution V for one estimand
# with per-unit identification weight b and unit weight u.
#
# NOTE: following the definition of \hat Sigma^K in Section 6 of main.tex, the
# contribution is V_i = u_i b_i (Y_i - mu) WITHOUT dividing by the Hajek
# denominator D = sum_i u_i b_i (which converges to 1; see the estimator
# definition in Section 7). The point estimate mu itself is the properly
# normalised Hajek ratio. TODO: if finite-sample D departs noticeably from 1,
# divide V by D for a sharper variance.
.lw_estimand <- function(b, Y, u) {
  den <- sum(u * b)
  if (!is.finite(den) || den == 0) {
    return(list(mu = NA_real_, V = rep(0, length(Y)), ok = FALSE))
  }
  mu <- sum(u * b * Y) / den
  list(mu = mu, V = u * b * (Y - mu), ok = TRUE)
}

# Canonical list of the six population average potential outcomes.
# regime: 1 -> phi_1 (weights1), 0 -> phi_0 (weights0).
# w: NA -> marginal Ybar(phi); 0/1 -> non-marginal Ybar(w, phi).
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

# Causal-effect contrasts over the potential-outcome keys (Section 4-5):
#   overall   CE^O          = Ybar(phi1)    - Ybar(phi0)
#   direct    CE^D(phi_c)   = Ybar(1,phi_c) - Ybar(0,phi_c)
#   indirect  CE^I(w)       = Ybar(w,phi1)  - Ybar(w,phi0)
#   total     CE^T          = Ybar(1,phi1)  - Ybar(0,phi0)
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

# Decide which potential outcomes are well-defined from the regime own-treatment
# probabilities. Marginal outcomes are always well-defined; Ybar(w, phi) needs
# P_phi(W = w) > 0 for every contributing unit.
.is_well_defined <- function(spec, p1_own, p0_own) {
  if (is.na(spec$w)) return(TRUE)
  p_own <- if (spec$regime == 1) p1_own else p0_own
  pw <- if (spec$w == 1L) p_own else (1 - p_own)
  all(pw > 0)
}

# Build the per-unit identification weight column b for one estimand:
#   marginal:      b_i = wts_i
#   non-marginal:  b_i = wts_i * 1(W_i = w) / P_phi(W_i = w)     (Theorem all_weights (b))
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

# Assemble the d x d HAC covariance for the supplied contribution matrix V,
# following the variance rule for (method, design):
#   - crn          : Sigma_2 = V^T K3+ V                     (Theorem var_2, any design)
#   - mrn / iptw,
#       bernoulli  : Sigma_1 = Sigma^{K_within} \vee Sigma^{K_overlap}   (Theorem var_1)
#       complete   : Sigma_1 - M_hat                          (bias-corrected, Section 6.3)
.assemble_Sigma <- function(V, A, C, W_C, method, design) {
  if (method == "crn") {
    K3 <- .kernel_K3_plus(A)
    return(list(Sigma = .sigma_from_kernel(V, K3), variance_type = "var_2 (K3+)"))
  }

  ker <- .cluster_kernels(A, C)
  S_within  <- .sigma_from_kernel(V, ker$K1)
  S_overlap <- .sigma_from_kernel(V, ker$K2)
  Sigma1 <- .lowner_max(S_within, S_overlap)

  if (design == "complete") {
    M <- .bias_correction_matrix(A, C, W_C, V)
    return(list(Sigma = Sigma1 - M, variance_type = "var_1 bias-corrected (complete rand.)"))
  }
  list(Sigma = Sigma1, variance_type = "var_1 (Loewner-max, Bernoulli)")
}

# Main engine: weights -> estimands -> Sigma -> tidy output.
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
