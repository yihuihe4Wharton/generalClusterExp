##############################################################
# Identification weights for the three estimators.
#
# Each constructor returns a list(weights1, weights0) of per-unit
# identification weights for the two counterfactual regimes phi_1 and phi_0:
#   - mrn : marginal Radon-Nikodym weights alpha_ij(phi_c)          (Prop. weights (i))
#   - crn : complete (joint) Radon-Nikodym weights alpha_ij^comp    (Prop. weights (ii))
#   - iptw: cluster-level inverse-probability weights               (Prop. weights (iv))
#
# weights1[i] estimates the contribution of unit i to Ybar(phi_1);
# weights0[i] the contribution to Ybar(phi_0). These are the columns of the
# identification matrix B_marg in Section 7 of main.tex.
##############################################################

# Probability of `treated_neighbors` treated out of `total_neighbors_in_cluster`
# draws, under the within-cluster (or cluster-level) assignment law.
.calc_treatment_probability <- function(dist, params,
                                        treated_neighbors,
                                        total_neighbors_in_cluster,
                                        total_size,
                                        cluster_id = NULL) {
  if (dist == "binom") {
    stats::dbinom(treated_neighbors,
                  size = total_neighbors_in_cluster,
                  prob = params$prob)
  } else if (dist == "hypergeom") {
    if (is.null(cluster_id) || length(params$m) == 1) {
      stats::dhyper(treated_neighbors,
                    m = params$m,
                    n = total_size - params$m,
                    k = total_neighbors_in_cluster)
    } else {
      stats::dhyper(treated_neighbors,
                    m = params$m[cluster_id],
                    n = total_size - params$m[cluster_id],
                    k = total_neighbors_in_cluster)
    }
  } else {
    stop("Unsupported distribution: ", dist)
  }
}

# Marginal Radon-Nikodym (MRN) weights (Theorem all_weights / Prop. weights (i)).
# Ported verbatim from our_helper::marg_Radon_weights; only the auxiliary
# bookkeeping needed for the leave-cluster-out jackknife is dropped.
.marg_rn_weights <- function(A, C, W_C, W_Y,
                             dist1 = "binom", params1 = list(prob = 0.5),
                             dist0 = "binom", params0 = list(prob = 0),
                             distC = "binom", paramsC = list(prob = 0.7),
                             unit_weights = NULL) {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  n  <- nrow(A)
  cK <- length(unique(C))
  if (is.null(unit_weights)) unit_weights <- rep(1 / n, n)

  weights1 <- rep(0, n)
  weights0 <- rep(0, n)

  adj <- lapply(seq_len(n), function(i) which(A[, i] != 0))
  active_individuals <- which(unit_weights != 0)
  N_total <- max(C)

  all_m_products_dp <- function(x, M = length(x)) {
    x <- as.numeric(x); K <- length(x); M <- as.integer(min(M, K))
    dp <- numeric(M + 1); dp[1] <- 1
    for (xi in x) {
      up_to <- min(M, K)
      for (k in seq(up_to, 1L, by = -1L)) dp[k + 1L] <- dp[k + 1L] + xi * dp[k]
    }
    dp
  }

  for (i in active_individuals) {
    nb <- adj[[i]]
    if (!length(nb)) next
    nb_clusters <- C[nb]
    totals <- tabulate(nb_clusters, nbins = cK)
    treated_rs <- rowsum(W_Y[nb], nb_clusters, reorder = FALSE)
    treated <- numeric(cK)
    treated[as.integer(rownames(treated_rs))] <- treated_rs[, 1]

    js <- which(totals > 0)
    weights1_factor <- rep(0, length(js))
    weights0_factor <- rep(0, length(js))
    for (j in seq_along(js)) {
      cid <- js[j]
      total_neighbors_in_cluster <- totals[cid]
      treated_neighbors <- treated[cid]
      prob1 <- .calc_treatment_probability(dist1, params1, treated_neighbors,
                                           total_neighbors_in_cluster, sum(C == cid), cid)
      prob0 <- .calc_treatment_probability(dist0, params0, treated_neighbors,
                                           total_neighbors_in_cluster, sum(C == cid), cid)
      weights0_factor[j] <- prob1 / prob0
      weights1_factor[j] <- prob0 / prob1
    }

    if (!all(is.finite(weights0_factor))) {
      weights0[i] <- 0
    } else {
      e0 <- all_m_products_dp(weights0_factor)
      terms0 <- vapply(seq_len(length(js) + 1), function(j) {
        e0[j] * .calc_treatment_probability(distC, paramsC, j - 1, length(js), N_total) /
          choose(length(js), j - 1)
      }, numeric(1))
      weights0[i] <- 1 / sum(terms0)
    }

    if (!all(is.finite(weights1_factor))) {
      weights1[i] <- 0
    } else {
      e1 <- all_m_products_dp(weights1_factor)
      terms1 <- vapply(seq_len(length(js) + 1), function(j) {
        e1[j] * .calc_treatment_probability(distC, paramsC, length(js) + 1 - j, length(js), N_total) /
          choose(length(js), length(js) + 1 - j)
      }, numeric(1))
      weights1[i] <- 1 / sum(terms1)
    }
  }

  list(weights1 = weights1, weights0 = weights0)
}

# Complete / joint Radon-Nikodym (CRN) weights (Prop. weights (ii)).
# Ported from our_helper::Radon_weights (cluster-agnostic weights).
.complete_rn_weights <- function(A, C, W_C, W_Y,
                                 dist1 = "binom", params1 = list(prob = 0.5),
                                 dist0 = "binom", params0 = list(prob = 0),
                                 unit_weights = NULL) {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  n  <- nrow(A)
  cK <- length(unique(C))
  if (is.null(unit_weights)) unit_weights <- rep(1 / n, n)

  weights1 <- rep(0, n)
  weights0 <- rep(0, n)

  adj <- lapply(seq_len(n), function(i) which(A[, i] != 0))
  active_individuals <- which(unit_weights != 0)

  for (i in active_individuals) {
    weights1[i] <- 1
    weights0[i] <- 1
    nb <- adj[[i]]
    if (!length(nb)) next
    nb_clusters <- C[nb]
    totals <- tabulate(nb_clusters, nbins = cK)
    treated_rs <- rowsum(W_Y[nb], nb_clusters, reorder = FALSE)
    treated <- numeric(cK)
    treated[as.integer(rownames(treated_rs))] <- treated_rs[, 1]

    for (j in which(totals > 0)) {
      total_neighbors_in_cluster <- totals[j]
      treated_neighbors <- treated[j]
      prob1 <- .calc_treatment_probability(dist1, params1, treated_neighbors,
                                           total_neighbors_in_cluster, sum(C == j), j)
      prob0 <- .calc_treatment_probability(dist0, params0, treated_neighbors,
                                           total_neighbors_in_cluster, sum(C == j), j)
      if (W_C[j] == 1) {
        weights0[i] <- weights0[i] * prob0 / prob1
      } else {
        weights1[i] <- weights1[i] * prob1 / prob0
      }
    }
  }

  list(weights1 = weights1, weights0 = weights0)
}

# IPTW weights (Leung 2025; Prop. weights (iv)):
#   beta_ij(phi_1) = 1{C_{N_ij} = 1} / P(C_{N_ij} = 1),
#   beta_ij(phi_0) = 1{C_{N_ij} = 0} / P(C_{N_ij} = 0),
# where the probabilities use the cluster-level assignment law (distC, paramsC)
# over the k = |N^{cl}_ij| distinct neighbouring clusters.
.iptw_weights <- function(A, C, W_C, distC = "binom", paramsC = list(prob = 0.7)) {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  n <- nrow(A)

  cluster_sets <- lapply(seq_len(n), function(i) unique(C[which(A[, i] != 0)]))
  k <- vapply(cluster_sets, length, integer(1))
  total_cl <- max(C)

  pi1 <- .calc_treatment_probability(distC, paramsC, k, k, total_cl)  # all neighbours treated
  pi0 <- .calc_treatment_probability(distC, paramsC, 0, k, total_cl)  # all neighbours control

  Tind <- vapply(cluster_sets, function(cl) {
    wc <- W_C[cl]
    if (all(wc == 1)) 1 else if (all(wc == 0)) 0 else -99
  }, numeric(1))

  weights1 <- ifelse(Tind == 1, 1 / pi1, 0)
  weights0 <- ifelse(Tind == 0, 1 / pi0, 0)
  weights1[!is.finite(weights1)] <- 0
  weights0[!is.finite(weights0)] <- 0

  list(weights1 = weights1, weights0 = weights0)
}

# Per-unit marginal probability of own treatment P_{phi}(W_ij = 1) under the
# within-cluster regime phi = (dist, params). Used both to form the conditional
# (non-marginal) weights and to decide which non-marginal estimands are
# well-defined.
.own_treatment_prob1 <- function(dist, params, C) {
  N <- length(C)
  if (dist == "binom") {
    rep(params$prob, N)
  } else if (dist == "hypergeom") {
    # TODO: confirm against Assumption (treatment mechanism), Section 4 of
    # main.tex. For a within-cluster hypergeometric law with m treated out of
    # the cluster size, the own-treatment marginal is m / cluster_size.
    sizes <- tabulate(C)
    cs <- sizes[C]
    m <- params$m
    if (length(m) == 1) m / cs else m[C] / cs
  } else {
    stop("Unsupported within-cluster regime distribution: ", dist)
  }
}
