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
# identification matrix B_marg in Section 7 of the companion manuscript.
##############################################################

#' Within-cluster / cluster-level assignment probability
#'
#' Probability of \code{treated_neighbors} treated out of
#' \code{total_neighbors_in_cluster} draws, under the within-cluster (or
#' cluster-level) assignment law \code{dist} with parameters \code{params}.
#'
#' @param dist \code{"binom"} or \code{"hypergeom"}.
#' @param params for \code{"binom"}, \code{list(prob = p)}; for
#'   \code{"hypergeom"}, \code{list(m = m)} where \code{m} may be a scalar or a
#'   per-cluster vector.
#' @param treated_neighbors number of treated draws (vectorised).
#' @param total_neighbors_in_cluster number of draws.
#' @param total_size population size of the cluster (or the number of clusters
#'   for a cluster-level law).
#' @param cluster_id optional cluster id used to pick the entry of a
#'   per-cluster \code{params$m}.
#' @return A numeric vector of probabilities.
#' @keywords internal
#' @examples
#' ## Binomial law: P(2 of 3 neighbours treated) with p = 0.5
#' generalClusterExp:::.calc_treatment_probability(
#'   "binom", list(prob = 0.5), 2, 3, total_size = 10)
#'
#' ## Hypergeometric law: 2 treated among 3 draws from a cluster of
#' ## size 10 containing m = 5 treated units
#' generalClusterExp:::.calc_treatment_probability(
#'   "hypergeom", list(m = 5), 2, 3, total_size = 10)
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

#' Marginal Radon--Nikodym (MRN) identification weights
#'
#' Per-unit marginal Radon--Nikodym weights \eqn{\alpha_{ij}(\phi_c)} (Theorem
#' all_weights / Prop. weights (i)) for the two counterfactual regimes. The
#' weights marginalise over the within-cluster assignment of each neighbouring
#' cluster, so a unit contributes whenever its cluster neighbourhood has
#' positive probability under both regimes.
#'
#' @param A n_unit x n_unit adjacency matrix (with self-loops).
#' @param C integer vector of cluster ids.
#' @param W_C 0/1 vector of cluster-level treatments, indexed by cluster id.
#' @param W_Y 0/1 vector of individual-level treatments.
#' @param dist1,params1 within-cluster assignment law of the treated regime phi_1.
#' @param dist0,params0 within-cluster assignment law of the control regime phi_0.
#' @param distC,paramsC cluster-level assignment law.
#' @param unit_weights optional per-unit weights (default: uniform); only units
#'   with nonzero weight are processed.
#' @return A list with numeric vectors \code{weights1} and \code{weights0}.
#' @keywords internal
#' @examples
#' set.seed(1)
#' A <- Matrix::Matrix(outer(1:12, 1:12, function(i, j) abs(i - j) <= 1) * 1,
#'                     sparse = TRUE)
#' C <- rep(1:3, each = 4)
#' W_C <- c(1, 0, 1)
#' W_Y <- rbinom(12, 1, ifelse(W_C[C] == 1, 0.5, 0.2))
#' w <- generalClusterExp:::.marg_rn_weights(A, C, W_C, W_Y,
#'   dist1 = "binom", params1 = list(prob = 0.5),
#'   dist0 = "binom", params0 = list(prob = 0.2),
#'   distC = "binom", paramsC = list(prob = 0.7))
#' cbind(phi1 = w$weights1, phi0 = w$weights0)
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

#' Complete / joint Radon--Nikodym (CRN) identification weights
#'
#' Per-unit complete (joint) Radon--Nikodym weights \eqn{\alpha_{ij}^{comp}}
#' (Prop. weights (ii)); cluster-agnostic. Unlike the marginal weights, these
#' condition on the realised cluster-level assignment, taking the product of
#' the per-cluster likelihood ratios over the neighbouring clusters.
#'
#' @inheritParams .marg_rn_weights
#' @return A list with numeric vectors \code{weights1} and \code{weights0}.
#' @keywords internal
#' @examples
#' set.seed(1)
#' A <- Matrix::Matrix(outer(1:12, 1:12, function(i, j) abs(i - j) <= 1) * 1,
#'                     sparse = TRUE)
#' C <- rep(1:3, each = 4)
#' W_C <- c(1, 0, 1)
#' W_Y <- rbinom(12, 1, ifelse(W_C[C] == 1, 0.5, 0.2))
#' w <- generalClusterExp:::.complete_rn_weights(A, C, W_C, W_Y,
#'   dist1 = "binom", params1 = list(prob = 0.5),
#'   dist0 = "binom", params0 = list(prob = 0.2))
#' cbind(phi1 = w$weights1, phi0 = w$weights0)
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

#' IPTW identification weights
#'
#' Cluster-level inverse-probability weights (Leung 2025; Prop. weights (iv)):
#' \deqn{\beta_{ij}(\phi_1) = 1\{C_{N_{ij}} = 1\} / P(C_{N_{ij}} = 1), \quad
#'       \beta_{ij}(\phi_0) = 1\{C_{N_{ij}} = 0\} / P(C_{N_{ij}} = 0),}
#' where the probabilities use the cluster-level assignment law
#' (\code{distC}, \code{paramsC}) over the \eqn{k = |N^{cl}_{ij}|} distinct
#' neighbouring clusters. Units whose neighbouring clusters are not all
#' treated (for phi_1) or all control (for phi_0) get weight 0.
#'
#' @inheritParams .marg_rn_weights
#' @return A list with numeric vectors \code{weights1} and \code{weights0}.
#' @keywords internal
#' @examples
#' A <- Matrix::Matrix(outer(1:12, 1:12, function(i, j) abs(i - j) <= 1) * 1,
#'                     sparse = TRUE)
#' C <- rep(1:3, each = 4)
#' W_C <- c(1, 0, 1)
#' w <- generalClusterExp:::.iptw_weights(A, C, W_C,
#'                                        distC = "binom",
#'                                        paramsC = list(prob = 0.7))
#' ## nonzero only where the cluster neighbourhood is all-treated / all-control
#' cbind(phi1 = w$weights1, phi0 = w$weights0)
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

#' Per-unit own-treatment probability under a within-cluster regime
#'
#' Marginal probability \eqn{P_{\phi}(W_{ij} = 1)} of a unit's own treatment
#' under the within-cluster regime \eqn{\phi} = (\code{dist}, \code{params}).
#' Used both to form the conditional (non-marginal) weights and to decide
#' which non-marginal estimands are well-defined.
#'
#' @param dist \code{"binom"} or \code{"hypergeom"}.
#' @param params law parameters; see \code{\link{.calc_treatment_probability}}.
#' @param C integer vector of cluster ids.
#' @return A numeric vector of length \code{length(C)}.
#' @keywords internal
#' @examples
#' C <- rep(1:3, each = 4)
#' generalClusterExp:::.own_treatment_prob1("binom", list(prob = 0.5), C)
#' ## hypergeometric: m = 2 treated in each cluster of size 4 -> 1/2 each
#' generalClusterExp:::.own_treatment_prob1("hypergeom", list(m = 2), C)
.own_treatment_prob1 <- function(dist, params, C) {
  N <- length(C)
  if (dist == "binom") {
    rep(params$prob, N)
  } else if (dist == "hypergeom") {
    # Under a within-cluster hypergeometric design with m treated out of
    # cluster_size, the own-treatment marginal probability is m / cluster_size.
    sizes <- tabulate(C)
    cs <- sizes[C]
    m <- params$m
    if (length(m) == 1) m / cs else m[C] / cs
  } else {
    stop("Unsupported within-cluster regime distribution: ", dist)
  }
}
