##############################################################
# Top-level estimators: mrn(), iptw(), crn().
#
# All three share the same interface and the same rich output (overall, direct,
# indirect and total causal effects; all well-defined marginal and non-marginal
# population average potential outcomes; and a variance estimate for every
# quantity). They differ only in (a) the identification weights and (b) the HAC
# variance rule, matching the MRN / IPTW / CRN estimators of Section 4.
##############################################################

#' Infer the cluster-level assignment design
#'
#' Infers the cluster-level assignment design from the cluster-level law:
#' \code{"binom"} maps to cluster-level independent Bernoulli randomization,
#' \code{"hypergeom"} to cluster-level (stratified) complete randomization.
#' An explicit \code{design} argument overrides the inference.
#'
#' @param distC cluster-level assignment law, \code{"binom"} or
#'   \code{"hypergeom"}.
#' @param design optional override, \code{"bernoulli"} or \code{"complete"};
#'   \code{NULL} to infer from \code{distC}.
#' @return \code{"bernoulli"} or \code{"complete"}.
#' @keywords internal
#' @examples
#' generalClusterExp:::.infer_design("binom", NULL)       # "bernoulli"
#' generalClusterExp:::.infer_design("hypergeom", NULL)   # "complete"
#' generalClusterExp:::.infer_design("binom", "complete") # explicit override
.infer_design <- function(distC, design) {
  if (!is.null(design)) {
    design <- match.arg(design, c("bernoulli", "complete"))
    return(design)
  }
  if (distC == "binom") "bernoulli"
  else if (distC == "hypergeom") "complete"
  else stop("Cannot infer design from distC = '", distC,
            "'. Pass design = 'bernoulli' or 'complete' explicitly.")
}

#' MRN estimator (marginal Radon--Nikodym weights)
#'
#' Linear-weighted Hajek estimator using the marginal Radon--Nikodym
#' identification weights (Theorem 3.6, Proposition 3.9(i)). Returns all
#' well-defined causal effects and population average potential outcomes with
#' variance estimates. The variance estimator is the Loewner-max HAC estimator
#' \eqn{\hat\Sigma_1} (Theorem 5.4) under cluster-level Bernoulli
#' randomization, and the bias-corrected estimator
#' \eqn{\hat\Sigma_3 = \hat\Sigma_1 - \hat M} (Theorem 5.11(ii)) under
#' cluster-level complete randomization.
#'
#' @param A n_unit x n_unit adjacency matrix of the interference network
#'   (with self-loops, as used throughout the project).
#' @param C integer vector of length n_unit giving each unit's cluster id.
#' @param Y numeric vector of observed outcomes (length n_unit).
#' @param W_C 0/1 vector of cluster-level treatments, indexed by cluster id.
#' @param W_Y 0/1 vector of individual-level treatments (length n_unit).
#' @param dist1,params1 within-cluster assignment law of the treated regime phi_1.
#' @param dist0,params0 within-cluster assignment law of the control regime phi_0.
#' @param distC,paramsC cluster-level assignment law.
#' @param unit_weights optional per-unit weights u_i = g_i / N_i
#'   (default: uniform 1 / n_unit).
#' @param design optional override, "bernoulli" or "complete"; inferred from
#'   \code{distC} when NULL.
#' @return A list with \code{causal_effects}, \code{potential_outcomes} (data
#'   frames of estimates / variances / standard errors), the estimand covariance
#'   matrix \code{Sigma}, and the chosen \code{variance_type}.
#' @examples
#' ## Small geometric cluster network (a scaled-down version of the
#' ## simulation pipeline in data_generation.R): units live around ten
#' ## cluster centres and interfere with all units within distance 1.5.
#' set.seed(1)
#' n_per <- 6; K <- 10; n <- n_per * K
#' C <- rep(seq_len(K), each = n_per)                      # cluster ids
#' centers <- cbind(rep(1:5, 2), rep(1:2, each = 5)) * 4   # cluster centres
#' loc <- centers[C, ] + matrix(runif(2 * n, -1, 1), n, 2) # unit locations
#' A <- Matrix::Matrix(as.matrix(dist(loc)) <= 1.5, sparse = TRUE) * 1
#'
#' ## Two-stage assignment: clusters are treated w.p. 0.7 (cluster-level
#' ## Bernoulli law distC); units in treated clusters get W = 1 w.p. 0.5
#' ## (regime phi_1), units in control clusters w.p. 0.2 (regime phi_0).
#' W_C <- rbinom(K, 1, 0.7)
#' W_Y <- rbinom(n, 1, ifelse(W_C[C] == 1, 0.5, 0.2))
#'
#' ## Outcomes with neighbourhood interference (spillover beta = 2 from
#' ## treated neighbours plus interaction gamma = 1 for treated units).
#' Y <- as.vector(-1 + A %*% (2 * W_Y) +
#'                  W_Y * as.vector(A %*% (1 * W_Y)) + rnorm(n))
#'
#' ## MRN estimator under cluster-level Bernoulli randomization:
#' ## variance is the Loewner-max HAC estimator Sigma1 (Theorem 5.4).
#' fit <- mrn(A, C, Y, W_C, W_Y,
#'            dist1 = "binom", params1 = list(prob = 0.5),
#'            dist0 = "binom", params0 = list(prob = 0.2),
#'            distC = "binom", paramsC = list(prob = 0.7))
#' fit$variance_type
#' fit$causal_effects
#' fit$potential_outcomes
#'
#' ## Cluster-level complete randomization (exactly 7 of 10 clusters
#' ## treated, distC = "hypergeom"): the bias-corrected variance
#' ## estimator Sigma3 of Theorem 5.11(ii) is used automatically.
#' W_C2 <- as.integer(seq_len(K) %in% sample.int(K, 7))
#' W_Y2 <- rbinom(n, 1, ifelse(W_C2[C] == 1, 0.5, 0.2))
#' Y2 <- as.vector(-1 + A %*% (2 * W_Y2) +
#'                   W_Y2 * as.vector(A %*% (1 * W_Y2)) + rnorm(n))
#' fit_bc <- mrn(A, C, Y2, W_C2, W_Y2,
#'               params1 = list(prob = 0.5), params0 = list(prob = 0.2),
#'               distC = "hypergeom", paramsC = list(m = 7))
#' fit_bc$variance_type
#' fit_bc$causal_effects
#' @export
mrn <- function(A, C, Y, W_C, W_Y,
                dist1 = "binom", params1 = list(prob = 0.5),
                dist0 = "binom", params0 = list(prob = 0),
                distC = "binom", paramsC = list(prob = 0.7),
                unit_weights = NULL, design = NULL) {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  des <- .infer_design(distC, design)
  w <- .marg_rn_weights(A, C, W_C, W_Y,
                        dist1 = dist1, params1 = params1,
                        dist0 = dist0, params0 = params0,
                        distC = distC, paramsC = paramsC,
                        unit_weights = unit_weights)
  .run_engine(A, C, Y, W_C, W_Y,
              weights1 = w$weights1, weights0 = w$weights0,
              p1_own = .own_treatment_prob1(dist1, params1, C),
              p0_own = .own_treatment_prob1(dist0, params0, C),
              unit_weights = unit_weights, method = "mrn", design = des)
}

#' IPTW estimator (inverse-probability-of-treatment weights)
#'
#' Linear-weighted Hajek estimator using the cluster-level IPTW identification
#' weights of Leung (2025) (Proposition 3.9(iv)). Same outputs and variance
#' rule as \code{\link{mrn}}: \eqn{\hat\Sigma_1} (Theorem 5.4) under Bernoulli
#' randomization and the bias-corrected \eqn{\hat\Sigma_3} (Theorem 5.11(ii))
#' under complete randomization.
#'
#' @inheritParams mrn
#' @return See \code{\link{mrn}}.
#' @examples
#' ## Same data-generating pipeline as in the mrn() example: geometric
#' ## cluster network, two-stage Bernoulli treatment, interference outcomes.
#' set.seed(1)
#' n_per <- 6; K <- 10; n <- n_per * K
#' C <- rep(seq_len(K), each = n_per)
#' centers <- cbind(rep(1:5, 2), rep(1:2, each = 5)) * 4
#' loc <- centers[C, ] + matrix(runif(2 * n, -1, 1), n, 2)
#' A <- Matrix::Matrix(as.matrix(dist(loc)) <= 1.5, sparse = TRUE) * 1
#'
#' W_C <- rbinom(K, 1, 0.7)
#' W_Y <- rbinom(n, 1, ifelse(W_C[C] == 1, 0.5, 0.2))
#' Y <- as.vector(-1 + A %*% (2 * W_Y) +
#'                  W_Y * as.vector(A %*% (1 * W_Y)) + rnorm(n))
#'
#' ## IPTW keeps only units whose neighbouring clusters are all treated
#' ## (for phi_1) or all control (for phi_0), reweighting by the inverse
#' ## exposure probability.
#' fit <- iptw(A, C, Y, W_C, W_Y,
#'             params1 = list(prob = 0.5), params0 = list(prob = 0.2),
#'             distC = "binom", paramsC = list(prob = 0.7))
#' fit$variance_type
#' fit$causal_effects
#' @export
iptw <- function(A, C, Y, W_C, W_Y,
                 dist1 = "binom", params1 = list(prob = 0.5),
                 dist0 = "binom", params0 = list(prob = 0),
                 distC = "binom", paramsC = list(prob = 0.7),
                 unit_weights = NULL, design = NULL) {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  des <- .infer_design(distC, design)
  w <- .iptw_weights(A, C, W_C, distC = distC, paramsC = paramsC)
  .run_engine(A, C, Y, W_C, W_Y,
              weights1 = w$weights1, weights0 = w$weights0,
              p1_own = .own_treatment_prob1(dist1, params1, C),
              p0_own = .own_treatment_prob1(dist0, params0, C),
              unit_weights = unit_weights, method = "iptw", design = des)
}

#' CRN estimator (complete Radon--Nikodym weights)
#'
#' Cluster-agnostic linear-weighted Hajek estimator using the complete (joint)
#' Radon--Nikodym identification weights (Theorem 3.7, Proposition 3.9(ii)).
#' Same outputs as \code{\link{mrn}}. The variance estimator is the
#' cluster-agnostic HAC estimator \eqn{\hat\Sigma_2} built from the PSD part
#' \eqn{K_3^+} of the interference kernel, and is used for every value of
#' \code{design}. Its conservativeness is established in Theorem 5.9 under
#' independent Bernoulli randomization; the companion manuscript does not
#' cover \eqn{\hat\Sigma_2} under complete randomization, so treat the
#' variance as heuristic in that case.
#'
#' @inheritParams mrn
#' @return See \code{\link{mrn}}.
#' @examples
#' ## Same data-generating pipeline as in the mrn() example: geometric
#' ## cluster network, two-stage Bernoulli treatment, interference outcomes.
#' set.seed(1)
#' n_per <- 6; K <- 10; n <- n_per * K
#' C <- rep(seq_len(K), each = n_per)
#' centers <- cbind(rep(1:5, 2), rep(1:2, each = 5)) * 4
#' loc <- centers[C, ] + matrix(runif(2 * n, -1, 1), n, 2)
#' A <- Matrix::Matrix(as.matrix(dist(loc)) <= 1.5, sparse = TRUE) * 1
#'
#' W_C <- rbinom(K, 1, 0.7)
#' W_Y <- rbinom(n, 1, ifelse(W_C[C] == 1, 0.5, 0.2))
#' Y <- as.vector(-1 + A %*% (2 * W_Y) +
#'                  W_Y * as.vector(A %*% (1 * W_Y)) + rnorm(n))
#'
#' ## CRN is cluster-agnostic: its variance uses the PSD-projected
#' ## interference kernel K3+ (Sigma2, Theorem 5.9).
#' fit <- crn(A, C, Y, W_C, W_Y,
#'            params1 = list(prob = 0.5), params0 = list(prob = 0.2),
#'            distC = "binom", paramsC = list(prob = 0.7))
#' fit$variance_type
#' fit$causal_effects
#' @export
crn <- function(A, C, Y, W_C, W_Y,
                dist1 = "binom", params1 = list(prob = 0.5),
                dist0 = "binom", params0 = list(prob = 0),
                distC = "binom", paramsC = list(prob = 0.7),
                unit_weights = NULL, design = NULL) {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  des <- .infer_design(distC, design)
  w <- .complete_rn_weights(A, C, W_C, W_Y,
                            dist1 = dist1, params1 = params1,
                            dist0 = dist0, params0 = params0,
                            unit_weights = unit_weights)
  .run_engine(A, C, Y, W_C, W_Y,
              weights1 = w$weights1, weights0 = w$weights0,
              p1_own = .own_treatment_prob1(dist1, params1, C),
              p0_own = .own_treatment_prob1(dist0, params0, C),
              unit_weights = unit_weights, method = "crn", design = des)
}
