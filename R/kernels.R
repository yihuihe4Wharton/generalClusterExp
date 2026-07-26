##############################################################
# HAC variance kernels and matrix operations (Section 6 of the companion
# manuscript).
#
# A kernel K is an n_unit x n_unit symmetric 0/1 (or PSD-projected) matrix that
# selects which pairs of units are treated as dependent. Given the n_unit x d
# matrix V of per-unit, per-estimand contributions \hat V_ij (one column per
# estimand), the corresponding HAC covariance estimator is
#       Sigma^K = V^T K V        (d x d).
##############################################################

# ---- PSD projection helpers -------------------------------------------------

#' PSD projection via dense eigendecomposition
#'
#' Projects a symmetric kernel matrix onto the positive semi-definite cone by
#' flooring its eigenvalues at \code{tol} (dense, exact PSD part).
#'
#' @param K symmetric kernel matrix (dense or sparse).
#' @param tol floor applied to the eigenvalues.
#' @return A dense symmetric PSD matrix.
#' @keywords internal
#' @examples
#' K <- matrix(c(1, 2, 2, 1), 2, 2)    # indefinite: eigenvalues 3 and -1
#' Kp <- generalClusterExp:::.psd_via_eigen(K)
#' eigen(Kp, symmetric = TRUE)$values  # all >= 0
.psd_via_eigen <- function(K, tol = 1e-12) {
  Kd <- as.matrix(K); Kd <- (Kd + t(Kd)) / 2
  eg <- eigen(Kd, symmetric = TRUE)
  lam_pos <- pmax(eg$values, tol)
  Kpsd <- eg$vectors %*% (lam_pos * t(eg$vectors))
  (Kpsd + t(Kpsd)) / 2
}

#' PSD correction via diagonal shift
#'
#' Makes a symmetric sparse kernel PSD by adding \code{(-lambda_min + eps)} to
#' the diagonal when its smallest eigenvalue \code{lambda_min} is negative.
#' Keeps the matrix sparse; uses \pkg{RSpectra} for the extreme eigenvalue when
#' available.
#'
#' @param K symmetric kernel matrix (coerced to sparse \code{dgCMatrix}).
#' @param eps extra diagonal mass added on top of \code{-lambda_min}.
#' @return A sparse symmetric PSD matrix.
#' @keywords internal
#' @examples
#' ## chain adjacency: indefinite (smallest eigenvalue is negative)
#' K <- Matrix::Matrix(abs(outer(1:4, 1:4, "-")) == 1, sparse = TRUE) * 1
#' Ks <- generalClusterExp:::.psd_via_shift(K)
#' min(eigen(as.matrix(Ks), symmetric = TRUE)$values) >= 0
.psd_via_shift <- function(K, eps = 1e-12) {
  if (!inherits(K, "dgCMatrix")) {
    K <- methods::as(methods::as(K, "generalMatrix"), "CsparseMatrix")
  }
  if (requireNamespace("RSpectra", quietly = TRUE)) {
    lam_min <- as.numeric(RSpectra::eigs_sym(K, k = 1, which = "SA")$values[1])
  } else {
    lam_min <- min(eigen(as.matrix(K), symmetric = TRUE, only.values = TRUE)$values)
  }
  if (is.finite(lam_min) && lam_min < 0) {
    K <- K + Matrix::Diagonal(nrow(K), x = (-lam_min) + eps)
  }
  K
}

#' Unit-level interference kernel K3+
#'
#' PSD part of the unit-level interference kernel
#' \deqn{(K_3)_{ij,i'j'} = 1\{ N_{ij} \cap N_{i'j'} \neq \emptyset \},}
#' i.e. two units are dependent when they share at least one interference
#' neighbour. Used by the cluster-agnostic CRN estimator (Theorem var_2).
#'
#' @param A n_unit x n_unit adjacency matrix (with self-loops).
#' @param psd_fun PSD projection: \code{"shift"} (default, sparse) or
#'   \code{"eigen"} (dense, exact PSD part).
#' @return An n_unit x n_unit PSD kernel matrix.
#' @keywords internal
#' @examples
#' ## chain network: each unit interferes with itself and adjacent units
#' A <- Matrix::Matrix(outer(1:6, 1:6, function(i, j) abs(i - j) <= 1) * 1,
#'                     sparse = TRUE)
#' K3 <- generalClusterExp:::.kernel_K3_plus(A)
#' round(as.matrix(K3), 2)  # overlap indicator (plus a small diagonal shift)
.kernel_K3_plus <- function(A, psd_fun = "shift") {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  K <- Matrix::tcrossprod(A)              # K_{ii'} = #shared neighbours
  if (length(K@x)) K@x[] <- 1             # binarise -> overlap indicator
  if (psd_fun == "shift") {
    .psd_via_shift(K)
  } else if (psd_fun == "eigen") {
    .psd_via_eigen(K)
  } else {
    stop("Unsupported psd_fun: ", psd_fun)
  }
}

#' Cluster kernels for the general (MRN / IPTW) variance estimator
#'
#' Builds the two cluster kernels of Theorem var_1:
#' \itemize{
#'   \item \code{K1} (within): \eqn{1\{i = i'\}} at the cluster level, i.e. two
#'     units are dependent when they belong to the same cluster (PSD);
#'   \item \code{K2} (overlap): \eqn{1\{ N^{cl}_{ij} \cap N^{cl}_{i'j'} \neq
#'     \emptyset \}}, i.e. their cluster-level interference neighbourhoods
#'     overlap.
#' }
#'
#' @param A n_unit x n_unit adjacency matrix (with self-loops).
#' @param C integer vector of cluster ids.
#' @return A list with sparse 0/1 kernels \code{K1} (within) and \code{K2}
#'   (overlap).
#' @keywords internal
#' @examples
#' A <- Matrix::Matrix(outer(1:6, 1:6, function(i, j) abs(i - j) <= 1) * 1,
#'                     sparse = TRUE)
#' C <- rep(1:3, each = 2)
#' ker <- generalClusterExp:::.cluster_kernels(A, C)
#' as.matrix(ker$K1)  # 1 iff same cluster
#' as.matrix(ker$K2)  # 1 iff overlapping cluster neighbourhoods
.cluster_kernels <- function(A, C) {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  n <- length(C)
  g <- as.integer(factor(C)); G <- max(g)

  # within-cluster kernel
  Zc <- Matrix::sparseMatrix(i = seq_len(n), j = g, x = 1, dims = c(n, G))
  K_within <- Matrix::tcrossprod(Zc)
  if (length(K_within@x)) K_within@x[] <- 1

  # overlapping cluster-neighbourhood kernel
  ii <- integer(0); jj <- integer(0)
  for (i in seq_len(n)) {
    clus <- unique(g[which(A[, i] != 0)])
    ii <- c(ii, rep.int(i, length(clus)))
    jj <- c(jj, clus)
  }
  B <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(n, G))
  K_overlap <- Matrix::tcrossprod(B)
  if (length(K_overlap@x)) K_overlap@x[] <- 1

  list(K1 = K_within, K2 = K_overlap)
}

# ---- matrix combination operators ------------------------------------------

#' HAC covariance from a kernel
#'
#' Computes the HAC covariance estimator \eqn{\Sigma^K = V^\top K V} (d x d,
#' symmetrised) from the n_unit x d contribution matrix \code{V} and a kernel
#' \code{K}.
#'
#' @param V n_unit x d matrix of per-unit, per-estimand contributions.
#' @param K n_unit x n_unit kernel matrix.
#' @return A d x d symmetric covariance matrix.
#' @keywords internal
#' @examples
#' set.seed(1)
#' V <- cbind(rnorm(6), rnorm(6)) / 6  # contributions of two estimands
#' K <- diag(6)                        # independence kernel
#' generalClusterExp:::.sigma_from_kernel(V, K)  # equals crossprod(V)
.sigma_from_kernel <- function(V, K) {
  S <- as.matrix(Matrix::crossprod(V, K %*% V))
  (S + t(S)) / 2
}

#' Loewner maximum of two symmetric matrices
#'
#' Computes \eqn{A \vee B = A + (B - A)_+}, the matrix generalisation of
#' \code{max()}: the result dominates both \code{A} and \code{B} in the PSD
#' (Loewner) order.
#'
#' @param A,B symmetric matrices of equal dimension.
#' @return A symmetric matrix, the Loewner maximum of \code{A} and \code{B}.
#' @keywords internal
#' @examples
#' A <- diag(c(2, 1)); B <- diag(c(1, 2))
#' generalClusterExp:::.lowner_max(A, B)  # here simply diag(2, 2)
.lowner_max <- function(A, B) {
  A <- (A + t(A)) / 2; B <- (B + t(B)) / 2
  D <- B - A
  eg <- eigen(D, symmetric = TRUE)
  Dp <- eg$vectors %*% (pmax(eg$values, 0) * t(eg$vectors))
  A + (Dp + t(Dp)) / 2
}

#' Bias-correction matrix under complete randomization
#'
#' Bias-correction matrix \eqn{\hat M} under cluster-level (stratified)
#' complete randomization (Section 6.3, bias-corrected estimator):
#' \deqn{\hat M = \sum_{k : p_k \in (p_{overlap}, 1 - p_{overlap})}
#'       \hat s_k \hat s_k^\top / (|I_k| p_k (1 - p_k)),}
#' \deqn{\hat s_k = \sum_{(i,j)} (T_{ij,k} - m_{ij,k} p_k) V_{ij},}
#' with \eqn{T_{ij,k}} the number of treated clusters of stratum k in
#' \eqn{N^{cl}_{ij}}, \eqn{m_{ij,k}} the number of clusters of stratum k in
#' \eqn{N^{cl}_{ij}}, and \eqn{p_k} the realised treated fraction of stratum k.
#'
#' @param A n_unit x n_unit adjacency matrix (with self-loops).
#' @param C integer vector of cluster ids.
#' @param W_C 0/1 vector of cluster-level treatments, indexed by cluster id.
#' @param V n_unit x d contribution matrix.
#' @param strata list of cluster-id vectors; \code{NULL} means one stratum.
#' @param p_overlap strata with \eqn{p_k} outside
#'   \eqn{(p_{overlap}, 1 - p_{overlap})} are skipped.
#' @return A d x d PSD correction matrix.
#' @keywords internal
#' @examples
#' set.seed(1)
#' A <- Matrix::Matrix(outer(1:12, 1:12, function(i, j) abs(i - j) <= 1) * 1,
#'                     sparse = TRUE)
#' C <- rep(1:4, each = 3)
#' W_C <- c(1, 0, 1, 0)  # exactly 2 of 4 clusters treated
#' V <- matrix(rnorm(12 * 2), 12, 2) / 12
#' generalClusterExp:::.bias_correction_matrix(A, C, W_C, V)
.bias_correction_matrix <- function(A, C, W_C, V, strata = NULL, p_overlap = 0) {
  if (!inherits(A, "sparseMatrix")) A <- Matrix::Matrix(A, sparse = TRUE)
  N <- length(C); d <- ncol(V)
  cluster_ids <- sort(unique(C))
  if (is.null(strata)) strata <- list(cluster_ids)
  K <- length(strata)

  max_cl <- max(cluster_ids)
  stratum_of <- integer(max_cl)
  for (k in seq_len(K)) stratum_of[strata[[k]]] <- k

  Ik        <- vapply(strata, length, integer(1))
  treated_k <- vapply(strata, function(s) sum(W_C[s] == 1), integer(1))
  p_k       <- treated_k / Ik
  active_k  <- p_k > p_overlap & p_k < 1 - p_overlap

  adj <- lapply(seq_len(N), function(i) which(A[, i] != 0))

  s_hat <- matrix(0, d, K)
  for (i in seq_len(N)) {
    vi <- V[i, ]
    if (all(vi == 0)) next
    clus <- unique(C[adj[[i]]])
    if (!length(clus)) next
    ks <- stratum_of[clus]; tc <- W_C[clus]
    in_k <- ks[ks > 0]
    if (!length(in_k)) next
    for (k in unique(in_k)) {
      sel  <- ks == k
      m_ik <- sum(sel)
      T_ik <- sum(tc[sel] == 1)
      s_hat[, k] <- s_hat[, k] + (T_ik - m_ik * p_k[k]) * vi
    }
  }

  M <- matrix(0, d, d)
  for (k in seq_len(K)) {
    if (!active_k[k]) next
    M <- M + (s_hat[, k] %o% s_hat[, k]) / (Ik[k] * p_k[k] * (1 - p_k[k]))
  }
  M
}
