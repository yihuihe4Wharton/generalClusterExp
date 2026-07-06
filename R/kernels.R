##############################################################
# HAC variance kernels and matrix operations (Section 6 of main.tex).
#
# A kernel K is an n_unit x n_unit symmetric 0/1 (or PSD-projected) matrix that
# selects which pairs of units are treated as dependent. Given the n_unit x d
# matrix V of per-unit, per-estimand contributions \hat V_ij (one column per
# estimand), the corresponding HAC covariance estimator is
#       Sigma^K = V^T K V        (d x d).
##############################################################

# ---- PSD projection helpers (from our_helper) -------------------------------

.psd_via_eigen <- function(K, tol = 1e-12) {
  Kd <- as.matrix(K); Kd <- (Kd + t(Kd)) / 2
  eg <- eigen(Kd, symmetric = TRUE)
  lam_pos <- pmax(eg$values, tol)
  Kpsd <- eg$vectors %*% (lam_pos * t(eg$vectors))
  (Kpsd + t(Kpsd)) / 2
}

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

# K3+ : PSD part of the unit-level interference kernel
#   (K_3)_{ij,i'j'} = 1{ N_ij  cap  N_i'j'  != empty },
# i.e. two units share at least one interference neighbour. Used by the
# cluster-agnostic CRN estimator (Theorem var_2). `psd_fun` selects the PSD
# projection: "shift" (default, sparse) or "eigen" (dense, exact PSD part).
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

# Cluster kernels for the general (mrn / iptw) estimators (Theorem var_1):
#   K_within   : 1{ i = i' }                      (same cluster; PSD)
#   K_overlap  : 1{ N^{cl}_ij  cap  N^{cl}_i'j' != empty }  (overlapping
#                cluster-level interference neighbourhoods)
# Returned with the manuscript's naming via $K1 (within) and $K2 (overlap).
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

# HAC covariance estimator Sigma^K = V^T K V (d x d, symmetrised).
.sigma_from_kernel <- function(V, K) {
  S <- as.matrix(Matrix::crossprod(V, K %*% V))
  (S + t(S)) / 2
}

# Loewner maximum  A \vee B = A + (B - A)_+  (matrix generalisation of max()).
# Guarantees A \vee B  >=  A  and  >=  B in the PSD order.
.lowner_max <- function(A, B) {
  A <- (A + t(A)) / 2; B <- (B + t(B)) / 2
  D <- B - A
  eg <- eigen(D, symmetric = TRUE)
  Dp <- eg$vectors %*% (pmax(eg$values, 0) * t(eg$vectors))
  A + (Dp + t(Dp)) / 2
}

# Bias-correction matrix \hat M under cluster-level complete randomization
# (Section 6.3, bias-corrected estimator):
#   M_hat   = sum_{k: p_k in (p_overlap, 1-p_overlap)} s_hat_k s_hat_k^T / (|I_k| p_k (1-p_k)),
#   s_hat_k = sum_{(i,j)} (T_ij,k - m_ij,k p_k) * V_ij,
# with T_ij,k = #treated clusters of stratum k in N^{cl}_ij,
#      m_ij,k = #clusters of stratum k in N^{cl}_ij,
#      p_k    = (#treated clusters in stratum k) / |I_k|.
# V is the n_unit x d contribution matrix, so s_hat_k and M are d-vectors /
# d x d matrices. `strata` is a list of cluster-id vectors; NULL = one stratum.
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
