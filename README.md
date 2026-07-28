# generalClusterExp

Estimators for cluster-randomized experiments with cross-cluster network
interference. The package implements the methods of

> Yihui He and Eric J. Tchetgen Tchetgen (2026).
> *A General Exposure-Mapping-Agnostic Framework for Causal Inference under
> Interference.* [arXiv:2607.04644](https://arxiv.org/abs/2607.04644)

as three top-level functions:

| function | weights | reference | variance estimator |
|----------|---------|-----------|--------------------|
| `mrn()`  | marginal Radon–Nikodym | Thm 3.6, Prop. 3.9(i)  | `Sigma1` Löwner-max (Thm 5.4) / `Sigma3` bias-corrected (Thm 5.11(ii)) |
| `iptw()` | cluster-level IPTW     | Prop. 3.9(iv)          | `Sigma1` Löwner-max (Thm 5.4) / `Sigma3` bias-corrected (Thm 5.11(ii)) |
| `crn()`  | complete Radon–Nikodym | Thm 3.7, Prop. 3.9(ii) | `Sigma2` (`K3+` kernel, Thm 5.9) |

## Installation

```r
# install.packages("remotes")
remotes::install_github("yihuihe4Wharton/generalClusterExp")
```

## Output

Each function returns a list with, for the supplied data,

* `causal_effects` — overall (`CE^O`), direct (`CE^D(phi_1)`, `CE^D(phi_0)`),
  indirect (`CE^I(1)`, `CE^I(0)`) and total (`CE^T`) effects, each with an
  estimate, variance and standard error. Effects that are not well-defined for
  the given regimes are returned as `NA` with `well_defined = FALSE`.
* `potential_outcomes` — all well-defined marginal `Ybar(phi)` and non-marginal
  `Ybar(w, phi)` population average potential outcomes, with variances.
* `Sigma` — the estimand covariance matrix used to derive every variance above.
* `variance_type` — which estimator was applied.

## Example

```r
library(generalClusterExp)

## Small geometric cluster network: units live around ten cluster centres
## and interfere with all units within distance 1.5.
set.seed(1)
n_per <- 6; K <- 10; n <- n_per * K
C <- rep(seq_len(K), each = n_per)                      # cluster ids
centers <- cbind(rep(1:5, 2), rep(1:2, each = 5)) * 4   # cluster centres
loc <- centers[C, ] + matrix(runif(2 * n, -1, 1), n, 2) # unit locations
A <- Matrix::Matrix(as.matrix(dist(loc)) <= 1.5, sparse = TRUE) * 1

## Two-stage assignment: clusters are treated w.p. 0.7; units in treated
## clusters get W = 1 w.p. 0.5 (regime phi_1), units in control clusters
## w.p. 0.2 (regime phi_0).
W_C <- rbinom(K, 1, 0.7)
W_Y <- rbinom(n, 1, ifelse(W_C[C] == 1, 0.5, 0.2))

## Outcomes with neighbourhood interference.
Y <- as.vector(-1 + A %*% (2 * W_Y) +
                 W_Y * as.vector(A %*% (1 * W_Y)) + rnorm(n))

fit <- mrn(A, C, Y, W_C, W_Y,
           dist1 = "binom", params1 = list(prob = 0.5),
           dist0 = "binom", params0 = list(prob = 0.2),
           distC = "binom", paramsC = list(prob = 0.7))
fit$causal_effects
fit$potential_outcomes

## Cluster-level complete randomization (exactly 7 of 10 clusters treated)
## automatically switches to the bias-corrected variance estimator:
W_C2 <- as.integer(seq_len(K) %in% sample.int(K, 7))
W_Y2 <- rbinom(n, 1, ifelse(W_C2[C] == 1, 0.5, 0.2))
Y2 <- as.vector(-1 + A %*% (2 * W_Y2) +
                  W_Y2 * as.vector(A %*% (1 * W_Y2)) + rnorm(n))
fit_bc <- mrn(A, C, Y2, W_C2, W_Y2,
              params1 = list(prob = 0.5), params0 = list(prob = 0.2),
              distC = "hypergeom", paramsC = list(m = 7))
fit_bc$variance_type
```

See `?mrn`, `?iptw` and `?crn` for full documentation and further examples.

## Reference

All section, theorem and proposition numbers below refer to
[arXiv:2607.04644](https://arxiv.org/abs/2607.04644).

* **Estimands** — direct, indirect, total and overall causal effects, and the
  marginal / non-marginal population average potential outcomes: Section 3
  ("Setup and identification").
* **Identification weights** — Theorem 3.6 characterizes all unbiased weights,
  Theorem 3.7 the cluster-agnostic ones, and Proposition 3.9 lists the
  examples implemented here: (i) marginal Radon–Nikodym → `mrn()`,
  (ii) complete Radon–Nikodym → `crn()`, (iv) cluster-level IPTW → `iptw()`.
* **Estimators** — the Hájek linear-weighted estimators and the weight
  matrices `B`, `B_marg`: Section 4.
* **Variance estimators** — Section 5. Theorem 5.4 gives `Sigma1`, the Löwner
  maximum of the within-cluster (`K1`) and cross-cluster (`K2`) HAC kernels
  under independent Bernoulli randomization. Theorem 5.9 gives `Sigma2`, the
  cluster-agnostic estimator built from `K3+`. Theorem 5.11(ii) gives
  `Sigma3 = Sigma1 - M`, the bias-corrected estimator under cluster-level
  complete randomization.

Theorem 5.9 establishes the conservativeness of `Sigma2` under independent
Bernoulli randomization only. Under cluster-level complete randomization
(`distC = "hypergeom"`, or `design = "complete"`) `crn()` still returns
`Sigma2`, but issues a warning and appends
`"; heuristic under complete randomization"` to `variance_type`. Prefer
`mrn()` or `iptw()` under that design — their bias-corrected `Sigma3`
(Theorem 5.11(ii)) is justified there.
