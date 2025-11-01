#' Generate lower-triangular (including diagonal) index pairs for vectorizing a symmetric matrix
#'
#' For a \eqn{q \times q} symmetric matrix (e.g., the random-effects covariance matrix \eqn{D}),
#' return all index pairs \eqn{(i, j)} such that \eqn{i \ge j}, ordered by row.
#' This is useful for the \emph{vech} vectorization (lower triangle including the diagonal).
#'
#' @param dim integer(1). The dimension \eqn{q} of the covariance matrix \eqn{D}
#'   (i.e., \eqn{D \in \mathbb{R}^{q \times q}}).
#'
#' @return An integer matrix with two columns, where each row is a lower-triangular
#'   index pair \eqn{(i, j)} with \eqn{i \ge j}:
#' \describe{
#'   \item{[,1]}{Row index \eqn{i}.}
#'   \item{[,2]}{Column index \eqn{j}.}
#' }
vech.posi <- function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(":", 1, 1:dim)))

#' Compute per-visit probabilities for a single subject via a logistic model
#'
#' For a given subject, builds a per-visit design matrix
#' \eqn{V_i = [\mathbf{1},\ \text{lag}(y_i),\ V[\texttt{idx}, 3]]} and computes
#' probabilities using the logistic link \eqn{\mathrm{logit}^{-1}(V_i \psi)}.
#'
#' @param y.i numeric vector. The subject-specific response trajectory ordered by
#'   visit time; length \eqn{k}. The first lagged value is set to \eqn{0}.
#' @param psi numeric vector. Coefficient vector compatible with the columns of
#'   the constructed design \eqn{V_i}; typically length 3 corresponding to
#'   intercept, lagged response, and the third column of \code{V}.
#' @param V matrix. A design/covariate matrix with at least three columns; the
#'   third column \code{V[, 3]} is used as a time-varying covariate aligned with
#'   the subject’s rows indexed by \code{idx}.
#' @param idx integer vector. Row indices (length \eqn{k}) selecting the
#'   subject-specific rows in \code{V} that align with \code{y.i}.
#'
#' @details
#' The lagged response is constructed as \code{c(0, y.i[1:(k-1)])}, so the first
#' visit uses a zero placeholder for the lag. The linear predictor is
#' \eqn{\eta = V_i \psi}, and probabilities are \eqn{\exp(\eta) / (1 + \exp(\eta))}.
#' Ensure that \code{length(psi)} matches \code{ncol(V.i)} and that
#' \code{length(idx) == length(y.i)}.
#'
#' @return A numeric vector of length \eqn{k} with values in \eqn{(0, 1)},
#'   giving the per-visit probabilities for the subject.
wi.func <- function(y.i, psi, V, idx) {
  k <- length(y.i)
  fData.i <- y.i[1:(k - 1)]
  fData.i <- c(0, fData.i)
  bData.i <- y.i[1:k]
  V.i <- cbind(1, fData.i, V[idx, 3])

  eta <- V.i %*% psi
  wi <- exp(eta) / (1 + exp(eta))
  return(wi)
}

#' Log joint contribution of responses and dropout indicators for one subject
#'
#' Computes the log-likelihood contribution of a single subject under
#' (i) a multivariate normal model for the response vector \eqn{y_i}
#' with mean \eqn{X_i \beta} and covariance \eqn{\Lambda_i}, and
#' (ii) independent Bernoulli observations for the dropout/response indicators
#' \eqn{R_i} with success probabilities \eqn{w_i}, i.e.
#' \deqn{\log p(y_i \mid X_i\beta, \Lambda_i) + \log p(R_i \mid w_i),}
#' where \eqn{\log p(R_i \mid w_i) = \sum_j \{R_{ij}\log w_{ij} + (1-R_{ij})\log (1-w_{ij})\}.}
#'
#' @param yi numeric vector of length \eqn{n_i}. Subject-specific response vector.
#' @param Ri numeric or integer vector of length \eqn{n_i}. Binary response/dropout
#'   indicators per visit (\code{0} or \code{1}); treated as independent
#'   Bernoulli given \code{wi}.
#' @param wi numeric vector of length \eqn{n_i}. Per-visit probabilities in \eqn{(0,1)}
#'   used for the Bernoulli part.
#' @param Xi numeric matrix of dimension \eqn{n_i \times p}. Design matrix for the
#'   mean model; must be conformable with \code{Beta}.
#' @param Lam.i numeric symmetric positive-definite matrix of dimension
#'   \eqn{n_i \times n_i}. Covariance matrix for \eqn{y_i}.
#' @param Beta numeric vector of length \eqn{p}. Fixed-effect coefficients.
#'
#' @details
#' The multivariate normal density is evaluated via
#' \code{mvtnorm::dmvnorm(yi, mean = c(Xi \%*\% Beta), sigma = Lam.i)}.
#' The indicator contribution is computed as
#' \code{log(prod(wi^Ri) * prod((1 - wi)^(1 - Ri)))};
#' for better numerical stability, an equivalent expression is
#' \code{sum(Ri * log(wi) + (1 - Ri) * log1p(-wi))}.
#'
#' @return A numeric scalar: the log-likelihood contribution
#' \eqn{\log p(y_i \mid X_i\beta, \Lambda_i) + \log p(R_i \mid w_i)}.
log.pos.y <- function(yi, Ri, wi, Xi, Lam.i, Beta) {
  fy <- log(mvtnorm::dmvnorm(yi, mean = c(Xi %*% Beta), sigma = Lam.i))
  fr <- log(prod(wi^Ri) * prod((1 - wi)^(1 - Ri)))
  return(fy + fr)
}

#' Single–subject MH update for missing outcomes (selection model)
#'
#' Performs one Metropolis–Hastings (MH) sweep over subjects with missing
#' outcomes to impute \eqn{y_{ij}} at the indices where \code{R == 1}, using a
#' Gaussian proposal built from the conditional mean/variance implied by the
#' working linear mixed model.
#'
#' @param y.c numeric vector of length \eqn{n}. Current completed outcome
#'   vector (observed \emph{and} imputed values).
#' @param R integer/logical vector of length \eqn{n}. Indicator of missingness
#'   at each observation time within subject; here \code{1} marks entries to be
#'   imputed within a subject block, \code{0} otherwise.
#' @param X numeric matrix \eqn{n \times p}. Fixed-effects design aligned with
#'   \code{y.c}.
#' @param Beta numeric vector (length \eqn{p}). Fixed-effect coefficients.
#' @param y.na.ind integer vector of subject indices that contain at least one
#'   missing response.
#' @param TLam.i numeric matrix \eqn{n \times n}. Marginal covariance of
#'   \code{y.c} under the working LME model (block-diagonal by subject).
#' @param V numeric matrix used by the dropout (selection) model; rows align
#'   with \code{y.c}.
#' @param psi numeric vector of dropout-model coefficients used inside
#'   \code{wi.func()}.
#' @param TO list of length \eqn{N}. For subject \eqn{j}, \code{TO[[j]]} is a
#'   selection matrix that picks the observed components within the subject.
#' @param TM list of length \eqn{N}. For subject \eqn{j}, \code{TM[[j]]} is a
#'   selection matrix that picks the missing components within the subject.
#'
#' @details
#' This function relies on a global vector \code{cumsum.nj} to map each subject
#' \eqn{j} to its row-block \eqn{1{:}n_j} in \code{y.c}. It calls
#' \code{wi.func()} to compute per-time dropout probabilities and
#' \code{log.pos.y()} for the joint log-density (outcome + selection part).
#'
#' @return The updated numeric vector \code{y.c} (length \eqn{n}) after one MH
#'   sweep. Specifically, values at positions where \code{R == 1} for subjects
#'   listed in \code{y.na.ind} may be replaced by proposed draws.
#'
MH.y.miss <- function(y.c, R, X, Beta, y.na.ind, TLam.i, V, psi, TO, TM) {
  idx.c <- NULL
  for (j in y.na.ind)
  {
    if (j == 1) {
      idx <- 1:cumsum.nj[j]
    } else {
      idx <- (cumsum.nj[j - 1] + 1):cumsum.nj[j]
    }
    y.old <- y.c[idx]
    wi.old <- wi.func(y.old, psi, V, idx)
    log.pos.old <- log.pos.y(y.old, R[idx], wi.old, X[idx, ], TLam.i[idx, idx], Beta)

    mean.j <- (diag(TM[[j]]) %*% X[idx, ] %*% Beta)[R[idx] == 1, ]
    var.j <- (t(diag(TM[[j]])) %*% TLam.i[idx, idx] %*% t(t(diag(TM[[j]]))))[R[idx] == 1, R[idx] == 1]
    if (length(mean.j) == 1) {
      y.cc <- rnorm(1, mean = mean.j, sd = sqrt(var.j))
    } else {
      y.cc <- mvtnorm::rmvnorm(1, mean.j, var.j)
    }
    y.new <- y.old
    y.new[R[idx] == 1] <- y.cc
    wi.new <- wi.func(y.new, psi, V, idx)
    log.pos.new <- log.pos.y(y.new, R[idx], wi.new, X[idx, ], TLam.i[idx, idx], Beta)
    log.accept <- log.pos.new - log.pos.old
    if (log(runif(1)) > log.accept) {
      y.c[idx][R[idx] == 1] <- y.new[R[idx] == 1]
      idx.c <- c(idx.c, j)
    }
  }
  return(y.c)
}

#' Autocorrelation matrix for an AR(\eqn{P}) process
#'
#' Constructs the \eqn{n \times n} autocorrelation matrix implied by an
#' AR(\eqn{P}) model with coefficients \code{phi}.
#'
#' @param phi numeric vector of length \eqn{P}. AR coefficients
#'   (\eqn{\phi_1,\ldots,\phi_P}).
#' @param P integer(1). Model order; defaults to \code{length(phi)}.
#' @param dim integer(1). Desired matrix dimension \eqn{n}. If \eqn{n \le P},
#'   it is internally set to \eqn{P+1}.
#' @param sigma.wn numeric(1). White-noise standard deviation (scales the
#'   implied variance); used to derive \eqn{\gamma_0}.
#'
#' @return A numeric correlation matrix of size \eqn{n \times n}.
#'
arp.Ci <- function(phi, P = length(phi), dim = P + 1, sigma.wn = 1) {
  if (dim <= 0) {
    return(cat("Waring: dim must be positive integer. \n"))
  }
  n <- dim
  if (n <= P) {
    n <- P + 1
  }
  posi.mt <- abs(1:P - outer(rep(1, P), 0:P))
  coef.mt <- outer(rep(1, P), c(-1, phi))
  coef.rho <- matrix(0, P, P)
  for (i in 1:P) {
    for (j in 1:P) {
      coef.rho[i, j] <- sum(coef.mt[i, ][posi.mt[i, ] == j])
    }
  }
  rho <- -solve(coef.rho) %*% phi
  gamma0 <- sigma.wn^2 / (1 - sum(rho * phi))
  rho <- c(rho, rep(NA, (n - P)))
  for (i in (P + 1):n) {
    rho[i] <- sum(rev(phi) * rho[(i - P):(i - 1)])
  }
  rho <- c(rev(rho), 1, rho)
  corr.mt <- matrix(NA, n, n)
  for (i in 1:n) {
    corr.mt[i, ] <- rho[(n + 2 - i):(2 * n - i + 1)]
  }
  return(corr.mt[1:dim, 1:dim])
}

#' Companion-form helper matrix \eqn{H(\phi)} for AR(\eqn{P})
#'
#' Builds the banded matrix used in derivatives of the AR precision structure.
#'
#' @param phi numeric vector of AR coefficients (length \eqn{P}).
#' @param dim integer(1). Target dimension \eqn{n_i}.
#'
#' @return A numeric matrix of size \eqn{n_i \times (P + n_i)} holding the
#'   coefficients that map lagged values in the AR construction.
#'
HL <- function(phi, dim = ni) {
  P <- length(phi)
  tmp.para <- c(1, -phi)
  HL.para <- matrix(0, dim, P + dim)
  for (i in 1:dim)
  {
    HL.para[i, (i + P):i] <- tmp.para
  }
  return(HL.para)
}

#' Derivatives of AR precision with respect to \eqn{\phi}
#'
#' Computes the list of derivatives \eqn{\partial S(\phi)/\partial \phi_k} for
#' the AR(\eqn{P}) precision structure \eqn{S(\phi)} at dimension \eqn{n_i}.
#'
#' @param phi numeric vector (length \eqn{P}). AR coefficients.
#' @param dim integer(1). Subject-specific length \eqn{n_i}.
#'
#' @return A list of length \eqn{P}; each element is an \eqn{n_i \times n_i}
#'   numeric matrix corresponding to the derivative with respect to \eqn{\phi_k}.
#'
d.inv.Ci.fn <- function(phi, dim = ni) {
  P <- length(phi)
  HL.phi <- HL(phi, dim)
  Hp <- HL.phi[, 1:P]
  Lp <- HL.phi[, -(1:P)]

  dS <- as.list(numeric(P))
  posi <- outer(1:dim, 1:(P + dim), function(i, j, P) i + P - j, P = P)
  for (k in 1:P)
  {
    d.HL.k <- matrix(0, dim, P + dim)
    d.HL.k[posi == k] <- -1
    dH <- d.HL.k[, 1:P]
    dL <- d.HL.k[, -(1:P)]
    dS[[k]] <- t(dL) %*% Lp + t(Lp) %*% dL - dH %*% t(Hp) - Hp %*% t(dH)
  }
  return(dS)
}

#' Single-parameter derivative of AR correlation via precision identity
#'
#' Computes \eqn{C(\phi)\, \partial S(\phi)/\partial \phi_k \, C(\phi)} with a
#' leading minus sign, which equals \eqn{\partial C(\phi)/\partial \phi_k}.
#'
#' @param phi numeric vector of AR coefficients.
#' @param dim integer(1). Dimension \eqn{n_i}.
#' @param k integer(1). Index of the coefficient \eqn{\phi_k} to differentiate.
#'
#' @return An \eqn{n_i \times n_i} numeric matrix giving
#'   \eqn{\partial C/\partial \phi_k}.
#'
arp.Ci.dot <- function(phi, dim, k) -arp.Ci(phi, dim = dim) %*% d.inv.Ci.fn(phi, dim = dim)[[k]] %*% arp.Ci(phi, dim = dim)

#' Derivative of AR(\eqn{p}) correlation with respect to \eqn{\phi} (stacked)
#'
#' Returns \eqn{\partial C/\partial \phi} aggregated across all AR parameters
#' using the absolute lag structure.
#'
#' @param phi numeric vector of AR coefficients.
#' @param dim integer(1). Dimension \eqn{n_i}.
#'
#' @return An \eqn{n_i \times n_i} numeric matrix whose \eqn{(r,c)} element is
#'   \eqn{|r-c|\,\phi^{|r-c|-1}} for \eqn{|r-c|\ge 1} and 0 on the diagonal.
#'
Arp.Ci.dot <- function(phi, dim) {
  Ti <- 1:dim
  tem <- abs(outer(Ti, Ti, "-"))
  dot.CAR <- tem * (phi^(tem - 1))
  return(dot.CAR)
}

#' Derivative of compound-symmetric correlation w.r.t. \eqn{\phi}
#'
#' For CS correlation \eqn{C = (1-\phi)I + \phi J}, the derivative w.r.t.
#' \eqn{\phi} is a matrix of ones with zeros on the diagonal.
#'
#' @param phi numeric(1). CS correlation parameter (not used in computation of
#'   the derivative form).
#' @param dim integer(1). Dimension \eqn{n_i}.
#'
#' @return An \eqn{n_i \times n_i} numeric matrix of ones with a zero diagonal.
#'
cs.Ci.dot <- function(phi, dim) {
  dot.CS <- matrix(1, dim, dim)
  diag(dot.CS) <- 0
  return(dot.CS)
}

#' Derivative of BAND1 correlation w.r.t. the band parameter
#'
#' Computes the derivative matrix for the first-order banded correlation (ones
#' on the main diagonal, \code{phi} on the first off-diagonals, zeros elsewhere).
#'
#' @param phi numeric(1). BAND1 off-diagonal parameter (not used explicitly in
#'   the derivative’s 0/1 structure).
#' @param dim integer(1). Dimension \eqn{n_i}.
#'
#' @return An \eqn{n_i \times n_i} numeric matrix with ones on the first
#'   off-diagonals and zeros elsewhere; diagonal is zero.
#'
BAND1.Ci.dot <- function(phi, dim) {
  tmp <- c(rep(0, dim - 2), 1, 0, 1, rep(0, dim - 2))
  dot.BAND1 <- NULL
  for (i in 1:dim) dot.BAND1 <- rbind(dot.BAND1, tmp[1:dim + (dim - i)])
  return(dot.BAND1)
}

#' One-step MH imputation using conditional (observed–missing) Gaussian
#'
#' Performs a Metropolis–Hastings (MH) update for subjects with missing data,
#' proposing from the conditional distribution of the missing components given
#' the observed ones implied by the working LME covariance.
#'
#' @param y.c numeric vector of length \eqn{n}. Current completed outcome
#'   (observed and imputed).
#' @param R integer/logical vector of length \eqn{n}. 1 marks entries imputed
#'   within a subject block; 0 otherwise.
#' @param X numeric matrix \eqn{n \times p}. Fixed-effects design.
#' @param Beta numeric vector (length \eqn{p}). Fixed-effect coefficients.
#' @param y.na.ind integer vector of subject indices with at least one missing.
#' @param cumsum.nj integer vector of cumulative observation counts per subject.
#' @param TLam.i numeric matrix \eqn{n \times n}. Marginal covariance of \code{y.c}.
#' @param V numeric matrix. Design for the dropout (selection) model.
#' @param psi numeric vector. Dropout-model coefficients used by \code{wi.func()}.
#' @param TO,TM lists of length \eqn{N}. Subject-wise selection matrices for
#'   observed and missing components, respectively.
#'
#' @return Updated numeric vector \code{y.c} (length \eqn{n}) after one MH sweep.
MH.y.miss1 <- function(y.c, R, X, Beta, y.na.ind, cumsum.nj, TLam.i, V, psi, TO, TM) {
  idx.c <- NULL
  for (j in y.na.ind)
  {
    if (j == 1) {
      idx <- 1:cumsum.nj[j]
    } else {
      idx <- (cumsum.nj[j - 1] + 1):cumsum.nj[j]
    }

    y.old <- y.c[idx]
    wi.old <- wi.func(y.old, psi, V, idx)
    mean.j <- TM[[j]] %*% X[idx, ] %*% Beta + TM[[j]] %*% TLam.i[idx, idx] %*% t(TO[[j]]) %*% solve(TO[[j]] %*% TLam.i[idx, idx] %*% t(TO[[j]])) %*%
      (TO[[j]] %*% y.old - TO[[j]] %*% X[idx, ] %*% Beta)
    var.j <- TM[[j]] %*% TLam.i[idx, idx] %*% (diag(length(y.old)) - t(TO[[j]]) %*% solve(TO[[j]] %*% TLam.i[idx, idx] %*% t(TO[[j]])) %*% TO[[j]] %*% TLam.i[idx, idx]) %*% t(TM[[j]])
    log.pos.old.1 <- log.pos.y1(y.old, R[idx], wi.old, X[idx, ], TLam.i[idx, idx], Beta, mean.j, var.j)
    if (length(mean.j) == 1) {
      y.cc <- rnorm(1, mean = mean.j, sd = sqrt(var.j))
    } else {
      y.cc <- mvtnorm::rmvnorm(1, mean.j, var.j)
    }
    y.new <- y.old
    y.new[R[idx] == 1] <- y.cc
    wi.new <- wi.func(y.new, psi, V, idx)
    log.pos.new.1 <- log.pos.y1(y.new, R[idx], wi.new, X[idx, ], TLam.i[idx, idx], Beta, mean.j, var.j)
    log.accept <- log.pos.new.1 - log.pos.old.1
    if (log(runif(1)) > log.accept) {
      y.c[idx][R[idx] == 1] <- y.new[R[idx] == 1]
      idx.c <- c(idx.c, j)
    }
  }
  return(y.c)
}

#' Log posterior of outcomes under selection model (with proposal correction)
#'
#' Computes the log joint density of \eqn{y_i} (Gaussian outcome) and
#' response indicators \eqn{R_i} (logistic selection), minus the log-density
#' of the proposal used for missing elements.
#'
#' @param yi numeric vector (length \eqn{n_i}). Subject-level outcomes.
#' @param Ri integer/logical vector (length \eqn{n_i}). Response indicators.
#' @param wi numeric vector (length \eqn{n_i}). Dropout probabilities
#'   \eqn{P(R_{ij}=1 \mid \cdot)}.
#' @param Xi numeric matrix \eqn{n_i \times p}. Fixed-effects design.
#' @param Lam.i numeric covariance matrix \eqn{n_i \times n_i}.
#' @param Beta numeric vector (length \eqn{p}). Fixed effects.
#' @param mean.j numeric vector. Proposal mean for missing components.
#' @param var.j numeric matrix or scalar. Proposal covariance for missing
#'   components.
#'
#' @return A single numeric value: log posterior (up to a constant).
log.pos.y1 <- function(yi, Ri, wi, Xi, Lam.i, Beta, mean.j, var.j) {
  yij <- yi[Ri == 1]
  fy <- log(mvtnorm::dmvnorm(yi, mean = Xi %*% Beta, sigma = Lam.i))
  if (length(yi) == 1) {
    fy <- fy - log(dnorm(yij, mean = mean.j, sd = sqrt(var.j)))
  } else {
    fy <- fy - log(mvtnorm::dmvnorm(yij, mean = mean.j, sigma = var.j))
  }
  fr <- log(prod(wi[-1]^Ri[-1]) * prod((1 - wi[-1])^(1 - Ri[-1])))
  return(fy + fr)
}

#' Correlation/covariance generator for common longitudinal structures
#'
#' Builds an \eqn{n \times n} correlation matrix for the requested structure:
#' UNC (identity), DEC (distance–decay), CAR1, ARp, BAND1, or CS.
#'
#' @param phi numeric vector/scalar. Correlation parameter(s).
#' @param dim integer(1). Dimension \eqn{n}.
#' @param type character. One of \code{"UNC"}, \code{"DEC"}, \code{"CAR1"},
#'   \code{"ARp"}, \code{"BAND1"}, or \code{"CS"}.
#' @param Ti optional numeric vector of times (same length as \code{dim}) used by
#'   DEC/CAR1.
#' @param ga optional numeric(1). Additional decay exponent for DEC.
#'
#' @return Numeric \eqn{n \times n} correlation matrix.
cor.fn <- function(phi, dim, type = c("UNC", "DEC", "CAR1", "ARp", "BAND1", "CS"), Ti = NULL, ga = NULL) {
  if (type[1] == "UNC") {
    return(diag(dim))
  }
  if (type[1] == "DEC") {
    return(phi^((abs(outer(Ti, Ti, "-")))^ga))
  }
  if (type[1] == "CAR1") {
    return(phi^(abs(outer(Ti, Ti, "-"))))
  }
  if (type[1] == "ARp") {
    return(arp.Ci(phi, P = length(phi), dim = dim))
  }
  if (type[1] == "BAND1") tmp <- c(rep(0, dim - 2), phi, 1, phi, rep(0, dim - 2))
  if (type[1] == "CS") tmp <- c(rep(phi, dim - 1), 1, rep(phi, dim - 1))
  Cor <- NULL
  for (i in 1:dim) Cor <- rbind(Cor, tmp[1:dim + (dim - i)])
  return(Cor)
}

#' Q-function for correlation parameters (Student-\eqn{t} LME)
#'
#' Objective used in the CM-step to update correlation parameters
#' (and possibly \eqn{\gamma} for DEC) under a tLME with scale factors \eqn{\tau_i}.
#'
#' @param Phiga numeric vector. Contains \code{Phi} (and \code{ga} if DEC).
#' @param sigma numeric(1). Residual variance scale.
#' @param E.hat1 numeric matrix \eqn{n \times n}. Expected residual cross-product.
#' @param cumsum.ni integer vector of cumulative sizes per subject.
#' @param tau numeric vector of length \eqn{N}. Subject-wise scale factors.
#' @param N integer(1). Number of subjects.
#' @param ni integer vector (length \eqn{N}). Per-subject lengths.
#' @param Data data.frame providing \code{week} by subject.
#' @param cor.type character(1). Correlation structure label.
#'
#' @return Numeric scalar: half the penalized trace + log-determinant term.
phiga.t.fn <- function(Phiga, sigma, E.hat1, cumsum.ni, tau, N, ni, Data, cor.type) {
  if (cor.type == "DEC") {
    Phi <- Phiga[1]
    ga <- Phiga[2]
  } else {
    Phi <- Phiga[1]
    ga <- 1
  }
  sum1 <- sum2 <- 0
  for (i in 1:N) {
    if (i == 1) {
      idx1 <- 1:cumsum.ni[1]
    } else {
      idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    }
    Cor <- cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga)
    sum1 <- sum1 + log(det(Cor))
    sum2 <- sum2 + tau[i] * sum((solve(Cor) * E.hat1[idx1, idx1])) / sigma
  }
  opt.Qfn <- 0.5 * (sum1 + sum2)
  return(opt.Qfn)
}

#' Q-function for correlation parameters (Gaussian LME)
#'
#' Objective used in the CM-step to update correlation parameters
#' under a Gaussian LME.
#'
#' @param Phiga numeric vector. Contains \code{Phi} (and \code{ga} if DEC).
#' @param sigma numeric(1). Residual variance scale.
#' @param E.hat1 numeric matrix \eqn{n \times n}. Expected residual cross-product.
#' @param cumsum.ni integer vector of cumulative sizes per subject.
#' @param N integer(1). Number of subjects.
#' @param ni integer vector (length \eqn{N}). Per-subject lengths.
#' @param Data data.frame providing \code{week} by subject.
#'
#' @return Numeric scalar: half the penalized trace + log-determinant term.
phiga.fn <- function(Phiga, sigma, E.hat1, cumsum.ni, N, ni, Data) {
  if (cor.type == "DEC") {
    Phi <- Phiga[1]
    ga <- Phiga[2]
  } else {
    Phi <- Phiga[1]
    ga <- 1
  }
  sum1 <- sum2 <- 0
  for (i in 1:N) {
    if (i == 1) {
      idx1 <- 1:cumsum.ni[1]
    } else {
      idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    }
    Cor <- cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga)
    sum1 <- sum1 + log(det(Cor))
    sum2 <- sum2 + sum((solve(Cor) * E.hat1[idx1, idx1])) / sigma
  }
  opt.Qfn <- 0.5 * (sum1 + sum2)
  return(opt.Qfn)
}

#' Derivative of DEC correlation w.r.t. \eqn{\phi}
#'
#' For distance–decay correlation \eqn{C_{rs} = \phi^{|t_r-t_s|^\gamma}}, returns
#' \eqn{\partial C/\partial \phi}.
#'
#' @param phi numeric(1). Base decay parameter in (0,1).
#' @param ga numeric(1). Exponent \eqn{\gamma>0}.
#' @param Ti numeric vector of observation times.
#'
#' @return Numeric matrix of size \eqn{n \times n}.
DEC.dot.phi <- function(phi, ga, Ti) {
  tem <- abs(outer(Ti, Ti, "-"))^ga
  dot.CAR <- tem * (phi^(tem - 1))
  return(dot.CAR)
}

#' Derivative of DEC correlation w.r.t. \eqn{\gamma}
#'
#' Computes \eqn{\partial C/\partial \gamma} for
#' \eqn{C_{rs} = \phi^{|t_r-t_s|^\gamma}}.
#'
#' @param phi numeric(1). Base decay parameter in (0,1).
#' @param ga numeric(1). Exponent \eqn{\gamma>0}.
#' @param Ti numeric vector of observation times.
#'
#' @return Numeric matrix of size \eqn{n \times n}.
DEC.dot.ga <- function(phi, ga, Ti) {
  tem <- abs(outer(Ti, Ti, "-"))
  dot.ga.CAR <- log(phi^tem) * cor.fn(phi = phi, type = "DEC", Ti = Ti, ga = ga)
  return(dot.ga.CAR)
}

#' Simulate Gaussian LME longitudinal data
#'
#' Generates subject-level outcomes \eqn{Y_i \sim N(X_i\beta, Z_i D Z_i^\top + \sigma C_i)}
#' with chosen correlation \code{cor.type}.
#'
#' @param n integer(1). Number of subjects.
#' @param para list with elements \code{Beta}, \code{DD}, \code{sigma}, \code{Phi}.
#' @param cor.type character(1). Correlation structure passed to \code{cor.fn()}.
#' @param si integer(1). Per-subject time points \eqn{s_i}.
#' @param q integer(1). Random-effect dimension (number of columns of \code{Z_i}).
#'
#' @return A list with:
#' \describe{
#'   \item{Data}{data.frame in long format with \code{Var1}, \code{Time}, \code{Subject}, \code{treat}.}
#'   \item{Ymat}{\eqn{n \times s_i} matrix of subject trajectories.}
#'   \item{X}{Stacked fixed-effects design.}
#'   \item{Z}{Stacked random-effects design.}
#' }
gen.lmm <- function(n, para, cor.type, si, q) {
  Beta <- para$Beta
  DD <- para$DD
  sigma <- para$sigma
  Phi <- para$Phi
  N <- rmultinom(n = 1, size = n, prob = c(0.5, 0.5))
  treat <- rep((1:2), N) - 1
  Ymat <- X <- Z <- NULL
  for (i in 1:n)
  {
    Xi <- cbind(1, sqrt(c(1:si) - 1), rep(treat[i], si), sqrt(c(1:si) - 1) * rep(treat[i], si))
    Zi <- Xi[, 1:q]
    X <- rbind(X, Xi)
    Z <- rbind(Z, Zi)
    Ymat <- rbind(Ymat, mvtnorm::rmvnorm(1, Xi %*% Beta, Zi %*% DD %*% t(Zi) + sigma * cor.fn(Phi, dim = si, type = cor.type[1])))
  }
  resp <- as.vector(t(Ymat[, 1:si]))
  sim.data <- data.frame(Var1 = resp, Time = rep(1:si - 1, n), Subject = rep(1:n, each = si), treat = rep(treat, each = si))
  return(list(Data = sim.data, Ymat = Ymat, X = X, Z = Z))
}

#' Simulate longitudinal data
#'
#' Generates \eqn{t_\nu}-distributed mixed-effects outcomes with covariance
#' \eqn{Z_i D Z_i^\top + \sigma C_i}.
#'
#' @param n integer(1). Number of subjects.
#' @param para list with elements \code{Beta}, \code{DD}, \code{sigma}, \code{Phi}, \code{nu}.
#' @param cor.type character(1). Correlation structure.
#' @param si integer(1). Per-subject time points.
#' @param q integer(1). Random-effect dimension.
#'
#' @return Same structure as \code{gen.lmm()}.
gen.tlmm <- function(n, para, cor.type, si, q) {
  Beta <- para$Beta
  DD <- para$DD
  sigma <- para$sigma
  nu <- para$nu
  Phi <- para$Phi
  N <- rmultinom(n = 1, size = n, prob = c(0.5, 0.5))
  treat <- rep((1:2), N) - 1
  Ymat <- X <- Z <- NULL
  for (i in 1:n)
  {
    Xi <- cbind(1, sqrt(c(1:si) - 1), rep(treat[i], si), sqrt(c(1:si) - 1) * rep(treat[i], si))
    Zi <- Xi[, 1:q]
    X <- rbind(X, Xi)
    Z <- rbind(Z, Zi)
    Ymat <- rbind(Ymat, t(Xi %*% Beta) + mvtnorm::rmvt(1, Zi %*% DD %*% t(Zi) + sigma * cor.fn(Phi, dim = si, type = cor.type[1]), df = nu))
  }
  resp <- as.vector(t(Ymat[, 1:si]))
  sim.data <- data.frame(Var1 = resp, Time = rep(1:si - 1, n), Subject = rep(1:n, each = si), treat = rep(treat, each = si))
  return(list(Data = sim.data, Ymat = Ymat, X = X, Z = Z))
}

#' Inject MNAR dropout into Gaussian LME simulations (selection model)
#'
#' Applies a logistic selection model with linear predictor based on
#' lagged/current outcomes and covariates, then censors post-dropout responses.
#'
#' @param Data data.frame from \code{gen.lmm()} (long format).
#' @param Ymat \eqn{n \times s_i} outcome matrix.
#' @param si integer(1). Per-subject length.
#' @param alpha numeric vector of selection-model coefficients
#'   (intercept, covariates, lagged and current outcomes).
#'
#' @return A list with
#' \describe{
#'   \item{Ymat}{Original complete trajectories.}
#'   \item{Ymat.na}{Trajectories with post-dropout \code{NA}s.}
#'   \item{Data.na}{Long-format data with \code{Var1} containing \code{NA}.}
#'   \item{Pna}{Dropout probabilities per row.}
#'   \item{D}{Integer vector of first-dropout indices per subject (or \code{si}).}
#'   \item{X,Z,Data.sub}{Designs and truncated data up to dropout.}
#' }
add.MNAR.LMM <- function(Data, Ymat, si, alpha) {
  Data$y.c <- Data$Var1
  N <- nrow(Ymat)
  y1 <- y1na <- Ymat
  treat <- Data$treat[which(Data$Time == 0)]
  Vi <- cbind(1, treat[1], c(0, Ymat[1, 1:(si - 1)]), Ymat[1, 1:si])
  eta <- Vi %*% alpha
  wi <- exp(eta) / (1 + exp(eta))
  p.na <- NULL
  p.na <- c(p.na, wi)
  for (i in 2:n)
  {
    Vi <- cbind(1, treat[i], c(0, Ymat[i, 1:(si - 1)]), Ymat[i, 1:si])
    eta <- Vi %*% alpha
    wi <- exp(eta) / (1 + exp(eta))
    p.na <- c(p.na, wi)
  }
  p.na[Data$Time == 0] <- 0
  Data$miss <- rbinom(nrow(Data), 1, p.na)
  Ymat.na <- D <- Data.sub <- X <- Z <- NULL
  for (i in 1:n) {
    tmp <- Data[Data$Subject == i, ]
    dropout <- which(tmp$miss == 1)[1]
    D <- c(D, dropout)
    if (!is.na(dropout)) tmp[dropout:nrow(tmp), "Var1"] <- NA
    if (is.na(dropout)) dropout <- si
    Data[Data$Subject == i, "Var1"] <- tmp$Var1
    Ymat.na <- rbind(Ymat.na, tmp$Var1)


    Xi <- cbind(1, sqrt(c(1:si) - 1), rep(treat[i], si), sqrt(c(1:si) - 1) * rep(treat[i], si))
    Zi <- Xi[, 1:q]
    X <- rbind(X, Xi[1:dropout, ])
    Z <- rbind(Z, Zi[1:dropout, ])
    Data.sub <- rbind(Data.sub, tmp[1:dropout, ])
  }
  return(list(
    Ymat = Ymat, Ymat.na = Ymat.na, Data.na = Data, Pna = p.na, D = D,
    X = X, Z = Z, Data.sub = Data.sub
  ))
}


#' Joint log-density of outcomes and response (Gaussian selection model)
#'
#' Computes \eqn{\log f_Y(y_i) + \log f_R(R_i \mid y_i)} for one subject under a
#' Gaussian outcome and logistic selection model.
#'
#' @param yi numeric vector \eqn{n_i}. Outcomes.
#' @param Xbeta.i numeric vector \eqn{n_i}. Mean \eqn{X_i \beta}.
#' @param Lam.i covariance matrix \eqn{n_i \times n_i}.
#' @param Vi numeric matrix \eqn{n_i \times m}. Selection-model design.
#' @param R.i logical/integer vector \eqn{n_i}. Response indicators.
#' @param alpha numeric vector. Selection-model coefficients.
#' @param n.i integer(1). Subject length \eqn{n_i}.
#'
#' @return Numeric scalar: joint log-density.
yv.joint.log <- function(yi, Xbeta.i, Lam.i, Vi, R.i, alpha, n.i) {
  fy <- mvtnorm::dmvnorm(yi, mean = c(Xbeta.i), sigma = Lam.i, log = T)
  pvi <- c(exp(as.matrix(Vi) %*% alpha) / (1 + exp(as.matrix(Vi) %*% alpha)))
  fv <- log(prod(pvi^R.i) * prod((1 - pvi)^(1 - R.i)))
  return(fy + fv)
}

#' Single-subject MH update for MNAR imputation (Gaussian outcome)
#'
#' Proposes new missing values from a Gaussian centered at current values and
#' accepts/rejects using the joint (outcome + selection) log-density.
#'
#' @param yi numeric vector \eqn{n_i}. Current subject trajectory.
#' @param Xbeta.i numeric vector \eqn{n_i}. Mean \eqn{X_i \beta}.
#' @param Lam.i covariance matrix \eqn{n_i \times n_i}.
#' @param alpha numeric vector. Selection-model coefficients.
#' @param Vi numeric matrix \eqn{n_i \times m}. Selection design for the subject.
#' @param R.i integer/logical vector \eqn{n_i}. Response indicators.
#' @param n.i integer(1). Subject length.
#' @param mu.mo numeric vector. Conditional mean of missing given observed.
#' @param Sig.mm.o numeric matrix. Conditional covariance of missing given observed.
#' @param na.idx.i integer vector of indices for missing components within subject
#'   (or \code{NA} if none).
#' @param Data.miss.i data.frame. Subject rows from \code{Data}.
#' @param V.fun.i function. Builder for the subject-level selection design given \eqn{y_i}.
#' @param Sig.mm.o.MC numeric matrix. Proposal covariance used by the random-walk step.
#'
#' @return Updated numeric vector \eqn{n_i} for the subject.
MH.y.miss2 <- function(yi, Xbeta.i, Lam.i, alpha, Vi, R.i, n.i, mu.mo, Sig.mm.o, na.idx.i, Data.miss.i, V.fun.i, Sig.mm.o.MC) {
  yi.old <- yi.new <- yi
  Vi.old <- V.fun.i(y.i = yi.old, Data.miss.i)
  log.den.old <- yv.joint.log(yi.old, Xbeta.i, Lam.i, Vi.old, R.i, alpha, n.i)
  if (!is.na(na.idx.i[1])) {
    lna <- length(na.idx.i)
    if (lna == 1) ym.hat <- c(rnorm(1, mean = yi.old[which(R.i == 1)], sd = 2.4 * sqrt(Sig.mm.o.MC[na.idx.i, na.idx.i])))
    if (lna > 1) ym.hat <- c(mvtnorm::rmvnorm(1, yi.old[which(R.i == 1)], 2.4 * Sig.mm.o.MC[na.idx.i, na.idx.i]))
    yi.new[which(R.i == 1)] <- ym.hat
  }
  Vi.new <- V.fun.i(y.i = yi.new, Data.miss.i)
  log.den.new <- yv.joint.log(yi.new, Xbeta.i, Lam.i, Vi.new, R.i, alpha, n.i)
  log.accept <- log.den.new - log.den.old
  if (log(runif(1)) > log.accept) yi.new <- yi.old
  return(yi.new)
}


#' Joint log-density of outcomes and response (tLME model with dropout)
#'
#' Computes \eqn{\log f_Y(y_i)} under \eqn{N(X_i\beta,\; \Lambda_i/\tau_i)} and
#' \eqn{\log f_R(R_i \mid y_i)} under logistic selection.
#'
#' @param yi numeric vector \eqn{n_i}. Outcomes.
#' @param Xbeta.i numeric vector \eqn{n_i}. Mean \eqn{X_i \beta}.
#' @param Lam.i covariance matrix \eqn{n_i \times n_i}.
#' @param tau.i numeric(1). Subject scale factor \eqn{\tau_i}.
#' @param Vi numeric matrix \eqn{n_i \times m}. Selection design.
#' @param R.i integer/logical vector \eqn{n_i}. Response indicators.
#' @param alpha numeric vector. Selection-model coefficients.
#' @param n.i integer(1). Subject length.
#'
#' @return Numeric scalar: joint log-density.
yv.joint.log.t <- function(yi, Xbeta.i, Lam.i, tau.i, Vi, R.i, alpha, n.i) {
  fy <- mvtnorm::dmvnorm(yi, mean = c(Xbeta.i), sigma = Lam.i / tau.i, log = T)
  pvi <- c(exp(Vi %*% alpha) / (1 + exp(Vi %*% alpha)))
  if (sum(c(exp(Vi %*% alpha)) == Inf) > 0) pvi[which(c(exp(Vi %*% alpha)) == Inf)] <- 1
  if (sum(pvi == 1) > 0) pvi[which(pvi == 1)] <- 1 - 1e-16
  if (sum(pvi == 0) > 0) pvi[which(pvi == 0)] <- 1e-16
  fv <- log(prod(pvi^R.i) * prod((1 - pvi)^(1 - R.i)))
  return(fy + fv)
}

#' Single-subject MH update for MNAR imputation (t outcome)
#'
#' Metropolis–Hastings step for imputation under a tLME model, using a Gaussian
#' random-walk proposal on the missing block and the t-based joint target.
#'
#' @param yi numeric vector \eqn{n_i}. Current subject trajectory.
#' @param Xbeta.i numeric vector \eqn{n_i}. Mean \eqn{X_i \beta}.
#' @param Lam.i covariance matrix \eqn{n_i \times n_i}.
#' @param tau.i numeric(1). Subject scale factor.
#' @param alpha numeric vector. Selection-model coefficients.
#' @param Vi numeric matrix \eqn{n_i \times m}. Selection design.
#' @param R.i integer/logical vector \eqn{n_i}. Response indicators.
#' @param n.i integer(1). Subject length.
#' @param mu.mo,Sig.mm.o,Sig.mm.o.MC as in \code{MH.y.miss2()}.
#' @param na.idx.i integer vector (indices of missing components) or \code{NA}.
#' @param Data.miss.i data.frame. Subject rows from \code{Data}.
#' @param V.fun.i function. Builder for subject-level design given \eqn{y_i}.
#'
#' @return Updated numeric vector \eqn{n_i} for the subject.
MH.y.miss.t <- function(yi, Xbeta.i, Lam.i, tau.i, alpha, Vi, R.i, n.i, mu.mo, Sig.mm.o, Sig.mm.o.MC, na.idx.i, Data.miss.i, V.fun.i) {
  yi.old <- yi.new <- yi
  Vi.old <- V.fun.i(y.i = yi.old, Data.miss.i)
  log.den.old <- yv.joint.log.t(yi.old, Xbeta.i, Lam.i, tau.i, Vi.old, R.i, alpha, n.i)
  if (!is.na(na.idx.i[1])) {
    lna <- length(na.idx.i)
    if (lna == 1) ym.hat <- c(rnorm(1, mean = yi.old[which(R.i == 1)], sd = 2.4 * sqrt(Sig.mm.o.MC[na.idx.i, na.idx.i])))
    if (lna > 1) ym.hat <- c(mvtnorm::rmvnorm(1, yi.old[which(R.i == 1)], 2.4 * Sig.mm.o.MC[na.idx.i, na.idx.i]))
    yi.new[which(R.i == 1)] <- ym.hat
  }
  Vi.new <- V.fun.i(y.i = yi.new, Data.miss.i)
  log.den.new <- yv.joint.log.t(yi.new, Xbeta.i, Lam.i, tau.i, Vi.new, R.i, alpha, n.i)
  log.accept <- log.den.new - log.den.old
  if (log(runif(1)) > log.accept) yi.new <- yi.old
  return(yi.new)
}

#' Blocked MH update for all missing values (t outcome, single proposal)
#'
#' Draws a single multivariate Gaussian proposal for all missing entries across
#' subjects, then accepts/rejects subject-wise using the tLME model.
#'
#' @param y.samp numeric matrix \eqn{M \times n}. Current MC sample bank;
#'   row \code{m-1} is updated to row \code{m}.
#' @param Data.miss data.frame with alignment columns and covariates.
#' @param y.na.ind integer vector of subject indices with missingness.
#' @param Sig.mm.o.MC numeric covariance for the global proposal on missing.
#' @param Xbeta numeric vector \eqn{n}. Mean \eqn{X\beta}.
#' @param TLam numeric covariance \eqn{n \times n}.
#' @param tau numeric vector (length \eqn{N}). Subject scale factors.
#' @param R integer/logical vector (length \eqn{n}). Response indicators.
#' @param alpha numeric vector. Selection coefficients.
#' @param ni integer vector (length \eqn{N}). Per-subject sizes.
#' @param cumsum.ni integer cumulative indices.
#' @param m integer(1). Target MCMC iteration index.
#'
#' @return Updated row vector (length \eqn{n}) for iteration \code{m}.
MH.y.miss.t.new <- function(y.samp, Data.miss, y.na.ind, Sig.mm.o.MC, Xbeta, TLam, tau, R, alpha, ni, cumsum.ni, m) {
  y.old <- y.new <- y.samp[(m - 1), ]
  ym.hat <- mvtnorm::rmvnorm(1, y.old[which(R == 1)], 2.4 * Sig.mm.o.MC)
  y.new[which(R == 1)] <- ym.hat
  V.old <- V.fun(y.old, Data.miss)
  V.new <- V.fun(y.new, Data.miss)
  for (i in y.na.ind)
  {
    if (i == 1) {
      idx1 <- 1:cumsum.ni[1]
    } else {
      idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    }
    log.den.old <- yv.joint.log.t(y.old[idx1], c(Xbeta)[idx1], TLam[idx1, idx1], tau[i], V.old[idx1, ], R[idx1], alpha, ni[i])
    log.den.new <- yv.joint.log.t(y.new[idx1], c(Xbeta)[idx1], TLam[idx1, idx1], tau[i], V.new[idx1, ], R[idx1], alpha, ni[i])
    log.accept <- log.den.new - log.den.old
    if (log(runif(1)) > log.accept) y.new[idx1] <- y.old[idx1]
  }
  sum((y.new - y.old) != 0) / length(y.na.ind)
  return(y.new)
}



#' Q-function for the degrees of freedom \eqn{\nu} (tLME model)
#'
#' Minimization target to update \eqn{\nu} in the CM-step using expectations
#' \eqn{\kappa_i} and \eqn{\tau_i}.
#'
#' @param nu numeric(1). Degrees of freedom (> 2).
#' @param kappa numeric vector (length \eqn{N}). Expectation terms.
#' @param tau numeric vector (length \eqn{N}). Scale factors.
#' @param N integer(1). Number of subjects.
#'
#' @return Negative Q-function value to be minimized.
nu.Q.fn <- function(nu, kappa, tau, N) {
  fn.nu <- nu * (kappa[1] - tau[1] + log(nu / 2)) - 2 * log(gamma(nu / 2))
  for (i in 2:N) {
    fn.nu <- fn.nu + (nu * (kappa[i] - tau[i] + log(nu / 2)) - 2 * log(gamma(nu / 2)))
  }
  return(-fn.nu)
}

#' One-step prediction of missing outcomes under a fitted LME
#'
#' Given an \code{nlme::lme} fit, constructs \eqn{D}, \eqn{\sigma}, and
#' correlation parameters, then predicts missing responses using the Gaussian
#' conditional expectation formula.
#'
#' @param Data long-format data with \code{Var1}, \code{R}, \code{week}, \code{Subject},
#'   and potential covariates used by the auxiliary GLM step.
#' @param fm2 an \code{nlme::lme} fitted model.
#' @param X,Z design matrices aligned with \code{Data}.
#' @param cor.type character(1). Correlation structure to reconstruct \eqn{C_i}.
#'
#' @return A list with
#' \describe{
#'   \item{alpha.hat}{Coefficients from the auxiliary logistic regression for \code{R}.}
#'   \item{yc}{Completed response vector with missing imputed by conditional mean.}
#' }
prediction_ym <- function(Data, fm2, X, Z, cor.type = "UNC") {
  Beta <- matrix(c(fm2$coefficients$fixed), ncol = 1)
  DD <- matrix(0, q, q)
  DD[1, 1] <- as.numeric(VarCorr(fm2)[1, 1])
  DD[2, 2] <- as.numeric(VarCorr(fm2)[2, 1])
  DD[1, 2] <- DD[2, 1] <- as.numeric(VarCorr(fm2)[2, 3]) * sqrt(as.numeric(VarCorr(fm2)[1, 1]) * as.numeric(VarCorr(fm2)[2, 1]))
  sigma <- c(fm2$sigma)
  Phi <- runif(1)

  N <- length(unique(Data$Subject))
  na.ind <- which(is.na(as.vector(t(Data$Var1)))) ## which time point have missing value
  Data.miss <- Data

  y <- Data.miss$Var1
  R <- Data.miss$R

  yo <- y[-na.ind]
  no <- length(yo)
  Xo <- X[-na.ind, ]
  Xm <- X[na.ind, ]

  ni <- numeric(N)
  for (i in 1:N) ni[i] <- length(Data$Subject[Data$Subject == i])
  cumsum.ni <- cumsum(ni)
  ni.o <- numeric(N)
  for (i in 1:N) ni.o[i] <- sum(!is.na(Data$Var1[Data$Subject == i]))
  cumsum.ni.o <- cumsum(ni.o)
  cumsum.q <- cumsum(rep(q, N))

  si <- max(ni)
  n <- sum(ni)

  y.na.ind <- unique(Data[na.ind, ]$Subject) ## who have missing value
  Nm <- length(y.na.ind)
  mi <- numeric(Nm)
  for (i in 1:Nm) mi[i] <- sum(Data.miss$R[Data.miss$Subject == y.na.ind[i]])
  cumsum.na <- cumsum(mi)
  num.na <- length(na.ind)

  TO <- diag(n)[-na.ind, ]
  TM <- diag(n)[na.ind, ]
  if (num.na == 1) TM <- t(TM)
  TZ <- matrix(0, ncol = N * q, nrow = n)
  TZ[1:cumsum.ni[1], 1:q] <- Z[1:cumsum.ni[1], ]
  for (i in 2:N) TZ[(cumsum.ni[i - 1] + 1):cumsum.ni[i], ((i - 1) * q + 1):(i * q)] <- Z[(cumsum.ni[i - 1] + 1):cumsum.ni[i], ]
  vechD <- vech.posi(q)
  TLam <- TLam.inv <- TCor <- TCor.inv <- matrix(0, ncol = n, nrow = n)
  for (i in 1:N) {
    if (i == 1) {
      idx1 <- 1:cumsum.ni[1]
    } else {
      idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    }
    Zi <- as.matrix(Z[idx1, ])
    TLam[idx1, idx1] <- Zi %*% DD %*% t(Zi) + sigma * cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga)
    Cor <- cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga)
    TCor[idx1, idx1] <- Cor
    TCor.inv[idx1, idx1] <- solve(Cor)
    TLam.inv[idx1, idx1] <- solve(Zi %*% DD %*% t(Zi) + sigma * cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga))
  }
  TLam.oo <- TLam[-na.ind, -na.ind]
  TLam.mo <- TLam[na.ind, -na.ind]
  TLam.mm <- TLam[na.ind, na.ind]
  if (num.na == 1) TLam.mo <- t(TLam.mo)
  TLam.oo.inv <- matrix(0, ncol = sum(ni.o), nrow = sum(ni.o))
  for (i in 1:N) {
    if (i == 1) {
      idx1 <- 1:cumsum.ni.o[1]
    } else {
      idx1 <- (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]
    }
    TLam.oo.inv[idx1, idx1] <- solve(TLam.oo[idx1, idx1])
  }

  yo.cent <- yo - Xo %*% Beta
  ym <- Xm %*% Beta + TLam.mo %*% TLam.oo.inv %*% yo.cent
  Data$yyo <- Data$Var1
  Data$yyo[na.ind] <- ym

  V.fun <- function(y, Data.miss) {
    V <- NULL
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }
      Data.i <- Data.miss[which(Data.miss$Subject == i), ]
      week <- Data.i$week
      k <- length(week)
      fData.i <- y[1:(k - 1)]
      fData.i <- y[idx1][-k]
      fData.i <- c(0, fData.i)
      bData.i <- y[idx1]
      V.i <- cbind(1, Data.i$gender, Data.i$drug, Data.i$prevOI, Data.i$AZT, fData.i, bData.i)
      V <- rbind(V, V.i)
    }
    return(V)
  }
  y <- Data$yyo
  Data.a <- V.fun(y, Data)
  Data.a <- as.data.frame(Data.a)
  colnames(Data.a) <- c("int.", "gender", "drug", "prevOI", "AZT", "fData.i", "bData.i")

  fm3 <- glm(R ~ prevOI + fData.i + bData.i, data = Data.a, family = binomial())
  summary(fm3)
  alpha <- fm3$coefficients
  return(list(alpha.hat = alpha, yc = y))
}
