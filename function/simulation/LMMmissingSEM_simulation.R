################################################################################
#                                                                                     
#   Filename    :    LMMmissingSEM.R  												  
#   Project     :    BiomJ article "Extending t linear mixed models for longitudinal 
#                    data with non-ignorable dropout applied to AIDS studies"                                                           
#   Authors     :    Yu-Chen Yang and Wan-Lun Wang and Luis M. Castro and Tsung-I Lin
#   Date        :    07.08.2025
#   Purpose     :    Function allowing to fit a LME model with non-ignorable dropout 
#                    under the MCAR, MAR, and MNAR mechanisms  
#
#   Input data files : A standard data frame (named 'Data') containing subject ID, 
#                      time points, observed responses, missingness indicators, 
#                      and additional covariates (if applicable). 
#   Output data files : an object colleting the fitting results of 
#                       linear mixed-effects (LME) model with  
#                       selection modeling-based MCAR, MAR, and MNAR mechanisms 
#
#   R Version   :    R-4.3.1                                                              
#   Required R packages : nlme; mvtnorm  
#
################################################################################ 
#' Fit an LME model with non-ignorable dropout via a selection model (MCAR/MAR/MNAR)
#'
#' Fits a Gaussian linear mixed-effects (LME) model for longitudinal outcomes
#' jointly with a selection-model dropout mechanism under MCAR, MAR, or MNAR.
#' The longitudinal covariance combines random-effect variation and a
#' within-subject correlation structure; the dropout model uses a logistic link
#' with a design \eqn{V} constructed internally according to the specified
#' mechanism.
#'
#' @param Data data.frame. Long-format longitudinal data. Must contain at least:
#'   \code{Subject} (ID), \code{Var1} (response; may contain \code{NA}),
#'   \code{R} (dropout/response indicator), \code{week} (time), and any
#'   covariates used by the dropout design (e.g., \code{prevOI}).
#' @param X matrix. Fixed-effects design matrix aligned with \code{Data}
#'   (\code{nrow(X) == nrow(Data)}). Columns correspond to fixed effects.
#' @param Z matrix. Random-effects design matrix aligned with \code{Data}.
#'   Its number of columns defines the random-effects dimension \eqn{q}.
#' @param init.para list. Initial parameters with components:
#'   \itemize{
#'     \item \code{Beta} (numeric length \eqn{p}): fixed effects.
#'     \item \code{DD} (numeric \eqn{q \times q} matrix): random-effects covariance \eqn{D}.
#'     \item \code{sigma} (numeric): residual scale multiplier.
#'     \item \code{alpha} (numeric): coefficients for the dropout logistic model.
#'     \item \code{Phi} (numeric): parameters of the within-subject correlation.
#'   }
#' @param cor.type character(1). Working correlation structure for within-subject
#'   errors. One of \code{"UNC"}, \code{"ARp"}, \code{"BAND1"}, \code{"CS"}.
#' @param M integer(1). Monte Carlo sample size for the data-augmentation chain.
#' @param M.LL integer(1). Monte Carlo sample size used when approximating the
#'   observed log-likelihood.
#' @param tol numeric(1). Tolerance used in convergence checks. Default \code{1e-6}.
#' @param max.iter integer(1). Maximum number of SEM iterations.
#' @param per integer(1). Print frequency (in iterations).
#' @param mechanism character(1). Dropout mechanism controlling the construction
#'   of the dropout-design \eqn{V}. One of \code{"MNAR"}, \code{"MAR"}, or \code{"MCAR"}.
#'
#' @details
#' For subject \eqn{i} with design matrices \eqn{X_i} and \eqn{Z_i}, the longitudinal
#' covariance is \eqn{Z_i D Z_i^\top + \sigma C_i(\Phi)}, where \eqn{C_i(\Phi)}
#' depends on \code{cor.type} (e.g., AR(p), CS, banded). The dropout probability
#' uses a logistic link with a per-row design \eqn{V}, built from \code{Data} (and
#' lagged responses under MNAR). Missing responses are imputed via an MCMC kernel
#' within a stochastic EM loop; parameters are updated by conditional maximization.
#' The function prints progress information and returns estimates, observed-information
#' components, and Monte Carlo diagnostics. Requires \pkg{mvtnorm}.
#'
#' @return A list with components:
#' \describe{
#'   \item{model.inf}{List with model diagnostics:
#'     \code{loglik} (final observed log-likelihood),
#'     \code{iter.lnL} (log-likelihood trace),
#'     \code{aic}, \code{bic}, and \code{time} (time to compute the information matrix).}
#'   \item{para.est}{List of parameter estimates:
#'     \code{Beta} (fixed effects), \code{sigma}, \code{D} (random-effects covariance),
#'     \code{Phi} (correlation parameters), \code{ga} (if applicable),
#'     \code{alpha} (dropout-model coefficients), and \code{b} (posterior means of random effects).}
#'   \item{iter}{Integer. Number of SEM iterations performed.}
#'   \item{y.c}{Numeric vector of length \eqn{n}. Completed/augmented responses.}
#'   \item{MSE.y}{Numeric scalar. Mean squared error against \code{Data$y.c} (if present).}
#'   \item{na.ind}{Integer vector. Row indices of missing responses in \code{Data$Var1}.}
#'   \item{Tpara}{Numeric matrix. Parameter trajectory across iterations (row \eqn{h} is iteration \eqn{h}).}
#'   \item{IM}{List. Output from \code{I.lmm.missing()} (observed/Fisher information estimators).}
#' }
LMM.miss.SEM <- function(Data, X, Z, V, g = g, init.para, cor.type = c("UNC", "ARp", "BAND1", "CS"), M = 100, M.LL = 1000, P = 1, tol = 1e-6, max.iter = max.iter, per = 1, mechanism = c("MNAR", "MAR", "MCAR")) {
  begin <- proc.time()[1]
  # initial values of parameter
  Beta <- init.para$Beta
  DD <- init.para$DD
  sigma <- init.para$sigma
  alpha <- init.para$alpha
  # if(mechanism=='MCAR') alpha = matrix(0)
  Phi <- init.para$Phi
  ga <- 1
  p <- ncol(X)
  q <- ncol(Z)
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

  dropout.idx <- NULL
  for (i in 1:N)
  {
    if (i == 1) {
      idx1 <- 1:cumsum.ni[1]
    } else {
      idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    }
    # which(R[idx1] == 1)
    dropout.idx <- c(dropout.idx, idx1[which(R[idx1] == 1)[-1]])
  }


  y.na.ind <- unique(Data[na.ind, ]$Subject) ## who have missing value
  Nm <- length(y.na.ind)
  mi <- numeric(Nm)
  for (i in 1:Nm) mi[i] <- sum(Data.miss$R[Data.miss$Subject == y.na.ind[i]])
  cumsum.na <- cumsum(mi)
  num.na <- length(na.ind)
  na.idx <- as.list(N)
  for (i in 1:N) na.idx[[i]] <- NA
  na.idx[[y.na.ind[[1]]]] <- 1:cumsum.na[1]
  for (i in 2:Nm) na.idx[[y.na.ind[i]]] <- (cumsum.na[i - 1] + 1):cumsum.na[i]

  if (mechanism == "MNAR") cat(rep("=", 25), "LMM (MNAR) with ", cor.type, " errors is fitted...; ", "missing = ", num.na / n * 100, "%", rep("=", 25), sep = "", "\n")
  if (mechanism == "MAR") cat(rep("=", 25), "LMM (MAR) with ", cor.type, " errors is fitted...; ", "missing = ", num.na / n * 100, "%", rep("=", 25), sep = "", "\n")
  if (mechanism == "MCAR") cat(rep("=", 25), "LMM (MCAR) with ", cor.type, " errors is fitted...; ", "missing = ", num.na / n * 100, "%", rep("=", 25), sep = "", "\n")

  TO <- diag(n)[-na.ind, ]
  TM <- diag(n)[na.ind, ]
  if (num.na == 1) TM <- t(TM)
  TZ <- matrix(0, ncol = N * q, nrow = n)
  TZ[1:cumsum.ni[1], 1:q] <- Z[1:cumsum.ni[1], ]
  for (i in 2:N) TZ[(cumsum.ni[i - 1] + 1):cumsum.ni[i], ((i - 1) * q + 1):(i * q)] <- Z[(cumsum.ni[i - 1] + 1):cumsum.ni[i], ]
  vechD <- vech.posi(q)
  TLam <- TCor <- TCor.inv <- matrix(0, ncol = n, nrow = n)
  for (i in 1:N) {
    if (i == 1) {
      idx1 <- 1:cumsum.ni[1]
    } else {
      idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    }
    Zi <- as.matrix(Z[idx1, ])
    Cor <- cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga)
    TLam[idx1, idx1] <- Zi %*% DD %*% t(Zi) + sigma * Cor
    TCor[idx1, idx1] <- Cor
    TCor.inv[idx1, idx1] <- solve(Cor)
  }
  # TLam.inv = solve(TLam)
  TLam.inv <- matrix(0, ncol = n, nrow = n)
  for (i in 1:N) {
    if (i == 1) {
      idx1 <- 1:cumsum.ni[1]
    } else {
      idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    }
    TLam.inv[idx1, idx1] <- solve(TLam[idx1, idx1])
  }

  TLam.oo <- TLam[-na.ind, -na.ind]
  TLam.mo <- TLam[na.ind, -na.ind]
  TLam.mm <- TLam[na.ind, na.ind]
  if (num.na == 1) TLam.mo <- t(TLam.mo)
  # TLam.oo.inv = solve(TLam.oo)
  TLam.oo.inv <- matrix(0, ncol = no, nrow = no)
  for (i in 1:N) {
    if (i == 1) {
      idx0 <- 1:cumsum.ni.o[1]
    } else {
      idx0 <- (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]
    }
    TLam.oo.inv[idx0, idx0] <- solve(TLam.oo[idx0, idx0])
  }
  Sig.mm.o.MC <- (TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo))


  # observed log-likelihood:
  Xbeta <- X %*% Beta
  y.samp <- matrix(rep(y, M), nrow = M, ncol = n, byrow = T)
  yo.hat <- yo
  y.hat <- t(TO) %*% yo + t(TM) %*% rep(mean(yo), num.na)
  y.samp[1, ] <- y.hat

  burn.in <- M / 2
  cho <- seq(1, burn.in, 5)
  mc.size <- length(cho)
  # mc.size = M/2

  ## A covance matrix for missing propabality
  if (mechanism == "MNAR") {
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
        # V.i = cbind(1, Data.i$treat[1:k], fData.i, bData.i)
        # V.i = cbind(1, Data.i$gender, Data.i$drug, Data.i$prevOI, Data.i$AZT, fData.i, bData.i)
        V.i <- cbind(1, Data.i$treat, fData.i, bData.i)
        V <- rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i <- function(y.i, Data.i) {
      k <- length(y.i)
      fData.i <- y.i[1:(k - 1)]
      fData.i <- y.i[-k]
      fData.i <- c(0, fData.i)
      bData.i <- y.i
      # V = cbind(rep(1, k), treat.i, fData.i, bData.i)
      # V = cbind(rep(1, k), Data.i$gender, Data.i$drug, Data.i$prevOI, Data.i$AZT, fData.i, bData.i)
      V <- cbind(rep(1, k), Data.i$treat, fData.i, bData.i)
      return(V)
    }
    # print(mechanism)
  }
  if (mechanism == "MAR") {
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
        V.i <- cbind(1, Data.i$treat, fData.i)
        V <- rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i <- function(y.i, Data.i) {
      k <- length(y.i)
      fData.i <- y.i[1:(k - 1)]
      fData.i <- y.i[-k]
      fData.i <- c(0, fData.i)
      bData.i <- y.i
      V <- cbind(rep(1, k), Data.i$treat, fData.i)
      return(V)
    }
  }

  if (mechanism == "MCAR") {
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
        V.i <- cbind(rep(1, length(idx1)))
        V <- rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i <- function(y.i, Data.i) {
      k <- length(y.i)
      V <- cbind(rep(1, k))
      return(V)
    }
  }
  V <- V.fun(y.hat, Data.miss)

  yo.cent <- yo.hat - Xo %*% Beta
  mu.mo <- Xm %*% Beta + TLam.mo %*% TLam.oo.inv %*% yo.cent
  Sig.mm.o <- TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo)

  #### log-likelihood
  wden <- numeric(N)
  for (i in 1:N)
  {
    if (i == 1) {
      idx1 <- 1:cumsum.ni[1]
    } else {
      idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    }

    if (sum(i == y.na.ind) == 0) {
      pvi <- c(exp(as.matrix(V[idx1, ]) %*% alpha) / (1 + exp(as.matrix(V[idx1, ]) %*% alpha)))
      wden[i] <- mvtnorm::dmvnorm(y[idx1], X[idx1, ] %*% Beta, TLam[idx1, idx1], log = F) * prod(1 - pvi)
    } else {
      na.i <- which(y.na.ind == i)
      y.hat.ll <- mu.mo[na.i] + mvtnorm::rmvnorm(M.LL, sigma = as.matrix(Sig.mm.o[na.i, na.i]))
      dt.ym <- mvtnorm::dmvnorm(matrix(y.hat.ll), mean = mu.mo[na.i], sigma = as.matrix(Sig.mm.o[na.i, na.i]), log = F)
      wden.i <- NULL
      for (jj in 1:(M.LL / 100))
      {
        y.hat.i <- c(y[idx1][R[idx1] == 0], y.hat.ll[jj])
        V.i <- V.fun.i(y.i = y.hat.i, Data.i = Data.miss[idx1, ])
        pvi <- c(exp(V.i %*% alpha) / (1 + exp(V.i %*% alpha)))
        fn <- mvtnorm::dmvnorm(y.hat.i, X[idx1, ] %*% Beta, TLam[idx1, idx1], log = F) * prod((1 - pvi)^(1 - R[idx1]) * pvi^(R[idx1]))
        wden.i <- c(wden.i, fn / dt.ym[jj])
      }
      wden[i] <- mean(wden.i)
    }
  }
  loglik.old <- iter.lnL <- sum(log(wden))

  theta.old <- c(Beta, DD[vechD], sigma, Phi, ga, alpha)
  iter <- 0
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("Linear mixed models with ", cor.type[1], " errors: ", "\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.old, sep = "", "\n")
  Tpara <- theta.old
  diff.lnL <- 10000
  diff <- 10000
  repeat  {
    iter <- iter + 1
    ######  E-Step:
    TD <- kronecker(diag(N), DD)
    # TSig.b = TD - TD %*% t(TZ) %*% TLam.inv %*% TZ %*% TD
    TSig.b <- matrix(0, N * q, N * q)
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
        idx2 <- 1:cumsum.q[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
        idx2 <- (cumsum.q[i - 1] + 1):cumsum.q[i]
      }
      TSig.b[idx2, idx2] <- DD - DD %*% t(TZ[idx1, idx2]) %*% TLam.inv[idx1, idx1] %*% TZ[idx1, idx2] %*% DD
    }

    ##### ##### ##### ##### ##### ##### ##### #####
    ##### genrate ym by MCMC
    yo.cent <- yo.hat - Xo %*% Beta
    mu.mo <- Xm %*% Beta + TLam.mo %*% TLam.oo.inv %*% yo.cent
    Sig.mm.o <- TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo)
    Xbeta <- X %*% Beta
    y.samp <- matrix(rep(y.hat, M), nrow = M, ncol = n, byrow = T)
    for (m in 2:M)
    {
      for (i in y.na.ind)
      {
        if (i == 1) {
          idx1 <- 1:cumsum.ni[1]
        } else {
          idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
        }
        y.samp[m, idx1] <- MH.y.miss2(
          yi = y.samp[(m - 1), idx1], Xbeta.i = c(Xbeta)[idx1], Lam.i = TLam[idx1, idx1], alpha = alpha,
          Vi = V[idx1, ], R.i = R[idx1], n.i = ni[i],
          mu.mo = mu.mo, Sig.mm.o = Sig.mm.o, na.idx.i = na.idx[[i]], Data.miss.i = Data.miss[idx1, ],
          V.fun.i = V.fun.i, Sig.mm.o.MC = Sig.mm.o.MC
        )
      }
      # print(m)
    }

    y.conv <- y.samp[-c(1:burn.in), ][cho, ]


    # ### 2
    b.hat <- matrix(0, nrow = N * q)
    sum.b2 <- matrix(0, N * q, N * q)
    sum.e2 <- matrix(0, n, n)
    sum.y2 <- 0
    sum.yb <- matrix(0, n, N * q)
    for (m in 1:mc.size)
    {
      Yhat.cent <- c(y.conv[m, ] - Xbeta)
      yy <- (y.conv[m, ] %*% t(y.conv[m, ]))
      sum.y2 <- sum.y2 + yy
      for (i in 1:N)
      {
        if (i == 1) {
          idx1 <- 1:cumsum.ni[1]
          idx2 <- 1:cumsum.q[1]
        } else {
          idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
          idx2 <- (cumsum.q[i - 1] + 1):cumsum.q[i]
        }
        b.hat[idx2] <- b.hat[idx2] + DD %*% t(Z[idx1, ]) %*% TLam.inv[idx1, idx1] %*% Yhat.cent[idx1]
        sum.b2[idx2, idx2] <- sum.b2[idx2, idx2] +
          (DD %*% t(Z[idx1, ]) %*% TLam.inv[idx1, idx1] %*% Yhat.cent[idx1] %*% t(DD %*% t(Z[idx1, ]) %*% TLam.inv[idx1, idx1] %*% (Yhat.cent[idx1])))
        sum.yb[idx1, idx2] <- sum.yb[idx1, idx2] + ((yy[idx1, idx1] - y.conv[m, ][idx1] %*% t(Xbeta[idx1, ])) %*% TLam.inv[idx1, idx1] %*% Z[idx1, ] %*% DD)
      }
    }
    y.hat <- colMeans(y.conv)
    b.hat <- b.hat / mc.size
    b2 <- sum.b2 / mc.size + TSig.b

    y2.hat <- sum.y2 / mc.size
    yb <- sum.yb / mc.size
    y_zbXbeta <- (y.hat - TZ %*% b.hat) %*% t(Xbeta)
    # E.hat1 = y2.hat - y_zbXbeta - t(y_zbXbeta) + Xbeta%*%t(Xbeta) - yb%*%t(TZ) - TZ%*%t(yb) + TZ%*%b2%*%t(TZ)
    E.hat1 <- matrix(0, n, n)
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
        idx2 <- 1:cumsum.q[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
        idx2 <- (cumsum.q[i - 1] + 1):cumsum.q[i]
      }
      E.hat1[idx1, idx1] <- y2.hat[idx1, idx1] - y_zbXbeta[idx1, idx1] - t(y_zbXbeta[idx1, idx1]) + Xbeta[idx1, ] %*% t(Xbeta[idx1, ]) -
        yb[idx1, idx2] %*% t(Z[idx1, ]) -
        Z[idx1, ] %*% t(yb[idx1, idx2]) + Z[idx1, ] %*% b2[idx2, idx2] %*% t(Z[idx1, ]) + Z[idx1, ] %*% TSig.b[idx2, idx2] %*% t(Z[idx1, ])
    }

    #####  CM-Step:
    ### Beta
    # Beta[,i] = solve(t(X) %*% TLam.inv[[i]] %*% X) %*% (t(X) %*% TLam.inv[[i]] %*% y)
    Beta <- solve(t(X) %*% TCor.inv %*% X) %*% (t(X) %*% TCor.inv %*% (y.hat - TZ %*% b.hat))

    ### D
    Nb2 <- 0
    for (i in 1:N) Nb2 <- Nb2 + b2[((i - 1) * q + 1):(i * q), ((i - 1) * q + 1):(i * q)]
    DD <- as.matrix(Nb2 / N)

    ### sigma
    Ce <- 0
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }
      Ce <- Ce + sum(TCor.inv[idx1, idx1] * E.hat1[idx1, idx1])
    }
    sigma <- Ce / sum(ni)

    ### Phi
    if (cor.type == "UNC") {
      Phi <- 1e-6
      ga <- 1
    }
    if (cor.type == "CAR1" | cor.type == "ARp" | cor.type == "CS") {
      ga <- 1
      Phi <- optim(par = Phi, fn = phiga.fn, method = "L-BFGS-B", lower = -1 + 1e-6, upper = 1 - 1e-6, sigma = sigma, E.hat1 = E.hat1, cumsum.ni = cumsum.ni, N = N, ni = ni, Data = Data)$par
    }
    if (cor.type == "DEC") {
      par.DEC <- optim(par = c(Phi, ga), fn = phiga.fn, method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1 - 1e-6, Inf), sigma = sigma, E.hat1 = E.hat1, cumsum.ni = cumsum.ni, N = N, ni = ni, Data = Data)$par
      Phi <- par.DEC[1]
      ga <- par.DEC[2]
    }
    if (cor.type == "BAND1") {
      ga <- 1
      Phi <- optim(par = Phi, fn = phiga.fn, method = "L-BFGS-B", lower = -1 / 2, upper = 1 / 2, sigma = sigma, E.hat1 = E.hat1, cumsum.ni = cumsum.ni, N = N, ni = ni, Data = Data)$par
    }


    ## alpha
    J.alpha <- S.alpha <- 0
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }
      treat.i <- Data$treat[idx1]
      R.i <- R[idx1]
      k3 <- k4 <- 0
      for (m in 1:mc.size)
      {
        a.i <- V.fun.i(y.conv[m, idx1], Data.miss[idx1, ])
        k3 <- k3 + t(a.i) %*% a.i
        pvi.i <- c(exp(a.i %*% alpha) / (1 + exp(a.i %*% alpha)))
        k4 <- k4 + t(a.i) %*% (R.i - pvi.i)
      }
      J.alpha <- J.alpha + k3 / mc.size
      S.alpha <- S.alpha + k4 / mc.size
    }
    alpha <- alpha + 4 * solve(J.alpha) %*% S.alpha

    # evaluate new log-likelihood
    TLam <- TCor <- TCor.inv <- matrix(0, ncol = n, nrow = n)
    for (i in 1:N) {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }
      Zi <- as.matrix(Z[idx1, ])
      Cor <- cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga)
      TLam[idx1, idx1] <- Zi %*% DD %*% t(Zi) + sigma * Cor
      TCor[idx1, idx1] <- Cor
      TCor.inv[idx1, idx1] <- solve(Cor)
    }
    # TLam.inv = solve(TLam)
    TLam.inv <- matrix(0, ncol = n, nrow = n)
    for (i in 1:N) {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }
      TLam.inv[idx1, idx1] <- solve(TLam[idx1, idx1])
    }
    TLam.oo <- TLam[-na.ind, -na.ind]
    TLam.mo <- TLam[na.ind, -na.ind]
    TLam.mm <- TLam[na.ind, na.ind]
    if (num.na == 1) TLam.mo <- t(TLam.mo)
    # TLam.oo.inv = solve(TLam.oo)
    TLam.oo.inv <- matrix(0, ncol = no, nrow = no)
    for (i in 1:N) {
      if (i == 1) {
        idx0 <- 1:cumsum.ni.o[1]
      } else {
        idx0 <- (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]
      }
      TLam.oo.inv[idx0, idx0] <- solve(TLam.oo[idx0, idx0])
    }


    #### log-likelihood 2
    V <- V.fun(y.hat, Data.miss)
    #### #### #### #### #### #### #### #### #### ####
    #### log-likelihood
    if (abs(diff.lnL) < tol || diff < tol || iter >= max.iter) {
      M.LL.i <- M.LL
    } else {
      M.LL.i <- M.LL / 100
    }
    wden <- numeric(N)
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }

      if (sum(i == y.na.ind) == 0) {
        pvi <- c(exp(as.matrix(V[idx1, ]) %*% alpha) / (1 + exp(as.matrix(V[idx1, ]) %*% alpha)))
        wden[i] <- mvtnorm::dmvnorm(y[idx1], X[idx1, ] %*% Beta, TLam[idx1, idx1], log = F) * prod(1 - pvi)
      } else {
        na.i <- which(y.na.ind == i)
        y.hat.ll <- mu.mo[na.i] + mvtnorm::rmvnorm(M.LL, sigma = as.matrix(Sig.mm.o[na.i, na.i]))
        dt.ym <- mvtnorm::dmvnorm(matrix(y.hat.ll), mean = mu.mo[na.i], sigma = as.matrix(Sig.mm.o[na.i, na.i]), log = F)
        wden.i <- NULL
        for (jj in 1:M.LL.i)
        {
          y.hat.i <- c(y[idx1][R[idx1] == 0], y.hat.ll[jj])
          V.i <- V.fun.i(y.i = y.hat.i, Data.i = Data.miss[idx1, ])
          pvi <- c(exp(V.i %*% alpha) / (1 + exp(V.i %*% alpha)))
          fn <- mvtnorm::dmvnorm(y.hat.i, X[idx1, ] %*% Beta, TLam[idx1, idx1], log = F) * prod((1 - pvi)^(1 - R[idx1]) * pvi^(R[idx1]))
          wden.i <- c(wden.i, fn / dt.ym[jj])
        }
        wden[i] <- mean(wden.i)
      }
    }
    loglik.new <- sum(log(wden))

    iter.lnL <- c(iter.lnL, loglik.new)
    diff.lnL <- (loglik.new - loglik.old)
    theta.new <- c(Beta, DD[vechD], sigma, Phi, ga, alpha)
    Tpara <- rbind(Tpara, theta.new)
    diff <- mean(abs((theta.new - theta.old) / theta.old)[-6])
    if (iter %% per == 0) cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, ",\t theta.diff = ", diff, sep = " ", "\n")
    if (abs(diff.lnL) < tol || diff < tol || iter >= max.iter) break
    loglik.old <- loglik.new
    theta.old <- theta.new
  }
  end <- proc.time()[1]
  # Parameter estimation
  cat(rep("=", 20), "Linear mixed models with ", cor.type[1], " errors", rep("=", 20), sep = "", "\n")
  cat("It took", end - begin, "seconds.\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, sep = "", "\n")
  cat("Beta =", Beta, "\n")
  cat("sigma =", sigma, "\n")
  cat("D =\n")
  print(DD)
  cat("Phi =", Phi, "\n")
  cat("ga =", ga, "\n")
  cat("alpha =", alpha, "\n")
  para.est <- list(Beta = Beta, sigma = sigma, D = DD, Phi = Phi, ga = ga, alpha = alpha, b = b.hat)

  # Model selection
  ma <- length(c(alpha))
  if (cor.type == "UNC") m <- g * (p + q * (q + 1) / 2 + 1) + ma
  if (cor.type == "ARp" | cor.type == "CS" | cor.type == "BAND1") m <- g * (p + q * (q + 1) / 2 + 1) + length(as.vector(Phi)) + ma
  if (cor.type == "DEC") m <- g * (p + q * (q + 1) / 2 + 1) + length(as.vector(Phi)) + length(ga) + ma
  aic <- 2 * m - 2 * loglik.new
  bic <- m * log(N) - 2 * loglik.new
  cat("aic =", aic, "\n")
  cat("bic =", bic, "\n")
  cat("MSE of missing response =", sum((y.hat - Data$y.c)^2) / sum(Data$R), "\n")
  cat(rep("=", 20), " Fisher information ", rep("=", 20), sep = "", "\n")
  begin <- proc.time()[1]
  IM <- I.lmm.missing(para.est,
    cor.type = cor.type, N = N, cumsum.ni = cumsum.ni, ni = ni, mc.size = mc.size, TLam.inv = TLam.inv,
    y.conv = y.conv, mechanism = mechanism, dropout.idx = dropout.idx, R = R, n = n, Data.miss = Data.miss
  )
  end <- proc.time()[1]
  cat("It took", end - begin, "seconds.\n")
  cat(rep("=", 50), sep = "", "\n")
  model.inf <- list(loglik = loglik.new, iter.lnL = iter.lnL, aic = aic, bic = bic, time = end - begin)
  return(list(
    model.inf = model.inf, para.est = para.est, iter = iter, y.c = y.hat, MSE.y = sum((y.hat - Data$y.c)^2) / sum(Data$R), MAE.y = sum(abs((y.hat - Data$y.c))) / sum(Data$R), MAPE.y = sum(abs((y.hat - Data$y.c) / Data$y.c)) / sum(Data$R),
    na.ind = na.ind, Tpara = Tpara, IM = IM
  ))
}

#' Observed-information (Fisher) components for the LME–dropout selection model
#'
#' Computes observed-information (or related sandwich/outer-product) matrices
#' for the fitted LME with selection-model dropout, using Monte Carlo averages
#' over completed-data draws obtained during estimation.
#'
#' @param para.est list. Parameter estimates from the fitted model:
#'   \code{Beta} (fixed effects), \code{D} (random-effects covariance),
#'   \code{sigma}, \code{Phi}, \code{ga} (if applicable), and \code{alpha}.
#' @param cor.type character(1). Working correlation structure (e.g., \code{"UNC"},
#'   \code{"CAR1"}, \code{"ARp"}, \code{"BAND1"}, \code{"CS"}, \code{"DEC"}).
#' @param N integer(1). Number of subjects.
#' @param cumsum.ni integer vector of length \eqn{N}. Cumulative sums of per-subject
#'   visit counts; used to slice block rows/columns.
#' @param ni integer vector of length \eqn{N}. Per-subject visit counts \eqn{n_i}.
#' @param mc.size integer(1). Monte Carlo sample size (number of completed-data draws)
#'   used for averaging the score/negative-Hessian components.
#' @param TLam.inv matrix. Block-diagonal inverse of the longitudinal covariance
#'   across subjects, aligned with the stacked data order.
#' @param y.conv matrix. Completed-response draws; typically of size
#'   \code{mc.size} \eqn{\times} \eqn{n} (stacked over all subjects).
#' @param mechanism character(1). Dropout mechanism controlling construction of
#'   the dropout-design \eqn{V} (\code{"MNAR"}, \code{"MAR"}, or \code{"MCAR"}).
#' @param dropout.idx integer vector. Indices of rows corresponding to dropout events
#'   used in internal diagnostics.
#' @param R numeric/integer vector of length \eqn{n}. Dropout/response indicators per row.
#' @param n integer(1). Total number of longitudinal rows across all subjects.
#' @param Data.miss data.frame. Long-format dataset used to rebuild the dropout
#'   design \eqn{V} (needs columns such as \code{Subject}, \code{Var1}, \code{week}, \code{prevOI}).
#'
#' @details
#' Let \eqn{p} be the number of fixed effects and \eqn{q} the random-effects dimension
#' (\eqn{q = ncol(D)}). The parameter vector stacks
#' \eqn{(\beta,\ \mathrm{vech}(D),\ \sigma,\ \Phi,\ \gamma,\ \alpha)}.
#' The function constructs per-parameter directional derivatives of subject-level
#' covariance blocks, accumulates Monte Carlo averages of score/Hessian terms,
#' and forms an observed-information approximation \eqn{I(\theta)} with optional
#' sandwich corrections. Internally, the function uses \code{vech.posi(q)} for indexing
#' the lower-triangular part of \eqn{D}; ensure that \eqn{q} is consistent with \code{para.est$D}.
#'
#' @return A list with components:
#' \describe{
#'   \item{out}{Matrix with two rows: the first row (\code{EST}) stacks the parameter
#'     estimates in the order \eqn{(\beta,\ \mathrm{vech}(D),\ \sigma,\ [\Phi,\ \gamma],\ \alpha)};
#'     the second row contains the corresponding standard errors. Column names
#'     are labeled by parameter blocks (e.g., \code{beta}, \code{d}, \code{sigma},
#'     \code{Phi}, \code{ga}, \code{alpha}) with appropriate repetition.}
#'   \item{se}{Numeric vector of standard errors (same order as \code{out} columns).}
#'   \item{I.theta}{Matrix. Estimated observed-information matrix for \eqn{\theta}.}
#'   \item{V.theta}{Matrix. Estimated variance–covariance matrix \eqn{I(\theta)^{-1}}
#'     (or its blockwise inverse variant depending on \code{cor.type}).}
#' }
I.lmm.missing <- function(para.est, cor.type = cor.type, N, cumsum.ni, ni, mc.size, TLam.inv, y.conv, mechanism, dropout.idx, R, n, Data.miss) {
  Beta <- para.est$Beta
  DD <- para.est$D
  sigma <- para.est$sigma
  Phi <- para.est$Phi
  ga <- para.est$ga
  alpha <- para.est$alpha
  vechD <- vech.posi(q)
  if (cor.type == "UNC") {
    ga <- NULL
    EST <- c(Beta, DD[vechD], sigma, alpha)
  }
  if (cor.type == "CAR1" | cor.type == "ARp" | cor.type == "CS" | cor.type == "BAND1") {
    ga <- 1
    EST <- c(Beta, DD[vechD], sigma, Phi, alpha)
  }
  if (cor.type == "DEC") {
    EST <- c(Beta, DD[vechD], sigma, Phi, ga, alpha)
  }

  ## A covance matrix for missing propabality
  if (mechanism == "MNAR") {
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
        # V.i = cbind(1, Data.i$treat[1:k], fData.i, bData.i)
        # V.i = cbind(1, Data.i$gender, Data.i$drug, Data.i$prevOI, Data.i$AZT, fData.i, bData.i)
        V.i <- cbind(1, Data.i$treat, fData.i, bData.i)
        V <- rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i <- function(y.i, Data.i) {
      k <- length(y.i)
      fData.i <- y.i[1:(k - 1)]
      fData.i <- y.i[-k]
      fData.i <- c(0, fData.i)
      bData.i <- y.i
      # V = cbind(rep(1, k), treat.i, fData.i, bData.i)
      # V = cbind(rep(1, k), Data.i$gender, Data.i$drug, Data.i$prevOI, Data.i$AZT, fData.i, bData.i)
      V <- cbind(rep(1, k), Data.i$treat, fData.i, bData.i)
      return(V)
    }
    # print(mechanism)
  }
  if (mechanism == "MAR") {
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
        V.i <- cbind(1, Data.i$treat, fData.i)
        V <- rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i <- function(y.i, Data.i) {
      k <- length(y.i)
      fData.i <- y.i[1:(k - 1)]
      fData.i <- y.i[-k]
      fData.i <- c(0, fData.i)
      bData.i <- y.i
      V <- cbind(rep(1, k), Data.i$treat, fData.i)
      return(V)
    }
  }

  if (mechanism == "MCAR") {
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
        V.i <- cbind(rep(1, length(idx1)))
        V <- rbind(V, V.i)
      }
      return(V)
    }
    V.fun.i <- function(y.i, Data.i) {
      k <- length(y.i)
      V <- cbind(rep(1, k))
      return(V)
    }
  }



  # Information matrix
  g1 <- q * (q + 1) / 2
  g2 <- r <- 1
  g3 <- length(c(Phi)) ## phi
  if (cor.type == "DEC") g3 <- 2
  g <- g1 + g2 + g3
  ma <- length(alpha)
  dot.L <- as.list(numeric(g))
  dot.L <- as.list(matrix(0, N, g))
  ### 0 ~ 300 dot D
  for (l in 1:g1)
  {
    dot.DD <- matrix(0, q, q)
    dot.DD[matrix(vechD[l, ], 1)] <- dot.DD[matrix(rev(vechD[l, ]), 1)] <- 1
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }
      dot.L[[(l - 1) * N + i]] <- Z[idx1, ] %*% dot.DD %*% t(Z[idx1, ])
    }
  }
  ### 301 ~ 400 dot sigma
  for (l in 1:g2)
  {
    dot.Sig <- matrix(0, r, r)
    dot.Sig <- 1
    for (i in 1:N) dot.L[[g1 * N + (l - 1) * N + i]] <- kronecker(dot.Sig, cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga))
  }

  dot.Sig <- matrix(0, 1, 1)
  dot.Sig <- 1
  for (i in 1:N) dot.L[[g1 * N + (l - 1) * N + i]] <- kronecker(dot.Sig, cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga))

  ### 400~ 500 dot phi
  if (cor.type == "UNC") {
    for (i in 1:N) dot.L[[(g1 + g2) * N + i]] <- diag(0, ni[i])
  }
  if (cor.type == "CAR1" | cor.type == "DEC") {
    for (i in 1:N) dot.L[[(g1 + g2) * N + i]] <- sigma * DEC.dot.phi(phi = Phi, ga = ga, Ti = Data$week[Data$Subject == i])
  }
  if (cor.type == "ARp") {
    # for(i in 1: N) dot.L[[(g1+g2)*N+i]] = arp.Ci.dot(Phi, dim=ni[i], l)
    for (i in 1:N) dot.L[[(g1 + g2) * N + i]] <- sigma * Arp.Ci.dot(Phi, dim = ni[i])
  }
  if (cor.type == "CS") {
    for (i in 1:N) dot.L[[(g1 + g2) * N + i]] <- sigma * cs.Ci.dot(Phi, dim = ni[i])
  }
  if (cor.type == "BAND1") {
    for (i in 1:N) dot.L[[(g1 + g2) * N + i]] <- sigma * BAND1.Ci.dot(Phi, dim = ni[i])
  }

  ### 500~ 600 dot ga
  if (cor.type == "DEC") for (i in 1:N) dot.L[[(g1 + g2 + 1) * N + i]] <- sigma * DEC.dot.ga(phi = Phi, ga = ga, Ti = Data$week[Data$Subject == i])

  H <- SS <- matrix(0, ncol = (p + g + ma), nrow = (p + g + ma))
  S.sum <- numeric((p + g + ma))
  H.alpha <- SS.alpha <- matrix(0, ma, ma)
  S.alpha <- numeric(ma)
  for (m in 1:mc.size)
  {
    if (m %% 100 == 0) cat(rep("=", 10), "Fisher information ", m / mc.size * 100, " % ", rep("=", 10), sep = "", "\n")
    S <- numeric((p + g + ma))
    # H.beta
    H[1:p, 1:p] <- H[1:p, 1:p] - t(X) %*% TLam.inv %*% X
    ys.cent <- y.conv[m, ] - X %*% Beta
    # S.beta
    sb <- t(X) %*% TLam.inv %*% ys.cent
    S[1:p] <- S[1:p] + sb

    # H.xi
    Linv.dotL <- as.list(numeric(N * g))
    for (s in 1:g) {
      for (i in 1:N) {
        if (i == 1) {
          idx1 <- 1:cumsum.ni[1]
        } else {
          idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
        }
        Linv.dotL[[(s - 1) * N + i]] <- TLam.inv[idx1, idx1] %*% dot.L[[(s - 1) * N + i]]
        H[1:p, (p + s)] <- H[1:p, (p + s)] - t(X[idx1, ]) %*% Linv.dotL[[(s - 1) * N + i]] %*% TLam.inv[idx1, idx1] %*% ys.cent[idx1]
        for (l in 1:s) {
          H[p + s, p + l] <- H[p + s, p + l] + 0.5 * sum(diag(Linv.dotL[[(s - 1) * N + i]] %*% Linv.dotL[[(l - 1) * N + i]])) -
            0.5 * sum(diag((ys.cent[idx1]) %*% t(ys.cent[idx1]) %*% (Linv.dotL[[(s - 1) * N + i]] %*% Linv.dotL[[(l - 1) * N + i]] %*% TLam.inv[idx1, idx1] + Linv.dotL[[(l - 1) * N + i]] %*% Linv.dotL[[(s - 1) * N + i]] %*% TLam.inv[idx1, idx1])))
        }
      }
    }
    for (l in (p + 1):(p + g - 1)) for (s in (l + 1):(p + g)) H[l, s] <- H[s, l]
    H[(p + 1):(p + g), 1:p] <- t(H[1:p, (p + 1):(p + g)])
    for (s in 1:g) {
      for (i in 1:N) {
        if (i == 1) {
          idx1 <- 1:cumsum.ni[1]
        } else {
          idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
        }
        Linv.dotL[[(s - 1) * N + i]] <- TLam.inv[idx1, idx1] %*% dot.L[[(s - 1) * N + i]]
        S[p + s] <- S[p + s] - 0.5 * sum(diag(Linv.dotL[[(s - 1) * N + i]])) + 0.5 * sum(diag(ys.cent[idx1] %*% t(ys.cent[idx1]) %*% Linv.dotL[[(s - 1) * N + i]] %*% TLam.inv[idx1, idx1]))
      }
    }
    V <- V.fun(y.conv[m, ], Data.miss)
    pvi <- c(exp(V %*% alpha) / (1 + exp(V %*% alpha)))
    ppvi <- pvi * (1 - pvi)
    S.alpha <- t(V) %*% (R - pvi)
    for (i in 1:n) H.alpha <- H.alpha + t(t(V[i, ])) %*% ppvi[i] %*% t(V[i, ])
    S[(p + g + 1):(p + g + ma)] <- S.alpha
    S.sum <- S.sum + S
    SS <- SS + S %*% t(S)
    # print(m)
  }
  H.mean <- H / mc.size
  H.alpha <- -H.alpha / mc.size
  SS.mean <- SS / mc.size
  S.mean <- S.sum / mc.size

  H.mean[(p + g + 1):(p + g + ma), (p + g + 1):(p + g + ma)] <- H.alpha

  I.theta <- -H.mean - SS.mean + S.mean %*% t(S.mean)
  if (cor.type == "UNC") {
    # V.theta = solve(I.theta[1:(p+g-1),1:(p+g-1)])
    V.theta <- solve(I.theta[-(p + g), -(p + g)])
    sd.theta <- c(sqrt(diag(V.theta)))
  }
  if (cor.type == "CAR1" | cor.type == "ARp" | cor.type == "CS" | cor.type == "BAND1") {
    V.theta <- solve(I.theta)
    # V.theta = solve(I.theta[1:(p+g-1),1:(p+g-1)])
    sd.theta <- c(sqrt(diag(V.theta)))
  }
  if (cor.type == "DEC") {
    V.theta <- solve(I.theta)
    # V.theta = solve(I.theta[1:9,1:9])
    sd.theta <- sqrt(diag(V.theta))
  }
  out <- rbind(EST, c(sd.theta))
  if (cor.type == "UNC") colnames(out) <- rep(c("beta", "d", "sigma", "alpha"), c(p, length(DD[vechD]), 1, ma))
  if (cor.type == "ARp" | cor.type == "CS" | cor.type == "BAND1") colnames(out) <- rep(c("beta", "d", "sigma", "Phi", "alpha"), c(p, length(DD[vechD]), 1, 1, ma))
  if (cor.type == "DEC") colnames(out) <- rep(c("beta", "d", "sigma", "Phi", "ga", "alpha"), c(p, length(DD[vechD]), 1, 1, 1, ma))
  SD <- list(out = out, se = c(sd.theta), I.theta = I.theta, V.theta = V.theta)
  # }
  return(SD)
}
