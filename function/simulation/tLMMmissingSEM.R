################################################################################
#                                                                                     
#   Filename    :    tLMMmissingSEM.R  												  
#   Project     :    BiomJ article "Extending t linear mixed models for longitudinal 
#                    data with non-ignorable dropout applied to AIDS studies"                                                           
#   Authors     :    Yu-Chen Yang and Wan-Lun Wang and Luis M. Castro and Tsung-I Lin
#   Date        :    07.08.2025
#   Purpose     :    Function allowing to fit a tLME model with non-ignorable dropout 
#                    under the MCAR, MAR, and MNAR mechanisms  
#
#   Input data files : A standard data frame (named 'Data') containing subject ID, 
#                      time points, observed responses, missingness indicators, 
#                      and additional covariates (if applicable). 
#   Output data files : an object colleting the fitting results of 
#                       t linear mixed-effects (tLME) model with  
#                       selection modeling-based MCAR, MAR, and MNAR mechanisms 
#
#   R Version   :    R-4.3.1                                                              
#   Required R packages : nlme; mvtnorm  
#
################################################################################ 
tLMM.miss.SEM <- function(Data, X, Z, V, g = g, init.para, cor.type = c("UNC", "ARp", "BAND1", "CS"), M = 100, M.LL = 1000, P = 1, tol = 1e-6, max.iter = max.iter, per = 1, mechanism = c("MNAR", "MAR", "MCAR")) {
  begin <- proc.time()[1]
  # initial values of parameter
  Beta <- init.para$Beta
  DD <- init.para$DD
  sigma <- init.para$sigma
  alpha <- init.para$alpha
  Phi <- init.para$Phi
  nu <- init.para$nu
  ga <- 1
  p <- ncol(X)
  q <- ncol(Z)
  N <- length(unique(Data$Subject))
  na.ind <- which(is.na(as.vector(t(Data$Var1)))) ## which time point have missing value
  Data.miss <- Data

  y <- Data.miss$Var1
  R <- Data.miss$R

  yo <- y[-na.ind]
  ym <- Data$yc[na.ind]
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
  na.idx <- as.list(N)
  for (i in 1:N) na.idx[[i]] <- NA
  na.idx[[y.na.ind[[1]]]] <- 1:cumsum.na[1]
  for (i in 2:Nm) na.idx[[y.na.ind[i]]] <- (cumsum.na[i - 1] + 1):cumsum.na[i]


  if (mechanism == "MNAR") cat(rep("=", 25), "tLMM (MNAR) with ", cor.type, " errors is fitted...; ", "missing = ", Nm / N * 100, "%", rep("=", 25), sep = "", "\n")
  if (mechanism == "MAR") cat(rep("=", 25), "tLMM (MAR) with ", cor.type, " errors is fitted...; ", "missing = ", Nm / N * 100, "%", rep("=", 25), sep = "", "\n")
  if (mechanism == "MCAR") cat(rep("=", 25), "tLMM (MCAR) with ", cor.type, " errors is fitted...; ", "missing = ", Nm / N * 100, "%", rep("=", 25), sep = "", "\n")

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
  Sig.mm.o.MC <- (TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo))

  # observed log-likelihood:
  Xbeta <- X %*% Beta
  y.samp <- matrix(rep(y, M), nrow = M, ncol = n, byrow = T)
  y.hat <- t(TO) %*% yo + t(TM) %*% ym
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

  yo.cent <- yo - Xo %*% Beta
  Delta <- t(yo.cent[1:cumsum.ni.o[1], ]) %*% TLam.oo.inv[1:cumsum.ni.o[1], 1:cumsum.ni.o[1]] %*% yo.cent[1:cumsum.ni.o[1], ]
  for (i in 2:N) Delta <- c(Delta, t(yo.cent[(cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i], ]) %*% TLam.oo.inv[(cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i], (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]] %*% yo.cent[(cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i], ])
  tau <- (nu + ni.o) / (nu + Delta)

  #### #### #### #### #### #### #### #### #### ####
  #### log-likelihood
  mu.mo <- Xm %*% Beta + TLam.mo %*% TLam.oo.inv %*% yo.cent
  Sig.mm.o <- (TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo)) / rep(tau, time = ni)[na.ind]
  wden <- numeric(N)
  for (i in 1:N)
  {
    if (i == 1) {
      idx0 <- 1:cumsum.ni.o[1]
      idx1 <- 1:cumsum.ni[1]
    } else {
      idx0 <- (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]
      idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
    }
    if (sum(i == y.na.ind) == 0) {
      pvi <- c(exp(as.matrix(V[idx1, ]) %*% alpha) / (1 + exp(as.matrix(V[idx1, ]) %*% alpha)))
      wden[i] <- mvtnorm::dmvt(y[idx1], X[idx1, ] %*% Beta, TLam[idx1, idx1], df = nu, log = F) * prod(1 - pvi)
    } else {
      na.i <- which(y.na.ind == i)
      y.hat.ll <- mu.mo[na.i, ] + mvtnorm::rmvt(M.LL, sigma = as.matrix(Sig.mm.o[na.i, na.i]), df = nu + ni.o[y.na.ind][na.i])
      dt.ym <- mvtnorm::dmvt(matrix(y.hat.ll), delta = mu.mo[na.i], sigma = as.matrix(Sig.mm.o[na.i, na.i]), df = nu + ni.o[y.na.ind][na.i], log = F)
      wden.i <- NULL
      for (jj in 1:(M.LL / 100))
      {
        y.hat.i <- c(y[idx1][R[idx1] == 0], y.hat.ll[jj])
        V.i <- V.fun.i(y.i = y.hat.i, Data.miss[idx1, ])
        pvi <- c(exp(V.i %*% alpha) / (1 + exp(V.i %*% alpha)))
        fn <- mvtnorm::dmvt(y.hat.i, X[idx1, ] %*% Beta, TLam[idx1, idx1], df = nu, log = F) * prod((1 - pvi)^(1 - R[idx1]) * pvi^(R[idx1]))
        wden.i <- c(wden.i, fn / dt.ym[jj])
      }
      wden[i] <- mean(wden.i)
    }
  }
  loglik.old <- iter.lnL <- sum(log(wden))

  theta.old <- c(Beta, DD[vechD], sigma, Phi, nu, alpha)
  iter <- 0
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  cat("t Linear mixed models with ", cor.type[1], " errors: ", "\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.old, sep = "", "\n")
  Tpara <- theta.old
  diff.lnL <- 10000
  diff <- 10000
  repeat  {
    iter <- iter + 1
    ##### ##### ##### ##### ##### ##### ##### #####
    ######  E-Step:
    yo.cent <- yo - Xo %*% Beta
    Delta <- t(yo.cent[1:cumsum.ni.o[1], ]) %*% TLam.oo.inv[1:cumsum.ni.o[1], 1:cumsum.ni.o[1]] %*% yo.cent[1:cumsum.ni.o[1], ]
    for (i in 2:N) Delta <- c(Delta, t(yo.cent[(cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i], ]) %*% TLam.oo.inv[(cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i], (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]] %*% yo.cent[(cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i], ])
    tau <- (nu + ni.o) / (nu + Delta)
    kappa <- digamma((nu + ni.o) / 2) - log((nu + Delta) / 2)

    TD <- kronecker(diag(N), DD)
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
      TSig.b[idx2, idx2] <- solve(t(TZ[idx1, idx2]) %*% TCor.inv[idx1, idx1] %*% TZ[idx1, idx2] / sigma + solve(DD)) / tau[i]
    }

    ##### ##### ##### ##### ##### ##### ##### #####
    ##### genrate ym by MCMC
    mu.mo <- Xm %*% Beta + TLam.mo %*% TLam.oo.inv %*% yo.cent
    Sig.mm.o <- (TLam.mm - TLam.mo %*% TLam.oo.inv %*% t(TLam.mo)) / rep(tau, time = ni)[na.ind]
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
        y.samp[m, idx1] <- MH.y.miss.t(
          yi = y.samp[(m - 1), idx1], Xbeta.i = c(Xbeta)[idx1], Lam.i = TLam[idx1, idx1], tau.i = tau[i], alpha = alpha,
          Vi = as.matrix(V[idx1, ]), R.i = R[idx1], n.i = ni[i], mu.mo = mu.mo, Sig.mm.o = Sig.mm.o,
          Sig.mm.o.MC = Sig.mm.o.MC, na.idx.i = na.idx[[i]], Data.miss.i = Data.miss[idx1, ],
          V.fun.i = V.fun.i
        )
      }
      # print(m)
    }

    y.conv <- y.samp[-c(1:burn.in), ][cho, ]

    ####
    Xbeta <- X %*% Beta
    b.hat <- matrix(0, nrow = N * q)
    sum.b2 <- matrix(0, N * q, N * q)
    sum.Omega2 <- matrix(0, n, n)

    for (m in 1:mc.size)
    {
      Yhat.cent <- c(y.conv[m, ] - Xbeta)
      for (i in 1:N)
      {
        if (i == 1) {
          idx1 <- 1:cumsum.ni[1]
          idx2 <- 1:cumsum.q[1]
        } else {
          idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
          idx2 <- (cumsum.q[i - 1] + 1):cumsum.q[i]
        }
        b.hat.i <- DD %*% t(Z[idx1, ]) %*% TLam.inv[idx1, idx1] %*% Yhat.cent[idx1]
        b.hat[idx2] <- b.hat[idx2] + b.hat.i
        sum.b2[idx2, idx2] <- sum.b2[idx2, idx2] +
          (DD %*% t(Z[idx1, ]) %*% TLam.inv[idx1, idx1] %*% Yhat.cent[idx1] %*% t(DD %*% t(Z[idx1, ]) %*% TLam.inv[idx1, idx1] %*% (Yhat.cent[idx1])))
        # sum.yb[idx1,idx2] = sum.yb[idx1,idx2] + ((yy[idx1,idx1] - y.conv[m, ][idx1]%*%t(Xbeta[idx1,])) %*% TLam.inv[idx1,idx1] %*% Z[idx1,] %*% DD)
        sum.Omega2[idx1, idx1] <- sum.Omega2[idx1, idx1] + (Yhat.cent[idx1] - Z[idx1, ] %*% b.hat.i) %*% t(Yhat.cent[idx1] - Z[idx1, ] %*% b.hat.i) + Z[idx1, ] %*% TSig.b[idx2, idx2] %*% t(Z[idx1, ])
      }
    }
    y.hat <- colMeans(y.conv)
    b.hat <- b.hat / mc.size
    b2 <- sum.b2 / mc.size + TSig.b
    E.hat <- sum.Omega2 / mc.size

    ##### ##### ##### ##### ##### ##### ##### #####
    #####  CM-Step:
    ### Beta
    # Beta = solve(t(rep(tau, times = ni) * X) %*% TCor.inv %*% X) %*% (t(rep(tau, times = ni) * X) %*% TCor.inv %*% (y.hat - TZ %*% b.hat))
    k1 <- k2 <- 0
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
        idx2 <- 1:cumsum.q[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
        idx2 <- (cumsum.q[i - 1] + 1):cumsum.q[i]
      }
      k1 <- k1 + t(tau[i] * X[idx1, ]) %*% TCor.inv[idx1, idx1] %*% X[idx1, ]
      k2 <- k2 + t(tau[i] * X[idx1, ]) %*% TCor.inv[idx1, idx1] %*% (y.hat[idx1] - Z[idx1, ] %*% b.hat[idx2])
    }
    Beta <- solve(k1) %*% k2


    ### D
    Nb2 <- 0
    for (i in 1:N) Nb2 <- Nb2 + tau[i] * b2[((i - 1) * q + 1):(i * q), ((i - 1) * q + 1):(i * q)]
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
      Ce <- Ce + sum(tau[i] * TCor.inv[idx1, idx1] * E.hat[idx1, idx1])
    }
    sigma <- Ce / sum(ni)

    ### Phi
    if (cor.type == "UNC") {
      Phi <- 1e-6
      ga <- 1
    }
    if (cor.type == "CAR1" | cor.type == "ARp" | cor.type == "CS") {
      ga <- 1
      Phi <- optim(
        par = Phi, fn = phiga.t.fn, method = "L-BFGS-B", lower = -1 + 1e-6, upper = 1 - 1e-6,
        sigma = sigma, E.hat1 = E.hat, cumsum.ni = cumsum.ni, tau = tau, N = N, ni = ni, Data = Data, cor.type = cor.type[1]
      )$par
    }
    if (cor.type == "DEC") {
      par.DEC <- optim(
        par = c(Phi, ga), fn = phiga.t.fn, method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(1 - 1e-6, Inf),
        sigma = sigma, E.hat1 = E.hat, cumsum.ni = cumsum.ni, tau = tau, N = N, ni = ni, Data = Data, cor.type = cor.type[1]
      )$par
      Phi <- par.DEC[1]
      ga <- par.DEC[2]
    }
    if (cor.type == "BAND1") {
      ga <- 1
      Phi <- optim(
        par = Phi, fn = phiga.t.fn, method = "L-BFGS-B", lower = -1 / 2, upper = 1 / 2,
        sigma = sigma, E.hat1 = E.hat, cumsum.ni = cumsum.ni, tau = tau, N = N, ni = ni, Data = Data, cor.type = cor.type[1]
      )$par
    }


    ### nu
    nu <- optim(par = nu, fn = nu.Q.fn, method = "L-BFGS-B", lower = 2, upper = 200, kappa = kappa, tau = tau, N = N)$par

    ### alpha
    J.alpha <- S.alpha <- 0
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }
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
        idx0 <- 1:cumsum.ni.o[1]
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx0 <- (cumsum.ni.o[i - 1] + 1):cumsum.ni.o[i]
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }
      if (sum(i == y.na.ind) == 0) {
        pvi <- c(exp(as.matrix(V[idx1, ]) %*% alpha) / (1 + exp(as.matrix(V[idx1, ]) %*% alpha)))
        wden[i] <- mvtnorm::dmvt(y[idx1], X[idx1, ] %*% Beta, TLam[idx1, idx1], df = nu, log = F) * prod(1 - pvi)
      } else {
        na.i <- which(y.na.ind == i)
        y.hat.ll <- mu.mo[na.i, ] + mvtnorm::rmvt(M.LL, sigma = as.matrix(Sig.mm.o[na.i, na.i]), df = nu + ni.o[y.na.ind][na.i])
        dt.ym <- mvtnorm::dmvt(matrix(y.hat.ll), delta = mu.mo[na.i], sigma = as.matrix(Sig.mm.o[na.i, na.i]), df = nu + ni.o[y.na.ind][na.i], log = F)
        wden.i <- NULL
        for (jj in 1:M.LL.i)
        {
          y.hat.i <- c(y[idx1][R[idx1] == 0], y.hat.ll[jj])
          # mvtnorm::dmvt(y.hat.i, X[idx1, ]%*%Beta, TLam[idx1,idx1], df = nu, log=F) * prod(1-pvi)
          V.i <- V.fun.i(y.i = y.hat.i, Data.miss[idx1, ])
          pvi <- c(exp(V.i %*% alpha) / (1 + exp(V.i %*% alpha)))
          fn <- mvtnorm::dmvt(y.hat.i, X[idx1, ] %*% Beta, TLam[idx1, idx1], df = nu, log = F) * prod((1 - pvi)^(1 - R[idx1]) * pvi^(R[idx1]))
          wden.i <- c(wden.i, fn / dt.ym[jj])
        }
        wden[i] <- mean(wden.i)
      }
    }
    loglik.new <- sum(log(wden))

    iter.lnL <- c(iter.lnL, loglik.new)
    diff.lnL <- (loglik.new - loglik.old)
    theta.new <- c(Beta, DD[vechD], sigma, Phi, nu, alpha)
    Tpara <- rbind(Tpara, theta.new)
    diff <- mean(abs((theta.new - theta.old) / theta.old)[-6])
    if (iter %% per == 0) cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, ",\t theta.diff = ", diff, sep = " ", "\n")
    if (abs(diff.lnL) < tol || diff < tol || iter >= max.iter) break
    loglik.old <- loglik.new
    theta.old <- theta.new
  }
  end <- proc.time()[1]
  # Parameter estimation
  cat(rep("=", 20), "t Linear mixed models with ", cor.type[1], " errors", rep("=", 20), sep = "", "\n")
  cat("It took", end - begin, "seconds.\n")
  cat("iter = ", iter, ",\t obs.loglik = ", loglik.new, sep = "", "\n")
  cat("Beta =", Beta, "\n")
  cat("sigma =", sigma, "\n")
  cat("D =\n")
  print(DD)
  cat("Phi =", Phi, "\n")
  cat("nu =", nu, "\n")
  cat("alpha =", alpha, "\n")
  para.est <- list(Beta = Beta, sigma = sigma, D = DD, Phi = Phi, nu = nu, alpha = alpha, b = b.hat)
  cat(rep("=", 20), " Fisher information ", rep("=", 20), sep = "", "\n")
  begin <- proc.time()[1]
  IM <- I.tlmm.missing(para.est,
    cor.type = cor.type, X = X, N = N, cumsum.ni = cumsum.ni, ni = ni, mc.size = mc.size, TLam.inv = TLam.inv,
    y.conv = y.conv, mechanism = mechanism, dropout.idx = dropout.idx, R = R, n = n, q = q, Data.miss = Data.miss
  )
  end <- proc.time()[1]
  cat("It took", end - begin, "seconds.\n")
  cat(rep("=", 50), sep = "", "\n")
  ma <- length(c(alpha))
  if (cor.type == "UNC") {
    m <- g * (p + 1 + q * (q + 1) / 2 + 1) + ma
  } else {
    m <- g * (p + 1 + q * (q + 1) / 2 + 1) + length(as.vector(Phi)) + ma
  }
  aic <- 2 * m - 2 * loglik.new
  bic <- m * log(N) - 2 * loglik.new
  cat("aic =", aic, "\n")
  cat("bic =", bic, "\n")
  cat("MSE of missing response =", sum((y.hat - Data$y.c)^2) / sum(Data$R), "\n")
  cat(paste(rep("=", 50), sep = "", collapse = ""), "\n")
  model.inf <- list(loglik = loglik.new, iter.lnL = iter.lnL, aic = aic, bic = bic, time = end - begin)
  return(list(
    model.inf = model.inf, para.est = para.est, iter = iter, y.c = y.hat, IM = IM,
    MSE.y = sum((y.hat - Data$y.c)^2) / sum(Data$R), MAE.y = sum(abs((y.hat - Data$y.c))) / sum(Data$R), MAPE.y = sum(abs((y.hat - Data$y.c) / Data$y.c)) / sum(Data$R),
    na.ind = na.ind, Tpara = Tpara
  ))
}


I.tlmm.missing <- function(para.est, cor.type = cor.type, X, N, cumsum.ni, ni, mc.size, TLam.inv, y.conv, mechanism, dropout.idx, R, n, q, Data.miss) {
  Beta <- para.est$Beta
  DD <- para.est$D
  sigma <- para.est$sigma
  Phi <- para.est$Phi
  nu <- para.est$nu
  alpha <- para.est$alpha
  vechD <- vech.posi(q)
  if (cor.type == "UNC") {
    EST <- c(Beta, DD[vechD], sigma, nu, alpha)
  } else {
    EST <- c(Beta, DD[vechD], sigma, Phi, nu, alpha)
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
  g1 <- q * (q + 1) / 2 ## D
  g2 <- r <- 1 ## sigma
  g3 <- length(Phi) ## phi
  g <- g1 + g2 + g3
  ma <- length(alpha)

  # dot alpha
  Tdot.L <- array(0, dim = c(n, n, g))
  for (l in 1:g1)
  {
    TZDZ <- matrix(0, ncol = n, nrow = n)
    dot.DD <- matrix(0, q, q)
    dot.DD[matrix(vechD[l, ], 1)] <- dot.DD[matrix(rev(vechD[l, ]), 1)] <- 1
    TZDZ[1:cumsum.ni[1], 1:cumsum.ni[1]] <- Z[1:cumsum.ni[1], ] %*% dot.DD %*% t(Z[1:cumsum.ni[1], ])
    for (i in 2:N)
    {
      rZj <- Z[(cumsum.ni[i - 1] + 1):cumsum.ni[i], ]
      TZDZ[(cumsum.ni[i - 1] + 1):cumsum.ni[i], (cumsum.ni[i - 1] + 1):cumsum.ni[i]] <- rZj %*% dot.DD %*% t(rZj)
    }
    Tdot.L[, , l] <- TZDZ
  }

  for (l in 1:g2)
  {
    Tdot.L[, , g1 + l][1:cumsum.ni[1], 1:cumsum.ni[1]] <- cor.fn(Phi, dim = ni[1], type = cor.type[1], Ti = Data$week[Data$Subject == 1], ga = 1)
    for (i in 2:N) Tdot.L[, , g1 + l][(cumsum.ni[i - 1] + 1):cumsum.ni[i], (cumsum.ni[i - 1] + 1):cumsum.ni[i]] <- cor.fn(Phi, dim = ni[i], type = cor.type[1], Ti = Data$week[Data$Subject == i], ga = ga)
  }

  if (cor.type == "UNC") {
    Tdot.L[, , (g1 + g2 + 1)] <- diag(0, sum(ni))
  }
  if (cor.type == "CAR1") {
    # for(i in 1: N) dot.L[[(g1+g2)*N+i]] = sigma*DEC.dot.phi(phi=Phi, ga=ga, Ti=Data$week[Data$Subject == i])
    for (l in 1:g3) {
      Tdot.L[, , (g1 + g2 + l)][1:cumsum.ni[1], 1:cumsum.ni[1]] <- sigma * DEC.dot.phi(phi = Phi, ga = 1, Ti = Data$week[Data$Subject == 1])
      for (i in 2:N) Tdot.L[, , (g1 + g2 + l)][(cumsum.ni[i - 1] + 1):cumsum.ni[i], (cumsum.ni[i - 1] + 1):cumsum.ni[i]] <- sigma * DEC.dot.phi(phi = Phi, ga = 1, Ti = Data$week[Data$Subject == i])
    }
  }
  if (cor.type == "ARp") {
    # for(i in 1: N) dot.L[[(g1+g2)*N+i]] = sigma*Arp.Ci.dot(Phi, dim=ni[i])
    for (l in 1:g3) {
      Tdot.L[, , (g1 + g2 + l)][1:cumsum.ni[1], 1:cumsum.ni[1]] <- sigma * Arp.Ci.dot(Phi, dim = ni[1])
      for (i in 2:N) Tdot.L[, , (g1 + g2 + l)][(cumsum.ni[i - 1] + 1):cumsum.ni[i], (cumsum.ni[i - 1] + 1):cumsum.ni[i]] <- sigma * Arp.Ci.dot(Phi, dim = ni[i])
    }
  }
  if (cor.type == "BAND1") {
    # for(i in 1: N) dot.L[[(g1+g2)*N+i]] = sigma*BAND1.Ci.dot(Phi, dim=ni[i])
    for (l in 1:g3) {
      Tdot.L[, , (g1 + g2 + l)][1:cumsum.ni[1], 1:cumsum.ni[1]] <- sigma * BAND1.Ci.dot(Phi, dim = ni[1])
      for (i in 2:N) Tdot.L[, , (g1 + g2 + l)][(cumsum.ni[i - 1] + 1):cumsum.ni[i], (cumsum.ni[i - 1] + 1):cumsum.ni[i]] <- sigma * BAND1.Ci.dot(Phi, dim = ni[i])
    }
  }



  H <- matrix(0, ncol = (p + g + 1), nrow = (p + g + 1))
  SS <- matrix(0, ncol = (p + g + 1 + ma), nrow = (p + g + 1 + ma))
  S <- numeric((p + g + 1 + ma))
  H.alpha <- matrix(0, ma, ma)
  H.gamma <- 0
  for (m in 1:mc.size)
  {
    # print(mc.size)
    if (m %% 100 == 0) cat(rep("=", 10), "Fisher information ", m / mc.size * 100, " % ", rep("=", 10), sep = "", "\n")
    ys.cent <- y.conv[m, ] - X %*% Beta
    Delta <- t(ys.cent[1:cumsum.ni[1], ]) %*% TLam.inv[1:cumsum.ni[1], 1:cumsum.ni[1]] %*% ys.cent[1:cumsum.ni[1], ]
    for (i in 2:N) Delta <- c(Delta, t(ys.cent[(cumsum.ni[i - 1] + 1):cumsum.ni[i], ]) %*% TLam.inv[(cumsum.ni[i - 1] + 1):cumsum.ni[i], (cumsum.ni[i - 1] + 1):cumsum.ni[i]] %*% ys.cent[(cumsum.ni[i - 1] + 1):cumsum.ni[i], ])
    tau <- (nu + ni) / (nu + Delta)
    ka1 <- 1 / (nu + Delta)
    ka2 <- (nu + ni) / (nu + Delta)^2
    ### Hesian matrix
    # H beta beta
    # H[1:p, 1:p] = H[1:p, 1:p] - t(rep(tau, ni) * X) %*% TLam.inv %*% X +
    #   t(rep(2*tau*ka1, ni) * X) %*% TLam.inv %*% ys.cent %*% t(ys.cent) %*% TLam.inv %*% X
    # H beta beta
    for (i in 1:N)
    {
      if (i == 1) {
        idx1 <- 1:cumsum.ni[1]
      } else {
        idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
      }
      H[1:p, 1:p] <- H[1:p, 1:p] - t(rep(tau[i], ni[i]) * X[idx1, ]) %*% TLam.inv[idx1, idx1] %*% X[idx1, ] +
        t(rep(2 * tau[i] * ka1[i], ni[i]) * X[idx1, ]) %*% TLam.inv[idx1, idx1] %*% ys.cent[idx1, ] %*% t(ys.cent[idx1, ]) %*% TLam.inv[idx1, idx1] %*% X[idx1, ]
    }
    # H beta nu
    H[1:p, p + g + 1] <- H[1:p, p + g + 1] + t(rep(ka1 - ka2, ni) * X) %*% TLam.inv %*% ys.cent
    # H nu nu
    H[p + g + 1, p + g + 1] <- H[p + g + 1, p + g + 1] + (1 / 4) *
      sum((trigamma((nu + ni) / 2) - trigamma(nu / 2) + 2 * ni / nu^2 + 4 * Delta * ka1 / nu - 2 * Delta * ka2 * (2 * nu + Delta) / (nu^2)))

    H1 <- matrix(0, ncol = (p + g + 1), nrow = (p + g + 1))
    TLinv.dotL <- as.list(numeric(g1 + g2 + g3))
    for (s in 1:(g1 + g2 + g3)) {
      TLinv.dotL[[s]] <- TLam.inv %*% Tdot.L[, , s]
      for (i in 1:N)
      {
        if (i == 1) {
          idx1 <- 1:cumsum.ni[1]
        } else {
          idx1 <- (cumsum.ni[i - 1] + 1):cumsum.ni[i]
        }
        # H beta gamma
        H1[1:p, (p + s)] <- H1[1:p, (p + s)] + t(rep(ka2[i], ni[i]) * X[idx1, ]) %*% TLam.inv[idx1, idx1] %*% ys.cent[idx1, ] %*% t(ys.cent[idx1, ]) %*% TLinv.dotL[[s]][idx1, idx1] %*% TLam.inv[idx1, idx1] %*% ys.cent[idx1, ] -
          t(rep(tau[i], ni[i]) * X[idx1, ]) %*% TLinv.dotL[[s]][idx1, idx1] %*% TLam.inv[idx1, idx1] %*% ys.cent[idx1, ]
        # H gamma gamma
        for (l in 1:s) {
          H1[(p + s), (p + l)] <- H1[(p + s), (p + l)] + sum(diag(TLinv.dotL[[l]][idx1, idx1] %*% TLinv.dotL[[s]][idx1, idx1])) / 2 +
            (t(rep(ka2[i], ni[i]) * ys.cent[idx1]) %*% (TLinv.dotL[[l]][idx1, idx1] %*% TLam.inv[idx1, idx1] %*% ys.cent[idx1, ] %*% t(ys.cent[idx1, ]) %*% TLinv.dotL[[s]][idx1, idx1] %*% TLam.inv[idx1, idx1] %*% ys.cent[idx1, ])) / 2 -
            (t(rep(tau[i], ni[i]) * ys.cent[idx1]) %*% (TLinv.dotL[[l]][idx1, idx1] %*% TLinv.dotL[[s]][idx1, idx1] %*% TLam.inv[idx1, idx1] + TLinv.dotL[[s]][idx1, idx1] %*% TLinv.dotL[[l]][idx1, idx1] %*% TLam.inv[idx1, idx1]) %*% ys.cent[idx1, ]) / 2
        }
        # H gamma nu
        H1[(p + s), p + g + 1] <- H1[(p + s), p + g + 1] + (t(rep((ka1[i] - ka2[i]), ni[i]) * ys.cent[idx1, ]) %*% TLinv.dotL[[s]][idx1, idx1] %*% TLam.inv[idx1, idx1] %*% ys.cent[idx1, ]) / 2
      }
    }
    H.gamma <- H.gamma + H1


    ### Score vector
    S1 <- numeric((p + g + 1))
    # S beta
    S1[1:p] <- S1[1:p] + t(rep(tau, ni) * X) %*% TLam.inv %*% ys.cent
    # S gamma
    for (s in 1:g)
    {
      S1[p + s] <- S1[p + s] - sum(diag(TLinv.dotL[[s]])) / 2 + (t(rep(tau, ni) * ys.cent) %*% TLinv.dotL[[s]] %*% TLam.inv %*% ys.cent) / 2
    }
    # S nu
    S1[p + g + 1] <- S1[p + g + 1] + sum((digamma((nu + ni) / 2) - digamma(nu / 2) - ni / nu - log(1 + Delta / (nu)) + tau * Delta / nu) / 2)
    # S alpha
    V <- V.fun(y.conv[m, ], Data.miss)
    pvi <- c(exp(V %*% alpha) / (1 + exp(V %*% alpha)))
    ppvi <- pvi * (1 - pvi)
    S.alpha <- t(V) %*% (R - pvi)
    for (i in 1:n)
    {
      H.alpha <- H.alpha + t(t(V[i, ])) %*% ppvi[i] %*% t(V[i, ])
    }

    S <- S + c(S1, c(S.alpha))
    ### SS
    # SS beta beta
    SS[1:p, 1:p] <- SS[1:p, 1:p] + S1[1:p] %*% t(S1[1:p])
    # SS beta gamma
    SS[1:p, (p + 1):(p + g)] <- SS[1:p, (p + 1):(p + g)] + S1[1:p] %*% t(S1[(p + 1):(p + g)])
    # SS beta nu
    SS[1:p, (p + g + 1)] <- SS[1:p, (p + g + 1)] + S1[1:p] %*% t(S1[p + g + 1])
    ## SS beta alpha
    SS[1:p, (p + g + 1 + 1):(p + g + 1 + ma)] <- SS[1:p, (p + g + 1 + 1):(p + g + 1 + ma)] + S1[1:p] %*% t(S.alpha)
    # SS gamma gamma
    SS[(p + 1):(p + g), (p + 1):(p + g)] <- SS[(p + 1):(p + g), (p + 1):(p + g)] + S1[(p + 1):(p + g)] %*% t(S1[(p + 1):(p + g)])
    # SS gamma nu
    SS[(p + 1):(p + g), (p + g + 1)] <- SS[(p + 1):(p + g), (p + g + 1)] + S1[(p + 1):(p + g)] %*% t(S1[p + g + 1])
    # SS gamma alpha
    SS[(p + 1):(p + g), (p + g + 1 + 1):(p + g + 1 + ma)] <- SS[(p + 1):(p + g), (p + g + 1 + 1):(p + g + 1 + ma)] + S1[(p + 1):(p + g)] %*% t(S.alpha)
    # SS nu nu
    SS[p + g + 1, p + g + 1] <- SS[p + g + 1, p + g + 1] + S1[p + g + 1] %*% t(S1[p + g + 1])
    # SS nu alpha
    SS[p + g + 1, (p + g + 1 + 1):(p + g + 1 + ma)] <- SS[p + g + 1, (p + g + 1 + 1):(p + g + 1 + ma)] + S1[p + g + 1] %*% t(S.alpha)
    # SS alpha
    SS[(p + g + 1 + 1):(p + g + 1 + ma), (p + g + 1 + 1):(p + g + 1 + ma)] <- SS[(p + g + 1 + 1):(p + g + 1 + ma), (p + g + 1 + 1):(p + g + 1 + ma)] + (S.alpha %*% t(S.alpha))
  }

  for (l in (p + 1):(p + g - 1)) for (s in (l + 1):(p + g)) H.gamma[l, s] <- H.gamma[s, l]
  H.gamma[(p + g + 1), (p + 1):(p + g)] <- t(H.gamma[(p + 1):(p + g), (p + g + 1)])
  H.gamma[(p + 1):(p + g), 1:p] <- t(H.gamma[1:p, (p + 1):(p + g)])

  H <- H + H.gamma
  H[(p + g + 1), 1:p] <- t(H[1:p, (p + g + 1)])
  H <- cbind(H, matrix(0, ncol = ma, nrow = (p + g + 1)))
  H <- rbind(H, matrix(0, ncol = p + g + 1 + ma, nrow = ma))
  H[(p + g + 1 + 1):(p + g + 1 + ma), (p + g + 1 + 1):(p + g + 1 + ma)] <- -H.alpha
  H.mean <- H / mc.size
  SS[(p + 1):(p + g), 1:p] <- t(SS[1:p, (p + 1):(p + g)])
  SS[(p + g + 1), 1:(p + g)] <- SS[1:(p + g), (p + g + 1)]
  SS[(p + g + 1 + 1):(p + g + 1 + ma), 1:(p + g + 1)] <- t(SS[1:(p + g + 1), (p + g + 1 + 1):(p + g + 1 + ma)])
  SS.mean <- SS / mc.size
  S.mean <- S / mc.size

  I.theta <- -H.mean - SS.mean + S.mean %*% t(S.mean)
  if (cor.type == "UNC") {
    # V.theta = solve(I.theta[1:(p+g-1),1:(p+g-1)])
    V.theta <- solve(I.theta[-(p + g), -(p + g)])
    sd.theta <- c(sqrt(diag(V.theta)))
  }
  if (cor.type == "CAR1" | cor.type == "ARp" | cor.type == "CS" | cor.type == "BAND1") {
    V.theta <- solve(I.theta)
    sd.theta <- c(sqrt(diag(V.theta)))
  }
  if (cor.type == "DEC") {
    V.theta <- solve(I.theta)
    sd.theta <- sqrt(diag(V.theta))
  }
  out <- rbind(EST, c(sd.theta))
  if (cor.type == "UNC") colnames(out) <- rep(c("beta", "d", "sigma", "nu", "alpha"), c(p, length(DD[vechD]), 1, 1, ma))
  if (cor.type == "ARp" | cor.type == "CS" | cor.type == "BAND1") colnames(out) <- rep(c("beta", "d", "sigma", "phi", "nu", "alpha"), c(p, length(DD[vechD]), 1, 1, 1, ma))
  SD <- list(out = out, se = c(sd.theta), I.theta = I.theta, V.theta = V.theta)
  # }
  return(SD)
}
