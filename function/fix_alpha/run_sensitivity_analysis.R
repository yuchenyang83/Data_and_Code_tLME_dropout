vech.posi <- function(dim) cbind(rep(1:dim, 1:dim), unlist(mapply(":", 1, 1:dim)))

OM.fun <- function(Data, n, nj, cumsum.nj, N) {
  TTo <- rep(1, N)
  TTo[is.na(Data$Var1)] <- 0
  TTm <- rep(0, N)
  TTm[is.na(Data$Var1)] <- 1
  TO <- TM <- list(n)
  for (j in 1:n)
  {
    if (j == 1) {
      idx1 <- 1:cumsum.nj[1]
    } else {
      idx1 <- (cumsum.nj[j - 1] + 1):cumsum.nj[j]
    }
    TO[[j]] <- TTo[idx1]
    TM[[j]] <- TTm[idx1]
  }
  return(list(TO = TO, TM = TM))
}

OM.fun1 <- function(Data, n, nj, cumsum.nj, N, R) {
  TTo <- rep(1, N)
  TTo[is.na(Data$Var1)] <- 0
  TTm <- rep(0, N)
  TTm[is.na(Data$Var1)] <- 1
  TO <- TM <- list(n)
  for (j in 1:n)
  {
    if (j == 1) {
      idx1 <- 1:cumsum.nj[1]
    } else {
      idx1 <- (cumsum.nj[j - 1] + 1):cumsum.nj[j]
    }
    TO[[j]] <- matrix(diag(TTo[idx1])[R[idx1] == 0, ], nrow = sum(R[idx1] == 0))
    TM[[j]] <- matrix(diag(TTm[idx1])[R[idx1] == 1, ], nrow = sum(R[idx1]))
  }
  return(list(TO = TO, TM = TM))
}

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

log.pos.y <- function(yi, Ri, wi, Xi, Lam.i, Beta) {
  fy <- log(mvtnorm::dmvnorm(yi, mean = c(Xi %*% Beta), sigma = Lam.i))
  fr <- log(prod(wi^Ri) * prod((1 - wi)^(1 - Ri)))
  return(fy + fr)
}

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

# Autocorrelation matrix:
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

arp.Ci.dot <- function(phi, dim, k) -arp.Ci(phi, dim = dim) %*% d.inv.Ci.fn(phi, dim = dim)[[k]] %*% arp.Ci(phi, dim = dim)

Arp.Ci.dot <- function(phi, dim) {
  Ti <- 1:dim
  tem <- abs(outer(Ti, Ti, "-"))
  dot.CAR <- tem * (phi^(tem - 1))
  return(dot.CAR)
}
cs.Ci.dot <- function(phi, dim) {
  dot.CS <- matrix(1, dim, dim)
  diag(dot.CS) <- 0
  return(dot.CS)
}
BAND1.Ci.dot <- function(phi, dim) {
  tmp <- c(rep(0, dim - 2), 1, 0, 1, rep(0, dim - 2))
  dot.BAND1 <- NULL
  for (i in 1:dim) dot.BAND1 <- rbind(dot.BAND1, tmp[1:dim + (dim - i)])
  return(dot.BAND1)
}

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

# Q t function
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

# Q normal function
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

# dot DEC:
DEC.dot.phi <- function(phi, ga, Ti) {
  tem <- abs(outer(Ti, Ti, "-"))^ga
  dot.CAR <- tem * (phi^(tem - 1))
  return(dot.CAR)
}
DEC.dot.ga <- function(phi, ga, Ti) {
  tem <- abs(outer(Ti, Ti, "-"))
  dot.ga.CAR <- log(phi^tem) * cor.fn(phi = phi, type = "DEC", Ti = Ti, ga = ga)
  return(dot.ga.CAR)
}

# Generate LMMs Data:
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

# Generate tLMMs Data:
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


# MH algorithm for MCEM
yv.joint.log <- function(yi, Xbeta.i, Lam.i, Vi, R.i, alpha, n.i) {
  fy <- mvtnorm::dmvnorm(yi, mean = c(Xbeta.i), sigma = Lam.i, log = T)
  pvi <- c(exp(as.matrix(Vi) %*% alpha) / (1 + exp(as.matrix(Vi) %*% alpha)))
  fv <- log(prod(pvi^R.i) * prod((1 - pvi)^(1 - R.i)))
  return(fy + fv)
}

MH.y.miss2 <- function(yi, Xbeta.i, Lam.i, alpha, Vi, R.i, n.i, mu.mo, Sig.mm.o, na.idx.i) {
  yi.old <- yi.new <- yi
  log.den.old <- yv.joint.log(yi.old, Xbeta.i, Lam.i, Vi, R.i, alpha, n.i)
  if (!is.na(na.idx.i[1])) {
    lna <- length(na.idx.i)
    if (lna == 1) ym.hat <- c(rnorm(1, mean = mu.mo[na.idx.i], sd = sqrt(Sig.mm.o[na.idx.i, na.idx.i])))
    if (lna > 1) ym.hat <- c(mvtnorm::rmvnorm(1, mu.mo[na.idx.i], Sig.mm.o[na.idx.i, na.idx.i]))
    yi.new[which(R.i == 1)] <- ym.hat
  }
  log.den.new <- yv.joint.log(yi.new, Xbeta.i, Lam.i, Vi, R.i, alpha, n.i)
  log.accept <- log.den.new - log.den.old
  if (log(runif(1)) > log.accept) yi.new <- yi.old
  return(yi.new)
}



yv.joint.log.t <- function(yi, Xbeta.i, Lam.i, tau.i, Vi, R.i, alpha, n.i) {
  fy <- mvtnorm::dmvnorm(yi, mean = c(Xbeta.i), sigma = Lam.i / tau.i, log = T)
  pvi <- c(exp(Vi %*% alpha) / (1 + exp(Vi %*% alpha)))
  if (sum(c(exp(Vi %*% alpha)) == Inf) > 0) pvi[which(c(exp(Vi %*% alpha)) == Inf)] <- 1
  if (sum(pvi == 1) > 0) pvi[which(pvi == 1)] <- 1 - 1e-16
  if (sum(pvi == 0) > 0) pvi[which(pvi == 0)] <- 1e-16
  fv <- log(prod(pvi^R.i) * prod((1 - pvi)^(1 - R.i)))
  return(fy + fv)
}

MH.y.miss.t <- function(yi, Xbeta.i, Lam.i, tau.i, alpha, Vi, R.i, n.i, mu.mo, Sig.mm.o, na.idx.i) {
  yi.old <- yi.new <- yi
  log.den.old <- yv.joint.log.t(yi.old, Xbeta.i, Lam.i, tau.i, Vi, R.i, alpha, n.i)
  if (!is.na(na.idx.i[1])) {
    lna <- length(na.idx.i)
    if (lna == 1) ym.hat <- c(rnorm(1, mean = mu.mo[na.idx.i], sd = sqrt(Sig.mm.o[na.idx.i, na.idx.i] / tau.i)))
    if (lna > 1) ym.hat <- c(mvtnorm::rmvnorm(1, mu.mo[na.idx.i], Sig.mm.o[na.idx.i, na.idx.i] / tau.i))
    yi.new[which(R.i == 1)] <- ym.hat
  }
  log.den.new <- yv.joint.log.t(yi.new, Xbeta.i, Lam.i, tau.i, Vi, R.i, alpha, n.i)
  log.accept <- log.den.new - log.den.old
  if (log(runif(1)) > log.accept) yi.new <- yi.old
  return(yi.new)
}

# Q-function for nu
nu.Q.fn <- function(nu, kappa, tau, N) {
  fn.nu <- nu * (kappa[1] - tau[1] + log(nu / 2)) - 2 * log(gamma(nu / 2))
  for (i in 2:N) {
    fn.nu <- fn.nu + (nu * (kappa[i] - tau[i] + log(nu / 2)) - 2 * log(gamma(nu / 2)))
  }
  return(-fn.nu)
}
