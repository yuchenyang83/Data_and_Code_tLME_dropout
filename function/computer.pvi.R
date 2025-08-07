#################################################################################
# 
# This R script computes the dropout probability for each subject at each 
# observation time point, based on the fitted dropout model.
#
#
#################################################################################
compute.pvi <- function(Data, X, Z, V, g = g, init.para, fit, cor.type = c("UNC", "ARp", "BAND1", "CS"), M = 100,
                        M.LL = 1000, P = 1, tol = 1e-6, max.iter = max.iter, per = 1, mechanism = c("MNAR", "MCAR")) {
  begin <- proc.time()[1]
  # initial values of parameter
  Beta <- init.para$Beta
  DD <- init.para$D
  sigma <- init.para$sigma
  alpha <- init.para$alpha
  Phi <- init.para$Phi
  ga <- init.para$ga
  nu <- init.para$nu
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


  if (mechanism == "MNAR") cat(rep("=", 25), "tLMM (MNAR) with ", cor.type, " errors is fitted...; ", "missing = ", num.na / N * 100, "%", rep("=", 25), sep = "", "\n")
  if (mechanism == "MAR") cat(rep("=", 25), "tLMM (MAR) with ", cor.type, " errors is fitted...; ", "missing = ", num.na / N * 100, "%", rep("=", 25), sep = "", "\n")
  if (mechanism == "MCAR") cat(rep("=", 25), "tLMM (MCAR) with ", cor.type, " errors is fitted...; ", "missing = ", num.na / N * 100, "%", rep("=", 25), sep = "", "\n")

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
        V.i <- cbind(1, Data.i$prevOI, fData.i, bData.i)
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
      V <- cbind(rep(1, k), Data.i$prevOI, fData.i, bData.i)
      return(V)
    }
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
        V.i <- cbind(1, Data.i$prevOI, fData.i)
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
      V <- cbind(rep(1, k), Data.i$prevOI, fData.i)
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
  yc <- fit$y.c
  V <- V.fun(yc, Data.miss)

  pvi <- c(exp(V %*% alpha) / (1 + exp(V %*% alpha)))
  return(pvi)
}