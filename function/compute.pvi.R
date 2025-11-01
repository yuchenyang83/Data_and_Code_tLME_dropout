#################################################################################
# 
# This R script computes the dropout probability for each subject at each 
# observation time point, based on the fitted dropout model.
#
#################################################################################
#' Compute per-visit dropout probabilities under a specified mechanism
#'
#' For each subject and each observation time, computes the dropout probability
#' based on a fitted dropout model and the chosen missingness mechanism
#' (`"MNAR"`, `"MAR"`, or `"MCAR"`). The function constructs the design matrix
#' \eqn{V} for the dropout model according to the mechanism and returns the
#' per-visit probabilities along with useful metadata.
#'
#' @param Data data.frame.
#'   Longitudinal data in **long format**. Must contain at least the columns:
#'   \itemize{
#'     \item \code{Subject} (integer or factor): subject identifier (1..N).
#'     \item \code{Var1} (numeric): response at each visit (may contain \code{NA}).
#'     \item \code{R} (integer/logical): dropout indicator per row (1 = missing/dropout at this visit, 0 = observed).
#'     \item \code{week} (numeric/integer): visit time index.
#'     \item \code{prevOI} (numeric/integer/logical): an example covariate used in the dropout model.
#'   }
#'   Additional columns are allowed and ignored here.
#'
#' @param X matrix.
#'   Fixed-effects design matrix aligned with \code{Data} row order; only its
#'   dimension (number of columns \eqn{p}) is used here. Must have \code{nrow(X) == nrow(Data)}.
#'
#' @param Z matrix.
#'   Random-effects design matrix aligned with \code{Data} row order; only its
#'   dimension (number of columns \eqn{q}) is referenced. Must have \code{nrow(Z) == nrow(Data)}.
#'
#' @param init.para list.
#'   Initial (or fitted) parameter container. The following components are
#'   expected (many are included for interface consistency; this function only
#'   needs \code{alpha}):
#'   \itemize{
#'     \item \code{Beta} (numeric): fixed-effect coefficients (unused here).
#'     \item \code{D} (matrix): random-effect covariance (unused here).
#'     \item \code{sigma} (numeric): error scale (unused here).
#'     \item \code{alpha} (numeric vector): coefficients for the dropout model
#'           using \eqn{V} (order must match the internally constructed \eqn{V}).
#'     \item \code{Phi}, \code{ga}, \code{nu}: additional parameters (unused here).
#'   }
#'
#' @param fit list.
#'   Object from a previous model fit that must contain at least
#'   \code{y.c} (numeric vector) — completed/pseudo responses aligned with
#'   \code{Data} row order; used to build lagged terms in \eqn{V}.
#'
#' @param cor.type character(1).
#'   Working correlation type label used only for printing. One of
#'   \code{"UNC"}, \code{"ARp"}, \code{"BAND1"}, \code{"CS"}. Default prints accordingly.
#'
#'
#'
#' @param mechanism character(1).
#'   Missingness mechanism controlling the construction of \eqn{V}. One of:
#'   \code{"MNAR"}, \code{"MAR"}, \code{"MCAR"}. Default includes all three in the matching.
#'   \itemize{
#'     \item \strong{MNAR}: \eqn{V = [1,\ \texttt{prevOI},\ y_{i,j-1},\ y_{ij}]}
#'     \item \strong{MAR}:  \eqn{V = [1,\ \texttt{prevOI},\ y_{i,j-1}]}
#'     \item \strong{MCAR}: \eqn{V = [1]}
#'   }
#'
#' @details
#' Let \eqn{\alpha} be the parameter vector for the dropout (selection) model.
#' Conditional dropout probability at each row is computed via the logistic link:
#' \deqn{p = \frac{\exp(V \alpha)}{1 + \exp(V \alpha)}.}
#' The design \eqn{V} is internally (re)built per subject using the completed
#' responses \code{fit$y.c} to form lagged terms. Subjects and rows must be
#' aligned across \code{Data}, \code{X}, \code{Z}, and \code{fit$y.c}.
#'
#' @return \item{pvi}{Numeric vector of length \code{nrow(Data)}. The per-row dropout
#'   probability \eqn{p} computed from the logistic model with the internally
#'   constructed \eqn{V}.}
compute.pvi <- function(Data, X, Z, init.para, fit, cor.type = c("UNC", "ARp", "BAND1", "CS"), mechanism = c("MNAR", "MCAR")) {
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
