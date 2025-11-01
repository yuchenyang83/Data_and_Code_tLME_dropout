#################################################################################
# 
# Load the AIDS dataset from 'Data_and_Code/Data/source/aids.RData'.
# You may replace this with your own data if desired.
#
# This R script allows you to perform model fitting for both tLME and LME models 
# under various missing data mechanisms: MNAR, MAR, and MCAR.
#
# The models can be fitted with one of the following within-subject error structures: UNC, CS, MA(1), or AR(1).
#
# The fitted model results will be saved to 'Data_and_Code/Data/fit_result.RData'.
#
#################################################################################


# rm(list = ls())
PATH <- getwd()
load(paste0(PATH, "/Data/source/aids.RData"))
aids$death.time <- aids$time
aids$time <- aids$obstime
aids$Y <- aids$CD4

n <- length(table(aids$time))
N <- length(unique(aids$id))
aids.mat <- aids.time <- aids.D <- NULL
kk <- 0
for (i in 1:N)
{
  Y <- aids[aids$id == unique(aids$id)[i], ]$Y
  time <- aids[aids$id == unique(aids$id)[i], ]$time
  drug <- unique(aids[aids$id == unique(aids$id)[i], ]$drug)
  aids.mat <- rbind(aids.mat, c(Y, rep(NA, n - length(Y)), drug))
  aids.time <- rbind(aids.time, c(time, rep(NA, n - length(time))))
  aids.D <- rbind(aids.D, c(max(time), drug))
  if (max(aids[aids$id == unique(aids$id)[i], ]$time) == 6) kk <- kk + 1
}
aids.mat[is.na(aids.mat[, 5]), ]
aids.time[is.na(aids.time[, 5]), ]

table(aids.D[, 2], aids.D[, 1])
table(aids.D[, 1]) / N
N - cumsum(colSums(table(aids.D[, 2], aids.D[, 1])))
table(aids$time) / N
1 - table(aids$time) / N
### ### ### ### ### ###
### ### ### ### ### ###  trajectory plot
### ### ### ### ### ###
id <- unique(aids$id)
n <- length(id)
nj <- numeric(n)
for (i in 1:n) nj[i] <- length(aids$time[aids$id == id[i]])
nj == table(aids$id)
aids$Subject <- rep(1:n, nj)
aids$Dropout <- rep(aids.D[, 1], nj)


######## ######## ########
######## ########
######## ########
library(nlme)
levels(aids$gender) <- c(0, 1)
levels(aids$prevOI) <- c(0, 1)
levels(aids$AZT) <- c(0, 1)
aids <- data.frame(aids$Subject, aids$time, aids$Y, aids$drug, aids$gender, aids$prevOI, aids$AZT)
colnames(aids) <- c("Subject", "week", "Var1", "drug", "gender", "prevOI", "AZT")
aids <- groupedData(Var1 ~ week | Subject, data = aids)
aids$drug <- as.character(aids$drug)
aids$drug[which(aids$drug == "ddI")] <- 0
aids$drug[which(aids$drug == "ddC")] <- 1
aids$drug <- as.numeric(aids$drug)
aids$D <- rep(aids.D[, 1], nj)
aids$gender <- as.character(aids$gender)
aids$gender <- as.numeric(aids$gender)
aids$prevOI <- as.character(aids$prevOI)
aids$prevOI <- as.numeric(aids$prevOI)
aids$AZT <- as.character(aids$AZT)
aids$AZT <- as.numeric(aids$AZT)
aids$week <- aids$week

table(aids$drug, aids$week)
table(aids$drug, aids$D)

cumsum.nj <- cumsum(nj)

###########################################################################
###########################################################################
##############################
##############################  Selection model
##############################
##############################
#############################
#############################
D.max <- max(aids$week)
aids.miss <- NULL
for (i in 1:n)
{
  if (unique(aids[aids$Subject == i, ]$D) != D.max) {
    aids.i <- aids[aids$Subject == i, ]
    if (aids.i$D[1] == 0) {
      week <- aids.i$week
      aids.i <- rbind(aids.i, aids.i[dim(aids.i)[1], ])
      aids.i$week <- c(week, 2)
      aids.i$Var1[which(aids.i$week == 2)] <- NA
    }
    if (aids.i$D[1] == 2) {
      week <- aids.i$week
      aids.i <- rbind(aids.i, aids.i[dim(aids.i)[1], ])
      aids.i$week <- c(week, 6)
      aids.i$Var1[which(aids.i$week == 6)] <- NA
    }
    if (aids.i$D[1] == 6) {
      week <- aids.i$week
      aids.i <- rbind(aids.i, aids.i[dim(aids.i)[1], ])
      aids.i$week <- c(week, 12)
      aids.i$Var1[which(aids.i$week == 12)] <- NA
    }
    if (aids.i$D[1] == 12) {
      week <- aids.i$week
      aids.i <- rbind(aids.i, aids.i[dim(aids.i)[1], ])
      aids.i$week <- c(week, 18)
      aids.i$Var1[which(aids.i$week == 18)] <- NA
    }
    aids.i <- cbind(aids.i, miss = 1)
    aids.miss <- rbind(aids.miss, aids.i)
  } else {
    aids.i <- aids[aids$Subject == i, ]
    aids.i <- cbind(aids.i, miss = 0)
    aids.miss <- rbind(aids.miss, aids.i)
  }
}

na.ind <- is.na(as.vector(t(aids.miss$Var1)))
aids.miss[aids.miss$Subject == 25, ]
aids.miss[aids.miss$Subject == 17, ]

aids.miss <- cbind(aids.miss, R = 0)
aids.miss$R[is.na(aids.miss$Var1)] <- 1


aids.miss$Subject <- as.numeric(as.character(aids.miss$Subject))
aids.miss$week <- aids.miss$week
aids.miss$Var1 <- aids.miss$Var1


Data <- aids.miss
Data <- groupedData(Var1 ~ week | Subject, data = Data)
fm1 <- lme(Var1 ~ week + prevOI + drug + gender + AZT, data = Data[which(Data$R != 1), ], random = ~ week | Subject, method = "ML")
summary(fm1)
intervals(fm1)
fm2 <- lme(Var1 ~ week + prevOI + drug + gender + AZT, data = Data[which(Data$R != 1), ], random = ~ week | Subject, correlation = corCAR1(), method = "ML")
summary(fm2)
intervals(fm2)

Data$Subject <- as.numeric(as.character(Data$Subject))
class(Data)
Data <- as.data.frame(Data)


#############################
#############################  Initialization
X <- cbind(1, Data$week, Data$prevOI, Data$drug, Data$gender, Data$AZT)
# Z = as.matrix(X[,1])     # RI
Z <- X[, c(1, 2)] # RIS
p <- ncol(X)
q <- ncol(Z)
g <- 1

Beta <- matrix(c(fm2$coefficients$fixed), ncol = 1)
DD <- matrix(0, q, q)
DD[1, 1] <- as.numeric(VarCorr(fm2)[1, 1])
DD[2, 2] <- as.numeric(VarCorr(fm2)[2, 1])
DD[1, 2] <- DD[2, 1] <- as.numeric(VarCorr(fm2)[2, 3]) * sqrt(as.numeric(VarCorr(fm2)[1, 1]) * as.numeric(VarCorr(fm2)[2, 1]))
sigma <- c(fm2$sigma)
Phi <- runif(g)
ga <- 1
cor.type <- "UNC"
tol <- 1e-6
max.iter <- 1000
per <- 100
M <- 1000 ## Metropolis algorithm
M.LL <- 500 ## importance sampling

## prediction of missing response
source(paste0(PATH, "/function/analyze_realdata_AIDS.r"))
alpha <- as.numeric(prediction_ym(Data, fm2, X, Z, "ARp")$alpha.hat)
Data$yc <- prediction_ym(Data, fm2, X, Z, "ARp")$yc


##### run MNAR
source(paste0(PATH, "/function/tLMMmissingSEM.R"))
source(paste0(PATH, "/function/LMMmissingSEM.R"))
############################################
###########  UNC
set.seed(20240102)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = alpha, nu = 20)
mechanism <- "MNAR"
# init.para = fit.MNAR$para.est
# init.para = list(Beta=fit.MNAR$para.est$Beta, DD=fit.MNAR$para.est$D, sigma=fit.MNAR$para.est$sigma, Phi=Phi, alpha=fit.MNAR$para.est$alpha)
fit.MNAR <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("UNC"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

#### run MAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = alpha[-length(alpha)])
fit.MAR <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("UNC"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

#### run MCAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = c(-5.3))
fit.MCAR <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("UNC"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


############################################
###########  ARp
set.seed(20240101)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = alpha)
mechanism <- "MNAR"
cor.type <- "ARp"
# init.para = fit.MNAR$para.est
# init.para = list(Beta=fit.MNAR$para.est$Beta, DD=fit.MNAR$para.est$D, sigma=fit.MNAR$para.est$sigma, Phi=Phi, alpha=fit.MNAR$para.est$alpha)
fit.MNAR.ARp <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("ARp"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

#### run MAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = alpha[-length(alpha)])
fit.MAR.ARp <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("ARp"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

#### run MCAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = c(-5.3))
fit.MCAR.ARp <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("ARp"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)

############################################
###########  "BAND1"
set.seed(20240101)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0.2, ga = ga, alpha = alpha, nu = 20)
mechanism <- "MNAR"
cor.type <- c("BAND1")
fit.MNAR.BAND1 <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("BAND1"), M = 1000, M.LL = M.LL, tol = tol, max.iter = 100, per = per,
  mechanism = "MNAR"
)

#### run MAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0, ga = ga, alpha = alpha[-length(alpha)])
fit.MAR.BAND1 <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("BAND1"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)
# fit.MAR.treat = LMM.AECM.miss(Data, X, Z, init.para, cor.type = c("UNC"), M = 1000, P=1, tol = 1e-2, max.iter=max.iter, per=per)

#### run MCAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0, ga = ga, alpha = c(-2.15), nu = 20)
fit.MCAR.BAND1 <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("BAND1"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)

############################################
###########  "CS"
set.seed(20240101)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0.2, ga = ga, alpha = alpha, nu = 20)
mechanism <- "MNAR"
cor.type <- c("CS")
fit.MNAR.CS <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("CS"), M = 1000, M.LL = M.LL, tol = tol, max.iter = 200, per = per,
  mechanism = "MNAR"
)
set.seed(20240101)
#### run MAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0, ga = ga, alpha = alpha[-length(alpha)])
fit.MAR.CS <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("CS"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)
# fit.MAR.treat = LMM.AECM.miss(Data, X, Z, init.para, cor.type = c("UNC"), M = 1000, P=1, tol = 1e-2, max.iter=max.iter, per=per)
set.seed(20240101)
#### run MCAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0, ga = ga, alpha = c(-2.15), nu = 20)
fit.MCAR.CS <- LMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("CS"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)

########################################################################################
##### run Student t
########################################################################################
##### run MNAR
############################################
###########  UNC
set.seed(20240101)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = alpha, nu = 20)
mechanism <- "MNAR"
est.MNAR <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("UNC"), M = 1000, M.LL = 1000, tol = tol, max.iter = 100, per = per,
  mechanism = "MNAR"
)

#### run MAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = alpha[-length(alpha)], nu = 20)
est.MAR <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("UNC"), M = 1000, M.LL = 1000, tol = tol, max.iter = 100, per = per,
  mechanism = "MAR"
)
# fit.MAR.drug = LMM.AECM.miss(Data, X, Z, init.para, cor.type = c("UNC"), M = 1000, P=1, tol = 1e-2, max.iter=max.iter, per=1)

#### run MCAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = c(-5.3), nu = 20)
est.MCAR <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("UNC"), M = 1000, M.LL = 1000, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


############################################
###########  ARp
set.seed(20240101)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = alpha, nu = 20)
mechanism <- "MNAR"
cor.type <- "ARp"
est.MNAR.ARp <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("ARp"), M = 1000, M.LL = 1000, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

#### run MAR
mechanism <- "MNAR"
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = alpha[-length(alpha)], nu = 20)
est.MAR.ARp <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("ARp"), M = 1000, M.LL = 1000, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

#### run MCAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, ga = ga, alpha = c(-5.3), nu = 20)
est.MCAR.ARp <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("ARp"), M = 1000, M.LL = 1000, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)

############################################
###########  "BAND1"
set.seed(20240101)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0.2, ga = ga, alpha = alpha, nu = 20)
mechanism <- "MNAR"
cor.type <- c("BAND1")
est.MNAR.BAND1 <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("BAND1"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)

#### run MAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0.2, ga = ga, alpha = alpha[-length(alpha)], nu = 20)
est.MAR.BAND1 <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("BAND1"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

#### run MCAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0.2, ga = ga, alpha = c(-2.15), nu = 20)
est.MCAR.BAND1 <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("BAND1"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)

############################################
###########  "CS"
set.seed(20240101)
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0.2, ga = ga, alpha = alpha, nu = 20)
mechanism <- "MNAR"
cor.type <- c("CS")
est.MNAR.CS <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("CS"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MNAR"
)
set.seed(20240101)
#### run MAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0.2, ga = ga, alpha = alpha[-length(alpha)], nu = 20)
est.MAR.CS <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("CS"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MAR"
)

set.seed(20240101)
#### run MCAR
init.para <- list(Beta = Beta, DD = DD, sigma = sigma, Phi = 0.2, ga = ga, alpha = c(-2.15), nu = 20)
est.MCAR.CS <- tLMM.miss.SEM(Data, X, Z, init.para,
  cor.type = c("CS"), M = 1000, M.LL = M.LL, tol = tol, max.iter = max.iter, per = per,
  mechanism = "MCAR"
)


c(fit.MNAR$model.inf$bic, fit.MAR$model.inf$bic, fit.MCAR$model.inf$bic)
c(fit.MNAR.ARp$model.inf$bic, fit.MAR.ARp$model.inf$bic, fit.MCAR.ARp$model.inf$bic)
c(fit.MNAR.CS$model.inf$bic, fit.MAR.CS$model.inf$bic, fit.MCAR.CS$model.inf$bic)
c(fit.MNAR.BAND1$model.inf$bic, fit.MAR.BAND1$model.inf$bic, fit.MCAR.BAND1$model.inf$bic)


c(est.MNAR$model.inf$bic, est.MAR$model.inf$bic, est.MCAR$model.inf$bic)
c(est.MNAR.ARp$model.inf$bic, est.MAR.ARp$model.inf$bic, est.MCAR.ARp$model.inf$bic)
c(est.MNAR.CS$model.inf$bic, est.MAR.CS$model.inf$bic, est.MCAR.CS$model.inf$bic)
c(est.MNAR.BAND1$model.inf$bic, est.MAR.BAND1$model.inf$bic, est.MCAR.BAND1$model.inf$bic)

save.image(paste0(PATH, "/Data/fit_result.RData"))

