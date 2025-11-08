################################################################################
#
#   Filename    :    simSEM25.R
#   Project    : BiomJ article "Extending t linear mixed models for longitudinal 
#                data with non-ignorable dropout applied to AIDS studies"                                                           
#   Authors     :    Yu-Chen Yang and Wan-Lun Wang and Luis M. Castro and Tsung-I Lin
#   Date        :    07.08.2025
#   Purpose     :    re-generate the results files for the simulation
#
#   Input code files  :  Data_and_Code/function/simulation/simulate_dropout_study.r;
#                        Data_and_Code/function/simulation/tLMMmissingSEM.r;
#                        Data_and_Code/function/simulation/LMMmissingSEM.r
#
#   Output .txt files :  Data_and_Code/Data/Simulation/SS-simulationSEM-t25/SIM1/...;
#                        Data_and_Code/Data/Simulation/SS-simulationSEM-t25/SIM2/...;
#                        Data_and_Code/Data/Simulation/SS-simulationSEM-t25/SIM3/...;
#                        Data_and_Code/Data/Simulation/SS-simulationSEM-t25/SIM4/...;
#                        Data_and_Code/Data/Simulation/SS-simulationSEM-t25/SIM5/...
#
#   R Version   :    R-4.3.1
#
###############################################################################
PATH11 <- paste0(PATH, "/function/simulation/")

# Load custom functions
source(paste0(PATH11, "simulate_dropout_study.r"))
source(paste0(PATH11, "tLMMmissingSEM_simulation.r"))
source(paste0(PATH11, "LMMmissingSEM_simulation.r"))

# Parameter setting:
p <- 4
si <- 6
q <- 2
g <- 1
Beta <- matrix(c(5.35, -0.3, -0.1, -0.65), ncol = g) ## intercept; time; treat; time x treat
DD <- matrix(0, q, q)
DD[1, 1] <- 0.4
DD[2, 2] <- 0.2
DD[2, 1] <- DD[1, 2] <- 0.1
sigma <- 0.5
nu <- 5
Phi <- 1e-6

### 25 dropout
alpha <- c(-2.51, -0.69, -0.30, 0.30) ## intercept; treat; y_si-1; y_si
### 50 dropout
# alpha = c(-1.56, -0.69, -0.30,  0.30) ## intercept; treat; y_si-1; y_si
### 75 dropout
# alpha = c(-0.72, -0.69, -0.30,  0.30) ## intercept; treat; y_si-1; y_si
para <- list(alpha = alpha, Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, nu = nu)

cor.type <- "UNC"
max.iter <- 100
if (n == 400) max.iter <- 85

# Repp <- 100
# set.seed(seednum)
for (Rep in Rep:Repp) {
  seednum <- n * 1000 + 1000000 + Rep
  set.seed(seednum)
  if (n == 25) PATH1 <- paste0(PATH, "/Data/Simulation/SS-simulationSEM-t25/SIM1/")
  if (n == 50) PATH1 <- paste0(PATH, "/Data/Simulation/SS-simulationSEM-t25/SIM2/")
  if (n == 100) PATH1 <- paste0(PATH, "/Data/Simulation/SS-simulationSEM-t25/SIM3/")
  if (n == 200) PATH1 <- paste0(PATH, "/Data/Simulation/SS-simulationSEM-t25/SIM4/")
  if (n == 400) PATH1 <- paste0(PATH, "/Data/Simulation/SS-simulationSEM-t25/SIM5/")

  if (n == 25 & sum(Rep == c(50, 110, 89, 77, 51, 137, 160, 14, 49, 15, 28, 42)) == 1) next
  if (n == 400 & Rep == 73) next

  cat(paste(c(rep("=", 15), rep(" ", 3), "The ", Rep, " time simulation: n = ", n, rep(" ", 3), rep("=", 15)), sep = "", collapse = ""), "\n")

  ################################
  ################   The new simulation Data
  ################
  cor.type <- "UNC"
  gen.Data <- gen.tlmm(n, para, cor.type = "UNC", si, q)
  Data <- gen.Data$Data
  Ymat <- gen.Data$Ymat

  gen.Data.miss <- add.MNAR.tLMM(gen.Data$Data, gen.Data$Ymat, si, para$alpha)
  Data.miss <- gen.Data.miss$Data.sub
  Data.miss$R <- rep(0, nrow(Data.miss))
  Data.miss$R[is.na(Data.miss$Var1)] <- 1
  Data.miss$week <- sqrt(Data.miss$Time)
  X <- gen.Data.miss$X
  Z <- gen.Data.miss$Z
  Data <- Data.miss

  fm2 <- try(lme(Var1 ~ week + treat * week, data = Data[which(Data$R != 1), ], random = ~ week | Subject, method = "ML"), silent = F)
  if (class(fm2) == "try-error") next
  alpha1 <- as.numeric(prediction_ym(Data, fm2, X, Z, "UNC")$alpha.hat)
  Data$yc <- prediction_ym(Data, fm2, X, Z, "UNC")$yc


  ###################################################
  ###################################################  Run MNAR
  ###################################################
  # mechanism='MNAR'
  cat(rep("=", 15), "Student's t of Linear mixed models with MNAR missing", cor.type[1], " errors: ", "\n")
  DD1 <- matrix(0, q, q)
  DD1[1, 1] <- as.numeric(VarCorr(fm2)[1, 1])
  DD1[2, 2] <- as.numeric(VarCorr(fm2)[2, 1])
  DD1[1, 2] <- DD1[2, 1] <- as.numeric(VarCorr(fm2)[2, 3]) * sqrt(as.numeric(VarCorr(fm2)[1, 1]) * as.numeric(VarCorr(fm2)[2, 1]))
  init.para <- list(Beta = matrix(c(fm2$coefficients$fixed), ncol = 1), DD = DD1, sigma = c(fm2$sigma), Phi = para$Phi, nu = para$nu, alpha = alpha1)
  est.MNAR <- try(tLMM.miss.SEM(Data, X, Z, V, g, init.para,
    cor.type = c("UNC"), M = 1000, M.LL = 1000, P = 1, tol = 1e-6, max.iter = max.iter, per = 500,
    mechanism = "MNAR"
  ), silent = F)
  if (class(est.MNAR) == "try-error") next

  fit.MNAR <- try(LMM.miss.SEM(Data, X, Z, V, g, init.para,
    cor.type = c("UNC"), M = 1000, M.LL = 1000, P = 1, tol = 1e-6, max.iter = max.iter, per = 500,
    mechanism = "MNAR"
  ), silent = F)
  if (class(fit.MNAR) == "try-error") next

  ###################################################
  ###################################################  Run MAR
  ###################################################
  cat(rep("=", 15), "Student's t of Linear mixed models with MAR missing", cor.type[1], " errors: ", "\n")
  init.para <- list(Beta = para$Beta, DD = para$DD, sigma = para$sigma, Phi = para$Phi, nu = para$nu, alpha = para$alpha)
  init.para$alpha <- c(-2.51, -0.69, -0.30)
  est.MAR <- try(tLMM.miss.SEM(Data, X, Z, V, g, init.para,
    cor.type = c("UNC"), M = 1000, M.LL = 1000, P = 1, tol = 1e-6, max.iter = max.iter, per = 500,
    mechanism = "MAR"
  ), silent = F)
  if (class(est.MAR) == "try-error") next
  fit.MAR <- try(LMM.miss.SEM(Data, X, Z, V, g, init.para,
    cor.type = c("UNC"), M = 1000, M.LL = 1000, P = 1, tol = 1e-6, max.iter = max.iter, per = 500,
    mechanism = "MAR"
  ), silent = F)
  if (class(fit.MAR) == "try-error") next

  ###################################################
  ###################################################  Run MCAR
  ###################################################
  cat(rep("=", 15), "Student's t of Linear mixed models with MCAR missing", cor.type[1], " errors: ", "\n")
  init.para$Beta <- para$Beta # - 0.1 * runif(1)
  init.para$alpha <- c(-2.51)
  est.MCAR <- try(tLMM.miss.SEM(Data, X, Z, V, g, init.para,
    cor.type = c("UNC"), M = 1000, M.LL = 1000, P = 1, tol = 1e-6, max.iter = max.iter, per = 500,
    mechanism = "MCAR"
  ), silent = F)
  if (class(est.MCAR) == "try-error") next
  fit.MCAR <- try(LMM.miss.SEM(Data, X, Z, V, g, init.para,
    cor.type = c("UNC"), M = 1000, M.LL = 1000, P = 1, tol = 1e-6, max.iter = max.iter, per = 500,
    mechanism = "MCAR"
  ), silent = F)
  if (class(fit.MCAR) == "try-error") next

  if (Rep == 1) {
    para.real <- c(para$Beta, para$sigma, para$DD[vech.posi(q)], para$nu, para$alpha)
    write(c(para.real), paste(PATH1, "realpara.txt", sep = ""), ncol = 1000, append = T)
  }
  if (Rep == 2) {
    para.real <- c(para$Beta, para$sigma, para$DD[vech.posi(q)], para$nu, para$alpha)
    write(c(para.real), paste(PATH1, "realpara.txt", sep = ""), ncol = 1000, append = T)
  }

  MNAR.para.est <- MAR.para.est <- MCAR.para.est <- NULL
  MNAR.para.est <- c(est.MNAR$para.est$Beta, est.MNAR$para.est$sigma, c(est.MNAR$para.est$D[vech.posi(q)]), est.MNAR$para.est$nu, c(est.MNAR$para.est$alpha))
  MAR.para.est <- c(est.MAR$para.est$Beta, est.MAR$para.est$sigma, c(est.MAR$para.est$D[vech.posi(q)]), est.MAR$para.est$nu, c(est.MAR$para.est$alpha))
  MCAR.para.est <- c(est.MCAR$para.est$Beta, est.MCAR$para.est$sigma, c(est.MCAR$para.est$D[vech.posi(q)]), est.MCAR$para.est$nu, c(est.MCAR$para.est$alpha))

  MNAR.para.fit <- MAR.para.fit <- MCAR.para.fit <- NULL
  MNAR.para.fit <- c(fit.MNAR$para.est$Beta, fit.MNAR$para.est$sigma, c(fit.MNAR$para.est$D[vech.posi(q)]), c(fit.MNAR$para.est$alpha))
  MAR.para.fit <- c(fit.MAR$para.est$Beta, fit.MAR$para.est$sigma, c(fit.MAR$para.est$D[vech.posi(q)]), c(fit.MAR$para.est$alpha))
  MCAR.para.fit <- c(fit.MCAR$para.est$Beta, fit.MCAR$para.est$sigma, c(fit.MCAR$para.est$D[vech.posi(q)]), c(fit.MCAR$para.est$alpha))

  loglik <- c(est.MNAR$model.inf$loglik, est.MAR$model.inf$loglik, est.MCAR$model.inf$loglik, fit.MNAR$model.inf$loglik, fit.MAR$model.inf$loglik, fit.MCAR$model.inf$loglik)
  aic <- c(est.MNAR$model.inf$aic, est.MAR$model.inf$aic, est.MCAR$model.inf$aic, fit.MNAR$model.inf$aic, fit.MAR$model.inf$aic, fit.MCAR$model.inf$aic)
  bic <- c(est.MNAR$model.inf$bic, est.MAR$model.inf$bic, est.MCAR$model.inf$bic, fit.MNAR$model.inf$bic, fit.MAR$model.inf$bic, fit.MCAR$model.inf$bic)
  MSE.y <- c(est.MNAR$MSE.y, est.MAR$MSE.y, est.MCAR$MSE.y, fit.MNAR$MSE.y, fit.MAR$MSE.y, fit.MCAR$MSE.y)
  MAE.y <- c(est.MNAR$MAE.y, est.MAR$MAE.y, est.MCAR$MAE.y, fit.MNAR$MAE.y, fit.MAR$MAE.y, fit.MCAR$MAE.y)
  MAPE.y <- c(est.MNAR$MAPE.y, est.MAR$MAPE.y, est.MCAR$MAPE.y, fit.MNAR$MAPE.y, fit.MAR$MAPE.y, fit.MCAR$MAPE.y)
  iter <- c(est.MNAR$iter, est.MAR$iter, est.MCAR$iter, fit.MNAR$iter, fit.MAR$iter, fit.MCAR$iter)

  MNAR.se.est <- est.MNAR$IM$se
  MAR.se.est <- est.MAR$IM$se
  MCAR.se.est <- est.MCAR$IM$se
  MNAR.se.fit <- fit.MNAR$IM$se
  MAR.se.fit <- fit.MAR$IM$se
  MCAR.se.fit <- fit.MCAR$IM$se

  write(c(Rep, loglik), paste(PATH1, "loglik.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MSE.y), paste(PATH1, "MSE.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, seednum), paste(PATH1, "seed.txt", sep = ""), ncol = 10000, append = T)

  write(c(Rep, MNAR.para.est), paste(PATH1, "MNAR.para.est.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MAR.para.est), paste(PATH1, "MAR.para.est.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MCAR.para.est), paste(PATH1, "MCAR.para.est.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MNAR.para.fit), paste(PATH1, "MNAR.para.fit.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MAR.para.fit), paste(PATH1, "MAR.para.fit.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MCAR.para.fit), paste(PATH1, "MCAR.para.fit.txt", sep = ""), ncol = 10000, append = T)

  write(c(Rep, MNAR.se.est), paste(PATH1, "MNAR.se.est.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MAR.se.est), paste(PATH1, "MAR.se.est.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MCAR.se.est), paste(PATH1, "MCAR.se.est.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MNAR.se.fit), paste(PATH1, "MNAR.se.fit.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MAR.se.fit), paste(PATH1, "MAR.se.fit.txt", sep = ""), ncol = 10000, append = T)
  write(c(Rep, MCAR.se.fit), paste(PATH1, "MCAR.se.fit.txt", sep = ""), ncol = 10000, append = T)
}
