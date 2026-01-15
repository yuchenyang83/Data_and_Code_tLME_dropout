#################################################################################
# 
# This R script provides a demonstration of the proposed method using a 
# simplified simulated dataset.
#
# It serves as a minimal working example illustrating the full workflow, 
# including data generation, model fitting, and result visualization, for 
# both tLME and LME models under MNAR mechanisms.
#
# For detailed explanation and interpretation, please refer to 
# Supplementary Appendix C, which demonstrates the practical implementation 
# of our method using this example.
#
#################################################################################

# Load custom functions
PATH <- getwd()
source(paste0(PATH, "/function/simulation/simulate_dropout_study.r"))
source(paste0(PATH, "/function/simulation/tLMMmissingSEM_simulation.r"))
source(paste0(PATH, "/function/simulation/LMMmissingSEM_simulation.r"))

# Set sample size
n <- 50

# True parameter settings:
p = 4; si = 6; q = 2; g = 1
Beta <- matrix(c(5.35 , -0.3, -0.1, -0.65), ncol = g) # Fixed effects
DD <- matrix(c(0.4, 0.1, 0.1, 0.2), nrow = q) # Covariance matrix
sigma <- 0.5 # Residual standard deviation
nu <- 5 # Degrees of freedom for t-distribution
Phi <- 1e-6 # UNC
alpha <- c(-0.72, -0.69, -0.30, 0.30) # Dropout process parameters
para <- list(alpha = alpha , Beta = Beta , DD = DD , sigma = sigma , Phi = Phi , nu = nu)

# Generate longitudinal data
set.seed (20250301)
gen.Data <- gen.tlmm(n, para , cor.type = "UNC", si , q)
Data <- gen.Data$Data
Ymat <- gen.Data$Ymat

dim(Ymat)
head(Ymat)

gen.Data.miss <- add.MNAR.tLMM(gen.Data$Data, gen.Data$Ymat, si, para$alpha)
gen.Data.miss$Ymat.na[1:10 ,]

Data.miss <- gen.Data.miss$Data.sub
Data.miss$R <- rep(0, nrow(Data.miss))
Data.miss$R[is.na(Data.miss$Var1)] <- 1 # Indicator for missingness
Data.miss$week <- sqrt(Data.miss$Time) # Transform time for better model stability
X <- gen.Data.miss$X
Z <- gen.Data.miss$Z
Data <- Data.miss
sum(Data$R)/n


k1 <- ggplot(data = Data[which(Data$R != 1), ], aes(x = Time , y = Var1 , col = as.factor(
  treat), shape = as.factor(treat), group = as.factor(Subject))) +
  geom_line() + geom_point(size = 3) + xlab("Time") +
  ylab(expression(paste("Measurements of ", y[ij]^"o"))) +
  theme(legend.position = "top") +
  scale_color_discrete(labels = c("Placebo", "Therapy")) +
  scale_shape_discrete(labels = c("Placebo", "Therapy")) +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit (1.5, "cm")) +
  theme(legend.title = element_blank ())
k1

k2 <- ggplot(data = Data[-which(is.na(Data$Var1)), ], aes(x = as.factor(Time), y = Var1 ,colour = as.factor(
  treat), shape = as.factor(treat), grop = as.factor(Time))) +
  geom_boxplot(na.rm = T) + xlab(NULL) + ylab(expression(paste("Measurements of ", y[ij]^"o"))) +
  theme(legend.position = "top") +
  scale_color_discrete(labels = c("Placebo", "Therapy")) +
  scale_shape_discrete(labels = c("Placebo", "Therapy")) +
  theme(text = element_text(size = 20)) +
  theme(legend.key.size = unit(1, "cm"),
        legend.key.width = unit (1.5, "cm")) +
  theme(legend.title = element_blank ())
k2

k3 <- ggplot(data = Data[-which(is.na(Data$Var1)), ], aes(x = as.factor(Time), fill = as.factor(treat))) +
  geom_bar(position = position_dodge ()) + xlab("Time") + ylab("Number of observed responses") +
  geom_text(stat = "count", aes(label = paste0(round (after_stat(count)/ n, 3) * 100, "%")),
            vjust = 1.6, color = "black", size = 5, position = position_dodge (0.9)) +
  theme(legend.position = "top") +
  scale_fill_discrete(name = NULL , labels = c("Placebo", "Therapy")) +
  theme(text = element_text(size = 20))
k3

# Fit LME model to observed data
fm2 <- lme(Var1 ~ week + treat * week, data = Data[which(Data$R != 1), ], 
           random = ~ week | Subject, method = "ML")

# Estimate dropout process parameters
alpha1 <- as.numeric(prediction_ym(Data, fm2, X, Z, "UNC")$alpha.hat)
Data$yc <- prediction_ym(Data, fm2, X, Z, "UNC")$yc

DD1 <- matrix(0, q, q)
DD1[1, 1] <- as.numeric(VarCorr(fm2)[1, 1])
DD1[2, 2] <- as.numeric(VarCorr(fm2)[2, 1])
DD1[1, 2] <- DD1[2, 1] <- as.numeric(VarCorr(fm2)[2, 3]) * sqrt(as.numeric(VarCorr(fm2)[1, 1]) * as.numeric(VarCorr(fm2)[2, 1]))
init.para <- list(Beta = matrix(c(fm2$coefficients$fixed), ncol = 1), DD = DD1, sigma = c(fm2$sigma), Phi = para$Phi, nu = para$nu, alpha = alpha1)

cor.type = "UNC"
# Fit t LME model with MNAR dropout
cat(rep("=", 15), "Student's t of Linear mixed models with MNAR missing", cor.type[1], " errors: ", "\n")
est.MNAR <- tLMM.miss.SEM(Data, X, Z, V, g, init.para = init.para, cor.type = c("UNC"), M = 1000, M.LL = 1000, tol = 1e-6, max.iter = 100, per = 20, mechanism = "MNAR")
est.MNAR$model.inf$bic; est.MNAR$MSE.y

# Fit LME model with MNAR dropout
cat(rep("=", 15), "Linear mixed models with MNAR missing", cor.type[1], " errors: ", "\n")
fit.MNAR <- LMM.miss.SEM(Data, X, Z, V, g, init.para = init.para, cor.type = c("UNC"), M = 1000, M.LL = 1000, tol = 1e-6, max.iter = 100, per = 20, mechanism = "MNAR")
fit.MNAR$model.inf$bic; fit.MNAR$MSE.y
