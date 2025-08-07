################################################################################
#                                                                                     
#   Filename    :    FigC1.R  												  
#   Project     :    BiomJ article "Extending t linear mixed models for longitudinal 
#                    data with non-ignorable dropout applied to AIDS studies"                                                           
#   Authors     :    Yu-Chen Yang and Wan-Lun Wang and Luis M. Castro and Tsung-I Lin
#   Date        :    07.08.2025
#   Purpose     :    produce Figure C1 in Appendix C
#
#   Input data files  :  None
#   Output data files :  Data_and_Code/results/FigureC1.pdf
#
#   R Version   :    R-4.3.1                                                              
#   Required R packages : ggplot2; ggtext
#
################################################################################ 

# Load custom functions
source(paste0(PATH, "/function/simulation/simulate_dropout_study.r"))
source(paste0(PATH, "/function/simulation/tLMMmissingSEM.r"))
source(paste0(PATH, "/function/simulation/LMMmissingSEM.r"))

n <- 50

# True parameter settings
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
alpha <- c(-0.72, -0.69, -0.30, 0.30)
para <- list(alpha = alpha, Beta = Beta, DD = DD, sigma = sigma, Phi = Phi, nu = nu)

set.seed(20250301)
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

k1 <- ggplot(data = Data[which(Data$R != 1), ], aes(x = Time, y = Var1, col = as.factor(treat), shape = as.factor(treat), grop = as.factor(Subject))) +
  geom_line() +
  geom_point(size = 3) +
  xlab("Time") +
  ylab(expression(paste("Measurements of ", y[ij]^"o"))) +
  theme(legend.position = "top") +
  scale_color_discrete(labels = c("Placebo", "Therapy")) +
  scale_shape_discrete(labels = c("Placebo", "Therapy")) +
  theme(text = element_text(size = 20)) +
  theme(
    legend.key.size = unit(1, "cm"),
    legend.key.width = unit(1.5, "cm")
  ) +
  theme(legend.title = element_blank())
k1

print(k1)

source(paste(PATH, "/function/multiplot.R", sep = ""))
layout <- matrix(c(1), nrow = 1, byrow = TRUE)
multiplot(plotlist = list(k1), layout = layout)


pdf(paste0(PATH, "/Result/FigureC1.pdf"), width = 12, height = 8)
multiplot(plotlist = list(k1), layout = layout)
dev.off()