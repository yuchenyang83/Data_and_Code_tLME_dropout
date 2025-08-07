################################################################################
#                                                                                     
#   Filename    :    FigC3.R  												  
#   Project     :    BiomJ article "Extending t linear mixed models for longitudinal 
#                    data with non-ignorable dropout applied to AIDS studies"                                                           
#   Authors     :    Yu-Chen Yang and Wan-Lun Wang and Luis M. Castro and Tsung-I Lin
#   Date        :    07.08.2025
#   Purpose     :    produce Figure C3 in Appendix C
#
#   Input data files  :  None
#   Output data files :  Data_and_Code/results/FigureC3.pdf
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

k3 <- ggplot(data = Data[-which(is.na(Data$Var1)), ], aes(x = as.factor(Time), fill = as.factor(treat))) +
  geom_bar(position = position_dodge()) +
  xlab("Time") +
  ylab("Number of observed responses") +
  geom_text(
    stat = "count", aes(label = paste0(round(after_stat(count) / n, 3) * 100, "%")), vjust = 1.6, color = "black", size = 5,
    position = position_dodge(0.9)
  ) +
  theme(legend.position = "top") +
  scale_fill_discrete(name = NULL, labels = c("Placebo", "Therapy")) +
  theme(text = element_text(size = 20))
k3

print(k3)

source(paste(PATH, "/function/multiplot.R", sep = ""))
layout <- matrix(c(1), nrow = 1, byrow = TRUE)
multiplot(plotlist = list(k3), layout = layout)


pdf(paste0(PATH, "/Result/FigureC3.pdf"), width = 12, height = 8)
multiplot(plotlist = list(k3), layout = layout)
dev.off()
