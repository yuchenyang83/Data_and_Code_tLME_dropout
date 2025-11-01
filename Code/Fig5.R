################################################################################
#                                                                                     
#   Filename    :    Fig5.R  												  
#   Project     :    BiomJ article "Extending t linear mixed models for longitudinal 
#                    data with non-ignorable dropout applied to AIDS studies"                                                           
#   Authors     :    Yu-Chen Yang and Wan-Lun Wang and Luis M. Castro and Tsung-I Lin
#   Date        :    07.08.2025
#   Purpose     :    produce Figure 5 for AIDS data
#
#   Input data files  :  Data_and_Code/Data/fit_result.RData
#   Output data files :  Data_and_Code/results/Figure5.pdf
#
#   R Version   :    R-4.3.1                                                              
#   Required R packages : ggplot2; ggtext; pROC; rlang; cowplot
#
################################################################################ 
load(paste0(PATH, "/Data/fit_result.RData"))
source(paste0(PATH, "/function/compute.pvi.R"))

est.MNAR.pvi <- compute.pvi(Data, X, Z, 
  init.para = est.MNAR.ARp$para.est, fit = est.MNAR.ARp,
  cor.type = c("ARp"), mechanism = c("MNAR")
)

est.MAR.pvi <- compute.pvi(Data, X, Z,
  init.para = est.MAR.ARp$para.est, fit = est.MAR.ARp,
  cor.type = c("ARp"), mechanism = c("MAR")
)

est.MCAR.pvi <- compute.pvi(Data, X, Z,
  init.para = est.MCAR.ARp$para.est, fit = est.MCAR.ARp,
  cor.type = c("ARp"), mechanism = c("MCAR")
)

est.MNAR.roc <- roc(Data$R, est.MNAR.pvi,
  plot = T, col = 3, cex.lab = 1.2, main = "ROC",
  cex.main = 1.5, cex.axis = 1.5, lwd = 1, print.auc = TRUE, auc.polygon = TRUE, max.auc.polygon = TRUE, print.thres = TRUE, levels = c(0, 1), direction = "<"
)

est.MAR.roc <- roc(Data$R, est.MAR.pvi,
  plot = T, col = 3, cex.lab = 1.2, main = "ROC",
  cex.main = 1.5, cex.axis = 1.5, lwd = 1, print.auc = TRUE, auc.polygon = TRUE, max.auc.polygon = TRUE, print.thres = TRUE, levels = c(0, 1), direction = "<"
)

est.MCAR.roc <- roc(Data$R, est.MCAR.pvi,
  plot = T, col = 3, cex.lab = 1.2, main = "ROC",
  cex.main = 1.5, cex.axis = 1.5, lwd = 1, print.auc = TRUE, auc.polygon = TRUE, max.auc.polygon = TRUE, print.thres = TRUE, levels = c(0, 1), direction = "<"
)


ya1 <- ggroc(list(MCAR = est.MCAR.roc, MAR = est.MAR.roc, MNAR = est.MNAR.roc), aes = c("linetype", "color"), size = 0.8) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color = "grey", linetype = "dashed") +
  guides(
    shape = guide_legend(title = NULL),
    colour = guide_legend(title = NULL),
    linetype = guide_legend(title = NULL)
  ) +
  scale_linetype_manual(
    values = c(1, 1, 1), breaks = c("MCAR", "MAR", "MNAR"),
    labels = c("tLME-MCAR-AR(1) (AUC=0.500)", "tLME-MAR-AR(1) (AUC=0.580)", "tLME-MNAR-AR(1) (AUC=0.985)")
  ) +
  scale_color_manual(
    values = c("#F8766D", "#00BA38", "#619CFF"),
    breaks = c("MCAR", "MAR", "MNAR"),
    labels = c("tLME-MCAR-AR(1) (AUC=0.500)", "tLME-MAR-AR(1) (AUC=0.580)", "tLME-MNAR-AR(1) (AUC=0.985)")
  ) +
  labs(title = "ROC curve") +
  theme(
    legend.position = c(.67, .18),
    text = element_text(size = 30),
    plot.title = element_text(hjust = 0.5)
  ) +
  xlab("Specificity") + ylab("Sensitivity")
print(ya1)

source(paste(PATH, "/function/multiplot.R", sep = ""))
layout <- matrix(c(1), nrow = 1, byrow = TRUE)
multiplot(plotlist = list(ya1), layout = layout)


pdf(paste0(PATH, "/Result/Figure5.pdf"), width = 12, height = 8)
multiplot(plotlist = list(ya1), layout = layout)
dev.off()
