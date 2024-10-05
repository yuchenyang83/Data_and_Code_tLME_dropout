load(paste0(PATH, "/Data/fit_result.RData"))
source(paste0(PATH, "/function/computer.pvi.R"))

est.MNAR.pvi = compute.pvi(Data, X, Z, V, g=g, init.para = est.MNAR.ARp$para.est, fit = est.MNAR.ARp,
                           cor.type = c("ARp"), M = 100, M.LL = 1000, P=1, tol = 1e-6, 
                           max.iter=max.iter, per=1, mechanism=c('MNAR'))

est.MAR.pvi = compute.pvi(Data, X, Z, V, g=g, init.para = est.MAR.ARp$para.est, fit = est.MAR.ARp,
                          cor.type = c("ARp"), M = 100, M.LL = 1000, P=1, tol = 1e-6, 
                          max.iter=max.iter, per=1, mechanism=c('MAR'))

est.MCAR.pvi = compute.pvi(Data, X, Z, V, g=g, init.para = est.MCAR.ARp$para.est, fit = est.MCAR.ARp, 
                           cor.type = c("ARp"), M = 100, M.LL = 1000, P=1, tol = 1e-6, 
                           max.iter=max.iter, per=1, mechanism=c('MCAR'))

est.MNAR.roc=roc(Data$R, est.MNAR.pvi,plot=T,col=3,cex.lab=1.2,main="ROC",
                 cex.main=1.5,cex.axis=1.5,lwd=1, print.auc=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, print.thres=TRUE)

est.MAR.roc=roc(Data$R, est.MAR.pvi,plot=T,col=3,cex.lab=1.2,main="ROC",
                cex.main=1.5,cex.axis=1.5,lwd=1, print.auc=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, print.thres=TRUE)

est.MCAR.roc=roc(Data$R, est.MCAR.pvi,plot=T,col=3,cex.lab=1.2,main="ROC",
                 cex.main=1.5,cex.axis=1.5,lwd=1, print.auc=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, print.thres=TRUE)


ya1 = ggroc(list(MCAR = est.MCAR.roc, MAR = est.MAR.roc, MNAR = est.MNAR.roc), aes=c("linetype", "color"), size = 0.8) + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +
  guides(shape=guide_legend(title=NULL), 
         colour=guide_legend(title=NULL), 
         linetype=guide_legend(title=NULL)) +
  scale_linetype_manual(values=c(1, 1, 1), breaks=c("MCAR", "MAR", "MNAR"),
                        labels=c("tLME-MCAR-ARp (AUC=0.500)", "tLME-MAR-ARp (AUC=0.580)", "tLME-MNAR-ARp (AUC=0.985)")) +
  scale_color_manual(values = c("#F8766D","#00BA38","#619CFF"), 
                     breaks=c("MCAR", "MAR", "MNAR"),
                     labels=c("tLME-MCAR-ARp (AUC=0.500)", "tLME-MAR-ARp (AUC=0.580)", "tLME-MNAR-ARp (AUC=0.985)")) +
  labs(title = 'ROC curve') +
  theme(legend.position=c(.67, .18), 
        text = element_text(size=30),
        plot.title = element_text(hjust = 0.5))+
  xlab("Specificity") + ylab("Sensitivity")
print(ya1)

