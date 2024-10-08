load(paste0(PATH, "/Data/fit_result.RData"))

2 * dim(est.MCAR$IM$out)[2] - 2 * est.MCAR$model.inf$loglik
dim(est.MCAR$IM$out)[2] * log(N) - 2 * est.MCAR$model.inf$loglik

c(dim(fit.MCAR$IM$out)[2], fit.MCAR$model.inf$loglik, fit.MCAR$model.inf$aic, fit.MCAR$model.inf$bic)
c(dim(fit.MAR$IM$out)[2], fit.MAR$model.inf$loglik, fit.MAR$model.inf$aic, fit.MAR$model.inf$bic)
c(dim(fit.MNAR$IM$out)[2], fit.MNAR$model.inf$loglik, fit.MNAR$model.inf$aic, fit.MNAR$model.inf$bic)

c(dim(est.MCAR$IM$out)[2], est.MCAR$model.inf$loglik, est.MCAR$model.inf$aic, est.MCAR$model.inf$bic)
c(dim(est.MAR$IM$out)[2], est.MAR$model.inf$loglik, est.MAR$model.inf$aic, est.MAR$model.inf$bic)
c(dim(est.MNAR$IM$out)[2], est.MNAR$model.inf$loglik, est.MNAR$model.inf$aic, est.MNAR$model.inf$bic)

c(dim(fit.MCAR.BAND1$IM$out)[2], fit.MCAR.BAND1$model.inf$loglik, fit.MCAR.BAND1$model.inf$aic, fit.MCAR.BAND1$model.inf$bic)
c(dim(fit.MAR.BAND1$IM$out)[2], fit.MAR.BAND1$model.inf$loglik, fit.MAR.BAND1$model.inf$aic, fit.MAR.BAND1$model.inf$bic)
c(dim(fit.MNAR.BAND1$IM$out)[2], fit.MNAR.BAND1$model.inf$loglik, fit.MNAR.BAND1$model.inf$aic, fit.MNAR.BAND1$model.inf$bic)

c(dim(est.MCAR.BAND1$IM$out)[2], est.MCAR.BAND1$model.inf$loglik, est.MCAR.BAND1$model.inf$aic, est.MCAR.BAND1$model.inf$bic)
c(dim(est.MAR.BAND1$IM$out)[2], est.MAR.BAND1$model.inf$loglik, est.MAR.BAND1$model.inf$aic, est.MAR.BAND1$model.inf$bic)
c(dim(est.MNAR.BAND1$IM$out)[2], est.MNAR.BAND1$model.inf$loglik, est.MNAR.BAND1$model.inf$aic, est.MNAR.BAND1$model.inf$bic)

c(dim(fit.MCAR.ARp$IM$out)[2], fit.MCAR.ARp$model.inf$loglik, fit.MCAR.ARp$model.inf$aic, fit.MCAR.ARp$model.inf$bic)
c(dim(fit.MAR.ARp$IM$out)[2], fit.MAR.ARp$model.inf$loglik, fit.MAR.ARp$model.inf$aic, fit.MAR.ARp$model.inf$bic)
c(dim(fit.MNAR.ARp$IM$out)[2], fit.MNAR.ARp$model.inf$loglik, fit.MNAR.ARp$model.inf$aic, fit.MNAR.ARp$model.inf$bic)

c(dim(est.MCAR.ARp$IM$out)[2], est.MCAR.ARp$model.inf$loglik, est.MCAR.ARp$model.inf$aic, est.MCAR.ARp$model.inf$bic)
c(dim(est.MAR.ARp$IM$out)[2], est.MAR.ARp$model.inf$loglik, est.MAR.ARp$model.inf$aic, est.MAR.ARp$model.inf$bic)
c(dim(est.MNAR.ARp$IM$out)[2], est.MNAR.ARp$model.inf$loglik, est.MNAR.ARp$model.inf$aic, est.MNAR.ARp$model.inf$bic)




print(matrix(c(c(dim(fit.MCAR$IM$out)[2], fit.MCAR$model.inf$loglik, fit.MCAR$model.inf$aic, fit.MCAR$model.inf$bic),
         c(dim(fit.MAR$IM$out)[2], fit.MAR$model.inf$loglik, fit.MAR$model.inf$aic, fit.MAR$model.inf$bic),
         c(dim(fit.MNAR$IM$out)[2], fit.MNAR$model.inf$loglik, fit.MNAR$model.inf$aic, fit.MNAR$model.inf$bic),
         
         c(dim(est.MCAR$IM$out)[2], est.MCAR$model.inf$loglik, est.MCAR$model.inf$aic, est.MCAR$model.inf$bic),
         c(dim(est.MAR$IM$out)[2], est.MAR$model.inf$loglik, est.MAR$model.inf$aic, est.MAR$model.inf$bic),
         c(dim(est.MNAR$IM$out)[2], est.MNAR$model.inf$loglik, est.MNAR$model.inf$aic, est.MNAR$model.inf$bic),
         
         c(12, fit.MCAR.CS$model.inf$loglik, fit.MCAR.CS$model.inf$aic, fit.MCAR.CS$model.inf$bic),
         c(14, fit.MAR.CS$model.inf$loglik, fit.MAR.CS$model.inf$aic, fit.MAR.CS$model.inf$bic),
         c(15, fit.MNAR.CS$model.inf$loglik, fit.MNAR.CS$model.inf$aic, fit.MNAR.CS$model.inf$bic),
         
         c(13, est.MCAR.CS$model.inf$loglik, est.MCAR.CS$model.inf$aic, est.MCAR.CS$model.inf$bic),
         c(15, est.MAR.CS$model.inf$loglik, est.MAR.CS$model.inf$aic, est.MAR.CS$model.inf$bic),
         c(16, est.MNAR.CS$model.inf$loglik, est.MNAR.CS$model.inf$aic, est.MNAR.CS$model.inf$bic),
         
         c(dim(fit.MCAR.BAND1$IM$out)[2], fit.MCAR.BAND1$model.inf$loglik, fit.MCAR.BAND1$model.inf$aic, fit.MCAR.BAND1$model.inf$bic),
         c(dim(fit.MAR.BAND1$IM$out)[2], fit.MAR.BAND1$model.inf$loglik, fit.MAR.BAND1$model.inf$aic, fit.MAR.BAND1$model.inf$bic),
         c(dim(fit.MNAR.BAND1$IM$out)[2], fit.MNAR.BAND1$model.inf$loglik, fit.MNAR.BAND1$model.inf$aic, fit.MNAR.BAND1$model.inf$bic),
         
         c(dim(est.MCAR.BAND1$IM$out)[2], est.MCAR.BAND1$model.inf$loglik, est.MCAR.BAND1$model.inf$aic, est.MCAR.BAND1$model.inf$bic),
         c(dim(est.MAR.BAND1$IM$out)[2], est.MAR.BAND1$model.inf$loglik, est.MAR.BAND1$model.inf$aic, est.MAR.BAND1$model.inf$bic),
         c(dim(est.MNAR.BAND1$IM$out)[2], est.MNAR.BAND1$model.inf$loglik, est.MNAR.BAND1$model.inf$aic, est.MNAR.BAND1$model.inf$bic),
         
         c(dim(fit.MCAR.ARp$IM$out)[2], fit.MCAR.ARp$model.inf$loglik, fit.MCAR.ARp$model.inf$aic, fit.MCAR.ARp$model.inf$bic),
         c(dim(fit.MAR.ARp$IM$out)[2], fit.MAR.ARp$model.inf$loglik, fit.MAR.ARp$model.inf$aic, fit.MAR.ARp$model.inf$bic),
         c(dim(fit.MNAR.ARp$IM$out)[2], fit.MNAR.ARp$model.inf$loglik, fit.MNAR.ARp$model.inf$aic, fit.MNAR.ARp$model.inf$bic),
         
         c(dim(est.MCAR.ARp$IM$out)[2], est.MCAR.ARp$model.inf$loglik, est.MCAR.ARp$model.inf$aic, est.MCAR.ARp$model.inf$bic),
         c(dim(est.MAR.ARp$IM$out)[2], est.MAR.ARp$model.inf$loglik, est.MAR.ARp$model.inf$aic, est.MAR.ARp$model.inf$bic),
         c(dim(est.MNAR.ARp$IM$out)[2], est.MNAR.ARp$model.inf$loglik, est.MNAR.ARp$model.inf$aic, est.MNAR.ARp$model.inf$bic)), nrow = 12))




