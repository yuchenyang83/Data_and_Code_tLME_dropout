load(paste0(PATH, "/Data/fit_result.RData"))

matrix(c(round(c(rbind(fit.MCAR$IM$out, abs(fit.MCAR$IM$out[1,]/fit.MCAR$IM$out[2,])), rep(NA, 6*3)), 3),
round(c(rbind(fit.MAR$IM$out, abs(fit.MAR$IM$out[1,]/fit.MAR$IM$out[2,])), rep(NA, 3)), 3),
round(c(rbind(fit.MNAR$IM$out, abs(fit.MNAR$IM$out[1,]/fit.MNAR$IM$out[2,]))), 3)), ncol=3)


xxx = matrix(c(round(c(rbind(est.MCAR.ARp$IM$out, abs(est.MCAR.ARp$IM$out[1,]/est.MCAR.ARp$IM$out[2,])), rep(NA, 3*3)), 4),
         round(c(rbind(est.MAR.ARp$IM$out, abs(est.MAR.ARp$IM$out[1,]/est.MAR.ARp$IM$out[2,])), rep(NA, 3)), 4),
         round(c(rbind(est.MNAR.ARp$IM$out, abs(est.MNAR.ARp$IM$out[1,]/est.MNAR.ARp$IM$out[2,]))), 4)), byrow = T, ncol=9)


print("MCAR")
print(round(t(cbind(rbind(est.MCAR.ARp$IM$out, abs(est.MCAR.ARp$IM$out[1,]/est.MCAR.ARp$IM$out[2,])), NA, NA, NA)), 3))
print("MAR")
print(round(t(cbind(rbind(est.MAR.ARp$IM$out, abs(est.MAR.ARp$IM$out[1,]/est.MAR.ARp$IM$out[2,])), NA)), 3))
print("MNAR")
print(round(t(cbind(rbind(est.MNAR.ARp$IM$out, abs(est.MNAR.ARp$IM$out[1,]/est.MNAR.ARp$IM$out[2,])))), 3))
