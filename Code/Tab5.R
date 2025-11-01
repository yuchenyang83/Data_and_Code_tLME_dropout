################################################################################
#                                                                                     
#   Filename    :    Tab5.R  												  
#   Project     :    BiomJ article "Extending t linear mixed models for longitudinal 
#                    data with non-ignorable dropout applied to AIDS studies"                                                           
#   Authors     :    Yu-Chen Yang and Wan-Lun Wang and Luis M. Castro and Tsung-I Lin
#   Date        :    07.08.2025
#   Purpose     :    produce Table 5 for AIDS
#
#   Input data files  :  Data_and_Code/Data/fit_result.RData
#   Output data files :  Data_and_Code/results/Table5.csv
#
#   R Version   :    R-4.3.1                                                              
#   Required R packages : None
#
################################################################################ 
load(paste0(PATH, "/Data/fit_result.RData"))

matrix(c(
  round(c(rbind(fit.MCAR$IM$out, abs(fit.MCAR$IM$out[1, ] / fit.MCAR$IM$out[2, ])), rep(NA, 6 * 3)), 3),
  round(c(rbind(fit.MAR$IM$out, abs(fit.MAR$IM$out[1, ] / fit.MAR$IM$out[2, ])), rep(NA, 3)), 3),
  round(c(rbind(fit.MNAR$IM$out, abs(fit.MNAR$IM$out[1, ] / fit.MNAR$IM$out[2, ]))), 3)
), ncol = 3)


xxx <- matrix(c(
  round(c(rbind(est.MCAR.ARp$IM$out, abs(est.MCAR.ARp$IM$out[1, ] / est.MCAR.ARp$IM$out[2, ])), rep(NA, 3 * 3)), 4),
  round(c(rbind(est.MAR.ARp$IM$out, abs(est.MAR.ARp$IM$out[1, ] / est.MAR.ARp$IM$out[2, ])), rep(NA, 3)), 4),
  round(c(rbind(est.MNAR.ARp$IM$out, abs(est.MNAR.ARp$IM$out[1, ] / est.MNAR.ARp$IM$out[2, ]))), 4)
), byrow = T, ncol = 9)


# print("MCAR")
# print(round(t(cbind(rbind(est.MCAR.ARp$IM$out, abs(est.MCAR.ARp$IM$out[1, ] / est.MCAR.ARp$IM$out[2, ])), NA, NA, NA)), 3))
# print("MAR")
# print(round(t(cbind(rbind(est.MAR.ARp$IM$out, abs(est.MAR.ARp$IM$out[1, ] / est.MAR.ARp$IM$out[2, ])), NA)), 3))
# print("MNAR")
# print(round(t(cbind(rbind(est.MNAR.ARp$IM$out, abs(est.MNAR.ARp$IM$out[1, ] / est.MNAR.ARp$IM$out[2, ])))), 3))

MCAR.r <- round(t(cbind(rbind(est.MCAR.ARp$IM$out, abs(est.MCAR.ARp$IM$out[1, ] / est.MCAR.ARp$IM$out[2, ])), NA, NA, NA)), 3)
MAR.r <- round(t(cbind(rbind(est.MAR.ARp$IM$out, abs(est.MAR.ARp$IM$out[1, ] / est.MAR.ARp$IM$out[2, ])), NA)), 3)
MNAR.r <- round(t(cbind(rbind(est.MNAR.ARp$IM$out, abs(est.MNAR.ARp$IM$out[1, ] / est.MNAR.ARp$IM$out[2, ])))), 3)

sum.table <- cbind(MCAR.r, MAR.r, MNAR.r)

colnames(sum.table) <- c("Est", "SE", "Est/SE", "Est", "SE", "Est/SE", "Est", "SE", "Est/SE")
row.names(sum.table) <- c("beta0", "beta1", "beta2", "beta3", "beta4", "beta5", "d11", "d12", "d22", "sigma^2", "phi", "nu", "alpha00", "alpha01", "alpha1", "alpha2")


write.csv(sum.table, paste0(PATH, "/Result/Table5.csv"))
